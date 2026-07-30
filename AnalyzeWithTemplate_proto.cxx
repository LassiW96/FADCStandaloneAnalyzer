// AnalyzeWithTemplate_proto.cxx
//
// Prototype implementation of a matched-filter + template-fit waveform
// analysis approach, adapted from the NPS calorimeter analysis in fit.C
// (Wassim Hamdi & Malek Mazouz) for the FADCStandaloneAnalyzer class structure.
//
// NOTHING here is wired into FADCStandaloneAnalyzer yet.
// This file is intentionally self-contained.
//
// Algorithm overview
// ──────────────────
//  1. Estimate pedestal from the first N samples.
//  2. Cross-correlate the pedestal-subtracted waveform with the reversed
//     central window of a normalised reference pulse (matched filter).
//  3. Find peaks in the MF output above a threshold (two complementary
//     strategies mirroring fit.C's npulsewf2).
//  4. For each MF peak, scan a small window of time offsets and solve for
//     amplitude + baseline analytically (closed-form linear least squares
//     at each candidate t0, then pick the minimum chi2).
//  5. Subtract each fitted pulse from the residual before fitting the next
//     (sequential peeling), so pile-up pulses don't bias each other.
//
// Correspondence to fit.C
// ──────────────────────────────────────────────────────────────────────────
//  fit.C concept                 │ This file
//  ──────────────────────────────┼──────────────────────────────────────────
//  ref_wf_N.txt                  │ FADCTemplate::fromWaveformParams()
//  mfyref[][], mfint[]           │ FADCTemplate::mfKernel, mfNorm
//  interpolation[i]->Eval()      │ FADCTemplate::eval() (linear interp.)
//  hsig_mf[], output[]           │ computeMatchedFilter()
//  npulsewf2()                   │ findMFPeaks()
//  Fitwf2() / TF1 + Minuit       │ fitSingleTemplate() (analytic scan)
//  TEST_2() event loop           │ AnalyzeWithTemplate()

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <string>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"
#include "waveformGenerator.h"   // for WaveformParams (tau_r, tau_f)

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TString.h>


// ─────────────────────────────────────────────────────────────────────────────
// Tunable parameters
// ─────────────────────────────────────────────────────────────────────────────

struct TemplateAnalysisConfig {
    // Pedestal
    Int_t    fNPedSamples  = 4;      // samples used to estimate the pedestal

    // Matched-filter window  (mirrors mfstart/mfend/mfleft/mfright in fit.C)
    Int_t    fMFLeft       = 5;      // samples left  of template peak in MF kernel
    Int_t    fMFRight      = 5;      // samples right of template peak in MF kernel
    Int_t    fMFStart      = 10;     // first sample of the MF search window
    Int_t    fMFEnd        = 100;    // last  sample of the MF search window

    // Peak finding  (mirrors specthres, mfthres in fit.C)
    Double_t fMFThreshold  = 0.5;   // MF amplitude threshold (normalised units)
    Int_t    fMFSeparation = 3;     // minimum distance between two MF peaks (samples)
    Int_t    fMaxPulses    = 12;    // cap on the number of pulses to fit

    // Template fit scan  (mirrors the +-4-bin window in Fitwf2's SetParLimits)
    Double_t fScanHalf     = 4.0;   // +- scan range around the MF candidate (samples)
    Double_t fScanStep     = 0.25;  // step size of the t0 scan (samples)

    // Noise model  (mirrors the sqrt(|content|*4.096/2)/4.096 in fit.C line 895)
    Double_t fAdcPerMV     = 4.096; // ADC counts per mV
    Double_t fMinSigma     = 0.5;   // floor for per-sample sigma (ADC counts)

    // Amplitude quality cut
    Double_t fAmpCutADC    = 1.0;   // discard fitted pulses below this amplitude
};


// ─────────────────────────────────────────────────────────────────────────────
// FADCTemplate
// ─────────────────────────────────────────────────────────────────────────────
// Stores the normalised reference pulse shape and the derived matched-filter
// kernel.  The shape is evaluated at integer sample positions; sub-sample
// values are obtained via linear interpolation (eval).
//
// In fit.C this comes from an external ref_wf_N.txt file and a ROOT cubic
// spline (ROOT::Math::Interpolator).  Here we build it directly from the
// WaveformParams already present in the project, so no file I/O is needed.
// Swap fromWaveformParams() for a file-reader if you need to load NPS files.

struct FADCTemplate {

    std::vector<Double_t> shape;      // normalised pulse shape, shape[i] in [0,1]
    std::vector<Double_t> mfKernel;   // reversed window used in the matched filter
    Double_t              mfNorm;     // sum(mfKernel) used as the normalisation divisor
    Int_t                 peakBin;    // sample index of the template peak (shape=1)

    // Build template from WaveformParams (SPE bi-exponential pulse shape).
    // Parameters:
    //   wp          - WaveformParams (provides tau_r and tau_f)
    //   sampleRate  - ns per sample (e.g. 4.0 for the NPS fADC)
    //   nSamples    - template length in samples (e.g. 110)
    //   peakTime_ns - where to place the pulse peak within the window (ns)
    //   mfLeft/Right- half-widths of the MF kernel in samples
    static FADCTemplate fromWaveformParams(const WaveformParams& wp,
                                           Double_t sampleRate  = 4.0,
                                           Int_t    nSamples    = 110,
                                           Double_t peakTime_ns = 20.0,
                                           Int_t    mfLeft      = 5,
                                           Int_t    mfRight     = 5);

    // Sub-sample evaluation via linear interpolation.
    // Returns 0 for t outside [0, shape.size()-2].
    Double_t eval(Double_t t) const;
};

// SPE pulse shape (bi-exponential), mirrors waveformGenerator::spePulse
static Double_t _spePulse(Double_t t, Double_t t0,
                           Double_t tau_r, Double_t tau_f) {
    if (t < t0) return 0.0;
    return std::exp(-(t - t0) / tau_f) - std::exp(-(t - t0) / tau_r);
}

FADCTemplate FADCTemplate::fromWaveformParams(const WaveformParams& wp,
                                               Double_t sampleRate,
                                               Int_t    nSamples,
                                               Double_t peakTime_ns,
                                               Int_t    mfLeft,
                                               Int_t    mfRight) {
    FADCTemplate tmpl;
    tmpl.shape.resize(nSamples, 0.0);

    // Evaluate the raw SPE shape
    for (Int_t i = 0; i < nSamples; ++i)
        tmpl.shape[i] = _spePulse(i * sampleRate, peakTime_ns, wp.tau_r, wp.tau_f);

    // Normalise to peak = 1 and record the peak bin
    Double_t peak = *std::max_element(tmpl.shape.begin(), tmpl.shape.end());
    if (peak <= 0.0)
        throw std::runtime_error("FADCTemplate: zero-amplitude shape -- check WaveformParams");

    tmpl.peakBin = -1;
    for (Int_t i = 0; i < nSamples; ++i) {
        tmpl.shape[i] /= peak;
        if (tmpl.shape[i] >= 1.0 - 1e-9) tmpl.peakBin = i;
    }

    // Build the matched-filter kernel:
    //   Take the window [peakBin - mfLeft, peakBin + mfRight] of the shape
    //   and reverse it.  Reversal turns cross-correlation into convolution,
    //   so a simple forward loop suffices in computeMatchedFilter().
    //   (Mirrors the mfyref[][] construction in fit.C lines 337-338.)
    tmpl.mfKernel.clear();
    tmpl.mfNorm = 0.0;
    for (Int_t i = tmpl.peakBin + mfRight; i >= tmpl.peakBin - mfLeft; --i) {
        Double_t v = (i >= 0 && i < nSamples) ? tmpl.shape[i] : 0.0;
        tmpl.mfKernel.push_back(v);
        tmpl.mfNorm += v;
    }
    if (tmpl.mfNorm <= 0.0)
        throw std::runtime_error("FADCTemplate: zero MF norm");

    return tmpl;
}

Double_t FADCTemplate::eval(Double_t t) const {
    if (t < 0.0 || t >= static_cast<Double_t>(shape.size() - 1))
        return 0.0;
    Int_t    lo   = static_cast<Int_t>(t);
    Double_t frac = t - lo;
    return shape[lo] * (1.0 - frac) + shape[lo + 1] * frac;
}


// ─────────────────────────────────────────────────────────────────────────────
// Output structs
// ─────────────────────────────────────────────────────────────────────────────

struct FitPulseResult {
    Double_t timeSample;    // fitted t0 in sample units
    Double_t amplitude;     // fitted amplitude A (pedestal-subtracted ADC counts)
    Double_t baseline;      // fitted baseline offset b (residual DC after ped sub)
    Double_t chi2ndf;       // chi2 / NDF at the best t0
    Double_t mfAmplitude;   // MF output value at the candidate peak
    Bool_t   valid;         // kFALSE if rejected (A < fAmpCutADC, det~0, etc.)
};

struct FADCPulseWF {
    Double_t                    pedestal;  // mean of the first fNPedSamples
    Int_t                       nPulses;   // number of valid fitted pulses
    std::vector<FitPulseResult> pulses;    // one entry per MF peak, sorted by time
};


// ─────────────────────────────────────────────────────────────────────────────
// Step 1 -- computeMatchedFilter
// ─────────────────────────────────────────────────────────────────────────────
// Cross-correlates the pedestal-subtracted waveform with the (reversed) MF
// kernel, then shifts the output so its minimum is zero.
//
// mfOut[i] = sum_j  (wfPS[i + j - mfRight] - wfMin) * mfKernel[j] / mfNorm
//
// The wfMin subtraction mirrors the `minsignal[bn]` removal in fit.C line 902:
//   output[it] += (hsig_i[i]->GetBinContent(...) - minsignal[i]) * mfyref[...] / mfint[i]
// It prevents a negative baseline from suppressing real peaks in the MF output.

static std::vector<Double_t> computeMatchedFilter(
        const std::vector<Double_t>& wfPS,
        const FADCTemplate&          tmpl,
        const TemplateAnalysisConfig& cfg)
{
    Int_t N  = static_cast<Int_t>(wfPS.size());
    Int_t kW = static_cast<Int_t>(tmpl.mfKernel.size()); // mfLeft + mfRight + 1

    Double_t wfMin = *std::min_element(wfPS.begin(), wfPS.end());

    std::vector<Double_t> mfOut(N, 0.0);
    Double_t mfOutMin = std::numeric_limits<Double_t>::max();

    for (Int_t i = cfg.fMFLeft; i < N - cfg.fMFRight; ++i) {
        Double_t val = 0.0;
        for (Int_t j = 0; j < kW; ++j) {
            Int_t wIdx = i + j - cfg.fMFRight;
            if (wIdx >= 0 && wIdx < N)
                val += (wfPS[wIdx] - wfMin) * tmpl.mfKernel[j];
        }
        mfOut[i] = val / tmpl.mfNorm;
        if (mfOut[i] < mfOutMin) mfOutMin = mfOut[i];
    }

    // Shift so the global minimum of the MF output = 0  (mirrors fit.C line 905)
    for (Int_t i = cfg.fMFLeft; i < N - cfg.fMFRight; ++i)
        mfOut[i] -= mfOutMin;

    return mfOut;
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 2 -- findMFPeaks
// ─────────────────────────────────────────────────────────────────────────────
// Finds pulse candidate positions in the MF output using two complementary
// strategies, mirroring npulsewf2() in fit.C:
//
//  Primary   -- simple local maximum: mfOut[i] > both neighbours, above threshold.
//               Sub-sample position is refined via parabolic interpolation.
//  Secondary -- "4 rising / 2 falling" pattern (fit.C lines 160-168).
//               Catches broad peaks where the strict local-max condition fails.
//               A secondary candidate is suppressed if it is within
//               fMFSeparation samples of an already-recorded primary peak.

struct MFPeak {
    Double_t position;    // sample-unit position (may be fractional)
    Double_t mfAmplitude; // MF output value at the peak
};

static std::vector<MFPeak> findMFPeaks(
        const std::vector<Double_t>& mfOut,
        const TemplateAnalysisConfig& cfg)
{
    Int_t N     = static_cast<Int_t>(mfOut.size());
    Int_t mfEnd = std::min(cfg.fMFEnd, N - 3);

    std::vector<MFPeak> peaks;

    auto isDuplicate = [&](Double_t pos) -> Bool_t {
        for (const auto& p : peaks)
            if (std::abs(p.position - pos) < static_cast<Double_t>(cfg.fMFSeparation))
                return kTRUE;
        return kFALSE;
    };

    // -- Primary: strict local maximum ----------------------------------------
    for (Int_t i = cfg.fMFStart + 1; i < mfEnd - 1; ++i) {
        if (mfOut[i] > cfg.fMFThreshold &&
            mfOut[i] > mfOut[i - 1]     &&
            mfOut[i] > mfOut[i + 1]) {

            // Parabolic sub-sample interpolation of the peak position:
            //   pos = i - (y[i+1] - y[i-1]) / (2 * (y[i-1] - 2*y[i] + y[i+1]))
            Double_t denom = 2.0 * (mfOut[i - 1] - 2.0 * mfOut[i] + mfOut[i + 1]);
            Double_t pos   = (std::abs(denom) > 1e-9)
                             ? i - (mfOut[i + 1] - mfOut[i - 1]) / denom
                             : static_cast<Double_t>(i);

            if (!isDuplicate(pos))
                peaks.push_back({pos, mfOut[i]});
        }
    }

    // -- Secondary: 4-rising / 2-falling pattern  (mirrors fit.C lines 160-168) -
    for (Int_t i = cfg.fMFStart + 3; i < mfEnd - 2; ++i) {
        if (mfOut[i]     > cfg.fMFThreshold  &&
            mfOut[i]     > mfOut[i - 1]      &&   //  --|
            mfOut[i - 1] > mfOut[i - 2]      &&   //    | 4 rising samples
            mfOut[i - 2] > mfOut[i - 3]      &&   //    |
            mfOut[i - 3] > 0.0               &&   //  --|
            mfOut[i + 1] < mfOut[i]          &&   //  --| 2 falling samples
            mfOut[i + 2] < mfOut[i + 1]) {        //  --|

            Double_t pos = static_cast<Double_t>(i);
            if (!isDuplicate(pos))
                peaks.push_back({pos, mfOut[i]});
        }
    }

    std::sort(peaks.begin(), peaks.end(),
              [](const MFPeak& a, const MFPeak& b) {
                  return a.position < b.position; });

    return peaks;
}


// ─────────────────────────────────────────────────────────────────────────────
// Step 3 -- fitSingleTemplate
// ─────────────────────────────────────────────────────────────────────────────
// Fits the model:    wf[i] ≈ baseline + A * tmpl.eval(i - t0)
// to the waveform samples in [fitLo, fitHi].
//
// Key insight: A and baseline enter *linearly*.  For a fixed t0 the optimal
// (A, baseline) pair is the closed-form solution of the 2x2 normal equations:
//
//   | sum(1/s^2)    sum(f/s^2)  | | b |   | sum(y/s^2)  |
//   | sum(f/s^2)   sum(f^2/s^2)| | A | = | sum(yf/s^2) |
//
// where f_i = tmpl.eval(i - t0), y_i = wfResidual[i], s_i = noise sigma.
//
// We scan t0 over [candidateBin +- fScanHalf] at fScanStep resolution,
// solve analytically at each t0, and keep the minimum-chi2 result.
// This replaces TF1 + Minuit from Fitwf2() in fit.C.
//
// Noise model (mirrors fit.C line 895):
//   sigma_i = sqrt(|wfPS[i]| * adcPerMV / 2) / adcPerMV,
//   with a floor of fMinSigma to avoid division by zero on near-zero samples.

static FitPulseResult fitSingleTemplate(
        const std::vector<Double_t>& wfResidual,
        const FADCTemplate&          tmpl,
        Double_t                     candidateBin,
        Int_t                        fitLo,
        Int_t                        fitHi,
        const TemplateAnalysisConfig& cfg)
{
    Int_t N = static_cast<Int_t>(wfResidual.size());
    fitLo = std::max(fitLo, 0);
    fitHi = std::min(fitHi, N - 1);
    Int_t nPts = fitHi - fitLo + 1;

    // Per-sample sigma
    std::vector<Double_t> sigma(N, cfg.fMinSigma);
    for (Int_t i = fitLo; i <= fitHi; ++i) {
        Double_t s = std::sqrt(std::abs(wfResidual[i]) * cfg.fAdcPerMV / 2.0)
                     / cfg.fAdcPerMV;
        sigma[i] = (s > cfg.fMinSigma) ? s : cfg.fMinSigma;
    }

    FitPulseResult best;
    best.valid       = kFALSE;
    best.chi2ndf     = std::numeric_limits<Double_t>::max();
    best.timeSample  = candidateBin;
    best.amplitude   = 0.0;
    best.baseline    = 0.0;
    best.mfAmplitude = 0.0;

    // Scan t0 window
    for (Double_t t0 = candidateBin - cfg.fScanHalf;
                  t0 <= candidateBin + cfg.fScanHalf;
                  t0 += cfg.fScanStep) {

        // Accumulate sums for the normal equations
        Double_t S1  = 0.0, Sf  = 0.0, Sf2 = 0.0;
        Double_t Sy  = 0.0, Syf = 0.0;

        for (Int_t i = fitLo; i <= fitHi; ++i) {
            Double_t f     = tmpl.eval(static_cast<Double_t>(i) - t0);
            Double_t invs2 = 1.0 / (sigma[i] * sigma[i]);
            S1  += invs2;
            Sf  += f * invs2;
            Sf2 += f * f * invs2;
            Sy  += wfResidual[i] * invs2;
            Syf += wfResidual[i] * f * invs2;
        }

        // Solve the 2x2 system
        Double_t det = S1 * Sf2 - Sf * Sf;
        if (std::abs(det) < 1e-12) continue;

        Double_t b = (Sy * Sf2 - Syf * Sf) / det;   // baseline
        Double_t A = (S1 * Syf - Sf  * Sy) / det;   // amplitude

        if (A <= 0.0) continue;  // unphysical -- skip

        // chi2
        Double_t chi2 = 0.0;
        for (Int_t i = fitLo; i <= fitHi; ++i) {
            Double_t f   = tmpl.eval(static_cast<Double_t>(i) - t0);
            Double_t res = wfResidual[i] - b - A * f;
            chi2 += (res * res) / (sigma[i] * sigma[i]);
        }
        Int_t    ndf     = std::max(nPts - 2, 1);   // 2 free parameters
        Double_t chi2ndf = chi2 / ndf;

        if (chi2ndf < best.chi2ndf) {
            best.chi2ndf   = chi2ndf;
            best.timeSample   = t0;
            best.amplitude    = A;
            best.baseline     = b;
            best.valid        = kTRUE;
        }
    }

    return best;
}


// ─────────────────────────────────────────────────────────────────────────────
// AnalyzeWithTemplate  -- top-level function
// ─────────────────────────────────────────────────────────────────────────────
// Chains the four steps above.
//
// Pile-up handling via sequential subtraction ("peeling"):
//   After each accepted template fit the fitted pulse is subtracted from the
//   working residual so subsequent fits see cleaner waveforms.  This is the
//   single-channel equivalent of the sequential fit approach in Fitwf2().
//   The 3x3 neighbour-sum trigger (Apply_thres in fit.C) is a calorimeter-
//   specific addition and is not reproduced here.

FADCPulseWF AnalyzeWithTemplate(
        const std::vector<Double_t>& samples,
        const FADCTemplate&          tmpl,
        const TemplateAnalysisConfig& cfg = TemplateAnalysisConfig{})
{
    FADCPulseWF result;
    result.nPulses = 0;

    Int_t N = static_cast<Int_t>(samples.size());
    if (N == 0) return result;

    // 1. Pedestal
    Int_t    nPed = std::min(cfg.fNPedSamples, N);
    Double_t ped  = 0.0;
    Double_t ped2 = 0.0;
    for (Int_t i = 0; i < nPed; ++i) ped += samples[i];
    ped /= nPed;
    for (Int_t i = N - nPed; i < N; i++) ped2 += samples[i];
    ped2 /= nPed;

    if(ped2 < ped) {result.pedestal = ped2;}
    else {result.pedestal = ped;}

    std::vector<Double_t> wfPS(N);
    for (Int_t i = 0; i < N; ++i) wfPS[i] = samples[i] - ped;

    // 2. Matched filter
    std::vector<Double_t> mfOut = computeMatchedFilter(wfPS, tmpl, cfg);

    // 3. Peak finding
    std::vector<MFPeak> mfPeaks = findMFPeaks(mfOut, cfg);
    if (static_cast<Int_t>(mfPeaks.size()) > cfg.fMaxPulses)
        mfPeaks.resize(cfg.fMaxPulses);

    // 4. Template fit with sequential subtraction
    std::vector<Double_t> wfResidual = wfPS;

    for (const auto& peak : mfPeaks) {

        FitPulseResult fit = fitSingleTemplate(
                wfResidual, tmpl,
                peak.position,
                0, N - 1,
                cfg);

        fit.mfAmplitude = peak.mfAmplitude;

        if (!fit.valid || fit.amplitude < cfg.fAmpCutADC) {
            fit.valid = kFALSE;
            result.pulses.push_back(fit);
            continue;
        }

        // Subtract the fitted pulse so the next iteration sees a cleaner residual
        for (Int_t i = 0; i < N; ++i) {
            Double_t f = tmpl.eval(static_cast<Double_t>(i) - fit.timeSample);
            wfResidual[i] -= fit.baseline + fit.amplitude * f;
        }

        result.pulses.push_back(fit);
        result.nPulses++;
    }

    // Final sort by time (already ordered from MF sort, but make it explicit)
    std::sort(result.pulses.begin(), result.pulses.end(),
              [](const FitPulseResult& a, const FitPulseResult& b) {
                  return a.timeSample < b.timeSample; });

    return result;
}

// Plotting function
void plotWaveformFit(const std::vector<Double_t>& samples,
                     const FADCTemplate& tmpl,
                     const FADCPulseWF& result) 
{
    Int_t N = static_cast<Int_t>(samples.size());
    if (N == 0) return;

    // 1. Setup the canvas
    TCanvas* c1 = new TCanvas("c1", "Waveform Fit", 800, 600);

    // 2. Plot the Raw Waveform as a Histogram
    TH1D* hWaveform = new TH1D("hWaveform", "Template Fit;Sample Number;ADC Counts", N, -0.5, N - 0.5);
    for (Int_t i = 0; i < N; ++i) {
        hWaveform->SetBinContent(i + 1, samples[i]);
    }
    hWaveform->SetLineColor(kBlack);
    hWaveform->SetLineWidth(2);
    hWaveform->SetStats(0); // Turn off the stat box for a cleaner look
    hWaveform->Draw("HIST");

    // 3. Create a TF1 for the TOTAL Fit (Pedestal + All Valid Pulses)
    // We capture 'tmpl' and 'result' by value [tmpl, result] so ROOT doesn't crash 
    // if the function redraws after leaving this scope.
    TF1* fTotal = new TF1("fTotal", 
        [tmpl, result](double* x, double* /*p*/) {
            double t = x[0];
            double val = result.pedestal;
            
            for (const auto& pulse : result.pulses) {
                if (pulse.valid) {
                    val += pulse.baseline + pulse.amplitude * tmpl.eval(t - pulse.timeSample);
                }
            }
            return val;
        }, -0.5, N - 0.5, 0); // 0 parameters because we use captures

    fTotal->SetLineColor(kRed);
    fTotal->SetLineWidth(2);
    fTotal->SetNpx(1000); // High resolution for smooth curves
    fTotal->Draw("SAME");

    // 4. (Optional) Draw each individual pulse separately in dashed blue
    int pIdx = 0;
    for (const auto& pulse : result.pulses) {
        if (!pulse.valid) continue;

        TF1* fSingle = new TF1(Form("fSingle_%d", pIdx++), 
            [tmpl, result, pulse](double* x, double* /*p*/) {
                double t = x[0];
                // An individual pulse sits on top of the pedestal + its own baseline
                return result.pedestal + pulse.baseline + pulse.amplitude * tmpl.eval(t - pulse.timeSample);
            }, -0.5, N - 0.5, 0);

        fSingle->SetLineColor(kBlue);
        fSingle->SetLineStyle(2); // Dashed line
        fSingle->SetLineWidth(2);
        fSingle->SetNpx(1000);
        fSingle->Draw("SAME");

        c1->Modified();
        c1->Update();
    }

    // Save the canvas to a file so it can be inspected without blocking.
    c1->SaveAs("waveform_fit.png");
    std::cout << "[plotWaveformFit] Canvas saved to waveform_fit.png\n";
}


// ─────────────────────────────────────────────────────────────────────────────
// Quick smoke-test  (compile with -DRUN_PROTO_TEST to enable main())
// ─────────────────────────────────────────────────────────────────────────────
#ifdef RUN_PROTO_TEST
#include <cassert>

int main() {
    TApplication app("app", nullptr, nullptr);

    WaveformParams wp;   // default: tau_r=1 ns, tau_f=3 ns
    FADCTemplate tmpl = FADCTemplate::fromWaveformParams(wp);

    std::cout << "Template peak bin : " << tmpl.peakBin          << "\n";
    std::cout << "MF kernel size    : " << tmpl.mfKernel.size()  << "\n";
    std::cout << "MF norm           : " << tmpl.mfNorm           << "\n";

    // Generate a synthetic waveform via waveformGenerator
    DigitizerParams dq;
    dq.num_samples = 110;
    dq.sample_rate = 4.0;
    dq.pedestal    = 100;

    waveformGenerator gen;
    gen.setDigitizerParams(dq);
    gen.setDetectorParams(wp);

    auto raw = gen.generateWaveform(PulseShape::SPE);
    std::vector<Double_t> samples(raw.begin(), raw.end());

    TemplateAnalysisConfig cfg;
    cfg.fMFStart = 5;
    cfg.fMFEnd   = 105;

    FADCPulseWF result = AnalyzeWithTemplate(samples, tmpl, cfg);

    std::cout << "Pedestal  : " << result.pedestal << "\n";
    std::cout << "nPulses   : " << result.nPulses  << "\n";
    for (int i = 0; i < (int)result.pulses.size(); ++i) {
        const auto& p = result.pulses[i];
        std::cout << "  Pulse " << i
                  << "  t0="        << p.timeSample
                  << "  A="         << p.amplitude
                  << "  b="         << p.baseline
                  << "  chi2/ndf="  << p.chi2ndf
                  << "  MF="        << p.mfAmplitude
                  << "  valid="     << p.valid << "\n";
    }

    // Run assertions before opening the canvas – app.Run() blocks until
    // the window is closed, so any code after it would never execute.
    assert(result.nPulses >= 1 && "Expected at least one fitted pulse");
    std::cout << "Smoke test passed.\n";

    // Save the fit plot to a PNG for inspection (non-blocking).
    plotWaveformFit(samples, tmpl, result);
    return 0;
}
#endif // RUN_PROTO_TEST
