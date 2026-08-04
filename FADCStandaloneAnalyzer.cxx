#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"

#include <algorithm>    // std::min, std::max, std::sort, std::min_element, std::max_element
#include <cmath>        // std::sqrt, std::abs, std::exp
#include <limits>       // std::numeric_limits
#include <stdexcept>    // std::runtime_error

FADCStandaloneAnalyzer::FADCStandaloneAnalyzer() 
{
    fConfig = {
      /*safety config parameters*/
      4, 50, 2, 5, 2, 3
    };

    fCalib = {
      /*safety calibration parameters*/
      0.2441, 4.0/50.0
    };
}

FADCPulse FADCStandaloneAnalyzer::Analyze(const std::vector<Double_t>& samples, 
                                          Double_t gain, Double_t tcal) const 
{
  const Int_t fNSamples = static_cast<Int_t>(samples.size());
  FADCPulse data;

  // Guard: nothing to do
  if (fNSamples <= 0 || fConfig.fNPedestalSamples <= 0) return data;

  data.fSampPulseInt.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseAmp.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseTime.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseTOT.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulsePileup.assign(fConfig.fMaxNPulses, 0);
  data.fSampPulseIntPedSub.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseIntMOLLERRaw.assign(fConfig.fMaxNPulses, 0.0);
  data.fSampPulseIntMOLLERVal.assign(fConfig.fMaxNPulses, 0.0);

  auto GetIntegral = [&](Int_t lo, Int_t hi) {
    Double_t sum = 0.0;
    lo = TMath::Max(lo, 0);
    hi = TMath::Min(hi, fNSamples - 1);
    for (Int_t i = lo; i <= hi; ++i) sum += samples[i];
    return sum;
  };

  // GetSample returns the pedestal-subtracted value of sample i
  auto GetSample = [&](Int_t i) {
    return samples[i] - data.fSampPed / fConfig.fNPedestalSamples;
  };

  data.fHasMulti = kTRUE;
  data.fSampPed = GetIntegral(0, fConfig.fNPedestalSamples - 1);

  Int_t NS = fConfig.fNPedestalSamples - 1;
  data.fNSampPulses = 0;

  if (fConfig.fSampThreshold == 0) {
    data.fSampPulseInt[data.fNSampPulses] = GetIntegral(TMath::Max(NS - fConfig.fNSB, 0),
        TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1));
    data.fNPeakSamples = TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1)
      - TMath::Max(NS - fConfig.fNSB, 0) + 1;
    data.fPeakPedestalRatio = 1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;
    data.fSampPulseIntPedSub[data.fNSampPulses] =
      data.fSampPulseInt[data.fNSampPulses] - data.fSampPed * data.fPeakPedestalRatio;
    data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
      data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV * fCalib.pC_Conv;
    data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
      data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV * fCalib.pC_Conv * gain;
    data.fSampPulseAmp[data.fNSampPulses]  = samples[NS];   // raw (no pedestal division)
    data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
    data.fNSampPulses = 1;
  } else {
    Bool_t    CheckSampBelowThres = kFALSE;
    Double_t  LastTailPed = GetIntegral(fNSamples - fConfig.fNPedestalSamples, fNSamples - 1);
    if (LastTailPed < data.fSampPed) {
      data.fSampPed = LastTailPed;
      CheckSampBelowThres = kTRUE;
    }

    while (NS < fNSamples && data.fNSampPulses < fConfig.fMaxNPulses) {
      if (CheckSampBelowThres) {
        if (GetSample(NS) < fConfig.fSampThreshold) CheckSampBelowThres = kFALSE;
      } else {
        Int_t ns_found = 0;
        for (Int_t nt = NS; nt < TMath::Min(NS + fConfig.fNSAT, fNSamples); ++nt) {
          if (GetSample(nt) > fConfig.fSampThreshold) ns_found++;
        }
        if (ns_found == fConfig.fNSAT) {
          data.fSampPulseInt[data.fNSampPulses] =
            GetIntegral(TMath::Max(NS - fConfig.fNSB, 0), TMath::Min(NS + fConfig.fNSA - 1, fNSamples - 1));
          data.fNPeakSamples = TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1)
            - TMath::Max(NS - fConfig.fNSB, 0) + 1;
          data.fPeakPedestalRatio = 1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;

          data.fSampPulseIntPedSub[data.fNSampPulses] = 
            data.fSampPulseInt[data.fNSampPulses] - data.fSampPed * data.fPeakPedestalRatio;

          data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
            data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV * fCalib.pC_Conv;

          data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
            data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV * fCalib.pC_Conv * gain;
          
          // Peak finding: scan NSA window for the highest sample
          // Use -1 as sentinel (0 is a valid bin index)
          Int_t    PeakBin = -1;
          Double_t PeakVal = GetSample(NS);
          for (Int_t nt = NS + 1; nt < TMath::Min((NS + fConfig.fNSA), fNSamples); ++nt) {
            if (GetSample(nt) < PeakVal && PeakBin == -1) {
              PeakBin = nt - 1;
            } else if (PeakBin == -1) {
              PeakVal = GetSample(nt);
            }
          }
          if (PeakBin > -1) {
            data.fSampPulseAmp[data.fNSampPulses] = samples[PeakBin]; // raw sample
            Int_t Time = NS * 64;
            Double_t VMid = (GetSample(PeakBin)) / 2.;
            for( Int_t nt = TMath::Max(NS - fConfig.fNSB, 0);
                nt < TMath::Min(PeakBin, int(fNSamples - 1)); nt++ ) {
              if( VMid >= GetSample(nt) && VMid < GetSample(nt + 1) ) {
                Time = 64 * nt + int(64 * (VMid - GetSample(nt)) /
                    (GetSample(nt + 1) - GetSample(nt)));
              }
            }
            data.fSampPulseTime[data.fNSampPulses] = Time;
          } else {
            data.fSampPulseAmp[data.fNSampPulses]  = samples[NS]; // raw sample
            data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
          }

          // Pileup calculations
          Int_t hi = TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1);
          Int_t tot = 0;
          for( Int_t nt = NS; nt <= hi; nt++ ) {
            if( GetSample(nt) > fConfig.fSampThreshold ) tot++;
          }
          data.fSampPulseTOT[data.fNSampPulses] = 64.0 * tot;

          Int_t pileup = 0;
          for( Int_t nt = NS + 1; nt < hi; nt++ ) {
            if( GetSample(nt) > fConfig.fSampThreshold &&
                GetSample(nt) > GetSample(nt - 1) &&
                GetSample(nt) > GetSample(nt + 1) ) {
              pileup = 1;
              break;
            }
          }
          data.fSampPulsePileup[data.fNSampPulses] = pileup;

          data.fNSampPulses++;
          NS = NS + fConfig.fNSA;
          CheckSampBelowThres = kTRUE;
        }
      }
      NS++;
    }
  }
  return data;
}

// Fit based functions
// FADC Template
FADCTemplate FADCTemplate::fromWaveformParams(const WaveformParams& wp,
                                               Double_t sampleRate,
                                               Int_t    nSamples,
                                               Double_t peakTime_ns,
                                               Int_t    mfLeft,
                                               Int_t    mfRight) {
    FADCTemplate tmpl;
    tmpl.shape.resize(nSamples, 0.0);

    // Hardcoded! Need to figureout the proper way.
    Double_t tau_r = 2.0;
    Double_t tau_f = 5.0;

    // Evaluate the raw SPE shape using WaveformUtil::spePulse in waveformGenerator
    for (Int_t i = 0; i < nSamples; ++i)
        tmpl.shape[i] = WaveformUtil::spePulse(i * sampleRate, peakTime_ns, tau_r, tau_f);

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

// Compute Matched Filter
std::vector<Double_t> FADCStandaloneAnalyzer::computeMatchedFilter(
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

// File-local struct shared between findMFPeaks and AnalyzeWithTemplate
struct MFPeak {
    Double_t position;    // sub-sample peak position (samples)
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

    // Per-sample sigma — sized to fit window only, accessed as sigma[i - fitLo]
    std::vector<Double_t> sigma(nPts, cfg.fMinSigma);
    for (Int_t i = fitLo; i <= fitHi; ++i) {
        Double_t s = std::sqrt(std::abs(wfResidual[i]) * cfg.fAdcPerMV / 2.0)
                     / cfg.fAdcPerMV;
        sigma[i - fitLo] = (s > cfg.fMinSigma) ? s : cfg.fMinSigma;
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

        // Accumulate normal-equation sums and chi2 in a single pass
        // (avoids calling tmpl.eval() twice per sample per t0 step)
        Double_t S1  = 0.0, Sf  = 0.0, Sf2 = 0.0;
        Double_t Sy  = 0.0, Syf = 0.0;

        for (Int_t i = fitLo; i <= fitHi; ++i) {
            Double_t f     = tmpl.eval(static_cast<Double_t>(i) - t0);
            Double_t invs2 = 1.0 / (sigma[i - fitLo] * sigma[i - fitLo]);
            S1  += invs2;
            Sf  += f * invs2;
            Sf2 += f * f * invs2;
            Sy  += wfResidual[i] * invs2;
            Syf += wfResidual[i] * f * invs2;
        }

        // Solve the 2x2 system for (baseline b, amplitude A)
        Double_t det = S1 * Sf2 - Sf * Sf;
        if (std::abs(det) < 1e-12) continue;

        Double_t b = (Sy * Sf2 - Syf * Sf) / det;   // baseline
        Double_t A = (S1 * Syf - Sf  * Sy) / det;   // amplitude

        if (A <= 0.0) continue;  // unphysical -- skip

        // chi2 — second pass eliminated: derive directly from the normal-equation sums
        //   chi2 = sum[ (y - b - A*f)^2 / sigma^2 ]
        //        = Sy_y2 - 2b*Sy - 2A*Syf + b^2*S1 + 2bA*Sf + A^2*Sf2
        // We still need sum(y^2/sigma^2), compute it here
        Double_t Sy2 = 0.0;
        for (Int_t i = fitLo; i <= fitHi; ++i) {
            Double_t invs2 = 1.0 / (sigma[i - fitLo] * sigma[i - fitLo]);
            Sy2 += wfResidual[i] * wfResidual[i] * invs2;
        }
        Double_t chi2    = Sy2 - 2.0*b*Sy - 2.0*A*Syf + b*b*S1 + 2.0*b*A*Sf + A*A*Sf2;
        Int_t    ndf     = std::max(nPts - 2, 1);   // 2 free parameters
        Double_t chi2ndf = chi2 / ndf;

        if (chi2ndf < best.chi2ndf) {
            best.chi2ndf    = chi2ndf;
            best.timeSample = t0;
            best.amplitude  = A;
            best.baseline   = b;
            best.valid      = kTRUE;
        }
    }

    return best;
}

// Top-level analysis
FADCPulseWF FADCStandaloneAnalyzer::AnalyzeWithTemplate(
        const std::vector<Double_t>& samples,
        const FADCTemplate&          tmpl,
        const TemplateAnalysisConfig& cfg) const
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
    for (Int_t i = N - nPed; i < N; ++i) ped2 += samples[i];
    ped2 /= nPed;

    // Use the lower of front/tail pedestals (same logic as Analyze())
    result.pedestal = (ped2 < ped) ? ped2 : ped;

    // Subtract the chosen pedestal consistently
    std::vector<Double_t> wfPS(N);
    for (Int_t i = 0; i < N; ++i) wfPS[i] = samples[i] - result.pedestal;

    // 2. Matched filter
    std::vector<Double_t> mfOut = computeMatchedFilter(wfPS, tmpl, cfg);

    // 3. Peak finding
    std::vector<MFPeak> mfPeaks = findMFPeaks(mfOut, cfg);

    // Sort by MF amplitude (strongest first) so dominant pulses are fitted
    // before noise peaks can corrupt the residual.
    std::sort(mfPeaks.begin(), mfPeaks.end(),
              [](const MFPeak& a, const MFPeak& b) {
                  return a.mfAmplitude > b.mfAmplitude;
              });

    // 4. Template fit with sequential subtraction (strongest-first peeling)
    std::vector<Double_t> wfResidual = wfPS;

    // Calculate the area (sum) of the normalized template shape
    Double_t templateArea = 0.0;
    for (Double_t val : tmpl.shape) {
        templateArea += val;
    }

    // Track positions of successfully fitted pulses to skip nearby noise peaks
    std::vector<Double_t> fittedPositions;

    for (const auto& peak : mfPeaks) {
        // Stop if we've fitted enough pulses
        if (result.nPulses >= cfg.fMaxPulses) break;

        // Skip if this peak is too close to an already-fitted pulse
        Bool_t tooClose = kFALSE;
        for (Double_t pos : fittedPositions) {
            if (std::abs(peak.position - pos) < static_cast<Double_t>(cfg.fMFSeparation)) {
                tooClose = kTRUE;
                break;
            }
        }
        if (tooClose) continue;

        // Local fit window around the candidate peak — prevents the baseline
        // from absorbing the signal when most of the waveform is flat noise.
        Int_t fitLo = std::max(0, static_cast<Int_t>(peak.position) - cfg.fFitMargin);
        Int_t fitHi = std::min(N - 1, static_cast<Int_t>(peak.position) + cfg.fFitMargin);

        FitPulseResult fit = fitSingleTemplate(
                wfResidual, tmpl,
                peak.position,
                fitLo, fitHi,
                cfg);

        fit.mfAmplitude = peak.mfAmplitude;

        // Calculate the pulse integral: Integral = amplitude * area(template)
        fit.integral = fit.amplitude * templateArea;

        if (!fit.valid || fit.amplitude < cfg.fAmpCutADC) {
            continue;   // skip noise — don't store or subtract
        }

        // Subtract ONLY the fitted pulse shape from the residual.
        // The baseline b is a local nuisance parameter (DC offset within the fit
        // window) — subtracting it globally would shift the entire waveform and
        // corrupt later fits.
        for (Int_t i = 0; i < N; ++i) {
            Double_t f = tmpl.eval(static_cast<Double_t>(i) - fit.timeSample);
            wfResidual[i] -= fit.amplitude * f;
        }

        fittedPositions.push_back(fit.timeSample);
        result.pulses.push_back(fit);
        result.nPulses++;
    }

    // Sort results chronologically for output
    std::sort(result.pulses.begin(), result.pulses.end(),
              [](const FitPulseResult& a, const FitPulseResult& b) {
                  return a.timeSample < b.timeSample; });

    return result;
}