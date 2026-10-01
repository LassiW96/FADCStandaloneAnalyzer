#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

FADCStandaloneAnalyzer::FADCStandaloneAnalyzer() {
  fConfig = {/*safety config parameters*/
             4, 50, 2, 5, 2, 3};

  fCalib = {/*safety calibration parameters*/
            0.2441, 4.0 / 50.0};
}

FADCPulse FADCStandaloneAnalyzer::Analyze(const std::vector<Double_t> &samples,
                                          Double_t gain, Double_t tcal) const {
  const Int_t fNSamples = static_cast<Int_t>(samples.size());
  FADCPulse data;

  // Guard: nothing to do
  if (fNSamples <= 0 || fConfig.fNPedestalSamples <= 0)
    return data;

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
    for (Int_t i = lo; i <= hi; ++i)
      sum += samples[i];
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
    data.fSampPulseInt[data.fNSampPulses] =
        GetIntegral(TMath::Max(NS - fConfig.fNSB, 0),
                    TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1));
    data.fNPeakSamples = TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1) -
                         TMath::Max(NS - fConfig.fNSB, 0) + 1;
    data.fPeakPedestalRatio =
        1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;
    data.fSampPulseIntPedSub[data.fNSampPulses] =
        data.fSampPulseInt[data.fNSampPulses] -
        data.fSampPed * data.fPeakPedestalRatio;
    data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
        data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV *
        fCalib.pC_Conv;
    data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
        data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV *
        fCalib.pC_Conv * gain;
    data.fSampPulseAmp[data.fNSampPulses] =
        samples[NS]; // raw (no pedestal division)
    data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
    data.fNSampPulses = 1;
  } else {
    Bool_t CheckSampBelowThres = kFALSE;
    Double_t LastTailPed =
        GetIntegral(fNSamples - fConfig.fNPedestalSamples, fNSamples - 1);
    if (LastTailPed < data.fSampPed) {
      data.fSampPed = LastTailPed;
      CheckSampBelowThres = kTRUE;
    }

    while (NS < fNSamples && data.fNSampPulses < fConfig.fMaxNPulses) {
      if (CheckSampBelowThres) {
        if (GetSample(NS) < fConfig.fSampThreshold)
          CheckSampBelowThres = kFALSE;
      } else {
        Int_t ns_found = 0;
        for (Int_t nt = NS; nt < TMath::Min(NS + fConfig.fNSAT, fNSamples);
             ++nt) {
          if (GetSample(nt) > fConfig.fSampThreshold)
            ns_found++;
        }
        if (ns_found == fConfig.fNSAT) {
          data.fSampPulseInt[data.fNSampPulses] =
              GetIntegral(TMath::Max(NS - fConfig.fNSB, 0),
                          TMath::Min(NS + fConfig.fNSA - 1, fNSamples - 1));
          data.fNPeakSamples =
              TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1) -
              TMath::Max(NS - fConfig.fNSB, 0) + 1;
          data.fPeakPedestalRatio =
              1.0 * data.fNPeakSamples / fConfig.fNPedestalSamples;

          data.fSampPulseIntPedSub[data.fNSampPulses] =
              data.fSampPulseInt[data.fNSampPulses] -
              data.fSampPed * data.fPeakPedestalRatio;

          data.fSampPulseIntMOLLERRaw[data.fNSampPulses] =
              data.fSampPulseInt[data.fNSampPulses] * fCalib.ChanTomV *
              fCalib.pC_Conv;

          data.fSampPulseIntMOLLERVal[data.fNSampPulses] =
              data.fSampPulseIntPedSub[data.fNSampPulses] * fCalib.ChanTomV *
              fCalib.pC_Conv * gain;

          // Peak finding: scan NSA window for the highest sample
          // Use -1 as sentinel (0 is a valid bin index)
          Int_t PeakBin = -1;
          Double_t PeakVal = GetSample(NS);
          for (Int_t nt = NS + 1;
               nt < TMath::Min((NS + fConfig.fNSA), fNSamples); ++nt) {
            if (GetSample(nt) < PeakVal && PeakBin == -1) {
              PeakBin = nt - 1;
            } else if (PeakBin == -1) {
              PeakVal = GetSample(nt);
            }
          }
          if (PeakBin > -1) {
            data.fSampPulseAmp[data.fNSampPulses] =
                samples[PeakBin]; // raw sample
            Int_t Time = NS * 64;
            Double_t VMid = (GetSample(PeakBin)) / 2.;
            for (Int_t nt = TMath::Max(NS - fConfig.fNSB, 0);
                 nt < TMath::Min(PeakBin, int(fNSamples - 1)); nt++) {
              if (VMid >= GetSample(nt) && VMid < GetSample(nt + 1)) {
                Time = 64 * nt + int(64 * (VMid - GetSample(nt)) /
                                     (GetSample(nt + 1) - GetSample(nt)));
              }
            }
            data.fSampPulseTime[data.fNSampPulses] = Time;
          } else {
            data.fSampPulseAmp[data.fNSampPulses] = samples[NS]; // raw sample
            data.fSampPulseTime[data.fNSampPulses] = 64 * NS;
          }

          // Pileup calculations
          Int_t hi = TMath::Min((NS + fConfig.fNSA - 1), fNSamples - 1);
          Int_t tot = 0;
          for (Int_t nt = NS; nt <= hi; nt++) {
            if (GetSample(nt) > fConfig.fSampThreshold)
              tot++;
          }
          data.fSampPulseTOT[data.fNSampPulses] = 64.0 * tot;

          Int_t pileup = 0;
          for (Int_t nt = NS + 1; nt < hi; nt++) {
            if (GetSample(nt) > fConfig.fSampThreshold &&
                GetSample(nt) > GetSample(nt - 1) &&
                GetSample(nt) > GetSample(nt + 1)) {
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

Double_t FADCTemplate::eval(Double_t t) const {
  if (t < 0.0 || t >= static_cast<Double_t>(shape.size()))
    return 0.0;
  Int_t lo = static_cast<Int_t>(t);
  Double_t frac = t - lo;
  if (lo >= static_cast<Int_t>(shape.size()) - 1)
    return shape.back();
  return shape[lo] * (1.0 - frac) + shape[lo + 1] * frac;
}

// Builds a template from a measured reference pulse (e.g., from a laser or
// LED calibration run). This path is model-independent: the measured shape is
// used directly, so no assumptions about sigma / tau are made.
//
// Usage pattern:
//   FADCTemplate tmpl = FADCTemplate::fromRawSamples(myCalibPulse, kTRUE);
FADCTemplate FADCTemplate::fromRawSamples(const std::vector<Double_t> &rawShape,
                                          Bool_t enable, Int_t mfLeft,
                                          Int_t mfRight) {
  FADCTemplate
      tmpl; // default-constructed: shape is empty, peakBin = -1 (undefined)
  if (!enable || rawShape.size() < 2)
    return tmpl;

  Int_t nSamples = static_cast<Int_t>(rawShape.size());
  tmpl.shape = rawShape;

  // Normalise to peak = 1 and record the peak bin
  Double_t peak = *std::max_element(tmpl.shape.begin(), tmpl.shape.end());
  if (peak <= 0.0)
    throw std::runtime_error(
        "FADCTemplate::fromRawSamples: zero-amplitude shape");

  tmpl.peakBin = -1;
  for (Int_t i = 0; i < nSamples; ++i) {
    tmpl.shape[i] /= peak;
    if (tmpl.shape[i] >= 1.0 - 1e-9)
      tmpl.peakBin = i;
  }
  if (tmpl.peakBin < 0) {
    // Fallback: take the bin with the highest value
    tmpl.peakBin = static_cast<Int_t>(
        std::max_element(tmpl.shape.begin(), tmpl.shape.end()) -
        tmpl.shape.begin());
  }

  // Build the matched-filter kernel
  tmpl.mfKernel.clear();
  tmpl.mfNorm = 0.0;
  for (Int_t i = tmpl.peakBin + mfRight; i >= tmpl.peakBin - mfLeft; --i) {
    Double_t v = (i >= 0 && i < nSamples) ? tmpl.shape[i] : 0.0;
    tmpl.mfKernel.push_back(v);
    tmpl.mfNorm += v;
  }
  if (tmpl.mfNorm <= 0.0)
    throw std::runtime_error("FADCTemplate::fromRawSamples: zero MF norm");

  return tmpl;
}

// ---------------------------------------------------------------------------
// fromWaveformSamples — build template from a raw waveform sample vector.
//
// Subtracts pedestal (mean of the leading nPedSamples) and delegates to
// fromRawSamples for peak normalisation and MF kernel construction.
// This is the primary factory for template-based analysis: pass in a clean
// single-pulse waveform (e.g. averaged calibration pulse).
// ---------------------------------------------------------------------------
FADCTemplate FADCTemplate::fromWaveformSamples(
    const std::vector<Double_t> &samples,
    Int_t nPedSamples, Int_t mfLeft, Int_t mfRight) {

  Int_t N = static_cast<Int_t>(samples.size());
  if (N < 2) return FADCTemplate{};  // empty template

  // Pedestal subtraction
  Int_t nPed = std::min(nPedSamples, N);
  Double_t ped = 0.0;
  for (Int_t i = 0; i < nPed; ++i) ped += samples[i];
  ped /= nPed;

  std::vector<Double_t> pedSubtracted(N);
  for (Int_t i = 0; i < N; ++i) pedSubtracted[i] = samples[i] - ped;

  // Delegate to fromRawSamples for normalisation and MF kernel
  return fromRawSamples(pedSubtracted, kTRUE, mfLeft, mfRight);
}

// Int_t convenience overload (e.g. from waveformGenerator or digitizer output)
FADCTemplate FADCTemplate::fromWaveformSamples(
    const std::vector<Int_t> &samples,
    Int_t nPedSamples, Int_t mfLeft, Int_t mfRight) {

  std::vector<Double_t> dSamples(samples.begin(), samples.end());
  return fromWaveformSamples(dSamples, nPedSamples, mfLeft, mfRight);
}

// Compute Matched Filter
std::vector<Double_t> FADCStandaloneAnalyzer::computeMatchedFilter(
    const std::vector<Double_t> &wfPS, const FADCTemplate &tmpl,
    const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(wfPS.size());
  Int_t kW = static_cast<Int_t>(tmpl.mfKernel.size()); // mfLeft + mfRight + 1

  Double_t wfMin = *std::min_element(wfPS.begin(), wfPS.end());

  std::vector<Double_t> mfOut(N, 0.0);
  Double_t mfOutMin = std::numeric_limits<Double_t>::max();

  for (Int_t i = 0; i < N; ++i) {
    Double_t val = 0.0;
    for (Int_t j = 0; j < kW; ++j) {
      Int_t wIdx = i + j - cfg.fMFRight;
      if (wIdx >= 0 && wIdx < N)
        val += (wfPS[wIdx] - wfMin) * tmpl.mfKernel[j];
    }
    mfOut[i] = val / tmpl.mfNorm;
    if (mfOut[i] < mfOutMin)
      mfOutMin = mfOut[i];
  }

  // Shift so the global minimum of the MF output = 0  (mirrors fit.C line 905)
  for (Int_t i = 0; i < N; ++i)
    mfOut[i] -= mfOutMin;

  return mfOut;
}

// ---------------------------------------------------------------------------
// File-local struct shared between findMFPeaks and AnalyzeWithTemplate
// ---------------------------------------------------------------------------
struct MFPeak {
  Double_t position;    // sub-sample peak position (samples)
  Double_t mfAmplitude; // MF output value at the peak
};

// ---------------------------------------------------------------------------
// findMFPeaks — locate peak candidates in the matched-filter trace
// ---------------------------------------------------------------------------
static std::vector<MFPeak> findMFPeaks(
    const std::vector<Double_t> &mfOut,
    const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(mfOut.size());
  Int_t mfEnd = std::min(cfg.fMFEnd, N - 3);

  std::vector<MFPeak> peaks;

  auto isDuplicate = [&](Double_t pos) -> Bool_t {
    for (const auto &p : peaks) {
      if (std::abs(p.position - pos) < static_cast<Double_t>(cfg.fMFSeparation))
        return kTRUE;
    }
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
      Double_t pos = (std::abs(denom) > 1e-9)
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
            [](const MFPeak &a, const MFPeak &b) {
              return a.position < b.position;
            });

  return peaks;
}


// ---------------------------------------------------------------------------
// evalRawEMG — file-local helper used by EMGFitFunctor and AnalyzeWithTemplate.
// Evaluates the raw (unnormalized) EMG pulse shape at time offset dt,
// with Gaussian sigma and exponential tau, both in sample units.
// ---------------------------------------------------------------------------
static Double_t evalRawEMG(Double_t dt, Double_t sigma, Double_t tau) {
  if (sigma <= 0.0 || tau <= 0.0)
    return 0.0;
  Double_t s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
  Double_t erf_arg = (dt - (sigma * sigma) / tau) / (sigma * std::sqrt(2.0));

  if (erf_arg > 700.0)
    return 0.0; // overflow guard

  Double_t val = std::exp(s2_2t2 - (dt / tau)) * (1.0 + std::erf(erf_arg));
  return (val > 1e-12) ? val : 0.0;
}

// ---------------------------------------------------------------------------
// EMGFitFunctor — chi2 callable for Stage-2 Minuit fit.
// Parameters: [0]=A (amplitude), [1]=b (baseline), [2]=t0,
//             [3]=sigma (sample units), [4]=tau (sample units).
// ---------------------------------------------------------------------------
struct EMGFitFunctor {
  const std::vector<Double_t> *wfResidual;
  Int_t fitLo;
  Int_t fitHi;
  Double_t minSigma;
  Double_t adcPerMV;
  Double_t noiseSigma;

  Double_t operator()(const Double_t *par) const {
    Double_t A = par[0];
    Double_t b = par[1];
    Double_t t0 = par[2];
    Double_t sigma = par[3];
    Double_t tau = par[4];

    Double_t chi2 = 0.0;
    for (Int_t i = fitLo; i <= fitHi; ++i) {
      Double_t dt = static_cast<Double_t>(i) - t0;
      Double_t f = evalRawEMG(dt, sigma, tau);
      Double_t expected = b + A * f;
      Double_t s_shot =
          std::sqrt(std::abs((*wfResidual)[i]) * adcPerMV / 2.0) / adcPerMV;
      Double_t sig = std::sqrt(s_shot * s_shot + noiseSigma * noiseSigma);
      if (sig < minSigma)
        sig = minSigma;
      Double_t residual = ((*wfResidual)[i] - expected) / sig;
      chi2 += residual * residual;
    }
    return chi2;
  }
};

// ---------------------------------------------------------------------------
// fitStage1Template — Stage-1 fast fit.
//
// Performs an analytical 2-parameter (A, b) scan over the fixed template shape.
// For each t0 in [candidateBin - fScanHalf, candidateBin + fScanHalf] at
// fScanStep resolution the weighted normal equations are solved exactly:
//   [ S1  Sf  ] [b]   [Sy ]
//   [ Sf  Sf2 ] [A] = [Syf]
// giving chi2 = sum((y - b - A*f)^2 / sigma_i^2) with NDF = nPts - 2.
// ---------------------------------------------------------------------------
static FitPulseResult fitStage1Template(
    const std::vector<Double_t> &wfResidual,
    const FADCTemplate &tmpl,
    Double_t candidateBin,
    Int_t fitLo,
    Int_t fitHi,
    const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(wfResidual.size());
  fitLo = std::max(fitLo, 0);
  fitHi = std::min(fitHi, N - 1);
  Int_t nPts = fitHi - fitLo + 1;

  // Precompute invSigma2 and Sy2 outside the t0 scan loop
  std::vector<Double_t> invSigma2(nPts);
  Double_t Sy2 = 0.0;
  for (Int_t i = fitLo; i <= fitHi; ++i) {
    Double_t s_shot = std::sqrt(std::abs(wfResidual[i]) * cfg.fAdcPerMV / 2.0)
                      / cfg.fAdcPerMV;
    Double_t s = std::sqrt(s_shot * s_shot + cfg.fNoiseSigmaADC * cfg.fNoiseSigmaADC);
    Double_t sig = (s > cfg.fMinSigma) ? s : cfg.fMinSigma;
    Double_t invs2 = 1.0 / (sig * sig);
    invSigma2[i - fitLo] = invs2;
    Sy2 += wfResidual[i] * wfResidual[i] * invs2;
  }

  FitPulseResult best;
  best.valid         = kFALSE;
  best.chi2ndf       = std::numeric_limits<Double_t>::max();
  best.chi2ndf_s1    = std::numeric_limits<Double_t>::max();
  best.timeSample    = candidateBin;
  best.timeHalfMax   = candidateBin;
  best.amplitude     = 0.0;
  best.integral      = 0.0;
  best.baseline      = 0.0;
  best.mfAmplitude   = 0.0;
  best.fitStage      = 1;
  best.fittedSigma   = 0.0;
  best.fittedTau     = 0.0;

  // Scan t0 window
  for (Double_t t0 = candidateBin - cfg.fScanHalf;
       t0 <= candidateBin + cfg.fScanHalf;
       t0 += cfg.fScanStep) {

    // Accumulate normal-equation sums in a single pass
    Double_t S1  = 0.0, Sf  = 0.0, Sf2 = 0.0;
    Double_t Sy  = 0.0, Syf = 0.0;

    for (Int_t i = fitLo; i <= fitHi; ++i) {
      Double_t f     = tmpl.eval(static_cast<Double_t>(i) - t0);
      Double_t invs2 = invSigma2[i - fitLo];
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

    if (A <= 0.0) continue;  // unphysical — skip

    // chi2 derived from normal-equation sums
    Double_t chi2    = Sy2 - 2.0*b*Sy - 2.0*A*Syf + b*b*S1 + 2.0*b*A*Sf + A*A*Sf2;
    Int_t    ndf     = std::max(nPts - 2, 1);   // 2 free parameters
    Double_t chi2ndf = chi2 / ndf;

    if (chi2ndf < best.chi2ndf) {
      best.chi2ndf    = chi2ndf;
      best.chi2ndf_s1 = chi2ndf;
      best.timeSample = t0;
      best.amplitude  = A;
      best.baseline   = b;
      best.valid      = kTRUE;
    }
  }
  return best;
}

// ---------------------------------------------------------------------------
// fitStage2Template — Stage-2 full Minuit2 fit (5 free parameters).
//
// wfPeakBin is the estimated waveform-peak sample (= Stage-1 t0 + peakBin,
// or original MF peak position).
// ---------------------------------------------------------------------------
static FitPulseResult fitStage2Template(
    const std::vector<Double_t> &wfResidual,
    const FADCTemplate &tmpl,
    Double_t wfPeakBin,
    Int_t fitLo,
    Int_t fitHi,
    const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(wfResidual.size());
  fitLo = std::max(fitLo, 0);
  fitHi = std::min(fitHi, N - 1);
  Int_t nPts = fitHi - fitLo + 1;

  FitPulseResult best;
  best.valid         = kFALSE;
  best.chi2ndf       = std::numeric_limits<Double_t>::max();
  best.chi2ndf_s1    = std::numeric_limits<Double_t>::max();
  best.timeSample    = wfPeakBin - static_cast<Double_t>(tmpl.peakBin);  // fallback
  best.timeHalfMax   = best.timeSample;
  best.amplitude     = 0.0;
  best.integral      = 0.0;
  best.baseline      = 0.0;
  best.mfAmplitude   = 0.0;
  best.fitStage      = 2;
  best.fittedSigma   = 0.0;
  best.fittedTau     = 0.0;

  std::unique_ptr<ROOT::Math::Minimizer> min(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!min) return best;

  min->SetMaxFunctionCalls(1000);
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);

  EMGFitFunctor chi2Functor{&wfResidual, fitLo, fitHi, cfg.fMinSigma, cfg.fAdcPerMV, cfg.fNoiseSigmaADC};
  ROOT::Math::Functor functor(chi2Functor, 5);
  min->SetFunction(functor);

  // Rough amplitude seed from the residual window
  Double_t maxVal = 0.0;
  for (Int_t i = fitLo; i <= fitHi; ++i) {
    if (wfResidual[i] > maxVal) maxVal = wfResidual[i];
  }

  // Nominal sigma and tau in sample units (sigma=5 ns, tau=10 ns at 4 ns/sample)
  const Double_t initSigma = 1.25;
  const Double_t initTau   = 2.5;

  // Seed t0 so initial EMG peaks near the known waveform peak.
  Double_t initT0 = wfPeakBin - (initSigma * initSigma / initTau);

  min->SetLimitedVariable(0, "Amplitude", maxVal,    maxVal * 0.1, 0.0,                    maxVal * 3.0);
  min->SetLimitedVariable(1, "Baseline",  0.0,       1.0,          -20.0,                  20.0);
  min->SetLimitedVariable(2, "t0",        initT0,    0.5,          initT0 - cfg.fScanHalf, initT0 + cfg.fScanHalf);
  min->SetLimitedVariable(3, "Sigma",     initSigma, 0.1,          0.1,                    5.0);
  min->SetLimitedVariable(4, "Tau",       initTau,   0.2,          0.5,                    15.0);

  min->Minimize();

  if (min->Status() == 0 || min->Status() == 1) {
    const Double_t *xs = min->X();
    Double_t A_raw  = xs[0];
    Double_t t0     = xs[2];
    Double_t sigma  = xs[3];
    Double_t tau    = xs[4];

    // Find the peak of the raw (unnormalized) EMG over the fit window to convert
    // to physical peak ADC count.
    Double_t emgPeakRaw = 0.0;
    for (Int_t i = fitLo; i <= fitHi; ++i) {
      Double_t f = evalRawEMG(static_cast<Double_t>(i) - t0, sigma, tau);
      if (f > emgPeakRaw) emgPeakRaw = f;
    }
    if (emgPeakRaw < 1e-12) emgPeakRaw = 1.0;

    best.amplitude   = A_raw * emgPeakRaw;
    best.baseline    = xs[1];
    best.timeSample  = t0;
    best.fittedSigma = sigma;
    best.fittedTau   = tau;

    Double_t chi2 = min->MinValue();
    Int_t ndf = std::max(nPts - 5, 1);
    best.chi2ndf = chi2 / ndf;
    best.valid   = kTRUE;
  }

  return best;
}

// ---------------------------------------------------------------------------
// fitSingleTemplate — two-stage dispatcher.
// ---------------------------------------------------------------------------
static FitPulseResult fitSingleTemplate(
    const std::vector<Double_t> &wfResidual,
    const FADCTemplate &tmpl,
    Double_t candidateBin,
    Int_t fitLo,
    Int_t fitHi,
    const TemplateAnalysisConfig &cfg) {
  // Stage 1: fast 2-parameter fixed-template scan
  FitPulseResult s1 = fitStage1Template(wfResidual, tmpl, candidateBin, fitLo, fitHi, cfg);

  if (s1.valid && s1.chi2ndf <= cfg.fStage2Chi2Threshold)
    return s1;  // Stage 1 accepted

  // Stage 2: full Minuit2 fit fallback
  Double_t wfPeakBin = s1.valid
      ? (s1.timeSample + static_cast<Double_t>(tmpl.peakBin))
      : (candidateBin  + static_cast<Double_t>(tmpl.peakBin));

  FitPulseResult s2 = fitStage2Template(wfResidual, tmpl, wfPeakBin, fitLo, fitHi, cfg);
  s2.chi2ndf_s1 = s1.chi2ndf;

  if (!s2.valid && s1.valid)
    return s1;

  return s2;
}

// ---------------------------------------------------------------------------
// addPulseContribution — add or subtract a fitted pulse's model from a waveform.
// sign = +1.0 to add back, -1.0 to subtract.
// ---------------------------------------------------------------------------
static void addPulseContribution(
    std::vector<Double_t> &wf,
    const FitPulseResult &pulse,
    const FADCTemplate &tmpl,
    Double_t sign) {
  Int_t N = static_cast<Int_t>(wf.size());
  if (pulse.fitStage == 1 || pulse.fittedSigma == 0.0) {
    for (Int_t i = 0; i < N; ++i) {
      Double_t f = tmpl.eval(static_cast<Double_t>(i) - pulse.timeSample);
      wf[i] += sign * pulse.amplitude * f;
    }
  } else {
    // Stage 2: normalized raw EMG
    Double_t emgPeakRaw = 0.0;
    for (Int_t i = 0; i < N; ++i) {
      Double_t f = evalRawEMG(
          static_cast<Double_t>(i) - pulse.timeSample,
          pulse.fittedSigma, pulse.fittedTau);
      if (f > emgPeakRaw) emgPeakRaw = f;
    }
    if (emgPeakRaw < 1e-12) emgPeakRaw = 1.0;
    for (Int_t i = 0; i < N; ++i) {
      Double_t f = evalRawEMG(
          static_cast<Double_t>(i) - pulse.timeSample,
          pulse.fittedSigma, pulse.fittedTau) / emgPeakRaw;
      wf[i] += sign * pulse.amplitude * f;
    }
  }
}

// ---------------------------------------------------------------------------
// computePulseDerivedQuantities — compute integral, timeHalfMax, and peak
// position for a FitPulseResult.
// ---------------------------------------------------------------------------
static void computePulseDerivedQuantities(
    FitPulseResult &fit,
    const std::vector<Double_t> &wfForTiming,
    const FADCTemplate &tmpl,
    Double_t templateArea,
    const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(wfForTiming.size());

  // Integral
  if (fit.fitStage == 1 || fit.fittedSigma == 0.0) {
    fit.integral = fit.amplitude * templateArea;
  } else {
    Double_t emgAreaRaw = 0.0;
    Double_t emgPeakRaw = 0.0;
    Int_t nTmpl = static_cast<Int_t>(tmpl.shape.size());
    for (Int_t i = 0; i < nTmpl; ++i) {
      Double_t dt = static_cast<Double_t>(i) - static_cast<Double_t>(tmpl.peakBin);
      Double_t f  = evalRawEMG(dt, fit.fittedSigma, fit.fittedTau);
      emgAreaRaw += f;
      if (f > emgPeakRaw) emgPeakRaw = f;
    }
    if (emgPeakRaw < 1e-12) emgPeakRaw = 1.0;
    fit.integral = fit.amplitude * emgAreaRaw / emgPeakRaw;
  }

  // Half-max rising-edge crossing (CFD-style timing)
  Int_t wfPeak;
  if (fit.fitStage == 1 || fit.fittedSigma == 0.0) {
    wfPeak = static_cast<Int_t>(fit.timeSample) + tmpl.peakBin;
  } else {
    wfPeak = static_cast<Int_t>(std::round(
        fit.timeSample + fit.fittedSigma * fit.fittedSigma / fit.fittedTau));
  }
  Int_t    hLo     = std::max(0,   wfPeak - cfg.fFitMargin);
  Int_t    hHi     = std::min(N-1, wfPeak + cfg.fFitMargin);
  Double_t locPeak = 0.0;
  Int_t    locBin  = wfPeak;
  for (Int_t i = hLo; i <= hHi; ++i) {
    if (wfForTiming[i] > locPeak) { locPeak = wfForTiming[i]; locBin = i; }
  }
  Double_t vMid    = locPeak / 2.0;
  Double_t halfMax = static_cast<Double_t>(locBin);  // fallback: peak bin
  for (Int_t i = hLo; i < locBin; ++i) {
    if (wfForTiming[i] <= vMid && wfForTiming[i + 1] > vMid) {
      halfMax = i + (vMid - wfForTiming[i])
                    / (wfForTiming[i + 1] - wfForTiming[i]);
      break;
    }
  }
  fit.timeHalfMax = halfMax;
}

// ---------------------------------------------------------------------------
// getPulsePeakPosition — compute the waveform-peak sample position for a pulse.
// ---------------------------------------------------------------------------
static Double_t getPulsePeakPosition(const FitPulseResult &fit,
                                     const FADCTemplate &tmpl) {
  if (fit.fitStage == 1 || fit.fittedSigma == 0.0)
    return fit.timeSample + static_cast<Double_t>(tmpl.peakBin);
  return fit.timeSample + fit.fittedSigma * fit.fittedSigma / fit.fittedTau;
}

// ---------------------------------------------------------------------------
// AnalyzeWithTemplate — Top-level template-based analysis method
// ---------------------------------------------------------------------------
FADCPulseWF FADCStandaloneAnalyzer::AnalyzeWithTemplate(
    const std::vector<Double_t> &samples,
    const FADCTemplate &tmpl,
    const TemplateAnalysisConfig &cfg) const {
  FADCPulseWF result;
  result.nPulses = 0;

  Int_t N = static_cast<Int_t>(samples.size());
  if (N == 0) return result;

  // 1. Pedestal — noise-aware front vs tail selection
  Int_t nPed = std::min(cfg.fNPedSamples, N);

  // Compute front pedestal mean and sample variance
  Double_t ped = 0.0;
  for (Int_t i = 0; i < nPed; ++i) ped += samples[i];
  ped /= nPed;
  Double_t pedVar = 0.0;
  for (Int_t i = 0; i < nPed; ++i) {
    Double_t d = samples[i] - ped;
    pedVar += d * d;
  }
  pedVar /= (nPed > 1) ? (nPed - 1) : 1;

  // Compute tail pedestal mean and sample variance
  Double_t ped2 = 0.0;
  for (Int_t i = N - nPed; i < N; ++i) ped2 += samples[i];
  ped2 /= nPed;
  Double_t ped2Var = 0.0;
  for (Int_t i = N - nPed; i < N; ++i) {
    Double_t d = samples[i] - ped2;
    ped2Var += d * d;
  }
  ped2Var /= (nPed > 1) ? (nPed - 1) : 1;

  Double_t pedNoiseSigma = std::sqrt((ped <= ped2) ? pedVar : ped2Var);

  if (pedNoiseSigma > 1e-6) {
    Double_t diffThreshold = 3.0 * pedNoiseSigma * std::sqrt(2.0 / nPed);
    if (std::abs(ped - ped2) > diffThreshold) {
      result.pedestal = std::min(ped, ped2);
    } else {
      result.pedestal = (ped + ped2) / 2.0;
    }
  } else {
    result.pedestal = std::min(ped, ped2);
  }

  // Subtract pedestal
  std::vector<Double_t> wfPS(N);
  for (Int_t i = 0; i < N; ++i) wfPS[i] = samples[i] - result.pedestal;

  TemplateAnalysisConfig localCfg = cfg;
  if (localCfg.fNoiseSigmaADC < 0.0) {
    localCfg.fNoiseSigmaADC = pedNoiseSigma;
  }

  // 2. Matched filter
  std::vector<Double_t> mfOut = computeMatchedFilter(wfPS, tmpl, cfg);

  // 3. Peak finding
  std::vector<MFPeak> mfPeaks = findMFPeaks(mfOut, cfg);

  // Sort by MF amplitude (strongest first)
  std::sort(mfPeaks.begin(), mfPeaks.end(),
            [](const MFPeak &a, const MFPeak &b) {
              return a.mfAmplitude > b.mfAmplitude;
            });

  // 4. Template fit with sequential subtraction (strongest-first peeling)
  std::vector<Double_t> wfResidual = wfPS;

  Double_t templateArea = 0.0;
  for (Double_t val : tmpl.shape) {
    templateArea += val;
  }

  std::vector<Double_t> fittedPositions;

  for (const auto &peak : mfPeaks) {
    if (result.nPulses >= cfg.fMaxPulses) break;

    Bool_t tooClose = kFALSE;
    for (Double_t pos : fittedPositions) {
      if (std::abs(peak.position - pos) < static_cast<Double_t>(cfg.fMFSeparation)) {
        tooClose = kTRUE;
        break;
      }
    }
    if (tooClose) continue;

    Int_t fitLo = std::max(0, static_cast<Int_t>(peak.position) - localCfg.fFitMargin);
    Int_t fitHi = std::min(N - 1, static_cast<Int_t>(peak.position) + localCfg.fFitMargin);

    Double_t candidateBin = peak.position - static_cast<Double_t>(tmpl.peakBin);

    FitPulseResult fit = fitSingleTemplate(
        wfResidual, tmpl,
        candidateBin,
        fitLo, fitHi,
        localCfg);

    fit.mfAmplitude = peak.mfAmplitude;

    if (!fit.valid || fit.amplitude < localCfg.fAmpCutADC) {
      continue;
    }

    computePulseDerivedQuantities(fit, wfResidual, tmpl, templateArea, localCfg);
    addPulseContribution(wfResidual, fit, tmpl, -1.0);

    fittedPositions.push_back(getPulsePeakPosition(fit, tmpl));
    result.pulses.push_back(fit);
    result.nPulses++;
  }

  // 5. Iterative re-fitting (coordinate descent for pileup resolution)
  if (localCfg.fRefitMaxIter > 0 && result.nPulses >= 2) {
    for (Int_t iter = 0; iter < localCfg.fRefitMaxIter; ++iter) {
      Double_t prevTotalChi2 = 0.0;
      for (const auto &p : result.pulses)
        prevTotalChi2 += p.chi2ndf;

      std::vector<Int_t> order(result.nPulses);
      for (Int_t i = 0; i < result.nPulses; ++i) order[i] = i;
      std::sort(order.begin(), order.end(),
                [&](Int_t a, Int_t b) {
                  return result.pulses[a].amplitude > result.pulses[b].amplitude;
                });

      for (Int_t idx : order) {
        FitPulseResult &pulse = result.pulses[idx];

        addPulseContribution(wfResidual, pulse, tmpl, +1.0);

        Double_t peakPos = getPulsePeakPosition(pulse, tmpl);
        Int_t fitLo = std::max(0,     static_cast<Int_t>(peakPos) - localCfg.fFitMargin);
        Int_t fitHi = std::min(N - 1, static_cast<Int_t>(peakPos) + localCfg.fFitMargin);

        Double_t candidateBin = peakPos - static_cast<Double_t>(tmpl.peakBin);

        FitPulseResult refit = fitSingleTemplate(
            wfResidual, tmpl,
            candidateBin,
            fitLo, fitHi,
            localCfg);

        refit.mfAmplitude = pulse.mfAmplitude;

        if (refit.valid && refit.amplitude >= localCfg.fAmpCutADC) {
          computePulseDerivedQuantities(refit, wfResidual, tmpl,
                                       templateArea, localCfg);
          pulse = refit;
        }

        addPulseContribution(wfResidual, pulse, tmpl, -1.0);
      }

      Double_t newTotalChi2 = 0.0;
      for (const auto &p : result.pulses)
        newTotalChi2 += p.chi2ndf;

      if (prevTotalChi2 > 0.0 &&
          std::abs(newTotalChi2 - prevTotalChi2) / prevTotalChi2
              < localCfg.fRefitConvergence) {
        break;
      }
    }
  }

  // 6. Post-fit quality cuts
  std::vector<FitPulseResult> cleanPulses;
  cleanPulses.reserve(result.pulses.size());
  for (const auto &p : result.pulses) {
    if (p.chi2ndf > localCfg.fMaxChi2NDF)
      continue;
    if (std::abs(p.baseline) > std::max(localCfg.fMaxBaselineFrac * p.amplitude,
                                        localCfg.fMaxBaselineAbs))
      continue;
    if (p.fitStage == 2 && (p.fittedSigma < 0.10 || p.fittedTau < 0.40))
      continue;
    cleanPulses.push_back(p);
  }
  result.pulses = cleanPulses;
  result.nPulses = static_cast<Int_t>(cleanPulses.size());

  // Sort results chronologically for output
  std::sort(result.pulses.begin(), result.pulses.end(),
            [](const FitPulseResult &a, const FitPulseResult &b) {
              return a.timeSample < b.timeSample;
            });

  // Diagnostics
  result.wfPS = std::move(wfPS);
  result.mfOut = std::move(mfOut);

  return result;
}

