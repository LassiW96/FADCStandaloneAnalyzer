#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <algorithm> // std::min, std::max, std::sort, std::min_element, std::max_element
#include <cmath>     // std::sqrt, std::abs, std::exp
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error

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
// evalRawEMG — file-local helper used by EMGFitFunctor and AnalyzeWithTemplate.
// Evaluates the raw (unnormalized) EMG pulse shape at time offset dt,
// with Gaussian sigma and exponential tau, both in sample units.
// Returns 0 for unphysical parameters or negligibly small values.
// ---------------------------------------------------------------------------
static Double_t evalRawEMG(Double_t dt, Double_t sigma, Double_t tau) {
  if (sigma <= 0.0 || tau <= 0.0)
    return 0.0;
  Double_t s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
  Double_t erf_arg = (dt - (sigma * sigma) / tau) / (sigma * std::sqrt(2.0));
  Double_t val = std::exp(s2_2t2 - (dt / tau)) * (1.0 + std::erf(erf_arg));
  return (val > 1e-12) ? val : 0.0;
}

// ---------------------------------------------------------------------------
// EMGFitFunctor — chi2 callable for Stage-2 Minuit fit.
// Parameters: [0]=A (scale, raw), [1]=b (baseline), [2]=t0,
//             [3]=sigma (sample units), [4]=tau (sample units).
// Note: A here is the scale of the *unnormalized* EMG; the physical peak
// amplitude is computed post-convergence in fitStage2Template.
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


// ═══════════════════════════════════════════════════════════════════════════════
// Second-derivative peak finder — fully template-free analysis path
// ═══════════════════════════════════════════════════════════════════════════════

// ---------------------------------------------------------------------------
// gaussianSmooth — 1-D Gaussian convolution with reflect boundary.
// Kernel sigma = halfWidth / 2.0 samples; kernel extends ±halfWidth samples.
// ---------------------------------------------------------------------------
std::vector<Double_t> FADCStandaloneAnalyzer::gaussianSmooth(
    const std::vector<Double_t> &wf, Int_t halfWidth) {
  Int_t N = static_cast<Int_t>(wf.size());
  if (halfWidth <= 0 || N <= 1)
    return wf;

  // Build normalised Gaussian kernel
  Double_t sigma = halfWidth / 2.0;
  Int_t kW = 2 * halfWidth + 1;
  std::vector<Double_t> kernel(kW);
  Double_t kSum = 0.0;
  for (Int_t j = 0; j < kW; ++j) {
    Double_t x = static_cast<Double_t>(j - halfWidth);
    kernel[j] = std::exp(-0.5 * x * x / (sigma * sigma));
    kSum += kernel[j];
  }
  for (Int_t j = 0; j < kW; ++j)
    kernel[j] /= kSum;

  // Convolve with reflect boundary
  std::vector<Double_t> out(N);
  for (Int_t i = 0; i < N; ++i) {
    Double_t val = 0.0;
    for (Int_t j = 0; j < kW; ++j) {
      Int_t idx = i + j - halfWidth;
      // Reflect at boundaries
      if (idx < 0)
        idx = -idx;
      if (idx >= N)
        idx = 2 * (N - 1) - idx;
      idx = std::max(0, std::min(idx, N - 1)); // safety clamp
      val += wf[idx] * kernel[j];
    }
    out[i] = val;
  }
  return out;
}

// ---------------------------------------------------------------------------
// computeSecondDerivative — public static method.
// Smooths the waveform with a Gaussian kernel, then computes the
// central-difference second derivative: d2[i] = y[i-1] - 2*y[i] + y[i+1].
// Edge values are copied from the nearest interior value.
// ---------------------------------------------------------------------------
std::vector<Double_t> FADCStandaloneAnalyzer::computeSecondDerivative(
    const std::vector<Double_t> &wfPS, Int_t smoothHalfWidth) {
  Int_t N = static_cast<Int_t>(wfPS.size());
  if (N < 3)
    return std::vector<Double_t>(N, 0.0);

  // Optional Gaussian pre-smoothing
  std::vector<Double_t> smoothed;
  const std::vector<Double_t> *yPtr = &wfPS;
  if (smoothHalfWidth > 0) {
    smoothed = gaussianSmooth(wfPS, smoothHalfWidth);
    yPtr = &smoothed;
  }
  const std::vector<Double_t> &y = *yPtr;

  // Central-difference second derivative
  std::vector<Double_t> d2(N, 0.0);
  for (Int_t i = 1; i < N - 1; ++i) {
    d2[i] = y[i - 1] - 2.0 * y[i] + y[i + 1];
  }
  // Edge values: copy from nearest interior value
  d2[0] = d2[1];
  d2[N - 1] = d2[N - 2];

  return d2;
}

struct EMGSeeds {
  Double_t sigma;
  Double_t tau;
};

// ---------------------------------------------------------------------------
// estimateEMGSeeds — estimate σ and τ from the smoothed waveform shape
// around a D2-identified peak, before Minuit runs.
//
//   σ ← 10%–90% rise time on the rising edge:  σ_seed = t_rise / 2.56
//   τ ← time from peak to 1/e amplitude on the falling edge: τ_seed = t_decay
//
// For pileup clusters (nPulsesInCluster > 1) the τ search window is halved
// to avoid measuring into the neighbouring pulse.  Falls back to the
// hardcoded NPS defaults (σ=1.25, τ=2.5 samples) when edge measurements
// are not possible.
// ---------------------------------------------------------------------------
static EMGSeeds estimateEMGSeeds(const std::vector<Double_t> &wfSmooth,
                                 Double_t peakPosition, Double_t peakAmplitude,
                                 Int_t nPulsesInCluster) {
  const Double_t defaultSigma = 1.25; // σ = 5 ns / 4 ns/sample
  const Double_t defaultTau = 2.5;    // τ = 10 ns / 4 ns/sample

  Int_t N = static_cast<Int_t>(wfSmooth.size());
  Int_t peakBin = static_cast<Int_t>(std::round(peakPosition));
  peakBin = std::max(1, std::min(peakBin, N - 2));

  EMGSeeds seeds = {defaultSigma, defaultTau};

  if (peakAmplitude < 1e-6)
    return seeds; // no signal — use defaults

  // ── σ estimate: 10%–90% rise time ────────────────────────────────────
  Double_t level10 = 0.10 * peakAmplitude;
  Double_t level90 = 0.90 * peakAmplitude;

  // Walk backward from the peak to find the crossings
  for (Int_t i = peakBin; i > 0; --i) {
    if (wfSmooth[i] <= level90 && wfSmooth[i + 1] > level90) {
      // Linear interpolation for sub-sample 90% crossing
      Double_t frac90 =
          (level90 - wfSmooth[i]) / (wfSmooth[i + 1] - wfSmooth[i]);
      Double_t t90 = static_cast<Double_t>(i) + frac90;

      // Continue walking backward for the 10% crossing
      for (Int_t j = i; j > 0; --j) {
        if (wfSmooth[j] <= level10 && wfSmooth[j + 1] > level10) {
          Double_t frac10 =
              (level10 - wfSmooth[j]) / (wfSmooth[j + 1] - wfSmooth[j]);
          Double_t t10 = static_cast<Double_t>(j) + frac10;
          Double_t tRise = t90 - t10;
          if (tRise > 0.0) {
            seeds.sigma = tRise / 2.56;
          }
          break;
        }
      }
      break;
    }
  }

  // Clamp σ to Minuit limits
  seeds.sigma = std::max(0.1, std::min(seeds.sigma, 5.0));

  // ── τ estimate: peak to 1/e decay ────────────────────────────────────
  Double_t level1e = peakAmplitude / std::exp(1.0); // ≈ 36.8% of peak

  // Limit the forward search to avoid measuring into a pileup neighbour
  Int_t maxForward = N - 1 - peakBin;
  if (nPulsesInCluster > 1) {
    maxForward = std::max(1, maxForward / 2);
  }
  Int_t searchEnd = std::min(N - 1, peakBin + maxForward);

  for (Int_t i = peakBin; i < searchEnd; ++i) {
    if (wfSmooth[i] >= level1e && wfSmooth[i + 1] < level1e) {
      // Linear interpolation for sub-sample 1/e crossing
      Double_t frac = (wfSmooth[i] - level1e) / (wfSmooth[i] - wfSmooth[i + 1]);
      Double_t tDecay =
          (static_cast<Double_t>(i) + frac) - static_cast<Double_t>(peakBin);
      if (tDecay > 0.0) {
        seeds.tau = tDecay;
      }
      break;
    }
  }

  // Clamp τ to Minuit limits
  seeds.tau = std::max(0.5, std::min(seeds.tau, 15.0));

  return seeds;
}

// File-local struct for D2 peak candidates (not exposed in header)
struct D2Peak {
  Double_t
      position; // sub-sample peak position (parabolic-refined on smoothed wf)
  Double_t
      d2Value; // d2 value at the minimum (negative = concave down = pulse peak)
  Int_t clusterID;           // which ROI cluster this peak belongs to
  Double_t seedSigma = 1.25; // per-pulse σ estimate from waveform edges
  Double_t seedTau = 2.5;    // per-pulse τ estimate from waveform edges
};

// Return type for findD2Peaks (avoids <utility> dependency for std::pair)
struct D2FindResult {
  std::vector<D2Peak> peaks;
  std::vector<D2Cluster> clusters;
};

// ---------------------------------------------------------------------------
// findD2Peaks — identify ROIs and pulse candidates from the d2 trace.
//
// Algorithm:
//   1. Scan d2 for contiguous regions where d2[i] < -d2Threshold (ROIs).
//   2. Within each ROI, find local minima of d2 separated by at least
//      minPeakSep samples.  Each minimum marks one pulse candidate.
//   3. For each d2 minimum, find the corresponding waveform peak (local max
//      of wfSmooth) and refine its position with parabolic interpolation.
//   4. Group peaks into D2Cluster structs with ROI boundaries and pulse count.
// ---------------------------------------------------------------------------
static D2FindResult findD2Peaks(const std::vector<Double_t> &d2,
                                const std::vector<Double_t> &wfSmooth,
                                Double_t d2Threshold, Int_t minPeakSep,
                                Int_t nPedSamples) {
  Int_t N = static_cast<Int_t>(d2.size());
  D2FindResult result;

  // 1. Find contiguous ROIs where d2 < -threshold
  struct ROI {
    Int_t start;
    Int_t end;
  };
  std::vector<ROI> rois;

  Bool_t inROI = kFALSE;
  Int_t roiStart = 0;

  for (Int_t i = 1; i < N; ++i) {
    if (!inROI && d2[i] < -d2Threshold) {
      inROI = kTRUE;
      roiStart = i;
    } else if (inROI && d2[i] >= -d2Threshold) {
      inROI = kFALSE;
      rois.push_back({roiStart, i - 1});
    }
  }
  // Close any open ROI at the end of the waveform
  if (inROI)
    rois.push_back({roiStart, N - 1});

  // 2. Process each ROI into a cluster
  for (size_t ci = 0; ci < rois.size(); ++ci) {
    const ROI &roi = rois[ci];
    D2Cluster cluster;
    cluster.roiStart = roi.start;
    cluster.roiEnd = roi.end;
    cluster.nPulses = 0;

    // Ensure the search range allows local-minimum detection (i-1, i+1 valid)
    Int_t searchStart = std::max(1, roi.start);
    Int_t searchEnd = std::min(N - 2, roi.end);

    if (searchStart > searchEnd) {
      // ROI too narrow for local-minimum detection — use the deepest bin
      Int_t bestBin = roi.start;
      for (Int_t i = roi.start + 1; i <= roi.end && i < N; ++i) {
        if (d2[i] < d2[bestBin])
          bestBin = i;
      }

      // Find corresponding waveform peak near the d2 minimum
      Int_t wfPeakBin = bestBin;
      Int_t sLo = std::max(0, bestBin - minPeakSep);
      Int_t sHi = std::min(N - 1, bestBin + minPeakSep);
      for (Int_t i = sLo; i <= sHi; ++i) {
        if (wfSmooth[i] > wfSmooth[wfPeakBin])
          wfPeakBin = i;
      }

      D2Peak peak;
      peak.position = static_cast<Double_t>(wfPeakBin);
      peak.d2Value = d2[bestBin];
      peak.clusterID = static_cast<Int_t>(ci);

      // Estimate seeds for a single isolated pulse
      EMGSeeds seeds =
          estimateEMGSeeds(wfSmooth, peak.position, wfSmooth[wfPeakBin], 1);
      peak.seedSigma = seeds.sigma;
      peak.seedTau = seeds.tau;

      result.peaks.push_back(peak);
      cluster.nPulses = 1;
      result.clusters.push_back(cluster);
      continue;
    }

    // Find all local minima of d2 within the ROI
    std::vector<Int_t> minBins;
    for (Int_t i = searchStart; i <= searchEnd; ++i) {
      if ((d2[i] < d2[i - 1] || i == searchStart) && d2[i] <= d2[i + 1]) {
        // Enforce minimum separation — keep the deeper minimum
        Bool_t merged = kFALSE;
        for (size_t k = 0; k < minBins.size(); ++k) {
          if (std::abs(i - minBins[k]) < minPeakSep) {
            if (d2[i] < d2[minBins[k]])
              minBins[k] = i;
            merged = kTRUE;
            break;
          }
        }
        if (!merged)
          minBins.push_back(i);
      }
    }

    // Fallback: if no strict local minimum was found, use the global minimum
    if (minBins.empty()) {
      Int_t bestBin = searchStart;
      for (Int_t i = searchStart + 1; i <= searchEnd; ++i) {
        if (d2[i] < d2[bestBin])
          bestBin = i;
      }
      minBins.push_back(bestBin);
    }

    // 3. For each d2 minimum, refine the peak position on the smoothed waveform
    for (Int_t mb : minBins) {
      // The d2 minimum corresponds to a waveform peak.
      // Find the local maximum of wfSmooth near this d2 minimum.
      Int_t sLo = std::max(0, mb - minPeakSep);
      Int_t sHi = std::min(N - 1, mb + minPeakSep);
      Int_t wfPeakBin = mb;
      Double_t wfPeakVal = wfSmooth[mb];
      for (Int_t i = sLo; i <= sHi; ++i) {
        if (wfSmooth[i] > wfPeakVal) {
          wfPeakVal = wfSmooth[i];
          wfPeakBin = i;
        }
      }

      // Parabolic sub-sample interpolation on the smoothed waveform:
      //   pos = peakBin - (y[peak+1] - y[peak-1]) /
      //                   (2 * (y[peak-1] - 2*y[peak] + y[peak+1]))
      Double_t pos = static_cast<Double_t>(wfPeakBin);
      if (wfPeakBin > 0 && wfPeakBin < N - 1) {
        Double_t denom =
            2.0 * (wfSmooth[wfPeakBin - 1] - 2.0 * wfSmooth[wfPeakBin] +
                   wfSmooth[wfPeakBin + 1]);
        if (std::abs(denom) > 1e-9) {
          pos = wfPeakBin -
                (wfSmooth[wfPeakBin + 1] - wfSmooth[wfPeakBin - 1]) / denom;
        }
      }

      D2Peak peak;
      peak.position = pos;
      peak.d2Value = d2[mb];
      peak.clusterID = static_cast<Int_t>(ci);

      // Determine how many pulses share this cluster to constrain the tau
      // search
      Int_t nPulsesInCluster = static_cast<Int_t>(minBins.size());
      EMGSeeds seeds = estimateEMGSeeds(wfSmooth, peak.position, wfPeakVal,
                                        nPulsesInCluster);
      peak.seedSigma = seeds.sigma;
      peak.seedTau = seeds.tau;

      result.peaks.push_back(peak);
      cluster.nPulses++;
    }

    result.clusters.push_back(cluster);
  }

  return result;
}

// ---------------------------------------------------------------------------
// fitEMGDirect — template-free EMG Minuit fit for the D2 path.
//
// wfPeakBin is the D2-identified waveform-peak sample position.
// The fit seeds t0 so that the initial EMG peaks near wfPeakBin:
//   t0_init = wfPeakBin - sigma_init²/tau_init
//
// Post-convergence, the raw Minuit amplitude A_raw is converted to the
// physical peak ADC count:  amplitude = A_raw * max(emgRaw) over the fit
// window.
// ---------------------------------------------------------------------------
static FitPulseResult fitEMGDirect(const std::vector<Double_t> &wfResidual,
                                   Double_t wfPeakBin, Int_t fitLo, Int_t fitHi,
                                   const TemplateAnalysisConfig &cfg,
                                   Double_t seedSigma = 1.25,
                                   Double_t seedTau = 2.5) {
  Int_t N = static_cast<Int_t>(wfResidual.size());
  fitLo = std::max(fitLo, 0);
  fitHi = std::min(fitHi, N - 1);
  Int_t nPts = fitHi - fitLo + 1;

  FitPulseResult best;
  best.valid = kFALSE;
  best.chi2ndf = std::numeric_limits<Double_t>::max();
  best.chi2ndf_s1 = std::numeric_limits<Double_t>::max();
  best.timeSample = wfPeakBin;
  best.timeHalfMax = wfPeakBin;
  best.amplitude = 0.0;
  best.integral = 0.0;
  best.baseline = 0.0;
  best.mfAmplitude = 0.0;
  best.fitStage = 2; // always EMG in D2 path
  best.fittedSigma = 0.0;
  best.fittedTau = 0.0;
  best.clusterID = -1;
  best.nPulsesInCluster = 1;

  std::unique_ptr<ROOT::Math::Minimizer> min(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!min)
    return best;

  min->SetMaxFunctionCalls(1000);
  min->SetTolerance(0.1);
  min->SetPrintLevel(0);

  EMGFitFunctor chi2Functor{&wfResidual,   fitLo,         fitHi,
                            cfg.fMinSigma, cfg.fAdcPerMV, cfg.fNoiseSigmaADC};
  ROOT::Math::Functor functor(chi2Functor, 5);
  min->SetFunction(functor);

  // Rough amplitude seed from the residual window
  Double_t maxVal = 0.0;
  for (Int_t i = fitLo; i <= fitHi; ++i)
    if (wfResidual[i] > maxVal)
      maxVal = wfResidual[i];

  // Dynamic seeds from waveform shape (or defaults if not estimated)
  const Double_t initSigma = seedSigma;
  const Double_t initTau = seedTau;

  // Seed t0 so the initial EMG peaks near the D2-identified waveform peak.
  // The raw EMG peaks at dt = sigma²/tau, so t0_init = wfPeakBin - sigma²/tau.
  Double_t initT0 = wfPeakBin - (initSigma * initSigma / initTau);

  min->SetLimitedVariable(0, "Amplitude", maxVal, maxVal * 0.1, 0.0,
                          maxVal * 20.0);
  min->SetLimitedVariable(1, "Baseline", 0.0, 1.0, -20.0, 20.0);
  min->SetLimitedVariable(2, "t0", initT0, 0.5, initT0 - cfg.fScanHalf,
                          initT0 + cfg.fScanHalf);
  min->SetLimitedVariable(3, "Sigma", initSigma, 0.1, 0.1, 5.0);
  min->SetLimitedVariable(4, "Tau", initTau, 0.2, 0.5, 15.0);

  min->Minimize();

  if (min->Status() == 0 || min->Status() == 1) {
    const Double_t *xs = min->X();
    Double_t A_raw = xs[0];
    Double_t t0 = xs[2];
    Double_t sigma = xs[3];
    Double_t tau = xs[4];

    // Convert raw amplitude to physical peak ADC count
    Double_t emgPeakRaw = 0.0;
    for (Int_t i = fitLo; i <= fitHi; ++i) {
      Double_t f = evalRawEMG(static_cast<Double_t>(i) - t0, sigma, tau);
      if (f > emgPeakRaw)
        emgPeakRaw = f;
    }
    if (emgPeakRaw < 1e-12)
      emgPeakRaw = 1.0;

    best.amplitude = A_raw * emgPeakRaw;
    best.baseline = xs[1];
    best.timeSample = t0;
    best.fittedSigma = sigma;
    best.fittedTau = tau;

    Double_t chi2 = min->MinValue();
    Int_t ndf = std::max(nPts - 5, 1);
    best.chi2ndf = chi2 / ndf;
    best.valid = kTRUE;
  }

  return best;
}

// ---------------------------------------------------------------------------
// addEMGContribution — add or subtract a fitted EMG pulse from a waveform.
// Template-free version of addPulseContribution for the D2 path.
// sign = +1.0 to add back (restore), -1.0 to subtract (peel).
// ---------------------------------------------------------------------------
static void addEMGContribution(std::vector<Double_t> &wf,
                               const FitPulseResult &pulse, Double_t sign) {
  Int_t N = static_cast<Int_t>(wf.size());
  // First pass: find the peak of the raw EMG for normalisation
  Double_t emgPeakRaw = 0.0;
  for (Int_t i = 0; i < N; ++i) {
    Double_t f = evalRawEMG(static_cast<Double_t>(i) - pulse.timeSample,
                            pulse.fittedSigma, pulse.fittedTau);
    if (f > emgPeakRaw)
      emgPeakRaw = f;
  }
  if (emgPeakRaw < 1e-12)
    emgPeakRaw = 1.0;
  // Second pass: add/subtract the normalised pulse shape
  for (Int_t i = 0; i < N; ++i) {
    Double_t f = evalRawEMG(static_cast<Double_t>(i) - pulse.timeSample,
                            pulse.fittedSigma, pulse.fittedTau) /
                 emgPeakRaw;
    wf[i] += sign * pulse.amplitude * f;
  }
}

// ---------------------------------------------------------------------------
// computePulseDerivedQuantitiesEMG — template-free derived quantities.
// Computes the pulse integral and half-max rising-edge timing for an
// EMG-fitted pulse.  The integration window is derived from the fitted
// sigma and tau rather than from a template shape.
// ---------------------------------------------------------------------------
static void
computePulseDerivedQuantitiesEMG(FitPulseResult &fit,
                                 const std::vector<Double_t> &wfForTiming,
                                 const TemplateAnalysisConfig &cfg) {
  Int_t N = static_cast<Int_t>(wfForTiming.size());

  // ── Integral ──────────────────────────────────────────────────────────
  // Integration window: [-3σ, peak + 5τ] centred on the EMG origin.
  Double_t sigma = (fit.fittedSigma > 0.0) ? fit.fittedSigma : 1.25;
  Double_t tau = (fit.fittedTau > 0.0) ? fit.fittedTau : 2.5;
  Double_t emgPeakOffset = sigma * sigma / tau;

  Int_t intLo = static_cast<Int_t>(std::floor(-3.0 * sigma));
  Int_t intHi = static_cast<Int_t>(std::ceil(emgPeakOffset + 5.0 * tau));

  Double_t emgAreaRaw = 0.0;
  Double_t emgPeakRaw = 0.0;
  for (Int_t i = intLo; i <= intHi; ++i) {
    Double_t f =
        evalRawEMG(static_cast<Double_t>(i), fit.fittedSigma, fit.fittedTau);
    emgAreaRaw += f;
    if (f > emgPeakRaw)
      emgPeakRaw = f;
  }
  if (emgPeakRaw < 1e-12)
    emgPeakRaw = 1.0;
  fit.integral = fit.amplitude * emgAreaRaw / emgPeakRaw;

  // ── Half-max rising-edge crossing (CFD-style timing) ──────────────────
  Int_t wfPeak = static_cast<Int_t>(std::round(
      fit.timeSample + fit.fittedSigma * fit.fittedSigma / fit.fittedTau));
  wfPeak = std::max(0, std::min(wfPeak, N - 1));
  Int_t hLo = std::max(0, wfPeak - cfg.fFitMargin);
  Int_t hHi = std::min(N - 1, wfPeak + cfg.fFitMargin);
  Double_t locPeak = 0.0;
  Int_t locBin = wfPeak;
  for (Int_t i = hLo; i <= hHi; ++i) {
    if (wfForTiming[i] > locPeak) {
      locPeak = wfForTiming[i];
      locBin = i;
    }
  }
  Double_t vMid = locPeak / 2.0;
  Double_t halfMax = static_cast<Double_t>(locBin); // fallback: peak bin
  for (Int_t i = hLo; i < locBin; ++i) {
    if (wfForTiming[i] <= vMid && wfForTiming[i + 1] > vMid) {
      halfMax =
          i + (vMid - wfForTiming[i]) / (wfForTiming[i + 1] - wfForTiming[i]);
      break;
    }
  }
  fit.timeHalfMax = halfMax;
}

// ---------------------------------------------------------------------------
// AnalyzeWithD2 — second-derivative based, fully template-free analysis.
//
// Pipeline:
//   1. Pedestal estimation (identical to AnalyzeWithTemplate)
//   2. Gaussian smoothing + second derivative computation
//   3. Auto-scaled d2 threshold from pedestal noise
//   4. ROI detection and peak identification via findD2Peaks
//   5. Direct EMG Minuit fit per peak (no template)
//   6. Iterative coordinate-descent re-fitting for pileup resolution
//   7. Chronological sort of output pulses
// ---------------------------------------------------------------------------
FADCPulseWF
FADCStandaloneAnalyzer::AnalyzeWithD2(const std::vector<Double_t> &samples,
                                      const TemplateAnalysisConfig &cfg) const {
  FADCPulseWF result;
  result.nPulses = 0;

  Int_t N = static_cast<Int_t>(samples.size());
  if (N == 0)
    return result;

  // ── 1. Pedestal — noise-aware front vs tail selection ─────────────────
  //       (identical logic to AnalyzeWithTemplate)
  Int_t nPed = std::min(cfg.fNPedSamples, N);

  // Front pedestal
  Double_t ped = 0.0;
  for (Int_t i = 0; i < nPed; ++i)
    ped += samples[i];
  ped /= nPed;
  Double_t pedVar = 0.0;
  for (Int_t i = 0; i < nPed; ++i) {
    Double_t d = samples[i] - ped;
    pedVar += d * d;
  }
  pedVar /= (nPed > 1) ? (nPed - 1) : 1;

  // Tail pedestal
  Double_t ped2 = 0.0;
  for (Int_t i = N - nPed; i < N; ++i)
    ped2 += samples[i];
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

  // Pedestal-subtracted waveform
  std::vector<Double_t> wfPS(N);
  for (Int_t i = 0; i < N; ++i)
    wfPS[i] = samples[i] - result.pedestal;
  result.wfPS = wfPS; // store for diagnostics

  // Auto-estimate electronic noise for the per-sample error model
  TemplateAnalysisConfig localCfg = cfg;
  if (localCfg.fNoiseSigmaADC < 0.0) {
    localCfg.fNoiseSigmaADC = pedNoiseSigma;
  }

  // ── 2. Smooth + second derivative ─────────────────────────────────────
  std::vector<Double_t> wfSmooth =
      gaussianSmooth(wfPS, localCfg.fSmoothHalfWidth);
  result.wfSmooth = wfSmooth;
  std::vector<Double_t> d2 =
      computeSecondDerivative(wfPS, localCfg.fSmoothHalfWidth);
  result.d2Out = d2;

  // ── 3. Auto-scale d2 threshold from pedestal noise ────────────────────
  //       Estimate σ_d2 from the pedestal region of the d2 trace.
  Int_t d2PedStart, d2PedEnd;
  if (ped <= ped2) {
    // Front is quieter (no early pulse)
    d2PedStart = 1; // skip edge copy
    d2PedEnd = std::min(d2PedStart + nPed, N - 1);
  } else {
    // Tail is quieter (early pulse detected)
    d2PedEnd = N - 1;
    d2PedStart = std::max(1, d2PedEnd - nPed);
  }
  Int_t d2PedN = d2PedEnd - d2PedStart;
  if (d2PedN < 1)
    d2PedN = 1;

  Double_t d2Mean = 0.0;
  for (Int_t i = d2PedStart; i < d2PedEnd; ++i)
    d2Mean += d2[i];
  d2Mean /= d2PedN;

  Double_t d2Var = 0.0;
  for (Int_t i = d2PedStart; i < d2PedEnd; ++i) {
    Double_t dd = d2[i] - d2Mean;
    d2Var += dd * dd;
  }
  d2Var /= (d2PedN > 1) ? (d2PedN - 1) : 1;
  Double_t d2Sigma = std::sqrt(d2Var);

  Double_t d2Threshold = localCfg.fD2ThresholdNSigma * d2Sigma;
  if (d2Threshold < 0.1)
    d2Threshold = 0.1; // floor to prevent noise triggers

  // ── 4. Find D2 peaks and clusters ─────────────────────────────────────
  D2FindResult d2Result =
      findD2Peaks(d2, wfSmooth, d2Threshold, localCfg.fD2MinPeakSep, nPed);
  result.clusters = d2Result.clusters;

  // Sort peaks by |d2| value (most negative = strongest peak first)
  std::sort(
      d2Result.peaks.begin(), d2Result.peaks.end(),
      [](const D2Peak &a, const D2Peak &b) { return a.d2Value < b.d2Value; });

  // ── 5. Fit each peak with direct EMG (template-free) ──────────────────
  std::vector<Double_t> wfResidual = wfPS;
  std::vector<Double_t> fittedPositions;

  for (const auto &peak : d2Result.peaks) {
    if (result.nPulses >= localCfg.fMaxPulses)
      break;

    // Skip if too close to an already-fitted pulse
    Bool_t tooClose = kFALSE;
    for (Double_t pos : fittedPositions) {
      if (std::abs(peak.position - pos) <
          static_cast<Double_t>(localCfg.fD2MinPeakSep)) {
        tooClose = kTRUE;
        break;
      }
    }
    if (tooClose)
      continue;

    // Fit window around the D2-identified peak
    Int_t fitLo =
        std::max(0, static_cast<Int_t>(peak.position) - localCfg.fFitMargin);
    Int_t fitHi = std::min(N - 1, static_cast<Int_t>(peak.position) +
                                      localCfg.fFitMargin);

    Int_t minPts = 2 * localCfg.fFitMargin;
    if (fitHi - fitLo < minPts) {
      if (fitLo == 0)
        fitHi = std::min(N - 1, fitLo + minPts);
      else if (fitHi == N - 1)
        fitLo = std::max(0, fitHi - minPts);
    }

    FitPulseResult fit = fitEMGDirect(wfResidual, peak.position, fitLo, fitHi,
                                      localCfg, peak.seedSigma, peak.seedTau);

    fit.mfAmplitude = std::abs(peak.d2Value); // store |d2| as diagnostic
    fit.clusterID = peak.clusterID;
    fit.nPulsesInCluster =
        (peak.clusterID >= 0 &&
         peak.clusterID < static_cast<Int_t>(d2Result.clusters.size()))
            ? d2Result.clusters[peak.clusterID].nPulses
            : 1;

    if (!fit.valid || fit.amplitude < localCfg.fAmpCutADC) {
      continue;
    }

    // Compute derived quantities (template-free)
    computePulseDerivedQuantitiesEMG(fit, wfResidual, localCfg);

    // Subtract this pulse from the residual for subsequent fits
    addEMGContribution(wfResidual, fit, -1.0);

    fittedPositions.push_back(
        fit.timeSample + fit.fittedSigma * fit.fittedSigma / fit.fittedTau);
    result.pulses.push_back(fit);
    result.nPulses++;
  }

  // ── 6. Iterative re-fitting (coordinate descent for pileup resolution)
  //       Each pulse is added back, re-fitted against the clean residual
  //       (all other pulses subtracted), then the new fit is subtracted.
  if (localCfg.fRefitMaxIter > 0 && result.nPulses >= 2) {
    for (Int_t iter = 0; iter < localCfg.fRefitMaxIter; ++iter) {
      Double_t prevTotalChi2 = 0.0;
      for (const auto &p : result.pulses)
        prevTotalChi2 += p.chi2ndf;

      // Process pulses strongest-first
      std::vector<Int_t> order(result.nPulses);
      for (Int_t i = 0; i < result.nPulses; ++i)
        order[i] = i;
      std::sort(order.begin(), order.end(), [&](Int_t a, Int_t b) {
        return result.pulses[a].amplitude > result.pulses[b].amplitude;
      });

      for (Int_t idx : order) {
        FitPulseResult &pulse = result.pulses[idx];

        // Add back this pulse's contribution to the residual
        addEMGContribution(wfResidual, pulse, +1.0);

        // Fit window from the current pulse's peak position
        Double_t peakPos = pulse.timeSample + pulse.fittedSigma *
                                                  pulse.fittedSigma /
                                                  pulse.fittedTau;
        Int_t fitLo =
            std::max(0, static_cast<Int_t>(peakPos) - localCfg.fFitMargin);
        Int_t fitHi =
            std::min(N - 1, static_cast<Int_t>(peakPos) + localCfg.fFitMargin);

        // Re-fit against the clean residual
        FitPulseResult refit =
            fitEMGDirect(wfResidual, peakPos, fitLo, fitHi, localCfg,
                         pulse.fittedSigma, pulse.fittedTau);

        refit.mfAmplitude = pulse.mfAmplitude;
        refit.clusterID = pulse.clusterID;
        refit.nPulsesInCluster = pulse.nPulsesInCluster;

        if (refit.valid && refit.amplitude >= localCfg.fAmpCutADC) {
          computePulseDerivedQuantitiesEMG(refit, wfResidual, localCfg);
          pulse = refit;
        }

        // Subtract the (possibly updated) pulse from the residual
        addEMGContribution(wfResidual, pulse, -1.0);
      }

      // Convergence check
      Double_t newTotalChi2 = 0.0;
      for (const auto &p : result.pulses)
        newTotalChi2 += p.chi2ndf;

      if (prevTotalChi2 > 0.0 &&
          std::abs(newTotalChi2 - prevTotalChi2) / prevTotalChi2 <
              localCfg.fRefitConvergence) {
        break;
      }
    }
  }

  // ── 7. Sort results chronologically for output ────────────────────────
  std::sort(result.pulses.begin(), result.pulses.end(),
            [](const FitPulseResult &a, const FitPulseResult &b) {
              return a.timeSample < b.timeSample;
            });

  return result;
}