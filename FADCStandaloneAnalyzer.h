#ifndef FADCSTANDALONEANALYZER_H
#define FADCSTANDALONEANALYZER_H

#include <cmath>
#include <vector>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"

// Data structures that change only per-module basis
// Configurations
struct FADCConfig {
  // configuration parameters - taken from db
  Int_t fNPedestalSamples;
  Double_t fSampThreshold;
  Int_t fNSB;
  Int_t fNSA;
  Int_t fNSAT;
  Int_t fMaxNPulses;
};

// Calibrations
struct FADCCalib {
  // Calibration parameters - these can be also taken from db at run time
  Double_t ChanTomV;
  Double_t pC_Conv;
};

// Output pulse data structure for hcana style analysis
struct FADCPulse {
  Double_t fSampPed;
  Bool_t fHasMulti;
  Int_t fNSampPulses;
  Int_t fNPeakSamples;
  Double_t fPeakPedestalRatio;

  std::vector<Double_t> fSampPulseInt;
  std::vector<Double_t> fSampPulseAmp;
  std::vector<Double_t> fSampPulseTime;
  std::vector<Double_t> fSampPulseIntPedSub;
  std::vector<Double_t> fSampPulseIntMOLLERRaw;
  std::vector<Double_t> fSampPulseIntMOLLERVal;
  std::vector<Double_t> fSampPulseTOT;
  std::vector<Int_t> fSampPulsePileup;
};

// Fit based analysis data structures
// Parameters
struct TemplateAnalysisConfig {
  // Pedestal
  Int_t fNPedSamples = 6; // samples used to estimate the pedestal

  // Matched-filter window (for computeMatchedFilter with raw template)
  Int_t fMFLeft = 5;  // samples left of template peak in MF kernel
  Int_t fMFRight = 5; // samples right of template peak in MF kernel

  // Noise model
  Double_t fAdcPerMV = 4.096; // ADC counts per mV
  Double_t fMinSigma = 0.5;   // floor for per-sample sigma (ADC counts)

  // Electronic noise sigma (ADC counts).
  // Negative (default) = auto-estimate from pedestal sample variance.
  // Zero = disable (pure shot-noise model, legacy behaviour).
  Double_t fNoiseSigmaADC = -1.0;

  // Amplitude quality cut
  Double_t fAmpCutADC = 1.0; // discard fitted pulses below this amplitude
  Int_t fMaxPulses = 12;     // cap on the number of pulses to fit

  // Local fit window and t0 search range around candidate peak
  Double_t fScanHalf = 4.0; // +- scan range around candidate (samples)
  Int_t fFitMargin = 20;    // fit over [peak - margin, peak + margin]

  // Iterative re-fitting (coordinate descent for pileup resolution).
  Int_t fRefitMaxIter = 3; // max refinement iterations (0 = disable)
  Double_t fRefitConvergence =
      0.01; // stop when |Δchi2_total/chi2_total| < this

  // Post-fit quality cuts
  Double_t fMaxChi2NDF = 25.0;     // reject fits with chi2/ndf above this
  Double_t fMaxBaselineFrac = 0.3; // reject if |baseline| > frac * amplitude
  Double_t fMaxBaselineAbs = 5.0;  // absolute cap on baseline (ADC)

  // Second-derivative peak finder
  // Gaussian pre-smoothing kernel half-width (samples).
  Int_t fSmoothHalfWidth = 3;

  // Significance level for the auto-scaled d2 threshold.
  // The threshold is computed as fD2ThresholdNSigma * sigma_d2, where sigma_d2
  // is the RMS of the second derivative in the pedestal region.
  Double_t fD2ThresholdNSigma = 5.0;

  Int_t fD2MinPeakSep = 3;
};

struct D2Cluster {
  Int_t roiStart; // first sample index of the ROI
  Int_t roiEnd;   // last  sample index of the ROI
  Int_t nPulses;  // number of d2 minima (= pulse candidates) in this cluster
};

// FADC Template (kernel for Matched Filter or raw calibration reference)
struct FADCTemplate {
  std::vector<Double_t> shape;    // normalised pulse shape, shape[i] in [0,1]
  std::vector<Double_t> mfKernel; // reversed window used in the matched filter
  Double_t mfNorm; // sum(mfKernel) used as the normalisation divisor
  Int_t peakBin;   // sample index of the template peak (shape=1)

  // Build template directly from a measured reference pulse (e.g., from a
  // laser or LED calibration run). rawShape is a vector of ADC-count samples
  // representing one clean pulse; it is normalised to peak = 1 internally.
  //
  // Parameters:
  //   rawShape         - measured pulse shape (at least 2 samples, peak > 0)
  //   enable           - kFALSE returns an empty FADCTemplate
  //   mfLeft / mfRight - MF kernel half-widths in samples
  static FADCTemplate fromRawSamples(const std::vector<Double_t> &rawShape,
                                     Bool_t enable = kTRUE, Int_t mfLeft = 5,
                                     Int_t mfRight = 5);

  // Sub-sample evaluation via linear interpolation.
  // Returns 0 for t outside [0, shape.size()-2].
  Double_t eval(Double_t t) const;
};

// Output data structures for fit based/template based analysis
struct FitPulseResult {
  Double_t timeSample; // fitted t0 in sample units s.t. EMG(i-t0,sigma,tau)
                       // peaks at t0+sigma²/tau
  Double_t
      timeHalfMax; // half-max rising-edge crossing in sample units (CFD-style)
  Double_t
      amplitude; // fitted peak ADC count (pedestal-subtracted, physical peak)
  Double_t integral;    // fitted integral of the pulse (ADC counts × samples)
  Double_t baseline;    // fitted baseline offset b (residual DC after ped sub)
  Double_t chi2ndf;     // chi2 / NDF of the final accepted fit
  Double_t chi2ndf_s1;  // diagnostic chi2 / NDF
  Double_t mfAmplitude; // |d2| magnitude or MF output at peak
  Double_t fittedSigma; // fitted Gaussian sigma in sample units
  Double_t fittedTau;   // fitted exponential tau in sample units
  Int_t fitStage;       // fit stage indicator (2 = EMG Minuit fit)
  Bool_t valid;         // kFALSE if rejected (A < fAmpCutADC, det~0, etc.)
  Int_t clusterID = -1; // -1 if non-cluster; >=0 for D2 cluster index
  Int_t nPulsesInCluster =
      1; // how many pulses share this cluster (1 = isolated)
};

struct FADCPulseWF {
  Double_t pedestal; // mean of the first fNPedSamples
  Int_t nPulses;     // number of valid fitted pulses
  std::vector<FitPulseResult>
      pulses; // one entry per fitted pulse, sorted by time

  // Diagnostics
  std::vector<Double_t> wfPS;
  std::vector<Double_t>
      wfSmooth; // smoothed trace (for diagnostics / intermediate steps)
  std::vector<Double_t> mfOut;
  std::vector<Double_t> d2Out;     // second-derivative trace (empty if non-D2)
  std::vector<D2Cluster> clusters; // D2 cluster info
};

class FADCStandaloneAnalyzer {
public:
  FADCStandaloneAnalyzer();
  ~FADCStandaloneAnalyzer() = default;

  // Configuration setters
  void setConfig(const FADCConfig &config) { fConfig = config; }
  void setCalib(const FADCCalib &calib) { fCalib = calib; }

  // Getters
  const FADCConfig &getConfig() const { return fConfig; }
  const FADCCalib &getCalib() const { return fCalib; }

  const FADCTemplate &getTemplate() const { return fTemplate; }

  static std::vector<Double_t> computeMatchedFilter(
      const std::vector<Double_t> &wfPS, const FADCTemplate &tmpl,
      const TemplateAnalysisConfig &cfg = TemplateAnalysisConfig{});

  static std::vector<Double_t> gaussianSmooth(const std::vector<Double_t> &wf,
                                              Int_t halfWidth);

  static std::vector<Double_t>
  computeSecondDerivative(const std::vector<Double_t> &wfPS,
                          Int_t smoothHalfWidth = 3);

  // Main analysis functions
  FADCPulse Analyze(const std::vector<Double_t> &samples, Double_t gain,
                    Double_t tcal) const;

  // Second-derivative based analysis — fully template-free.
  FADCPulseWF AnalyzeWithD2(
      const std::vector<Double_t> &samples,
      const TemplateAnalysisConfig &cfg = TemplateAnalysisConfig{}) const;

private:
  FADCConfig fConfig;
  FADCCalib fCalib;

  FADCTemplate fTemplate;
}; // end FADCStandaloneAnalyzer

#endif // FADCSTANDALONEANALYZER_H