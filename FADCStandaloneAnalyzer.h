#ifndef FADCSTANDALONEANALYZER_H
#define FADCSTANDALONEANALYZER_H

#include <vector>
#include <cmath>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"

// EMG (exponential-modified Gaussian) model parameters — self-contained, no waveformGenerator dependency.
// The EMG pulse shape f(t) = exp(s2/2t2 - t/tau) * (1 + erf((t - s2/tau) / (sigma*sqrt(2))))
// has a naturally smooth, rounded peak; no external smearing is needed.
// Defaults match the standard NPS detector response; override in test scripts as needed.
struct SPEParams {
    Double_t sigma = 5.0;   // Gaussian sigma controlling the rise time (ns)
    Double_t tau   = 10.0;  // Exponential decay constant controlling the fall time (ns)
};

//Data structures that change only per-module basis
//Configurations
struct FADCConfig {
    // configuration parameters - taken from db
    Int_t fNPedestalSamples;
    Double_t fSampThreshold;
    Int_t fNSB;
    Int_t fNSA;
    Int_t fNSAT;
    Int_t fMaxNPulses;
};

//Calibrations
struct FADCCalib {
    // Calibration parameters - these can be also taken from db at run time
    Double_t ChanTomV;
    Double_t pC_Conv;

    // These change per channel basis
    // Double_t tcal;
    // Double_t gain;
};

//Output Pulse Data
struct FADCPulse {
    Double_t fSampPed;
    Bool_t fHasMulti;
    Int_t fNSampPulses;
    Int_t fNPeakSamples;
    Double_t fPeakPedestalRatio;

    // Output Pulse Data
    std::vector<Double_t> fSampPulseInt;
    std::vector<Double_t> fSampPulseAmp;
    std::vector<Double_t> fSampPulseTime;
    std::vector<Double_t> fSampPulseIntPedSub;
    std::vector<Double_t> fSampPulseIntMOLLERRaw;
    std::vector<Double_t> fSampPulseIntMOLLERVal;
    std::vector<Double_t> fSampPulseTOT;
    std::vector<Int_t>    fSampPulsePileup;
};

// Fit based analysis data structures
// Parameters
struct TemplateAnalysisConfig {
    // Pedestal
    Int_t    fNPedSamples  = 4;      // samples used to estimate the pedestal

    // Matched-filter window  (mirrors mfstart/mfend/mfleft/mfright in fit.C)
    Int_t    fMFLeft       = 5;      // samples left  of template peak in MF kernel
    Int_t    fMFRight      = 5;      // samples right of template peak in MF kernel
    Int_t    fMFStart      = fMFLeft;     // first sample of the MF search window
    Int_t    fMFEnd        = 250;    // last  sample of the MF search window

    // Peak finding  (mirrors specthres, mfthres in fit.C)
    Double_t fMFThreshold  = 0.1;   // MF amplitude threshold (normalised units)
    Int_t    fMFSeparation = 3;     // minimum distance between two MF peaks (samples)
    Int_t    fMaxPulses    = 12;    // cap on the number of pulses to fit

    // Template fit scan  (mirrors the +-4-bin window in Fitwf2's SetParLimits)
    Double_t fScanHalf     = 4.0;   // +- scan range around the MF candidate (samples)
    Double_t fScanStep     = 0.25;  // step size of the t0 scan (samples)

    // Local fit window: how many samples around the MF candidate to include in the fit.
    // Should be large enough to contain the full pulse shape (~2–3 decay time constants)
    // but small enough to exclude other pulses and avoid baseline absorption.
    Int_t    fFitMargin    = 20;    // fit over [peak - margin, peak + margin]

    // Noise model  (mirrors the sqrt(|content|*4.096/2)/4.096 in fit.C line 895)
    Double_t fAdcPerMV     = 4.096; // ADC counts per mV
    Double_t fMinSigma     = 0.5;   // floor for per-sample sigma (ADC counts)

    // Electronic noise sigma (ADC counts).  Added in quadrature with the
    // shot-noise term in the per-sample error model:  σ² = σ_shot² + σ_elec².
    // Negative (default) = auto-estimate from pedestal sample variance.
    // Zero = disable (pure shot-noise model, legacy behaviour).
    Double_t fNoiseSigmaADC = -1.0;

    // Amplitude quality cut
    Double_t fAmpCutADC    = 1.0;   // discard fitted pulses below this amplitude

    // Two-stage fitting threshold
    // Stage 1 is the fast 2-parameter (A, b) fixed-template scan.
    // If its chi2/ndf exceeds this value Stage 2 (Minuit with all 5 parameters) is invoked.
    Double_t fStage2Chi2Threshold = 5.0;

    // Iterative re-fitting (coordinate descent for pileup resolution).
    // After the initial sequential peeling pass, each pulse is re-fitted
    // against a residual with all *other* pulses subtracted, cycling until
    // convergence.  Set fRefitMaxIter = 0 to disable (legacy single-pass).
    Int_t    fRefitMaxIter      = 3;     // max refinement iterations (0 = disable)
    Double_t fRefitConvergence  = 0.01;  // stop when |Δchi2_total/chi2_total| < this

    // ── Second-derivative peak finder ────────────────────────────────────────────
    // These fields are used by AnalyzeWithD2(). They have no effect on
    // the legacy AnalyzeWithTemplate() matched-filter path.

    // Gaussian pre-smoothing kernel half-width (samples).
    // The smoothing sigma is fSmoothHalfWidth / 2.0 samples.
    Int_t    fSmoothHalfWidth      = 3;

    // Significance level for the auto-scaled d2 threshold.
    // The threshold is computed as fD2ThresholdNSigma * sigma_d2, where sigma_d2
    // is the RMS of the second derivative in the pedestal region.
    Double_t fD2ThresholdNSigma    = 5.0;

    // Minimum separation between second-derivative minima within a cluster
    // to count as distinct pulses (samples).
    Int_t    fD2MinPeakSep         = 3;

};

// Second-derivative cluster: a contiguous region of interest (ROI) where the
// second derivative falls below the noise threshold, indicating one or more
// overlapping pulses.
struct D2Cluster {
    Int_t    roiStart;      // first sample index of the ROI
    Int_t    roiEnd;        // last  sample index of the ROI
    Int_t    nPulses;       // number of d2 minima (= pulse candidates) in this cluster
};

// FADC Template (kernal for Matched Filter)
struct FADCTemplate {

    std::vector<Double_t> shape;      // normalised pulse shape, shape[i] in [0,1]
    std::vector<Double_t> mfKernel;   // reversed window used in the matched filter
    Double_t              mfNorm;     // sum(mfKernel) used as the normalisation divisor
    Int_t                 peakBin;    // sample index of the template peak (shape=1)

    // Build template from SPEParams using the EMG (exponential-modified Gaussian) pulse shape.
    // Parameters:
    //   sp          - SPEParams (sigma, tau; set directly in test scripts)
    //   sampleRate  - ns per sample (e.g. 4.0 for the NPS fADC)
    //   nSamples    - template length in samples (e.g. 110)
    //   peakTime_ns - time of the EMG peak within the window (ns).
    //                 Pass -1.0 (default) to auto-compute: sigma^2/tau + 3*sigma.
    //   mfLeft      - MF kernel samples left  of the peak.
    //                 Pass -1 (default) to auto-compute: ceil(2*sigma/sampleRate).
    //   mfRight     - MF kernel samples right of the peak.
    //                 Pass -1 (default) to auto-compute: ceil(3*tau/sampleRate).
    static FADCTemplate fromWaveformParams(const SPEParams& sp,
                                           Double_t sampleRate  = 4.0,
                                           Int_t    nSamples    = 110,
                                           Double_t peakTime_ns = -1.0,  // -1 = auto: sigma^2/tau + 3*sigma
                                           Int_t    mfLeft      = -1,    // -1 = auto: ceil(2*sigma/sampleRate)
                                           Int_t    mfRight     = -1);   // -1 = auto: ceil(3*tau/sampleRate)

    // Build template directly from a measured reference pulse (e.g., from a
    // laser or LED calibration run). rawShape is a vector of ADC-count samples
    // representing one clean pulse; it is normalised to peak = 1 internally.
    // Use this factory when calibration data are available. Fall back to
    // fromWaveformParams() (EMG model) when they are not.
    //
    // Parameters:
    //   rawShape         - measured pulse shape (at least 2 samples, peak > 0)
    //   enable           - kFALSE returns an empty FADCTemplate (caller uses EMG fallback)
    //   mfLeft / mfRight - MF kernel half-widths in samples
    static FADCTemplate fromRawSamples(const std::vector<Double_t>& rawShape,
                                       Bool_t enable   = kTRUE,
                                       Int_t  mfLeft   = 5,
                                       Int_t  mfRight  = 5);

    // Sub-sample evaluation via linear interpolation.
    // Returns 0 for t outside [0, shape.size()-2].
    Double_t eval(Double_t t) const;
};

// Outputs
struct FitPulseResult {
    Double_t timeSample;    // fitted t0 in sample units
                            //   Stage 1: t0 s.t. tmpl.eval(i-t0) peaks at (t0+peakBin)
                            //   Stage 2: t0 s.t. EMG(i-t0,sigma,tau) peaks at (t0+sigma²/tau)
    Double_t timeHalfMax;   // half-max rising-edge crossing in sample units (CFD-style)
    Double_t amplitude;     // fitted peak ADC count (pedestal-subtracted)
                            //   Stage 1: directly = A from the 2-parameter scan
                            //   Stage 2: A_raw * emgPeakRaw (normalised to physical peak)
    Double_t integral;      // fitted integral of the pulse (ADC counts × samples)
    Double_t baseline;      // fitted baseline offset b (residual DC after ped sub)
    Double_t chi2ndf;       // chi2 / NDF of the final accepted fit
    Double_t chi2ndf_s1;    // chi2 / NDF of the Stage-1 fixed-template fit (diagnostic;
                            //   equals chi2ndf when Stage 1 was accepted)
    Double_t mfAmplitude;   // MF output value at the candidate peak
    Double_t fittedSigma;   // Stage-2 fitted Gaussian sigma in sample units (0.0 for Stage-1)
    Double_t fittedTau;     // Stage-2 fitted exponential tau in sample units  (0.0 for Stage-1)
    Int_t    fitStage;      // 1 = fixed-template scan accepted, 2 = Minuit fallback used
    Bool_t   valid;         // kFALSE if rejected (A < fAmpCutADC, det~0, etc.)
    Int_t    clusterID        = -1;  // -1 if MF path; >=0 for D2 cluster index
    Int_t    nPulsesInCluster = 1;   // how many pulses share this cluster (1 = isolated)
};

struct FADCPulseWF {
    Double_t                    pedestal;  // mean of the first fNPedSamples
    Int_t                       nPulses;   // number of valid fitted pulses
    std::vector<FitPulseResult> pulses;    // one entry per MF peak, sorted by time

    // Diagnostics
    std::vector<Double_t> wfPS;
    std::vector<Double_t> mfOut;
    std::vector<Double_t> d2Out;              // second-derivative trace (empty if MF path)
    std::vector<D2Cluster> clusters;          // D2 cluster info (empty if MF path)
};

class FADCStandaloneAnalyzer {
public:
  FADCStandaloneAnalyzer();
  ~FADCStandaloneAnalyzer() = default;

  // Configuration setters
  void setConfig(const FADCConfig& config) {fConfig = config;}
  void setCalib(const FADCCalib& calib) {fCalib = calib;}

  // Getters
  const FADCConfig& getConfig() const { return fConfig; }
  const FADCCalib&  getCalib()  const { return fCalib; }

  const FADCTemplate& getTemplate() const { return fTemplate; }

  static std::vector<Double_t> computeMatchedFilter(const std::vector<Double_t>&  wfPS,
                                                    const FADCTemplate&           tmpl,
                                                    const TemplateAnalysisConfig& cfg = TemplateAnalysisConfig{});

  // Main analysis function
  FADCPulse Analyze(const std::vector<Double_t>& samples,
                                 Double_t gain, Double_t tcal) const;

  
  // Fit based functions
  // Note: use FADCTemplate::fromWaveformParams(...) static factory to build a template.
  // Main Template Analysis Function
  FADCPulseWF AnalyzeWithTemplate(const std::vector<Double_t>& samples,
                                  const FADCTemplate&          tmpl,
                                  const TemplateAnalysisConfig& cfg = TemplateAnalysisConfig{}) const;

  // Second-derivative based analysis — fully template-free.
  // Peak finding via second derivative of the (optionally smoothed) waveform.
  // Fitting via direct EMG model (Minuit, no template shape required).
  static std::vector<Double_t> computeSecondDerivative(
          const std::vector<Double_t>& wfPS,
          Int_t smoothHalfWidth = 3);

  FADCPulseWF AnalyzeWithD2(const std::vector<Double_t>& samples,
                            const TemplateAnalysisConfig& cfg = TemplateAnalysisConfig{}) const;

private:
  FADCConfig fConfig;
  FADCCalib  fCalib;

  FADCTemplate fTemplate;
  // Helper functions (computeMatchedFilter, findMFPeaks, fitSingleTemplate) are
  // implemented as file-static free functions in FADCStandaloneAnalyzer.cxx.
};    // end FADCStandaloneAnalyzer

#endif // FADCSTANDALONEANALYZER_H