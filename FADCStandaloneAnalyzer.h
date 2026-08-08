#ifndef FADCSTANDALONEANALYZER_H
#define FADCSTANDALONEANALYZER_H

#include <vector>
#include <cmath>

#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"

#include "waveformGenerator.h"  // for WaveformParams (tau_r, tau_f)

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

    // Amplitude quality cut
    Double_t fAmpCutADC    = 1.0;   // discard fitted pulses below this amplitude

    // ── Peak Optimization (added) ────────────────────────────────────────────
    // Parabolic sub-sample t0 refinement (Improvement 2).
    // After the chi² scan identifies the best grid point, a 3-point parabola is
    // fitted to the surrounding chi² values to pin down the true minimum between
    // scan steps. This reduces timing bias and the amplitude inflation that comes
    // from a slight t0 mismatch.
    Bool_t   fRefineT0 = kTRUE;   // kTRUE = apply parabolic refinement; kFALSE = scan-only
    // ─────────────────────────────────────────────────────────────────────────
};

// FADC Template (kernal for Matched Filter)
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
    // Build template from WaveformParams (SPE bi-exponential pulse shape).
    // smearSigma  - Gaussian convolution width in samples (Improvement 1);
    //               rounds the sharp mathematical peak to model bandwidth limits.
    //               Safe starting value for unknown bandwidth at 4 ns/sample: 0.5.
    //               Set to 0.0 to disable smearing and recover the original behaviour.
    static FADCTemplate fromWaveformParams(const WaveformParams& wp,
                                           Double_t sampleRate  = 4.0,
                                           Int_t    nSamples    = 110,
                                           Double_t peakTime_ns = 20.0,
                                           Int_t    mfLeft      = 5,
                                           Int_t    mfRight     = 5,
                                           Double_t smearSigma  = 0.5);  // ── Peak Optimization (added)

    // ── Peak Optimization: calibration-data factory (added) ──────────────────
    // Build template directly from a measured reference pulse (e.g., from a
    // laser or LED calibration run). rawShape is a vector of ADC-count samples
    // representing one clean pulse; it is normalised to peak = 1 internally.
    // Use this factory when calibration data are available. Fall back to
    // fromWaveformParams() (bi-exponential model) when they are not.
    //
    // Parameters:
    //   rawShape - measured pulse shape (at least 2 samples, peak > 0)
    //   enable   - kFALSE returns an empty FADCTemplate; the caller should then
    //              use fromWaveformParams() as backup.
    //   mfLeft / mfRight - MF kernel half-widths (same meaning as above)
    static FADCTemplate fromRawSamples(const std::vector<Double_t>& rawShape,
                                       Bool_t enable   = kTRUE,
                                       Int_t  mfLeft   = 5,
                                       Int_t  mfRight  = 5);
    // ─────────────────────────────────────────────────────────────────────────

    // Sub-sample evaluation via linear interpolation.
    // Returns 0 for t outside [0, shape.size()-2].
    Double_t eval(Double_t t) const;
};

// Outputs
struct FitPulseResult {
    Double_t timeSample;    // fitted t0 in sample units (template-fit χ² scan)
    Double_t timeHalfMax;   // half-max rising-edge crossing in sample units (CFD-style)
    Double_t amplitude;     // fitted amplitude A (pedestal-subtracted ADC counts)
    Double_t integral;      // fitted integral of the pulse
    Double_t baseline;      // fitted baseline offset b (residual DC after ped sub)
    Double_t chi2ndf;       // chi2 / NDF at the best t0
    Double_t mfAmplitude;   // MF output value at the candidate peak
    Bool_t   valid;         // kFALSE if rejected (A < fAmpCutADC, det~0, etc.)

    // ── Peak Optimization (added) ────────────────────────────────────────────
    // Observed local maximum of the residual waveform inside the fit window,
    // stored for diagnostic comparison with the fitted amplitude.
    // If amplitude >> peakAmplitude, the template peak is sharper than the data.
    Double_t peakAmplitude; // observed waveform max in [fitLo, fitHi] (diagnostic)
    // ─────────────────────────────────────────────────────────────────────────
};

struct FADCPulseWF {
    Double_t                    pedestal;  // mean of the first fNPedSamples
    Int_t                       nPulses;   // number of valid fitted pulses
    std::vector<FitPulseResult> pulses;    // one entry per MF peak, sorted by time

    // Diagnostics
    std::vector<Double_t> wfPS;
    std::vector<Double_t> mfOut;
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

private:
  FADCConfig fConfig;
  FADCCalib  fCalib;

  FADCTemplate fTemplate;
  // Helper functions (computeMatchedFilter, findMFPeaks, fitSingleTemplate) are
  // implemented as file-static free functions in FADCStandaloneAnalyzer.cxx.
};    // end FADCStandaloneAnalyzer

#endif // FADCSTANDALONEANALYZER_H