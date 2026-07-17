#ifndef WAVEFORMGENERATOR_H
#define WAVEFORMGENERATOR_H

#include <vector>
#include <random>
#include "Rtypes.h"
#include "RtypesCore.h"
#include "TMath.h"

#include "FADCStandaloneAnalyzer.h"

struct WaveformData {
    Int_t TotNSamples;
    Int_t InitPedSamples;

    Int_t Pedestal;
    Int_t PedRange;

    Int_t NPulses;
    Int_t PulseAmp;
    Int_t PulseWidth;
    Int_t PulseSeparation;
};

struct WaveformInfo {
    std::vector<Double_t> PedSubIntegral; // Double to match analyzer's floating-point formula
    std::vector<Int_t>    PulseAmp;
    std::vector<Int_t>    PulseTime;
};

// Structs for the waveform generator algorithm v2 - keep the safety values
struct WaveformParams {
    Double_t tau_r = 2.0;           // Rise time (ns)
    Double_t tau_f = 10.0;          // Decay time (ns)
    Double_t tts_sigma = 1.5;       // Transit time spread (ns)
    Double_t mean_gain = 5.0;       // Mean amplitude per PE (mV)
    Double_t gain_sigma = 1.0;      // Gain variation (mV)

    // After pulse parameters
    Double_t ap_prob = 0.03;         // After pulse generation probability
    Double_t ap_delay_mean = 80.0;   // Mean delay for afterpulse (ns)
    Double_t ap_delay_sigma = 5.0;   // Spread in afterpulse delay (ns)

};

struct DigitizerParams {
    Double_t sample_rate = 4.0;      // Sampling rate of the digitizer (250 MHz)
    Int_t num_samples = 256;         // #of samples in the sampling window
    Double_t pedestal = 100.0;       // Pedestal in ADC counts
    Double_t noise_sigma = 2.5;      // Electronic noise (mV)
    Double_t v_lsb = 1.0;            // Voltage resolution (mV per ADC count)
    Int_t max_adc = 4095;            // Maximum ADC value (12-bit)
    Int_t max_n_pulses = 4;          // Max # of pulses per window
};

class waveformGenerator {
    public:
        waveformGenerator(); // Seed the random number generator. New seed everytime run.
        virtual ~waveformGenerator() = default;

        // Waveform Data Setters
        void setWaveformData(const WaveformData& waveformData) {fWaveformData = waveformData;}
        void setAnalyzer(const FADCStandaloneAnalyzer& analyzer) {saAnalyzer = analyzer;}

        // NEW: Detector + Digitizer Setter
        void setDetectorParams(const WaveformParams& params) {fWf = params;}
        void setDigitizerParams(const DigitizerParams& params) {fDaq = params;}

        // Waveform generators
        // Simple generators
        std::vector<Int_t> noiseWaveform() const;
        std::vector<Int_t> singlePulseWaveform() const;
        std::vector<Int_t> multiplePulseWaveform() const;

        // Generator v2
        std::vector<Int_t> generateWaveform() const;

        // Calc pulse parameters
        //std::vector<Int_t> calcPulseIntegral(const std::vector<Int_t> wave) const;
        void calcPulseParameters(const std::vector<Int_t> wave);

        // Accessors for pulse parameters computed by calcPulseParameters()
        const WaveformInfo& getWaveformInfo() const { return fWaveformInfo; }
        void clearWaveformInfo() { fWaveformInfo = {}; }
        
    private:
        WaveformData fWaveformData;
        WaveformInfo fWaveformInfo;

        // Members for waveform generator v2
        DigitizerParams fDaq;
        WaveformParams fWf;

        FADCStandaloneAnalyzer saAnalyzer;

        // Per-instance RNG — thread-safe (no shared global state)
        mutable std::mt19937 fRng;

        // Helpers for waveform generator v2
        // SPE template
        Double_t spePulse(Double_t t, Double_t t0, Double_t amp) const;

        // Hit times + after pulse generator
        std::vector<Double_t> generatePETimes() const;
};

#endif // WAVEFORMGENERATOR_H