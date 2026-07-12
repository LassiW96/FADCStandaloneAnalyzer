#ifndef WAVEFORMGENERATOR_H
#define WAVEFORMGENERATOR_H

#include <vector>
#include <random>
#include <TRandom.h>
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
    std::vector<Int_t> PedSubIntegral;
    std::vector<Int_t> PulseAmp;
    std::vector<Int_t> PulseTime;
};

class waveformGenerator {
    public:
        waveformGenerator(); // Seed the random number generator. New seed everytime run.
        virtual ~waveformGenerator() = default;

        // Waveform Data Setters
        void setWaveformData(const WaveformData& waveformData) {fWaveformData = waveformData;}
        void setAnalyzer(const FADCStandaloneAnalyzer& analyzer) {saAnalyzer = analyzer;}

        // Waveform generators
        std::vector<Int_t> noiseWaveform() const;
        std::vector<Int_t> singlePulseWaveform() const;
        std::vector<Int_t> multiplePulseWaveform() const;

        // Calc pulse parameters
        //std::vector<Int_t> calcPulseIntegral(const std::vector<Int_t> wave) const;
        void calcPulseParameters(const std::vector<Int_t> wave);

        // Accessors for pulse parameters computed by calcPulseParameters()
        const WaveformInfo& getWaveformInfo() const { return fWaveformInfo; }
        void clearWaveformInfo() { fWaveformInfo = {}; }
        
    private:
        WaveformData fWaveformData;
        WaveformInfo fWaveformInfo;

        FADCStandaloneAnalyzer saAnalyzer;
};

#endif // WAVEFORMGENERATOR_H