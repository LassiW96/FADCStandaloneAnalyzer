#include <iostream>
#include "waveformGenerator.h"

int main() {
    // Define WaveformData
    WaveformData waveData;
    waveData.TotNSamples = 20;
    waveData.InitPedSamples = 4;
    waveData.Pedestal = 100;
    waveData.PedRange = 4;
    waveData.NPulses = 1;
    waveData.PulseAmp = 100;
    waveData.PulseWidth = 10;
    waveData.PulseSeparation = 10;

    waveformGenerator generator;
    generator.setWaveformData(waveData);
    auto noise = generator.noiseWaveform();
    for( int i = 0; i < noise.size(); i++ ) {
        std::cout << noise[i] << std::endl;
    }
    return 0;
}