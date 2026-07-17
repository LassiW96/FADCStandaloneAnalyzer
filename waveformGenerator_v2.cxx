#include "waveformGenerator_v2.h"
#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

waveformGenerator_v2::waveformGenerator_v2() : fRng(std::random_device{}()) {}

std::vector<Int_t> waveformGenerator_v2::noiseWaveform() const {
    int min = fDaq.pedestal - fDaq.ped_sigma / 2.0;
    int max = fDaq.pedestal + fDaq.ped_sigma / 2.0 + 1;
    std::uniform_real_distribution<double> rndm(0.0, 1.0);
    std::vector<Int_t> wave;
    for( int i = 0; i < fDaq.num_samples; i++ ) {
    wave.push_back(min + (max - min) * rndm(fRng));
    }
    return wave;
}

std::vector<Int_t> waveformGenerator_v2::singlePulseWaveform() const {
    int min_ped = fDaq.pedestal - fDaq.ped_sigma / 2.0;
    int max_ped = fDaq.pedestal + fDaq.ped_sigma / 2.0 + 1;
    std::uniform_real_distribution<double> rndm(0.0, 1.0);
    std::vector<Int_t> wave;

    int pulse_min = fDaq.init_ped_samples;
    int pulse_max = fDaq.num_samples - 1;
    int random_idx = pulse_min + int(pulse_max - pulse_min) * rndm(fRng);
    for( int i = 0; i < fDaq.num_samples; i++ ) {
        if(i >= random_idx && i < random_idx + fWf.pulse_width) {
            int j = i - random_idx;
            wave.push_back(min_ped + fWf.mean_gain/TMath::Power(2, j));
        } else {
            wave.push_back(min_ped + (max_ped - min_ped) * rndm(fRng));
        }
    }
    return wave;
}

std::vector<Int_t> waveformGenerator_v2::multiplePulseWaveform() const {
    int min_ped = fDaq.pedestal - fDaq.ped_sigma / 2.0;
    int max_ped = fDaq.pedestal + fDaq.ped_sigma / 2.0 + 1;
    std::uniform_real_distribution<double> rndm(0.0, 1.0);
    std::vector<Int_t> wave;
    std::vector<std::pair<Int_t, Int_t>> pulse_info;

    int num_pulses = 2 + int((fDaq.max_n_pulses - 2) * rndm(fRng));

    // Maximum pulse space required
    int max_space = (num_pulses * fWf.pulse_width) + 
                ((num_pulses - 1) * fWf.pulse_separation);

    // Maximum available space for pulse start
    int max_first_start = fDaq.num_samples - max_space;

    if(max_first_start < fDaq.init_ped_samples) {
    max_first_start = fDaq.init_ped_samples;
    }

    // Current starting point
    int current_start = fDaq.init_ped_samples + 
                    int((max_first_start - fDaq.init_ped_samples) * rndm(fRng));

    // Random amps and separations
    for(int p = 0; p < num_pulses; p++) {
    int min_amp = fWf.mean_gain / 2;
    int current_amp = min_amp + int((fWf.mean_gain - min_amp) * rndm(fRng));

    pulse_info.push_back({current_start, current_amp});

    int current_separation = fWf.pulse_separation / 2 + 
                                int((fWf.pulse_separation - fWf.pulse_separation / 2.) * rndm(fRng));
    current_start += fWf.pulse_width + current_separation;
    }

    // Build the waveform
    for( int i = 0; i < fDaq.num_samples; i++) {
    Bool_t in_pulse = false;
    Int_t pulse_relative_idx = 0;
    Int_t active_amp = 0;

    for(auto& p : pulse_info) {
        if(i >= p.first && i < p.first + fWf.pulse_width) {
            in_pulse = true;
            pulse_relative_idx = i - p.first;
            active_amp = p.second;
            break;
        }
    }
    if(in_pulse) {
        wave.push_back(min_ped + active_amp/TMath::Power(2, pulse_relative_idx));
    } else {
        wave.push_back(min_ped + (max_ped - min_ped) * rndm(fRng));
    }
    }
    return wave;
}

void waveformGenerator_v2::calcPulseParameters(const std::vector<Int_t> wave) {
    FADCConfig fConfig = saAnalyzer.getConfig();
    Int_t NPedSamps = fConfig.fNPedestalSamples;
    Int_t NSB       = fConfig.fNSB;
    Int_t NSA       = fConfig.fNSA;
    Int_t Threshold = fConfig.fSampThreshold;
    Int_t NSamples  = (Int_t)wave.size();

    std::vector<Int_t> AvgIntegrals = {};

    // Pedestal: raw sum of first NPedSamps samples (matches analyzer's fSampPed)
    Double_t AvgPed = 0.0;
    Double_t PedSum = 0.0;
    for (int i = 0; i < NPedSamps; i++)
        PedSum += wave[i];
    // Use the same formula as the analyzer:
    //   AvgPed (float) = PedSum / NPedSamps  — NOT integer division

    // Calculate the pedestal using last 4 samples as well
    Double_t PedSumLastFour = 0.0;
    for( int i = NSamples - NPedSamps; i < NSamples; i++) {
        PedSumLastFour += wave[i];
    }
    
    if(PedSumLastFour < PedSum) {
        AvgPed = PedSumLastFour / NPedSamps;
    } else {
        AvgPed = PedSum / NPedSamps;
    }

    Bool_t PulseFound = false;

    for(int i = 0; i < wave.size(); i++) {
        if(wave[i] - AvgPed > Threshold) {
            if(!PulseFound) {
                // Integration window: [FirstSamp - NSB, FirstSamp + NSA - 1], clamped to valid range
                Int_t lo = TMath::Max(i - NSB, 0);
                Int_t hi = TMath::Min(i + NSA - 1, NSamples - 1);
                // Pedestal-subtracted integral, using the same formula as the analyzer:
                //   PedSubInt = RawInt - PedSum * (window_width / NPedSamps)
                //             = RawInt - AvgPed * window_width
                Int_t WindowWidth = hi - lo + 1;
                Double_t AvgIntegral = 0.0;
                Int_t RawIntegral = 0;
                for (int j = lo; j <= hi; j++) {
                    RawIntegral += wave[j];
                }
                // Match analyzer: PedSubInt = RawInt - PedSum * (WindowWidth / NPedSamps)
                AvgIntegral = RawIntegral - PedSum * (static_cast<Double_t>(WindowWidth) / NPedSamps);
                PulseFound = true;
                fWaveformInfo.PedSubIntegral.push_back(AvgIntegral);

                Int_t  pulse_amp = 0;
                Int_t pulse_Time = 0;

                Int_t PeakBin = 0;
                Double_t PeakVal = wave[i] - AvgPed;
                for( Int_t nt = i + 1; nt < TMath::Min((i + NSA), int(NSamples)); nt++ ) {
                if( (wave[nt] - AvgPed) < PeakVal && PeakBin == 0 ) {
                    PeakBin = nt - 1;
                } else if( PeakBin == 0 ) {
                    PeakVal = wave[nt] - AvgPed;
                }
                }
                if( PeakBin > 0 ) {
                pulse_amp = wave[PeakBin];
                Int_t Time = i * 64;
                Double_t VMid = (wave[PeakBin] - AvgPed) / 2.;
                for( Int_t nt = TMath::Max(i - fConfig.fNSB, 0);
                    nt < TMath::Min(PeakBin, int(NSamples - 1)); nt++ ) {
                    if( VMid >= (wave[nt] - AvgPed) && VMid < (wave[nt + 1] - AvgPed) ) {
                    Time = 64 * nt + int(64 * (VMid - (wave[nt] - AvgPed)) /
                        ((wave[nt + 1] - AvgPed) - (wave[nt] - AvgPed)));
                    }
                }
                pulse_Time = Time;
                } else {
                pulse_amp = wave[i];
                pulse_Time = 64 * i;
                }

                fWaveformInfo.PulseAmp.push_back(pulse_amp);
                fWaveformInfo.PulseTime.push_back(pulse_Time);
            }
        } else {
            PulseFound = false;
        }
    }
}

// Waveform generator v2 related functions - now the parameters are read from the object state
Double_t waveformGenerator_v2::spePulse(Double_t t, Double_t t0, Double_t amp) const {
    if(t < t0) return 0.0;
    return amp * (std::exp(-(t - t0) / fWf.tau_f) - std::exp(-(t - t0) / fWf.tau_r));
}

std::vector<Double_t> waveformGenerator_v2::generatePETimes() const {
    // Distributions for event and after pulse gen
    std::uniform_int_distribution<Int_t> num_pulses_dist(1, fDaq.max_n_pulses);
    std::uniform_real_distribution<Double_t> window_dist(10.0, fDaq.num_samples * fDaq.sample_rate - 50.0);
    std::poisson_distribution<Int_t> pe_per_hit_dist(15); // Avg 15 primary PEs per hit
    std::uniform_real_distribution<Double_t> ap_chance(0.0, 1.0); // Prob for after pulse
    std::normal_distribution<Double_t> ap_time_dist(fWf.ap_delay_mean, fWf.ap_delay_sigma);

    std::vector<Double_t> all_pe_times;
    Int_t num_pulses = num_pulses_dist(fRng);

    for (Int_t p = 0; p < num_pulses; p++) {
        Double_t hit_time = window_dist(fRng);
        Int_t num_pes = pe_per_hit_dist(fRng);

        for (Int_t i = 0; i < num_pes; i++) {
            all_pe_times.push_back(hit_time);

            // After pulse probability roll
            if(ap_chance(fRng) < fWf.ap_prob) {
                Double_t delayed_time = hit_time + ap_time_dist(fRng);
                all_pe_times.push_back(delayed_time);
            }
        }
    }
    // Sort the PE times
    std::sort(all_pe_times.begin(), all_pe_times.end());
    return all_pe_times;
}

std::vector<Int_t> waveformGenerator_v2::generateWaveform() const {
    std::vector<Double_t> pe_times = generatePETimes();

    // Distributions
    std::normal_distribution<Double_t> tts_dist(0.0, fWf.tts_sigma);
    std::normal_distribution<Double_t> gain_dist(fWf.mean_gain, fWf.gain_sigma);
    std::normal_distribution<Double_t> noise_dist(0.0, fDaq.noise_sigma);

    std::vector<Double_t> analog_wf(fDaq.num_samples, 0.0);

    for(Double_t true_times : pe_times) {
        Double_t smeared_time = true_times + tts_dist(fRng);
        Double_t smeared_gain = gain_dist(fRng);
        if(smeared_gain < 0.0) smeared_gain = 0.0;

        for(Int_t i = 0; i < fDaq.num_samples; i++) {
            Double_t current_time = i * fDaq.sample_rate;
            // Optimization
            if(current_time > smeared_time - 5.0 && current_time < smeared_time + 50) {
                analog_wf[i] += spePulse(current_time, smeared_time, smeared_gain);
            }
        }
    }

    // Add noise and digitize
    std::vector<Int_t> adc_wf(fDaq.num_samples, 0);
    for(Int_t i = 0; i < fDaq.num_samples; i++) {
        Double_t noisy_voilage = analog_wf[i] + noise_dist(fRng);
        Int_t adc_count = static_cast<Int_t>(std::floor(noisy_voilage/fDaq.v_lsb) + fDaq.pedestal);
        adc_count = std::max(0, std::min(adc_count, fDaq.max_adc));
        adc_wf[i] = adc_count;
    }
    return adc_wf;
}