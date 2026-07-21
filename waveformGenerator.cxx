#include "waveformGenerator.h"
#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

waveformGenerator::waveformGenerator() : fRng(std::random_device{}()) {}

std::vector<Int_t> waveformGenerator::noiseWaveform() const {
    int min = fDaq.pedestal - fDaq.ped_sigma / 2.0;
    int max = fDaq.pedestal + fDaq.ped_sigma / 2.0 + 1;
    std::uniform_real_distribution<double> rndm(0.0, 1.0);
    std::vector<Int_t> wave;
    for( int i = 0; i < fDaq.num_samples; i++ ) {
    wave.push_back(min + (max - min) * rndm(fRng));
    }
    return wave;
}

std::vector<Int_t> waveformGenerator::singlePulseWaveform() const {
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

std::vector<Int_t> waveformGenerator::multiplePulseWaveform() const {
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

void waveformGenerator::calcPulseParameters(const std::vector<Int_t> wave) {
    FADCConfig fConfig = saAnalyzer.getConfig();
    Int_t NPedSamps = fConfig.fNPedestalSamples;
    Int_t NSB       = fConfig.fNSB;
    Int_t NSA       = fConfig.fNSA;
    Int_t NSAT      = fConfig.fNSAT;
    Int_t Threshold = fConfig.fSampThreshold;
    Int_t NSamples  = (Int_t)wave.size();
    Int_t NS = NPedSamps - 1;

    Int_t PedFromDaq = fDaq.pedestal;

    auto GetIntegral = [&](Int_t lo, Int_t hi) {
        Double_t sum = 0.0;
        lo = TMath::Max(lo, 0);
        hi = TMath::Min(hi, NSamples - 1);
        for( Int_t i = lo; i <= hi; i++ ) sum += wave[i];
        return sum;
    };

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
    
    if((PedSum < (PedFromDaq + 8) * 4 && PedSum > (PedFromDaq - 8) * 4) &&
        (PedSumLastFour < (PedFromDaq + 8) * 4 && PedSumLastFour > (PedFromDaq - 8) * 4)) {
            if(PedSumLastFour < PedSum) {
                AvgPed = PedSumLastFour / NPedSamps;
            } else {
                AvgPed = PedSum / NPedSamps;
            }
        }
    else AvgPed = PedFromDaq;
    
    //Bool_t PulseFound = false;
    Int_t NSampPulses = 0;
    Bool_t CheckSampBelowThres = false;

    std::vector<Int_t> AvgIntegrals = {};

    while(NS < int(NSamples) && NSampPulses < fConfig.fMaxNPulses) {
      if(CheckSampBelowThres) {
        if( wave[NS] - AvgPed < Threshold ) CheckSampBelowThres = kFALSE;
      } else {
        Int_t ns_found = 0;
        for( Int_t nt = NS; nt < TMath::Min(NS + NSAT, int(NSamples)); nt++) {
          if( wave[nt] - AvgPed > Threshold) ns_found++;
        }
        if(ns_found == fConfig.fNSAT) {
            Int_t lo = TMath::Max(NS - fConfig.fNSB, 0);
            Int_t hi = TMath::Min(NS + fConfig.fNSA - 1, int(NSamples - 1));
            Int_t rawInt = GetIntegral(lo, hi);
            Int_t WindowWidth = hi - lo + 1;
            fWaveformInfo.PedSubIntegral.push_back(rawInt - AvgPed * WindowWidth);
          
          // Peak and time calculation
          Int_t PeakBin = 0;
          Double_t PeakVal = wave[NS] - AvgPed;
          for( Int_t nt = NS + 1; nt < TMath::Min((NS + fConfig.fNSA), int(NSamples)); nt++ ) {
            if( wave[nt] - AvgPed < PeakVal && PeakBin == 0 ) {
              PeakBin = nt - 1;
            } else if( PeakBin == 0 ) {
              PeakVal = wave[nt] - AvgPed;
            }
          }
          if( PeakBin > 0 ) {
            fWaveformInfo.PulseAmp.push_back(wave[PeakBin]);
            Int_t Time = NS * 64;
            Double_t VMid = (wave[PeakBin] - AvgPed) / 2.;
            for( Int_t nt = TMath::Max(NS - fConfig.fNSB, 0);
                nt < TMath::Min(PeakBin, int(NSamples - 1)); nt++ ) {
              if( VMid >= wave[nt] - AvgPed && VMid < wave[nt + 1] - AvgPed ) {
                Time = 64 * nt + int(64 * (VMid - (wave[nt] - AvgPed)) /
                    (wave[nt + 1] - AvgPed - (wave[nt] - AvgPed)));
              }
            }
            fWaveformInfo.PulseTime.push_back(Time);
          } else {
            fWaveformInfo.PulseAmp.push_back(wave[NS]);
            fWaveformInfo.PulseTime.push_back(64 * NS);
          }

          // Pileup calculations
        //   Int_t hi = TMath::Min((NS + fConfig.fNSA - 1), int(NSamples - 1));
        //   Int_t tot = 0;
        //   for( Int_t nt = NS; nt <= hi; nt++ ) {
        //     if( wave[nt] - AvgPed > fConfig.fSampThreshold ) tot++;
        //   }
        //   fWaveformInfo.PulseTOT.push_back(64.0 * tot);

        //   Int_t pileup = 0;
        //   for( Int_t nt = NS + 1; nt < hi; nt++ ) {
        //     if( wave[nt] - AvgPed > fConfig.fSampThreshold &&
        //         wave[nt] - AvgPed > wave[nt - 1] - AvgPed &&
        //         wave[nt] - AvgPed > wave[nt + 1] - AvgPed ) {
        //       pileup = 1;
        //       break;
        //     }
        //   }
        //   fWaveformInfo.PulsePileup.push_back(pileup);

        NSampPulses++;
        NS = NS + fConfig.fNSA;
        CheckSampBelowThres = kTRUE;
        }
      }
      NS++;
    }
}

void waveformGenerator::simpleParameterCalc(const std::vector<Int_t> wave) {
    // No need of this for now
}

// Waveform generator v2 related functions - now the parameters are read from the object state
Double_t waveformGenerator::spePulse(Double_t t, Double_t t0, Double_t amp) const {
    if(t < t0) return 0.0;
    return amp * (std::exp(-(t - t0) / fWf.tau_f) - std::exp(-(t - t0) / fWf.tau_r));
}

Double_t waveformGenerator::squarePulse(Double_t t, Double_t t0, Double_t amp) const {
    Double_t width_ns = fWf.pulse_width * fDaq.sample_rate; // Convert samples -> ns
    if(t < t0 || t > t0 + width_ns) return 0.0;
    return amp;
}

Double_t waveformGenerator::trianglePulse(Double_t t, Double_t t0, Double_t amp) const {
    Double_t width_ns = fWf.pulse_width * fDaq.sample_rate;
    Int_t dt = fDaq.sample_rate;
    t0 = round(t0 / dt) * dt;
    Double_t half_w = width_ns / 2.0;
    if(t < t0 - half_w || t > t0 + half_w) return 0.0;

    if(t < t0) {
        return amp * (t - (t0 - half_w)) / (half_w);
    } else {
        return amp * ((t0 + half_w) - t) / (half_w);
    }
}

std::vector<Double_t> waveformGenerator::generatePETimes() const {
    // Distributions for event and after pulse gen
    std::uniform_int_distribution<Int_t> num_pulses_dist(1, fDaq.max_n_pulses);
    std::uniform_real_distribution<Double_t> window_dist(10.0, fDaq.num_samples * fDaq.sample_rate - 50.0);
    std::poisson_distribution<Int_t> pe_per_hit_dist(15); // Avg 15 primary PEs per hit
    std::uniform_real_distribution<Double_t> ap_chance(0.0, 1.0); // Prob for after pulse
    std::normal_distribution<Double_t> ap_time_dist(fWf.ap_delay_mean, fWf.ap_delay_sigma); // After pulse time dist

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

std::vector<Int_t> waveformGenerator::generateWaveform() const {
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

std::vector<Int_t> waveformGenerator::squareWaveform() const {
    std::vector<Double_t> pe_times = generatePETimes();

    std::vector<Double_t> square_wf(fDaq.num_samples, 0.0);

    for(Double_t true_times : pe_times) {
        for(Int_t i = 0; i < fDaq.num_samples; i++) {
            Double_t current_time = i * fDaq.sample_rate;
            // Optimization
            if(current_time > true_times - 5.0 && current_time < true_times + 50) {
                square_wf[i] += squarePulse(current_time, true_times, fWf.mean_gain);
            }
        }
    }

    // Digitize
    std::vector<Int_t> adc_wf(fDaq.num_samples, 0);
    for(Int_t i = 0; i < fDaq.num_samples; i++) {
        Int_t adc_count = static_cast<Int_t>(std::floor(square_wf[i]/fDaq.v_lsb) + fDaq.pedestal);
        adc_count = std::max(0, std::min(adc_count, fDaq.max_adc));
        adc_wf[i] = adc_count;
    }
    return adc_wf;
}

std::vector<Int_t> waveformGenerator::triangleWaveform() const {
    std::vector<Double_t> pe_times = generatePETimes();

    std::vector<Double_t> square_wf(fDaq.num_samples, 0.0);

    for(Double_t true_times : pe_times) {
        for(Int_t i = 0; i < fDaq.num_samples; i++) {
            Double_t current_time = i * fDaq.sample_rate;
            // Optimization
            if(current_time > true_times - 5.0 && current_time < true_times + 50) {
                square_wf[i] += trianglePulse(current_time, true_times, fWf.mean_gain);
            }
        }
    }

    // Digitize
    std::vector<Int_t> adc_wf(fDaq.num_samples, 0);
    for(Int_t i = 0; i < fDaq.num_samples; i++) {
        Int_t adc_count = static_cast<Int_t>(std::floor(square_wf[i]/fDaq.v_lsb) + fDaq.pedestal);
        adc_count = std::max(0, std::min(adc_count, fDaq.max_adc));
        adc_wf[i] = adc_count;
    }
    return adc_wf;
}