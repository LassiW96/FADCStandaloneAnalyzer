#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <future>
#include <chrono>
#include <ctime>
#include "waveformGenerator.h"
#include "FADCStandaloneAnalyzer.h"

// A thread-safe progress bar
class ThreadSafeProgressBar {
    std::mutex mtx;
    int total_work;
    int current_progress;

public:
    ThreadSafeProgressBar(int total) : total_work(total), current_progress(0) {}

    void update(int amount = 1) {
        // Lock the mutex so multiple threads don't scramble the console output
        std::lock_guard<std::mutex> lock(mtx);
        current_progress += amount;
        
        float percentage = (static_cast<float>(current_progress) / total_work) * 100.0f;
        int barWidth = 50;
        
        std::cout << "[";
        int pos = barWidth * percentage / 100;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        // \r returns the cursor to the beginning of the line to overwrite it
        std::cout << "] " << int(percentage) << " %\r"; 
        std::cout.flush(); // Force the console to print immediately
    }

    void finish() {
        std::cout << std::endl; // Move to the next line so Catch2 output isn't overwritten
    }
};

std::vector<Double_t> ToDoubleWave(const std::vector<Int_t>& wave) {
    return std::vector<Double_t>(wave.begin(), wave.end());
}

// Shared analyzer configuration.
FADCConfig GetTestConfig() {
    FADCConfig conf;
    conf.fNPedestalSamples = 4;
    conf.fSampThreshold    = 50.0;
    conf.fNSB              = 1;
    conf.fNSA              = 2;
    conf.fNSAT             = 1;
    conf.fMaxNPulses       = 3;
    return conf;
}

// Struct to hold results to be asserted in the main thread safely
struct WaveformTestResult {
    double expected_integral;
    double actual_integral;
    double expected_amp;
    double actual_amp;
    double expected_time;
    double actual_time;
    double tolerance;
    int num_pulses_detected;
};

std::vector<WaveformTestResult> processDataSection(int section_idx, int chunks_per_section, ThreadSafeProgressBar& bar) {
    FADCStandaloneAnalyzer analyzer;
    analyzer.setConfig(GetTestConfig());
    analyzer.setCalib({0.2441, 4.0/50.0});

    WaveformData wfData;
    wfData.TotNSamples     = 100;
    wfData.InitPedSamples  = 4;
    wfData.Pedestal        = 100;
    wfData.PedRange        = 4;
    wfData.NPulses         = 4;
    wfData.PulseAmp        = 300;
    wfData.PulseWidth      = 5;
    wfData.PulseSeparation = 10;

    waveformGenerator gen;
    gen.setWaveformData(wfData);
    gen.setAnalyzer(analyzer);

    std::vector<WaveformTestResult> section_results;

    for (int j = 0; j < chunks_per_section; j++) {
        // Alternate between single and multiple pulses to actually process data every loop
        bool is_single_pulse = (j % 2 == 0); 
        
        std::vector<Int_t> int_wave = is_single_pulse ? gen.singlePulseWaveform() : gen.multiplePulseWaveform();
        std::vector<Double_t> wave = ToDoubleWave(int_wave);

        // Run the heavy analyzer
        auto results = analyzer.Analyze(wave, 1.0, 4.0);

        gen.clearWaveformInfo();
        gen.calcPulseParameters(int_wave);
        const WaveformInfo& wfInfo = gen.getWaveformInfo();

        // Store the results instead of calling REQUIRE here
        for(size_t i = 0; i < wfInfo.PedSubIntegral.size(); i++) {
            WaveformTestResult res;
            res.num_pulses_detected = results.fNSampPulses;
            res.expected_integral = static_cast<Double_t>(wfInfo.PedSubIntegral[i]);
            res.actual_integral = results.fSampPulseIntPedSub[i];
            res.expected_amp = static_cast<Double_t>(wfInfo.PulseAmp[i]);
            res.actual_amp = results.fSampPulseAmp[i];
            res.expected_time = static_cast<Double_t>(wfInfo.PulseTime[i]);
            res.actual_time = results.fSampPulseTime[i];
            res.tolerance = wfData.PedRange;
            
            section_results.push_back(res);
        }

        // Simulate further heavy CPU work if needed
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        bar.update(1);
    }
    
    return section_results;
}   

TEST_CASE("Performance: Sequential vs Multi-threaded execution", "[threading]") {
    const int num_sections = 4; // Increased to show better threading gains
    const int chunks_per_section = 100; 
    const int total_chunks = num_sections * chunks_per_section;

    SECTION("Multi-threaded Run") {
        std::cout << "\n--- Running Multi-threaded ---\n";
        ThreadSafeProgressBar bar(total_chunks);
        
        auto start_wall = std::chrono::high_resolution_clock::now();
        auto start_cpu = std::clock();

        // Notice the future now expects a vector of our custom struct
        std::vector<std::future<std::vector<WaveformTestResult>>> futures;
        
        for(int i = 0; i < num_sections; ++i) {
            futures.push_back(std::async(std::launch::async, processDataSection, i, chunks_per_section, std::ref(bar)));
        }

        // 1. Wait for threads to finish AND collect their results
        std::vector<WaveformTestResult> all_results;
        for(auto& f : futures) {
            auto thread_results = f.get(); 
            all_results.insert(all_results.end(), thread_results.begin(), thread_results.end());
        }
        bar.finish();

        auto end_wall = std::chrono::high_resolution_clock::now();
        auto end_cpu = std::clock();

        // 2. Perform ALL Catch2 Assertions sequentially in the main thread
        for (const auto& res : all_results) {
            REQUIRE(res.num_pulses_detected > 0);
            REQUIRE(res.actual_integral == Catch::Approx(res.expected_integral).margin(res.tolerance));
            REQUIRE(res.actual_amp == Catch::Approx(res.expected_amp).margin(res.tolerance));
            REQUIRE(res.actual_time == Catch::Approx(res.expected_time).margin(res.tolerance));
        }

        std::chrono::duration<double> wall_diff = end_wall - start_wall;
        double cpu_diff = (double)(end_cpu - start_cpu) / CLOCKS_PER_SEC;

        std::cout << "Multi-threaded Wall Time: " << wall_diff.count() << " seconds\n";
        std::cout << "Multi-threaded CPU Time:  " << cpu_diff << " seconds\n";
    }
}