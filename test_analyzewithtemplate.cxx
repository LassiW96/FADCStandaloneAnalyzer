#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include "waveformGenerator.h"

// ROOT Headers
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <TMultiGraph.h>

TEST_CASE("Generate Waveform Fit Plot", "[Plotting]") {
    waveformGenerator gen;
    DigitizerParams daq;
    daq.noise_sigma = 1.0;
    daq.num_samples = 120;

    WaveformParams wf;
    wf.tau_r = 2;
    wf.tau_f = 7;
    wf.mean_gain = 3.0;

    gen.setDigitizerParams(daq);
    gen.setDetectorParams(wf);

    FADCStandaloneAnalyzer analyzer;
    FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
    analyzer.setConfig(config);

    TemplateAnalysisConfig tmplCfg;
    tmplCfg.fMFSeparation = 8;
    tmplCfg.fScanHalf = 10;
    tmplCfg.fMFLeft = 8;
    tmplCfg.fMFRight = 8;
    tmplCfg.fAmpCutADC = 4.0;
    

    FADCTemplate tmpl = FADCTemplate::fromWaveformParams(WaveformParams{});
    
    SECTION("Multiple iterations"){
        const int num_wfs = 1000;
        int passed_count = 0;
        std::vector<int> failed_runs;
        
        Double_t int_margin = wf.mean_gain * wf.pulse_width * 5.0;
        Double_t amp_margin = wf.mean_gain * 5.0;
        Double_t time_margin = daq.sample_rate * 10.0;
        
        // Create the ROOT File
        TFile outfile("test_waveform_fits2.root", "RECREATE");

        TH1F h_int_diff("h_int_diff", "Integral Difference", 100, -100, 100);
        TH1F h_amp_diff("h_amp_diff", "Amplitude Difference", 100, -100, 100);
        TH1F h_time_diff("h_time_diff", "Time Difference", 100, -100, 100);
        
        for(int i = 0; i < num_wfs; i++){
            std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
            std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
            FADCPulse standardResult = analyzer.Analyze(samples, 1.0, 0.0);
            const TruthPulseInfo& wfInfo = gen.getTruthInfo();

            FADCPulseWF templateResult = analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

            Int_t N = static_cast<Int_t>(samples.size());
            // TH1D hWaveform("hWaveform" + TString::Itoa(i, 10), 
            //     "Template Fit (k = " + TString::Itoa(i, 10) + ");Sample;ADC", N, -0.5, N - 0.5);
            // for (Int_t j = 0; j < N; ++j) {
            //     hWaveform.SetBinContent(j + 1, samples[j]);
            // }
            // hWaveform.SetStats(0);

            TGraph* gSamples = new TGraph(N);
            for (Int_t j = 0; j < N; ++j) {
                gSamples->SetPoint(j, (double)j, samples[j]);
            }
            gSamples->SetMarkerStyle(3);
            gSamples->SetMarkerSize(0.5);
            gSamples->SetMarkerColor(kBlue);

            TF1 fTotal("fTotal", 
                [tmpl, templateResult](double* x, double* /*p*/) {
                    double t = x[0];
                    double val = templateResult.pedestal;
                    for (const auto& pulse : templateResult.pulses) {
                        if (pulse.valid) {
                            val += pulse.amplitude * tmpl.eval(t - pulse.timeSample);
                        }
                    }
                    return val;
                }, -0.5, N - 0.5, 0);

            fTotal.SetNpx(1000);

            TGraph* gFit = new TGraph(&fTotal);
            gFit->SetLineColor(kRed);
            gFit->SetLineWidth(1);

            // hWaveform.GetListOfFunctions()->Add(gFit, "L");

            TMultiGraph mg;
            mg.Add(gSamples, "P"); 
            mg.Add(gFit, "L");     

            // ── Validate results ──────────────────────────────────────────────
            bool run_passed = true;

            std::string s1 = "";
            std::string s2 = "";
            std::string s3 = "";

            if(templateResult.nPulses < 1) {
                run_passed = false;
            }

            if(templateResult.pulses.size() < wfInfo.pulseIntegral.size()) {
                run_passed = false;
            }
            else {
                for(size_t k = 0; k < wfInfo.pulseIntegral.size(); k++) {
                    if(std::abs(templateResult.pulses[k].integral - wfInfo.pulseIntegral[k]) > int_margin) {
                        run_passed = false;
                        s1 = "Integral mismatch at pulse " + std::to_string(k);
                    }

                    if(std::abs(templateResult.pulses[k].amplitude - (wfInfo.pulseAmplitude[k] - daq.pedestal)) > amp_margin) {
                        run_passed = false;
                        s2 = "Amplitude mismatch at pulse " + std::to_string(k);
                    }

                    if(std::abs(templateResult.pulses[k].timeSample * daq.tdc_subsamples - (wfInfo.pulseTime[k] - daq.tdc_subsamples*5)) > time_margin) {
                        // run_passed = false;
                        s3 = "Time mismatch at pulse " + std::to_string(k);
                    }

                    // Fill hists
                    h_int_diff.Fill(templateResult.pulses[k].integral - wfInfo.pulseIntegral[k]);
                    h_amp_diff.Fill(templateResult.pulses[k].amplitude - (wfInfo.pulseAmplitude[k] - daq.pedestal));
                    h_time_diff.Fill(templateResult.pulses[k].timeSample * daq.tdc_subsamples - (wfInfo.pulseTime[k] - daq.tdc_subsamples*5));
                }
            }

            if(run_passed) {
                passed_count++;
            }
            else {
                failed_runs.push_back(i);

                if(s1.empty() && s2.empty() && s3.empty()) {
                    mg.SetName(Form("Waveform_Run_%d", i));
                    mg.SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
                }
                else {
                    mg.SetName(Form("Waveform_Run_%d_%s_%s_%s", i, s1.c_str(), s2.c_str(), s3.c_str()));
                    mg.SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
                }
                // hWaveform.Write();
                mg.Write();
            }
        }

        h_int_diff.SetTitle("Integral Difference");
        h_int_diff.Write();

        h_amp_diff.SetTitle("Amplitude Difference");
        h_amp_diff.Write();

        h_time_diff.SetTitle("Time Difference");
        h_time_diff.Write();

        outfile.Close();

        double pass_probability = (static_cast<double>(passed_count) / num_wfs) * 100.0;

        std::cout << "\n======================================================\n";
        std::cout << "           TEMPLATE WAVEFORM TEST STATISTICS          \n";
        std::cout << "======================================================\n";
        std::cout << "Total Runs Executed : " << num_wfs << "\n";
        std::cout << "Passed Runs         : " << passed_count << "\n";
        std::cout << "Failed Runs         : " << failed_runs.size() << "\n";
        std::cout << "Pass Probability    : " << std::fixed << std::setprecision(2) << pass_probability << "%\n";
        std::cout << "Waveforms Saved To  : test_waveform_fits.root\n\n";

        if (!failed_runs.empty()) {
            std::cout << "Failed Run Indices  : ";
            for (size_t i = 0; i < failed_runs.size(); ++i) {
                std::cout << failed_runs[i] << (i < failed_runs.size() - 1 ? ", " : "");
            }
            std::cout << "\n";
        }
        std::cout << "======================================================\n\n";

        REQUIRE(pass_probability >= 75.0);
    }
}