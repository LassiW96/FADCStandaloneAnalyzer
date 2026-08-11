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
    daq.noise_sigma = 0.5;
    daq.num_samples = 120;

    WaveformParams wf;
    wf.tau_r = 2;
    wf.tau_f = 7;
    wf.mean_gain = 10.0;

    gen.setDigitizerParams(daq);
    gen.setDetectorParams(wf);

    FADCStandaloneAnalyzer analyzer;
    FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
    analyzer.setConfig(config);

    TemplateAnalysisConfig tmplCfg;
    tmplCfg.fMFSeparation = 8;
    tmplCfg.fScanHalf = 6;     // ±6 from corrected candidateBin = peak - peakBin
    tmplCfg.fMFLeft = 8;
    tmplCfg.fMFRight = 8;
    tmplCfg.fMFStart = 0;     // search from sample 0 (MF loop now starts at 0)
    tmplCfg.fAmpCutADC = 4.0;

    SPEParams sp;
    sp.sigma = 1;
    sp.tau = 8;

    FADCTemplate tmpl = FADCTemplate::fromWaveformParams(sp);
    
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
        TH1F h_time_halfmax_diff("h_time_halfmax_diff", "Time Difference (half-max CFD)", 100, -100, 100);
        
        struct SavedPlot {
            int index;
            bool passed;
            TMultiGraph* mg;
        };
        std::vector<SavedPlot> saved_plots;

        for(int i = 0; i < num_wfs; i++){
            std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
            std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
            FADCPulse standardResult = analyzer.Analyze(samples, 1.0, 0.0);
            const TruthPulseInfo& wfInfo = gen.getTruthInfo();

            FADCPulseWF templateResult = analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

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

                    if(std::abs(templateResult.pulses[k].timeHalfMax * daq.tdc_subsamples - (wfInfo.pulseTime[k])) > time_margin) {
                        // run_passed = false;
                        s3 = "Time mismatch at pulse " + std::to_string(k);
                    }

                    // Fill hists
                    h_int_diff.Fill(templateResult.pulses[k].integral - wfInfo.pulseIntegral[k]);
                    h_amp_diff.Fill(templateResult.pulses[k].amplitude - (wfInfo.pulseAmplitude[k] - daq.pedestal));
                    // Half-max CFD time vs. truth (same unit conversion)
                    h_time_halfmax_diff.Fill(templateResult.pulses[k].timeHalfMax * daq.tdc_subsamples - (wfInfo.pulseTime[k]));
                }
            }

            if(run_passed) {
                passed_count++;
            }
            else {
                failed_runs.push_back(i);
            }

            if(!run_passed || (i % 100 == 0)) {
                Int_t N = static_cast<Int_t>(samples.size());

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

                // Build fit graph at 0.25-sample steps
                const double fitStep = 0.25;
                const double fitXMin = -0.5;
                const double fitXMax = N - 0.5;
                int nFitPts = static_cast<int>((fitXMax - fitXMin) / fitStep) + 1;
                TGraph* gFit = new TGraph(nFitPts);
                for (int k = 0; k < nFitPts; ++k) {
                    double xv = fitXMin + k * fitStep;
                    gFit->SetPoint(k, xv, fTotal.Eval(xv));
                }
                gFit->SetLineColor(kRed);
                gFit->SetLineWidth(1);

                TMultiGraph* mg = new TMultiGraph();
                mg->Add(gSamples, "P"); 
                mg->Add(gFit, "L");     

                if(!run_passed) {
                    if(s1.empty() && s2.empty() && s3.empty()) {
                        mg->SetName(Form("Waveform_Run_%d", i));
                        mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
                    }
                    else {
                        mg->SetName(Form("Waveform_Run_%d_%s_%s_%s", i, s1.c_str(), s2.c_str(), s3.c_str()));
                        mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
                    }
                }
                else {
                    mg->SetName(Form("Waveform_Run_%d", i));
                    mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
                }

                saved_plots.push_back({i, run_passed, mg});
            }
        }

        double pass_probability = (static_cast<double>(passed_count) / num_wfs) * 100.0;

        for (const auto& plot : saved_plots) {
            if (!plot.passed) {
                plot.mg->Write();
            }
            else if (pass_probability > 99.0 && (plot.index % 100 == 0)) {
                plot.mg->Write();
            }
        }

        h_int_diff.SetTitle("Integral Difference");
        h_int_diff.Write();

        h_amp_diff.SetTitle("Amplitude Difference");
        h_amp_diff.Write();

        h_time_halfmax_diff.SetTitle("Time Difference (half-max CFD)");
        h_time_halfmax_diff.Write();

        outfile.Close();

        std::cout << "\n======================================================\n";
        std::cout << "           TEMPLATE WAVEFORM TEST STATISTICS          \n";
        std::cout << "======================================================\n";
        std::cout << "Total Runs Executed : " << num_wfs << "\n";
        std::cout << "Passed Runs         : " << passed_count << "\n";
        std::cout << "Failed Runs         : " << failed_runs.size() << "\n";
        std::cout << "Pass Probability    : " << std::fixed << std::setprecision(2) << pass_probability << "%\n";
        std::cout << "Waveforms Saved To  : test_waveform_fits2.root\n\n";

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