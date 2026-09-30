#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include "waveformGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

// ROOT Headers
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>

TEST_CASE("Generate Waveform Fit Plot D2 v2", "[Plotting_v2][D2_v2]") {
  waveformGenerator gen;
  DigitizerParams daq;
  daq.noise_sigma = 1.5;
  daq.num_samples = 74;
  daq.max_n_pulses = 1;

  WaveformParams wf;
  wf.sigma_g = 1.5;
  wf.tau = 50;
  wf.mean_gain = 10.0;

  // Secondary pulse (pileup simulation)
  wf.sp_enable = false;
  wf.sp_delay = 25.0;
  wf.sp_amplitude_frac = 0.5;

  gen.setDigitizerParams(daq);
  gen.setDetectorParams(wf);

  FADCStandaloneAnalyzer analyzer;
  FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
  analyzer.setConfig(config);

  TemplateAnalysisConfig tmplCfg;
  // D2 specific configurations
  tmplCfg.fSmoothHalfWidth = 1; // Gaussian pre-smoothing kernel half-width
  tmplCfg.fD2ThresholdNSigma =
      5.0; // Significance level for the auto-scaled d2 threshold
  tmplCfg.fD2MinPeakSep = 1; // Minimum separation between pulse candidates
  tmplCfg.fAmpCutADC =
      5.0; // lowered from 8.0 to keep broad low-amplitude pulses
  tmplCfg.fFitMargin =
      5; // increased from 3 for better baseline in long-tail fits
  tmplCfg.fMaxBaselineAbs = 8.0; // relaxed from 5.0 for long-tail baseline bias
  tmplCfg.fScanHalf = 2.0;

  SECTION("Multiple iterations") {
    const int num_wfs = 1000;
    int passed_count = 0;
    std::vector<int> failed_runs;

    // Integral comparison parameters derived from the EMG model & DAQ specs
    Double_t tauSamples = (wf.tau / daq.sample_rate);
    Double_t sigmaSamples = (wf.sigma_g / daq.sample_rate);
    Double_t nIntSamples = 3.0 * sigmaSamples + 5.0 * tauSamples;
    Double_t abs_noise_margin = 3.0 * daq.noise_sigma * std::sqrt(nIntSamples);

    Double_t time_margin = daq.sample_rate * 5.0 * daq.tdc_subsamples;

    // Create the ROOT File
    TString name;
    if (wf.sp_enable) {
      name = TString::Format(
          "../root_files/d2_v2_test1/"
          "test_waveform_fits_d2_v2_pileup_delay_%g_sigma_%g_tau_%g_set2.root",
          wf.sp_delay, wf.sigma_g, wf.tau);
    } else {
      name =
          TString::Format("../root_files/d2_v2_single_pulse/"
                          "test_waveform_fits_d2_v2_sigma_%g_tau_%g_set2.root",
                          wf.sigma_g, wf.tau);
    }
    TFile outfile(name, "RECREATE");

    TH1F h_int_diff("h_int_diff",
                    "Integral Difference (Fit - Truth);ADC counts;Entries", 120,
                    -150, 150);
    TH1F h_int_rel_diff(
        "h_int_rel_diff",
        "Relative Integral Difference (Fit - Truth)/Truth;Fraction;Entries",
        100, -0.5, 0.5);
    TH1F h_amp_diff("h_amp_diff",
                    "Amplitude Difference (Fit - Truth);ADC counts;Entries",
                    100, -50, 50);
    TH1F h_time_halfmax_diff(
        "h_time_halfmax_diff",
        "Time Difference (half-max CFD);1/64 sample units;Entries", 100, -50,
        50);

    struct SavedPlot {
      int index;
      bool passed;
      TMultiGraph *mg;
    };
    std::vector<SavedPlot> saved_plots;

    for (int i = 0; i < num_wfs; i++) {
      std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
      std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
      const TruthPulseInfo &wfInfo = gen.getTruthInfo();

      FADCPulseWF templateResult = analyzer.AnalyzeWithEMG(samples, tmplCfg);

      // ── Validate results ──────────────────────────────────────────────
      bool run_passed = true;
      std::string s1 = "";
      std::string s2 = "";
      std::string s3 = "";

      if (templateResult.nPulses < 1) {
        run_passed = false;
        s3 = "No_pulses_found";
      } else if (templateResult.pulses.size() < wfInfo.pulseIntegral.size()) {
        run_passed = false;
        s3 = Form("Missed_pulse_found_%zu_expected_%zu",
                  templateResult.pulses.size(), wfInfo.pulseIntegral.size());
      } else {
        std::vector<bool> matched(templateResult.pulses.size(), false);

        for (size_t k = 0; k < wfInfo.pulseIntegral.size(); k++) {
          Double_t truth_time_sample = wfInfo.pulseTime[k] / 64.0;

          // 1. Time-association matching within +/- 4 samples
          int best_idx = -1;
          Double_t min_dt = 1e9;
          for (size_t p = 0; p < templateResult.pulses.size(); p++) {
            if (matched[p])
              continue;
            Double_t dt = std::abs(templateResult.pulses[p].timeSample -
                                   truth_time_sample);
            if (dt < min_dt) {
              min_dt = dt;
              best_idx = p;
            }
          }

          if (best_idx == -1 || min_dt > 4.0) {
            run_passed = false;
            s3 = "Unmatched_pulse_" + std::to_string(k);
            break;
          }
          matched[best_idx] = true;

          const auto &fitPulse = templateResult.pulses[best_idx];

          // 2. Dynamic Integral Margin:
          // Primary pulse: 12% relative tolerance
          // Secondary pulse (sp_enable): 20% relative tolerance due to tail
          // pileup
          Double_t rel_tol = (k == 0) ? 0.12 : 0.20;

          // Boundary handling: If secondary pulse tail is truncated by the
          // finite window, waveformGenerator's Riemann sum cuts off early while
          // AnalyzeWithD2 integrates the model.
          bool is_near_boundary =
              (truth_time_sample + 3.0 * tauSamples > daq.num_samples);
          if (is_near_boundary) {
            rel_tol += 0.15;
          }

          Double_t dyn_int_margin =
              abs_noise_margin + rel_tol * wfInfo.pulseIntegral[k];
          Double_t diff_int = fitPulse.integral - wfInfo.pulseIntegral[k];

          if (std::abs(diff_int) > dyn_int_margin) {
            run_passed = false;
            s1 = "Integral_mismatch_pulse_" + std::to_string(k);
          }

          // 3. Dynamic Amplitude Margin: 3*sigma noise floor + 10% relative
          Double_t truth_amp = wfInfo.pulseAmplitude[k] - daq.pedestal;
          Double_t dyn_amp_margin = 3.0 * daq.noise_sigma + 0.10 * truth_amp;
          Double_t diff_amp = fitPulse.amplitude - truth_amp;
          if (std::abs(diff_amp) > dyn_amp_margin) {
            run_passed = false;
            s2 = "Amp_mismatch_pulse_" + std::to_string(k);
          }

          // 4. Time verification
          Double_t time_diff =
              fitPulse.timeHalfMax * daq.tdc_subsamples - wfInfo.pulseTime[k];
          if (std::abs(time_diff) > time_margin) {
            s3 = "Time_mismatch_pulse_" + std::to_string(k);
          }

          // Fill diagnostic histograms
          h_int_diff.Fill(diff_int);
          if (wfInfo.pulseIntegral[k] > 0.0) {
            h_int_rel_diff.Fill(diff_int / wfInfo.pulseIntegral[k]);
          }
          h_amp_diff.Fill(diff_amp);
          h_time_halfmax_diff.Fill(time_diff);
        }
      }

      if (run_passed) {
        passed_count++;
      } else {
        failed_runs.push_back(i);
      }

      if (!run_passed || (i % 100 == 0)) {
        Int_t N = static_cast<Int_t>(samples.size());

        TGraph *gSamples = new TGraph(N);
        for (Int_t j = 0; j < N; ++j) {
          gSamples->SetPoint(j, (double)j, samples[j]);
        }
        gSamples->SetMarkerStyle(3);
        gSamples->SetMarkerSize(0.5);
        gSamples->SetMarkerColor(kBlue);

        TF1 fTotal(
            "fTotal",
            [templateResult](double *x, double * /*p*/) {
              double t = x[0];
              double val = templateResult.pedestal;
              for (const auto &pulse : templateResult.pulses) {
                if (pulse.valid) {
                  // D2 path always uses the fitted EMG model (Stage 2)
                  double dt = t - pulse.timeSample;
                  double sigma = pulse.fittedSigma;
                  double tau = pulse.fittedTau;
                  if (sigma > 0 && tau > 0) {
                    double s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
                    double erf_arg =
                        (dt - sigma * sigma / tau) / (sigma * std::sqrt(2.0));
                    double raw =
                        std::exp(s2_2t2 - dt / tau) * (1.0 + std::erf(erf_arg));
                    // Normalize so peak = 1, then scale by fitted amplitude
                    double dt_peak = sigma * sigma / tau;
                    double raw_peak = std::exp(s2_2t2 - dt_peak / tau) *
                                      (1.0 + std::erf(0.0));
                    if (raw_peak > 1e-12)
                      val += pulse.amplitude * raw / raw_peak;
                  }
                }
              }
              return val;
            },
            -0.5, N - 0.5, 0);

        // Build fit graph at 0.25-sample steps
        const double fitStep = 0.25;
        const double fitXMin = -0.5;
        const double fitXMax = N - 0.5;
        int nFitPts = static_cast<int>((fitXMax - fitXMin) / fitStep) + 1;
        TGraph *gFit = new TGraph(nFitPts);
        for (int k = 0; k < nFitPts; ++k) {
          double xv = fitXMin + k * fitStep;
          gFit->SetPoint(k, xv, fTotal.Eval(xv));
        }
        gFit->SetLineColor(kRed);
        gFit->SetLineWidth(1);

        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gSamples, "P");
        mg->Add(gFit, "L");

        if (!run_passed) {
          if (s1.empty() && s2.empty() && s3.empty()) {
            mg->SetName(Form("Waveform_Run_%d", i));
            mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
          } else {
            mg->SetName(Form("Waveform_Run_%d_%s_%s_%s", i, s1.c_str(),
                             s2.c_str(), s3.c_str()));
            mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
          }
        } else {
          mg->SetName(Form("Waveform_Run_%d", i));
          mg->SetTitle(Form("Template Fit Waveform (Run %d);Sample;ADC", i));
        }

        saved_plots.push_back({i, run_passed, mg});
      }
    }

    double pass_probability =
        (static_cast<double>(passed_count) / num_wfs) * 100.0;

    for (const auto &plot : saved_plots) {
      if (!plot.passed) {
        plot.mg->Write();
      } else if (pass_probability > 98.0 && (plot.index % 100 == 0)) {
        plot.mg->Write();
      }
    }

    h_int_diff.Write();
    h_int_rel_diff.Write();
    h_amp_diff.Write();
    h_time_halfmax_diff.Write();

    outfile.Close();

    std::cout << "\n======================================================\n";
    std::cout << "        TEMPLATE WAVEFORM D2 V2 TEST STATISTICS       \n";
    std::cout << "======================================================\n";
    std::cout << "Total Runs Executed : " << num_wfs << "\n";
    std::cout << "Passed Runs         : " << passed_count << "\n";
    std::cout << "Failed Runs         : " << failed_runs.size() << "\n";
    std::cout << "Pass Probability    : " << std::fixed << std::setprecision(2)
              << pass_probability << "%\n";
    std::cout << "Waveforms Saved To  : " << name.Data() << "\n\n";

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
