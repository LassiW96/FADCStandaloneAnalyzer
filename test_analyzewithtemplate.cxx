#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include "waveformGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <iomanip>
#include <iostream>

// ROOT Headers
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>

TEST_CASE("Generate Waveform Fit Plot", "[Plotting]") {
  waveformGenerator gen;
  DigitizerParams daq;
  daq.noise_sigma = 1.5;
  daq.num_samples = 74;
  daq.max_n_pulses = 1;

  WaveformParams wf;
  wf.sigma_g = 2.6;
  wf.tau = 7;
  wf.mean_gain = 10.0;

  // Secondary pulse
  wf.sp_enable = false;
  wf.sp_delay = 65.0;
  wf.sp_amplitude_frac = 0.5;

  gen.setDigitizerParams(daq);
  gen.setDetectorParams(wf);

  // ── Generate a clean (noiseless) template waveform ────────────────────────
  // Build the template analytically from WaveformUtil::spePulse using the same
  // sigma_g and tau as the detector model.  The pulse is placed at a known t0
  // (20 ns, after the 4-sample pedestal window) so the template baseline is 0
  // and peakBin is deterministic.
  const Int_t tmplNSamples = 50;   // compact template window
  const Double_t tmplT0_ns = 20.0; // pulse onset at 20 ns (after 16 ns pedestal window)
  const Int_t tmplNPed = daq.init_ped_samples;

  // Evaluate the EMG pulse at each sample; find the raw peak for normalisation
  std::vector<Double_t> tmplAnalog(tmplNSamples, 0.0);
  Double_t rawPeak = 0.0;
  for (Int_t i = 0; i < tmplNSamples; ++i) {
    Double_t t = i * daq.sample_rate; // time in ns
    tmplAnalog[i] = WaveformUtil::spePulse(t, tmplT0_ns, wf.sigma_g, wf.tau);
    if (tmplAnalog[i] > rawPeak)
      rawPeak = tmplAnalog[i];
  }

  // Scale to the mean_gain amplitude and add the pedestal, then digitise
  std::vector<Double_t> cleanWaveform(tmplNSamples);
  for (Int_t i = 0; i < tmplNSamples; ++i) {
    Double_t amplitude_mV =
        (rawPeak > 1e-12) ? wf.mean_gain * tmplAnalog[i] / rawPeak : 0.0;
    cleanWaveform[i] = std::round(amplitude_mV / daq.v_lsb + daq.pedestal);
  }

  FADCTemplate tmpl =
      FADCTemplate::fromWaveformSamples(cleanWaveform, tmplNPed,
                                        /* mfLeft */ 5, /* mfRight */ 5);

  // Verify the template was built successfully
  REQUIRE(!tmpl.shape.empty());
  REQUIRE(tmpl.peakBin >= 0);

  FADCStandaloneAnalyzer analyzer;
  FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
  analyzer.setConfig(config);

  TemplateAnalysisConfig tmplCfg;
  // Template / matched-filter specific configurations
  tmplCfg.fMFLeft = 5;                  // MF kernel half-width left
  tmplCfg.fMFRight = 5;                 // MF kernel half-width right
  tmplCfg.fMFStart = 1;                 // first sample of MF search window (include early pulses)
  tmplCfg.fMFEnd = daq.num_samples - 1; // last sample of MF search window
  tmplCfg.fMFThreshold = 3.0 * daq.noise_sigma; // 3-sigma noise threshold
  tmplCfg.fMFSeparation = 3;            // minimum distance between MF peaks
  tmplCfg.fAmpCutADC = 8.0;
  tmplCfg.fFitMargin = 6;
  tmplCfg.fScanHalf = 15.0;
  tmplCfg.fRefitMaxIter = 3;

  SECTION("Multiple iterations") {
    const int num_wfs = 1000;
    int passed_count = 0;
    std::vector<int> failed_runs;

    Double_t int_margin = wf.mean_gain * 20.0;
    Double_t amp_margin = wf.mean_gain * 5.0;
    Double_t time_margin = daq.sample_rate * 5.0;

    // Create the ROOT File
    TString name;
    if (wf.sp_enable) {
      name = TString::Format(
          "../root_files/template_an_test1/"
          "test_waveform_fits_pileup_delay_%g_noise_%g_updated.root",
          wf.sp_delay, daq.noise_sigma);
    } else {
      name = TString::Format("../root_files/template_an_test1/"
                             "test_waveform_fits_noise_%g_updated.root",
                             daq.noise_sigma);
    }
    TFile outfile(name, "RECREATE");

    TH1F h_int_diff("h_int_diff", "Integral Difference", 100, -100, 100);
    TH1F h_amp_diff("h_amp_diff", "Amplitude Difference", 100, -50, 50);
    TH1F h_time_halfmax_diff("h_time_halfmax_diff",
                             "Time Difference (half-max CFD)", 100, -50, 50);

    struct SavedPlot {
      int index;
      bool passed;
      TMultiGraph *mg;
    };
    std::vector<SavedPlot> saved_plots;

    for (int i = 0; i < num_wfs; i++) {
      std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
      std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
      FADCPulse standardResult = analyzer.Analyze(samples, 1.0, 0.0);
      const TruthPulseInfo &wfInfo = gen.getTruthInfo();

      FADCPulseWF templateResult =
          analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

      // ── Validate results ──────────────────────────────────────────────
      bool run_passed = true;

      std::string s1 = "";
      std::string s2 = "";
      std::string s3 = "";

      if (templateResult.nPulses < 1) {
        run_passed = false;
      }

      if (templateResult.pulses.size() < wfInfo.pulseIntegral.size()) {
        run_passed = false;
      } else {
        for (size_t k = 0; k < wfInfo.pulseIntegral.size(); k++) {
          if (std::abs(templateResult.pulses[k].integral -
                       wfInfo.pulseIntegral[k]) > int_margin) {
            run_passed = false;
            s1 = "Integral mismatch at pulse " + std::to_string(k);
          }

          if (std::abs(templateResult.pulses[k].amplitude -
                       (wfInfo.pulseAmplitude[k] - daq.pedestal)) >
              amp_margin) {
            run_passed = false;
            s2 = "Amplitude mismatch at pulse " + std::to_string(k);
          }

          if (std::abs(templateResult.pulses[k].timeHalfMax *
                           daq.tdc_subsamples -
                       (wfInfo.pulseTime[k])) > time_margin) {
            // run_passed = false;
            s3 = "Time mismatch at pulse " + std::to_string(k);
          }

          // Fill hists
          h_int_diff.Fill(templateResult.pulses[k].integral -
                          wfInfo.pulseIntegral[k]);
          h_amp_diff.Fill(templateResult.pulses[k].amplitude -
                          (wfInfo.pulseAmplitude[k] - daq.pedestal));
          // Half-max CFD time vs. truth (same unit conversion)
          h_time_halfmax_diff.Fill(templateResult.pulses[k].timeHalfMax *
                                       daq.tdc_subsamples -
                                   (wfInfo.pulseTime[k]));
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

        // Build fit graph using stage-aware model evaluation
        auto evalPulseModel = [&](Double_t t, const FitPulseResult &pulse) -> Double_t {
          if (!pulse.valid) return 0.0;
          if (pulse.fitStage == 1 || pulse.fittedSigma == 0.0) {
            return pulse.amplitude * tmpl.eval(t - pulse.timeSample);
          } else {
            Double_t dt = t - pulse.timeSample;
            Double_t sigma = pulse.fittedSigma;
            Double_t tau = pulse.fittedTau;
            if (sigma <= 0.0 || tau <= 0.0) return 0.0;
            Double_t s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
            Double_t erf_arg = (dt - (sigma * sigma) / tau) / (sigma * std::sqrt(2.0));
            if (erf_arg > 700.0) return 0.0;
            Double_t f = std::exp(s2_2t2 - (dt / tau)) * (1.0 + std::erf(erf_arg));
            Double_t dt_peak = (sigma * sigma) / tau;
            Double_t f_peak = std::exp(s2_2t2 - (dt_peak / tau)) * (1.0 + std::erf(0.0));
            return (f_peak > 1e-12 && f > 1e-12) ? (pulse.amplitude * f / f_peak) : 0.0;
          }
        };

        const double fitStep = 0.25;
        const double fitXMin = -0.5;
        const double fitXMax = N - 0.5;
        int nFitPts = static_cast<int>((fitXMax - fitXMin) / fitStep) + 1;
        TGraph *gFit = new TGraph(nFitPts);
        for (int k = 0; k < nFitPts; ++k) {
          double xv = fitXMin + k * fitStep;
          double val = templateResult.pedestal;
          for (const auto &pulse : templateResult.pulses) {
            val += evalPulseModel(xv, pulse);
          }
          gFit->SetPoint(k, xv, val);
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
    std::cout << "Pass Probability    : " << std::fixed << std::setprecision(2)
              << pass_probability << "%\n";
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