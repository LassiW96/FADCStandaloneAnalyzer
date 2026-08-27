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
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

TEST_CASE("AnalyzeWithD2 Intermediate Steps", "[D2_Intermediate_Steps]") {
  // ── 0. Test Configuration ────────────────────────────────────────────────
  // Number of waveforms to generate and visualize (adjustable)
  const int NUM_WAVEFORMS = 5;

  waveformGenerator gen;
  DigitizerParams daq;
  daq.noise_sigma = 2.5;
  daq.num_samples = 74;
  daq.max_n_pulses = 1;

  WaveformParams wf;
  wf.sigma_g = 3.0;
  wf.tau = 12.0;
  wf.mean_gain = 10.0;

  // Secondary pulse (pileup) configuration
  wf.sp_enable = true;
  wf.sp_delay = 20.0;
  wf.sp_amplitude_frac = 0.5;

  gen.setDigitizerParams(daq);
  gen.setDetectorParams(wf);

  FADCStandaloneAnalyzer analyzer;
  FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
  analyzer.setConfig(config);

  TemplateAnalysisConfig tmplCfg;
  // D2 specific configurations
  tmplCfg.fSmoothHalfWidth = 2; // Gaussian pre-smoothing kernel half-width
  tmplCfg.fD2ThresholdNSigma =
      5.0;                   // Significance level for auto-scaled d2 threshold
  tmplCfg.fD2MinPeakSep = 1; // Minimum separation between pulse candidates
  tmplCfg.fAmpCutADC = 8.0;
  tmplCfg.fFitMargin = 8;
  tmplCfg.fScanHalf = 15.0;

  // Ensure output directory exists
  gSystem->mkdir("../root_files", kTRUE);

  TString rootFileName = "../root_files/test_d2_intermediate_steps.root";
  TFile outfile(rootFileName, "RECREATE");
  if (outfile.IsZombie()) {
    // Fallback to current directory if parent path is not writable
    rootFileName = "test_d2_intermediate_steps.root";
    outfile.Open(rootFileName, "RECREATE");
  }

  std::cout << "\n======================================================\n";
  std::cout << "    FADC AnalyzeWithD2 INTERMEDIATE STEPS TEST        \n";
  std::cout << "======================================================\n";
  std::cout << "Generating and analyzing " << NUM_WAVEFORMS
            << " waveforms...\n";
  std::cout << "Output ROOT file: " << rootFileName << "\n\n";

  // Diagnostic histograms
  TH1F h_int_diff("h_int_diff",
                  "Integral Diff (Fitted - Truth);ADC#timesSample;Entries", 100,
                  -100, 100);
  TH1F h_amp_diff("h_amp_diff",
                  "Amplitude Diff (Fitted - Truth);ADC counts;Entries", 100,
                  -50, 50);
  TH1F h_time_diff("h_time_diff",
                   "Time Diff (Fitted CFD - Truth);Subsamples;Entries", 100,
                   -50, 50);

  int total_pulses_found = 0;

  for (int i = 0; i < NUM_WAVEFORMS; ++i) {
    std::cout << "── Processing Waveform " << i
              << " ──────────────────────────────\n";

    // 1. Generate waveform with truth info
    std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
    std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
    const TruthPulseInfo &wfInfo = gen.getTruthInfo();
    const Int_t N = static_cast<Int_t>(samples.size());

    // 2. Run AnalyzeWithD2
    FADCPulseWF d2Result = analyzer.AnalyzeWithD2(samples, tmplCfg);
    total_pulses_found += d2Result.nPulses;

    std::cout << "  Pedestal: " << std::fixed << std::setprecision(2)
              << d2Result.pedestal << " ADC counts\n";
    std::cout << "  Truth Pulses: " << wfInfo.pulseIntegral.size()
              << " | Fitted Pulses: " << d2Result.nPulses << "\n";

    for (size_t p = 0; p < d2Result.pulses.size(); ++p) {
      const auto &pulse = d2Result.pulses[p];
      std::cout << "    Pulse [" << p
                << "]: Amplitude = " << std::setprecision(2) << pulse.amplitude
                << " ADC, t0 = " << pulse.timeSample
                << " samples, Integral = " << pulse.integral
                << ", sigma = " << pulse.fittedSigma
                << ", tau = " << pulse.fittedTau
                << ", chi2/ndf = " << pulse.chi2ndf << "\n";

      if (p < wfInfo.pulseIntegral.size()) {
        h_int_diff.Fill(pulse.integral - wfInfo.pulseIntegral[p]);
        h_amp_diff.Fill(pulse.amplitude -
                        (wfInfo.pulseAmplitude[p] - daq.pedestal));
        h_time_diff.Fill(pulse.timeHalfMax * daq.tdc_subsamples -
                         wfInfo.pulseTime[p]);
      }
    }

    // ── Build Plots for the 4 Intermediate Steps ──────────────────────────

    // ── Step 1: Raw Waveform ──────────────────────────────────────────────
    TGraph *gRaw = new TGraph(N);
    gRaw->SetName(Form("g_raw_wf_%d", i));
    gRaw->SetTitle(
        Form("Step 1: Raw Waveform (Waveform %d);Sample Index;ADC Counts", i));
    for (Int_t j = 0; j < N; ++j) {
      gRaw->SetPoint(j, static_cast<Double_t>(j), samples[j]);
    }
    gRaw->SetMarkerStyle(20);
    gRaw->SetMarkerSize(0.6);
    gRaw->SetMarkerColor(kBlue + 1);
    gRaw->SetLineColor(kBlue - 7);
    gRaw->SetLineWidth(1);

    // Pedestal Line for Step 1
    TLine *linePedStep1 =
        new TLine(0.0, d2Result.pedestal, N - 1.0, d2Result.pedestal);
    linePedStep1->SetLineColor(kGray + 2);
    linePedStep1->SetLineStyle(2);
    linePedStep1->SetLineWidth(2);

    // ── Step 2: Smoothed Waveform ─────────────────────────────────────────
    TGraph *gSmooth = new TGraph(N);
    gSmooth->SetName(Form("g_smooth_wf_%d", i));
    gSmooth->SetTitle(
        Form("Step 2: Gaussian Smoothed Waveform (#sigma_{kernel} = %.1f) "
             "(Waveform %d);Sample Index;ADC Counts",
             tmplCfg.fSmoothHalfWidth / 2.0, i));

    // Plot smoothed waveform with pedestal baseline restored for direct
    // comparison
    for (Int_t j = 0; j < N; ++j) {
      Double_t smoothVal = (j < static_cast<Int_t>(d2Result.wfSmooth.size()))
                               ? (d2Result.wfSmooth[j] + d2Result.pedestal)
                               : samples[j];
      gSmooth->SetPoint(j, static_cast<Double_t>(j), smoothVal);
    }
    gSmooth->SetMarkerStyle(21);
    gSmooth->SetMarkerSize(0.5);
    gSmooth->SetMarkerColor(kTeal + 2);
    gSmooth->SetLineColor(kTeal + 2);
    gSmooth->SetLineWidth(2);

    TLine *linePedStep2 =
        new TLine(0.0, d2Result.pedestal, N - 1.0, d2Result.pedestal);
    linePedStep2->SetLineColor(kGray + 2);
    linePedStep2->SetLineStyle(2);
    linePedStep2->SetLineWidth(2);

    // ── Step 3: D2 Result (Second-Derivative) ─────────────────────────────
    TGraph *gD2 = new TGraph(N);
    gD2->SetName(Form("g_d2_wf_%d", i));
    gD2->SetTitle(
        Form("Step 3: Second Derivative Trace d^{2}y/dt^{2} (Waveform "
             "%d);Sample Index;d^{2}y/dt^{2} (ADC/sample^{2})",
             i));
    for (Int_t j = 0; j < N; ++j) {
      Double_t d2Val = (j < static_cast<Int_t>(d2Result.d2Out.size()))
                           ? d2Result.d2Out[j]
                           : 0.0;
      gD2->SetPoint(j, static_cast<Double_t>(j), d2Val);
    }
    gD2->SetMarkerStyle(24);
    gD2->SetMarkerSize(0.6);
    gD2->SetMarkerColor(kMagenta + 2);
    gD2->SetLineColor(kMagenta + 1);
    gD2->SetLineWidth(2);

    // Zero-baseline line for D2
    TLine *lineZeroD2 = new TLine(0.0, 0.0, N - 1.0, 0.0);
    lineZeroD2->SetLineColor(kBlack);
    lineZeroD2->SetLineStyle(3);
    lineZeroD2->SetLineWidth(1);

    // D2 Threshold line estimation
    Int_t nPed = config.fNPedestalSamples;
    Int_t d2PedN = std::min(nPed, N - 1);
    Double_t d2Mean = 0.0, d2Var = 0.0;
    for (Int_t j = 1; j <= d2PedN; ++j) {
      d2Mean += d2Result.d2Out[j];
    }
    d2Mean /= (d2PedN > 0 ? d2PedN : 1);
    for (Int_t j = 1; j <= d2PedN; ++j) {
      Double_t diff = d2Result.d2Out[j] - d2Mean;
      d2Var += diff * diff;
    }
    d2Var /= (d2PedN > 1 ? (d2PedN - 1) : 1);
    Double_t d2Sigma = std::sqrt(d2Var);
    Double_t d2ThreshVal = tmplCfg.fD2ThresholdNSigma * d2Sigma;
    if (d2ThreshVal < 0.1)
      d2ThreshVal = 0.1;

    TLine *lineD2Thresh = new TLine(0.0, -d2ThreshVal, N - 1.0, -d2ThreshVal);
    lineD2Thresh->SetLineColor(kOrange + 7);
    lineD2Thresh->SetLineStyle(2);
    lineD2Thresh->SetLineWidth(2);

    // ── Step 4: EMG Model Fitting Overlay ─────────────────────────────────
    TGraph *gSamplesForFit = new TGraph(N);
    for (Int_t j = 0; j < N; ++j) {
      gSamplesForFit->SetPoint(j, static_cast<Double_t>(j), samples[j]);
    }
    gSamplesForFit->SetMarkerStyle(20);
    gSamplesForFit->SetMarkerSize(0.6);
    gSamplesForFit->SetMarkerColor(kBlue + 1);
    gSamplesForFit->SetLineColor(kBlue - 9);

    // Total fit curve
    TF1 fTotal(
        Form("fTotal_wf_%d", i),
        [d2Result](double *x, double * /*p*/) {
          double t = x[0];
          double val = d2Result.pedestal;
          for (const auto &pulse : d2Result.pulses) {
            if (pulse.valid && pulse.fittedSigma > 0 && pulse.fittedTau > 0) {
              double dt = t - pulse.timeSample;
              double sigma = pulse.fittedSigma;
              double tau = pulse.fittedTau;
              double s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
              double erf_arg =
                  (dt - sigma * sigma / tau) / (sigma * std::sqrt(2.0));
              double raw =
                  std::exp(s2_2t2 - dt / tau) * (1.0 + std::erf(erf_arg));
              double dt_peak = sigma * sigma / tau;
              double raw_peak =
                  std::exp(s2_2t2 - dt_peak / tau) * (1.0 + std::erf(0.0));
              if (raw_peak > 1e-12) {
                val += pulse.amplitude * raw / raw_peak;
              }
            }
          }
          return val;
        },
        -0.5, N - 0.5, 0);

    const double fitStep = 0.2;
    const double fitXMin = -0.5;
    const double fitXMax = N - 0.5;
    int nFitPts = static_cast<int>((fitXMax - fitXMin) / fitStep) + 1;
    TGraph *gFitTotal = new TGraph(nFitPts);
    gFitTotal->SetName(Form("g_fit_total_wf_%d", i));
    for (int k = 0; k < nFitPts; ++k) {
      double xv = fitXMin + k * fitStep;
      gFitTotal->SetPoint(k, xv, fTotal.Eval(xv));
    }
    gFitTotal->SetLineColor(kRed + 1);
    gFitTotal->SetLineWidth(2);

    TMultiGraph *mgFit = new TMultiGraph();
    mgFit->SetName(Form("mg_fit_wf_%d", i));
    mgFit->SetTitle(Form(
        "Step 4: EMG Model Fit Overlay (Waveform %d);Sample Index;ADC Counts",
        i));
    mgFit->Add(gSamplesForFit, "P");
    mgFit->Add(gFitTotal, "L");

    // Add individual pulse components if multiple pulses are fitted
    std::vector<TGraph *> gPulses;
    const Color_t pulseColors[] = {kGreen + 2, kViolet + 1, kOrange + 2,
                                   kCyan + 2};
    for (size_t p = 0; p < d2Result.pulses.size(); ++p) {
      const auto &pulse = d2Result.pulses[p];
      if (!pulse.valid)
        continue;

      TF1 fSingle(
          Form("fSingle_wf_%d_p_%zu", i, p),
          [d2Result, pulse](double *x, double * /*p*/) {
            double t = x[0];
            double val = d2Result.pedestal;
            if (pulse.fittedSigma > 0 && pulse.fittedTau > 0) {
              double dt = t - pulse.timeSample;
              double sigma = pulse.fittedSigma;
              double tau = pulse.fittedTau;
              double s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
              double erf_arg =
                  (dt - sigma * sigma / tau) / (sigma * std::sqrt(2.0));
              double raw =
                  std::exp(s2_2t2 - dt / tau) * (1.0 + std::erf(erf_arg));
              double dt_peak = sigma * sigma / tau;
              double raw_peak =
                  std::exp(s2_2t2 - dt_peak / tau) * (1.0 + std::erf(0.0));
              if (raw_peak > 1e-12) {
                val += pulse.amplitude * raw / raw_peak;
              }
            }
            return val;
          },
          -0.5, N - 0.5, 0);

      TGraph *gPulse = new TGraph(nFitPts);
      gPulse->SetName(Form("g_pulse_wf_%d_%zu", i, p));
      for (int k = 0; k < nFitPts; ++k) {
        double xv = fitXMin + k * fitStep;
        gPulse->SetPoint(k, xv, fSingle.Eval(xv));
      }
      Color_t pCol = pulseColors[p % 4];
      gPulse->SetLineColor(pCol);
      gPulse->SetLineStyle(7); // dashed
      gPulse->SetLineWidth(2);
      gPulses.push_back(gPulse);
      mgFit->Add(gPulse, "L");
    }

    // ── Create 2x2 Canvas for all 4 Intermediate Steps ───────────────────
    TCanvas *cSteps = new TCanvas(
        Form("c_steps_wf_%d", i),
        Form("AnalyzeWithD2 Intermediate Steps (Waveform %d)", i), 1000, 800);
    cSteps->Divide(2, 2);

    // Pad 1: Step 1 Raw Waveform
    cSteps->cd(1);
    gPad->SetGrid(1, 1);
    gRaw->Draw("ALP");
    linePedStep1->Draw("SAME");
    TLegend *leg1 = new TLegend(0.60, 0.72, 0.88, 0.88);
    leg1->AddEntry(gRaw, "Raw ADC Samples", "lp");
    leg1->AddEntry(linePedStep1, Form("Pedestal (%.1f)", d2Result.pedestal),
                   "l");
    leg1->Draw();

    // Pad 2: Step 2 Smoothed Waveform
    cSteps->cd(2);
    gPad->SetGrid(1, 1);
    gSmooth->Draw("ALP");
    linePedStep2->Draw("SAME");
    TLegend *leg2 = new TLegend(0.60, 0.72, 0.88, 0.88);
    leg2->AddEntry(gSmooth, "Smoothed Trace", "lp");
    leg2->AddEntry(linePedStep2, Form("Pedestal (%.1f)", d2Result.pedestal),
                   "l");
    leg2->Draw();

    // Pad 3: Step 3 D2 Result
    cSteps->cd(3);
    gPad->SetGrid(1, 1);
    gD2->Draw("ALP");
    lineZeroD2->Draw("SAME");
    lineD2Thresh->Draw("SAME");
    TLegend *leg3 = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg3->AddEntry(gD2, "d^{2}y/dt^{2} Trace", "lp");
    leg3->AddEntry(lineD2Thresh, Form("Threshold (-%.2f)", d2ThreshVal), "l");
    leg3->Draw();

    // Draw candidate peak markers on Pad 3
    for (const auto &pulse : d2Result.pulses) {
      if (pulse.valid && pulse.fittedSigma > 0 && pulse.fittedTau > 0) {
        Double_t peakPos =
            pulse.timeSample +
            (pulse.fittedSigma * pulse.fittedSigma) / pulse.fittedTau;
        Int_t pIdx = std::max(
            0, std::min(N - 1, static_cast<Int_t>(std::round(peakPos))));
        Double_t peakD2 = d2Result.d2Out[pIdx];
        TMarker *m = new TMarker(peakPos, peakD2, 29); // Star marker
        m->SetMarkerColor(kRed);
        m->SetMarkerSize(1.5);
        m->Draw("SAME");
      }
    }

    // Pad 4: Step 4 EMG Fit Overlay
    cSteps->cd(4);
    gPad->SetGrid(1, 1);
    mgFit->Draw("A");
    TLegend *leg4 = new TLegend(0.55, 0.65, 0.88, 0.88);
    leg4->AddEntry(gSamplesForFit, "Raw Samples", "p");
    leg4->AddEntry(gFitTotal, "Total EMG Fit", "l");
    for (size_t p = 0; p < gPulses.size(); ++p) {
      leg4->AddEntry(gPulses[p], Form("Pulse %zu Component", p + 1), "l");
    }
    leg4->Draw();

    // ── Write objects to ROOT File ────────────────────────────────────────
    TDirectory *dir = outfile.mkdir(Form("Waveform_%d", i));
    dir->cd();
    gRaw->Write();
    gSmooth->Write();
    gD2->Write();
    mgFit->Write();
    cSteps->Write("c_steps");

    outfile.cd();
    cSteps->Write(Form("c_steps_wf_%d", i));

    // Validations for each run
    REQUIRE(d2Result.nPulses >= 1);
    REQUIRE(d2Result.pulses[0].valid == kTRUE);
  }

  outfile.cd();
  h_int_diff.Write();
  h_amp_diff.Write();
  h_time_diff.Write();
  outfile.Close();

  std::cout << "\n======================================================\n";
  std::cout << "                 SUMMARY & RESULTS                    \n";
  std::cout << "======================================================\n";
  std::cout << "Total Waveforms Analyzed : " << NUM_WAVEFORMS << "\n";
  std::cout << "Total Pulses Detected    : " << total_pulses_found << "\n";
  std::cout << "Saved ROOT File          : " << rootFileName << "\n";
  std::cout << "======================================================\n\n";

  REQUIRE(total_pulses_found >= NUM_WAVEFORMS);
}
