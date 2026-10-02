#include "FADCStandaloneAnalyzer.h"
#include "waveformGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

// ROOT includes
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMultiGraph.h>
#include <TSystem.h>

TEST_CASE("AnalyzeWithTemplate 6 Intermediate Stages", "[Template_Stages]") {
  // ── Setup Digitizer and Detector Parameters ──
  waveformGenerator gen;
  DigitizerParams daq;
  daq.noise_sigma = 1.5;
  daq.num_samples = 74;
  daq.sample_rate = 4.0; // ns/sample
  daq.pedestal = 100.0;
  daq.max_n_pulses = 1;  // Exactly 1 primary pulse (sp_enable adds overlapping secondary)

  WaveformParams wf;
  wf.sigma_g = 2.6;
  wf.tau = 7.0;
  wf.mean_gain = 15.0;
  wf.sp_enable = true; // Enable secondary pulse for pileup testing
  wf.sp_delay = 20.0;  // 20 ns separation (overlapping pileup)
  wf.sp_amplitude_frac = 0.5;

  gen.setDigitizerParams(daq);
  gen.setDetectorParams(wf);

  // ── Build Template from analytical SPE pulse ──
  const Int_t tmplNSamples = 50;
  const Double_t tmplT0_ns = 20.0; // Place onset well after pedestal window (samples 0..3 are pure baseline)
  const Int_t tmplNPed = daq.init_ped_samples;

  std::vector<Double_t> tmplAnalog(tmplNSamples, 0.0);
  Double_t rawPeak = 0.0;
  for (Int_t i = 0; i < tmplNSamples; ++i) {
    Double_t t = i * daq.sample_rate;
    tmplAnalog[i] = WaveformUtil::spePulse(t, tmplT0_ns, wf.sigma_g, wf.tau);
    if (tmplAnalog[i] > rawPeak)
      rawPeak = tmplAnalog[i];
  }

  std::vector<Double_t> cleanWaveform(tmplNSamples);
  for (Int_t i = 0; i < tmplNSamples; ++i) {
    Double_t amplitude_mV =
        (rawPeak > 1e-12) ? wf.mean_gain * tmplAnalog[i] / rawPeak : 0.0;
    cleanWaveform[i] = std::round(amplitude_mV / daq.v_lsb + daq.pedestal);
  }
  FADCTemplate tmpl = FADCTemplate::fromWaveformSamples(
      cleanWaveform, tmplNPed, 5, 5);

  FADCStandaloneAnalyzer analyzer;
  TemplateAnalysisConfig cfg;
  cfg.fMFLeft = 5;
  cfg.fMFRight = 5;
  cfg.fMFStart = 5;
  cfg.fMFEnd = daq.num_samples - 1;
  cfg.fMFThreshold = 0.1;
  cfg.fMFSeparation = 3;
  cfg.fAmpCutADC = 5.0;
  cfg.fFitMargin = 6;
  cfg.fScanHalf = 15.0;
  cfg.fRefitMaxIter = 3;

  gSystem->mkdir("../root_files/template_an_test1", kTRUE);
  TFile outfile("../root_files/template_an_test1/test_template_6_stages.root",
                "RECREATE");

  // ── Generate waveform and run analysis ──
  std::vector<Int_t> raw = gen.generateWaveform(PulseShape::SPE);
  std::vector<Double_t> samples(raw.begin(), raw.end());
  const TruthPulseInfo &truth = gen.getTruthInfo();
  Int_t N = samples.size();

  // Run full analysis
  FADCPulseWF result = analyzer.AnalyzeWithTemplate(samples, tmpl, cfg);

  // Helper to evaluate a fitted pulse model (handles Stage 1 template & Stage 2 EMG)
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

  // -------------------------------------------------------------------------
  // VALIDATION & CANVAS PLOTTING FOR ALL 6 STAGES
  // -------------------------------------------------------------------------
  TCanvas *c6 =
      new TCanvas("c_6_stages", "AnalyzeWithTemplate 6 Stages", 1200, 800);
  c6->Divide(3, 2);

  // ── STAGE 1: Pedestal & Baseline Subtraction ─────────────────────────────
  INFO("Testing Stage 1: Pedestal Estimation");
  REQUIRE(result.pedestal ==
          Catch::Approx(daq.pedestal).margin(3.0 * daq.noise_sigma));
  REQUIRE(result.wfPS.size() == samples.size());

  c6->cd(1);
  gPad->SetGrid(1, 1);
  TGraph *gRaw = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gRaw->SetPoint(j, j, samples[j]);
  gRaw->SetTitle("Stage 1: Raw Waveform & Pedestal;Sample;ADC");
  gRaw->SetMarkerStyle(20);
  gRaw->Draw("ALP");

  TLine *lPed = new TLine(0, result.pedestal, N - 1, result.pedestal);
  lPed->SetLineColor(kRed);
  lPed->SetLineStyle(2);
  lPed->Draw("SAME");

  // ── STAGE 2: Matched Filter Output ───────────────────────────────────────
  INFO("Testing Stage 2: Matched Filter Output");
  REQUIRE(result.mfOut.size() == samples.size());

  c6->cd(2);
  gPad->SetGrid(1, 1);
  TGraph *gMF = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gMF->SetPoint(j, j, result.mfOut[j]);
  gMF->SetTitle("Stage 2: Matched Filter Output;Sample;MF Response");
  gMF->SetLineColor(kAzure + 2);
  gMF->Draw("AL");

  TLine *lThresh = new TLine(0, cfg.fMFThreshold, N - 1, cfg.fMFThreshold);
  lThresh->SetLineColor(kOrange + 2);
  lThresh->SetLineStyle(2);
  lThresh->Draw("SAME");

  // ── STAGE 3: Peak Finding in MF ──────────────────────────────────────────
  INFO("Testing Stage 3: Peak Detection in MF");
  REQUIRE(result.nPulses >= 1);

  c6->cd(3);
  gPad->SetGrid(1, 1);
  TGraph *gMFPeaks = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gMFPeaks->SetPoint(j, j, result.mfOut[j]);
  gMFPeaks->SetTitle("Stage 3: MF Peaks Identified;Sample;MF Response");
  gMFPeaks->Draw("AL");

  for (const auto &p : result.pulses) {
    Double_t peakPos = (p.fitStage == 1 || p.fittedSigma == 0.0)
                           ? (p.timeSample + tmpl.peakBin)
                           : (p.timeSample + p.fittedSigma * p.fittedSigma / p.fittedTau);
    Int_t idx = std::max(0, std::min(N - 1, static_cast<Int_t>(std::round(peakPos))));
    Double_t mfVal = result.mfOut[idx];
    TMarker *mk = new TMarker(peakPos, mfVal, 29); // Star marker
    mk->SetMarkerColor(kRed);
    mk->SetMarkerSize(1.5);
    mk->Draw("SAME");
  }

  // ── STAGE 4: Sequential Peeling & Residual ───────────────────────────────
  INFO("Testing Stage 4: Sequential Peeling & Residual");
  c6->cd(4);
  gPad->SetGrid(1, 1);
  std::vector<Double_t> wfResidual = result.wfPS;
  if (!result.pulses.empty()) {
    const auto &p1 = result.pulses[0];
    for (int j = 0; j < N; ++j) {
      wfResidual[j] -= evalPulseModel(j, p1);
    }
  }
  TGraph *gRes = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gRes->SetPoint(j, j, wfResidual[j]);
  gRes->SetTitle(
      "Stage 4: Residual After 1st Pulse Subtraction;Sample;ADC");
  gRes->SetLineColor(kGreen + 2);
  gRes->Draw("AL");

  // ── STAGE 5: Refit Comparison (With vs Without Iterations) ───────────────
  INFO("Testing Stage 5: Coordinate Descent Convergence");
  TemplateAnalysisConfig cfgNoRefit = cfg;
  cfgNoRefit.fRefitMaxIter = 0; // Disable stage 5
  FADCPulseWF resultNoRefit =
      analyzer.AnalyzeWithTemplate(samples, tmpl, cfgNoRefit);

  if (result.nPulses >= 2 && resultNoRefit.nPulses >= 2) {
    Double_t chi2_with_refit = 0.0, chi2_no_refit = 0.0;
    for (const auto &p : result.pulses)
      chi2_with_refit += p.chi2ndf;
    for (const auto &p : resultNoRefit.pulses)
      chi2_no_refit += p.chi2ndf;

    CHECK(chi2_with_refit <= chi2_no_refit + 0.1);
  }

  c6->cd(5);
  gPad->SetGrid(1, 1);
  TGraph *gPS = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gPS->SetPoint(j, j, result.wfPS[j]);
  gPS->SetTitle("Stage 5: Refitted Pulse Components;Sample;ADC");
  gPS->Draw("AP");

  for (size_t p = 0; p < result.pulses.size(); ++p) {
    TGraph *gComp = new TGraph(N);
    for (int j = 0; j < N; ++j) {
      gComp->SetPoint(j, j, evalPulseModel(j, result.pulses[p]));
    }
    gComp->SetLineColor(p == 0 ? kRed : kMagenta);
    gComp->SetLineWidth(2);
    gComp->Draw("L SAME");
  }

  // ── STAGE 6: Chronological Sorting & Quality Cuts ────────────────────────
  INFO("Testing Stage 6: Chronological Ordering & Quality Cuts");
  for (size_t p = 1; p < result.pulses.size(); ++p) {
    REQUIRE(result.pulses[p].timeSample >= result.pulses[p - 1].timeSample);
    REQUIRE(result.pulses[p].chi2ndf <= cfg.fMaxChi2NDF);
  }

  c6->cd(6);
  gPad->SetGrid(1, 1);
  TMultiGraph *mgFinal = new TMultiGraph();
  mgFinal->SetTitle("Stage 6: Final Fitted Overlay;Sample;ADC");

  TGraph *gFinalRaw = new TGraph(N);
  for (int j = 0; j < N; ++j)
    gFinalRaw->SetPoint(j, j, samples[j]);
  mgFinal->Add(gFinalRaw, "P");

  TGraph *gFinalFit = new TGraph(N);
  for (int j = 0; j < N; ++j) {
    Double_t val = result.pedestal;
    for (const auto &p : result.pulses) {
      val += evalPulseModel(j, p);
    }
    gFinalFit->SetPoint(j, j, val);
  }
  gFinalFit->SetLineColor(kRed);
  gFinalFit->SetLineWidth(2);
  mgFinal->Add(gFinalFit, "L");
  mgFinal->Draw("A");

  c6->Write();
  outfile.Close();
}
