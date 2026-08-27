#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include "waveformGenerator.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <iomanip>
#include <iostream>

// ROOT Headers
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

TEST_CASE("Single pulse without noise", "[basic1]") {
  waveformGenerator gen;
  DigitizerParams daq;
  daq.noise_sigma = 2.5;
  daq.num_samples = 74;
  daq.max_n_pulses = 1;

  WaveformParams wf;
  wf.sigma_g = 3;
  wf.tau = 10;
  wf.mean_gain = 10.0;

  // Secondary pulse
  wf.sp_enable = true;
  wf.sp_delay = 30.0;
  wf.sp_amplitude_frac = 0.5;

  gen.setDigitizerParams(daq);
  gen.setDetectorParams(wf);

  FADCStandaloneAnalyzer analyzer;
  FADCConfig config = {4, 50, 2, 5, 2, 12}; // Match standard params
  analyzer.setConfig(config);

  TemplateAnalysisConfig tmplCfg;
  tmplCfg.fMFSeparation = 5;
  tmplCfg.fScanHalf = 6; // ±6 from corrected candidateBin = peak - peakBin
  tmplCfg.fMFLeft = 8;
  tmplCfg.fMFRight = 8;
  tmplCfg.fMFStart = 0; // search from sample 0 (MF loop now starts at 0)
  tmplCfg.fAmpCutADC = 8.0;
  tmplCfg.fFitMargin = 8;

  SPEParams sp;
  sp.sigma = 3;
  sp.tau = 10;

  FADCTemplate tmpl = FADCTemplate::fromWaveformParams(sp);

  std::vector<Int_t> rawWaveform = gen.generateWaveform(PulseShape::SPE);
  std::vector<Double_t> samples(rawWaveform.begin(), rawWaveform.end());
  const TruthPulseInfo &wfInfo = gen.getTruthInfo();

  FADCPulseWF templateResult =
      analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

  // Error margins
  Double_t int_margin = wf.mean_gain * 3.0; // Sould be within N # of samples

  std::cout << "wfGen integral: " << wfInfo.pulseIntegral[0] << std::endl;
  std::cout << "templateResult integral: " << templateResult.pulses[0].integral
            << std::endl;
  std::cout << "int_margin: " << int_margin << std::endl;
  std::cout << "fit stage: " << templateResult.pulses[0].fitStage << std::endl;
  std::cout << "========================================\n\n";

  std::cout << "wfGen amplitude: " << wfInfo.pulseAmplitude[0] - daq.pedestal
            << std::endl;
  std::cout << "templateResult amplitude: "
            << templateResult.pulses[0].amplitude << std::endl;
  std::cout << "========================================\n\n";

  std::cout << "wfGen time: " << wfInfo.pulseTime[0] << std::endl;
  std::cout << "templateResult time: "
            << templateResult.pulses[0].timeHalfMax * daq.tdc_subsamples
            << std::endl;
  std::cout << "========================================\n\n";

  REQUIRE(templateResult.pulses.size() >= 1);
  REQUIRE(templateResult.pulses[0].valid);
  // REQUIRE(std::abs(templateResult.pulses[0].integral -
  // wfInfo.pulseIntegral[0]) < int_margin);

  // Create TApplication if ROOT hasn't made one yet (needed for the GUI event
  // loop).
  int argc = 0;
  char **argv = nullptr;
  if (!gApplication) {
    int argc = 1;
    static char appName[] = "TestApp";
    static char *argv[] = {appName, nullptr};
    new TApplication("TestApp", &argc, argv);
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Single pulse template fit;Sample;ADC counts");
  TCanvas *c = new TCanvas("c", "basic1 fit", 800, 500);
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
      [tmpl, templateResult](double *x, double * /*p*/) {
        double t = x[0];
        double val = templateResult.pedestal;
        for (const auto &pulse : templateResult.pulses) {
          if (pulse.valid) {
            if (pulse.fitStage == 1 || pulse.fittedSigma == 0.0) {
              // Stage 1: fixed template shape
              val += pulse.amplitude * tmpl.eval(t - pulse.timeSample);
            } else {
              // Stage 2: fitted EMG with custom sigma/tau
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
                double raw_peak =
                    std::exp(s2_2t2 - dt_peak / tau) * (1.0 + std::erf(0.0));
                if (raw_peak > 1e-12)
                  val += pulse.amplitude * raw / raw_peak;
              }
            }
          }
        }
        return val;
      },
      -0.5, N - 0.5, 0);

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

  mg->Add(gSamples, "P");
  mg->Add(gFit, "L");

  c->cd();
  mg->Draw("AP");
  c->Update(); // flush to the display before Run() takes over the event loop

  // Poll the event loop until the canvas is closed by the user.
  // When the window closes, ROOT removes it from gROOT->GetListOfCanvases(),
  // so the loop exits cleanly and the test binary terminates normally.
  while (gROOT->GetListOfCanvases()->FindObject("c")) {
    gSystem->ProcessEvents(); // dispatch X11/GUI events
    gSystem->Sleep(20);       // 20 ms sleep to avoid busy-waiting
  }
}