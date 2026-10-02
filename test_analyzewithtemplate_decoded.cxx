/*
 * test_analyzewithtemplate_decoded.cxx
 *
 * Analog of test_analyzewithtemplate.cxx but reads real decoded FADC
 * waveforms from:
 *   /MyFiles/StandAloneAnalyzer/wm_data/decoded_data/
 * instead of using the waveformGenerator.
 *
 * Because there is no truth information available for real data, the
 * per-waveform validation section (integral/amplitude/time comparison to
 * truth) is replaced by:
 *   - Statistics histograms of the reconstructed quantities
 *       h_nPulses      - number of pulses found per waveform
 *       h_amplitude    - fitted amplitude of the first valid pulse
 *       h_pedestal     - estimated pedestal value per waveform
 *       h_chi2ndf      - chi2/ndf of the best fit per waveform
 *   - Per-event waveform overlays written to the ROOT file for ALL
 *     waveforms (not just failures), named "Waveform_ev<N>_ch<C>"
 *
 * Branch layout follows plot_waveforms.C:
 *   TTree "T"  in the decoded ROOT file
 *     Ndata.moller.scint.samps        (Int_t)
 *     moller.scint.samps              (Double_t[])
 *     Ndata.moller.scint.samps_idx    (Int_t)
 *     moller.scint.samps_idx          (Double_t[])
 *     Ndata.moller.scint.samps_elemID (Int_t)
 *     moller.scint.samps_elemID       (Double_t[])
 *
 * Compile & run (from the test1/ directory):
 *   g++ -std=c++17 test_main.cxx test_analyzewithtemplate_decoded.cxx \
 *       FADCStandaloneAnalyzer.cxx waveformGenerator.cxx \
 *       -lCatch2Main -lCatch2 $(root-config --cflags --libs) -lMinuit2 \
 *       -o run_decoded_template && ./run_decoded_template
 *
 * Output: ../root_files/decoded_analyzewithtemplate_chan5.root
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
// waveformGenerator.h not needed here: SPEParams is defined in
// FADCStandaloneAnalyzer.h

// ROOT Headers
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

// ── Decoded data file
// ─────────────────────────────────────────────────────────
static const char *kInputFile =
    "/home/lasitha/MyFiles/StandAloneAnalyzer/wm_data/decoded_data/"
    "moller_fadc_replayed_329_seg0_1_08_10_26.root";

// Maximum array sizes (same as plot_waveforms.C)
static const int kMaxSamps = 20000;
static const int kMaxChans = 1000;

// ═══════════════════════════════════════════════════════════════════════════════
// Catch2 test case
// ═══════════════════════════════════════════════════════════════════════════════

TEST_CASE("AnalyzeWithTemplate on Decoded FADC Data",
          "[DecodedData][Statistics]") {

  // ── 1.  Open the decoded ROOT file ───────────────────────────────────────
  TFile *inFile = TFile::Open(kInputFile);
  REQUIRE(inFile != nullptr);
  REQUIRE(!inFile->IsZombie());

  TTree *T = dynamic_cast<TTree *>(inFile->Get("T"));
  REQUIRE(T != nullptr);

  // ── 2.  Branch setup (mirrors plot_waveforms.C) ───────────────────────────
  T->SetBranchStatus("*", 0);
  T->SetBranchStatus("Ndata.moller.scint.samps", 1);
  T->SetBranchStatus("moller.scint.samps", 1);
  T->SetBranchStatus("Ndata.moller.scint.samps_idx", 1);
  T->SetBranchStatus("moller.scint.samps_idx", 1);
  T->SetBranchStatus("Ndata.moller.scint.samps_elemID", 1);
  T->SetBranchStatus("moller.scint.samps_elemID", 1);

  Int_t nsamps = 0;
  Int_t nsamps_idx = 0;
  Int_t nsamps_elemID = 0;
  Double_t samps[kMaxSamps];
  Double_t samps_idx[kMaxChans];
  Double_t samps_elemID[kMaxChans];

  T->SetBranchAddress("Ndata.moller.scint.samps", &nsamps);
  T->SetBranchAddress("moller.scint.samps", samps);
  T->SetBranchAddress("Ndata.moller.scint.samps_idx", &nsamps_idx);
  T->SetBranchAddress("moller.scint.samps_idx", samps_idx);
  T->SetBranchAddress("Ndata.moller.scint.samps_elemID", &nsamps_elemID);
  T->SetBranchAddress("moller.scint.samps_elemID", samps_elemID);

  // ── 3.  Template & analysis configuration ────────────────────────────────
  // Construct a 4-sample rising and 10-sample lowering pulse template
  const int nRise = 3;
  const int nFall = 8;
  std::vector<Double_t> rawShape;
  rawShape.reserve(nRise + nFall + 1);

  // 4 samples rising to peak (0.0 to 1.0)
  for (int i = 0; i <= nRise; ++i) {
    rawShape.push_back(static_cast<Double_t>(i) / nRise);
  }
  // 10 samples lowering from peak to zero (1.0 to 0.0)
  for (int i = 1; i <= nFall; ++i) {
    rawShape.push_back(1.0 - static_cast<Double_t>(i) / nFall);
  }

  // Create template using fromRawSamples:
  // peak is at index 4 (tmpl.peakBin = 4)
  // mfLeft = 4 (covers the 4-sample rising edge), mfRight = 5
  FADCTemplate tmpl = FADCTemplate::fromRawSamples(
      rawShape, kTRUE, /* mfLeft */ 4, /* mfRight */ 5);
  REQUIRE(!tmpl.shape.empty());
  REQUIRE(tmpl.peakBin >= 0);

  FADCStandaloneAnalyzer analyzer;
  FADCConfig config = {4, 50, 2, 5, 2, 12};
  analyzer.setConfig(config);

  TemplateAnalysisConfig tmplCfg;
  tmplCfg.fMFLeft = 4;         // samples left of peak in MF kernel
  tmplCfg.fMFRight = 5;        // samples right of peak in MF kernel
  tmplCfg.fMFStart = 1;        // first sample of search window
  tmplCfg.fMFThreshold = 10.0; // MF amplitude threshold (ADC counts)
  tmplCfg.fMFSeparation = 4;   // min separation between pulse candidates
  tmplCfg.fAmpCutADC = 8.0;    // amplitude cut
  tmplCfg.fFitMargin = 10;     // fit window margin around peak
  tmplCfg.fScanHalf = 15.0;    // scan range
  tmplCfg.fRefitMaxIter = 3;   // iterative refitting

  SECTION("Multiple iterations on decoded data") {

    // ── How many waveforms to process ────────────────────────────────────
    //   Set num_wfs to a large number (e.g. 10000) to sweep the whole file,
    //   or a smaller value for a quick sanity check.
    const int num_wfs = 1000;

    // Which channel index to follow across events (-1 = all channels in
    // each event, same semantics as plot_waveforms.C cycleChannels mode).
    const int targetChannel = 5;
    bool cycleChannels = (targetChannel < 0);

    // ── Output file ───────────────────────────────────────────────────────
    gSystem->mkdir("../root_files/wm_decoded/template_an", kTRUE);
    TString name = TString::Format("../root_files/wm_decoded/template_an/"
                                   "decoded_analyzewithtemplate_chan%d.root",
                                   targetChannel);
    TFile outfile(name, "RECREATE");
    REQUIRE(outfile.IsOpen());

    // ── Statistics histograms ─────────────────────────────────────────────
    // These replace the truth-residual histograms in
    // test_analyzewithtemplate.cxx since real data has no truth information.
    TH1F h_nPulses("h_nPulses", "Pulses Found per Waveform;N_{pulses};Counts",
                   13, -0.5, 12.5);
    TH1F h_amplitude(
        "h_amplitude",
        "Fitted Amplitude (first valid pulse);A [ADC counts];Counts", 100, 0,
        500);
    TH1F h_pedestal("h_pedestal",
                    "Pedestal per Waveform;Pedestal [ADC counts];Counts", 100,
                    0, 500);
    TH1F h_chi2ndf("h_chi2ndf",
                   "chi^{2}/NDF (first valid pulse);#chi^{2}/NDF;Counts", 100,
                   0, 20);
    TH1F h_integral(
        "h_integral",
        "Fitted Integral (first valid pulse);Integral [ADC counts];Counts", 100,
        0, 5000);

    // ── Event loop ────────────────────────────────────────────────────────
    long long nentries = T->GetEntries();
    int wfCount = 0;   // waveforms processed so far
    int zeroCount = 0; // events where no pulse was found

    for (long long iev = 0; iev < nentries && wfCount < num_wfs; ++iev) {
      T->GetEntry(iev);

      if (nsamps_idx <= 0)
        continue; // skip empty events

      int nch = nsamps_idx;
      int nsamps_per_ch = (nch > 0) ? (nsamps / nch) : 0;
      if (nsamps_per_ch <= 0)
        continue;

      int ch_start = cycleChannels ? 0 : targetChannel;
      int ch_end = cycleChannels ? (nch - 1) : targetChannel;

      if (ch_start >= nch) {
        printf("Warning: channel %d not present in event %lld (nch=%d). "
               "Skipping.\n",
               targetChannel, iev, nch);
        continue;
      }
      ch_end = TMath::Min(ch_end, nch - 1);

      for (int ci = ch_start; ci <= ch_end && wfCount < num_wfs; ++ci) {
        int startIdx = static_cast<int>(samps_idx[ci]);
        int elemID = static_cast<int>(samps_elemID[ci]);

        // ── Build sample vector ───────────────────────────────────────
        std::vector<Double_t> samples(nsamps_per_ch);
        for (int s = 0; s < nsamps_per_ch; ++s)
          samples[s] = samps[startIdx + s];

        // ── Run AnalyzeWithTemplate ───────────────────────────────────
        tmplCfg.fMFEnd = static_cast<Int_t>(samples.size()) - 1;
        FADCPulseWF templateResult =
            analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

        // ── Fill statistics histograms ────────────────────────────────
        h_nPulses.Fill(templateResult.nPulses);
        h_pedestal.Fill(templateResult.pedestal);

        if (templateResult.nPulses == 0) {
          ++zeroCount;
        }

        for (const auto &p : templateResult.pulses) {
          if (!p.valid)
            continue;
          // Only the first valid pulse goes into the per-waveform hists
          h_amplitude.Fill(p.amplitude);
          h_chi2ndf.Fill(p.chi2ndf);
          h_integral.Fill(p.integral);
          break;
        }

        // ── Build waveform + fit overlay (only when at least one pulse found)
        // ──
        if (templateResult.nPulses > 0) {
          Int_t N = static_cast<Int_t>(samples.size());

          TGraph *gSamples = new TGraph(N);
          gSamples->SetName(Form("wf_raw_ev%lld_ch%d", iev, elemID));
          for (Int_t j = 0; j < N; ++j)
            gSamples->SetPoint(j, (double)j, samples[j]);
          gSamples->SetMarkerStyle(20);
          gSamples->SetMarkerSize(0.8);
          gSamples->SetMarkerColor(kBlue + 1);
          gSamples->SetLineColor(kBlue + 1);

          // Build fit graph using stage-aware model evaluation
          auto evalPulseModel = [&](Double_t t,
                                    const FitPulseResult &pulse) -> Double_t {
            if (!pulse.valid)
              return 0.0;
            if (pulse.fitStage == 1 || pulse.fittedSigma == 0.0) {
              return pulse.amplitude * tmpl.eval(t - pulse.timeSample);
            } else {
              Double_t dt = t - pulse.timeSample;
              Double_t sigma = pulse.fittedSigma;
              Double_t tau = pulse.fittedTau;
              if (sigma <= 0.0 || tau <= 0.0)
                return 0.0;
              Double_t s2_2t2 = (sigma * sigma) / (2.0 * tau * tau);
              Double_t erf_arg =
                  (dt - (sigma * sigma) / tau) / (sigma * std::sqrt(2.0));
              if (erf_arg > 700.0)
                return 0.0;
              Double_t f =
                  std::exp(s2_2t2 - (dt / tau)) * (1.0 + std::erf(erf_arg));
              Double_t dt_peak = (sigma * sigma) / tau;
              Double_t f_peak =
                  std::exp(s2_2t2 - (dt_peak / tau)) * (1.0 + std::erf(0.0));
              return (f_peak > 1e-12 && f > 1e-12)
                         ? (pulse.amplitude * f / f_peak)
                         : 0.0;
            }
          };

          // Build fit graph at 0.25-sample steps
          const double fitStep = 0.25;
          const double fitXMin = 0.0;
          const double fitXMax = static_cast<double>(N - 1);
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
          gFit->SetName(Form("fit_ev%lld_ch%d", iev, elemID));
          gFit->SetLineColor(kRed);
          gFit->SetLineWidth(2);

          TMultiGraph *mg = new TMultiGraph();
          mg->SetName(Form("Waveform_ev%lld_ch%d_%dpulse(s)", iev, elemID,
                           templateResult.nPulses));
          mg->SetTitle(Form("Template Fit  ev %lld  ch %d;Sample;ADC counts",
                            iev, elemID));
          mg->Add(gSamples, "P");
          mg->Add(gFit, "L");
          mg->Write();
        }

        ++wfCount;
      }

      // In cycle-all-channels mode, only use the first non-empty event
      if (cycleChannels && wfCount > 0)
        break;
    }

    // ── Write histograms ──────────────────────────────────────────────────
    h_nPulses.Write();
    h_amplitude.Write();
    h_pedestal.Write();
    h_chi2ndf.Write();
    h_integral.Write();

    outfile.Close();

    // ── Summary printout (mirrors test_analyzewithtemplate.cxx style) ─────
    double frac_with_pulse =
        (wfCount > 0) ? (1.0 - static_cast<double>(zeroCount) / wfCount) * 100.0
                      : 0.0;

    std::cout << "\n======================================================\n";
    std::cout << "       DECODED DATA TEMPLATE ANALYSIS STATISTICS      \n";
    std::cout << "======================================================\n";
    std::cout << "Total Waveforms Processed : " << wfCount << "\n";
    std::cout << "Waveforms With >= 1 Pulse : " << (wfCount - zeroCount)
              << "\n";
    std::cout << "Waveforms With No Pulse   : " << zeroCount << "\n";
    std::cout << "Pulse-Finding Efficiency  : " << std::fixed
              << std::setprecision(2) << frac_with_pulse << "%\n";
    std::cout << "Output File               : " << name.Data() << "\n";
    std::cout << "======================================================\n\n";

    // Basic sanity check: we should have processed at least one waveform
    REQUIRE(wfCount > 0);
  }

  inFile->Close();
}
