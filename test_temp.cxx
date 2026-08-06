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

TGraph* makeGraph(const std::vector<Double_t>& y,
                  const std::string& name,
                  const std::string& title,
                  int color = kBlue,
                  int markerStyle = 20)
{
    int N = static_cast<int>(y.size());
    TGraph* g = new TGraph(N);
    g->SetName(name.c_str());
    g->SetTitle((title + ";Sample;Value").c_str());
    for (int i = 0; i < N; ++i) g->SetPoint(i, i, y[i]);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(markerStyle);
    return g;
}

void writeWaveformDiagnostic(TFile& outfile,
                             int runIdx,
                             const std::vector<Double_t>& samples,
                             const std::vector<Double_t>& wfPS,
                             const std::vector<Double_t>& mfOut,
                             const FADCTemplate&           tmpl,
                             const FADCPulseWF&            result)
{
    outfile.mkdir(Form("run_%d", runIdx));
    outfile.cd(Form("run_%d", runIdx));
    // ── template shape (once, but written per-run for easy browsing) ──
    TGraph* gShape = makeGraph(tmpl.shape,   "shape",    "Template Shape",   kBlack, 20);
    TGraph* gKernel = makeGraph(tmpl.mfKernel,"mfKernel","MF Kernel (reversed)", kGreen+2, 21);
    // ── per-waveform curves ───────────────────────────────────────────
    int N = static_cast<int>(samples.size());
    TGraph* gSamples = new TGraph(N);
    gSamples->SetName("wf_raw");
    for (int j = 0; j < N; ++j) gSamples->SetPoint(j, j, samples[j]);
    gSamples->SetMarkerStyle(3);
    gSamples->SetMarkerColor(kBlue);
    TGraph* gWfPS  = makeGraph(wfPS,  "wfPS",  "Ped-Subtracted Waveform", kCyan+1, 24);
    TGraph* gMFOut = makeGraph(mfOut, "mfOut", "Matched Filter Output",   kRed,    25);
    // ── fitted-pulses overlay (TF1 sampled into a TGraph) ───────────
    TF1 fTotal("fTotal",
        [&tmpl, &result](double* x, double* /*p*/) {
            double val = result.pedestal;
            for (const auto& p : result.pulses)
                if (p.valid)
                    val += p.amplitude * tmpl.eval(x[0] - p.timeSample);
            return val;
        }, -0.5, N - 0.5, 0);
    fTotal.SetNpx(1000);
    TGraph* gFit = new TGraph(&fTotal);
    gFit->SetName("fit_overlay");
    gFit->SetLineColor(kRed+1);
    gFit->SetLineWidth(2);
    // ── Write all objects ─────────────────────────────────────────────
    gShape->Write();
    gKernel->Write();
    gSamples->Write();
    gWfPS->Write();
    gMFOut->Write();
    gFit->Write();
    outfile.cd();   // back to root of file
}

TEST_CASE("Generate Waveform Fit Plot", "[Plotting]") {
    // ── Setup ─────────────────────────────────────────────────────────
    waveformGenerator gen;
    DigitizerParams daq;  daq.noise_sigma = 1.0;  daq.num_samples = 120;
    WaveformParams  wf;   wf.tau_r = 2;  wf.tau_f = 7;  wf.mean_gain = 3.0;
    gen.setDigitizerParams(daq);
    gen.setDetectorParams(wf);
    FADCTemplate tmpl = FADCTemplate::fromWaveformParams(wf);
    TemplateAnalysisConfig tmplCfg;
    tmplCfg.fMFSeparation = 8;
    tmplCfg.fScanHalf     = 10;
    tmplCfg.fMFLeft       = 8;
    tmplCfg.fMFRight      = 8;
    tmplCfg.fAmpCutADC    = 4.0;
    FADCStandaloneAnalyzer analyzer;
    SECTION("Write template diagnostic ROOT file") {
        TFile outfile("template_fit_debug.root", "RECREATE");
        REQUIRE(outfile.IsOpen());
        // ── Write template graphs to top-level directory ───────────────
        TGraph* gShape  = makeGraph(tmpl.shape,     "template_shape",  "Template Shape",        kBlack, 20);
        TGraph* gKernel = makeGraph(tmpl.mfKernel,  "template_kernel", "MF Kernel (reversed)",  kGreen+2, 21);
        gShape->Write();
        gKernel->Write();
        // ── Per-waveform loop ─────────────────────────────────────────
        const int N_WF = 20;   // small number for quick inspection
        for (int i = 0; i < N_WF; ++i) {
            std::vector<Int_t>    raw = gen.generateWaveform(PulseShape::SPE);
            std::vector<Double_t> samples(raw.begin(), raw.end());
            // Reconstruct wfPS (ped subtraction — mirrors AnalyzeWithTemplate)
            int    nPed = std::min(tmplCfg.fNPedSamples, (int)samples.size());
            double ped  = 0.0;
            for (int j = 0; j < nPed; ++j) ped += samples[j];
            ped /= nPed;
            std::vector<Double_t> wfPS(samples.size());
            for (size_t j = 0; j < samples.size(); ++j) wfPS[j] = samples[j] - ped;
            // Call the now-public ComputeMatchedFilter
            std::vector<Double_t> mfOut =
                FADCStandaloneAnalyzer::computeMatchedFilter(wfPS, tmpl, tmplCfg);
            FADCPulseWF result = analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);
            writeWaveformDiagnostic(outfile, i, samples, wfPS, mfOut, tmpl, result);
        }
        outfile.Close();
    }
}