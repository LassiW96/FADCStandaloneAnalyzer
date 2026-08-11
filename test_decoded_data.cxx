/*
 * test_decoded_data.cxx
 *
 * Similar to test_temp.cxx but reads real decoded FADC waveforms from
 *   /MyFiles/StandAloneAnalyzer/wm_data/decoded_data/
 * instead of using the waveformGenerator.
 *
 * The waveform data layout follows plot_waveforms.C:
 *   - TTree "T" inside the ROOT file
 *   - Branches used:
 *       Ndata.moller.scint.samps        (Int_t)   – total number of samples in event
 *       moller.scint.samps              (Double_t[]) – flat sample array
 *       Ndata.moller.scint.samps_idx    (Int_t)   – number of channels in event
 *       moller.scint.samps_idx          (Double_t[]) – per-channel start index into samps[]
 *       Ndata.moller.scint.samps_elemID (Int_t)   – number of channels with element IDs
 *       moller.scint.samps_elemID       (Double_t[]) – hardware channel (element) ID
 *
 * Compile & run (from the test1/ directory):
 *   g++ -std=c++17 test_main.cxx test_decoded_data.cxx \
 *       FADCStandaloneAnalyzer.cxx waveformGenerator.cxx \
 *       -lCatch2Main -lCatch2 $(root-config --cflags --libs) \
 *       -o run_tests_decoded && ./run_tests_decoded
 *
 * Output ROOT file: decoded_fit_debug.root
 *   Structure mirrors template_fit_debug.root produced by test_temp.cxx.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "FADCStandaloneAnalyzer.h"
#include "RtypesCore.h"
#include "waveformGenerator.h"   // for WaveformParams (used to build the template)

// ROOT Headers
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>
#include <TMultiGraph.h>
#include <TMath.h>

// ── Decoded data file ─────────────────────────────────────────────────────────
static const char* kInputFile =
    "/home/lasitha/MyFiles/StandAloneAnalyzer/wm_data/decoded_data/"
    "moller_fadc_replayed_329_seg0_1_08_10_26.root";

// Maximum array sizes (must be at least as large as the decoder fills)
static const int kMaxSamps    = 20000;
static const int kMaxChans    = 1000;

// ── Utility: make a TGraph from a std::vector<Double_t> ──────────────────────
static TGraph* makeGraph(const std::vector<Double_t>& y,
                         const std::string& name,
                         const std::string& title,
                         int color       = kBlue,
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

// ── Utility: write a full diagnostic sub-directory for one waveform ───────────
/*
 * Mirrors writeWaveformDiagnostic() in test_temp.cxx.
 * Adds an "event/channel" label to the directory name so it is easy to
 * identify which decoded event / hardware channel it came from.
 */
static void writeWaveformDiagnostic(TFile&                         outfile,
                                    int                            runIdx,
                                    long long                      eventNum,
                                    int                            elemID,
                                    const std::vector<Double_t>&   samples,
                                    const std::vector<Double_t>&   wfPS,
                                    const std::vector<Double_t>&   mfOut,
                                    const FADCTemplate&            tmpl,
                                    const FADCPulseWF&             result)
{
    // Directory name encodes the running index, event number, and channel ID
    TString dirName = Form("wf_%03d_ev%lld_ch%d", runIdx, eventNum, elemID);
    outfile.mkdir(dirName.Data());
    outfile.cd(dirName.Data());

    // ── template graphs (written once per run for easy browsing) ─────────────
    TGraph* gShape  = makeGraph(tmpl.shape,    "shape",    "Template Shape",       kBlack,   20);
    TGraph* gKernel = makeGraph(tmpl.mfKernel, "mfKernel","MF Kernel (reversed)", kGreen+2, 21);

    // ── raw samples ──────────────────────────────────────────────────────────
    int N = static_cast<int>(samples.size());
    TGraph* gSamples = new TGraph(N);
    gSamples->SetName("wf_raw");
    gSamples->SetTitle(Form("Raw ADC  ev%lld ch%d;Sample;ADC counts", eventNum, elemID));
    for (int j = 0; j < N; ++j) gSamples->SetPoint(j, j, samples[j]);
    gSamples->SetMarkerStyle(3);
    gSamples->SetMarkerColor(kBlue);
    gSamples->SetLineColor(kBlue);

    // ── ped-subtracted & matched-filter output ───────────────────────────────
    TGraph* gWfPS  = makeGraph(wfPS,  "wfPS",  "Ped-Subtracted Waveform", kCyan+1, 24);
    TGraph* gMFOut = makeGraph(mfOut, "mfOut", "Matched Filter Output",   kRed,    25);

    // ── fitted-pulse overlay sampled at 0.25-step intervals ─────────────────
    TF1 fTotal("fTotal",
        [&tmpl, &result](double* x, double* /*p*/) {
            double val = result.pedestal;
            for (const auto& p : result.pulses)
                if (p.valid)
                    val += p.amplitude * tmpl.eval(x[0] - p.timeSample);
            return val;
        }, -0.5, N - 0.5, 0);
    const double fitStep = 0.25;
    const double fitXMin = -0.5;
    const double fitXMax = N - 0.5;
    int nFitPts = static_cast<int>((fitXMax - fitXMin) / fitStep) + 1;
    TGraph* gFit = new TGraph(nFitPts);
    for (int k = 0; k < nFitPts; ++k) {
        double xv = fitXMin + k * fitStep;
        gFit->SetPoint(k, xv, fTotal.Eval(xv));
    }
    gFit->SetName("fit_overlay");
    gFit->SetTitle("Fit Overlay;Sample;ADC counts");
    gFit->SetLineColor(kRed+1);
    gFit->SetLineWidth(2);

    // ── Write all ─────────────────────────────────────────────────────────────
    gShape->Write();
    gKernel->Write();
    gSamples->Write();
    gWfPS->Write();
    gMFOut->Write();
    gFit->Write();

    outfile.cd();   // return to file root
}

// ═══════════════════════════════════════════════════════════════════════════════
// Catch2 test case
// ═══════════════════════════════════════════════════════════════════════════════

TEST_CASE("Analyze Decoded FADC Waveforms", "[DecodedData][Plotting]") {

    // ── 1.  Open the decoded ROOT file ───────────────────────────────────────
    TFile* inFile = TFile::Open(kInputFile);
    REQUIRE(inFile != nullptr);
    REQUIRE(!inFile->IsZombie());

    TTree* T = dynamic_cast<TTree*>(inFile->Get("T"));
    REQUIRE(T != nullptr);

    // ── 2.  Set up branch reading (mirrors plot_waveforms.C) ─────────────────
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("Ndata.moller.scint.samps",        1);
    T->SetBranchStatus("moller.scint.samps",               1);
    T->SetBranchStatus("Ndata.moller.scint.samps_idx",    1);
    T->SetBranchStatus("moller.scint.samps_idx",           1);
    T->SetBranchStatus("Ndata.moller.scint.samps_elemID", 1);
    T->SetBranchStatus("moller.scint.samps_elemID",        1);

    Int_t    nsamps        = 0;
    Int_t    nsamps_idx    = 0;
    Int_t    nsamps_elemID = 0;
    Double_t samps        [kMaxSamps];
    Double_t samps_idx    [kMaxChans];
    Double_t samps_elemID [kMaxChans];

    T->SetBranchAddress("Ndata.moller.scint.samps",        &nsamps);
    T->SetBranchAddress("moller.scint.samps",               samps);
    T->SetBranchAddress("Ndata.moller.scint.samps_idx",    &nsamps_idx);
    T->SetBranchAddress("moller.scint.samps_idx",           samps_idx);
    T->SetBranchAddress("Ndata.moller.scint.samps_elemID", &nsamps_elemID);
    T->SetBranchAddress("moller.scint.samps_elemID",        samps_elemID);

    // ── 3.  Build the template from bi-exponential parameters ────────────────
    //        Adjust tau_r / tau_f to match the real detector response if known.
    WaveformParams wf;
    wf.tau_r     = 3.0;
    wf.tau_f     = 10.0;
    wf.mean_gain = 10.0;

    FADCTemplate tmpl = FADCTemplate::fromWaveformParams(wf);

    // ── 4.  Configure the template analysis ──────────────────────────────────
    TemplateAnalysisConfig tmplCfg;
    tmplCfg.fMFSeparation = 8;
    tmplCfg.fScanHalf     = 10;
    tmplCfg.fMFLeft       = 2;
    tmplCfg.fMFRight      = 2;
    tmplCfg.fAmpCutADC    = 4.0;

    FADCStandaloneAnalyzer analyzer;

    SECTION("Write decoded waveform diagnostic ROOT file") {

        TFile outfile("decoded_fit_debug.root", "RECREATE");
        REQUIRE(outfile.IsOpen());

        // Write the template graphs to the top level for reference
        TGraph* gShape  = makeGraph(tmpl.shape,    "template_shape",  "Template Shape",       kBlack,   20);
        TGraph* gKernel = makeGraph(tmpl.mfKernel, "template_kernel", "MF Kernel (reversed)", kGreen+2, 21);
        gShape->Write();
        gKernel->Write();

        // ── 5.  Loop over events and channels ────────────────────────────────
        //        Collect up to N_WF waveforms from channel index 0 (targetChannel).
        //        Set targetChannel = -1 to cycle through all channels in the
        //        first non-empty event (same semantics as plot_waveforms.C).
        const int  N_WF          = 20;   // number of waveforms to process
        const int  targetChannel = 3;    // channel index to follow (-1 = all)
        bool cycleChannels       = (targetChannel < 0);

        long long nentries = T->GetEntries();
        int wfCount = 0;

        for (long long iev = 0; iev < nentries && wfCount < N_WF; ++iev) {
            T->GetEntry(iev);

            if (nsamps_idx <= 0) continue;          // skip empty events

            int nch           = nsamps_idx;
            int nsamps_per_ch = (nch > 0) ? (nsamps / nch) : 0;
            if (nsamps_per_ch <= 0) continue;

            int ch_start = cycleChannels ? 0           : targetChannel;
            int ch_end   = cycleChannels ? (nch - 1)   : targetChannel;

            if (ch_start >= nch) {
                printf("Warning: channel %d not present in event %lld (nch=%d). Skipping.\n",
                       targetChannel, iev, nch);
                continue;
            }
            ch_end = TMath::Min(ch_end, nch - 1);

            for (int ci = ch_start; ci <= ch_end && wfCount < N_WF; ++ci) {
                int startIdx = static_cast<int>(samps_idx[ci]);
                int elemID   = static_cast<int>(samps_elemID[ci]);

                // ── Build sample vector for this channel ──────────────────────
                std::vector<Double_t> samples(nsamps_per_ch);
                for (int s = 0; s < nsamps_per_ch; ++s)
                    samples[s] = samps[startIdx + s];

                // ── Pedestal subtraction ──────────────────────────────────────
                int    nPed = std::min(tmplCfg.fNPedSamples, (int)samples.size());
                double ped  = 0.0;
                for (int j = 0; j < nPed; ++j) ped += samples[j];
                ped /= nPed;

                std::vector<Double_t> wfPS(samples.size());
                for (size_t j = 0; j < samples.size(); ++j)
                    wfPS[j] = samples[j] - ped;

                // ── Matched filter ────────────────────────────────────────────
                std::vector<Double_t> mfOut =
                    FADCStandaloneAnalyzer::computeMatchedFilter(wfPS, tmpl, tmplCfg);

                // ── Full template fit ─────────────────────────────────────────
                FADCPulseWF result = analyzer.AnalyzeWithTemplate(samples, tmpl, tmplCfg);

                // ── Console summary ───────────────────────────────────────────
                printf("wf %3d | ev %6lld | ch %3d | ped %.1f | nPulses %d\n",
                       wfCount, iev, elemID, result.pedestal, result.nPulses);
                for (int pi = 0; pi < (int)result.pulses.size(); ++pi) {
                    const auto& p = result.pulses[pi];
                    if (!p.valid) continue;
                    printf("         pulse %d: t0=%.2f  A=%.2f  chi2/ndf=%.2f\n",
                           pi, p.timeSample, p.amplitude, p.chi2ndf);
                }

                // ── Write diagnostic graphs to output file ────────────────────
                writeWaveformDiagnostic(outfile, wfCount, iev, elemID,
                                        samples, wfPS, mfOut, tmpl, result);
                ++wfCount;
            }

            // In cycle-all-channels mode, only process the first non-empty event
            if (cycleChannels && wfCount > 0) break;
        }

        printf("\nProcessed %d waveform(s) from:\n  %s\n", wfCount, kInputFile);

        outfile.Close();

        // Sanity-check: at least one waveform was processed
        REQUIRE(wfCount > 0);
    }

    inFile->Close();
}
