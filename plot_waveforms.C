/*
 * plot_waveforms.C
 *
 * Plots the first 10 waveforms from the moller.scint.samps branch
 * in the decoded FADC replay output.
 *
 * Usage (from ROOT prompt or shell):
 *   root -l plot_waveforms.C
 *   root -l 'plot_waveforms.C("/path/to/file.root", 10, 0)'
 *
 * Arguments:
 *   fname      - path to the decoded ROOT file
 *   nwaveforms - number of waveforms to plot (default: 10)
 *   channel    - which channel (elemID index) to plot, 0-based (default: 0)
 *                pass -1 to cycle through all available channels in event 0
 *   startEvent - first event to start from (default: 0)
 */

void plot_waveforms(
    const char* fname      = "/home/lasitha/MyFiles/StandAloneAnalyzer/wm_data/decoded_data/moller_fadc_replayed_329_seg0_1_08_10_26.root",
    int nwaveforms         = 10,
    int channel            = 0,
    int startEvent         = 0
)
{
    // --- Open file and tree ---------------------------------------------------
    TFile *f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        printf("Error: cannot open file %s\n", fname);
        return;
    }

    TTree *T = (TTree*)f->Get("T");
    if (!T) {
        printf("Error: TTree 'T' not found in %s\n", fname);
        return;
    }

    // Disable all branches; enable only what we need to avoid missing-dict errors
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("Ndata.moller.scint.samps",    1);
    T->SetBranchStatus("moller.scint.samps",           1);
    T->SetBranchStatus("Ndata.moller.scint.samps_idx", 1);
    T->SetBranchStatus("moller.scint.samps_idx",       1);
    T->SetBranchStatus("Ndata.moller.scint.samps_elemID", 1);
    T->SetBranchStatus("moller.scint.samps_elemID",    1);

    // --- Branch variables -----------------------------------------------------
    const int MAXSAMPS  = 20000;
    const int MAXCHANS  = 1000;

    Int_t    nsamps, nsamps_idx, nsamps_elemID;
    Double_t samps[MAXSAMPS];
    Double_t samps_idx[MAXCHANS];
    Double_t samps_elemID[MAXCHANS];

    T->SetBranchAddress("Ndata.moller.scint.samps",       &nsamps);
    T->SetBranchAddress("moller.scint.samps",              samps);
    T->SetBranchAddress("Ndata.moller.scint.samps_idx",   &nsamps_idx);
    T->SetBranchAddress("moller.scint.samps_idx",          samps_idx);
    T->SetBranchAddress("Ndata.moller.scint.samps_elemID",&nsamps_elemID);
    T->SetBranchAddress("moller.scint.samps_elemID",       samps_elemID);

    // --- Canvas setup ---------------------------------------------------------
    // Arrange nwaveforms in a roughly square grid
    int ncols = (int)ceil(sqrt((double)nwaveforms));
    int nrows = (int)ceil((double)nwaveforms / ncols);

    TCanvas *c = new TCanvas("c_waveforms", "FADC Waveform Samples", 300*ncols, 250*nrows);
    c->Divide(ncols, nrows);

    // --- Plot loop ------------------------------------------------------------
    int plotted = 0;
    long nentries = T->GetEntries();

    // channel == -1: show all available channels from the first non-empty event
    // channel >= 0 : show that specific channel index across successive events
    bool cycleChannels = (channel < 0);

    for (long iev = startEvent; iev < nentries && plotted < nwaveforms; iev++) {
        T->GetEntry(iev);

        if (nsamps_idx <= 0) continue;   // skip empty events

        // Number of samples per channel (uniform in this decoder)
        int nch           = nsamps_idx;
        int nsamps_per_ch = (nch > 0) ? (nsamps / nch) : 0;

        if (nsamps_per_ch <= 0) continue;

        // Determine which channel indices to draw from this event
        int ch_start = cycleChannels ? 0          : channel;
        int ch_end   = cycleChannels ? (nch - 1)  : channel;

        if (ch_start >= nch) {
            printf("Warning: channel %d does not exist in event %lld (nch=%d). Skipping.\n",
                   channel, (long long)iev, nch);
            continue;
        }
        ch_end = TMath::Min(ch_end, nch - 1);

        for (int ci = ch_start; ci <= ch_end && plotted < nwaveforms; ci++) {
            int start  = (int)samps_idx[ci];
            int elemID = (int)samps_elemID[ci];

            // Build x (sample number) and y (ADC counts) arrays
            std::vector<double> x(nsamps_per_ch), y(nsamps_per_ch);
            for (int s = 0; s < nsamps_per_ch; s++) {
                x[s] = s;
                y[s] = samps[start + s];
            }

            // Draw
            c->cd(plotted + 1);
            gPad->SetLeftMargin(0.15);
            gPad->SetBottomMargin(0.15);

            TGraph *gr = new TGraph(nsamps_per_ch, x.data(), y.data());
            gr->SetTitle(Form("Ev %lld  ch %d;Sample;ADC counts", (long long)iev, elemID));
            gr->SetLineColor(kBlue + 1);
            gr->SetLineWidth(2);
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.5);
            gr->SetMarkerColor(kBlue + 1);
            gr->Draw("APL");

            plotted++;
        }

        // In single-channel mode each event gives one plot; in cycle mode
        // we only need the first non-empty event.
        if (cycleChannels && plotted > 0) break;
    }

    if (plotted == 0) {
        printf("Warning: no waveforms found in the requested range.\n");
        return;
    }

    c->Update();

    printf("Plotted %d waveform(s).\n", plotted);
    if (cycleChannels)
        printf("Mode: all channels from the first non-empty event (startEvent=%d).\n", startEvent);
    else
        printf("Mode: channel %d across events starting from event %d.\n", channel, startEvent);
    printf("Samples per channel: %d\n",
           (nsamps_idx > 0) ? (nsamps / nsamps_idx) : 0);
}
