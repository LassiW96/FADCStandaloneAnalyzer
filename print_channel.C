#include <iostream>
#include <iomanip>
#include <vector>
#include <TFile.h>
#include <TTree.h>

void print_channel(long long targetEvent = 11, int targetChannel = 0,
                   const char* filename = "/home/lasitha/MyFiles/StandAloneAnalyzer/wm_data/decoded_data/moller_fadc_replayed_329_seg0_1_08_10_26.root") {
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file: " << filename << std::endl;
        return;
    }

    TTree *T = dynamic_cast<TTree*>(f->Get("T"));
    if (!T) {
        std::cerr << "Error: Tree 'T' not found in file." << std::endl;
        f->Close();
        return;
    }

    const int kMaxSamps = 20000;
    const int kMaxChans = 1000;

    Int_t nsamps = 0;
    Int_t nsamps_idx = 0;
    Int_t nsamps_elemID = 0;
    Double_t samps[kMaxSamps];
    Double_t samps_idx[kMaxChans];
    Double_t samps_elemID[kMaxChans];

    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("Ndata.moller.scint.samps", 1);
    T->SetBranchStatus("moller.scint.samps", 1);
    T->SetBranchStatus("Ndata.moller.scint.samps_idx", 1);
    T->SetBranchStatus("moller.scint.samps_idx", 1);
    T->SetBranchStatus("Ndata.moller.scint.samps_elemID", 1);
    T->SetBranchStatus("moller.scint.samps_elemID", 1);

    T->SetBranchAddress("Ndata.moller.scint.samps", &nsamps);
    T->SetBranchAddress("moller.scint.samps", samps);
    T->SetBranchAddress("Ndata.moller.scint.samps_idx", &nsamps_idx);
    T->SetBranchAddress("moller.scint.samps_idx", samps_idx);
    T->SetBranchAddress("Ndata.moller.scint.samps_elemID", &nsamps_elemID);
    T->SetBranchAddress("moller.scint.samps_elemID", samps_elemID);

    if (targetEvent < 0 || targetEvent >= T->GetEntries()) {
        std::cerr << "Error: Event " << targetEvent << " is out of range [0, " << T->GetEntries() - 1 << "]" << std::endl;
        f->Close();
        return;
    }

    T->GetEntry(targetEvent);

    int nch = nsamps_idx;
    if (nch <= 0) {
        std::cout << "Event " << targetEvent << " contains no channel data." << std::endl;
        f->Close();
        return;
    }

    int nsamps_per_ch = nsamps / nch;
    int foundIdx = -1;

    for (int ci = 0; ci < nch; ++ci) {
        if (static_cast<int>(samps_elemID[ci]) == targetChannel) {
            foundIdx = ci;
            break;
        }
    }

    if (foundIdx < 0) {
        std::cout << "Channel " << targetChannel << " not found in event " << targetEvent << "." << std::endl;
        f->Close();
        return;
    }

    int startIdx = static_cast<int>(samps_idx[foundIdx]);
    std::cout << "========================================================\n";
    std::cout << "File:            " << filename << "\n";
    std::cout << "Event:           " << targetEvent << "\n";
    std::cout << "Channel (elemID):" << targetChannel << "\n";
    std::cout << "Samples/Channel: " << nsamps_per_ch << " (start index in flat buffer: " << startIdx << ")\n";
    std::cout << "========================================================\n";
    for (int s = 0; s < nsamps_per_ch; ++s) {
        std::cout << "  sample[" << std::setw(2) << s << "] = " << std::fixed << std::setprecision(5) << samps[startIdx + s] << "\n";
    }

    f->Close();
}
