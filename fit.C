// Script wrote, tested and maintaned by Wassim Hamdi and Malek Mazouz

#include "TCanvas.h"
 #include "TROOT.h"
 #include "TGraphErrors.h"
 #include "TF1.h"
 #include "TH1F.h"
 #include "TH2F.h"
 #include "TLegend.h"
 #include "TMath.h"
 #include "TLorentzVector.h"
 #include "TChain.h"
 #include "TTree.h"
 #include "TStyle.h"
 #include "TSystem.h"
 #include "TRandom.h"
 #include "TFile.h"
 #include "TROOT.h"
 #include "TMatrixD.h"
 #include "TLatex.h"
 #include "TLine.h"
 #include "TVirtualFFT.h"
 #include "TStopwatch.h"
 #include "TSpectrum.h"
 #include "TObjectTable.h"
 #include "TApplication.h"
 #include <Math/Interpolator.h>
 #include <iostream>
 #include <fstream>

const int ntime=110; // number of time samples for each fADC channel (100 for runs 55 and 56)
const int nfADC=7; //number of fADC modules (16 channels each) used in this run
const int nchannel=16; //number of channels in each fADC module
const int ncol=30; //number of calo colomns tested in this run (3 for runs 55 and 56)
const int nlin=36; //number of calo blocks in each colomn
const int nblocks=ncol*nlin; //number of tested calo blocks
const int maxpulses=4; //nb maximal de pulses dans la fentre ntime 
const int nsignals=nblocks*maxpulses;
Double_t dt=4.; //time bin (sample) width (4 ns), the total time window is then ntime*dt
const int ntemp=56; //nb of temp sensors

ROOT::Math::Interpolator *interpolation[nblocks];
TH1F *hsig_i[nblocks];//histogram showing the waveform for block i and a given event
Int_t wfnpulse[nblocks],wfnpulsei[nblocks],wfnpulsesup[nblocks];
const int maxwfpulses=12; //nb maximal de pulses que la wfa peut trouver
Double_t wfampl[nblocks][maxwfpulses]; 
Double_t wftime[nblocks][maxwfpulses]; 
Int_t wfnpulse_sev[nblocks][maxwfpulses];
TF1 *finter[nblocks]; 
Double_t timeref[nblocks];
Double_t calo_dist=9.5;// calo distance to the target (m)
Double_t timerefacc=(calo_dist-9.5)/(3.e8*1.e-9*4);//time offset (4ns units) due to the change of calo_dist
Int_t cpulse,highnbpulse;
Double_t corr_time_HMS;
Float_t cortime[nblocks];

//Matched Filter parameters
const int mfleft=5; // start bin of the refwf used as a filter = max of the refwf - mfleft
const int mfright=5; // end bin of the refwf used as a filter = max of the refwf + mfright
const int mfwidth=mfleft+mfright+1; //width of the refwf used as a filter (4ns units)
const int mfstart=10; // Search for peaks in the wf between bins mfstart and mfend (4ns units)
const int mfend=100; // Search for peaks in the wf between bins mfstart and mfend (4ns units)
const double specthres=0.01; //peaks with amplitude less than specthres*highest_peak are discarded when TSpectrum::Search() is called
const double mfthres=0.5; //peaks with amplitude less than mfthres in the mf of the wf are discarded 
const double trig_thres=10; //threshold (mV) on sum of 3x3 wf to allow the fit of the central wf 
const int coinc_width=20; // (4ns units) the trig_thres will be applied on sum of 3x3 wf in the region {expected coinc time +/- coinc_width} instead of all wf range 
Double_t mfyref[nblocks][mfwidth]; // part of the refwf that will be used in the MF
Double_t mfint[nblocks];
Float_t mfmin;
Int_t peaksup;
TH1F *hsig_mf[nblocks];//histogram showing the matched filter of the waveform for block i and a given event
TH1F *hsig_trig[nblocks];//histogram showing the sum of 3x3 waveforms around block i and a given event (will be used as a trigger for the wf fit)
Int_t trig[nblocks]; //must be equal 1 to allow the fit of the wf of block i
TSpectrum *spec;
Int_t npeaks;
Double_t minsignal[nblocks];
Int_t pres[nblocks];
Double_t dum,minwf,maxwf;
int end_event, start_event;

const int evts_per_job = 100000;

//////////////////////////////////////////////////////////////////////////////
// La fonction du fit de la wf
Double_t func(Double_t *x, Double_t *par){ 
 Int_t j= (int)(par[0]);
 Double_t val=0;
 for(Int_t p=0;p<maxwfpulses;p++){
 if(x[0]-par[2+2*p]>1&&x[0]-par[2+2*p]<109)val+=par[3+2*p]*interpolation[j]->Eval(x[0]-par[2+2*p]);
 }
 return val+par[1];
}

//////////////////////////////////////////////////////////////////////////////
Int_t voisin8(Int_t i, Int_t j){//neighbour  number j of block i
Int_t voisin=-1;
Int_t ilin=(int)(i/ncol);
Int_t icol=i%ncol;
Int_t col=-100;
Int_t lin=-100;
Int_t good=0;
if(j==0){if(icol<ncol-1){col=icol+1;lin=ilin;good=1;}}
if(j==1){if(icol>0){col=icol-1;lin=ilin;good=1;}}
if(j==2){if(ilin<nlin-1){col=icol;lin=ilin+1;good=1;}}
if(j==3){if(ilin>0){col=icol;lin=ilin-1;good=1;}}
if(j==4){if(icol<ncol-1&&ilin<nlin-1){col=icol+1;lin=ilin+1;good=1;}}
if(j==5){if(icol>0&&ilin<nlin-1){col=icol-1;lin=ilin+1;good=1;}}
if(j==6){if(icol<ncol-1&&ilin>0){col=icol+1;lin=ilin-1;good=1;}}
if(j==7){if(icol>0&&ilin>0){col=icol-1;lin=ilin-1;good=1;}}

if(j>-1&&j<8&&good==1){voisin=col+lin*ncol;}

return voisin;
}

//////////////////////////////////////////////////////////////////////////////
int getMemoryUsageKB() {
    std::ifstream statusFile("/proc/self/status");
    if (!statusFile.is_open()) {
        std::cerr << "Error: Could not open /proc/self/status" << std::endl;
        return -1;
    }

    std::string line;
    while (std::getline(statusFile, line)) {
        size_t pos = line.find("VmRSS:");
        if (pos != std::string::npos) {
            std::stringstream ss(line.substr(pos + 7)); // Skip "VmRSS: "
            long memoryKB;
            if (ss >> memoryKB) {
                std::cout << "VmRSS: " << memoryKB << " KB" << std::endl;
                return 0; // Success
            } else {
                std::cerr << "Error: Could not parse VmRSS value" << std::endl;
                return 1;
            }
        }
    }

    std::cerr << "Error: VmRSS line not found" << std::endl;
    return 1;
}

//////////////////////////////////////////////////////////////////////////////

void npulsewf2(Int_t bn, Int_t minsamp,Int_t maxsamp){

wfnpulse[bn]=0;
for(Int_t p=0;p<maxwfpulses;p++){wfampl[bn][p]=-1;}

//Detect the pulses with TSpectrum::Search() applied to the matched filter transform of the wf  
    npeaks = spec->Search(hsig_mf[bn], 2,"nobackground", specthres);
    for(Int_t ip=0;ip<npeaks;ip++){
    if(spec->GetPositionX()[ip]-2>TMath::Max(mfstart,minsamp)&&spec->GetPositionX()[ip]-2<TMath::Min(mfend,maxsamp)&&spec->GetPositionY()[ip]>mfthres){
    wfampl[bn][wfnpulse[bn]]=TMath::Abs(hsig_i[bn]->GetBinContent((int)(spec->GetPositionX()[ip]-2))-minsignal[bn]);
    wftime[bn][wfnpulse[bn]]=spec->GetPositionX()[ip]-2;
    wfnpulse[bn]++;
    }}
// Check if additional pulses can be identified in the matched filter transform of the wf    
    for(Int_t it=mfstart;it<mfend;it++){if(hsig_mf[bn]->GetBinContent(it+1)>mfthres){    if(hsig_mf[bn]->GetBinContent(it+1)>hsig_mf[bn]->GetBinContent(it)&&hsig_mf[bn]->GetBinContent(it)>hsig_mf[bn]->GetBinContent(it-1)&&hsig_mf[bn]->GetBinContent(it-1)>hsig_mf[bn]->GetBinContent(it-2)/*&&hsig_mf[bn]->GetBinContent(it-2)>hsig_mf[bn]->GetBinContent(it-3)*/){ //4 increasing samples
    if(hsig_mf[bn]->GetBinContent(it+2)<hsig_mf[bn]->GetBinContent(it+1)&&hsig_mf[bn]->GetBinContent(it+3)<hsig_mf[bn]->GetBinContent(it+2)){//2 decreasing samples
    peaksup=1;
    for(Int_t ip=0;ip<npeaks;ip++){if(TMath::Abs(hsig_mf[bn]->GetBinCenter(it+1)-spec->GetPositionX()[ip])<3.1) peaksup=0;}
    if(peaksup==1){
    wfampl[bn][wfnpulse[bn]]=TMath::Abs(hsig_i[bn]->GetBinContent(it+1-2));
    wftime[bn][wfnpulse[bn]]=hsig_mf[bn]->GetBinCenter(it+1)-2;
    wfnpulse[bn]++;
    }}}}}


if(wfnpulse[bn]>maxwfpulses-2){cout<<"Warning : High nb of pulses in wf of block "<<bn<<endl;}
}

//////////////////////////////////////////////////////////////////////////////

void Fitwf2(Int_t bn, Int_t minsamp,Int_t maxsamp){

//Adjust the parameters of the fit function
for(Int_t p=0;p<maxwfpulses;p++){finter[bn]->FixParameter(2+2*p,0.);finter[bn]->FixParameter(3+2*p,0.);}
if(wfnpulse[bn]>0){
for(Int_t p=0;p<TMath::Min(maxwfpulses,wfnpulse[bn]);p++){
finter[bn]->ReleaseParameter(2+2*p);finter[bn]->ReleaseParameter(3+2*p);
finter[bn]->SetParameter(2+2*p,wftime[bn][p]-timeref[bn]/*+(corr_time_HMS-cortime[bn])/dt*/);finter[bn]->SetParameter(3+2*p,wfampl[bn][p]);
finter[bn]->SetParLimits(2+2*p,wftime[bn][p]-timeref[bn]/*+(corr_time_HMS-cortime[bn])/dt*/-4,wftime[bn][p]-timeref[bn]/*+(corr_time_HMS-cortime[bn])/dt*/+4);finter[bn]->SetParLimits(3+2*p,wfampl[bn][p]*0.2,wfampl[bn][p]*5);
}}

finter[bn]->SetParameter(1,0.);finter[bn]->SetParLimits(1,-100,100.);
if(wfnpulse[bn]>0)hsig_i[bn]->Fit(Form("finter_%d",bn),"Q","",minsamp,maxsamp);
}

//////////////////////////////////////////////////////////////////////////////

void Apply_thres(Int_t bn){
trig[bn]=0;
minwf=1.e6;
maxwf=-1.e6;
//compute the sum of 3x3 wf around block bn
for(Int_t it=0;it<ntime;it++){
dum=hsig_i[bn]->GetBinContent(it+1);
for(Int_t j=0;j<8;j++){if(voisin8(bn,j)>-1){if(pres[voisin8(bn,j)]){dum+=hsig_i[voisin8(bn,j)]->GetBinContent(it+1);}}}
hsig_trig[bn]->SetBinContent(it+1,dum);
minwf=TMath::Min(dum,minwf);
if(TMath::Abs(hsig_i[bn]->GetBinCenter(it+1)-(timeref[bn]+timerefacc))<coinc_width)maxwf=TMath::Max(dum,maxwf);
}
if((maxwf-minwf)>trig_thres)trig[bn]=1;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////


void TEST_2(Int_t runnumber, Int_t seg, Int_t job_number){
const int nslots=1104; //nb maximal de slots dans tous les fADC
const int Ndata=nslots*(ntime+1+1); //(1104 slots fADC au total mais pas tous utilises y compris 2 PM scintillateurs?) should correspond to Ndata.NPS.cal.fly.adcSampWaveform variable in the root file
Double_t ADCtomV=1000./4096; //1000 mV correspond a 4096 canaux
Double_t integtopC= (1./4096)*(dt*1.e3/50); // (1 V / 4096 adc channels) * (4000 ps time sample / 50 ohms input resistance) = 0.020 pc/channel 

TStopwatch t;
t.Start();
Double_t calodist=9.5;

//if(runnumber>1500 && runnumber<5000){calodist=9.5;} //calib elastique

if(runnumber>1571 && runnumber<3667){calodist=3.;}
if(runnumber>3666 && runnumber<4632){calodist=4.;}
if(runnumber>4635 && runnumber<4953){calodist=6.;}
if(runnumber>4965 && runnumber<5344){calodist=4.;}
if(runnumber>5354 && runnumber<5464){calodist=3.;}
if(runnumber>5523 && runnumber<7013){calodist=3.5;}

 
calo_dist=calodist;
timerefacc=(calo_dist-9.5)/(3.e8*1.e-9*4);
 
long memoryKB = getMemoryUsageKB();

cout<<" timerefacc = "<<timerefacc<<endl;
// where the elastic pulse is expected to be
Int_t binmin=30;
Int_t binmax=109;

// Input root file
TChain *tree=new TChain("T");
// tree->Add(Form("/cache/hallc/c-nps/analysis/pass1/replays/production/elastics/nps_hms_coin_%d_%d_1_-1.root", runnumber, seg)); //elastique

 tree->Add(Form("../nps_hms_coin_%d_%d_1_-1.root", runnumber, seg)); //production data

 //if(runnumber==5427)tree->Add(Form("rootfiles/run_%d_0_1000_ps4.root",runnumber)); //production 5138 ps6 et 5427 ps4
 //if(runnumber==5138)tree->Add(Form("rootfiles/run_%d_0_1000_ps6.root",runnumber)); //production 5138 ps6 et 5427 ps4


 Int_t nentries=tree->GetEntries();
//nentries=10000;

//Matched Filter variables and histos
Float_t output[ntime];
for(Int_t i=0;i<nblocks;i++){
hsig_mf[i] = new TH1F(Form("hsig_mf%d",i),Form("hsig_mf%d : matched filter of the wf",i),ntime,0,ntime); 
hsig_mf[i]->SetLineColor(1);
hsig_mf[i]->SetLineWidth(2);
hsig_mf[i]->GetXaxis()->SetTitle("Time (4 ns)");
hsig_mf[i]->GetYaxis()->SetTitle("(arbitrary units)");
hsig_mf[i]->GetYaxis()->SetLabelSize(0.05);
hsig_trig[i] = new TH1F(Form("hsig_trig%d",i),Form("hsig_trig%d : Sum of 3x3 wf around this block (fit performed : NO)",i),ntime,0,ntime); 
hsig_trig[i]->SetLineColor(1);
hsig_trig[i]->SetLineWidth(2);
hsig_trig[i]->GetXaxis()->SetTitle("Time (4 ns)");
hsig_trig[i]->GetYaxis()->SetTitle("(mV)");
hsig_trig[i]->GetYaxis()->SetLabelSize(0.05);
}
spec = new TSpectrum();


//Read reference waveforms
const int nbparameters=2+2*maxwfpulses;
cout<<"nbparameters = "<<nbparameters<<endl;
Double_t x[ntime],y[ntime];
ifstream filewf;
Double_t dum1,ymax;
Int_t preswf[nblocks];
for(Int_t i=0;i<nblocks;i++){
  //filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/ref_wf_%d.txt",i));
 
  // filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));

  

  if (runnumber>6183 && runnumber<7500) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6171-6183/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>6168 && runnumber<6171) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>5236 && runnumber<6151) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5217-5236/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>5208 && runnumber<5217) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/5183-5208/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>3898 && runnumber<5183) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/3883-3898/fit_e_runs/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2920 && runnumber<3883) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2900-2920/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2885 && runnumber<2900) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2875-2885/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>2871 && runnumber<2875) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/2855-2871/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>1982 && runnumber<2855) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1969-1982/RWF/ref_wf_%d.txt",i));
    } 
    else if (runnumber>1560 && runnumber<1961) {
        filewf.open(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/1423-1511/RWF/ref_wf_%d.txt",i));
    }

 //filewf.open(Form("refwf/final_refwf_block_%d.txt",i));
timeref[i]=-1.e6;

 preswf[i]=0;
 
if(filewf.is_open()){

  filewf>>timeref[i]>>dum1;

  for (Int_t it=0;it<ntime;it++){filewf>>x[it]>>y[it];}

  ymax=0;

 for (Int_t it=0;it<ntime;it++){if(y[it]>ymax){ymax=y[it];timeref[i]=x[it];}}

 mfint[i]=0;

 for (Int_t it=0;it<ntime;it++){if(TMath::Abs(timeref[i]-x[it])<0.001){for(Int_t jt=0;jt<mfwidth;jt++){mfyref[i][jt]=y[it+jt-mfleft];mfint[i]+=mfyref[i][jt];}}}

 interpolation[i] = new ROOT::Math::Interpolator(ntime, ROOT::Math::Interpolation::kCSPLINE);

 interpolation[i]->SetData(ntime, x, y);

 finter[i]= new TF1(Form("finter_%d",i),func,-200,ntime+200,nbparameters);

 finter[i]->FixParameter(0,i); //numero de bloc
//for(Int_t p=0;p<maxwfpulses;p++){cout<<"p ="<<p<<endl;finter[i]->FixParameter(2+2*p,0.);finter[i]->FixParameter(3+2*p,0.);}
finter[i]->SetLineColor(4);
finter[i]->SetNpx(1100);
preswf[i]=1;
}
filewf.close();
//cout<<i<<" "<<timeref[i]<<endl;
}

//Check if all ref wf are present
Int_t probrefwf=0;
for(Int_t i=0;i<nblocks;i++){if(preswf[i]==0){probrefwf++;cout<<"REF WF for block "<<i<<" is missing !"<<endl;}}


//chunck added by wassim

Float_t tdcoffset[nblocks];


// Load base offsets
 
ifstream filetdc("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/Work_Analysis/WF/BK_TEST/TEST_BOOM/6151-6168/fit_e_runs/RWF/tdc_offset_param.txt");
 
if (filetdc.is_open()) {
    for (Int_t i = 0; i < nblocks; i++)
        filetdc >> tdcoffset[i];
    filetdc.close();
} else {
    std::cerr << "Warning: Could not open base tdc_offset_param.txt" << std::endl;
}
 

ifstream filetdc2;

// choose the right offset file depending on runnumber
if (runnumber >= 1572 && runnumber < 1622) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1572.txt");
}
else if (runnumber >= 1622 && runnumber < 1690) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1622.txt");
}
else if (runnumber >= 1690 && runnumber < 1723) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1690.txt");
}
else if (runnumber >= 1723 && runnumber < 1801) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1723.txt");
}
else if (runnumber >= 1801 && runnumber < 1853) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1801.txt");
}
else if (runnumber >= 1853 && runnumber < 1976) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1853.txt");
}
else if (runnumber >= 1976 && runnumber < 2079) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_1976.txt");
}
else if (runnumber >= 2079 && runnumber < 2162) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2079.txt");
}
else if (runnumber >= 2162 && runnumber < 2257) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2162.txt");
}
else if (runnumber >= 2257 && runnumber < 2352) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2257.txt");
}
else if (runnumber >= 2352 && runnumber < 2466) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2352.txt");
}
else if (runnumber >= 2466 && runnumber < 2569) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2466.txt");
}
else if (runnumber >= 2569 && runnumber < 2663) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2569.txt");
}
else if (runnumber >= 2663 && runnumber < 2769) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2663.txt");
}
else if (runnumber >= 2769 && runnumber < 2914) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2769.txt");
}
else if (runnumber >= 2914 && runnumber < 3058) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_2914.txt");
}
else if (runnumber >= 3058 && runnumber < 3164) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3058.txt");
}
else if (runnumber >= 3164 && runnumber < 3269) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3164.txt");
}
else if (runnumber >= 3269 && runnumber < 3383) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3269.txt");
}
else if (runnumber >= 3383 && runnumber < 3486) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3383.txt");
}
else if (runnumber >= 3486 && runnumber < 3593) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3486.txt");
}
else if (runnumber >= 3593 && runnumber < 3728) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3593.txt");
}
else if (runnumber >= 3728 && runnumber < 3776) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3728.txt");
}
else if (runnumber >= 3776 && runnumber < 3855) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3776.txt");
}
else if (runnumber >= 3855 && runnumber < 3955) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3855.txt");
}
else if (runnumber >= 3955 && runnumber < 4058) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_3955.txt");
}
else if (runnumber >= 4058 && runnumber < 4156) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4058.txt");
}
else if (runnumber >= 4156 && runnumber < 4260) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4156.txt");
}
else if (runnumber >= 4260 && runnumber < 4362) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4260.txt");
}
else if (runnumber >= 4362 && runnumber < 4456) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4362.txt");
}
else if (runnumber >= 4456 && runnumber < 4557) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4456.txt");
}
else if (runnumber >= 4557 && runnumber < 4636) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4557.txt");
}
else if (runnumber >= 4636 && runnumber < 4718) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4636.txt");
}
else if (runnumber == 4718) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4718_test.txt");
}
else if (runnumber >= 4760 && runnumber < 4851) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4760.txt");
}
else if (runnumber >= 4851 && runnumber < 4967) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4851.txt");
}
else if (runnumber >= 4967 && runnumber < 5054) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_4967.txt");
}
else if (runnumber >= 5054 && runnumber < 5172) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5054.txt");
}
else if (runnumber >= 5172 && runnumber < 5272) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5172.txt");
}
else if (runnumber >= 5272 && runnumber < 5356) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5272.txt");
}
else if (runnumber >= 5356 && runnumber < 5524) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5356.txt");
}
else if (runnumber >= 5524 && runnumber < 5577) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5524.txt");
}
else if (runnumber >= 5577 && runnumber < 5657) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5577.txt");
}
else if (runnumber >= 5657 && runnumber < 5751) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5657.txt");
}
else if (runnumber >= 5751 && runnumber < 5864) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5751.txt");
}
else if (runnumber >= 5864 && runnumber < 5967) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5864.txt");
}
else if (runnumber >= 5967 && runnumber < 6055) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_5967.txt");
}
else if (runnumber >= 6055 && runnumber < 6153) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6055.txt");
}
else if (runnumber >= 6153 && runnumber < 6253) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6153.txt");
}
else if (runnumber >= 6253 && runnumber < 6352) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6253.txt");
}
else if (runnumber >= 6352 && runnumber < 6457) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6352.txt");
}
else if (runnumber >= 6457 && runnumber < 6557) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6457.txt");
}
else if (runnumber >= 6557 && runnumber < 6661) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6557.txt");
}
else if (runnumber >= 6661 && runnumber < 6758) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6661.txt");
}
else if (runnumber >= 6758 && runnumber < 6853) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6758.txt");
}
else if (runnumber >= 6853 && runnumber < 6942) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6853.txt");
}
else if (runnumber >= 6942) {
    filetdc2.open("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/TEST/offsets/Offsets_6942.txt");
}
else {
    std::cerr << "No offset file match for run " << runnumber << std::endl;
}


if (filetdc2.is_open()) {
    double dum = 0;
    for (Int_t i = 0; i < nblocks; i++) {
        filetdc2 >> dum;
        tdcoffset[i] += dum;
    }
    filetdc2.close();
} else {
    std::cerr << "Warning: Could not open secondary offset file Offsets_2914.txt" << std::endl;
}


 
 
//Lecture des corrections timing (obtenu au step 0 de analyse.C du setting elastique 1974)
ifstream filetime;
Float_t rmstime[nblocks];
filetime.open("filetime_step_i.txt");
for(Int_t i=0;i<nblocks;i++){filetime>>dum>>cortime[i]>>dum>>rmstime[i]>>dum;}


// Input variables
// T->Print(); to see the type (int, double, ...) of each variable in the rootfile
tree->SetBranchStatus("*", false);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampWaveform",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampWaveform",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcCounter",true);
tree->SetBranchStatus("NPS.cal.fly.adcCounter",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseAmp",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseAmp",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseInt",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseInt",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPed",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPed",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTime",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTime",true);
tree->SetBranchStatus("Ndata.NPS.cal.fly.adcSampPulseTimeRaw",true);
tree->SetBranchStatus("NPS.cal.fly.adcSampPulseTimeRaw",true);
tree->SetBranchStatus("T.hms.hT1_tdcTime",true);
tree->SetBranchStatus("T.hms.hT2_tdcTime",true);
tree->SetBranchStatus("T.hms.hT3_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG1_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG2_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG3_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG4_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG5_tdcTime",true);
tree->SetBranchStatus("T.hms.hTRIG6_tdcTime",true);
tree->SetBranchStatus("H.react.x",true);
tree->SetBranchStatus("H.react.y",true);
tree->SetBranchStatus("H.react.z",true);
tree->SetBranchStatus("H.gtr.dp",true);
tree->SetBranchStatus("H.gtr.th",true);
tree->SetBranchStatus("H.gtr.ph",true);
tree->SetBranchStatus("H.gtr.px",true);
tree->SetBranchStatus("H.gtr.py",true);
tree->SetBranchStatus("H.gtr.pz",true);
tree->SetBranchStatus("H.gtr.beta",true);
tree->SetBranchStatus("H.hod.beta",true);
tree->SetBranchStatus("H.cer.npeSum",true);
tree->SetBranchStatus("H.cal.etottracknorm",true);
tree->SetBranchStatus("g.evnum",true);


/////
Int_t NSampWaveForm; tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampWaveform",&NSampWaveForm);//should correspond to Ndata
Double_t SampWaveForm[Ndata]; tree->SetBranchAddress("NPS.cal.fly.adcSampWaveform",&SampWaveForm); // Contain the waveforms for the tested calo blocks

Int_t  NadcCounter;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcCounter",&NadcCounter) ;// total number of pulses above theshold (could have up to 4 pulses per block)
Double_t adcCounter[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcCounter",&adcCounter) ;

Int_t NadcSampPulseAmp;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseAmp",&NadcSampPulseAmp) ;
Double_t adcSampPulseAmp[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseAmp",&adcSampPulseAmp) ;

Int_t NadcSampPulseInt;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseInt",&NadcSampPulseInt) ;
Double_t adcSampPulseInt[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseInt",&adcSampPulseInt) ;

Int_t NadcSampPulsePed;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPed",&NadcSampPulsePed) ;
Double_t adcSampPulsePed[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPed",&adcSampPulsePed) ;

Int_t NadcSampPulseTime;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseTime",&NadcSampPulseTime) ;
Double_t adcSampPulseTime[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseTime",&adcSampPulseTime) ;

Int_t NadcSampPulseTimeRaw;tree->SetBranchAddress("Ndata.NPS.cal.fly.adcSampPulseTimeRaw",&NadcSampPulseTimeRaw) ;
Double_t adcSampPulseTimeRaw[nsignals];tree->SetBranchAddress("NPS.cal.fly.adcSampPulseTimeRaw",&adcSampPulseTimeRaw) ;

Double_t hT1_tdcTime;tree->SetBranchAddress("T.hms.hT1_tdcTime",&hT1_tdcTime) ;
Double_t hT2_tdcTime;tree->SetBranchAddress("T.hms.hT2_tdcTime",&hT2_tdcTime) ;
Double_t hT3_tdcTime;tree->SetBranchAddress("T.hms.hT3_tdcTime",&hT3_tdcTime) ;
Double_t hTRIG1_tdcTime;tree->SetBranchAddress("T.hms.hTRIG1_tdcTime",&hTRIG1_tdcTime) ;
Double_t hTRIG2_tdcTime;tree->SetBranchAddress("T.hms.hTRIG2_tdcTime",&hTRIG2_tdcTime) ;
Double_t hTRIG3_tdcTime;tree->SetBranchAddress("T.hms.hTRIG3_tdcTime",&hTRIG3_tdcTime) ;
Double_t hTRIG4_tdcTime;tree->SetBranchAddress("T.hms.hTRIG4_tdcTime",&hTRIG4_tdcTime) ;
Double_t hTRIG5_tdcTime;tree->SetBranchAddress("T.hms.hTRIG5_tdcTime",&hTRIG5_tdcTime) ;
Double_t hTRIG6_tdcTime;tree->SetBranchAddress("T.hms.hTRIG6_tdcTime",&hTRIG6_tdcTime) ;
Double_t vx;tree->SetBranchAddress("H.react.x",&vx) ;
Double_t vy;tree->SetBranchAddress("H.react.y",&vy) ;
Double_t vz;tree->SetBranchAddress("H.react.z",&vz) ;
Double_t dp;tree->SetBranchAddress("H.gtr.dp",&dp) ;
Double_t th;tree->SetBranchAddress("H.gtr.th",&th) ;
Double_t ph;tree->SetBranchAddress("H.gtr.ph",&ph) ;
Double_t px;tree->SetBranchAddress("H.gtr.px",&px) ;
Double_t py;tree->SetBranchAddress("H.gtr.py",&py) ;
Double_t pz;tree->SetBranchAddress("H.gtr.pz",&pz) ;
Double_t gtrbeta;tree->SetBranchAddress("H.gtr.beta",&gtrbeta) ;
Double_t hodbeta;tree->SetBranchAddress("H.hod.beta",&hodbeta) ;
Double_t cernpe;tree->SetBranchAddress("H.cer.npeSum",&cernpe) ;
Double_t caltracknorm;tree->SetBranchAddress("H.cal.etottracknorm",&caltracknorm) ;
Double_t hcnps_intlk_cz_t_back[ntemp];
Double_t hcnps_intlk_cz_t_front[ntemp];
Double_t g_event;tree->SetBranchAddress("g.evnum", &g_event); 
// Enable only the required branches
for (int i = 0; i < ntemp; ++i) {
    tree->SetBranchStatus(Form("hcnps_intlk_cz_t_back_%i", i+1), 1);
    tree->SetBranchStatus(Form("hcnps_intlk_cz_t_front_%i", i+1), 1);
    tree->SetBranchAddress(Form("hcnps_intlk_cz_t_back_%i", i+1), &hcnps_intlk_cz_t_back[i]);
    tree->SetBranchAddress(Form("hcnps_intlk_cz_t_front_%i", i+1), &hcnps_intlk_cz_t_front[i]);
}
 
 cout<<"Nb of entries in rootfiles = "<<nentries<<endl;

Int_t tracage = 0; // 1 to see the waveforms event by event
TFile *fout = nullptr;
TTree *treeout = nullptr;

if (tracage == 0) {
  //TString rootfilePath = Form("/volatile/hallc/nps/wassim/ELASTICS_HMS_CALIBRATED/TEST_MFILTER/nps_production_%d_%d_%d.root", runnumber, seg, job_number);
  // TString rootfilePath = Form("/volatile/hallc/nps/wassim/nps_production_pass0_%d_%d_%d.root", runnumber, seg, job_number);
  //TString rootfilePath = Form("/cache/hallc/c-nps/users/wassim/KINC_x60_3/nps_production_%d_%d_%d.root", runnumber, seg, job_number);
  
 TString rootfilePath = Form("nps_production_MF_%d_%d_%d.root", runnumber, seg, job_number);

  
    std::cout << "Rootfile location: " << rootfilePath << std::endl;

    fout = new TFile(rootfilePath, "RECREATE");  //  Open file first
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error opening ROOT file: " << rootfilePath << std::endl;
        return;
    }

    treeout = new TTree("T", "Tree organized");  //  Create tree after file
    treeout->SetAutoFlush(-300000000);
 }



/*
// Output rootfile
Int_t tracage=1; //1 to see the waveforms event by event
TFile *fout;
if(tracage==0) fout=new TFile(Form("rootfiles/nps_elaswf4t_%d_%d_%d.root",runnumber,(int)(trig_thres),coinc_width),"recreate");
if(tracage==1) fout=new TFile(Form("rootfiles/nps_elaswf4t_%d_tracage.root",runnumber),"recreate");
TTree *treeout=new TTree("T","Tree organized");
treeout->SetAutoFlush(-300000000);
*/
// Output variables
Double_t signal[nblocks][ntime];//treeout->Branch("signal",&signal,Form("signal[%d][%d]/D",nblocks,ntime));//if we want to save the raw waveforms in the output rootfile
Double_t ampl2[nblocks];//treeout->Branch("ampl2",&ampl2,Form("ampl2[%d]/D",nblocks));//amplitude of the pulse relatively to background
Double_t ampl[nblocks];treeout->Branch("ampl",&ampl,Form("ampl[%d]/D",nblocks));//amplitude of the pulse (a comparer avec Sampampl)
Double_t amplwf[nblocks];treeout->Branch("amplwf",&amplwf,Form("amplwf[%d]/D",nblocks));//amplwfitude of the pulse (a comparer avec Sampamplwf)
Double_t baselinewf[nblocks];treeout->Branch("baselinewf",&baselinewf,Form("baselinewf[%d]/D",nblocks));//baselinewfitude of the pulse (a comparer avec Sampbaselinewf)
treeout->Branch("wfampl",&wfampl,Form("wfampl[%d][%d]/D",nblocks,maxwfpulses));
treeout->Branch("wftime",&wftime,Form("wftime[%d][%d]/D",nblocks,maxwfpulses));
//Double_t wftimei[nblocks][maxwfpulses]; treeout->Branch("wftimei",&wftimei,Form("wftimei[%d][%d]/D",nblocks,maxwfpulses));
Double_t Npulse[nblocks];treeout->Branch("Npulse",&Npulse,Form("Npulse[%d]/D",nblocks));//number of pulses in one block 
treeout->Branch("wfnpulse",&wfnpulse,Form("wfnpulse[%d]/I",nblocks));//number of pulses in one block 
//treeout->Branch("wfnpulse_sev",&wfnpulse_sev,Form("wfnpulse_sev[%d][%d]/I",nblocks,maxwfpulses)); // 1 if the pulse can be detected with a more severe condition
//treeout->Branch("wfnpulsesup",&wfnpulsesup,Form("wfnpulsesup[%d]/I",nblocks));//1 if a supl fit was performed in the pulse region, 0 otherwise 
//treeout->Branch("wfnpulsei",&wfnpulsei,Form("wfnpulsei[%d]/I",nblocks));//number of initial pulses in one block found by the algo 
Int_t wfnpulse10[nblocks];//treeout->Branch("wfnpulse10",&wfnpulse10,Form("wfnpulse10[%d]/I",nblocks));//number of pulses (>10mV) in one block 
Int_t wfnpulse1[nblocks];//treeout->Branch("wfnpulse1",&wfnpulse1,Form("wfnpulse1[%d]/I",nblocks));//number of pulses (>1mV) in one block 
Double_t Sampampl[nblocks];treeout->Branch("Sampampl",&Sampampl,Form("Sampampl[%d]/D",nblocks));//amplitude of the 1st pulse 
Double_t Samptime[nblocks];treeout->Branch("Samptime",&Samptime,Form("Samptime[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t Sampped[nblocks];treeout->Branch("Sampped",&Sampped,Form("Sampped[%d]/D",nblocks));//ped of the 1st pulse 'ns)
Double_t Sampener[nblocks];//treeout->Branch("Sampener",&Sampener,Form("Sampener[%d]/D",nblocks));//integral of the 1st pulse 
Double_t ener[nblocks];treeout->Branch("ener",&ener,Form("ener[%d]/D",nblocks));//integral of waveforms in the pulse window - integral of waveform in the backgroung window
Double_t time[nblocks];//treeout->Branch("time",&time,Form("time[%d]/D",nblocks));//time position of the pulse maximum
Double_t timewf[nblocks];treeout->Branch("timewf",&timewf,Form("timewf[%d]/D",nblocks));//timewf position of the pulse 
Double_t chi2[nblocks];treeout->Branch("chi2",&chi2,Form("chi2[%d]/D",nblocks));//time of the 1st pulse 'ns)
Double_t bkg[nblocks];treeout->Branch("bkg",&bkg,Form("bkg[%d]/D",nblocks));//integral of waveforms in the backgroung window
Double_t enertot;treeout->Branch("enertot",&enertot,"enertot/D");// sum of all ener[i] : total energy deposited in the calo (!! energy calibration is not done yet)
Double_t integ[nblocks];//treeout->Branch("integ",&integ,Form("integ[%d]/D",nblocks));//integral of waveforms in [0,dt*ntime ns] window
Double_t integtot;//treeout->Branch("integtot",&integtot,"integtot/D");// sum of all integ[i]
Double_t noise[nblocks];treeout->Branch("noise",&noise,Form("noise[%d]/D",nblocks));//RMS of waveforms relatively to zero
Double_t larg50[nblocks];//treeout->Branch("larg50",&larg50,Form("larg50[%d]/D",nblocks));//RMS of waveforms relatively to zero
Double_t larg90[nblocks];//treeout->Branch("larg90",&larg90,Form("larg90[%d]/D",nblocks));//RMS of waveforms relatively to zero
treeout->Branch("trig",&trig,Form("trig[%d]/I",nblocks));//1 if the waveform data is available
treeout->Branch("pres",&pres,Form("pres[%d]/I",nblocks));//1 if the waveform data is available
treeout->Branch("preswf",&preswf,Form("preswf[%d]/I",nblocks));//1 if the ref shapes is available
Int_t nrun;treeout->Branch("nrun",&nrun,"nrun/I");// run number
treeout->Branch("job_number",&job_number,"job_number/I");// run number  
Int_t event;treeout->Branch("event",&event,"event/I");// event number in input rootfiles
treeout->Branch("dp",&dp,"dp/D");
treeout->Branch("th",&th,"th/D");
treeout->Branch("ph",&ph,"ph/D");
treeout->Branch("px",&px,"px/D");
treeout->Branch("py",&py,"py/D");
treeout->Branch("pz",&pz,"pz/D");
treeout->Branch("vx",&vx,"vx/D");
treeout->Branch("vy",&vy,"vy/D");
treeout->Branch("vz",&vz,"vz/D");
treeout->Branch("cernpe",&cernpe,"cernpe/D");
treeout->Branch("caltracknorm",&caltracknorm,"caltracknorm/D");
treeout->Branch("gtrbeta",&gtrbeta,"gtrbeta/D");
treeout->Branch("hodbeta",&hodbeta,"hodbeta/D");
treeout->Branch("corr_time_HMS",&corr_time_HMS,"corr_time_HMS/D");
//treeout->Branch("timerefacc",&timerefacc,"timerefacc/D");
Double_t cputime;treeout->Branch("cputime",&cputime,"cputime/D");
treeout->Branch("g_event", &g_event, "g_event/D");

//Other variables
Double_t sigmax[nblocks];
Int_t bloc,nsamp,ns,ilin,icol,ilinc,icolc,inp;
Double_t max50,max90,min50,min90;
Int_t nbfit=0;
Int_t nsamp2[nblocks];

//Problem variables
Int_t nsampwf=0;
Int_t ndataprob=0;
Int_t probnsamp=0;
Int_t probslot=0;
Int_t probHMStime=0;

//Histograms
for(Int_t i=0;i<nblocks;i++){
hsig_i[i] = new TH1F(Form("hsig_i%d",i),Form("hsig_i%d",i),ntime,0,ntime); 
hsig_i[i]->SetLineColor(1);
hsig_i[i]->SetLineWidth(2);
hsig_i[i]->GetXaxis()->SetTitle("Time (4 ns)");
hsig_i[i]->GetYaxis()->SetTitle("(mV)");
hsig_i[i]->GetYaxis()->SetLabelSize(0.05);
}

// Variables du tracage
TCanvas *c1 = new TCanvas("c1","c1",1300,1000);c1->Divide(ncol,nlin,0,0);
TCanvas *c11 = new TCanvas("c11","c11",700,500);
TCanvas *c12 = new TCanvas("c12","c12",700,500);
TCanvas *c13 = new TCanvas("c13","c13",700,500);

TLatex* tex = new TLatex();
tex->SetTextSize(0.015);
TH1::AddDirectory(kFALSE);
TLine *ln = new TLine();ln->SetLineStyle(2);ln->SetLineColor(2);
TLine *ln2 = new TLine();ln2->SetLineStyle(2);ln2->SetLineColor(3);ln2->SetLineWidth(2);

Float_t timemean2[nblocks];
for(Int_t ii=0;ii<nblocks;ii++){timemean2[ii]=170+timerefacc*dt;} //where the coincidence pulse is expected (ns units)

TH1F *h1time = new TH1F("h1time","h1time",200,-50,50);
TH1F *h2time = new TH1F("h2time","h2time",200,-100,100);


/////// ANALYZE //////////////////////////////////////////////////////////
 
start_event = job_number * evts_per_job;
 
end_event = (job_number+1)* evts_per_job;
if(end_event > nentries) end_event = nentries;

for(Int_t evt=start_event; evt<end_event; evt++){

    tree->GetEntry(evt);
    event=evt;
    nrun=runnumber;
    
    if (evt%1000==0) {cout << " Entry = " << evt <<"/"<<nentries<<"  run="<<nrun<<"  cpu time="<<t.RealTime()<<endl;t.Continue();}
    
    if(TMath::Abs(th)<0.08&&TMath::Abs(ph)<0.04&&TMath::Abs(dp)<10&&TMath::Abs(vz)<4){
    // if(TMath::Abs(th)<0.08&&TMath::Abs(ph)<0.04&&TMath::Abs(dp)<8&&dp>250*ph-2&&dp<250*ph+3&&TMath::Abs(vz)<4){
      
    if(NSampWaveForm>Ndata){nsampwf++; cout<<"!!!! NSampWaveForm problem  "<<evt<<"  "<<NSampWaveForm<<" "<<Ndata<<endl;}

    if(NSampWaveForm<=Ndata){ //NSampWaveForm must be <= Ndata (otherwise correct Ndata value)
    for(Int_t i=0;i<nblocks;i++){pres[i]=0;minsignal[i]=1000;for(Int_t j=0;j<ntime;j++){signal[i][j]=0;hsig_mf[i]->SetBinContent(j+1,0.);}}//liste de presence des blocs! pres=0 si bloc absent, pres=1 s'il est present
 
//Extract the data from the complex variable SampWaveForm[] (NPS.cal.fly.adcSampWaveform)
    ns=0; //ns represent for a given event the element number of the NPS.cal.fly.adcSampWaveform variable 
    while (ns < NSampWaveForm) {
    bloc=SampWaveForm[ns];ns++; //bloc number (actually the slot number)
    nsamp=SampWaveForm[ns];ns++; //time samples (should be equal to ntime (100))
if(bloc==2000)bloc=1080; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
if(bloc==2001)bloc=1081; //modification of the bloc number because the fADC (16 slots) corresponding to slot 720-736 is used for the scintillators
if(bloc<0||bloc>nslots-0.5){cout<<"slot number problem "<<evt<<" "<<bloc<<endl;ns=NSampWaveForm+1;}//to exit the while() 
if(nsamp!=ntime){cout<<"nsamp warning "<<evt<<" "<<nsamp<<" "<<bloc<<" (Fit will be performed only for "<<nsamp<<" samples instead of 110)"<<endl;probnsamp++;} //just a check!
if(bloc>-0.5&&bloc<nslots){// that's what we expect!
//cout<<"lecture"<<endl;
nsamp2[bloc]=nsamp;
pres[bloc]=1;
for (Int_t it=0;it<nsamp;it++){
if(bloc>-0.5&&bloc<nblocks){signal[bloc][it]=SampWaveForm[ns];hsig_i[bloc]->SetBinContent(it+1,signal[bloc][it]);minsignal[bloc]=TMath::Min(minsignal[bloc],signal[bloc][it]);}
ns++;
}
}//fin if(bloc number is good)
    }//fin while()


//initialisation
enertot=0;integtot=0;
for(Int_t i=0;i<nblocks;i++){timewf[i]=-100;amplwf[i]=-100;;baselinewf[i]=-300;chi2[i]=-1;ener[i]=0;integ[i]=0;noise[i]=0;bkg[i]=0;sigmax[i]=-100;Sampampl[i]=-100;Sampped[i]=-100;Samptime[i]=-100;Sampener[i]=-100;Npulse[i]=0;trig[i]=-100;}
for(Int_t i=0;i<nblocks;i++){wfnpulse[i]=-100;wfnpulsesup[i]=-100;wfnpulsei[i]=-100;wfnpulse10[i]=-100;wfnpulse1[i]=-100;for(Int_t p=0;p<maxwfpulses;p++){wftime[i][p]=-300;wfampl[i][p]=-100;wfnpulse_sev[i][p]=-100;}}
for(Int_t i=0;i<nblocks;i++){hsig_i[i]->SetLineColor(1);}

//Read the hcana calculated variables
if(!(NadcCounter==NadcSampPulseAmp&&NadcSampPulseInt==NadcSampPulsePed&&NadcSampPulseInt==NadcSampPulseTime)){cout<<"!!!!! Problem Ndata !!!!!! "<<evt<<endl;ndataprob++;}
corr_time_HMS=0;
for(Int_t iNdata = 0; iNdata < NadcCounter; iNdata++){
if(adcCounter[iNdata] ==2000) adcCounter[iNdata]=1080;// nouveau numero du scintillateur attribue par Malek
if(adcCounter[iNdata] ==2001) adcCounter[iNdata]=1081;// nouveau numero du scintillateur attribue par Malek


//determine HMS time correction
if(iNdata==0)corr_time_HMS=adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])];
if(iNdata!=0){if(TMath::Abs(corr_time_HMS-(adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])]))>0.001){cout<<"problem HMS time correction event "<<evt<<endl;probHMStime++;}}

if(!(adcCounter[iNdata] >= 0 && adcCounter[iNdata] <nblocks+2)/*&&adcCounter[iNdata]!=196*/)cout<<"****** Problem adcCounter ******* "<<evt<<" "<<iNdata<<" "<<adcCounter[iNdata]<<endl;
if(adcCounter[iNdata] >= 0 && adcCounter[iNdata] <nblocks){
Npulse[(int)(adcCounter[iNdata])]+=1;
//hsig_i[(int)(adcCounter[iNdata])]->SetLineColor(2);
if(Npulse[(int)(adcCounter[iNdata])]==1){
Sampampl[(int)(adcCounter[iNdata])]=adcSampPulseAmp[iNdata];
Samptime[(int)(adcCounter[iNdata])]=adcSampPulseTime[iNdata];
Sampener[(int)(adcCounter[iNdata])]=adcSampPulseInt[iNdata];
Sampped[(int)(adcCounter[iNdata])]=adcSampPulsePed[iNdata];
}
if(Npulse[(int)(adcCounter[iNdata])]>1){
if(TMath::Abs(Samptime[(int)(adcCounter[iNdata])]-timemean2[(int)(adcCounter[iNdata])])>TMath::Abs(adcSampPulseTime[iNdata]-timemean2[(int)(adcCounter[iNdata])])){//si le 2eme pulse est plus proche du temps de reference que le 1er pulse alors on prend le 2eme
Sampampl[(int)(adcCounter[iNdata])]=adcSampPulseAmp[iNdata];
Samptime[(int)(adcCounter[iNdata])]=adcSampPulseTime[iNdata];
Sampener[(int)(adcCounter[iNdata])]=adcSampPulseInt[iNdata];
Sampped[(int)(adcCounter[iNdata])]=adcSampPulsePed[iNdata];
}}
}}

//Apply matched filter + threshold on 3x3 blocks
for(Int_t i=0;i<nblocks;i++){
if(pres[i]==1&&preswf[i]==1){ 
for(Int_t it=0;it<nsamp;it++){
hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(hsig_i[i]->GetBinContent(it+1)*4.096/2.))/4.096);
if(hsig_i[i]->GetBinContent(it+1)<1.)hsig_i[i]->SetBinError(it+1,TMath::Sqrt(TMath::Abs(1.*4.096/2.))/4.096);
//hsig_i[i]->SetBinError(it+1,0.5); //Charles suggestion
}
mfmin=1.e6;
for(Int_t it=mfleft;it<nsamp-mfright;it++){
output[it]=0.;
for(Int_t jt=0;jt<mfwidth;jt++){output[it]+=(hsig_i[i]->GetBinContent(it+jt-mfright+1)-minsignal[i])*mfyref[i][mfwidth-1-jt]/mfint[i];}
mfmin=TMath::Min(mfmin,output[it]);
}
for(Int_t it=mfleft;it<nsamp-mfright;it++){hsig_mf[i]->SetBinContent(it+1,output[it]-mfmin);}
npulsewf2(i,0,TMath::Min(ntime,nsamp2[i]));
Apply_thres(i);
}}

// Fit of the wf
for(Int_t i=0;i<nblocks;i++){if(pres[i]==1&&preswf[i]==1&&trig[i]==1&&nsamp==ntime){Fitwf2(i,0,TMath::Min(ntime,nsamp2[i]));nbfit++;hsig_i[i]->SetLineColor(2);}}
// Get the fit results
for(Int_t i=0;i<nblocks;i++){
for(Int_t p=0;p<maxwfpulses;p++){wftime[i][p]=-300;wfampl[i][p]=-100;/*wftimei[i][p]=-300;*/}//important!!
if(pres[i]==1&&preswf[i]==1&&trig[i]==1&&nsamp==ntime){ 
wfnpulsei[i]=wfnpulse[i];
wfnpulse1[i]=0;
cpulse=-100;
for(Int_t p=0;p<TMath::Max(wfnpulse[i],1);p++){if(wfampl[i][p]>1&&p!=cpulse)wfnpulse1[i]++;}
if(tracage==1)cout<<i<<"  cpulse ="<<cpulse<<" wfnpulse1= "<<wfnpulse1[i]<<endl;
if(wfnpulse[i]>0){
chi2[i]=finter[i]->GetChisquare()/finter[i]->GetNDF();
wfnpulse10[i]=0;
baselinewf[i]=finter[i]->GetParameter(1);
for(Int_t p=0;p<TMath::Max(wfnpulse[i],1);p++){
//wftimei[i][p]=finter[i]->GetParameter(2+2*p);
wftime[i][p]=finter[i]->GetParameter(2+2*p)*dt+corr_time_HMS-cortime[i]-timerefacc*dt;//temps du pulse en ns
wfampl[i][p]=finter[i]->GetParameter(3+2*p); //amplitude du pulse en mV
if(wfampl[i][p]>10)wfnpulse10[i]++;
if(wfampl[i][p]>20){h2time->Fill(wftime[i][p],1.);h1time->Fill(finter[i]->GetParameter(2+2*p)-timerefacc+corr_time_HMS/dt,1.);}
if(p==0){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}
if(p>0){if(TMath::Abs(wftime[i][p])<TMath::Abs(timewf[i])){timewf[i]=wftime[i][p];amplwf[i]=wfampl[i][p];}}//pour prendre le pulse dont le temps est le plus proche de timerefacc
}}}}

//Calculation
for(Int_t i=0;i<nblocks;i++){
if(pres[i]==1){
for(Int_t it=0;it<nsamp2[i];it++){ 
integ[i]+=signal[i][it];
integtot+=signal[i][it];
if(it>binmin&&it<binmax){//cosmic pulse window
//if(it>47&&it<84){//cosmic pulse window
ener[i]+=signal[i][it];
enertot+=signal[i][it];
}
if(!(it>binmin&&it<binmax)){bkg[i]+=signal[i][it];}//background window
if(signal[i][it]>sigmax[i]){time[i]=it;sigmax[i]=signal[i][it];ampl[i]=signal[i][it];}//determine the pulse maximum
}

ener[i]-=bkg[i]*(binmax-binmin-1)/(nsamp2[i]-(binmax-binmin-1));//subtract the bkg contribution (we have to normalize since the window widths are not the same)
bkg[i]=bkg[i]/(nsamp2[i]-(binmax-binmin-1));//mean value of bkg
for(Int_t it=0;it<nsamp2[i];it++){ 
if(!(it>binmin&&it<binmax)){noise[i]+=(signal[i][it]-bkg[i])*(signal[i][it]-bkg[i])/(nsamp2[i]-(binmax-binmin-1));}// RMS of the bkg
}
noise[i]=TMath::Sqrt(noise[i]);
}}

// Calculation of signal widths
for(Int_t i=0;i<nblocks;i++){
if(pres[i]==1){
ampl2[i]=ampl[i]-bkg[i];//amplitude of the pulse relatively to bkg
max50=0;max90=50;min50=100;min90=100;
for(Int_t it=time[i];it<nsamp2[i];it++){//aller vers la droite du maximum
if((signal[i][it]-bkg[i])>=ampl2[i]*0.5)max50=it;
if((signal[i][it]-bkg[i])>=ampl2[i]*0.1)max90=it;
}
for(Int_t it=time[i];it>-0.5;it--){//aller vers la gauche du maximum
if((signal[i][it]-bkg[i])>=ampl2[i]*0.5)min50=it;
if((signal[i][it]-bkg[i])>=ampl2[i]*0.1)min90=it;
}
larg50[i]=max50-min50;
larg90[i]=max90-min90;
}}

	      
 
//// Tracage
if(tracage==1){ //draw the waveforms event by event
    gStyle->SetOptStat(0);    
    c1->Update();
    c1->Clear();
    c1->Divide(ncol,nlin,0,0);
    for(Int_t i=0;i<nblocks;i++){
    if(pres[i]==1){//si le bloc est present

    icol=i%ncol;//numero de la colonne a laquelle le bloc appartient
    ilin=(int)(i/ncol);//numero de la ligne a laquelle le bloc appartient
    ilinc=nlin-ilin-1; //changement de numerotation pour la rendre compatible avec le numero de pad dans la canvas
//    ilinc=ilin+1; //changement de numerotation pour la rendre compatible avec le numero de pad dans la canvas
//    icolc=ncol-icol; //changement de numerotation pour la rendre compatible avec le numero de pad dans la canvas
    icolc=icol+1; //changement de numerotation pour la rendre compatible avec le numero de pad dans la canvas
    c1->cd(icolc+ilinc*ncol);
    tex->SetNDC(1);
    tex->SetTextColor(4);
    tex->SetTextSize(0.6);
    //if(evt==2586){hsig_i[i]->SetMinimum(-4.99);hsig_i[i]->SetMaximum(15.01);}
    hsig_i[i]->Draw();if(preswf[i]==1&&trig[i]==1)finter[i]->DrawF1(0,TMath::Min(ntime,nsamp2[i]),"same");
    //tex->DrawLatex(5,hsig_i[i]->GetYaxis()->GetXmin(),Form("%d",i));
    //if(i==804||i==901)
    tex->DrawLatex(0.1,0.1,Form("%d",i));
    }}
    c1->cd();
 //   tex->SetTextColor(2);tex->DrawLatex(0.1,0.989,"raw waveforms");
    c1->Update();
   
    cout <<"current event = "<<evt<<"   NadcCounter="<<NadcCounter<<endl;
    cout<<"blocs touches : "<<endl;
    for(Int_t iNdata = 0; iNdata < NadcCounter; iNdata++){
    cout<<"num="<<adcCounter[iNdata]<<" col="<<(int)(adcCounter[iNdata])%ncol<<" lin="<<(int)(adcCounter[iNdata]/ncol)<<" ";
    cout<<"  amp(mV)="<<adcSampPulseAmp[iNdata]<<"  integ(pC)="<<adcSampPulseInt[iNdata]<<"  integ(adc)="<<adcSampPulseInt[iNdata]/integtopC<<" ";
    cout<<"  integ(mV)="<<ener[(int)(adcCounter[iNdata])]<<"  integ(can)="<<ener[(int)(adcCounter[iNdata])]/ADCtomV<<"  ";
    cout<<"  ped(can)="<<adcSampPulsePed[iNdata]/ADCtomV<<endl;
    cout<<"  t="<<adcSampPulseTime[iNdata]/dt<<" tr="<<Samptime[(int)(adcCounter[iNdata])]<<endl;
    }
    cout<<"--------------------------------"<<endl;
    for(Int_t iNdata = 0; iNdata < NadcCounter; iNdata++){
    cout<<"num="<<adcCounter[iNdata]<<"  corr="<<adcSampPulseTime[iNdata]-adcSampPulseTimeRaw[iNdata]/16.-tdcoffset[(int)(adcCounter[iNdata])]<<endl;
    }
   cout<<endl;
    cout << " Enter 1 for next event, event number, 0 to quit, -1 to print, -2 to see a particular block" << endl;
    cin>> inp;
    if (inp >1) evt=inp-1;
    if (inp ==0) evt=nentries; //pour finir la boucle for(evt)
    if (inp ==-1) {c1->Print(Form("figures/signals_run%d_evt%d.png",nrun,evt));}
    if (inp==-2){
    while(inp!=-1){
    cout<<"Enter bloc number to zoom (-1 to skip)"<<endl;
    cin>>inp;    
    if(inp>-1&&pres[inp]==1&&inp%ncol>0){
    npeaks = spec->Search(hsig_mf[inp], 2,"nobackground", specthres);
    c11->cd();hsig_i[inp]->Draw();
    if(preswf[inp]==1&&trig[inp]==1) finter[inp]->DrawF1(1,TMath::Min(ntime,nsamp2[inp]),"same");
 //   npeaks = spec->Search(hsig_i[inp], 2.5*0.5, "", 0.01);
    for(Int_t ip=0;ip<npeaks;ip++){
    if(spec->GetPositionX()[ip]-2>mfstart&&spec->GetPositionX()[ip]-2<mfend&&spec->GetPositionY()[ip]>mfthres){
    ln->DrawLine(spec->GetPositionX()[ip]-2,hsig_i[inp]->GetMinimum(),spec->GetPositionX()[ip]-2,hsig_i[inp]->GetMaximum());
    }}
    
    for(Int_t it=mfstart;it<mfend;it++){if(hsig_mf[inp]->GetBinContent(it+1)>mfthres){    if(hsig_mf[inp]->GetBinContent(it+1)>hsig_mf[inp]->GetBinContent(it)&&hsig_mf[inp]->GetBinContent(it)>hsig_mf[inp]->GetBinContent(it-1)&&hsig_mf[inp]->GetBinContent(it-1)>hsig_mf[inp]->GetBinContent(it-2)/*&&hsig_mf[inp]->GetBinContent(it-2)>hsig_mf[inp]->GetBinContent(it-3)*/){ //4 increasing sample
    if(hsig_mf[inp]->GetBinContent(it+2)<hsig_mf[inp]->GetBinContent(it+1)&&hsig_mf[inp]->GetBinContent(it+3)<hsig_mf[inp]->GetBinContent(it+2)){//2 increasing sample
    peaksup=1;
    for(Int_t ip=0;ip<npeaks;ip++){if(TMath::Abs(hsig_mf[inp]->GetBinCenter(it+1)-spec->GetPositionX()[ip])<3.1) peaksup=0;}
    if(peaksup==1)ln->DrawLine(hsig_mf[inp]->GetBinCenter(it+1)-2,hsig_i[inp]->GetMinimum(),hsig_mf[inp]->GetBinCenter(it+1)-2,hsig_i[inp]->GetMaximum());
    }}}}
    c11->Update();
    
    c12->cd();hsig_mf[inp]->Draw();
    for(Int_t ip=0;ip<npeaks;ip++){
    if(spec->GetPositionX()[ip]-2>10&&spec->GetPositionX()[ip]-2<110&&spec->GetPositionY()[ip]>0.5){
    ln->DrawLine(spec->GetPositionX()[ip],hsig_mf[inp]->GetMinimum(),spec->GetPositionX()[ip],hsig_mf[inp]->GetMaximum());
    }}
    ln2->DrawLine(0,mfthres,nsamp,mfthres);
    //spec->Background(hsig_mf[inp],20,"same");
    c12->Update();
    cout<<"minsignal = "<<minsignal[inp]<<endl; 
    cout<<"Nb de pulses ="<<wfnpulse[inp]<<"   chi2/ndf = "<<chi2[inp]<<endl;
if(wfnpulse[inp]>0){
    for(Int_t p=0;p<TMath::Max(wfnpulse[inp],1);p++){
    cout<<"pulse n:"<<p<<" time initial="<<(wftime[inp][p]-corr_time_HMS+cortime[inp]+timerefacc*dt)/dt+timeref[inp]<<"(4ns) time corrected="<<wftime[inp][p]<<"(ns)  ampl="<<wfampl[inp][p]<<" (mV)"<<endl;
  //wftime[i][p]=finter[i]->GetParameter(2+2*p)*dt+corr_time_HMS-cortime[i]
//  wftime[i][p]=finter[i]->GetParameter(2+2*p)*dt+corr_time_HMS-cortime[i]-timerefacc*dt;//temps du pulse en ns
  }}
    cout<<"pulse retenu :  time corrected="<<timewf[inp]<<"   ampl="<<amplwf[inp]<<endl;
    
    c11->Update();
    
    c13->cd();
    if(trig[inp]==1)hsig_trig[inp]->SetTitle(Form("hsig_trig%d : Sum of 3x3 wfs around this block (fit performed : YES)",inp));
    hsig_trig[inp]->Draw();
    ln2->DrawLine(timeref[inp]+timerefacc-coinc_width,hsig_trig[inp]->GetMinimum(),timeref[inp]+timerefacc-coinc_width,hsig_trig[inp]->GetMaximum());
    ln2->DrawLine(timeref[inp]+timerefacc+coinc_width,hsig_trig[inp]->GetMinimum(),timeref[inp]+timerefacc+coinc_width,hsig_trig[inp]->GetMaximum());
    c13->Update();
    }
    }//end while
    cout<<"Enter 1 for next event"<<endl;
    cin>>inp;
    }    
 }//end if(tracage=1)
    cputime=t.RealTime();
    fout->cd();
    treeout->Fill();
    
    }// end if(NSampWaveForm<=Ndata)
    }// end if HMS cut
    if(evt%10000==0){cout<<"memory : "<<getMemoryUsageKB()<<endl;gObjectTable->Print();treeout->AutoSave("SaveSelf");
    }
    
 } //end for(evt)

/////////////////////////////////////////////////////////////////////////////////////////
Double_t h1timemax=h1time->GetBinCenter(h1time->FindFirstBinAbove(h1time->GetMaximum()*0.99));
Double_t h2timemax=h2time->GetBinCenter(h2time->FindFirstBinAbove(h2time->GetMaximum()*0.99));



/////// Write the output files //////////////////////////////////////////////////////////

cout<<"fin de l'analyse"<<endl;
ofstream fileprob(Form("/w/hallc-scshelf2102/nps/wassim/ANALYSIS/err_out_files/KINC_X36_5/fileprob4_%d_%d_%d.txt",runnumber,(int)(trig_thres),coinc_width));
cout<<probrefwf<<" "<<nsampwf<<" "<<ndataprob<<" "<<probnsamp<<" "<<probslot<<" "<<probHMStime<<" "<<h1timemax<<" "<<h2timemax<<" "<<highnbpulse<<endl;
cout<<cputime<<endl;
cout<<"total nb of performed fits="<<nbfit<<endl;
fileprob<<probrefwf<<" "<<nsampwf<<" "<<ndataprob<<" "<<probnsamp<<" "<<probslot<<" "<<probHMStime<<" "<<h1timemax<<" "<<h2timemax<<" "<<highnbpulse<<endl;
fileprob<<cputime<<endl;
fileprob<<"total nb of performed fits="<<nbfit<<endl;

  fout->cd();
  fout->Write();
  h1time->Write("h1time");
  h2time->Write("h2time");
  fout->Close();



}
