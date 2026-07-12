#ifndef MOLLERGenericDetector_h
#define MOLLERGenericDetector_h

////////////////////////////////////////////////////////////////////////////////
//
//MOLLERGenericDetector
//
// A generic detector class which could contain the following types
// of data:
// - ADC
//  - Can be single valued OR
//  - Can contain pulse information such as integral, amplitude and time
//  - Can be using Waveform data
//  - Can be single valued (leading edge) OR:
//  - Can contain both leading and trailing edge and Time-Over-Threshold
//
//  Organized into row, col, and layer structure (but makes no assumption about
//  number of rows, cols or layers are constant throughout).
//
////////////////////////////////////////////////////////////////////////////////

#include "THaNonTrackingDetector.h"
#include "TRotation.h"
#include "TVector3.h"
#include "THaDetMap.h"
#include "MOLLERElement.h"
#include "Helper.h"
#include "TFile.h"
#include "TTree.h"

namespace MOLLERModeADC {
  enum Mode{
    kNone,
    kADC,       //< Contains pulse information also
    kADCSimple, //< Does not contain pulse information (nor reference info)
    kWaveform,  //< Contains waveform data
  };
}

struct MOLLERGenericOutputData {
  std::vector<Int_t> ADCrow;
  std::vector<Int_t> ADCcol;
  std::vector<Int_t> ADClayer;
  std::vector<Int_t> ADCelemID;

  std::vector<Int_t> ped;
  std::vector<Int_t> a_mult;
  std::vector<Double_t> a;
  std::vector<Double_t> a_p;
  std::vector<Double_t> a_c;
  std::vector<Double_t> a_amp;
  std::vector<Double_t> a_amp_p;
  std::vector<Double_t> a_amp_c;
  std::vector<Double_t> a_amptrig_p;
  std::vector<Double_t> a_amptrig_c;
  std::vector<Double_t> a_time;

  std::vector<Int_t> samps_elemID;
  std::vector<Int_t> nsamps;
  std::vector<Int_t> sidx;
  std::vector<Double_t> samps;

  void clear() {
    ADCrow.clear();
    ADCcol.clear();
    ADCelemID.clear();
    ADClayer.clear();
    ped.clear();
    a_mult.clear();
    a.clear();
    a_p.clear();
    a_c.clear();
    a_amp.clear();
    a_amp_p.clear();
    a_amp_c.clear();
    a_amptrig_p.clear();
    a_amptrig_c.clear();
    a_time.clear();
    nsamps.clear();
    sidx.clear();
    samps.clear();
    samps_elemID.clear();
  }
};

struct MOLLERTHcSampData {
  Bool_t fHasMulti;
  Int_t fNSampPulses;
  Double_t fSampPed;
  Int_t fNPeakSamples;
  Double_t fPeakPedestalRatio;

  std::vector<Double_t> fSampPulseInt;
  std::vector<Double_t> fSampPulseAmp;
  std::vector<Double_t> fSampPulseTime;
  std::vector<Double_t> fSampPulseIntPedSub;
  std::vector<Double_t> fSampPulseIntMOLLERRaw;
  std::vector<Double_t> fSampPulseIntMOLLERVal;
    // ADD THESE
  std::vector<Double_t> fSampPulseTOT;
  std::vector<Int_t>    fSampPulsePileup;

  MOLLERTHcSampData() :
    fHasMulti(kFALSE),
    fNSampPulses(0),
    fSampPed(0.0),
    fNPeakSamples(0),
    fPeakPedestalRatio(0.0)
  {}
};

class MOLLERGenericDetector : public THaNonTrackingDetector {

public:
  explicit MOLLERGenericDetector( const char* name,
                                  const char* description = "",
                                  THaApparatus* a = nullptr );
  virtual ~MOLLERGenericDetector();

  virtual void Clear( Option_t* opt="" );

  void InitMultiPeakTree();

  void SetModeADC(MOLLERModeADC::Mode mode);
  void SetDisableRefADC(Bool_t b) { fDisableRefADC = b; }
  void SetStoreRawHits(Bool_t var) { fStoreRawHits = var; }
  void SetStoreEmptyElements(Bool_t b) { fStoreEmptyElements = b; }

  Bool_t WithADC() { return fModeADC != MOLLERModeADC::kNone; };

  virtual Int_t Decode( const THaEvData& );
  virtual Int_t CoarseProcess(TClonesArray& tracks);
  virtual Int_t FineProcess(TClonesArray& tracks);

  virtual Int_t DecodeADC( const THaEvData&, MOLLERElement *blk,
                           THaDetMap::Module *d, Int_t chan, Bool_t IsRef );

  virtual MOLLERElement* MakeElement(Double_t x, Double_t y, Double_t z,
                                     Int_t row, Int_t col, Int_t layer,
                                     Int_t id = 0);

  Double_t SizeRow() const { return fSizeRow; };
  Double_t SizeCol() const { return fSizeCol; };

protected:

  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );
  virtual Int_t FindGoodHit(MOLLERElement *);

  MOLLERTHcSampData SetSampIntTimePedestalPeakStandalone(
      const std::vector<Double_t>& samples,
      Int_t fNPedestalSamples,
      Double_t fSampThreshold,
      Int_t fNSB,
      Int_t fNSA,
      Int_t fNSAT,
      Int_t fMaxNPulses ) const;

  void InitTHcStandaloneTree();

  void FillTHcStandaloneTree(
      const THaEvData& evdata,
      MOLLERElement* blk,
      THaDetMap::Module* d,
      Int_t chan,
      const std::vector<Double_t>& samples,
      const MOLLERTHcSampData& data);

  // Configuration
  Int_t fNRefElem;
  Int_t fNrows;
  std::vector<Int_t> fNcols;
  Double_t fSizeRow;
  Double_t fSizeCol;
  Int_t fNcolsMax;
  Int_t fNlayers;
  MOLLERModeADC::Mode fModeADC;
  Bool_t fDisableRefADC;
  Bool_t fStoreEmptyElements;
  Bool_t fIsMC;
  Int_t fF1_RollOver;
  Int_t fF1_TimeWindow;

  // Mapping
  UShort_t fChanMapStart;
  std::vector<std::vector<Int_t> > fChanMap;
  std::vector<std::vector<Int_t> > fRefChanMap;
  std::vector<Bool_t> fModuleRefTimeFlag;
  std::vector<Int_t> fRefChanLo;
  std::vector<Int_t> fRefChanHi;

  // Output variable containers
  MOLLERGenericOutputData fGood;
  MOLLERGenericOutputData fRaw;
  MOLLERGenericOutputData fRefGood;
  MOLLERGenericOutputData fRefRaw;

  // Standalone hcana-style waveform diagnostic output
  TFile* fTHcStandaloneFile;
  TTree* fTHcStandaloneTree;

  Int_t fSA_evnum;
  Int_t fSA_elemID;
  Int_t fSA_crate;
  Int_t fSA_slot;
  Int_t fSA_chan;

  Double_t fSA_ped;
  Int_t fSA_a_mult;

  std::vector<Double_t> fSA_a;
  std::vector<Double_t> fSA_a_amp;
  std::vector<Double_t> fSA_a_time;
  std::vector<Double_t> fSA_tot;
  std::vector<Int_t> fSA_pileup;
  std::vector<Double_t> fSA_samps;

  std::vector<MOLLERElement*> fElements;
  std::vector<MOLLERElement*> fRefElements;
  std::vector<std::vector<std::vector<MOLLERElement*> > > fElementGrid;

  Bool_t fCoarseProcessed;
  Bool_t fFineProcessed;

  Double_t fConst;
  Double_t fSlope;
  Double_t fAccCharge;

  Int_t fNhits;
  Int_t fNRefhits;
  Int_t fNGoodADChits;

  Bool_t fStoreRawHits;

  TFile *fMultiPeakFile;
TTree *fMultiPeakTree;

Int_t fMulti_evnum;
Int_t fMulti_elemID;
Int_t fMulti_npulses;

std::vector<Double_t> fMulti_wf_raw;
std::vector<Double_t> fMulti_wf_cal;
std::vector<Double_t> fMulti_pulse_time;
std::vector<Double_t> fMulti_pulse_amp;
std::vector<Double_t> fMulti_pulse_int;

  /* TFile *fMultiPeakFile;
TTree *fMultiPeakTree;

Int_t fMulti_evnum;
Int_t fMulti_elemID;
Int_t fMulti_npulses;

std::vector<Double_t> fMulti_wf_raw;
std::vector<Double_t> fMulti_wf_cal;
std::vector<Double_t> fMulti_pulse_time;
std::vector<Double_t> fMulti_pulse_amp;
std::vector<Double_t> fMulti_pulse_int;
*/
  //  void InitMultiPeakTree();

private:
  void ClearOutputVariables();

  ClassDef(MOLLERGenericDetector,0)
};

#endif // MOLLERGenericDetector_h
