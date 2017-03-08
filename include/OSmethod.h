#ifndef OSmethod_h
#define OSmethod_h

// C++
#include <iostream>
#include <fstream>
#include <iomanip> // For setprecision
#include <vector>
#include <map>

// ROOT
#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TSystem.h"

// Include classes
#include "Tree.h"
#include "Settings.h"
#include "Category.h"
#include "FinalStates.h"
#include "FakeRates.h"
#include "bitops.h"

using namespace std;

const int num_of_processes         = Settings::num_of_processes;
const int num_of_flavours          = Settings::num_of_flavours;
const int num_of_final_states      = Settings::num_of_final_states;
const int num_of_categories        = Settings::num_of_categories;
const int num_of_eta_bins          = Settings::num_of_eta_bins;

class OSmethod: public Tree
{

public:
	
	OSmethod();
	~OSmethod();
   
   void FillHistos( TString );
   void MakeHistogramsZX( TString, TString );
   void DeclareHistos();
   void SaveHistos( TString );
   void GetHistos( TString );
   void SubtractWZ( bool );
   void ProduceFakeRates( TString );
   void RemoveNegativeBins( TH2F* );
   void Set_pT_binning( int, float* );
   void SetLumi( float );
   int find_current_process( TString );
   int FindFinalState();
   int FindFinalStateZX();

   
private:

   TFile *input_file, *input_file_data;
   TFile *fOutHistos;
   TTree *input_tree, *input_tree_data;

   TH1F *hCounters;
   
   Long64_t n_gen_events;
   
   vector<string> _s_process, _s_flavour;
   TString _histo_name;
   
   float jetPt[99];
   float jetEta[99];
   float jetPhi[99];
   float jetMass[99];
   float jetQGL[99];
   float jetPgOverPq[99];
   
   float _pT_bins[99];
   
   int _current_process, _current_final_state, _current_category, _n_pT_bins;
   float _lumi, _yield_SR;
   double gen_sum_weights, _event_weight;
   vector< vector <float> > _expected_yield_SR, _number_of_events_CR;

   
   TH2F *passing[num_of_processes][num_of_flavours], *failing[num_of_processes][num_of_flavours];
   
   TGraphErrors *FR_OS_electron_EB, *FR_OS_electron_EE, *FR_OS_muon_EB, *FR_OS_muon_EE;
   
   vector<Float_t> vector_X[num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_Y[num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_EX[num_of_eta_bins][num_of_flavours];
   vector<Float_t> vector_EY[num_of_eta_bins][num_of_flavours];
   
};
#endif
