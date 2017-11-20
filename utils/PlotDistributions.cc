// Standalone script to plot distributions in control regions
// Running with: root -l -b utils/PlotDistributions.cc++()

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "TROOT.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "THStack.h"
#include "TTree.h"
#include "TGraphErrors.h"

using namespace std;

enum LeptonFlavours
{
	NoMatch = 0,
	Conversion = 1,
	Lepton = 2,
	LightJet = 3,
	HeavyJet = 4,
	NUM_OF_FLAVOURS
};

float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};
int _n_pT_bins = 7;
TH1D* h_LepPT[LeptonFlavours::NUM_OF_FLAVOURS][2][2][2];
TH1D* h_LepMissingHit[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
TH1D* h_JetbTagger[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
TH1D* h_LepSIP[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
TH1D* h_LepBDT[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
TH1D* h_LepCounter[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
int n_events[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
int n_events_afterCuts[LeptonFlavours::NUM_OF_FLAVOURS][2];
vector <TString> _s_CR, _s_LepFlav, _s_EEorEB, _s_MatchFlavour;


int matchFlavour(int MatchJetPartonFlavour, int GenMCTruthMatchId, int GenMCTruthMatchMotherId)
{
	int flavour = 0;
	
	
	//no match to jets but match to photon
	if ( MatchJetPartonFlavour == 0 && (fabs(GenMCTruthMatchId) == 22 ||  fabs(GenMCTruthMatchMotherId) == 22) ) flavour = LeptonFlavours::Conversion;
	
	//no match to jets but match to lepton
	else if ( MatchJetPartonFlavour == 0 && (fabs(GenMCTruthMatchId) == 11 ||  fabs(GenMCTruthMatchId) == 13) ) flavour = LeptonFlavours::Lepton;
	
	//no match to jets
	else if ( MatchJetPartonFlavour == 0 ) flavour = LeptonFlavours::NoMatch;
	
	//match to light (u,d,s & g) jets
	else if ( fabs(MatchJetPartonFlavour) == 1 || fabs(MatchJetPartonFlavour) == 2 || fabs(MatchJetPartonFlavour) == 3 || fabs(MatchJetPartonFlavour) == 21) flavour = LeptonFlavours::LightJet;
	
	//match to heavy jets (c & b)
	else if ( fabs(MatchJetPartonFlavour) == 4 || fabs(MatchJetPartonFlavour) == 5 ) flavour = LeptonFlavours::HeavyJet;
	
	//match to gluon jets
	//else if ( fabs(MatchJetPartonFlavour) == 21 ) flavour = LeptonFlavours::GluonJet;
	
	return flavour;
}



void fillDistributionsCRZLL(TString input_file_name, bool isData = false)
{
	TFile* _file = TFile::Open(input_file_name);
	TTree* _Tree = (TTree*) _file->Get("CRZLLTree/candTree");
	int  nentries = _Tree->GetEntries();
	TH1F* hCounters = (TH1F*)_file->Get("CRZLLTree/Counters");
	Long64_t gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
	Float_t PFMET;
	Float_t ZZMass;
	Float_t Z1Mass;
	Short_t Z1Flav;
	Float_t Z2Mass;
	Short_t Z2Flav;
	vector<float>   *LepPt = 0;
	vector<float>   *LepEta = 0;
	vector<float>   *LepPhi = 0;
	vector<short>   *LepLepId = 0;
	vector<float>   *LepSIP = 0;
	vector<float>   *LepBDT = 0;
	vector<bool>    *LepisID = 0;
	vector<float>   *LepCombRelIsoPF = 0;
	vector<int>     *LepMissingHit = 0;
	Bool_t          passIsoPreFSR;
	vector<float>   *JetPt = 0;
	vector<float>   *JetEta = 0;
	vector<float>   *JetPhi = 0;
	vector<float>   *JetMass = 0;
	vector<float>   *JetBTagger = 0;
	vector<float>   *JetIsBtagged = 0;
	vector<float>   *JetIsBtaggedWithSF = 0;
	vector<float>   *JetQGLikelihood = 0;
	vector<float>   *JetAxis2 = 0;
	vector<float>   *JetMult = 0;
	vector<float>   *JetPtD = 0;
	vector<float>   *JetSigma = 0;
	vector<short>   *MatchJetPartonFlavour = 0;
	vector<float>   *MatchJetbTagger = 0;
	Float_t         genHEPMCweight;
	Float_t         PUWeight;
	Float_t         dataMCWeight;
	Float_t         trigEffWeight;
	Float_t         overallEventWeight;
	Float_t         HqTMCweight;
	Float_t         xsec;
	vector<short> 	 *GenMCTruthMatchId = 0;
	vector<short> 	 *GenMCTruthMatchMotherId = 0;
	
	_Tree->SetBranchAddress("PFMET",&PFMET);
	_Tree->SetBranchAddress("ZZMass",&ZZMass);
	_Tree->SetBranchAddress("Z1Mass",&Z1Mass);
	_Tree->SetBranchAddress("Z1Flav",&Z1Flav);
	_Tree->SetBranchAddress("Z2Mass",&Z2Mass);
	_Tree->SetBranchAddress("Z2Flav",&Z2Flav);
	_Tree->SetBranchAddress("LepPt",&LepPt);
	_Tree->SetBranchAddress("LepEta",&LepEta);
	_Tree->SetBranchAddress("LepPhi",&LepPhi);
	_Tree->SetBranchAddress("LepLepId",&LepLepId);
	_Tree->SetBranchAddress("LepSIP",&LepSIP);
	_Tree->SetBranchAddress("LepBDT",&LepBDT);
	_Tree->SetBranchAddress("LepisID",&LepisID);
	_Tree->SetBranchAddress("LepCombRelIsoPF",&LepCombRelIsoPF);
	_Tree->SetBranchAddress("LepMissingHit", &LepMissingHit);
	_Tree->SetBranchAddress("passIsoPreFSR",&passIsoPreFSR);
	_Tree->SetBranchAddress("JetPt",&JetPt);
	_Tree->SetBranchAddress("JetEta",&JetEta);
	_Tree->SetBranchAddress("JetPhi",&JetPhi);
	_Tree->SetBranchAddress("JetMass",&JetMass);
	_Tree->SetBranchAddress("JetBTagger",&JetBTagger);
	_Tree->SetBranchAddress("JetIsBtagged",&JetIsBtagged);
	_Tree->SetBranchAddress("JetIsBtaggedWithSF",&JetIsBtaggedWithSF);
	_Tree->SetBranchAddress("JetQGLikelihood",&JetQGLikelihood);
	_Tree->SetBranchAddress("JetAxis2",&JetAxis2);
	_Tree->SetBranchAddress("JetMult",&JetMult);
	_Tree->SetBranchAddress("JetPtD",&JetPtD);
	_Tree->SetBranchAddress("JetSigma",&JetSigma);
	_Tree->SetBranchAddress("MatchJetPartonFlavour",&MatchJetPartonFlavour);
	_Tree->SetBranchAddress("MatchJetbTagger",&MatchJetbTagger);
	_Tree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	_Tree->SetBranchAddress("PUWeight",&PUWeight);
	_Tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
	_Tree->SetBranchAddress("trigEffWeight",&trigEffWeight);
	_Tree->SetBranchAddress("overallEventWeight",&overallEventWeight);
	_Tree->SetBranchAddress("HqTMCweight",&HqTMCweight);
	_Tree->SetBranchAddress("xsec",&xsec);
	_Tree->SetBranchAddress("GenMCTruthMatchId",&GenMCTruthMatchId);
	_Tree->SetBranchAddress("GenMCTruthMatchMotherId",&GenMCTruthMatchMotherId);
	
	int jet1_category = -1;
	int lep1_flavour = -1;
	int lep1_EBorEE = -1;
	int jet2_category = -1;
	int lep2_flavour = -1;
	int lep2_EBorEE = -1;
	float _overallEventWeight;
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		for (int j = 0; j < 2 ; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				n_events[i_jf][j][k] = 0;
			}
		}
	}
	
	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		_overallEventWeight = 1.;
		if(!isData) _overallEventWeight = overallEventWeight;
		//Determine lepton origin
		jet1_category = 0;
		if(!isData) jet1_category = matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		lep1_flavour = (fabs(LepLepId->at(2)) == 11) ? 0 : 1;
		if(lep1_flavour == 0 && (abs(LepEta->at(2)) < 1.479) )  lep1_EBorEE = 0;
		if(lep1_flavour == 0 && (abs(LepEta->at(2)) >= 1.479) ) lep1_EBorEE = 1;
		if(lep1_flavour == 1 && (abs(LepEta->at(2)) < 1.2) )    lep1_EBorEE = 0;
		if(lep1_flavour == 1 && (abs(LepEta->at(2)) >= 1.2) )   lep1_EBorEE = 1;
		if(Z2Flav > 0) n_events[jet1_category][lep1_flavour][0]++;
		if(Z2Flav < 0) n_events[jet1_category][lep1_flavour][1]++;
		
		jet2_category = 0;
		if(!isData) jet2_category = matchFlavour(MatchJetPartonFlavour->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3));
		lep2_flavour = (fabs(LepLepId->at(3)) == 11) ? 0 : 1;
		if(lep2_flavour == 0 && (abs(LepEta->at(3)) < 1.479) )  lep2_EBorEE = 0;
		if(lep2_flavour == 0 && (abs(LepEta->at(3)) >= 1.479) ) lep2_EBorEE = 1;
		if(lep2_flavour == 1 && (abs(LepEta->at(3)) < 1.2) )    lep2_EBorEE = 0;
		if(lep2_flavour == 1 && (abs(LepEta->at(3)) >= 1.2) )   lep2_EBorEE = 1;
		if(Z2Flav > 0) n_events[jet2_category][lep2_flavour][0]++;
		if(Z2Flav < 0) n_events[jet2_category][lep2_flavour][1]++;
		if(Z2Flav > 0) h_LepPT[jet1_category][lep1_flavour][lep1_EBorEE][0]->Fill(LepPt->at(2), _overallEventWeight);
		if(Z2Flav > 0) h_LepPT[jet2_category][lep2_flavour][lep1_EBorEE][0]->Fill(LepPt->at(3), _overallEventWeight);
		if(Z2Flav < 0) h_LepPT[jet1_category][lep1_flavour][lep1_EBorEE][1]->Fill(LepPt->at(2), _overallEventWeight);
		if(Z2Flav < 0) h_LepPT[jet2_category][lep2_flavour][lep1_EBorEE][1]->Fill(LepPt->at(3), _overallEventWeight);
		if(Z2Flav > 0) h_LepMissingHit[lep1_flavour][lep1_EBorEE][0]->Fill(LepPt->at(2), LepMissingHit->at(2));
		if(Z2Flav > 0) h_LepMissingHit[lep2_flavour][lep2_EBorEE][0]->Fill(LepPt->at(3), LepMissingHit->at(3));
		if(Z2Flav < 0) h_LepMissingHit[lep1_flavour][lep1_EBorEE][1]->Fill(LepPt->at(2), LepMissingHit->at(2));
		if(Z2Flav < 0) h_LepMissingHit[lep2_flavour][lep2_EBorEE][1]->Fill(LepPt->at(3), LepMissingHit->at(3));
		if(Z2Flav > 0) h_LepCounter[lep1_flavour][lep1_EBorEE][0]->Fill(LepPt->at(2));
		if(Z2Flav > 0) h_LepCounter[lep2_flavour][lep2_EBorEE][0]->Fill(LepPt->at(3));
		if(Z2Flav < 0) h_LepCounter[lep1_flavour][lep1_EBorEE][1]->Fill(LepPt->at(2));
		if(Z2Flav < 0) h_LepCounter[lep2_flavour][lep2_EBorEE][1]->Fill(LepPt->at(3));
		
	}
	
	if(!isData)
	{
		cout << "==========================================================================================" << endl;
		cout << "[INFO] Control printout for " << input_file_name << endl;
		cout << "==========================================================================================" << endl;
		cout << "[INFO] Total number of electrons in SS ZLL CR = " << n_events[LeptonFlavours::Conversion][0][0] + n_events[LeptonFlavours::NoMatch][0][0] + n_events[LeptonFlavours::LightJet][0][0] + n_events[LeptonFlavours::HeavyJet][0][0] + n_events[LeptonFlavours::Lepton][0][0] << endl;
		cout << "[INFO] Number of electrons matched to leptons = " << n_events[LeptonFlavours::Lepton][0][0] << endl;
		cout << "[INFO] Number of electrons matched to photons = " << n_events[LeptonFlavours::Conversion][0][0] << endl;
		cout << "[INFO] Number of electrons not matched to jet = " << n_events[LeptonFlavours::NoMatch][0][0] << endl;
		cout << "[INFO] Number of electrons matched to light jet = " << n_events[LeptonFlavours::LightJet][0][0] << endl;
		cout << "[INFO] Number of electrons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][0][0] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of muons in SS ZLL CR = " << n_events[LeptonFlavours::Conversion][1][0] + n_events[LeptonFlavours::NoMatch][1][0] + n_events[LeptonFlavours::LightJet][1][0] + n_events[LeptonFlavours::HeavyJet][1][0] + n_events[LeptonFlavours::Lepton][1][0] << endl;
		cout << "[INFO] Number of muons matched to leptons = " << n_events[LeptonFlavours::Lepton][1][0] << endl;
		cout << "[INFO] Number of muons matched to photons = " << n_events[LeptonFlavours::Conversion][1][0] << endl;
		cout << "[INFO] Number of muons not matched to jet = " << n_events[LeptonFlavours::NoMatch][1][0] << endl;
		cout << "[INFO] Number of muons matched to light jet = " << n_events[LeptonFlavours::LightJet][1][0] << endl;
		cout << "[INFO] Number of muons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][1][0] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of electrons in OS ZLL CR = " << n_events[LeptonFlavours::Conversion][0][1] + n_events[LeptonFlavours::NoMatch][0][1] + n_events[LeptonFlavours::LightJet][0][1] + n_events[LeptonFlavours::HeavyJet][0][1] + n_events[LeptonFlavours::Lepton][0][1] << endl;
		cout << "[INFO] Number of electrons matched to leptons = " << n_events[LeptonFlavours::Lepton][0][1] << endl;
		cout << "[INFO] Number of electrons matched to photons = " << n_events[LeptonFlavours::Conversion][0][1] << endl;
		cout << "[INFO] Number of electrons not matched to jet = " << n_events[LeptonFlavours::NoMatch][0][1] << endl;
		cout << "[INFO] Number of electrons matched to light jet = " << n_events[LeptonFlavours::LightJet][0][1] << endl;
		cout << "[INFO] Number of electrons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][0][1] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of muons in OS ZLL CR = " << n_events[LeptonFlavours::Conversion][1][1] + n_events[LeptonFlavours::NoMatch][1][1] + n_events[LeptonFlavours::LightJet][1][1] + n_events[LeptonFlavours::HeavyJet][1][1] + n_events[LeptonFlavours::Lepton][1][1] << endl;
		cout << "[INFO] Number of muons matched to leptons = " << n_events[LeptonFlavours::Lepton][1][1] << endl;
		cout << "[INFO] Number of muons matched to photons = " << n_events[LeptonFlavours::Conversion][1][1] << endl;
		cout << "[INFO] Number of muons not matched to jet = " << n_events[LeptonFlavours::NoMatch][1][1] << endl;
		cout << "[INFO] Number of muons matched to light jet = " << n_events[LeptonFlavours::LightJet][1][1] << endl;
		cout << "[INFO] Number of muons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][1][1] << endl;
		cout << "============================================================" << endl;
	}
}


void fillDistributionsCRZL(TString input_file_name, bool isData=false)
{
	TFile* _file = TFile::Open(input_file_name);
	TTree* _Tree = (TTree*)_file->Get("CRZLTree/candTree");
	int  nentries = _Tree->GetEntries();
	TH1F* hCounters = (TH1F*)_file->Get("CRZLTree/Counters");
	Long64_t gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
	
	Float_t PFMET;
	Float_t Z1Mass;
	Short_t Z1Flav;
	vector<float>   *LepPt = 0;
	vector<float>   *LepEta = 0;
	vector<float>   *LepPhi = 0;
	vector<short>   *LepLepId = 0;
	vector<float>   *LepSIP = 0;
	vector<float>   *LepBDT = 0;
	vector<bool>    *LepisID = 0;
	vector<float>   *LepCombRelIsoPF = 0;
	vector<int>     *LepMissingHit = 0;
	Bool_t          passIsoPreFSR;
	vector<float>   *JetPt = 0;
	vector<float>   *JetEta = 0;
	vector<float>   *JetPhi = 0;
	vector<float>   *JetMass = 0;
	vector<float>   *JetBTagger = 0;
	vector<float>   *JetIsBtagged = 0;
	vector<float>   *JetIsBtaggedWithSF = 0;
	vector<float>   *JetQGLikelihood = 0;
	vector<float>   *JetAxis2 = 0;
	vector<float>   *JetMult = 0;
	vector<float>   *JetPtD = 0;
	vector<float>   *JetSigma = 0;
	vector<short>   *MatchJetPartonFlavour = 0;
	vector<float>   *MatchJetbTagger = 0;
	Float_t         genHEPMCweight;
	Float_t         PUWeight;
	Float_t         dataMCWeight;
	Float_t         trigEffWeight;
	Float_t         overallEventWeight;
	Float_t         HqTMCweight;
	Float_t         xsec;
	vector<short> 	 *GenMCTruthMatchId = 0;
	vector<short> 	 *GenMCTruthMatchMotherId = 0;
	
	_Tree->SetBranchAddress("PFMET",&PFMET);
	_Tree->SetBranchAddress("Z1Mass",&Z1Mass);
	_Tree->SetBranchAddress("Z1Flav",&Z1Flav);
	_Tree->SetBranchAddress("LepPt",&LepPt);
	_Tree->SetBranchAddress("LepEta",&LepEta);
	_Tree->SetBranchAddress("LepPhi",&LepPhi);
	_Tree->SetBranchAddress("LepLepId",&LepLepId);
	_Tree->SetBranchAddress("LepSIP",&LepSIP);
	_Tree->SetBranchAddress("LepBDT",&LepBDT);
	_Tree->SetBranchAddress("LepisID",&LepisID);
	_Tree->SetBranchAddress("LepCombRelIsoPF",&LepCombRelIsoPF);
	_Tree->SetBranchAddress("LepMissingHit", &LepMissingHit);
	_Tree->SetBranchAddress("passIsoPreFSR",&passIsoPreFSR);
	_Tree->SetBranchAddress("JetPt",&JetPt);
	_Tree->SetBranchAddress("JetEta",&JetEta);
	_Tree->SetBranchAddress("JetPhi",&JetPhi);
	_Tree->SetBranchAddress("JetMass",&JetMass);
	_Tree->SetBranchAddress("JetBTagger",&JetBTagger);
	_Tree->SetBranchAddress("JetIsBtagged",&JetIsBtagged);
	_Tree->SetBranchAddress("JetIsBtaggedWithSF",&JetIsBtaggedWithSF);
	_Tree->SetBranchAddress("JetQGLikelihood",&JetQGLikelihood);
	_Tree->SetBranchAddress("JetAxis2",&JetAxis2);
	_Tree->SetBranchAddress("JetMult",&JetMult);
	_Tree->SetBranchAddress("JetPtD",&JetPtD);
	_Tree->SetBranchAddress("JetSigma",&JetSigma);
	_Tree->SetBranchAddress("MatchJetPartonFlavour",&MatchJetPartonFlavour);
	_Tree->SetBranchAddress("MatchJetbTagger",&MatchJetbTagger);
	_Tree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	_Tree->SetBranchAddress("PUWeight",&PUWeight);
	_Tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
	_Tree->SetBranchAddress("trigEffWeight",&trigEffWeight);
	_Tree->SetBranchAddress("overallEventWeight",&overallEventWeight);
	_Tree->SetBranchAddress("HqTMCweight",&HqTMCweight);
	_Tree->SetBranchAddress("xsec",&xsec);
	_Tree->SetBranchAddress("GenMCTruthMatchId",&GenMCTruthMatchId);
	_Tree->SetBranchAddress("GenMCTruthMatchMotherId",&GenMCTruthMatchMotherId);
	
	int jet_category = -1;
	int lep_flavour = -1;
	int lep_EBorEE = -1;
	float _overallEventWeight;
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		n_events[i_jf][0][0] = 0;
		n_events_afterCuts[i_jf][0] = 0;
		n_events[i_jf][1][0] = 0;
		n_events_afterCuts[i_jf][1] = 0;
	}
	
	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		//Determine lepton origin
		jet_category = 0;
		if(!isData) jet_category = matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		lep_flavour = (fabs(LepLepId->at(2)) == 11) ? 0 : 1;
		if(lep_flavour == 0 && (abs(LepEta->at(2)) < 1.479) )  lep_EBorEE = 0;
		if(lep_flavour == 0 && (abs(LepEta->at(2)) >= 1.479) ) lep_EBorEE = 1;
		if(lep_flavour == 1 && (abs(LepEta->at(2)) < 1.2) )    lep_EBorEE = 0;
		if(lep_flavour == 1 && (abs(LepEta->at(2)) >= 1.2) )   lep_EBorEE = 1;
		
		n_events[jet_category][lep_flavour][0]++;
		
		// Calculate fake rates
		if ( Z1Mass < 40. ) {continue;}
		if ( Z1Mass > 120. ) {continue;}
		if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {continue;}
		if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {continue;}
		if ( LepSIP->at(2) > 4.) {continue;}
		if ( PFMET > 25. ) {continue;}
		
		// Fill distributions
		_overallEventWeight = 1.;
		//if (!input_file_name.Contains("Data")) _overallEventWeight = overallEventWeight;
		
		h_LepPT[jet_category][lep_flavour][lep_EBorEE][0]->Fill(LepPt->at(2), _overallEventWeight);
		if (!isData) h_JetbTagger[jet_category][lep_flavour][lep_EBorEE]->Fill(MatchJetbTagger->at(2), _overallEventWeight);
		h_LepBDT[jet_category][lep_flavour][lep_EBorEE]->Fill(LepBDT->at(2), _overallEventWeight);
		h_LepSIP[jet_category][lep_flavour][lep_EBorEE]->Fill(LepSIP->at(2), _overallEventWeight);
		h_LepMissingHit[lep_flavour][lep_EBorEE][0]->Fill(LepPt->at(2) ,LepMissingHit->at(2));
		h_LepCounter[lep_flavour][lep_EBorEE][0]->Fill(LepPt->at(2));
		
		n_events_afterCuts[jet_category][lep_flavour]++;
	}
	if(!isData)
	{
		cout << "==========================================================================================" << endl;
		cout << "[INFO] Control printout for " << input_file_name << endl;
		cout << "==========================================================================================" << endl;
		cout << "[INFO] Total number of electrons BEFORE cuts = " << n_events[LeptonFlavours::Conversion][0][0] + n_events[LeptonFlavours::NoMatch][0][0] + n_events[LeptonFlavours::LightJet][0][0] + n_events[LeptonFlavours::HeavyJet][0][0] + n_events[LeptonFlavours::Lepton][0][0] << endl;
		cout << "[INFO] Number of electrons matched to leptons = " << n_events[LeptonFlavours::Lepton][0][0] << endl;
		cout << "[INFO] Number of electrons matched to photons = " << n_events[LeptonFlavours::Conversion][0][0] << endl;
		cout << "[INFO] Number of electrons not matched to jet = " << n_events[LeptonFlavours::NoMatch][0][0] << endl;
		cout << "[INFO] Number of electrons matched to light jet = " << n_events[LeptonFlavours::LightJet][0][0] << endl;
		cout << "[INFO] Number of electrons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][0][0] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of muons BEFORE cuts = " << n_events[LeptonFlavours::Conversion][1][0] + n_events[LeptonFlavours::NoMatch][1][0] + n_events[LeptonFlavours::LightJet][1][0] + n_events[LeptonFlavours::HeavyJet][1][0] + n_events[LeptonFlavours::Lepton][1][0] << endl;
		cout << "[INFO] Number of muons matched to leptons = " << n_events[LeptonFlavours::Lepton][1][0] << endl;
		cout << "[INFO] Number of muons matched to photons = " << n_events[LeptonFlavours::Conversion][1][0] << endl;
		cout << "[INFO] Number of muons not matched to jet = " << n_events[LeptonFlavours::NoMatch][1][0] << endl;
		cout << "[INFO] Number of muons matched to light jet = " << n_events[LeptonFlavours::LightJet][1][0] << endl;
		cout << "[INFO] Number of muons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][1][0] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of electrons AFTER cuts = " << n_events_afterCuts[LeptonFlavours::Conversion][0] + n_events_afterCuts[LeptonFlavours::NoMatch][0] + n_events_afterCuts[LeptonFlavours::LightJet][0] + n_events_afterCuts[LeptonFlavours::HeavyJet][0] + n_events_afterCuts[LeptonFlavours::Lepton][0] << endl;
		cout << "[INFO] Number of electrons matched to leptons = " << n_events_afterCuts[LeptonFlavours::Lepton][0] << endl;
		cout << "[INFO] Number of electrons matched to photons = " << n_events_afterCuts[LeptonFlavours::Conversion][0] << endl;
		cout << "[INFO] Number of electrons not matched to jet = " << n_events_afterCuts[LeptonFlavours::NoMatch][0] << endl;
		cout << "[INFO] Number of electrons matched to light jet = " << n_events_afterCuts[LeptonFlavours::LightJet][0] << endl;
		cout << "[INFO] Number of electrons matched to heavy jet = " << n_events_afterCuts[LeptonFlavours::HeavyJet][0] << endl;
		cout << "============================================================" << endl;
		cout << "[INFO] Total number of muons AFTER cuts = " << n_events_afterCuts[LeptonFlavours::Conversion][1] + n_events_afterCuts[LeptonFlavours::NoMatch][1] + n_events_afterCuts[LeptonFlavours::LightJet][1] + n_events_afterCuts[LeptonFlavours::HeavyJet][1] + n_events_afterCuts[LeptonFlavours::Lepton][1] << endl;
		cout << "[INFO] Number of muons matched to leptons = " << n_events_afterCuts[LeptonFlavours::Lepton][1] << endl;
		cout << "[INFO] Number of muons matched to photons = " << n_events_afterCuts[LeptonFlavours::Conversion][1] << endl;
		cout << "[INFO] Number of muons not matched to jet = " << n_events_afterCuts[LeptonFlavours::NoMatch][1] << endl;
		cout << "[INFO] Number of muons matched to light jet = " << n_events_afterCuts[LeptonFlavours::LightJet][1] << endl;
		cout << "[INFO] Number of muons matched to heavy jet = " << n_events_afterCuts[LeptonFlavours::HeavyJet][1] << endl;
		cout << "============================================================" << endl;
	}
}

void Reset()
{
	for ( int i = 0; i < LeptonFlavours::NUM_OF_FLAVOURS ; i++)
	{
		for( int j = 0; j < 2; j++)
		{
			for ( int k = 0; k < 2; k++)
			{
				for ( int l = 0; l < 2; l++)
				{
					h_LepPT[i][j][k][l]->Reset();
					h_LepMissingHit[j][k][l]->Reset();
					h_LepCounter[j][k][l]->Reset();
				}
			}
		}
	}
}

void PlotDistributions()
{
	gROOT->SetBatch();
	
	_s_CR.push_back("SS");
	_s_CR.push_back("OS");
	
	_s_LepFlav.push_back("ele");
	_s_LepFlav.push_back("mu");
	
	_s_EEorEB.push_back("EB");
	_s_EEorEB.push_back("EE");
	
	_s_MatchFlavour.push_back("NoMatch");
	_s_MatchFlavour.push_back("Conversion");
	_s_MatchFlavour.push_back("Lepton");
	_s_MatchFlavour.push_back("LightJet");
	_s_MatchFlavour.push_back("HeavyJet");
	
	float lumi = 35.9;
	
	TString path = "NewData/";
	TString file_name = "/ZZ4lAnalysis.root";
	
	TString DY       = path + "DYJetsToLL_M50" + file_name;
	TString TTJets   = path + "TTJets_DiLept_ext1" + file_name;
	TString WZTo3LNu = path + "WZTo3LNu" + file_name;
	TString Data     = "Moriond_2017/Data/ZZ4lAnalysis.root";
	
	TString histo_name;
	
	histo_name ="h_LepMissingHit_SS_ele_EB";
	h_LepMissingHit[0][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepMissingHit_SS_ele_EE";
	h_LepMissingHit[0][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_SS_ele_EB";
	h_LepCounter[0][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_SS_ele_EE";
	h_LepCounter[0][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	
	histo_name ="h_LepMissingHit_OS_ele_EB";
	h_LepMissingHit[0][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepMissingHit_OS_ele_EE";
	h_LepMissingHit[0][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_OS_ele_EB";
	h_LepCounter[0][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_OS_ele_EE";
	h_LepCounter[0][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	
	histo_name ="h_LepMissingHit_SS_mu_EB";
	h_LepMissingHit[1][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepMissingHit_SS_mu_EE";
	h_LepMissingHit[1][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_SS_mu_EB";
	h_LepCounter[1][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_SS_mu_EE";
	h_LepCounter[1][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	
	histo_name ="h_LepMissingHit_OS_mu_EB";
	h_LepMissingHit[1][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepMissingHit_OS_mu_EE";
	h_LepMissingHit[1][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_OS_mu_EB";
	h_LepCounter[1][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	histo_name ="h_LepCounter_OS_mu_EE";
	h_LepCounter[1][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		histo_name ="h_LepPT_SS_ele_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][0][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_SS_ele_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][0][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_SS_mu_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][1][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_SS_mu_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][1][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		
		histo_name ="h_LepPT_OS_ele_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][0][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_OS_ele_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][0][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_OS_mu_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][1][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_OS_mu_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepPT[i_jf][1][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	}
	
	fillDistributionsCRZLL(DY);
	fillDistributionsCRZLL(TTJets);
	fillDistributionsCRZLL(WZTo3LNu);
	
	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	TString canvas_name;
	c1->cd();
	TH1D* sum_pT;
	
	for(int i_lep_flav = 0; i_lep_flav < 2; i_lep_flav++)
	{
		for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
		{
			for (int i_CR = 0; i_CR < 2; i_CR++)
			{
				sum_pT = (TH1D*)h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->Clone();
				sum_pT->Add(h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][i_CR]);
				sum_pT->Add(h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][i_CR]);
				sum_pT->Add(h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][i_CR]);
				sum_pT->Add(h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][i_CR]);
				
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->Divide(sum_pT);
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->SetLineColor(kBlue);
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->SetMarkerColor(kBlue);
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->SetMaximum(1.);
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->SetMinimum(0.);
				h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][i_CR]->Draw("HIST ");
				h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][i_CR]->Divide(sum_pT);
				h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][i_CR]->SetLineColor(kBlack);
				h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][i_CR]->SetMarkerColor(kBlack);
				h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][i_CR]->Draw("HIST SAME");
				h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][i_CR]->Divide(sum_pT);
				h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][i_CR]->SetLineColor(kRed);
				h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][i_CR]->SetMarkerColor(kRed);
				h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][i_CR]->Draw("HIST SAME");
				h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][i_CR]->Divide(sum_pT);
				h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][i_CR]->SetLineColor(kGreen);
				h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][i_CR]->SetMarkerColor(kGreen);
				h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][i_CR]->Draw("HIST SAME");
				h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][i_CR]->Divide(sum_pT);
				h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][i_CR]->SetLineColor(kViolet);
				h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][i_CR]->SetMarkerColor(kViolet);
				h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][i_CR]->Draw("HIST SAME");
				
				canvas_name = "./MCTruthStudy/CRZLL_"+_s_CR.at(i_CR)+"_PT_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".pdf";
				c1->SaveAs(canvas_name);
				canvas_name = "./MCTruthStudy/CRZLL_"+_s_CR.at(i_CR)+"_PT_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".png";
				c1->SaveAs(canvas_name);
				
				sum_pT->Reset();
			}
		}
	}
	
	TH1D *data_MissingHit_EB_SS, *data_MissingHit_EE_SS, *MC_MissingHit_EB_SS, *MC_MissingHit_EE_SS;
	TH1D *data_MissingHit_EB_OS, *data_MissingHit_EE_OS, *MC_MissingHit_EB_OS, *MC_MissingHit_EE_OS;
	
	h_LepMissingHit[0][0][0]->Divide(h_LepCounter[0][0][0]);
	h_LepMissingHit[0][1][0]->Divide(h_LepCounter[0][1][0]);
	h_LepMissingHit[0][0][1]->Divide(h_LepCounter[0][0][1]);
	h_LepMissingHit[0][1][1]->Divide(h_LepCounter[0][1][1]);
	
	MC_MissingHit_EB_SS = (TH1D*)h_LepMissingHit[0][0][0]->Clone();
	MC_MissingHit_EE_SS = (TH1D*)h_LepMissingHit[0][1][0]->Clone();
	MC_MissingHit_EB_OS = (TH1D*)h_LepMissingHit[0][0][1]->Clone();
	MC_MissingHit_EE_OS = (TH1D*)h_LepMissingHit[0][1][1]->Clone();
	
	Reset();
	fillDistributionsCRZLL(Data, true);
	
	h_LepMissingHit[0][0][0]->Divide(h_LepCounter[0][0][0]);
	h_LepMissingHit[0][1][0]->Divide(h_LepCounter[0][1][0]);
	h_LepMissingHit[0][0][1]->Divide(h_LepCounter[0][0][1]);
	h_LepMissingHit[0][1][1]->Divide(h_LepCounter[0][1][1]);
	
	data_MissingHit_EB_SS = (TH1D*)h_LepMissingHit[0][0][0]->Clone();
	data_MissingHit_EE_SS = (TH1D*)h_LepMissingHit[0][1][0]->Clone();
	data_MissingHit_EB_OS = (TH1D*)h_LepMissingHit[0][0][1]->Clone();
	data_MissingHit_EE_OS = (TH1D*)h_LepMissingHit[0][1][1]->Clone();
	
	MC_MissingHit_EB_SS->DrawNormalized("HIST");
	data_MissingHit_EB_SS->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZLL_SS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZLL_SS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".png";
	c1->SaveAs(canvas_name);
	MC_MissingHit_EE_SS->DrawNormalized("HIST");
	data_MissingHit_EE_SS->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZLL_SS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZLL_SS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".png";
	c1->SaveAs(canvas_name);
	
	MC_MissingHit_EB_OS->DrawNormalized("HIST");
	data_MissingHit_EB_OS->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZLL_OS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZLL_OS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".png";
	c1->SaveAs(canvas_name);
	MC_MissingHit_EE_OS->DrawNormalized("HIST");
	data_MissingHit_EE_OS->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZLL_OS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZLL_OS_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".png";
	c1->SaveAs(canvas_name);
	
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		histo_name ="h_LepbTag_ele_EE_"+_s_MatchFlavour.at(i_jf);
		h_JetbTagger[i_jf][0][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_ele_EB_"+_s_MatchFlavour.at(i_jf);
		h_JetbTagger[i_jf][0][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_mu_EE_"+_s_MatchFlavour.at(i_jf);
		h_JetbTagger[i_jf][1][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_mu_EB_"+_s_MatchFlavour.at(i_jf);
		h_JetbTagger[i_jf][1][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		
		histo_name ="h_LepSIP_ele_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepSIP[i_jf][0][0] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_ele_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepSIP[i_jf][0][1] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_mu_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepSIP[i_jf][1][0] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_mu_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepSIP[i_jf][1][1] = new TH1D(histo_name,histo_name,40., 0.,4.);
		
		histo_name ="h_LepBDT_ele_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepBDT[i_jf][0][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_ele_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepBDT[i_jf][0][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_mu_EE_"+_s_MatchFlavour.at(i_jf);
		h_LepBDT[i_jf][1][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_mu_EB_"+_s_MatchFlavour.at(i_jf);
		h_LepBDT[i_jf][1][1] = new TH1D(histo_name,histo_name,20, -1., 1.);

		
	}
	
	Reset();
	fillDistributionsCRZL(DY);
	fillDistributionsCRZL(TTJets);
	fillDistributionsCRZL(WZTo3LNu);

	c1->cd();
	for(int i_lep_flav = 0; i_lep_flav < 2; i_lep_flav++)
	{
		for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
		{
			sum_pT = (TH1D*)h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->Clone();
			sum_pT->Add(h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][0]);
			sum_pT->Add(h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][0]);
			sum_pT->Add(h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][0]);
			sum_pT->Add(h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][0]);
			
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->SetLineColor(kBlue);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->SetMarkerColor(kBlue);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->SetMaximum(1.0);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->SetMinimum(-0.01);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE][0]->Draw("HIST ");
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][0]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][0]->SetLineColor(kBlack);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][0]->SetMarkerColor(kBlack);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE][0]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][0]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][0]->SetLineColor(kRed);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][0]->SetMarkerColor(kRed);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE][0]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][0]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][0]->SetLineColor(kGreen);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][0]->SetMarkerColor(kGreen);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE][0]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][0]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][0]->SetLineColor(kViolet);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][0]->SetMarkerColor(kViolet);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE][0]->Draw("HIST SAME");
			
			canvas_name = "./MCTruthStudy/CRZL_PT_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/CRZL_PT_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".png";
			c1->SaveAs(canvas_name);
			
			sum_pT->Reset();
		}
	}
	
	for(int i_lep_flav = 0; i_lep_flav < 2; i_lep_flav++)
	{
		for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
		{
			h_LepSIP[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetLineColor(kBlue);
			h_LepSIP[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kBlue);
			h_LepSIP[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->DrawNormalized("HIST");
			h_LepSIP[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->SetLineColor(kBlack);
			h_LepSIP[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->SetMarkerColor(kBlack);
			h_LepSIP[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->DrawNormalized("HIST SAME");
			h_LepSIP[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetLineColor(kRed);
			h_LepSIP[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kRed);
			h_LepSIP[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->DrawNormalized("HIST SAME");
			h_LepSIP[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->SetLineColor(kGreen);
			h_LepSIP[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->SetMarkerColor(kGreen);
			h_LepSIP[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->DrawNormalized("HIST SAME");
			h_LepSIP[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->SetLineColor(kViolet);
			h_LepSIP[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->SetMarkerColor(kViolet);
			h_LepSIP[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->DrawNormalized("HIST SAME");
			
			canvas_name = "./MCTruthStudy/CRZL_SIP_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/CRZL_SIP_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".png";
			c1->SaveAs(canvas_name);
		}
	}
	
	for(int i_lep_flav = 0; i_lep_flav < 2; i_lep_flav++)
	{
		for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
		{
			h_JetbTagger[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetLineColor(kBlue);
			h_JetbTagger[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kBlue);
			h_JetbTagger[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->DrawNormalized("HIST");
			h_JetbTagger[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetLineColor(kRed);
			h_JetbTagger[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kRed);
			h_JetbTagger[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->DrawNormalized("HIST SAME");
			
			canvas_name = "./MCTruthStudy/CRZL_bTaggScore_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/CRZL_bTaggScore_distribution_" + _s_LepFlav.at(i_lep_flav) + "_" + _s_EEorEB.at(i_EBorEE) + ".png";
			c1->SaveAs(canvas_name);
		}
	}
	
	for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
	{
		h_LepBDT[LeptonFlavours::LightJet][0][i_EBorEE]->SetLineColor(kBlue);
		h_LepBDT[LeptonFlavours::LightJet][0][i_EBorEE]->SetMarkerColor(kBlue);
		h_LepBDT[LeptonFlavours::LightJet][0][i_EBorEE]->DrawNormalized("HIST");
		h_LepBDT[LeptonFlavours::NoMatch][0][i_EBorEE]->SetLineColor(kBlack);
		h_LepBDT[LeptonFlavours::NoMatch][0][i_EBorEE]->SetMarkerColor(kBlack);
		h_LepBDT[LeptonFlavours::NoMatch][0][i_EBorEE]->DrawNormalized("HIST SAME");
		h_LepBDT[LeptonFlavours::HeavyJet][0][i_EBorEE]->SetLineColor(kRed);
		h_LepBDT[LeptonFlavours::HeavyJet][0][i_EBorEE]->SetMarkerColor(kRed);
		h_LepBDT[LeptonFlavours::HeavyJet][0][i_EBorEE]->DrawNormalized("HIST SAME");
		h_LepBDT[LeptonFlavours::Conversion][0][i_EBorEE]->SetLineColor(kGreen);
		h_LepBDT[LeptonFlavours::Conversion][0][i_EBorEE]->SetMarkerColor(kGreen);
		h_LepBDT[LeptonFlavours::Conversion][0][i_EBorEE]->DrawNormalized("HIST SAME");
		h_LepBDT[LeptonFlavours::Lepton][0][i_EBorEE]->SetLineColor(kViolet);
		h_LepBDT[LeptonFlavours::Lepton][0][i_EBorEE]->SetMarkerColor(kViolet);
		h_LepBDT[LeptonFlavours::Lepton][0][i_EBorEE]->DrawNormalized("HIST SAME");
		
		canvas_name = "./MCTruthStudy/CRZL_BDT_distribution_ele_" + _s_EEorEB.at(i_EBorEE) + ".pdf";
		c1->SaveAs(canvas_name);
		canvas_name = "./MCTruthStudy/CRZL_BDT_distribution_ele_" + _s_EEorEB.at(i_EBorEE) + ".png";
		c1->SaveAs(canvas_name);
	}
	TH1D *data_MissingHit_EB, *data_MissingHit_EE, *MC_MissingHit_EB, *MC_MissingHit_EE;
	
	h_LepMissingHit[0][0][0]->Divide(h_LepCounter[0][0][0]);
	h_LepMissingHit[0][1][0]->Divide(h_LepCounter[0][1][0]);
	h_LepMissingHit[1][0][0]->Divide(h_LepCounter[1][0][0]);
	h_LepMissingHit[1][1][0]->Divide(h_LepCounter[1][1][0]);
	
	MC_MissingHit_EB = (TH1D*)h_LepMissingHit[0][0][0]->Clone();
	MC_MissingHit_EE = (TH1D*)h_LepMissingHit[0][1][0]->Clone();

	Reset();
	fillDistributionsCRZL(Data, true);
	
	h_LepMissingHit[0][0][0]->Divide(h_LepCounter[0][0][0]);
	h_LepMissingHit[0][1][0]->Divide(h_LepCounter[0][1][0]);
	h_LepMissingHit[1][0][0]->Divide(h_LepCounter[1][0][0]);
	h_LepMissingHit[1][1][0]->Divide(h_LepCounter[1][1][0]);
	
	data_MissingHit_EB = (TH1D*)h_LepMissingHit[0][0][0]->Clone();
	data_MissingHit_EE = (TH1D*)h_LepMissingHit[0][1][0]->Clone();
	
	MC_MissingHit_EB->DrawNormalized("HIST");
	data_MissingHit_EB->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZL_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZL_avgMissingHit_distribution_ele_" + _s_EEorEB.at(0) + ".png";
	c1->SaveAs(canvas_name);
	MC_MissingHit_EE->DrawNormalized("HIST");
	data_MissingHit_EE->DrawNormalized("P SAME");
	canvas_name = "./MCTruthStudy/CRZL_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".pdf";
	c1->SaveAs(canvas_name);
	canvas_name = "./MCTruthStudy/CRZL_avgMissingHit_distribution_ele_" + _s_EEorEB.at(1) + ".png";
	c1->SaveAs(canvas_name);
}
