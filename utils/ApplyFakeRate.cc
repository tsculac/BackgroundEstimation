// Standalone script to apply simple FR coming from MCTruthAnalyzer on MC Z+LL control region to estimate the Z+X yields
// Running with: root -l -b utils/ApplyFakeRate.cc++()

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
TGraphErrors *fr_ele_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_ele_EE[LeptonFlavours::NUM_OF_FLAVOURS];
TGraphErrors *fr_mu_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_mu_EE[LeptonFlavours::NUM_OF_FLAVOURS];
TGraphErrors *g_FR_mu_EB,*g_FR_mu_EE,*g_FR_e_EB,*g_FR_e_EE;
TH1D* h_LepPT[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
int n_events[LeptonFlavours::NUM_OF_FLAVOURS][2];

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

float GetFakeRate_WithFlavour(float lep_Pt, float lep_eta, int lep_ID, int flavour)
{
	float fr = -1.;
	
	float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
	int   my_lep_ID = abs(lep_ID);
	
	int bin = 0;
	if ( my_lep_Pt > 5 && my_lep_Pt <= 7 ) bin = 0;
	else if ( my_lep_Pt >  7 && my_lep_Pt <= 10 ) bin = 1;
	else if ( my_lep_Pt > 10 && my_lep_Pt <= 20 ) bin = 2;
	else if ( my_lep_Pt > 20 && my_lep_Pt <= 30 ) bin = 3;
	else if ( my_lep_Pt > 30 && my_lep_Pt <= 40 ) bin = 4;
	else if ( my_lep_Pt > 40 && my_lep_Pt <= 50 ) bin = 5;
	else if ( my_lep_Pt > 50 && my_lep_Pt <= 80 ) bin = 6;


	if ( fabs(my_lep_ID) == 11 ) bin = bin-1; // there is no [5, 7] bin in the electron fake rate
	
	if ( my_lep_ID == 11 )
	{
		if ( fabs(lep_eta) < 1.479 )
			fr = (fr_ele_EB[flavour]->GetY())[bin];
		else
			fr = (fr_ele_EE[flavour]->GetY())[bin];
	}
	else if ( my_lep_ID == 13 )
	{
		if ( fabs(lep_eta) < 1.2 )
			fr = (fr_mu_EB[flavour]->GetY())[bin];
		else
			fr = (fr_mu_EE[flavour]->GetY())[bin];
	}
	else
	{
		cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
		fr = 0;
	}
	
	return fr;
}

float GetFakeRate(float lep_Pt, float lep_eta, int lep_ID)
{
	float my_lep_Pt = lep_Pt >= 80. ? 79. : lep_Pt;
	int   my_lep_ID = abs(lep_ID);
	
	int bin = 0;
	if ( my_lep_Pt > 5 && my_lep_Pt <= 7 ) bin = 0;
	else if ( my_lep_Pt >  7 && my_lep_Pt <= 10 ) bin = 1;
	else if ( my_lep_Pt > 10 && my_lep_Pt <= 20 ) bin = 2;
	else if ( my_lep_Pt > 20 && my_lep_Pt <= 30 ) bin = 3;
	else if ( my_lep_Pt > 30 && my_lep_Pt <= 40 ) bin = 4;
	else if ( my_lep_Pt > 40 && my_lep_Pt <= 50 ) bin = 5;
	else if ( my_lep_Pt > 50 && my_lep_Pt <= 80 ) bin = 6;
	
	if ( fabs(my_lep_ID) == 11 ) bin = bin-1; // there is no [5, 7] bin in the electron fake rate
	
	if ( my_lep_ID == 11 )
	{
		if ( fabs(lep_eta) < 1.479 )
			return (g_FR_e_EB->GetY())[bin];
		else
			return (g_FR_e_EE->GetY())[bin];
	}
	else if ( my_lep_ID == 13 )
	{
		if ( fabs(lep_eta) < 1.2 )
			return (g_FR_mu_EB->GetY())[bin];
		else
			return (g_FR_mu_EE->GetY())[bin];
	}
	else
	{
		cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
		return 0;
	}
}

void fillZXHisto_WithFlavourInfo(TString input_file_name, float lumi, TH1F* m4l_ZX[3][4])
{
	vector<float> _fs_ROS_SS;
//	_fs_ROS_SS.push_back(1.22);//4mu
//	_fs_ROS_SS.push_back(0.97);//4e
//	_fs_ROS_SS.push_back(1.30);//2e2mu
//	_fs_ROS_SS.push_back(0.98);//2mu2e
	
	_fs_ROS_SS.push_back(1.);//4mu
	_fs_ROS_SS.push_back(1.);//4e
	_fs_ROS_SS.push_back(1.);//2e2mu
	_fs_ROS_SS.push_back(1.);//2mu2e
	
	TFile *fr_file = TFile::Open("output.root");

	TString fr_name;
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		fr_name = "FR_ele_EB_"+to_string(i_jf);
		fr_ele_EB[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_ele_EE_"+to_string(i_jf);
		fr_ele_EE[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EB_"+to_string(i_jf);
		fr_mu_EB[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EE_"+to_string(i_jf);
		fr_mu_EE[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);

	}
	
	
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
	int jet2_category = -1;
	int _current_final_state = -1;
	int _current_proces = -1;

	if(input_file_name.Contains("DY"))     _current_proces = 0;
	if(input_file_name.Contains("TTJets")) _current_proces = 1;
	if(input_file_name.Contains("WZ"))     _current_proces = 2;
	

	float _event_weight;
	float _yield_SR;

	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		if(Z2Flav < 0) continue;

		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 169) _current_final_state = 0;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 121) _current_final_state = 1;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 169) _current_final_state = 2;
		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 121) _current_final_state = 3;
		
		//Determine lepton origin
		jet1_category = matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		
		jet2_category = matchFlavour(MatchJetPartonFlavour->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3));
		
		_event_weight = (lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;
		
		// Calculate yield
		_yield_SR = _fs_ROS_SS.at(_current_final_state)*GetFakeRate_WithFlavour(LepPt->at(2),LepEta->at(2),LepLepId->at(2),jet1_category)*GetFakeRate_WithFlavour(LepPt->at(3),LepEta->at(3),LepLepId->at(3),jet2_category);
		m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _yield_SR*_event_weight);
		
//		if (m4l_ZX[0][2]->Integral() != m4l_ZX[0][2]->Integral())
//		{
//			cout << "YIELD = " << _yield_SR << endl;
//			cout << "WEIGHT = " << _event_weight << endl;
//			cout << "ZZMass = " << ZZMass << endl;
//			cout << LepPt->at(2) << " " << LepPt->at(3) << endl;
//			cout << LepEta->at(2) << " " << LepEta->at(3) << endl;
//			cout << LepLepId->at(2) << " " << LepLepId->at(3) << endl;
//			cout << jet1_category << " " << jet2_category << endl;
//			cout << GetFakeRate_WithFlavour(LepPt->at(2),LepEta->at(2),LepLepId->at(2),jet1_category) << " " << GetFakeRate_WithFlavour(LepPt->at(3),LepEta->at(3),LepLepId->at(3),jet2_category) << endl;
//			cout << "==================" << endl;
//			break;
//		}

	}

cout << "[INFO] Z+X histograms filled for file " + input_file_name + "!" << endl;
}


void fillZXHisto(TString input_file_name, float lumi, TH1F* m4l_ZX[3][4])
{
	vector<float> _fs_ROS_SS;
	//	_fs_ROS_SS.push_back(1.22);//4mu
	//	_fs_ROS_SS.push_back(0.97);//4e
	//	_fs_ROS_SS.push_back(1.30);//2e2mu
	//	_fs_ROS_SS.push_back(0.98);//2mu2e
	
	_fs_ROS_SS.push_back(1.);//4mu
	_fs_ROS_SS.push_back(1.);//4e
	_fs_ROS_SS.push_back(1.);//2e2mu
	_fs_ROS_SS.push_back(1.);//2mu2e
	
	TFile *fr_file = TFile::Open("FakeRates_SS_Moriond17.root");
	
	g_FR_mu_EB = (TGraphErrors*)fr_file->Get("FR_SS_muon_EB");
	g_FR_mu_EE = (TGraphErrors*)fr_file->Get("FR_SS_muon_EE");
	g_FR_e_EB  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EB");
	g_FR_e_EE  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EE");
	
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
	
	int _current_final_state = -1;
	int _current_proces = -1;
	
	if(input_file_name.Contains("DY"))     _current_proces = 0;
	if(input_file_name.Contains("TTJets")) _current_proces = 1;
	if(input_file_name.Contains("WZ"))     _current_proces = 2;
	
	float _event_weight;
	float _yield_SR;
	
	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		if(Z2Flav < 0) continue;
		
		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 169) _current_final_state = 0;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 121) _current_final_state = 1;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 169) _current_final_state = 2;
		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 121) _current_final_state = 3;
		
		_event_weight = (lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;
		
		// Calculate yield
		_yield_SR = _fs_ROS_SS.at(_current_final_state)*GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2))*GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
		//		cout << "YIELD = " << _yield_SR << endl;
		//		cout << "WEIGHT = " << _event_weight << endl;
		//		cout << "==================" << endl;
		m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _yield_SR*_event_weight);
		
		
	}
	
	cout << "[INFO] Z+X histograms filled for file " + input_file_name + "!" << endl;
}

void fillZXHisto_CutBased(TString input_file_name, float lumi, TH1F* m4l_ZX[3][4])
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
	
	int _current_final_state = -1;
	int _current_proces = -1;
	
	if(input_file_name.Contains("DY"))     _current_proces = 0;
	if(input_file_name.Contains("TTJets")) _current_proces = 1;
	if(input_file_name.Contains("WZ"))     _current_proces = 2;
	
	float _event_weight;
	float _yield_SR;
	
	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		if(Z2Flav < 0) continue;
		
		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 169) _current_final_state = 0;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 121) _current_final_state = 1;
		if(abs(Z1Flav) == 121 && abs(Z2Flav) == 169) _current_final_state = 2;
		if(abs(Z1Flav) == 169 && abs(Z2Flav) == 121) _current_final_state = 3;
		
		_event_weight = (lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;
		
		// Calculate yield
		if(LepisID->at(2) && LepisID->at(3) && LepCombRelIsoPF->at(2) < 0.35 && LepCombRelIsoPF->at(3) < 0.35) m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _event_weight);
		
		
	}
	
	cout << "[INFO] Z+X histograms filled for file " + input_file_name + "!" << endl;
}

void fillDistributions(TString input_file_name)
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
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		n_events[i_jf][0] = 0;
		n_events[i_jf][1] = 0;
	}
	
	for(int i_event = 0; i_event < nentries; i_event++)
	{
		_Tree->GetEvent(i_event);
		
		if(Z2Flav < 0) continue;
		
		//Determine lepton origin
		jet1_category = matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		lep1_flavour = (fabs(LepLepId->at(2)) == 11) ? 0 : 1;
		if(lep1_flavour == 0 && (abs(LepEta->at(2)) < 1.479) )  lep1_EBorEE = 0;
		if(lep1_flavour == 0 && (abs(LepEta->at(2)) >= 1.479) ) lep1_EBorEE = 1;
		if(lep1_flavour == 1 && (abs(LepEta->at(2)) < 1.2) )    lep1_EBorEE = 0;
		if(lep1_flavour == 1 && (abs(LepEta->at(2)) >= 1.2) )   lep1_EBorEE = 1;
		n_events[jet1_category][lep1_flavour]++;
		
		jet2_category = matchFlavour(MatchJetPartonFlavour->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3));
		lep2_flavour = (fabs(LepLepId->at(3)) == 11) ? 0 : 1;
		if(lep2_flavour == 0 && (abs(LepEta->at(3)) < 1.479) )  lep2_EBorEE = 0;
		if(lep2_flavour == 0 && (abs(LepEta->at(3)) >= 1.479) ) lep2_EBorEE = 1;
		if(lep2_flavour == 1 && (abs(LepEta->at(3)) < 1.2) )    lep2_EBorEE = 0;
		if(lep2_flavour == 1 && (abs(LepEta->at(3)) >= 1.2) )   lep2_EBorEE = 1;
		n_events[jet2_category][lep2_flavour]++;
		
		h_LepPT[jet1_category][lep1_flavour][lep1_EBorEE]->Fill(LepPt->at(2), overallEventWeight);
		h_LepPT[jet2_category][lep2_flavour][lep1_EBorEE]->Fill(LepPt->at(3), overallEventWeight);
	}
	
	cout << "==========================================================================================" << endl;
	cout << "[INFO] Control printout for " << input_file_name << endl;
	cout << "==========================================================================================" << endl;
	cout << "[INFO] Total number of electrons in ZLL CR = " << n_events[LeptonFlavours::Conversion][0] + n_events[LeptonFlavours::NoMatch][0] + n_events[LeptonFlavours::LightJet][0] + n_events[LeptonFlavours::HeavyJet][0] + n_events[LeptonFlavours::Lepton][0] << endl;
	cout << "[INFO] Number of electrons matched to leptons = " << n_events[LeptonFlavours::Lepton][0] << endl;
	cout << "[INFO] Number of electrons matched to photons = " << n_events[LeptonFlavours::Conversion][0] << endl;
	cout << "[INFO] Number of electrons not matched to jet = " << n_events[LeptonFlavours::NoMatch][0] << endl;
	cout << "[INFO] Number of electrons matched to light jet = " << n_events[LeptonFlavours::LightJet][0] << endl;
	cout << "[INFO] Number of electrons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][0] << endl;
	cout << "============================================================" << endl;
	cout << "[INFO] Total number of muons in ZLL CR = " << n_events[LeptonFlavours::Conversion][1] + n_events[LeptonFlavours::NoMatch][1] + n_events[LeptonFlavours::LightJet][1] + n_events[LeptonFlavours::HeavyJet][1] + n_events[LeptonFlavours::Lepton][1] << endl;
	cout << "[INFO] Number of muons matched to leptons = " << n_events[LeptonFlavours::Lepton][1] << endl;
	cout << "[INFO] Number of muons matched to photons = " << n_events[LeptonFlavours::Conversion][1] << endl;
	cout << "[INFO] Number of muons not matched to jet = " << n_events[LeptonFlavours::NoMatch][1] << endl;
	cout << "[INFO] Number of muons matched to light jet = " << n_events[LeptonFlavours::LightJet][1] << endl;
	cout << "[INFO] Number of muons matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet][1] << endl;
	cout << "============================================================" << endl;
}

void DrawHistos(TH1F* m4l_ZX[3][4], TString FRType)
{
	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	THStack *hs_4mu = new THStack("hs_4mu","test stacked histograms");
	c1->cd();
	m4l_ZX[0][0]->SetFillColor(kGreen);
	m4l_ZX[1][0]->SetFillColor(kBlue);
	m4l_ZX[2][0]->SetFillColor(kViolet);
	hs_4mu->Add(m4l_ZX[2][0]);
	hs_4mu->Add(m4l_ZX[1][0]);
	hs_4mu->Add(m4l_ZX[0][0]);
	hs_4mu->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4mu.png");
	
	c1->cd();
	THStack *hs_4e = new THStack("hs_4e","test stacked histograms");
	m4l_ZX[0][1]->SetFillColor(kGreen);
	m4l_ZX[1][1]->SetFillColor(kBlue);
	m4l_ZX[2][1]->SetFillColor(kViolet);
	hs_4e->Add(m4l_ZX[2][1]);
	hs_4e->Add(m4l_ZX[1][1]);
	hs_4e->Add(m4l_ZX[0][1]);
	hs_4e->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4e.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4e.png");
	
	c1->cd();
	THStack *hs_2e2mu = new THStack("hs_2e2mu","test stacked histograms");
	m4l_ZX[0][2]->SetFillColor(kGreen);
	m4l_ZX[1][2]->SetFillColor(kBlue);
	m4l_ZX[2][2]->SetFillColor(kViolet);
	hs_2e2mu->Add(m4l_ZX[2][2]);
	hs_2e2mu->Add(m4l_ZX[1][2]);
	hs_2e2mu->Add(m4l_ZX[0][2]);
	hs_2e2mu->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2e2mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2e2mu.png");
	
	c1->cd();
	THStack *hs_2mu2e = new THStack("hs_2mu2e","test stacked histograms");
	m4l_ZX[0][3]->SetFillColor(kGreen);
	m4l_ZX[1][3]->SetFillColor(kBlue);
	m4l_ZX[2][3]->SetFillColor(kViolet);
	hs_2mu2e->Add(m4l_ZX[2][3]);
	hs_2mu2e->Add(m4l_ZX[1][3]);
	hs_2mu2e->Add(m4l_ZX[0][3]);
	hs_2mu2e->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2mu2e.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2mu2e.png");
	
	cout << "==========================================" << endl;
	cout << "[INFO] Total yields " << FRType << endl;
	cout << "==========================================" << endl;
	cout << "[INFO] DYJets " << endl;
	cout << "[INFO] 4mu   = " << m4l_ZX[0][0]->Integral() << endl;
	cout << "[INFO] 4e    = " << m4l_ZX[0][1]->Integral() << endl;
	cout << "[INFO] 2e2mu = " << m4l_ZX[0][2]->Integral() << endl;
	cout << "[INFO] 2mu2e = " << m4l_ZX[0][3]->Integral() << endl;
	cout << "=====================" << endl;
	cout << "[INFO] TTJets " << endl;
	cout << "[INFO] 4mu   = " << m4l_ZX[1][0]->Integral() << endl;
	cout << "[INFO] 4e    = " << m4l_ZX[1][1]->Integral() << endl;
	cout << "[INFO] 2e2mu = " << m4l_ZX[1][2]->Integral() << endl;
	cout << "[INFO] 2mu2e = " << m4l_ZX[1][3]->Integral() << endl;
	cout << "=====================" << endl;
	cout << "[INFO] WZTo3LNu " << endl;
	cout << "[INFO] 4mu   = " << m4l_ZX[2][0]->Integral() << endl;
	cout << "[INFO] 4e    = " << m4l_ZX[2][1]->Integral() << endl;
	cout << "[INFO] 2e2mu = " << m4l_ZX[2][2]->Integral() << endl;
	cout << "[INFO] 2mu2e = " << m4l_ZX[2][3]->Integral() << endl;
	cout << "=====================" << endl;
	cout << "[INFO] INCLUSIVE " << endl;
	cout << "[INFO] 4mu   = " << m4l_ZX[0][0]->Integral() + m4l_ZX[1][0]->Integral() + m4l_ZX[2][0]->Integral() << endl;
	cout << "[INFO] 4e    = " << m4l_ZX[0][1]->Integral() + m4l_ZX[1][1]->Integral() + m4l_ZX[2][1]->Integral() << endl;
	cout << "[INFO] 2e2mu = " << m4l_ZX[0][2]->Integral() + m4l_ZX[1][2]->Integral() + m4l_ZX[2][2]->Integral() << endl;
	cout << "[INFO] 2mu2e = " << m4l_ZX[0][3]->Integral() + m4l_ZX[1][3]->Integral() + m4l_ZX[2][3]->Integral() << endl;
	cout << "=====================" << endl;
	cout << "[INFO] 4l   = " <<  m4l_ZX[0][0]->Integral() + m4l_ZX[1][0]->Integral() + m4l_ZX[2][0]->Integral()
	+ m4l_ZX[0][1]->Integral() + m4l_ZX[1][1]->Integral() + m4l_ZX[2][1]->Integral()
	+ m4l_ZX[0][2]->Integral() + m4l_ZX[1][2]->Integral() + m4l_ZX[2][2]->Integral()
	+ m4l_ZX[0][3]->Integral() + m4l_ZX[1][3]->Integral() + m4l_ZX[2][3]->Integral() << endl;
}

void ResetHistos(TH1F* m4l_ZX[3][4])
{
	m4l_ZX[0][0]->Reset();
	m4l_ZX[0][1]->Reset();
	m4l_ZX[0][2]->Reset();
	m4l_ZX[0][3]->Reset();
	m4l_ZX[1][0]->Reset();
	m4l_ZX[1][1]->Reset();
	m4l_ZX[1][2]->Reset();
	m4l_ZX[1][3]->Reset();
	m4l_ZX[2][0]->Reset();
	m4l_ZX[2][1]->Reset();
	m4l_ZX[2][2]->Reset();
	m4l_ZX[2][3]->Reset();
}

void ApplyFakeRate()
{
	gROOT->SetBatch();
	
	float lumi = 35.9;
	
	TString path = "NewData/";
	TString file_name = "/ZZ4lAnalysis.root";
	
	TString DY       = path + "DYJetsToLL_M50" + file_name;
	TString TTJets   = path + "TTJets_DiLept_ext1" + file_name;
	TString WZTo3LNu = path + "WZTo3LNu" + file_name;
	
	TString histo_name;
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		histo_name ="h_LepPT_ele_EE_"+to_string(i_jf);
		h_LepPT[i_jf][0][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_ele_EB_"+to_string(i_jf);
		h_LepPT[i_jf][0][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_mu_EE_"+to_string(i_jf);
		h_LepPT[i_jf][1][0] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
		histo_name ="h_LepPT_mu_EB_"+to_string(i_jf);
		h_LepPT[i_jf][1][1] = new TH1D(histo_name,histo_name,_n_pT_bins, pT_bins);
	}
	
	fillDistributions(DY);
	fillDistributions(TTJets);
	fillDistributions(WZTo3LNu);
	
	TCanvas *c1 = new TCanvas("c1","c1",900,900);
	TString canvas_name;
	c1->cd();
	TH1D* sum_pT;
	
	for(int i_lep_flav = 0; i_lep_flav < 2; i_lep_flav++)
	{
		for(int i_EBorEE = 0; i_EBorEE < 2; i_EBorEE++)
		{
			sum_pT = (TH1D*)h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->Clone();
			sum_pT->Add(h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]);
			sum_pT->Add(h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]);
			sum_pT->Add(h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]);
			sum_pT->Add(h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]);
			
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetLineColor(kBlue);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kBlue);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetMaximum(1.0);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->SetMinimum(-0.01);
			h_LepPT[LeptonFlavours::LightJet][i_lep_flav][i_EBorEE]->Draw("HIST ");
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->SetLineColor(kBlack);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->SetMarkerColor(kBlack);
			h_LepPT[LeptonFlavours::NoMatch][i_lep_flav][i_EBorEE]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetLineColor(kRed);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->SetMarkerColor(kRed);
			h_LepPT[LeptonFlavours::HeavyJet][i_lep_flav][i_EBorEE]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->SetLineColor(kGreen);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->SetMarkerColor(kGreen);
			h_LepPT[LeptonFlavours::Conversion][i_lep_flav][i_EBorEE]->Draw("HIST SAME");
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->Divide(sum_pT);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->SetLineColor(kViolet);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->SetMarkerColor(kViolet);
			h_LepPT[LeptonFlavours::Lepton][i_lep_flav][i_EBorEE]->Draw("HIST SAME");
			
			canvas_name = "./MCTruthStudy/CRZLL_PT_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/CRZLL_PT_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".png";
			c1->SaveAs(canvas_name);
			
			sum_pT->Reset();
		}
	}
	
	TH1F* m4l_ZX[3][4];
	m4l_ZX[0][0] = new TH1F("m4l_ZX_DY_4mu","m4l_ZX_4mu",40,70.,870.);
	m4l_ZX[0][1] = new TH1F("m4l_ZX_DY_4e","m4l_ZX_4e",40,70.,870.);
	m4l_ZX[0][2] = new TH1F("m4l_ZX_DY_2e2mu","m4l_ZX_2e2mu",40,70.,870.);
	m4l_ZX[0][3] = new TH1F("m4l_ZX_DY_2mu2e","m4l_ZX_2mu2e",40,70.,870.);
	m4l_ZX[1][0] = new TH1F("m4l_ZX_TT_4mu","m4l_ZX_4mu",40,70.,870.);
	m4l_ZX[1][1] = new TH1F("m4l_ZX_TT_4e","m4l_ZX_4e",40,70.,870.);
	m4l_ZX[1][2] = new TH1F("m4l_ZX_TT_2e2mu","m4l_ZX_2e2mu",40,70.,870.);
	m4l_ZX[1][3] = new TH1F("m4l_ZX_TT_2mu2e","m4l_ZX_2mu2e",40,70.,870.);
	m4l_ZX[2][0] = new TH1F("m4l_ZX_WZ_4mu","m4l_ZX_4mu",40,70.,870.);
	m4l_ZX[2][1] = new TH1F("m4l_ZX_WZ_4e","m4l_ZX_4e",40,70.,870.);
	m4l_ZX[2][2] = new TH1F("m4l_ZX_WZ_2e2mu","m4l_ZX_2e2mu",40,70.,870.);
	m4l_ZX[2][3] = new TH1F("m4l_ZX_WZ_2mu2e","m4l_ZX_2mu2e",40,70.,870.);
	
	fillZXHisto_WithFlavourInfo(DY, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(TTJets, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "FlavourInfo");
	ResetHistos(m4l_ZX);
	
	fillZXHisto(DY, lumi, m4l_ZX);
	fillZXHisto(TTJets, lumi, m4l_ZX);
	fillZXHisto(WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "Average");
	ResetHistos(m4l_ZX);
	
	fillZXHisto_CutBased(DY, lumi, m4l_ZX);
	fillZXHisto_CutBased(TTJets, lumi, m4l_ZX);
	fillZXHisto_CutBased(WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "CutBased");
	ResetHistos(m4l_ZX);
}
