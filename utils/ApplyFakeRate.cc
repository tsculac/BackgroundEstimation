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
int n_events[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
vector <TString> _s_CR, _s_LepFlav, _s_EEorEB, _s_MatchFlavour;


int matchFlavour(int MatchJetPartonFlavour, float dR_lep_RecoJet, int GenMCTruthMatchId, int GenMCTruthMatchMotherId)
{
	int flavour = 0;
	
	
	//no match to jets but match to photon
	if ( (MatchJetPartonFlavour == 0 || fabs(dR_lep_RecoJet) > 0.3) && (fabs(GenMCTruthMatchId) == 22 ||  fabs(GenMCTruthMatchMotherId) == 22) ) flavour = LeptonFlavours::Conversion;
	
	//no match to jets but match to lepton
	else if ( (MatchJetPartonFlavour == 0 || fabs(dR_lep_RecoJet) > 0.3) && (fabs(GenMCTruthMatchId) == 11 ||  fabs(GenMCTruthMatchId) == 13) ) flavour = LeptonFlavours::Lepton;
	
	//no match to jets
	else if ( MatchJetPartonFlavour == 0 || fabs(dR_lep_RecoJet) > 0.3) flavour = LeptonFlavours::NoMatch;
	
	//match to light (u,d,s & g) jets
	else if ( fabs(dR_lep_RecoJet)<=0.3 && (fabs(MatchJetPartonFlavour) == 1 || fabs(MatchJetPartonFlavour) == 2 || fabs(MatchJetPartonFlavour) == 3 || fabs(MatchJetPartonFlavour) == 21) ) flavour = LeptonFlavours::LightJet;
	
	//match to heavy jets (c & b)
	else if ( fabs(dR_lep_RecoJet)<=0.3 && (fabs(MatchJetPartonFlavour) == 4 || fabs(MatchJetPartonFlavour) == 5 ) ) flavour = LeptonFlavours::HeavyJet;
	
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
	
	TFile *fr_file = TFile::Open("output_DY.root");

	TString fr_name;
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		fr_name = "FR_ele_EB_"+_s_MatchFlavour.at(i_jf);
		fr_ele_EB[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_ele_EE_"+_s_MatchFlavour.at(i_jf);
		fr_ele_EE[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EB_"+_s_MatchFlavour.at(i_jf);
		fr_mu_EB[i_jf] = (TGraphErrors*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EE_"+_s_MatchFlavour.at(i_jf);
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
	vector<float>   *dR_lep_RecoJet = 0;
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
	_Tree->SetBranchAddress("dR_lep_RecoJet",&dR_lep_RecoJet);
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
		jet1_category = matchFlavour(MatchJetPartonFlavour->at(2), dR_lep_RecoJet->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		jet2_category = matchFlavour(MatchJetPartonFlavour->at(3), dR_lep_RecoJet->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3));
		
		_event_weight = (lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;
		if (_current_proces == 0) _event_weight*=0.33;
		
		// Calculate yield
		_yield_SR = _fs_ROS_SS.at(_current_final_state)*GetFakeRate_WithFlavour(LepPt->at(2),LepEta->at(2),LepLepId->at(2),jet1_category)*GetFakeRate_WithFlavour(LepPt->at(3),LepEta->at(3),LepLepId->at(3),jet2_category);
		m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _yield_SR*_event_weight);
		
//		if (jet1_category == LeptonFlavours::Conversion || jet2_category == LeptonFlavours::Conversion)
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


void fillZXHisto(bool corrected, TString input_file_name, float lumi, TH1F* m4l_ZX[3][4])
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
	
	TFile *fr_file = TFile::Open("FakeRates_SS_Moriond17_MC.root");
	
	if (corrected)
	{
		g_FR_mu_EB = (TGraphErrors*)fr_file->Get("FR_SS_muon_EB");
		g_FR_mu_EE = (TGraphErrors*)fr_file->Get("FR_SS_muon_EE");
		g_FR_e_EB  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EB");
		g_FR_e_EE  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EE");
	}
	else
	{
		g_FR_mu_EB = (TGraphErrors*)fr_file->Get("FR_SS_muon_EB_unc");
		g_FR_mu_EE = (TGraphErrors*)fr_file->Get("FR_SS_muon_EE_unc");
		g_FR_e_EB  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EB_unc");
		g_FR_e_EE  = (TGraphErrors*)fr_file->Get("FR_SS_electron_EE_unc");
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
	vector<float>   *dR_lep_RecoJet = 0;
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
	_Tree->SetBranchAddress("dR_lep_RecoJet",&dR_lep_RecoJet);
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
		if (_current_proces == 0) _event_weight*=0.33;
		
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
	vector<float>   *dR_lep_RecoJet = 0;
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
	_Tree->SetBranchAddress("dR_lep_RecoJet",&dR_lep_RecoJet);
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
		if (_current_proces == 0) _event_weight*=0.33;
		
		//if (matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2)) == LeptonFlavours::Conversion || matchFlavour(MatchJetPartonFlavour->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3)) == LeptonFlavours::Conversion) continue;
		// Calculate yield
		if(LepisID->at(2) && LepisID->at(3) && LepCombRelIsoPF->at(2) < 0.35 && LepCombRelIsoPF->at(3) < 0.35) m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _event_weight);
		
		
	}
	
	cout << "[INFO] Z+X histograms filled for file " + input_file_name + "!" << endl;
}



void DrawHistos(TH1F* m4l_ZX[3][4], TString FRType)
{
	TCanvas *c1 = new TCanvas("c1"+FRType,"c1"+FRType,900,900);
	THStack *hs_4mu = new THStack("hs_4mu","test stacked histograms");
	c1->cd();
	m4l_ZX[0][0]->SetFillColor(kGreen-1);
	m4l_ZX[1][0]->SetFillColor(kBlue);
	m4l_ZX[2][0]->SetFillColor(kViolet);
	hs_4mu->Add(m4l_ZX[2][0]);
	hs_4mu->Add(m4l_ZX[1][0]);
	hs_4mu->Add(m4l_ZX[0][0]);
	hs_4mu->Draw("HIST");
	hs_4mu->GetXaxis()->SetTitle("m_{4#mu} [GeV]");
	hs_4mu->GetYaxis()->SetTitle("Events");
	if(FRType != "CutBased") hs_4mu->SetMinimum(0.);
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4mu.png");
	
	c1->cd();
	THStack *hs_4e = new THStack("hs_4e","test stacked histograms");
	m4l_ZX[0][1]->SetFillColor(kGreen-1);
	m4l_ZX[1][1]->SetFillColor(kBlue);
	m4l_ZX[2][1]->SetFillColor(kViolet);
	hs_4e->Add(m4l_ZX[2][1]);
	hs_4e->Add(m4l_ZX[1][1]);
	hs_4e->Add(m4l_ZX[0][1]);
	hs_4e->Draw("HIST");
	hs_4e->GetXaxis()->SetTitle("m_{4e} [GeV]");
	hs_4e->GetYaxis()->SetTitle("Events");
	if(FRType != "CutBased") hs_4e->SetMinimum(0.);
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4e.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_4e.png");
	
	c1->cd();
	THStack *hs_2e2mu = new THStack("hs_2e2mu","test stacked histograms");
	m4l_ZX[0][2]->SetFillColor(kGreen-1);
	m4l_ZX[1][2]->SetFillColor(kBlue);
	m4l_ZX[2][2]->SetFillColor(kViolet);
	hs_2e2mu->Add(m4l_ZX[2][2]);
	hs_2e2mu->Add(m4l_ZX[1][2]);
	hs_2e2mu->Add(m4l_ZX[0][2]);
	hs_2e2mu->Draw("HIST");
	hs_2e2mu->GetXaxis()->SetTitle("m_{2e2#mu} [GeV]");
	hs_2e2mu->GetYaxis()->SetTitle("Events");
	if(FRType != "CutBased") hs_2e2mu->SetMinimum(0.);
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2e2mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_"+FRType+"_FR_2e2mu.png");
	
	c1->cd();
	THStack *hs_2mu2e = new THStack("hs_2mu2e","test stacked histograms");
	m4l_ZX[0][3]->SetFillColor(kGreen-1);
	m4l_ZX[1][3]->SetFillColor(kBlue);
	m4l_ZX[2][3]->SetFillColor(kViolet);
	hs_2mu2e->Add(m4l_ZX[2][3]);
	hs_2mu2e->Add(m4l_ZX[1][3]);
	hs_2mu2e->Add(m4l_ZX[0][3]);
	hs_2mu2e->Draw("HIST");
	hs_2mu2e->GetXaxis()->SetTitle("m_{2#mu2e} [GeV]");
	hs_2mu2e->GetYaxis()->SetTitle("Events");
	if(FRType != "CutBased") hs_2mu2e->SetMinimum(0.);
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
	gStyle->SetOptStat(0);
	
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
	
	TString path = "NewMC/";
	TString path_ext = "ImprovedMatching_MC/";
	TString file_name = "/ZZ4lAnalysis.root";
	
	TString DY       = path + "DYJetsToLL_M50" + file_name;
	TString TTJets   = path + "TTJets_DiLept_ext1" + file_name;
	TString WZTo3LNu = path + "WZTo3LNu" + file_name;
	
	TString DY_mix = path_ext + "DYJetsToLL_M50" + file_name;
	TString DY_0j  = path_ext + "DYToLL_0J" + file_name;
	TString DY_0jb = path_ext + "DYToLL_0J_backup" + file_name;
	TString DY_1j  = path_ext + "DYToLL_1J" + file_name;
	TString DY_1jb = path_ext + "DYToLL_1J_backup" + file_name;
	TString DY_2j  = path_ext + "DYToLL_2J" + file_name;
	TString DY_2jb = path_ext + "DYToLL_2J_backup" + file_name;
	
	TH1F* m4l_ZX[3][4];
	int n_bins = 20;
	m4l_ZX[0][0] = new TH1F("m4l_ZX_DY_4mu","m4l_ZX_4mu",n_bins,70.,870.);
	m4l_ZX[0][1] = new TH1F("m4l_ZX_DY_4e","m4l_ZX_4e",n_bins,70.,870.);
	m4l_ZX[0][2] = new TH1F("m4l_ZX_DY_2e2mu","m4l_ZX_2e2mu",n_bins,70.,870.);
	m4l_ZX[0][3] = new TH1F("m4l_ZX_DY_2mu2e","m4l_ZX_2mu2e",n_bins,70.,870.);
	m4l_ZX[1][0] = new TH1F("m4l_ZX_TT_4mu","m4l_ZX_4mu",n_bins,70.,870.);
	m4l_ZX[1][1] = new TH1F("m4l_ZX_TT_4e","m4l_ZX_4e",n_bins,70.,870.);
	m4l_ZX[1][2] = new TH1F("m4l_ZX_TT_2e2mu","m4l_ZX_2e2mu",n_bins,70.,870.);
	m4l_ZX[1][3] = new TH1F("m4l_ZX_TT_2mu2e","m4l_ZX_2mu2e",n_bins,70.,870.);
	m4l_ZX[2][0] = new TH1F("m4l_ZX_WZ_4mu","m4l_ZX_4mu",n_bins,70.,870.);
	m4l_ZX[2][1] = new TH1F("m4l_ZX_WZ_4e","m4l_ZX_4e",n_bins,70.,870.);
	m4l_ZX[2][2] = new TH1F("m4l_ZX_WZ_2e2mu","m4l_ZX_2e2mu",n_bins,70.,870.);
	m4l_ZX[2][3] = new TH1F("m4l_ZX_WZ_2mu2e","m4l_ZX_2mu2e",n_bins,70.,870.);
	
	TH1F *m4l_4e_cutBased,*m4l_4e_Flav,*m4l_4e_avgUnc,*m4l_4e_avg;
	TH1F *m4l_4mu_cutBased,*m4l_4mu_Flav,*m4l_4mu_avgUnc,*m4l_4mu_avg;
	
	m4l_4e_cutBased = new TH1F("m4l_ZX_DY_4e_cutBased","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4e_Flav     = new TH1F("m4l_ZX_DY_4e_Flav","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4e_avgUnc   = new TH1F("m4l_ZX_DY_4e_avgUnc","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4e_avg      = new TH1F("m4l_ZX_DY_4e_avg","m4l_ZX_4e",n_bins,70.,870.);
	
	m4l_4mu_cutBased = new TH1F("m4l_ZX_DY_4mu_cutBased","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4mu_Flav     = new TH1F("m4l_ZX_DY_4mu_Flav","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4mu_avgUnc   = new TH1F("m4l_ZX_DY_4mu_avgUnc","m4l_ZX_4e",n_bins,70.,870.);
	m4l_4mu_avg      = new TH1F("m4l_ZX_DY_4mu_avg","m4l_ZX_4e",n_bins,70.,870.);
	
	
	fillZXHisto_WithFlavourInfo(DY_mix, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_0j, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_0jb, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_1j, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_1jb, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_2j, lumi, m4l_ZX);
	fillZXHisto_WithFlavourInfo(DY_2jb, lumi, m4l_ZX);
	//fillZXHisto_WithFlavourInfo(TTJets, lumi, m4l_ZX);
	//fillZXHisto_WithFlavourInfo(WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "FlavourInfo");
	m4l_4e_Flav = (TH1F*)m4l_ZX[0][1]->Clone();
	m4l_4e_Flav->Add(m4l_ZX[1][1]);
	m4l_4mu_Flav = (TH1F*)m4l_ZX[0][0]->Clone();
	m4l_4mu_Flav->Add(m4l_ZX[1][0]);
	ResetHistos(m4l_ZX);
	
	fillZXHisto(true, DY_mix, lumi, m4l_ZX);
	fillZXHisto(true, DY_0j, lumi, m4l_ZX);
	fillZXHisto(true, DY_0jb, lumi, m4l_ZX);
	fillZXHisto(true, DY_1j, lumi, m4l_ZX);
	fillZXHisto(true, DY_1jb, lumi, m4l_ZX);
	fillZXHisto(true, DY_2j, lumi, m4l_ZX);
	fillZXHisto(true, DY_2jb, lumi, m4l_ZX);
	//fillZXHisto(true, TTJets, lumi, m4l_ZX);
	//fillZXHisto(true, WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "Average_ConversionCorrected");
	m4l_4e_avg = (TH1F*)m4l_ZX[0][1]->Clone();
	m4l_4e_avg->Add(m4l_ZX[1][1]);
	m4l_4mu_avg = (TH1F*)m4l_ZX[0][0]->Clone();
	m4l_4mu_avg->Add(m4l_ZX[1][0]);
	ResetHistos(m4l_ZX);
	
	fillZXHisto(false, DY_mix, lumi, m4l_ZX);
	fillZXHisto(false, DY_0j, lumi, m4l_ZX);
	fillZXHisto(false, DY_0jb, lumi, m4l_ZX);
	fillZXHisto(false, DY_1j, lumi, m4l_ZX);
	fillZXHisto(false, DY_1jb, lumi, m4l_ZX);
	fillZXHisto(false, DY_2j, lumi, m4l_ZX);
	fillZXHisto(false, DY_2jb, lumi, m4l_ZX);
	//fillZXHisto(false, TTJets, lumi, m4l_ZX);
	//fillZXHisto(false, WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "Average_UnCorrected");
	m4l_4e_avgUnc = (TH1F*)m4l_ZX[0][1]->Clone();
	m4l_4e_avgUnc->Add(m4l_ZX[1][1]);
	m4l_4mu_avgUnc = (TH1F*)m4l_ZX[0][0]->Clone();
	m4l_4mu_avgUnc->Add(m4l_ZX[1][0]);
	ResetHistos(m4l_ZX);
	
	fillZXHisto_CutBased(DY_mix, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_0j, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_0jb, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_1j, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_1jb, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_2j, lumi, m4l_ZX);
	fillZXHisto_CutBased(DY_2jb, lumi, m4l_ZX);
	//fillZXHisto_CutBased(TTJets, lumi, m4l_ZX);
	//fillZXHisto_CutBased(WZTo3LNu, lumi, m4l_ZX);
	DrawHistos(m4l_ZX, "CutBased");
	m4l_4e_cutBased = (TH1F*)m4l_ZX[0][1]->Clone();
	m4l_4e_cutBased->Add(m4l_ZX[1][1]);
	m4l_4mu_cutBased = (TH1F*)m4l_ZX[0][0]->Clone();
	m4l_4mu_cutBased->Add(m4l_ZX[1][0]);
	ResetHistos(m4l_ZX);
	
	TCanvas *c = new TCanvas("c","c",900,900);
	c->cd();
	m4l_4e_avgUnc->SetLineColor(kBlack);
	m4l_4e_avgUnc->SetLineWidth(3);
	m4l_4e_avgUnc->SetLineStyle(2);
	m4l_4e_avgUnc->SetFillColorAlpha(kWhite, 0.0);
	m4l_4e_avgUnc->Draw("HIST");
	m4l_4e_cutBased->SetLineColor(kGreen-1);
	m4l_4e_cutBased->SetLineWidth(3);
	m4l_4e_cutBased->SetFillColorAlpha(kWhite, 0.0);
	m4l_4e_cutBased->Draw("HIST SAME");
	m4l_4e_avg->SetLineColor(kBlack);
	m4l_4e_avg->SetLineWidth(3);
	m4l_4e_avg->SetFillColorAlpha(kWhite, 0.0);
	m4l_4e_avg->Draw("HIST SAME");
	m4l_4e_Flav->SetLineColor(kBlue);
	m4l_4e_Flav->SetLineWidth(3);
	m4l_4e_Flav->SetFillColorAlpha(kWhite, 0.0);
	m4l_4e_Flav->Draw("HIST SAME");
	
	c->SaveAs("./MCTruthStudy/ZX_CompareMethods_4e.pdf");
	c->SaveAs("./MCTruthStudy/ZX_CompareMethods_4e.png");
	
	c->Clear();
	m4l_4mu_avgUnc->SetLineColor(kBlack);
	m4l_4mu_avgUnc->SetLineWidth(3);
	m4l_4mu_avgUnc->SetLineStyle(2);
	m4l_4mu_avgUnc->SetFillColorAlpha(kWhite, 0.0);
	m4l_4mu_avgUnc->Draw("HIST");
	m4l_4mu_cutBased->SetLineColor(kGreen-1);
	m4l_4mu_cutBased->SetLineWidth(3);
	m4l_4mu_cutBased->SetFillColorAlpha(kWhite, 0.0);
	m4l_4mu_cutBased->Draw("HIST SAME");
	m4l_4mu_avg->SetLineColor(kBlack);
	m4l_4mu_avg->SetLineWidth(3);
	m4l_4mu_avg->SetFillColorAlpha(kWhite, 0.0);
	m4l_4mu_avg->Draw("HIST SAME");
	m4l_4mu_Flav->SetLineColor(kBlue);
	m4l_4mu_Flav->SetLineWidth(3);
	m4l_4mu_Flav->SetFillColorAlpha(kWhite, 0.0);
	m4l_4mu_Flav->Draw("HIST SAME");
	
	c->SaveAs("./MCTruthStudy/ZX_CompareMethods_4mu.pdf");
	c->SaveAs("./MCTruthStudy/ZX_CompareMethods_4mu.png");
	
}
