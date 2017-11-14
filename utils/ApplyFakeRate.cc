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

TH1D *fr_ele_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_ele_EE[LeptonFlavours::NUM_OF_FLAVOURS];
TH1D *fr_mu_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_mu_EE[LeptonFlavours::NUM_OF_FLAVOURS];


int matchFlavourOnlyJets(int MatchJetPartonFlavour, int GenMCTruthMatchId, int GenMCTruthMatchMotherId)
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

float GetFakeRate(float LepPt, float LepEta, int LepId, int flavour)
{
	float fr = -1.;
	
	float my_lep_Pt = LepPt >= 55. ? 54. : LepPt;
	int   my_lep_ID = abs(LepId);
	
	int bin = 0;
	if ( my_lep_Pt > 5 && my_lep_Pt <= 15 ) bin = 1;
	else if ( my_lep_Pt > 15 && my_lep_Pt <= 25 ) bin = 2;
	else if ( my_lep_Pt > 25 && my_lep_Pt <= 35 ) bin = 3;
	else if ( my_lep_Pt > 35 && my_lep_Pt <= 45 ) bin = 4;
	else if ( my_lep_Pt > 45 && my_lep_Pt <= 55 ) bin = 5;
	else if ( my_lep_Pt > 55 && my_lep_Pt <= 65 ) bin = 6;
	else if ( my_lep_Pt > 65 && my_lep_Pt <= 75 ) bin = 7;
	else if ( my_lep_Pt > 75 && my_lep_Pt <= 85 ) bin = 8;

	if ( my_lep_ID == 11 )
	{
		if ( fabs(LepEta) < 1.479 )
			fr = fr_ele_EB[flavour]->GetBinContent(bin);
		else
			fr = fr_ele_EE[flavour]->GetBinContent(bin);
	}
	else if ( my_lep_ID == 13 )
	{
		if ( fabs(LepEta) < 1.2 )
			fr = fr_mu_EB[flavour]->GetBinContent(bin);
		else
			fr = fr_mu_EE[flavour]->GetBinContent(bin);
	}
	else
	{
		cout << "[ERROR] Wrong lepton ID: " << my_lep_ID << endl;
		fr = 0;
	}
	return fr;
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
		fr_ele_EB[i_jf] = (TH1D*)fr_file->Get(fr_name);
		
		fr_name = "FR_ele_EE_"+to_string(i_jf);
		fr_ele_EE[i_jf] = (TH1D*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EB_"+to_string(i_jf);
		fr_mu_EB[i_jf] = (TH1D*)fr_file->Get(fr_name);
		
		fr_name = "FR_mu_EE_"+to_string(i_jf);
		fr_mu_EE[i_jf] = (TH1D*)fr_file->Get(fr_name);

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
	int lep1_flavour = -1;
	int jet2_category = -1;
	int lep2_flavour = -1;
	int _current_final_state = -1;
	int _current_proces = -1;

	if(input_file_name.Contains("DY"))     _current_proces = 0;
	if(input_file_name.Contains("TTJets")) _current_proces = 1;
	if(input_file_name.Contains("WZ"))     _current_proces = 2;
	
	TFile* fOutHistos = new TFile("./MCTruthStudy/ZXyields.root", "recreate");
	fOutHistos->cd();

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
		jet1_category = matchFlavourOnlyJets(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		lep1_flavour = (fabs(LepLepId->at(2)) == 11) ? 0 : 1;
		
		jet2_category = matchFlavourOnlyJets(MatchJetPartonFlavour->at(3), GenMCTruthMatchId->at(3), GenMCTruthMatchMotherId->at(3));
		lep2_flavour = (fabs(LepLepId->at(3)) == 11) ? 0 : 1;

		
		_event_weight = (lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;
		
		// Calculate yield
		_yield_SR = _fs_ROS_SS.at(_current_final_state)*GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2),jet1_category)*GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3),jet2_category);
//		cout << "YIELD = " << _yield_SR << endl;
//		cout << "WEIGHT = " << _event_weight << endl;
//		cout << "==================" << endl;
		m4l_ZX[_current_proces][_current_final_state]->Fill(ZZMass, _yield_SR*_event_weight);

		
	}

cout << "[INFO] Z+X histograms filled for file " + input_file_name + "!" << endl;
}

void DrawHistos(TH1F* m4l_ZX[3][4])
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
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_4mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_4mu.png");
	
	c1->cd();
	THStack *hs_4e = new THStack("hs_4e","test stacked histograms");
	m4l_ZX[0][1]->SetFillColor(kGreen);
	m4l_ZX[1][1]->SetFillColor(kBlue);
	m4l_ZX[2][1]->SetFillColor(kViolet);
	hs_4e->Add(m4l_ZX[2][1]);
	hs_4e->Add(m4l_ZX[1][1]);
	hs_4e->Add(m4l_ZX[0][1]);
	hs_4e->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_4e.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_4e.png");
	
	c1->cd();
	THStack *hs_2e2mu = new THStack("hs_2e2mu","test stacked histograms");
	m4l_ZX[0][2]->SetFillColor(kGreen);
	m4l_ZX[1][2]->SetFillColor(kBlue);
	m4l_ZX[2][2]->SetFillColor(kViolet);
	hs_2e2mu->Add(m4l_ZX[2][2]);
	hs_2e2mu->Add(m4l_ZX[1][2]);
	hs_2e2mu->Add(m4l_ZX[0][2]);
	hs_2e2mu->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_2e2mu.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_2e2mu.png");
	
	c1->cd();
	THStack *hs_2mu2e = new THStack("hs_2mu2e","test stacked histograms");
	m4l_ZX[0][3]->SetFillColor(kGreen);
	m4l_ZX[1][3]->SetFillColor(kBlue);
	m4l_ZX[2][3]->SetFillColor(kViolet);
	hs_2mu2e->Add(m4l_ZX[2][3]);
	hs_2mu2e->Add(m4l_ZX[1][3]);
	hs_2mu2e->Add(m4l_ZX[0][3]);
	hs_2mu2e->Draw("HIST");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_2mu2e.pdf");
	c1->SaveAs("./MCTruthStudy/ZX_flavourInfoinFR_2mu2e.png");
	
	cout << "=====================" << endl;
	cout << "[INFO] Total yields" << endl;
	cout << "=====================" << endl;
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
	
	float lumi = 38.0;
	
	TString path = "NewData/";
	TString file_name = "/ZZ4lAnalysis.root";
	
	TString DY       = path + "DYJetsToLL_M50" + file_name;
	TString TTJets   = path + "TTJets_DiLept_ext1" + file_name;
	TString WZTo3LNu = path + "WZTo3LNu" + file_name;
	
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
	
	DrawHistos(m4l_ZX);
	ResetHistos(m4l_ZX);
	
	
}
