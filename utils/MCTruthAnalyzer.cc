// Standalone script to analyze composition of fakes based on MC truth matching
// Running with: root -l -b utils/MCTruthAnalyzer.cc++()

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
int _n_pT_bins = 8;

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

void RemoveNegativeBins2D(TH2F *h)
{
	for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
	{
		for (int i_bin_y = 1; i_bin_y <= h->GetYaxis()->GetNbins(); i_bin_y++)
		{
			if( h->GetBinContent(i_bin_x,i_bin_y) < 0.) h->SetBinContent(i_bin_x,i_bin_y,0);
		}
		
	}
	
}

TGraphErrors* makeFRgraph(int flav, TString eta, TH2F *passing, TH2F *failing)
{
	TGraphErrors *fr;
	vector<float> vector_X, vector_Y, vector_EX, vector_EY;
	int eta_bin = 0;
	if(eta == "EB") eta_bin = 1;
	if(eta == "EE") eta_bin = 2;
	
	RemoveNegativeBins2D(passing);
	RemoveNegativeBins2D(failing);
	
	for(int i_pT_bin = 0; i_pT_bin < _n_pT_bins - 1; i_pT_bin++ )
	{
		double temp_NP = 0;
		double temp_NF = 0;
		double temp_error_x = 0;
		double temp_error_NP = 0;
		double temp_error_NF = 0;
		
		if ( flav == 0 && i_pT_bin == 0) continue; // electrons do not have 5 - 7 GeV bin
		
		temp_NP = passing->IntegralAndError(passing->GetXaxis()->FindBin(pT_bins[i_pT_bin]),passing->GetXaxis()->FindBin(pT_bins[i_pT_bin+1]) - 1, eta_bin, eta_bin, temp_error_NP);
		temp_NF = failing->IntegralAndError(failing->GetXaxis()->FindBin(pT_bins[i_pT_bin]),failing->GetXaxis()->FindBin(pT_bins[i_pT_bin+1]) - 1, eta_bin, eta_bin, temp_error_NF);
		
		vector_X.push_back((pT_bins[i_pT_bin] + pT_bins[i_pT_bin + 1])/2);
		if ( temp_NP == 0) vector_Y.push_back(0.);
		else if ( temp_NF == 0) vector_Y.push_back(1.);
		else vector_Y.push_back(temp_NP/(temp_NP+temp_NF));
		
		vector_EX.push_back((pT_bins[i_pT_bin + 1] - pT_bins[i_pT_bin])/2);
		if ( temp_NP == 0 || temp_NF== 0) vector_EY.push_back(0.);
		else vector_EY.push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
	}
	
	fr = new TGraphErrors(vector_X.size(), &(vector_X[0]),&(vector_Y[0]),&(vector_EX[0]),&(vector_EY[0]));
	return fr;
	
}

void MCTruthAnalyzer()
{
	gROOT->ProcessLine(".L ./ext/setTDRStyle_cpp.so");
	gROOT->ProcessLine("setTDRStyle();");
	gROOT->SetBatch();
	
	TString path = "NewData/";
	TString file_name = "/ZZ4lAnalysis.root";
	
	TString DY      = path + "DYJetsToLL_M50" + file_name;
	
	TFile* DY_file = TFile::Open(DY);
	TTree* DY_Tree = (TTree*) DY_file->Get("CRZLTree/candTree");
	int  nentries = DY_Tree->GetEntries();
	TH1F* hCounters = (TH1F*)DY_file->Get("CRZLTree/Counters");
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
	
	DY_Tree->SetBranchAddress("PFMET",&PFMET);
	DY_Tree->SetBranchAddress("Z1Mass",&Z1Mass);
	DY_Tree->SetBranchAddress("Z1Flav",&Z1Flav);
	DY_Tree->SetBranchAddress("LepPt",&LepPt);
	DY_Tree->SetBranchAddress("LepEta",&LepEta);
	DY_Tree->SetBranchAddress("LepPhi",&LepPhi);
	DY_Tree->SetBranchAddress("LepLepId",&LepLepId);
	DY_Tree->SetBranchAddress("LepSIP",&LepSIP);
	DY_Tree->SetBranchAddress("LepBDT",&LepBDT);
	DY_Tree->SetBranchAddress("LepisID",&LepisID);
	DY_Tree->SetBranchAddress("LepCombRelIsoPF",&LepCombRelIsoPF);
	DY_Tree->SetBranchAddress("LepMissingHit", &LepMissingHit);
	DY_Tree->SetBranchAddress("passIsoPreFSR",&passIsoPreFSR);
	DY_Tree->SetBranchAddress("JetPt",&JetPt);
	DY_Tree->SetBranchAddress("JetEta",&JetEta);
	DY_Tree->SetBranchAddress("JetPhi",&JetPhi);
	DY_Tree->SetBranchAddress("JetMass",&JetMass);
	DY_Tree->SetBranchAddress("JetBTagger",&JetBTagger);
	DY_Tree->SetBranchAddress("JetIsBtagged",&JetIsBtagged);
	DY_Tree->SetBranchAddress("JetIsBtaggedWithSF",&JetIsBtaggedWithSF);
	DY_Tree->SetBranchAddress("JetQGLikelihood",&JetQGLikelihood);
	DY_Tree->SetBranchAddress("JetAxis2",&JetAxis2);
	DY_Tree->SetBranchAddress("JetMult",&JetMult);
	DY_Tree->SetBranchAddress("JetPtD",&JetPtD);
	DY_Tree->SetBranchAddress("JetSigma",&JetSigma);
	DY_Tree->SetBranchAddress("MatchJetPartonFlavour",&MatchJetPartonFlavour);
	DY_Tree->SetBranchAddress("MatchJetbTagger",&MatchJetbTagger);
	DY_Tree->SetBranchAddress("genHEPMCweight",&genHEPMCweight);
	DY_Tree->SetBranchAddress("PUWeight",&PUWeight);
	DY_Tree->SetBranchAddress("dataMCWeight",&dataMCWeight);
	DY_Tree->SetBranchAddress("trigEffWeight",&trigEffWeight);
	DY_Tree->SetBranchAddress("overallEventWeight",&overallEventWeight);
	DY_Tree->SetBranchAddress("HqTMCweight",&HqTMCweight);
	DY_Tree->SetBranchAddress("xsec",&xsec);
	DY_Tree->SetBranchAddress("GenMCTruthMatchId",&GenMCTruthMatchId);
	DY_Tree->SetBranchAddress("GenMCTruthMatchMotherId",&GenMCTruthMatchMotherId);
	
	int jet_category = -1;
	int lep_flavour = -1;
	int lep_EBorEE = -1;
	
	TH1D* h_LepPT[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
	TH1D* h_JetbTagger[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
	TH1D* h_LepSIP[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
	TH1D* h_LepBDT[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
	TH1D* h_LepMissingHit[LeptonFlavours::NUM_OF_FLAVOURS][2][2];
	TString histo_name;
	
	TH2F* passingSEL[LeptonFlavours::NUM_OF_FLAVOURS][2];
	TH2F* faillingSEL[LeptonFlavours::NUM_OF_FLAVOURS][2];
	
	int n_events[LeptonFlavours::NUM_OF_FLAVOURS];
	int n_events_afterCuts[LeptonFlavours::NUM_OF_FLAVOURS];
	
	TFile* fOutHistos = new TFile("output.root", "recreate");
	fOutHistos->cd();
	

	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		histo_name ="h_LepPT_ele_EE_"+to_string(i_jf);
		h_LepPT[i_jf][0][0] = new TH1D(histo_name,histo_name,8, 5.,85.);
		histo_name ="h_LepPT_ele_EB_"+to_string(i_jf);
		h_LepPT[i_jf][0][1] = new TH1D(histo_name,histo_name,8, 5.,85.);
		histo_name ="h_LepPT_mu_EE_"+to_string(i_jf);
		h_LepPT[i_jf][1][0] = new TH1D(histo_name,histo_name,8, 5.,85.);
		histo_name ="h_LepPT_mu_EB_"+to_string(i_jf);
		h_LepPT[i_jf][1][1] = new TH1D(histo_name,histo_name,8, 5.,85.);
		
		histo_name ="h_LepbTag_ele_EE_"+to_string(i_jf);
		h_JetbTagger[i_jf][0][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_ele_EB_"+to_string(i_jf);
		h_JetbTagger[i_jf][0][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_mu_EE_"+to_string(i_jf);
		h_JetbTagger[i_jf][1][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepbTag_mu_EB_"+to_string(i_jf);
		h_JetbTagger[i_jf][1][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		
		histo_name ="h_LepSIP_ele_EE_"+to_string(i_jf);
		h_LepSIP[i_jf][0][0] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_ele_EB_"+to_string(i_jf);
		h_LepSIP[i_jf][0][1] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_mu_EE_"+to_string(i_jf);
		h_LepSIP[i_jf][1][0] = new TH1D(histo_name,histo_name,40., 0.,4.);
		histo_name ="h_LepSIP_mu_EB_"+to_string(i_jf);
		h_LepSIP[i_jf][1][1] = new TH1D(histo_name,histo_name,40., 0.,4.);
		
		histo_name ="h_LepBDT_ele_EE_"+to_string(i_jf);
		h_LepBDT[i_jf][0][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_ele_EB_"+to_string(i_jf);
		h_LepBDT[i_jf][0][1] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_mu_EE_"+to_string(i_jf);
		h_LepBDT[i_jf][1][0] = new TH1D(histo_name,histo_name,20, -1., 1.);
		histo_name ="h_LepBDT_mu_EB_"+to_string(i_jf);
		h_LepBDT[i_jf][1][1] = new TH1D(histo_name,histo_name,20, -1., 1.);

		
		histo_name = "Passing_ele_" + to_string(i_jf);
		passingSEL[i_jf][0] = new TH2F(histo_name,"", 80, 0, 80, 2, 0, 2);
		histo_name = "Failing_ele_" + to_string(i_jf);
		faillingSEL[i_jf][0] = new TH2F(histo_name,"", 80, 0, 80, 2, 0, 2);
		
		histo_name = "Passing_mu_" + to_string(i_jf);
		passingSEL[i_jf][1] = new TH2F(histo_name,"", 80, 0, 80, 2, 0, 2);
		histo_name = "Failing_mu_" + to_string(i_jf);
		faillingSEL[i_jf][1] = new TH2F(histo_name,"", 80, 0, 80, 2, 0, 2);
		
		n_events[i_jf] = 0;
		n_events_afterCuts[i_jf] = 0;
	}

	for(int i_event = 0; i_event < nentries; i_event++)
	{
		DY_Tree->GetEvent(i_event);
		
		//Determine lepton origin
		jet_category = matchFlavour(MatchJetPartonFlavour->at(2), GenMCTruthMatchId->at(2), GenMCTruthMatchMotherId->at(2));
		lep_flavour = (fabs(LepLepId->at(2)) == 11) ? 0 : 1;
		if(lep_flavour == 0 && (abs(LepEta->at(2)) < 1.479) )  lep_EBorEE = 0;
		if(lep_flavour == 0 && (abs(LepEta->at(2)) >= 1.479) ) lep_EBorEE = 1;
		if(lep_flavour == 1 && (abs(LepEta->at(2)) < 1.2) )    lep_EBorEE = 0;
		if(lep_flavour == 1 && (abs(LepEta->at(2)) >= 1.2) )   lep_EBorEE = 1;
		

		n_events[jet_category]++;
		
		// Calculate fake rates
		if ( Z1Mass < 40. ) {continue;}
		if ( Z1Mass > 120. ) {continue;}
		if ( (LepPt->at(0) > LepPt->at(1)) && (LepPt->at(0) < 20. || LepPt->at(1) < 10.) ) {continue;}
		if ( (LepPt->at(1) > LepPt->at(0)) && (LepPt->at(1) < 20. || LepPt->at(0) < 10.) ) {continue;}
		if ( LepSIP->at(2) > 4.) {continue;}
		if ( PFMET > 25. ) {continue;}
		if ( LepisID->at(2) && LepCombRelIsoPF->at(2) < 0.35 )
		{
			if (fabs(LepLepId->at(2)) == 11) passingSEL[jet_category][0]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5, overallEventWeight);
			if (fabs(LepLepId->at(2)) == 13) passingSEL[jet_category][1]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5, overallEventWeight);
		}
		else
		{
			if (fabs(LepLepId->at(2)) == 11) faillingSEL[jet_category][0]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5, overallEventWeight);
			if (fabs(LepLepId->at(2)) == 13) faillingSEL[jet_category][1]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5, overallEventWeight);
		}
		
		// Fill distributions
		h_LepPT[jet_category][lep_flavour][lep_EBorEE]->Fill(LepPt->at(2), overallEventWeight);
		h_JetbTagger[jet_category][lep_flavour][lep_EBorEE]->Fill(MatchJetbTagger->at(2), overallEventWeight);
		h_LepBDT[jet_category][lep_flavour][lep_EBorEE]->Fill(LepBDT->at(2), overallEventWeight);
		h_LepSIP[jet_category][lep_flavour][lep_EBorEE]->Fill(LepSIP->at(2), overallEventWeight);
		
		n_events_afterCuts[jet_category]++;
	}
	
	cout << "==============================" << endl;
	cout << "[INFO] Control printout " << endl;
	cout << "==============================" << endl;
	cout << "[INFO] Total number of events BEFORE cuts = " << n_events[LeptonFlavours::Conversion] + n_events[LeptonFlavours::NoMatch] + n_events[LeptonFlavours::LightJet] + n_events[LeptonFlavours::HeavyJet] + n_events[LeptonFlavours::Lepton] << endl;
	cout << "[INFO] Number of events matched to leptons = " << n_events[LeptonFlavours::Lepton] << endl;
	cout << "[INFO] Number of events matched to photons = " << n_events[LeptonFlavours::Conversion] << endl;
	cout << "[INFO] Number of events not matched to jet = " << n_events[LeptonFlavours::NoMatch] << endl;
	cout << "[INFO] Number of events matched to light jet = " << n_events[LeptonFlavours::LightJet] << endl;
	cout << "[INFO] Number of events matched to heavy jet = " << n_events[LeptonFlavours::HeavyJet] << endl;
	cout << "============================================================" << endl;
	cout << "[INFO] Total number of events AFTER cuts = " << n_events_afterCuts[LeptonFlavours::Conversion] + n_events_afterCuts[LeptonFlavours::NoMatch] + n_events_afterCuts[LeptonFlavours::LightJet] + n_events_afterCuts[LeptonFlavours::HeavyJet] + n_events_afterCuts[LeptonFlavours::Lepton] << endl;
	cout << "[INFO] Number of events matched to leptons = " << n_events_afterCuts[LeptonFlavours::Lepton] << endl;
	cout << "[INFO] Number of events matched to photons = " << n_events_afterCuts[LeptonFlavours::Conversion] << endl;
	cout << "[INFO] Number of events not matched to jet = " << n_events_afterCuts[LeptonFlavours::NoMatch] << endl;
	cout << "[INFO] Number of events matched to light jet = " << n_events_afterCuts[LeptonFlavours::LightJet] << endl;
	cout << "[INFO] Number of events matched to heavy jet = " << n_events_afterCuts[LeptonFlavours::HeavyJet] << endl;
	cout << "============================================================" << endl;
	
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

			canvas_name = "./MCTruthStudy/PT_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/PT_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".png";
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
			
			canvas_name = "./MCTruthStudy/SIP_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/SIP_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".png";
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

			canvas_name = "./MCTruthStudy/bTaggScore_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".pdf";
			c1->SaveAs(canvas_name);
			canvas_name = "./MCTruthStudy/bTaggScore_distribution_" + to_string(i_lep_flav) + "_" + to_string(i_EBorEE) + ".png";
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
		
		canvas_name = "./MCTruthStudy/BDT_distribution_ele_" + to_string(i_EBorEE) + ".pdf";
		c1->SaveAs(canvas_name);
		canvas_name = "./MCTruthStudy/BDT_distribution_ele_" + to_string(i_EBorEE) + ".png";
		c1->SaveAs(canvas_name);
	}

	
	TGraphErrors *fr_ele_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_ele_EE[LeptonFlavours::NUM_OF_FLAVOURS];
	TGraphErrors *fr_mu_EB[LeptonFlavours::NUM_OF_FLAVOURS],*fr_mu_EE[LeptonFlavours::NUM_OF_FLAVOURS];
	TString fr_name;
	
	for (int i_jf = 0; i_jf < LeptonFlavours::NUM_OF_FLAVOURS; i_jf++)
	{
		fr_name = "FR_ele_EB_"+to_string(i_jf);
		fr_ele_EB[i_jf] = makeFRgraph(0, "EB" , passingSEL[i_jf][0], faillingSEL[i_jf][0]);
		fr_ele_EB[i_jf]->SetName(fr_name);
		fr_ele_EB[i_jf]->Write();
		
		fr_name = "FR_ele_EE_"+to_string(i_jf);
		fr_ele_EE[i_jf] = makeFRgraph(0, "EE", passingSEL[i_jf][0], faillingSEL[i_jf][0]);
		fr_ele_EE[i_jf]->SetName(fr_name);
		fr_ele_EE[i_jf]->Write();
		
		fr_name = "FR_mu_EB_"+to_string(i_jf);
		fr_mu_EB[i_jf] = makeFRgraph(1, "EB" , passingSEL[i_jf][1], faillingSEL[i_jf][1]);
		fr_mu_EB[i_jf]->SetName(fr_name);
		fr_mu_EB[i_jf]->Write();
		
		fr_name = "FR_mu_EE_"+to_string(i_jf);
		fr_mu_EE[i_jf] = makeFRgraph(1, "EE", passingSEL[i_jf][1], faillingSEL[i_jf][1]);
		fr_mu_EE[i_jf]->SetName(fr_name);
		fr_mu_EE[i_jf]->Write();
		
	}
	
//	auto c0 = new TCanvas("c0","multigraph L3",900,900);
//	fr_ele_EB[LeptonFlavours::LightJet]->SetLineColor(kBlue);
//	fr_ele_EB[LeptonFlavours::LightJet]->SetMarkerColor(kBlue);
//	fr_ele_EB[LeptonFlavours::LightJet]->Draw("EP0");
//	fr_ele_EB[LeptonFlavours::NoMatch]->SetLineColor(kBlack);
//	fr_ele_EB[LeptonFlavours::NoMatch]->SetMarkerColor(kBlack);
//	fr_ele_EB[LeptonFlavours::NoMatch]->Draw("EP0 SAME");
//	fr_ele_EB[LeptonFlavours::HeavyJet]->SetLineColor(kRed);
//	fr_ele_EB[LeptonFlavours::HeavyJet]->SetMarkerColor(kRed);
//	fr_ele_EB[LeptonFlavours::HeavyJet]->Draw("EP0 SAME");
//	fr_ele_EB[LeptonFlavours::Conversion]->SetLineColor(kGreen);
//	fr_ele_EB[LeptonFlavours::Conversion]->SetMarkerColor(kGreen);
//	fr_ele_EB[LeptonFlavours::Conversion]->Draw("EP0 SAME");
//	fr_ele_EB[LeptonFlavours::Lepton]->SetLineColor(kViolet);
//	fr_ele_EB[LeptonFlavours::Lepton]->SetMarkerColor(kViolet);
//	fr_ele_EB[LeptonFlavours::Lepton]->Draw("EP0 SAME");
//	fr_ele_EB[LeptonFlavours::LightJet]->SetMaximum(0.95);
//	c0->SaveAs("./MCTruthStudy/FakeRates_ele_EB.pdf");
//	c0->SaveAs("./MCTruthStudy/FakeRates_ele_EB.png");
//
//	fr_ele_EE[LeptonFlavours::LightJet]->SetLineColor(kBlue);
//	fr_ele_EE[LeptonFlavours::LightJet]->SetMarkerColor(kBlue);
//	fr_ele_EE[LeptonFlavours::LightJet]->Draw("EP0");
//	fr_ele_EE[LeptonFlavours::NoMatch]->SetLineColor(kBlack);
//	fr_ele_EE[LeptonFlavours::NoMatch]->SetMarkerColor(kBlack);
//	fr_ele_EE[LeptonFlavours::NoMatch]->Draw("EP0 SAME");
//	fr_ele_EE[LeptonFlavours::HeavyJet]->SetLineColor(kRed);
//	fr_ele_EE[LeptonFlavours::HeavyJet]->SetMarkerColor(kRed);
//	fr_ele_EE[LeptonFlavours::HeavyJet]->Draw("EP0 SAME");
//	fr_ele_EE[LeptonFlavours::Conversion]->SetLineColor(kGreen);
//	fr_ele_EE[LeptonFlavours::Conversion]->SetMarkerColor(kGreen);
//	fr_ele_EE[LeptonFlavours::Conversion]->Draw("EP0 SAME");
//	fr_ele_EE[LeptonFlavours::Lepton]->SetLineColor(kViolet);
//	fr_ele_EE[LeptonFlavours::Lepton]->SetMarkerColor(kViolet);
//	fr_ele_EE[LeptonFlavours::Lepton]->Draw("EP0 SAME");
//	fr_ele_EE[LeptonFlavours::LightJet]->SetMaximum(0.95);
//	c0->SaveAs("./MCTruthStudy/FakeRates_ele_EE.pdf");
//	c0->SaveAs("./MCTruthStudy/FakeRates_ele_EE.png");
//
//	fr_mu_EB[LeptonFlavours::LightJet]->SetLineColor(kBlue);
//	fr_mu_EB[LeptonFlavours::LightJet]->SetMarkerColor(kBlue);
//	fr_mu_EB[LeptonFlavours::LightJet]->Draw("EP0");
//	fr_mu_EB[LeptonFlavours::NoMatch]->SetLineColor(kBlack);
//	fr_mu_EB[LeptonFlavours::NoMatch]->SetMarkerColor(kBlack);
//	fr_mu_EB[LeptonFlavours::NoMatch]->Draw("EP0 SAME");
//	fr_mu_EB[LeptonFlavours::HeavyJet]->SetLineColor(kRed);
//	fr_mu_EB[LeptonFlavours::HeavyJet]->SetMarkerColor(kRed);
//	fr_mu_EB[LeptonFlavours::HeavyJet]->Draw("EP0 SAME");
//	fr_mu_EB[LeptonFlavours::Conversion]->SetLineColor(kGreen);
//	fr_mu_EB[LeptonFlavours::Conversion]->SetMarkerColor(kGreen);
//	fr_mu_EB[LeptonFlavours::Conversion]->Draw("EP0 SAME");
//	fr_mu_EB[LeptonFlavours::Lepton]->SetLineColor(kViolet);
//	fr_mu_EB[LeptonFlavours::Lepton]->SetMarkerColor(kViolet);
//	fr_mu_EB[LeptonFlavours::Lepton]->Draw("EP0 SAME");
//	fr_mu_EB[LeptonFlavours::LightJet]->SetMaximum(0.95);
//	c0->SaveAs("./MCTruthStudy/FakeRates_mu_EB.pdf");
//	c0->SaveAs("./MCTruthStudy/FakeRates_mu_EB.png");
//
//	fr_mu_EE[LeptonFlavours::LightJet]->SetLineColor(kBlue);
//	fr_mu_EE[LeptonFlavours::LightJet]->SetMarkerColor(kBlue);
//	fr_mu_EE[LeptonFlavours::LightJet]->Draw("EP0");
//	fr_mu_EE[LeptonFlavours::NoMatch]->SetLineColor(kBlack);
//	fr_mu_EE[LeptonFlavours::NoMatch]->SetMarkerColor(kBlack);
//	fr_mu_EE[LeptonFlavours::NoMatch]->Draw("EP0 SAME");
//	fr_mu_EE[LeptonFlavours::HeavyJet]->SetLineColor(kRed);
//	fr_mu_EE[LeptonFlavours::HeavyJet]->SetMarkerColor(kRed);
//	fr_mu_EE[LeptonFlavours::HeavyJet]->Draw("EP0 SAME");
//	fr_mu_EE[LeptonFlavours::Conversion]->SetLineColor(kGreen);
//	fr_mu_EE[LeptonFlavours::Conversion]->SetMarkerColor(kGreen);
//	fr_mu_EE[LeptonFlavours::Conversion]->Draw("EP0 SAME");
//	fr_mu_EE[LeptonFlavours::Lepton]->SetLineColor(kViolet);
//	fr_mu_EE[LeptonFlavours::Lepton]->SetMarkerColor(kViolet);
//	fr_mu_EE[LeptonFlavours::Lepton]->Draw("EP0 SAME");
//	fr_mu_EE[LeptonFlavours::LightJet]->SetMaximum(0.95);
//	c0->SaveAs("./MCTruthStudy/FakeRates_mu_EE.pdf");
//	c0->SaveAs("./MCTruthStudy/FakeRates_mu_EE.png");

	fOutHistos->Close();
	delete fOutHistos;
	
}
