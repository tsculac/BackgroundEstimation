//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 21 14:52:28 2016 by ROOT version 6.02/05
// from TTree candTree/Event Summary
// found on file: ICHEP_2016/ggH125/ZZ4lAnalysis.root
//////////////////////////////////////////////////////////

#ifndef Tree_h
#define Tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class Tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNumber;
   Long64_t        EventNumber;
   Int_t           LumiNumber;
   Short_t         NRecoMu;
   Short_t         NRecoEle;
   Short_t         Nvtx;
   Short_t         NObsInt;
   Float_t         NTrueInt;
   Float_t         PFMET;
   Float_t         PFMETPhi;
   Float_t         PFMETNoHF;
   Float_t         PFMETNoHFPhi;
   Short_t         nCleanedJets;
   Short_t         nCleanedJetsPt30;
   Short_t         nCleanedJetsPt30_jecUp;
   Short_t         nCleanedJetsPt30_jecDn;
   Short_t         nCleanedJetsPt30BTagged;
   Short_t         nCleanedJetsPt30BTagged_bTagSF;
   Short_t         trigWord;
   Float_t         ZZMass;
   Float_t         ZZMassErr;
   Float_t         ZZMassErrCorr;
   Float_t         ZZMassPreFSR;
   Short_t         ZZsel;
   Float_t         ZZPt;
   Float_t         ZZEta;
   Float_t         ZZPhi;
   Int_t           CRflag;
   Float_t         Z1Mass;
   Float_t         Z1Pt;
   Short_t         Z1Flav;
   Float_t         ZZMassRefit;
   Float_t         ZZMassRefitErr;
   Float_t         ZZMassUnrefitErr;
   Float_t         Z2Mass;
   Float_t         Z2Pt;
   Short_t         Z2Flav;
   Float_t         costhetastar;
   Float_t         helphi;
   Float_t         helcosthetaZ1;
   Float_t         helcosthetaZ2;
   Float_t         phistarZ1;
   Float_t         phistarZ2;
   Float_t         xi;
   Float_t         xistar;
   vector<float>   *LepPt;
   vector<float>   *LepEta;
   vector<float>   *LepPhi;
   vector<short>   *LepLepId;
   vector<float>   *LepSIP;
   vector<float>   *LepTime;
   vector<bool>    *LepisID;
   vector<short>   *LepisLoose;
   vector<float>   *LepBDT;
   vector<char>    *LepMissingHit;
   vector<float>   *LepCombRelIsoPF;
   vector<float>   *fsrPt;
   vector<float>   *fsrEta;
   vector<float>   *fsrPhi;
   vector<short>   *fsrLept;
   Bool_t          passIsoPreFSR;
   Float_t         p0plus_VAJHU;
   Float_t         p0_g1prime2_VAJHU;
   Float_t         p0hplus_VAJHU;
   Float_t         p0minus_VAJHU;
   Float_t         p0_g1prime2_zgs_VAJHU;
   Float_t         p0hplus_zgs_VAJHU;
   Float_t         p0minus_zgs_VAJHU;
   Float_t         p0hplus_gsgs_VAJHU;
   Float_t         p0minus_gsgs_VAJHU;
   Float_t         pg1g1prime2_VAJHU;
   Float_t         pg1g2_VAJHU;
   Float_t         pg1g2_pi2_VAJHU;
   Float_t         pg1g4_VAJHU;
   Float_t         pg1g4_pi2_VAJHU;
   Float_t         p0plus_zz_g1prime2_zgs_VAJHU;
   Float_t         p0plus_zz_g1prime2_zgs_pi2_VAJHU;
   Float_t         p0plus_zz_0hplus_zgs_VAJHU;
   Float_t         p0plus_zz_0minus_zgs_VAJHU;
   Float_t         p0plus_zz_0hplus_gsgs_VAJHU;
   Float_t         p0plus_zz_0minus_gsgs_VAJHU;
   Float_t         p1_VAJHU;
   Float_t         p1_prodIndep_VAJHU;
   Float_t         p1plus_VAJHU;
   Float_t         p1plus_prodIndep_VAJHU;
   Float_t         p2plus_gg_VAJHU;
   Float_t         p2plus_prodIndep_VAJHU;
   Float_t         p2plus_qqb_VAJHU;
   Float_t         p2h2plus_gg_VAJHU;
   Float_t         p2h2plus_qqb_VAJHU;
   Float_t         p2h2plus_prodIndep_VAJHU;
   Float_t         p2h3plus_gg_VAJHU;
   Float_t         p2h3plus_qqb_VAJHU;
   Float_t         p2h3plus_prodIndep_VAJHU;
   Float_t         p2h4plus_gg_VAJHU;
   Float_t         p2h4plus_qqb_VAJHU;
   Float_t         p2h4plus_prodIndep_VAJHU;
   Float_t         p2bplus_gg_VAJHU;
   Float_t         p2bplus_qqb_VAJHU;
   Float_t         p2bplus_prodIndep_VAJHU;
   Float_t         p2h6plus_gg_VAJHU;
   Float_t         p2h6plus_qqb_VAJHU;
   Float_t         p2h6plus_prodIndep_VAJHU;
   Float_t         p2h7plus_gg_VAJHU;
   Float_t         p2h7plus_qqb_VAJHU;
   Float_t         p2h7plus_prodIndep_VAJHU;
   Float_t         p2hminus_gg_VAJHU;
   Float_t         p2hminus_qqb_VAJHU;
   Float_t         p2hminus_prodIndep_VAJHU;
   Float_t         p2h9minus_gg_VAJHU;
   Float_t         p2h9minus_qqb_VAJHU;
   Float_t         p2h9minus_prodIndep_VAJHU;
   Float_t         p2h10minus_gg_VAJHU;
   Float_t         p2h10minus_qqb_VAJHU;
   Float_t         p2h10minus_prodIndep_VAJHU;
   Float_t         p0plus_VAMCFM;
   Float_t         ggzz_VAMCFM;
   Float_t         ggzz_p0plus_VAMCFM;
   Float_t         bkg_VAMCFM;
   Float_t         bkg_prodIndep_VAMCFM;
   Float_t         pZJJ_VAMCFM;
   Float_t         Dgg10_VAMCFM;
   Float_t         p0plus_m4l;
   Float_t         p0plus_m4l_ScaleUp;
   Float_t         p0plus_m4l_ScaleDown;
   Float_t         p0plus_m4l_ResUp;
   Float_t         p0plus_m4l_ResDown;
   Float_t         bkg_m4l;
   Float_t         bkg_m4l_ScaleUp;
   Float_t         bkg_m4l_ScaleDown;
   Float_t         bkg_m4l_ResUp;
   Float_t         bkg_m4l_ResDown;
   Float_t         pwh_leptonic_VAJHU;
   Float_t         pzh_leptonic_VAJHU;
   Float_t         phjj_VAJHU_highestPTJets;
   Float_t         pvbf_VAJHU_highestPTJets;
   Float_t         phjj_VAJHU_highestPTJets_up;
   Float_t         pvbf_VAJHU_highestPTJets_up;
   Float_t         phjj_VAJHU_highestPTJets_dn;
   Float_t         pvbf_VAJHU_highestPTJets_dn;
   Float_t         phjj_VAJHU_bestDjet;
   Float_t         pvbf_VAJHU_bestDjet;
   Float_t         phjj_VAJHU_bestDjet_up;
   Float_t         pvbf_VAJHU_bestDjet_up;
   Float_t         phjj_VAJHU_bestDjet_dn;
   Float_t         pvbf_VAJHU_bestDjet_dn;
   Float_t         pAux_vbf_VAJHU;
   Float_t         pAux_vbf_VAJHU_up;
   Float_t         pAux_vbf_VAJHU_dn;
   Float_t         phj_VAJHU;
   Float_t         phj_VAJHU_up;
   Float_t         phj_VAJHU_dn;
   Float_t         pwh_hadronic_VAJHU;
   Float_t         pwh_hadronic_VAJHU_up;
   Float_t         pwh_hadronic_VAJHU_dn;
   Float_t         pzh_hadronic_VAJHU;
   Float_t         pzh_hadronic_VAJHU_up;
   Float_t         pzh_hadronic_VAJHU_dn;
   Float_t         ptth_VAJHU;
   Float_t         ptth_VAJHU_up;
   Float_t         ptth_VAJHU_dn;
   Float_t         pbbh_VAJHU;
   Float_t         pbbh_VAJHU_up;
   Float_t         pbbh_VAJHU_dn;
   vector<float>   *JetPt;
   vector<float>   *JetEta;
   vector<float>   *JetPhi;
   vector<float>   *JetMass;
   vector<float>   *JetBTagger;
   vector<float>   *JetIsBtagged;
   vector<float>   *JetQGLikelihood;
   vector<float>   *JetAxis2;
   vector<float>   *JetMult;
   vector<float>   *JetPtD;
   vector<float>   *JetSigma;
   Float_t         DiJetMass;
   Float_t         DiJetDEta;
   Float_t         DiJetFisher;
   Short_t         nExtraLep;
   Short_t         nExtraZ;
   vector<float>   *ExtraLepPt;
   vector<float>   *ExtraLepEta;
   vector<float>   *ExtraLepPhi;
   vector<short>   *ExtraLepLepId;
   Float_t         ZXFakeweight;
   Float_t         KFactor_QCD_ggZZ_Nominal;
   Float_t         KFactor_QCD_ggZZ_PDFScaleDn;
   Float_t         KFactor_QCD_ggZZ_PDFScaleUp;
   Float_t         KFactor_QCD_ggZZ_QCDScaleDn;
   Float_t         KFactor_QCD_ggZZ_QCDScaleUp;
   Float_t         KFactor_QCD_ggZZ_AsDn;
   Float_t         KFactor_QCD_ggZZ_AsUp;
   Float_t         KFactor_QCD_ggZZ_PDFReplicaDn;
   Float_t         KFactor_QCD_ggZZ_PDFReplicaUp;
   Float_t         KFactor_EW_qqZZ;
   Float_t         KFactor_EW_qqZZ_unc;
   Float_t         KFactor_QCD_qqZZ_dPhi;
   Float_t         KFactor_QCD_qqZZ_M;
   Float_t         KFactor_QCD_qqZZ_Pt;
   Short_t         genFinalState;
   Int_t           genProcessId;
   Float_t         genHEPMCweight;
   Float_t         PUWeight;
   Float_t         dataMCWeight;
   Float_t         trigEffWeight;
   Float_t         overallEventWeight;
   Float_t         HqTMCweight;
   Float_t         xsec;
   Short_t         genExtInfo;
   Float_t         GenHMass;
   Float_t         GenHPt;
   Float_t         GenHRapidity;
   Float_t         GenZ1Mass;
   Float_t         GenZ1Pt;
   Float_t         GenZ1Phi;
   Float_t         GenZ1Flav;
   Float_t         GenZ2Mass;
   Float_t         GenZ2Pt;
   Float_t         GenZ2Phi;
   Float_t         GenZ2Flav;
   Float_t         GenLep1Pt;
   Float_t         GenLep1Eta;
   Float_t         GenLep1Phi;
   Short_t         GenLep1Id;
   Float_t         GenLep2Pt;
   Float_t         GenLep2Eta;
   Float_t         GenLep2Phi;
   Short_t         GenLep2Id;
   Float_t         GenLep3Pt;
   Float_t         GenLep3Eta;
   Float_t         GenLep3Phi;
   Short_t         GenLep3Id;
   Float_t         GenLep4Pt;
   Float_t         GenLep4Eta;
   Float_t         GenLep4Phi;
   Short_t         GenLep4Id;
   Float_t         GenAssocLep1Pt;
   Float_t         GenAssocLep1Eta;
   Float_t         GenAssocLep1Phi;
   Short_t         GenAssocLep1Id;
   Float_t         GenAssocLep2Pt;
   Float_t         GenAssocLep2Eta;
   Float_t         GenAssocLep2Phi;
   Short_t         GenAssocLep2Id;
   vector<float>   *reweightingweights;
   Float_t         LHEPDFScale;
   Float_t         LHEweight_QCDscale_muR1_muF1;
   Float_t         LHEweight_QCDscale_muR1_muF2;
   Float_t         LHEweight_QCDscale_muR1_muF0p5;
   Float_t         LHEweight_QCDscale_muR2_muF1;
   Float_t         LHEweight_QCDscale_muR2_muF2;
   Float_t         LHEweight_QCDscale_muR2_muF0p5;
   Float_t         LHEweight_QCDscale_muR0p5_muF1;
   Float_t         LHEweight_QCDscale_muR0p5_muF2;
   Float_t         LHEweight_QCDscale_muR0p5_muF0p5;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_NRecoMu;   //!
   TBranch        *b_NRecoEle;   //!
   TBranch        *b_Nvtx;   //!
   TBranch        *b_NObsInt;   //!
   TBranch        *b_NTrueInt;   //!
   TBranch        *b_PFMET;   //!
   TBranch        *b_PFMETPhi;   //!
   TBranch        *b_PFMETNoHF;   //!
   TBranch        *b_PFMETNoHFPhi;   //!
   TBranch        *b_nCleanedJets;   //!
   TBranch        *b_nCleanedJetsPt30;   //!
   TBranch        *b_nCleanedJetsPt30_jecUp;   //!
   TBranch        *b_nCleanedJetsPt30_jecDn;   //!
   TBranch        *b_nCleanedJetsPt30BTagged;   //!
   TBranch        *b_nCleanedJetsPt30BTagged_bTagSF;   //!
   TBranch        *b_trigWord;   //!
   TBranch        *b_ZZMass;   //!
   TBranch        *b_ZZMassErr;   //!
   TBranch        *b_ZZMassErrCorr;   //!
   TBranch        *b_ZZMassPreFSR;   //!
   TBranch        *b_ZZsel;   //!
   TBranch        *b_ZZPt;   //!
   TBranch        *b_ZZEta;   //!
   TBranch        *b_ZZPhi;   //!
   TBranch        *b_CRflag;   //!
   TBranch        *b_Z1Mass;   //!
   TBranch        *b_Z1Pt;   //!
   TBranch        *b_Z1Flav;   //!
   TBranch        *b_ZZMassRefit;   //!
   TBranch        *b_ZZMassRefitErr;   //!
   TBranch        *b_ZZMassUnrefitErr;   //!
   TBranch        *b_Z2Mass;   //!
   TBranch        *b_Z2Pt;   //!
   TBranch        *b_Z2Flav;   //!
   TBranch        *b_costhetastar;   //!
   TBranch        *b_helphi;   //!
   TBranch        *b_helcosthetaZ1;   //!
   TBranch        *b_helcosthetaZ2;   //!
   TBranch        *b_phistarZ1;   //!
   TBranch        *b_phistarZ2;   //!
   TBranch        *b_xi;   //!
   TBranch        *b_xistar;   //!
   TBranch        *b_LepPt;   //!
   TBranch        *b_LepEta;   //!
   TBranch        *b_LepPhi;   //!
   TBranch        *b_LepLepId;   //!
   TBranch        *b_LepSIP;   //!
   TBranch        *b_LepTime;   //!
   TBranch        *b_LepisID;   //!
   TBranch        *b_LepisLoose;   //!
   TBranch        *b_LepBDT;   //!
   TBranch        *b_LepMissingHit;   //!
   TBranch        *b_LepCombRelIsoPF;   //!
   TBranch        *b_fsrPt;   //!
   TBranch        *b_fsrEta;   //!
   TBranch        *b_fsrPhi;   //!
   TBranch        *b_fsrLept;   //!
   TBranch        *b_passIsoPreFSR;   //!
   TBranch        *b_p0plus_VAJHU;   //!
   TBranch        *b_p0_g1prime2_VAJHU;   //!
   TBranch        *b_p0hplus_VAJHU;   //!
   TBranch        *b_p0minus_VAJHU;   //!
   TBranch        *b_p0_g1prime2_zgs_VAJHU;   //!
   TBranch        *b_p0hplus_zgs_VAJHU;   //!
   TBranch        *b_p0minus_zgs_VAJHU;   //!
   TBranch        *b_p0hplus_gsgs_VAJHU;   //!
   TBranch        *b_p0minus_gsgs_VAJHU;   //!
   TBranch        *b_pg1g1prime2_VAJHU;   //!
   TBranch        *b_pg1g2_VAJHU;   //!
   TBranch        *b_pg1g2_pi2_VAJHU;   //!
   TBranch        *b_pg1g4_VAJHU;   //!
   TBranch        *b_pg1g4_pi2_VAJHU;   //!
   TBranch        *b_p0plus_zz_g1prime2_zgs_VAJHU;   //!
   TBranch        *b_p0plus_zz_g1prime2_zgs_pi2_VAJHU;   //!
   TBranch        *b_p0plus_zz_0hplus_zgs_VAJHU;   //!
   TBranch        *b_p0plus_zz_0minus_zgs_VAJHU;   //!
   TBranch        *b_p0plus_zz_0hplus_gsgs_VAJHU;   //!
   TBranch        *b_p0plus_zz_0minus_gsgs_VAJHU;   //!
   TBranch        *b_p1_VAJHU;   //!
   TBranch        *b_p1_prodIndep_VAJHU;   //!
   TBranch        *b_p1plus_VAJHU;   //!
   TBranch        *b_p1plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2plus_gg_VAJHU;   //!
   TBranch        *b_p2plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2plus_qqb_VAJHU;   //!
   TBranch        *b_p2h2plus_gg_VAJHU;   //!
   TBranch        *b_p2h2plus_qqb_VAJHU;   //!
   TBranch        *b_p2h2plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h3plus_gg_VAJHU;   //!
   TBranch        *b_p2h3plus_qqb_VAJHU;   //!
   TBranch        *b_p2h3plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h4plus_gg_VAJHU;   //!
   TBranch        *b_p2h4plus_qqb_VAJHU;   //!
   TBranch        *b_p2h4plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2bplus_gg_VAJHU;   //!
   TBranch        *b_p2bplus_qqb_VAJHU;   //!
   TBranch        *b_p2bplus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h6plus_gg_VAJHU;   //!
   TBranch        *b_p2h6plus_qqb_VAJHU;   //!
   TBranch        *b_p2h6plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h7plus_gg_VAJHU;   //!
   TBranch        *b_p2h7plus_qqb_VAJHU;   //!
   TBranch        *b_p2h7plus_prodIndep_VAJHU;   //!
   TBranch        *b_p2hminus_gg_VAJHU;   //!
   TBranch        *b_p2hminus_qqb_VAJHU;   //!
   TBranch        *b_p2hminus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h9minus_gg_VAJHU;   //!
   TBranch        *b_p2h9minus_qqb_VAJHU;   //!
   TBranch        *b_p2h9minus_prodIndep_VAJHU;   //!
   TBranch        *b_p2h10minus_gg_VAJHU;   //!
   TBranch        *b_p2h10minus_qqb_VAJHU;   //!
   TBranch        *b_p2h10minus_prodIndep_VAJHU;   //!
   TBranch        *b_p0plus_VAMCFM;   //!
   TBranch        *b_ggzz_VAMCFM;   //!
   TBranch        *b_ggzz_p0plus_VAMCFM;   //!
   TBranch        *b_bkg_VAMCFM;   //!
   TBranch        *b_bkg_prodIndep_VAMCFM;   //!
   TBranch        *b_pZJJ_VAMCFM;   //!
   TBranch        *b_Dgg10_VAMCFM;   //!
   TBranch        *b_p0plus_m4l;   //!
   TBranch        *b_p0plus_m4l_ScaleUp;   //!
   TBranch        *b_p0plus_m4l_ScaleDown;   //!
   TBranch        *b_p0plus_m4l_ResUp;   //!
   TBranch        *b_p0plus_m4l_ResDown;   //!
   TBranch        *b_bkg_m4l;   //!
   TBranch        *b_bkg_m4l_ScaleUp;   //!
   TBranch        *b_bkg_m4l_ScaleDown;   //!
   TBranch        *b_bkg_m4l_ResUp;   //!
   TBranch        *b_bkg_m4l_ResDown;   //!
   TBranch        *b_pwh_leptonic_VAJHU;   //!
   TBranch        *b_pzh_leptonic_VAJHU;   //!
   TBranch        *b_phjj_VAJHU_highestPTJets;   //!
   TBranch        *b_pvbf_VAJHU_highestPTJets;   //!
   TBranch        *b_phjj_VAJHU_highestPTJets_up;   //!
   TBranch        *b_pvbf_VAJHU_highestPTJets_up;   //!
   TBranch        *b_phjj_VAJHU_highestPTJets_dn;   //!
   TBranch        *b_pvbf_VAJHU_highestPTJets_dn;   //!
   TBranch        *b_phjj_VAJHU_bestDjet;   //!
   TBranch        *b_pvbf_VAJHU_bestDjet;   //!
   TBranch        *b_phjj_VAJHU_bestDjet_up;   //!
   TBranch        *b_pvbf_VAJHU_bestDjet_up;   //!
   TBranch        *b_phjj_VAJHU_bestDjet_dn;   //!
   TBranch        *b_pvbf_VAJHU_bestDjet_dn;   //!
   TBranch        *b_pAux_vbf_VAJHU;   //!
   TBranch        *b_pAux_vbf_VAJHU_up;   //!
   TBranch        *b_pAux_vbf_VAJHU_dn;   //!
   TBranch        *b_phj_VAJHU;   //!
   TBranch        *b_phj_VAJHU_up;   //!
   TBranch        *b_phj_VAJHU_dn;   //!
   TBranch        *b_pwh_hadronic_VAJHU;   //!
   TBranch        *b_pwh_hadronic_VAJHU_up;   //!
   TBranch        *b_pwh_hadronic_VAJHU_dn;   //!
   TBranch        *b_pzh_hadronic_VAJHU;   //!
   TBranch        *b_pzh_hadronic_VAJHU_up;   //!
   TBranch        *b_pzh_hadronic_VAJHU_dn;   //!
   TBranch        *b_ptth_VAJHU;   //!
   TBranch        *b_ptth_VAJHU_up;   //!
   TBranch        *b_ptth_VAJHU_dn;   //!
   TBranch        *b_pbbh_VAJHU;   //!
   TBranch        *b_pbbh_VAJHU_up;   //!
   TBranch        *b_pbbh_VAJHU_dn;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetMass;   //!
   TBranch        *b_JetBTagger;   //!
   TBranch        *b_JetIsBtagged;   //!
   TBranch        *b_JetQGLikelihood;   //!
   TBranch        *b_JetAxis2;   //!
   TBranch        *b_JetMult;   //!
   TBranch        *b_JetPtD;   //!
   TBranch        *b_JetSigma;   //!
   TBranch        *b_DiJetMass;   //!
   TBranch        *b_DiJetDEta;   //!
   TBranch        *b_DiJetFisher;   //!
   TBranch        *b_nExtraLep;   //!
   TBranch        *b_nExtraZ;   //!
   TBranch        *b_ExtraLepPt;   //!
   TBranch        *b_ExtraLepEta;   //!
   TBranch        *b_ExtraLepPhi;   //!
   TBranch        *b_ExtraLepLepId;   //!
   TBranch        *b_ZXFakeweight;   //!
   TBranch        *b_KFactor_QCD_ggZZ_Nominal;   //!
   TBranch        *b_KFactor_QCD_ggZZ_PDFScaleDn;   //!
   TBranch        *b_KFactor_QCD_ggZZ_PDFScaleUp;   //!
   TBranch        *b_KFactor_QCD_ggZZ_QCDScaleDn;   //!
   TBranch        *b_KFactor_QCD_ggZZ_QCDScaleUp;   //!
   TBranch        *b_KFactor_QCD_ggZZ_AsDn;   //!
   TBranch        *b_KFactor_QCD_ggZZ_AsUp;   //!
   TBranch        *b_KFactor_QCD_ggZZ_PDFReplicaDn;   //!
   TBranch        *b_KFactor_QCD_ggZZ_PDFReplicaUp;   //!
   TBranch        *b_KFactor_EW_qqZZ;   //!
   TBranch        *b_KFactor_EW_qqZZ_unc;   //!
   TBranch        *b_KFactor_QCD_qqZZ_dPhi;   //!
   TBranch        *b_KFactor_QCD_qqZZ_M;   //!
   TBranch        *b_KFactor_QCD_qqZZ_Pt;   //!
   TBranch        *b_genFinalState;   //!
   TBranch        *b_genProcessId;   //!
   TBranch        *b_genHEPMCweight;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_dataMCWeight;   //!
   TBranch        *b_trigEffWeight;   //!
   TBranch        *b_overallEventWeight;   //!
   TBranch        *b_HqTMCweight;   //!
   TBranch        *b_xsec;   //!
   TBranch        *b_genExtInfo;   //!
   TBranch        *b_GenHMass;   //!
   TBranch        *b_GenHPt;   //!
   TBranch        *b_GenHRapidity;   //!
   TBranch        *b_GenZ1Mass;   //!
   TBranch        *b_GenZ1Pt;   //!
   TBranch        *b_GenZ1Phi;   //!
   TBranch        *b_GenZ1Flav;   //!
   TBranch        *b_GenZ2Mass;   //!
   TBranch        *b_GenZ2Pt;   //!
   TBranch        *b_GenZ2Phi;   //!
   TBranch        *b_GenZ2Flav;   //!
   TBranch        *b_GenLep1Pt;   //!
   TBranch        *b_GenLep1Eta;   //!
   TBranch        *b_GenLep1Phi;   //!
   TBranch        *b_GenLep1Id;   //!
   TBranch        *b_GenLep2Pt;   //!
   TBranch        *b_GenLep2Eta;   //!
   TBranch        *b_GenLep2Phi;   //!
   TBranch        *b_GenLep2Id;   //!
   TBranch        *b_GenLep3Pt;   //!
   TBranch        *b_GenLep3Eta;   //!
   TBranch        *b_GenLep3Phi;   //!
   TBranch        *b_GenLep3Id;   //!
   TBranch        *b_GenLep4Pt;   //!
   TBranch        *b_GenLep4Eta;   //!
   TBranch        *b_GenLep4Phi;   //!
   TBranch        *b_GenLep4Id;   //!
   TBranch        *b_GenAssocLep1Pt;   //!
   TBranch        *b_GenAssocLep1Eta;   //!
   TBranch        *b_GenAssocLep1Phi;   //!
   TBranch        *b_GenAssocLep1Id;   //!
   TBranch        *b_GenAssocLep2Pt;   //!
   TBranch        *b_GenAssocLep2Eta;   //!
   TBranch        *b_GenAssocLep2Phi;   //!
   TBranch        *b_GenAssocLep2Id;   //!
   TBranch        *b_reweightingweights;   //!
   TBranch        *b_LHEPDFScale;   //!
   TBranch        *b_LHEweight_QCDscale_muR1_muF1;   //!
   TBranch        *b_LHEweight_QCDscale_muR1_muF2;   //!
   TBranch        *b_LHEweight_QCDscale_muR1_muF0p5;   //!
   TBranch        *b_LHEweight_QCDscale_muR2_muF1;   //!
   TBranch        *b_LHEweight_QCDscale_muR2_muF2;   //!
   TBranch        *b_LHEweight_QCDscale_muR2_muF0p5;   //!
   TBranch        *b_LHEweight_QCDscale_muR0p5_muF1;   //!
   TBranch        *b_LHEweight_QCDscale_muR0p5_muF2;   //!
   TBranch        *b_LHEweight_QCDscale_muR0p5_muF0p5;   //!

   Tree(TTree *tree=0);
   virtual ~Tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TString input_file_name);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Tree_cxx
Tree::Tree(TTree *tree) : fChain(0) 
{
}

Tree::~Tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Tree::Init(TTree *tree, TString input_file_name)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LepPt = 0;
   LepEta = 0;
   LepPhi = 0;
   LepLepId = 0;
   LepSIP = 0;
   LepTime = 0;
   LepisID = 0;
   LepisLoose = 0;
   LepBDT = 0;
   LepMissingHit = 0;
   LepCombRelIsoPF = 0;
   fsrPt = 0;
   fsrEta = 0;
   fsrPhi = 0;
   fsrLept = 0;
   JetPt = 0;
   JetEta = 0;
   JetPhi = 0;
   JetMass = 0;
   JetBTagger = 0;
   JetIsBtagged = 0;
   JetQGLikelihood = 0;
   JetAxis2 = 0;
   JetMult = 0;
   JetPtD = 0;
   JetSigma = 0;
   ExtraLepPt = 0;
   ExtraLepEta = 0;
   ExtraLepPhi = 0;
   ExtraLepLepId = 0;
   reweightingweights = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("NRecoMu", &NRecoMu, &b_NRecoMu);
   fChain->SetBranchAddress("NRecoEle", &NRecoEle, &b_NRecoEle);
   fChain->SetBranchAddress("Nvtx", &Nvtx, &b_Nvtx);
   fChain->SetBranchAddress("NObsInt", &NObsInt, &b_NObsInt);
   fChain->SetBranchAddress("NTrueInt", &NTrueInt, &b_NTrueInt);
   fChain->SetBranchAddress("PFMET", &PFMET, &b_PFMET);
   fChain->SetBranchAddress("PFMETPhi", &PFMETPhi, &b_PFMETPhi);
   fChain->SetBranchAddress("PFMETNoHF", &PFMETNoHF, &b_PFMETNoHF);
   fChain->SetBranchAddress("PFMETNoHFPhi", &PFMETNoHFPhi, &b_PFMETNoHFPhi);
   fChain->SetBranchAddress("nCleanedJets", &nCleanedJets, &b_nCleanedJets);
   fChain->SetBranchAddress("nCleanedJetsPt30", &nCleanedJetsPt30, &b_nCleanedJetsPt30);
   fChain->SetBranchAddress("nCleanedJetsPt30_jecUp", &nCleanedJetsPt30_jecUp, &b_nCleanedJetsPt30_jecUp);
   fChain->SetBranchAddress("nCleanedJetsPt30_jecDn", &nCleanedJetsPt30_jecDn, &b_nCleanedJetsPt30_jecDn);
   fChain->SetBranchAddress("nCleanedJetsPt30BTagged", &nCleanedJetsPt30BTagged, &b_nCleanedJetsPt30BTagged);
   fChain->SetBranchAddress("nCleanedJetsPt30BTagged_bTagSF", &nCleanedJetsPt30BTagged_bTagSF, &b_nCleanedJetsPt30BTagged_bTagSF);
   fChain->SetBranchAddress("trigWord", &trigWord, &b_trigWord);
   fChain->SetBranchAddress("ZZMass", &ZZMass, &b_ZZMass);
   fChain->SetBranchAddress("ZZMassErr", &ZZMassErr, &b_ZZMassErr);
   fChain->SetBranchAddress("ZZMassErrCorr", &ZZMassErrCorr, &b_ZZMassErrCorr);
   fChain->SetBranchAddress("ZZMassPreFSR", &ZZMassPreFSR, &b_ZZMassPreFSR);
   fChain->SetBranchAddress("ZZsel", &ZZsel, &b_ZZsel);
   fChain->SetBranchAddress("ZZPt", &ZZPt, &b_ZZPt);
   fChain->SetBranchAddress("ZZEta", &ZZEta, &b_ZZEta);
   fChain->SetBranchAddress("ZZPhi", &ZZPhi, &b_ZZPhi);
   fChain->SetBranchAddress("CRflag", &CRflag, &b_CRflag);
   fChain->SetBranchAddress("Z1Mass", &Z1Mass, &b_Z1Mass);
   fChain->SetBranchAddress("Z1Pt", &Z1Pt, &b_Z1Pt);
   fChain->SetBranchAddress("Z1Flav", &Z1Flav, &b_Z1Flav);
   fChain->SetBranchAddress("ZZMassRefit", &ZZMassRefit, &b_ZZMassRefit);
   fChain->SetBranchAddress("ZZMassRefitErr", &ZZMassRefitErr, &b_ZZMassRefitErr);
   fChain->SetBranchAddress("ZZMassUnrefitErr", &ZZMassUnrefitErr, &b_ZZMassUnrefitErr);
   fChain->SetBranchAddress("Z2Mass", &Z2Mass, &b_Z2Mass);
   fChain->SetBranchAddress("Z2Pt", &Z2Pt, &b_Z2Pt);
   fChain->SetBranchAddress("Z2Flav", &Z2Flav, &b_Z2Flav);
   fChain->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   fChain->SetBranchAddress("helphi", &helphi, &b_helphi);
   fChain->SetBranchAddress("helcosthetaZ1", &helcosthetaZ1, &b_helcosthetaZ1);
   fChain->SetBranchAddress("helcosthetaZ2", &helcosthetaZ2, &b_helcosthetaZ2);
   fChain->SetBranchAddress("phistarZ1", &phistarZ1, &b_phistarZ1);
   fChain->SetBranchAddress("phistarZ2", &phistarZ2, &b_phistarZ2);
   fChain->SetBranchAddress("xi", &xi, &b_xi);
   fChain->SetBranchAddress("xistar", &xistar, &b_xistar);
   fChain->SetBranchAddress("LepPt", &LepPt, &b_LepPt);
   fChain->SetBranchAddress("LepEta", &LepEta, &b_LepEta);
   fChain->SetBranchAddress("LepPhi", &LepPhi, &b_LepPhi);
   fChain->SetBranchAddress("LepLepId", &LepLepId, &b_LepLepId);
   fChain->SetBranchAddress("LepSIP", &LepSIP, &b_LepSIP);
   fChain->SetBranchAddress("LepTime", &LepTime, &b_LepTime);
   fChain->SetBranchAddress("LepisID", &LepisID, &b_LepisID);
   fChain->SetBranchAddress("LepisLoose", &LepisLoose, &b_LepisLoose);
   fChain->SetBranchAddress("LepBDT", &LepBDT, &b_LepBDT);
   fChain->SetBranchAddress("LepMissingHit", &LepMissingHit, &b_LepMissingHit);
   fChain->SetBranchAddress("LepCombRelIsoPF", &LepCombRelIsoPF, &b_LepCombRelIsoPF);
   fChain->SetBranchAddress("fsrPt", &fsrPt, &b_fsrPt);
   fChain->SetBranchAddress("fsrEta", &fsrEta, &b_fsrEta);
   fChain->SetBranchAddress("fsrPhi", &fsrPhi, &b_fsrPhi);
   fChain->SetBranchAddress("fsrLept", &fsrLept, &b_fsrLept);
   fChain->SetBranchAddress("passIsoPreFSR", &passIsoPreFSR, &b_passIsoPreFSR);
   fChain->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU, &b_p0plus_VAJHU);
   fChain->SetBranchAddress("p0_g1prime2_VAJHU", &p0_g1prime2_VAJHU, &b_p0_g1prime2_VAJHU);
   fChain->SetBranchAddress("p0hplus_VAJHU", &p0hplus_VAJHU, &b_p0hplus_VAJHU);
   fChain->SetBranchAddress("p0minus_VAJHU", &p0minus_VAJHU, &b_p0minus_VAJHU);
   fChain->SetBranchAddress("p0_g1prime2_zgs_VAJHU", &p0_g1prime2_zgs_VAJHU, &b_p0_g1prime2_zgs_VAJHU);
   fChain->SetBranchAddress("p0hplus_zgs_VAJHU", &p0hplus_zgs_VAJHU, &b_p0hplus_zgs_VAJHU);
   fChain->SetBranchAddress("p0minus_zgs_VAJHU", &p0minus_zgs_VAJHU, &b_p0minus_zgs_VAJHU);
   fChain->SetBranchAddress("p0hplus_gsgs_VAJHU", &p0hplus_gsgs_VAJHU, &b_p0hplus_gsgs_VAJHU);
   fChain->SetBranchAddress("p0minus_gsgs_VAJHU", &p0minus_gsgs_VAJHU, &b_p0minus_gsgs_VAJHU);
   fChain->SetBranchAddress("pg1g1prime2_VAJHU", &pg1g1prime2_VAJHU, &b_pg1g1prime2_VAJHU);
   fChain->SetBranchAddress("pg1g2_VAJHU", &pg1g2_VAJHU, &b_pg1g2_VAJHU);
   fChain->SetBranchAddress("pg1g2_pi2_VAJHU", &pg1g2_pi2_VAJHU, &b_pg1g2_pi2_VAJHU);
   fChain->SetBranchAddress("pg1g4_VAJHU", &pg1g4_VAJHU, &b_pg1g4_VAJHU);
   fChain->SetBranchAddress("pg1g4_pi2_VAJHU", &pg1g4_pi2_VAJHU, &b_pg1g4_pi2_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_g1prime2_zgs_VAJHU", &p0plus_zz_g1prime2_zgs_VAJHU, &b_p0plus_zz_g1prime2_zgs_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_g1prime2_zgs_pi2_VAJHU", &p0plus_zz_g1prime2_zgs_pi2_VAJHU, &b_p0plus_zz_g1prime2_zgs_pi2_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_0hplus_zgs_VAJHU", &p0plus_zz_0hplus_zgs_VAJHU, &b_p0plus_zz_0hplus_zgs_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_0minus_zgs_VAJHU", &p0plus_zz_0minus_zgs_VAJHU, &b_p0plus_zz_0minus_zgs_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_0hplus_gsgs_VAJHU", &p0plus_zz_0hplus_gsgs_VAJHU, &b_p0plus_zz_0hplus_gsgs_VAJHU);
   fChain->SetBranchAddress("p0plus_zz_0minus_gsgs_VAJHU", &p0plus_zz_0minus_gsgs_VAJHU, &b_p0plus_zz_0minus_gsgs_VAJHU);
   fChain->SetBranchAddress("p1_VAJHU", &p1_VAJHU, &b_p1_VAJHU);
   fChain->SetBranchAddress("p1_prodIndep_VAJHU", &p1_prodIndep_VAJHU, &b_p1_prodIndep_VAJHU);
   fChain->SetBranchAddress("p1plus_VAJHU", &p1plus_VAJHU, &b_p1plus_VAJHU);
   fChain->SetBranchAddress("p1plus_prodIndep_VAJHU", &p1plus_prodIndep_VAJHU, &b_p1plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2plus_gg_VAJHU", &p2plus_gg_VAJHU, &b_p2plus_gg_VAJHU);
   fChain->SetBranchAddress("p2plus_prodIndep_VAJHU", &p2plus_prodIndep_VAJHU, &b_p2plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2plus_qqb_VAJHU", &p2plus_qqb_VAJHU, &b_p2plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h2plus_gg_VAJHU", &p2h2plus_gg_VAJHU, &b_p2h2plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h2plus_qqb_VAJHU", &p2h2plus_qqb_VAJHU, &b_p2h2plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h2plus_prodIndep_VAJHU", &p2h2plus_prodIndep_VAJHU, &b_p2h2plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h3plus_gg_VAJHU", &p2h3plus_gg_VAJHU, &b_p2h3plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h3plus_qqb_VAJHU", &p2h3plus_qqb_VAJHU, &b_p2h3plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h3plus_prodIndep_VAJHU", &p2h3plus_prodIndep_VAJHU, &b_p2h3plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h4plus_gg_VAJHU", &p2h4plus_gg_VAJHU, &b_p2h4plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h4plus_qqb_VAJHU", &p2h4plus_qqb_VAJHU, &b_p2h4plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h4plus_prodIndep_VAJHU", &p2h4plus_prodIndep_VAJHU, &b_p2h4plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2bplus_gg_VAJHU", &p2bplus_gg_VAJHU, &b_p2bplus_gg_VAJHU);
   fChain->SetBranchAddress("p2bplus_qqb_VAJHU", &p2bplus_qqb_VAJHU, &b_p2bplus_qqb_VAJHU);
   fChain->SetBranchAddress("p2bplus_prodIndep_VAJHU", &p2bplus_prodIndep_VAJHU, &b_p2bplus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h6plus_gg_VAJHU", &p2h6plus_gg_VAJHU, &b_p2h6plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h6plus_qqb_VAJHU", &p2h6plus_qqb_VAJHU, &b_p2h6plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h6plus_prodIndep_VAJHU", &p2h6plus_prodIndep_VAJHU, &b_p2h6plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h7plus_gg_VAJHU", &p2h7plus_gg_VAJHU, &b_p2h7plus_gg_VAJHU);
   fChain->SetBranchAddress("p2h7plus_qqb_VAJHU", &p2h7plus_qqb_VAJHU, &b_p2h7plus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h7plus_prodIndep_VAJHU", &p2h7plus_prodIndep_VAJHU, &b_p2h7plus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2hminus_gg_VAJHU", &p2hminus_gg_VAJHU, &b_p2hminus_gg_VAJHU);
   fChain->SetBranchAddress("p2hminus_qqb_VAJHU", &p2hminus_qqb_VAJHU, &b_p2hminus_qqb_VAJHU);
   fChain->SetBranchAddress("p2hminus_prodIndep_VAJHU", &p2hminus_prodIndep_VAJHU, &b_p2hminus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h9minus_gg_VAJHU", &p2h9minus_gg_VAJHU, &b_p2h9minus_gg_VAJHU);
   fChain->SetBranchAddress("p2h9minus_qqb_VAJHU", &p2h9minus_qqb_VAJHU, &b_p2h9minus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h9minus_prodIndep_VAJHU", &p2h9minus_prodIndep_VAJHU, &b_p2h9minus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p2h10minus_gg_VAJHU", &p2h10minus_gg_VAJHU, &b_p2h10minus_gg_VAJHU);
   fChain->SetBranchAddress("p2h10minus_qqb_VAJHU", &p2h10minus_qqb_VAJHU, &b_p2h10minus_qqb_VAJHU);
   fChain->SetBranchAddress("p2h10minus_prodIndep_VAJHU", &p2h10minus_prodIndep_VAJHU, &b_p2h10minus_prodIndep_VAJHU);
   fChain->SetBranchAddress("p0plus_VAMCFM", &p0plus_VAMCFM, &b_p0plus_VAMCFM);
   fChain->SetBranchAddress("ggzz_VAMCFM", &ggzz_VAMCFM, &b_ggzz_VAMCFM);
   fChain->SetBranchAddress("ggzz_p0plus_VAMCFM", &ggzz_p0plus_VAMCFM, &b_ggzz_p0plus_VAMCFM);
   fChain->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM, &b_bkg_VAMCFM);
   fChain->SetBranchAddress("bkg_prodIndep_VAMCFM", &bkg_prodIndep_VAMCFM, &b_bkg_prodIndep_VAMCFM);
   fChain->SetBranchAddress("pZJJ_VAMCFM", &pZJJ_VAMCFM, &b_pZJJ_VAMCFM);
   fChain->SetBranchAddress("Dgg10_VAMCFM", &Dgg10_VAMCFM, &b_Dgg10_VAMCFM);
   fChain->SetBranchAddress("p0plus_m4l", &p0plus_m4l, &b_p0plus_m4l);
   fChain->SetBranchAddress("p0plus_m4l_ScaleUp", &p0plus_m4l_ScaleUp, &b_p0plus_m4l_ScaleUp);
   fChain->SetBranchAddress("p0plus_m4l_ScaleDown", &p0plus_m4l_ScaleDown, &b_p0plus_m4l_ScaleDown);
   fChain->SetBranchAddress("p0plus_m4l_ResUp", &p0plus_m4l_ResUp, &b_p0plus_m4l_ResUp);
   fChain->SetBranchAddress("p0plus_m4l_ResDown", &p0plus_m4l_ResDown, &b_p0plus_m4l_ResDown);
   fChain->SetBranchAddress("bkg_m4l", &bkg_m4l, &b_bkg_m4l);
   fChain->SetBranchAddress("bkg_m4l_ScaleUp", &bkg_m4l_ScaleUp, &b_bkg_m4l_ScaleUp);
   fChain->SetBranchAddress("bkg_m4l_ScaleDown", &bkg_m4l_ScaleDown, &b_bkg_m4l_ScaleDown);
   fChain->SetBranchAddress("bkg_m4l_ResUp", &bkg_m4l_ResUp, &b_bkg_m4l_ResUp);
   fChain->SetBranchAddress("bkg_m4l_ResDown", &bkg_m4l_ResDown, &b_bkg_m4l_ResDown);
   fChain->SetBranchAddress("pwh_leptonic_VAJHU", &pwh_leptonic_VAJHU, &b_pwh_leptonic_VAJHU);
   fChain->SetBranchAddress("pzh_leptonic_VAJHU", &pzh_leptonic_VAJHU, &b_pzh_leptonic_VAJHU);
   fChain->SetBranchAddress("phjj_VAJHU_highestPTJets", &phjj_VAJHU_highestPTJets, &b_phjj_VAJHU_highestPTJets);
   fChain->SetBranchAddress("pvbf_VAJHU_highestPTJets", &pvbf_VAJHU_highestPTJets, &b_pvbf_VAJHU_highestPTJets);
   fChain->SetBranchAddress("phjj_VAJHU_highestPTJets_up", &phjj_VAJHU_highestPTJets_up, &b_phjj_VAJHU_highestPTJets_up);
   fChain->SetBranchAddress("pvbf_VAJHU_highestPTJets_up", &pvbf_VAJHU_highestPTJets_up, &b_pvbf_VAJHU_highestPTJets_up);
   fChain->SetBranchAddress("phjj_VAJHU_highestPTJets_dn", &phjj_VAJHU_highestPTJets_dn, &b_phjj_VAJHU_highestPTJets_dn);
   fChain->SetBranchAddress("pvbf_VAJHU_highestPTJets_dn", &pvbf_VAJHU_highestPTJets_dn, &b_pvbf_VAJHU_highestPTJets_dn);
   fChain->SetBranchAddress("phjj_VAJHU_bestDjet", &phjj_VAJHU_bestDjet, &b_phjj_VAJHU_bestDjet);
   fChain->SetBranchAddress("pvbf_VAJHU_bestDjet", &pvbf_VAJHU_bestDjet, &b_pvbf_VAJHU_bestDjet);
   fChain->SetBranchAddress("phjj_VAJHU_bestDjet_up", &phjj_VAJHU_bestDjet_up, &b_phjj_VAJHU_bestDjet_up);
   fChain->SetBranchAddress("pvbf_VAJHU_bestDjet_up", &pvbf_VAJHU_bestDjet_up, &b_pvbf_VAJHU_bestDjet_up);
   fChain->SetBranchAddress("phjj_VAJHU_bestDjet_dn", &phjj_VAJHU_bestDjet_dn, &b_phjj_VAJHU_bestDjet_dn);
   fChain->SetBranchAddress("pvbf_VAJHU_bestDjet_dn", &pvbf_VAJHU_bestDjet_dn, &b_pvbf_VAJHU_bestDjet_dn);
   fChain->SetBranchAddress("pAux_vbf_VAJHU", &pAux_vbf_VAJHU, &b_pAux_vbf_VAJHU);
   fChain->SetBranchAddress("pAux_vbf_VAJHU_up", &pAux_vbf_VAJHU_up, &b_pAux_vbf_VAJHU_up);
   fChain->SetBranchAddress("pAux_vbf_VAJHU_dn", &pAux_vbf_VAJHU_dn, &b_pAux_vbf_VAJHU_dn);
   fChain->SetBranchAddress("phj_VAJHU", &phj_VAJHU, &b_phj_VAJHU);
   fChain->SetBranchAddress("phj_VAJHU_up", &phj_VAJHU_up, &b_phj_VAJHU_up);
   fChain->SetBranchAddress("phj_VAJHU_dn", &phj_VAJHU_dn, &b_phj_VAJHU_dn);
   fChain->SetBranchAddress("pwh_hadronic_VAJHU", &pwh_hadronic_VAJHU, &b_pwh_hadronic_VAJHU);
   fChain->SetBranchAddress("pwh_hadronic_VAJHU_up", &pwh_hadronic_VAJHU_up, &b_pwh_hadronic_VAJHU_up);
   fChain->SetBranchAddress("pwh_hadronic_VAJHU_dn", &pwh_hadronic_VAJHU_dn, &b_pwh_hadronic_VAJHU_dn);
   fChain->SetBranchAddress("pzh_hadronic_VAJHU", &pzh_hadronic_VAJHU, &b_pzh_hadronic_VAJHU);
   fChain->SetBranchAddress("pzh_hadronic_VAJHU_up", &pzh_hadronic_VAJHU_up, &b_pzh_hadronic_VAJHU_up);
   fChain->SetBranchAddress("pzh_hadronic_VAJHU_dn", &pzh_hadronic_VAJHU_dn, &b_pzh_hadronic_VAJHU_dn);
   fChain->SetBranchAddress("ptth_VAJHU", &ptth_VAJHU, &b_ptth_VAJHU);
   fChain->SetBranchAddress("ptth_VAJHU_up", &ptth_VAJHU_up, &b_ptth_VAJHU_up);
   fChain->SetBranchAddress("ptth_VAJHU_dn", &ptth_VAJHU_dn, &b_ptth_VAJHU_dn);
   fChain->SetBranchAddress("pbbh_VAJHU", &pbbh_VAJHU, &b_pbbh_VAJHU);
   fChain->SetBranchAddress("pbbh_VAJHU_up", &pbbh_VAJHU_up, &b_pbbh_VAJHU_up);
   fChain->SetBranchAddress("pbbh_VAJHU_dn", &pbbh_VAJHU_dn, &b_pbbh_VAJHU_dn);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetMass", &JetMass, &b_JetMass);
   fChain->SetBranchAddress("JetBTagger", &JetBTagger, &b_JetBTagger);
   fChain->SetBranchAddress("JetIsBtagged", &JetIsBtagged, &b_JetIsBtagged);
   fChain->SetBranchAddress("JetQGLikelihood", &JetQGLikelihood, &b_JetQGLikelihood);
   fChain->SetBranchAddress("JetAxis2", &JetAxis2, &b_JetAxis2);
   fChain->SetBranchAddress("JetMult", &JetMult, &b_JetMult);
   fChain->SetBranchAddress("JetPtD", &JetPtD, &b_JetPtD);
   fChain->SetBranchAddress("JetSigma", &JetSigma, &b_JetSigma);
   fChain->SetBranchAddress("DiJetMass", &DiJetMass, &b_DiJetMass);
   fChain->SetBranchAddress("DiJetDEta", &DiJetDEta, &b_DiJetDEta);
   fChain->SetBranchAddress("DiJetFisher", &DiJetFisher, &b_DiJetFisher);
   fChain->SetBranchAddress("nExtraLep", &nExtraLep, &b_nExtraLep);
   fChain->SetBranchAddress("nExtraZ", &nExtraZ, &b_nExtraZ);
   fChain->SetBranchAddress("ExtraLepPt", &ExtraLepPt, &b_ExtraLepPt);
   fChain->SetBranchAddress("ExtraLepEta", &ExtraLepEta, &b_ExtraLepEta);
   fChain->SetBranchAddress("ExtraLepPhi", &ExtraLepPhi, &b_ExtraLepPhi);
   fChain->SetBranchAddress("ExtraLepLepId", &ExtraLepLepId, &b_ExtraLepLepId);
   fChain->SetBranchAddress("ZXFakeweight", &ZXFakeweight, &b_ZXFakeweight);
   if ( !(input_file_name.Contains("Data")) )
   {
      if ( input_file_name.Contains("ggT") )
      {
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_Nominal", &KFactor_QCD_ggZZ_Nominal, &b_KFactor_QCD_ggZZ_Nominal);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleDn", &KFactor_QCD_ggZZ_PDFScaleDn, &b_KFactor_QCD_ggZZ_PDFScaleDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFScaleUp", &KFactor_QCD_ggZZ_PDFScaleUp, &b_KFactor_QCD_ggZZ_PDFScaleUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleDn", &KFactor_QCD_ggZZ_QCDScaleDn, &b_KFactor_QCD_ggZZ_QCDScaleDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_QCDScaleUp", &KFactor_QCD_ggZZ_QCDScaleUp, &b_KFactor_QCD_ggZZ_QCDScaleUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_AsDn", &KFactor_QCD_ggZZ_AsDn, &b_KFactor_QCD_ggZZ_AsDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_AsUp", &KFactor_QCD_ggZZ_AsUp, &b_KFactor_QCD_ggZZ_AsUp);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaDn", &KFactor_QCD_ggZZ_PDFReplicaDn, &b_KFactor_QCD_ggZZ_PDFReplicaDn);
         fChain->SetBranchAddress("KFactor_QCD_ggZZ_PDFReplicaUp", &KFactor_QCD_ggZZ_PDFReplicaUp, &b_KFactor_QCD_ggZZ_PDFReplicaUp);
      }
      
      if ( input_file_name.Contains("ZZTo4l") )
      {
         fChain->SetBranchAddress("KFactor_EW_qqZZ", &KFactor_EW_qqZZ, &b_KFactor_EW_qqZZ);
         fChain->SetBranchAddress("KFactor_EW_qqZZ_unc", &KFactor_EW_qqZZ_unc, &b_KFactor_EW_qqZZ_unc);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_dPhi", &KFactor_QCD_qqZZ_dPhi, &b_KFactor_QCD_qqZZ_dPhi);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_M", &KFactor_QCD_qqZZ_M, &b_KFactor_QCD_qqZZ_M);
         fChain->SetBranchAddress("KFactor_QCD_qqZZ_Pt", &KFactor_QCD_qqZZ_Pt, &b_KFactor_QCD_qqZZ_Pt);
      }
      
      fChain->SetBranchAddress("genFinalState", &genFinalState, &b_genFinalState);
      fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
      fChain->SetBranchAddress("genHEPMCweight", &genHEPMCweight, &b_genHEPMCweight);
      fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
      fChain->SetBranchAddress("dataMCWeight", &dataMCWeight, &b_dataMCWeight);
      fChain->SetBranchAddress("trigEffWeight", &trigEffWeight, &b_trigEffWeight);
      fChain->SetBranchAddress("overallEventWeight", &overallEventWeight, &b_overallEventWeight);
      fChain->SetBranchAddress("HqTMCweight", &HqTMCweight, &b_HqTMCweight);
      fChain->SetBranchAddress("xsec", &xsec, &b_xsec);
      fChain->SetBranchAddress("genExtInfo", &genExtInfo, &b_genExtInfo);
      fChain->SetBranchAddress("GenHMass", &GenHMass, &b_GenHMass);
      fChain->SetBranchAddress("GenHPt", &GenHPt, &b_GenHPt);
      fChain->SetBranchAddress("GenHRapidity", &GenHRapidity, &b_GenHRapidity);
      fChain->SetBranchAddress("GenZ1Mass", &GenZ1Mass, &b_GenZ1Mass);
      fChain->SetBranchAddress("GenZ1Pt", &GenZ1Pt, &b_GenZ1Pt);
      fChain->SetBranchAddress("GenZ1Phi", &GenZ1Phi, &b_GenZ1Phi);
      fChain->SetBranchAddress("GenZ1Flav", &GenZ1Flav, &b_GenZ1Flav);
      fChain->SetBranchAddress("GenZ2Mass", &GenZ2Mass, &b_GenZ2Mass);
      fChain->SetBranchAddress("GenZ2Pt", &GenZ2Pt, &b_GenZ2Pt);
      fChain->SetBranchAddress("GenZ2Phi", &GenZ2Phi, &b_GenZ2Phi);
      fChain->SetBranchAddress("GenZ2Flav", &GenZ2Flav, &b_GenZ2Flav);
      fChain->SetBranchAddress("GenLep1Pt", &GenLep1Pt, &b_GenLep1Pt);
      fChain->SetBranchAddress("GenLep1Eta", &GenLep1Eta, &b_GenLep1Eta);
      fChain->SetBranchAddress("GenLep1Phi", &GenLep1Phi, &b_GenLep1Phi);
      fChain->SetBranchAddress("GenLep1Id", &GenLep1Id, &b_GenLep1Id);
      fChain->SetBranchAddress("GenLep2Pt", &GenLep2Pt, &b_GenLep2Pt);
      fChain->SetBranchAddress("GenLep2Eta", &GenLep2Eta, &b_GenLep2Eta);
      fChain->SetBranchAddress("GenLep2Phi", &GenLep2Phi, &b_GenLep2Phi);
      fChain->SetBranchAddress("GenLep2Id", &GenLep2Id, &b_GenLep2Id);
      fChain->SetBranchAddress("GenLep3Pt", &GenLep3Pt, &b_GenLep3Pt);
      fChain->SetBranchAddress("GenLep3Eta", &GenLep3Eta, &b_GenLep3Eta);
      fChain->SetBranchAddress("GenLep3Phi", &GenLep3Phi, &b_GenLep3Phi);
      fChain->SetBranchAddress("GenLep3Id", &GenLep3Id, &b_GenLep3Id);
      fChain->SetBranchAddress("GenLep4Pt", &GenLep4Pt, &b_GenLep4Pt);
      fChain->SetBranchAddress("GenLep4Eta", &GenLep4Eta, &b_GenLep4Eta);
      fChain->SetBranchAddress("GenLep4Phi", &GenLep4Phi, &b_GenLep4Phi);
      fChain->SetBranchAddress("GenLep4Id", &GenLep4Id, &b_GenLep4Id);
      fChain->SetBranchAddress("GenAssocLep1Pt", &GenAssocLep1Pt, &b_GenAssocLep1Pt);
      fChain->SetBranchAddress("GenAssocLep1Eta", &GenAssocLep1Eta, &b_GenAssocLep1Eta);
      fChain->SetBranchAddress("GenAssocLep1Phi", &GenAssocLep1Phi, &b_GenAssocLep1Phi);
      fChain->SetBranchAddress("GenAssocLep1Id", &GenAssocLep1Id, &b_GenAssocLep1Id);
      fChain->SetBranchAddress("GenAssocLep2Pt", &GenAssocLep2Pt, &b_GenAssocLep2Pt);
      fChain->SetBranchAddress("GenAssocLep2Eta", &GenAssocLep2Eta, &b_GenAssocLep2Eta);
      fChain->SetBranchAddress("GenAssocLep2Phi", &GenAssocLep2Phi, &b_GenAssocLep2Phi);
      fChain->SetBranchAddress("GenAssocLep2Id", &GenAssocLep2Id, &b_GenAssocLep2Id);
      //   fChain->SetBranchAddress("reweightingweights", &reweightingweights, &b_reweightingweights);
      fChain->SetBranchAddress("LHEPDFScale", &LHEPDFScale, &b_LHEPDFScale);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF1", &LHEweight_QCDscale_muR1_muF1, &b_LHEweight_QCDscale_muR1_muF1);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF2", &LHEweight_QCDscale_muR1_muF2, &b_LHEweight_QCDscale_muR1_muF2);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR1_muF0p5", &LHEweight_QCDscale_muR1_muF0p5, &b_LHEweight_QCDscale_muR1_muF0p5);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF1", &LHEweight_QCDscale_muR2_muF1, &b_LHEweight_QCDscale_muR2_muF1);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF2", &LHEweight_QCDscale_muR2_muF2, &b_LHEweight_QCDscale_muR2_muF2);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR2_muF0p5", &LHEweight_QCDscale_muR2_muF0p5, &b_LHEweight_QCDscale_muR2_muF0p5);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF1", &LHEweight_QCDscale_muR0p5_muF1, &b_LHEweight_QCDscale_muR0p5_muF1);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF2", &LHEweight_QCDscale_muR0p5_muF2, &b_LHEweight_QCDscale_muR0p5_muF2);
      fChain->SetBranchAddress("LHEweight_QCDscale_muR0p5_muF0p5", &LHEweight_QCDscale_muR0p5_muF0p5, &b_LHEweight_QCDscale_muR0p5_muF0p5);
   }
   
   Notify();
}

Bool_t Tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Tree_cxx
