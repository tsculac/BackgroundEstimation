// Include classes
#include "OSmethod.h"

// Constructor
//============================================================
OSmethod::OSmethod():Tree()
{
   _current_process = -999;
   _current_final_state = -999;
   _current_category = -999;
   
   _s_process.push_back("Data");
   _s_process.push_back("WZ");
   _s_process.push_back("qqZZ");
   _s_process.push_back("DY");
   _s_process.push_back("ttbar");
   
   _s_flavour.push_back("ele");
   _s_flavour.push_back("mu");
   
   _s_final_state.push_back("4e");
   _s_final_state.push_back("4mu");
   _s_final_state.push_back("2e2mu");
   _s_final_state.push_back("2mu2e");
   _s_final_state.push_back("4l");
   
   _s_category.push_back("UnTagged");
   _s_category.push_back("VBF1jTagged");
   _s_category.push_back("VBF2jTagged");
   _s_category.push_back("VHLeptTagged");
   _s_category.push_back("VHHadrTagged");
   _s_category.push_back("ttHTagged");
   _s_category.push_back("VHMETTagged");
   _s_category.push_back("Inclusive");
   
   _s_region.push_back("2P2F");
   _s_region.push_back("3P1F");
   
   DeclareFRHistos();
   DeclareDataMCHistos();
   DeclareZXHistos();
}
//============================================================



// Destructor
//====================
OSmethod::~OSmethod()
{
}
//====================


//===============================================================================
void OSmethod::FillFRHistos( TString input_file_data_name )
{
   input_file_data = new TFile("./" + input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      //if ( !test_bit(CRflag, ZL) ) continue;
      if( fabs(Z1Mass - 91.2) < 7 && PFMET < 25.)
      {
         // Final event weight
         _k_factor = calculate_K_factor(input_file_data_name);
         _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;

         if(LepisID->at(2) && LepCombRelIsoPF->at(2) < 0.35)
         {
            if(fabs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
         {
            if(fabs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.479) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.2) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
      }
   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::FillDataMCPlots( TString input_file_data_name )
{
   input_file_data = new TFile("./" + input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , true);
   
   _current_process = find_current_process(input_file_data_name);
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      
      _current_category = categoryMor17(nExtraLep, nExtraZ, nCleanedJetsPt30, nCleanedJetsPt30BTagged_bTagSF, jetQGL,
                                        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
                                        p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
                                        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, jetPhi, ZZMass, PFMET, true, false);
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      
      if ( test_bit(CRflag, CRZLLos_2P2F) ) histos_1D[Settings::reg2P2F][_current_process][_current_final_state][_current_category]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
      if ( test_bit(CRflag, CRZLLos_3P1F) ) histos_1D[Settings::reg3P1F][_current_process][_current_final_state][_current_category]->Fill(ZZMass, (_current_process == Settings::Data) ? 1 :  _event_weight);
   
   } // END events loop
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================



//===============================================================================
void OSmethod::MakeHistogramsZX( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   
   input_file_data = new TFile("./" + input_file_data_name);
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , false);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_2P2F)) && !(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      
      _current_category = categoryMor17(nExtraLep, nExtraZ, nCleanedJetsPt30, nCleanedJetsPt30BTagged_bTagSF, jetQGL,
                                        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
                                        p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
                                        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, jetPhi, ZZMass, PFMET, true, false);

      if ( test_bit(CRflag, CRZLLos_2P2F) )
      {
         _f3 = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         _f4 = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
//         cout << "===============" << endl;
//         cout << "f3 = " << _f3 << endl;
//         cout << "f4 = " << _f4 << endl;
//         cout << "weight = " << (_f3/(1-_f3))*(_f4/(1-_f4)) << endl;
         h_from2P2F_SR[_current_final_state][_current_category]->Fill(ZZMass, (_f3/(1-_f3))*(_f4/(1-_f4)) );
         h_from2P2F_3P1F[_current_final_state][_current_category]->Fill(ZZMass, (_f3/(1-_f3))+(_f4/(1-_f4)) );
      }
      if ( test_bit(CRflag, CRZLLos_3P1F) )
      {
         if(LepisID->at(3) && LepCombRelIsoPF->at(3) < 0.35) _f4 = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
         else _f4 = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));
         
         h_from3P1F_SR[_current_final_state][_current_category]->Fill(ZZMass,_f4/(1-_f4) );
      }
      
   }
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================

//===============================================================================
void OSmethod::MakeZXMCContribution( TString input_file_data_name, TString  input_file_FR_name )
{
   
   FakeRates *FR = new FakeRates( input_file_FR_name );
   input_file_data = new TFile("./" + input_file_data_name);
   
   hCounters = (TH1F*)input_file_data->Get("CRZLLTree/Counters");
   gen_sum_weights = (Long64_t)hCounters->GetBinContent(40);
   
   input_tree_data = (TTree*)input_file_data->Get("CRZLLTree/candTree");
   Init( input_tree_data, input_file_data_name , false);
   
   
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (!(test_bit(CRflag, CRZLLos_3P1F))) continue;
      
      _current_final_state = FindFinalState();
      
      _current_category = categoryMor17(nExtraLep, nExtraZ, nCleanedJetsPt30, nCleanedJetsPt30BTagged_bTagSF, jetQGL,
                                        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
                                        p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
                                        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, jetPhi, ZZMass, PFMET, true, false);
      
      _k_factor = calculate_K_factor(input_file_data_name);
      _event_weight = (_lumi * 1000 * xsec * _k_factor * overallEventWeight) / gen_sum_weights;
      
      if(LepisID->at(3) && LepCombRelIsoPF->at(3) < 0.35) _f4 = FR->GetFakeRate(LepPt->at(3),LepEta->at(3),LepLepId->at(3));
      else _f4 = FR->GetFakeRate(LepPt->at(2),LepEta->at(2),LepLepId->at(2));

      h_from3P1F_SR_ZZonly[_current_final_state][_current_category]->Fill(ZZMass, _event_weight * (_f4/(1-_f4)) );
      
   }
   
   cout << "[INFO] Processing of " << input_file_data_name << " done." << endl;
}
//===============================================================================




//===============================================================
void OSmethod::DeclareFRHistos()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

   }

}
//===============================================================

//===============================================================
void OSmethod::DeclareDataMCHistos()
{
   for (int i_reg = 0; i_reg < num_of_regions; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
               histos_1D[i_reg][i_proc][i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
            }
         }
      }
   }
   
}
//===============================================================

//===============================================================
void OSmethod::DeclareZXHistos()
{
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
      {
         _histo_name = "h_from2P2F_SR_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         h_from2P2F_SR[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
         
         _histo_name = "h_from2P2F_3P1F_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         h_from2P2F_3P1F[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
         
         _histo_name = "h_from3P1F_SR_final_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         h_from3P1F_SR_final[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
         
         _histo_name = "h_from3P1F_SR_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         h_from3P1F_SR[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
         
         _histo_name = "h_from3P1F_SR_ZZonly_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         h_from3P1F_SR_ZZonly[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
         
         _histo_name = "ZX_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         _histo_labels = ";" + Plots::M4l().var_X_label + ";" + Plots::M4l().var_Y_label;
         histos_ZX[i_fs][i_cat] = new TH1F(_histo_name, _histo_labels, Plots::M4l().var_N_bin, Plots::M4l().var_min, Plots::M4l().var_max);
      }
   }
}
//===============================================================

//===============================================================
void OSmethod::SaveFRHistos( TString file_name,  bool remove_negative_bins)
{
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   // Copy data histos to total histos, if there is no WZ subtraction this is the final histo for fake rate calculation
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::Data][i_flav], 1.);
   }
   
   if ( remove_negative_bins ) // Set negative bins to zero
   {
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         RemoveNegativeBins2D( passing[Settings::Total][i_flav] );
         RemoveNegativeBins2D( passing[Settings::Total][i_flav] );
         
         RemoveNegativeBins2D( failing[Settings::Total][i_flav] );
         RemoveNegativeBins2D( failing[Settings::Total][i_flav] );
      }
      cout << "[INFO] Negative bins removed." << endl;
   }

   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         passing[i_proc][i_flav]->Write();
         failing[i_proc][i_flav]->Write();
      }
      
   passing[Settings::Total][i_flav]->Write();
   failing[Settings::Total][i_flav]->Write();
      
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All FakeRate histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveDataMCHistos( TString file_name )
{
   FillDataMCInclusive();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_reg = 0; i_reg < num_of_regions; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               histos_1D[i_reg][i_proc][i_fs][i_cat]->Write();
            }
         }
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Data/MC histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillDataMCInclusive( )
{
   for (int i_reg = 0; i_reg < num_of_regions; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
            {
               histos_1D[i_reg][i_proc][i_fs][Settings::inclusive]->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
               histos_1D[i_reg][i_proc][Settings::fs4l][i_cat]    ->Add(histos_1D[i_reg][i_proc][i_fs][i_cat]);
            }
         }
      }
   }
   
   for (int i_reg = 0; i_reg < num_of_regions; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
         {
            histos_1D[i_reg][i_proc][Settings::fs4l][Settings::inclusive]->Add(histos_1D[i_reg][i_proc][i_fs][Settings::inclusive]);
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SaveZXHistos( TString file_name )
{
   FillZXInclusive();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
      {
         h_from2P2F_SR[i_fs][i_cat]->Write();
         h_from2P2F_3P1F[i_fs][i_cat]->Write();
         h_from3P1F_SR_final[i_fs][i_cat]->Write();
         h_from3P1F_SR[i_fs][i_cat]->Write();
         h_from3P1F_SR_ZZonly[i_fs][i_cat]->Write();
         histos_ZX[i_fs][i_cat]->Write();
      }
   }
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] All Z+X histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::FillZXInclusive( )
{
   for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat < Settings::inclusive; i_cat++)
      {
         h_from2P2F_SR[i_fs][Settings::inclusive]->Add(h_from2P2F_SR[i_fs][i_cat]);
         h_from2P2F_SR[Settings::fs4l][i_cat]    ->Add(h_from2P2F_SR[i_fs][i_cat]);
         
         h_from2P2F_3P1F[i_fs][Settings::inclusive]->Add(h_from2P2F_3P1F[i_fs][i_cat]);
         h_from2P2F_3P1F[Settings::fs4l][i_cat]    ->Add(h_from2P2F_3P1F[i_fs][i_cat]);
         
         h_from3P1F_SR_final[i_fs][Settings::inclusive]->Add(h_from3P1F_SR_final[i_fs][i_cat]);
         h_from3P1F_SR_final[Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR_final[i_fs][i_cat]);
         
         h_from3P1F_SR[i_fs][Settings::inclusive]->Add(h_from3P1F_SR[i_fs][i_cat]);
         h_from3P1F_SR[Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR[i_fs][i_cat]);
         
         h_from3P1F_SR_ZZonly[i_fs][Settings::inclusive]->Add(h_from3P1F_SR_ZZonly[i_fs][i_cat]);
         h_from3P1F_SR_ZZonly[Settings::fs4l][i_cat]    ->Add(h_from3P1F_SR_ZZonly[i_fs][i_cat]);
      }
   }

   for (int i_fs = 0; i_fs < Settings::fs4l; i_fs++)
   {
      h_from2P2F_SR[Settings::fs4l][Settings::inclusive]->Add(h_from2P2F_SR[i_fs][Settings::inclusive]);
      h_from2P2F_3P1F[Settings::fs4l][Settings::inclusive]->Add(h_from2P2F_3P1F[i_fs][Settings::inclusive]);
      h_from3P1F_SR_final[Settings::fs4l][Settings::inclusive]->Add(h_from3P1F_SR_final[i_fs][Settings::inclusive]);
      h_from3P1F_SR[Settings::fs4l][Settings::inclusive]->Add(h_from3P1F_SR[i_fs][Settings::inclusive]);
      h_from3P1F_SR_ZZonly[Settings::fs4l][Settings::inclusive]->Add(h_from3P1F_SR_ZZonly[i_fs][Settings::inclusive]);
      
   }
   
   for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat <= Settings::inclusive; i_cat++)
      {
         h_from3P1F_SR_final[i_fs][i_cat]->Add(h_from3P1F_SR[i_fs][i_cat], 1.);
         h_from3P1F_SR_final[i_fs][i_cat]->Add(h_from3P1F_SR_ZZonly[i_fs][i_cat], -1.);
         h_from3P1F_SR_final[i_fs][i_cat]->Add(h_from2P2F_SR[i_fs][i_cat], -2.);
      }
   }
   
   for (int i_fs = 0; i_fs <= Settings::fs4l; i_fs++)
   {
      for (int i_cat = 0; i_cat <= Settings::inclusive; i_cat++)
      {
         histos_ZX[i_fs][i_cat]->Add(h_from3P1F_SR_final[i_fs][i_cat], 1.);
         histos_ZX[i_fs][i_cat]->Add(h_from2P2F_SR[i_fs][i_cat], 1.);
      }
   }
   

   
   cout << "[INFO] All Z+X histograms summed." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetFRHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = (TH2F*)histo_file->Get(_histo_name);
         
      }
      
      _histo_name = "Passing_Total_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      _histo_name = "Failing_Total_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      
   }
   
   cout << "[INFO] All FakeRate histograms retrieved from file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetDataMCHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_reg = 0; i_reg < num_of_regions; i_reg ++)
   {
      for (int i_proc = 0; i_proc < Settings::Total; i_proc++)
      {
         for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
         {
            for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
            {
               _histo_name = "M4l_" + _s_region.at(i_reg) + "_" + _s_process.at(i_proc) + "_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
               histos_1D[i_reg][i_proc][i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
            }
         }
      }
   }
   
   cout << "[INFO] All Data/MC histograms retrieved from file." << endl;
}

//===============================================================

//===============================================================
void OSmethod::GetZXHistos( TString file_name)
{
   TFile* histo_file = TFile::Open(file_name);
   
   for (int i_fs = 0; i_fs < num_of_final_states; i_fs++)
   {
      for (int i_cat = 0; i_cat < num_of_categories; i_cat++)
      {
         _histo_name = "h_from2P2F_SR_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         h_from2P2F_SR[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         
         _histo_name = "h_from2P2F_3P1F_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         h_from2P2F_3P1F[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         
         _histo_name = "h_from3P1F_SR_final_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         h_from3P1F_SR_final[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         
         _histo_name = "h_from3P1F_SR_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         h_from3P1F_SR[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         
         _histo_name = "h_from3P1F_SR_ZZonly_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         h_from3P1F_SR_ZZonly[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
         
         _histo_name = "ZX_" + _s_final_state.at(i_fs) + "_" + _s_category.at(i_cat);
         histos_ZX[i_fs][i_cat] = (TH1F*)histo_file->Get(_histo_name);
      }
   }
   
   cout << "[INFO] All Z+X histograms retrieved from file." << endl;
}

//===============================================================

//========================================================================================================
void OSmethod::PlotDataMC_2P2F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("2P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive]   ->SetFillColor(kMagenta);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive]   ->SetFillColor(kGreen-1);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive]->SetFillColor(kBlue-2);
      
      histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive]   ->SetLineColor(kMagenta);
      histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive]   ->SetLineColor(kGreen-1);
      histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive]->SetLineColor(kBlue-2);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->SetMarkerSize(0.8);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->SetMarkerStyle(20);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
      stack->Add(histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive]);
   
      stack->Draw("HIST");
      
      float data_max = histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->GetBinContent(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->GetBinErrorUp(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.1);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_2P2F("right",histos_1D[Settings::reg2P2F][Settings::Data][i_fs][Settings::inclusive],histos_1D[Settings::reg2P2F][Settings::WZ][i_fs][Settings::inclusive],histos_1D[Settings::reg2P2F][Settings::qqZZ][i_fs][Settings::inclusive],histos_1D[Settings::reg2P2F][Settings::DY][i_fs][Settings::inclusive],histos_1D[Settings::reg2P2F][Settings::ttbar][i_fs][Settings::inclusive]);
      legend->Draw();

      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_2P2F_" + _s_final_state.at(i_fs) + "_" + _s_category.at(Settings::inclusive);
      SavePlots(c, _out_file_name);

   }
}
//========================================================================================================


//========================================================================================================
void OSmethod::PlotDataMC_3P1F( TString variable_name, TString folder )
{
   TCanvas *c;
   c = new TCanvas("3P2F", variable_name, 600, 600);
   
   if ( GetVarLogX( variable_name) ) c->SetLogx();
   if ( GetVarLogY( variable_name) ) c->SetLogy();
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive]   ->SetFillColor(kMagenta);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive] ->SetFillColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive]   ->SetFillColor(kGreen-1);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive]->SetFillColor(kBlue-2);
      
      histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive]   ->SetLineColor(kMagenta);
      histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive] ->SetLineColor(kCyan+1);
      histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive]   ->SetLineColor(kGreen-1);
      histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive]->SetLineColor(kBlue-2);
      
      h_from2P2F_3P1F[i_fs][Settings::inclusive]->SetLineColor(kRed);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->SetMarkerSize(0.8);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->SetMarkerStyle(20);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->SetBinErrorOption(TH1::kPoisson);
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->SetLineColor(kBlack);
      
      THStack *stack = new THStack( "stack", "stack" );
      stack->Add(histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive]);
      stack->Add(histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive]);
      
      stack->Draw("HIST");
      
      h_from2P2F_3P1F[i_fs][Settings::inclusive]->Draw("HIST SAME");
      
      float data_max = histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->GetBinContent(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->GetMaximumBin());
      float data_max_error = histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->GetBinErrorUp(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->GetMaximumBin());
      
      stack->SetMinimum(1e-5);
      stack->SetMaximum((data_max + data_max_error)*1.35);
      
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      stack->GetXaxis()->SetTitle(_fs_label);
      stack->GetXaxis()->SetTitleSize(0.04);
      stack->GetXaxis()->SetLabelSize(0.04);
      stack->GetYaxis()->SetTitle(histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->GetYaxis()->GetTitle());
      stack->GetYaxis()->SetTitleSize(0.04);
      stack->GetYaxis()->SetLabelSize(0.04);
      
      stack->GetXaxis()->SetTitleOffset(1.2);
      stack->GetYaxis()->SetTitleOffset(1.25);
      
      histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive]->Draw("SAME p E1 X0");
      
      TLegend *legend;
      legend  = CreateLegend_3P1F("right",histos_1D[Settings::reg3P1F][Settings::Data][i_fs][Settings::inclusive],h_from2P2F_3P1F[i_fs][Settings::inclusive],histos_1D[Settings::reg3P1F][Settings::WZ][i_fs][Settings::inclusive],histos_1D[Settings::reg3P1F][Settings::qqZZ][i_fs][Settings::inclusive],histos_1D[Settings::reg3P1F][Settings::DY][i_fs][Settings::inclusive],histos_1D[Settings::reg3P1F][Settings::ttbar][i_fs][Settings::inclusive]);
      legend->Draw();
      
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + variable_name + "_3P1F_" + _s_final_state.at(i_fs) + "_" + _s_category.at(Settings::inclusive);
      SavePlots(c, _out_file_name);
      
   }
}
//========================================================================================================

//========================================================================================================
void OSmethod::PlotZXContributions( TString folder )
{
   TCanvas *c;
   c = new TCanvas("c", "c", 600, 600);
   
   for( int i_fs = 0; i_fs < Settings::fs4l ; i_fs++ )
   {
      h_from3P1F_SR[i_fs][Settings::inclusive]       ->SetLineColor(kBlue);
      h_from2P2F_SR[i_fs][Settings::inclusive]       ->SetLineColor(kYellow);
      h_from3P1F_SR_final[i_fs][Settings::inclusive] ->SetLineColor(kBlack);
      h_from3P1F_SR_ZZonly[i_fs][Settings::inclusive]->SetLineColor(kRed);
      histos_ZX[i_fs][Settings::inclusive]           ->SetLineColor(kGreen);
      
      h_from3P1F_SR[i_fs][Settings::inclusive]->SetMinimum(-1); // For visualisation of negative bins
      
      h_from3P1F_SR[i_fs][Settings::inclusive]       ->Draw("HIST");
      h_from2P2F_SR[i_fs][Settings::inclusive]       ->Draw("HIST SAME");
      h_from3P1F_SR_final[i_fs][Settings::inclusive] ->Draw("HIST SAME");
      h_from3P1F_SR_ZZonly[i_fs][Settings::inclusive]->Draw("HIST SAME");
      histos_ZX[i_fs][Settings::inclusive]           ->Draw("HIST SAME");
      
      TString _fs_label;
      if ( i_fs == Settings::fs4e) _fs_label = "m_{4#font[12]{e}} (GeV)";
      if ( i_fs == Settings::fs4mu) _fs_label = "m_{4#font[12]{#mu}} (GeV)";
      if ( i_fs == Settings::fs2e2mu) _fs_label = "m_{2#font[12]{e}2#font[12]{#mu}} (GeV)";
            if ( i_fs == Settings::fs2mu2e) _fs_label = "m_{2#font[12]{#mu}2#font[12]{e}} (GeV)";
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetXaxis()->SetTitle(_fs_label);
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetXaxis()->SetTitleSize(0.04);
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetXaxis()->SetLabelSize(0.04);
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetYaxis()->SetTitle(h_from2P2F_SR[i_fs][Settings::inclusive]->GetYaxis()->GetTitle());
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetYaxis()->SetTitleSize(0.04);
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetYaxis()->SetLabelSize(0.04);
      
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetXaxis()->SetTitleOffset(1.2);
      h_from2P2F_SR[i_fs][Settings::inclusive]->GetYaxis()->SetTitleOffset(1.25);
      
      TLegend *legend;
      legend  = CreateLegend_ZXcontr( "right", h_from2P2F_SR[i_fs][Settings::inclusive], h_from3P1F_SR[i_fs][Settings::inclusive],h_from3P1F_SR_ZZonly[i_fs][Settings::inclusive],h_from3P1F_SR_final[i_fs][Settings::inclusive],histos_ZX[i_fs][Settings::inclusive] );
      legend->Draw();
      
      // Draw lumi
      CMS_lumi *lumi = new CMS_lumi;
      lumi->set_lumi(c, _lumi, 0);
      
      TString _out_file_name;
      _out_file_name = folder + "/" + "ZX_Contributions_" + _s_final_state.at(i_fs) + "_" + _s_category.at(Settings::inclusive);
      SavePlots(c, _out_file_name);
      
   }
}
//========================================================================================================


//===============================================================
void OSmethod::SubtractWZ()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::WZ][i_flav], -1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::WZ][i_flav], -1.);
   }

   cout << "[INFO] WZ contribution subtracted." << endl;
   
}
//===============================================================

//===============================================================
void OSmethod::ProduceFakeRates( TString file_name )
{
   for(int i_pT_bin = 0; i_pT_bin < _n_pT_bins - 1; i_pT_bin++ )
   {
      double temp_NP = 0;
      double temp_NF = 0;
      double temp_error_x = 0;
      double temp_error_NP = 0;
      double temp_error_NF = 0;
      
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         if ( i_flav == Settings::ele && i_pT_bin == 0) continue; // electrons do not have 5 - 7 GeV bin
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
//         cout << "========================================" << endl;
//         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
//         cout << "NP = " << temp_NP << endl;
//         cout << "error NP = " << temp_error_NP << endl;
//         cout << "NF = " << temp_NF << endl;
//         cout << "error NF = " << temp_error_NF << endl;
//         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
//         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
//         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
//         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::corrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::corrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::corrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         // Just for fake rate plots calculate the same for histograms without WZ subtraction
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 1, 1, temp_error_NF);
         
         //         cout << "========================================" << endl;
         //         cout << "pT bin = " << _pT_bins[i_pT_bin] << endl;
         //         cout << "NP = " << temp_NP << endl;
         //         cout << "error NP = " << temp_error_NP << endl;
         //         cout << "NF = " << temp_NF << endl;
         //         cout << "error NF = " << temp_error_NF << endl;
         //         cout << "X = " << (_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2 << endl;
         //         cout << "error X = " << (_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2 << endl;
         //         cout << "Y = " << temp_NP/(temp_NP+temp_NF) << endl;
         //         cout << "error Y = " << sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)) << endl;
         
         vector_X[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EB][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));
         
         temp_NP = passing[Settings::Data][i_flav]->IntegralAndError(passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Data][i_flav]->IntegralAndError(failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Data][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::uncorrected][Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::uncorrected][Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::uncorrected][Settings::EE][i_flav].push_back(sqrt(pow((temp_NF/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NP,2) + pow((temp_NP/pow(temp_NF+temp_NP,2)),2)*pow(temp_error_NF,2)));

      }
   }
   
   FR_OS_electron_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB->SetName("FR_OS_electron_EB");
   
   FR_OS_electron_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::corrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::corrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE->SetName("FR_OS_electron_EE");
   
   FR_OS_muon_EB = new TGraphErrors (vector_X[Settings::corrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB->SetName("FR_OS_muon_EB");
   
   FR_OS_muon_EE = new TGraphErrors (vector_X[Settings::corrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::corrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::corrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE->SetName("FR_OS_muon_EE");
   
   FR_OS_electron_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB_unc->SetName("FR_OS_electron_EB_unc");
   
   FR_OS_electron_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::uncorrected][Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::uncorrected][Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE_unc->SetName("FR_OS_electron_EE_unc");
   
   FR_OS_muon_EB_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB_unc->SetName("FR_OS_muon_EB_unc");
   
   FR_OS_muon_EE_unc = new TGraphErrors (vector_X[Settings::uncorrected][Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::uncorrected][Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::uncorrected][Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE_unc->SetName("FR_OS_muon_EE_unc");
   
   PlotFR();
   
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   FR_OS_electron_EB->Write();
   FR_OS_electron_EE->Write();
   FR_OS_muon_EB->Write();
   FR_OS_muon_EE->Write();
   
   fOutHistos->Close();
   delete fOutHistos;
   
   cout << "[INFO] Fake rates produced and stored in a file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::PlotFR()
{
   TCanvas *c_ele, *c_mu;
   c_ele = new TCanvas("FR_ele", "FR_ele", 600, 600);
   c_mu  = new TCanvas("FR_mu", "FR_mu", 600, 600);
   
   mg_electrons = new TMultiGraph();
   mg_muons = new TMultiGraph();
   
   mg_electrons->Add(FR_OS_electron_EB);
   FR_OS_electron_EB->SetLineColor(kBlue);
   FR_OS_electron_EB->SetLineStyle(2);
   FR_OS_electron_EB->SetMarkerSize(0);
   FR_OS_electron_EB->SetTitle("barel corrected");
   mg_electrons->Add(FR_OS_electron_EE);
   FR_OS_electron_EE->SetLineColor(kRed);
   FR_OS_electron_EE->SetLineStyle(2);
   FR_OS_electron_EE->SetMarkerSize(0);
   FR_OS_electron_EE->SetTitle("endcap corrected");
   mg_electrons->Add(FR_OS_electron_EB_unc);
   FR_OS_electron_EB_unc->SetLineColor(kBlue);
   FR_OS_electron_EB_unc->SetLineStyle(1);
   FR_OS_electron_EB_unc->SetMarkerSize(0);
   FR_OS_electron_EB_unc->SetTitle("barel uncorrected");
   mg_electrons->Add(FR_OS_electron_EE_unc);
   FR_OS_electron_EE_unc->SetLineColor(kRed);
   FR_OS_electron_EE_unc->SetLineStyle(1);
   FR_OS_electron_EE_unc->SetMarkerSize(0);
   FR_OS_electron_EE_unc->SetTitle("endcap uncorrected");
   
   mg_muons->Add(FR_OS_muon_EB);
   FR_OS_muon_EB->SetLineColor(kBlue);
   FR_OS_muon_EB->SetLineStyle(2);
   FR_OS_muon_EB->SetMarkerSize(0);
   FR_OS_muon_EB->SetTitle("barel corrected");
   mg_muons->Add(FR_OS_muon_EE);
   FR_OS_muon_EE->SetLineColor(kRed);
   FR_OS_muon_EE->SetLineStyle(2);
   FR_OS_muon_EE->SetMarkerSize(0);
   FR_OS_muon_EE->SetTitle("endcap corrected");
   mg_muons->Add(FR_OS_muon_EB_unc);
   FR_OS_muon_EB_unc->SetLineColor(kBlue);
   FR_OS_muon_EB_unc->SetLineStyle(1);
   FR_OS_muon_EB_unc->SetMarkerSize(0);
   FR_OS_muon_EB_unc->SetTitle("barel uncorrected");
   mg_muons->Add(FR_OS_muon_EE_unc);
   FR_OS_muon_EE_unc->SetLineColor(kRed);
   FR_OS_muon_EE_unc->SetLineStyle(1);
   FR_OS_muon_EE_unc->SetMarkerSize(0);
   FR_OS_muon_EE_unc->SetTitle("endcap uncorrected");
   
   
   gStyle->SetEndErrorSize(0);
   
   TLegend *leg_ele,*leg_mu;

   c_ele->cd();
   mg_electrons->Draw("AP");
   mg_electrons->SetMaximum(0.35);
   leg_ele = CreateLegend_FR("left",FR_OS_electron_EB_unc,FR_OS_electron_EB,FR_OS_electron_EE_unc,FR_OS_electron_EE);
   leg_ele->Draw();
   SavePlots(c_ele, "Plots/FR_electrons");
   
   c_mu->cd();
   mg_muons->Draw("AP");
   mg_muons->SetMaximum(0.35);
   leg_mu = CreateLegend_FR("left",FR_OS_muon_EB_unc,FR_OS_muon_EB,FR_OS_muon_EE_unc,FR_OS_muon_EE);
   leg_mu->Draw();
   SavePlots(c_mu, "Plots/FR_muons");
   
}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins1D(TH1F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      if( h->GetBinContent(i_bin_x) < 0.) h->SetBinContent(i_bin_x, 0);
   }
   
}
//===============================================================

//===============================================================
void OSmethod::RemoveNegativeBins2D(TH2F *h)
{
   for (int i_bin_x = 1; i_bin_x <= h->GetXaxis()->GetNbins(); i_bin_x++)
   {
      for (int i_bin_y = 1; i_bin_y <= h->GetYaxis()->GetNbins(); i_bin_y++)
      {
         if( h->GetBinContent(i_bin_x,i_bin_y) < 0.) h->SetBinContent(i_bin_x,i_bin_y,0);
      }
      
   }
   
}
//===============================================================

//===============================================================
void OSmethod::Set_pT_binning(int size, float *bins)
{
   _n_pT_bins = size;

   for (int i = 0; i < size; i++)
   {
      _pT_bins[i] = bins[i];
   }
}
//===============================================================

//===============================================================
void OSmethod::SetLumi(float lumi)
{
   _lumi = lumi;
}
//===============================================================


//==========================================================
int OSmethod::find_current_process( TString input_file_name )
{
   
   int current_process = -999;
   
   // Assign dataset to correct process
   if ( input_file_name.Contains("Data") )           current_process = Settings::Data;
   if ( input_file_name.Contains("WZ") )             current_process = Settings::WZ;
   if ( input_file_name.Contains("ZZTo4l") )         current_process = Settings::qqZZ;
   if ( input_file_name.Contains("DYJetsToLL") )     current_process = Settings::DY;
   if ( input_file_name.Contains("TTJets") )         current_process = Settings::ttbar;
   if ( input_file_name.Contains("TTTo2L2Nu") )      current_process = Settings::ttbar;
   
   return current_process;
}
//==========================================================


//=============================
int OSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == -121 )
         final_state = Settings::fs4e;
      else if ( Z2Flav == -169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == -121 )
         final_state = Settings::fs2mu2e;
      else if ( Z2Flav == -169 )
         final_state = Settings::fs4mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else
   {
      cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z1Flav = " << Z1Flav << endl;
   }
   
   return final_state;
}
//=============================


//=================================
float OSmethod::calculate_K_factor(TString input_file_name)
{
   
   float k_factor = 1;
   
   if ( input_file_name.Contains("ZZTo4l"))
   {
      k_factor = KFactor_EW_qqZZ * KFactor_QCD_qqZZ_M; // As of Moriond2016
   }
   else if ( input_file_name.Contains("ggTo"))
   {
      k_factor = KFactor_QCD_ggZZ_Nominal; // as of Moriond2016
   }
   return k_factor;
}
//=================================

//===================================================
bool OSmethod::GetVarLogX ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_x);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort;
      return bool(Plots::M4l().var_log_x);
   }
}
//===================================================



//===================================================
bool OSmethod::GetVarLogY ( TString variable_name )
{
   //=============
   // M4l
   //=============
   if(variable_name == "M4l")                return bool(Plots::M4l().var_log_y);

   else
   {
      cout << "[ERROR] Wrong variable name choosen!" << endl;
      abort;
      return bool(Plots::M4l().var_log_y);
   }
}
//===================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_FR( string position, TGraphErrors *EB_unc, TGraphErrors *EB_cor,TGraphErrors *EE_unc,TGraphErrors *EE_cor )
{
   TLegend *leg;
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( EB_unc, "barrel uncorrected", "l" );
   leg->AddEntry( EB_cor, "barrel corrected","l");
   leg->AddEntry( EE_unc, "endcap uncorrected", "l" );
   leg->AddEntry( EE_cor, "endcap corrected", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_ZXcontr( string position, TH1F *h_2P2F_SR, TH1F *h_3P1F_SR,TH1F *h_3P1F_ZZ,TH1F *h_3P1F_SR_final,TH1F *total )
{
   TLegend *leg;
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   
   leg->AddEntry( h_2P2F_SR, "2P2F", "l" );
   leg->AddEntry( h_3P1F_SR, "3P1F w/o removal","l");
   leg->AddEntry( h_3P1F_ZZ, "3P1F ZZ contr.", "l" );
   leg->AddEntry( h_3P1F_SR_final, "3P1F final", "l" );
   leg->AddEntry( total, "Z+X final", "l" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_2P2F( string position, TH1F *data, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================

//=========================================================================================================
TLegend* OSmethod::CreateLegend_3P1F( string position, TH1F *data, TH1F *h_2P2F, TH1F *WZ,TH1F *qqZZ,TH1F *DY,TH1F *ttbar )
{
   TLegend *leg;
   if(position == "right") leg = new TLegend( .64, .65, .97, .9 );
   else if(position == "left") leg = new TLegend(.18,.65,.51,.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   
   leg->AddEntry( data, "Data", "p" );
   leg->AddEntry( h_2P2F, "2P2F extr.", "l" );
   leg->AddEntry( WZ,"WZ","f");
   leg->AddEntry( qqZZ, "Z#gamma*, ZZ", "f" );
   leg->AddEntry( DY, "Z + jets", "f" );
   leg->AddEntry( ttbar, "t#bar{t} + jets", "f" );
   
   return leg;
}
//=========================================================================================================



//=======================================
void OSmethod::SavePlots( TCanvas *c, TString name)
{
   c->SaveAs(name + ".pdf");
   c->SaveAs(name + ".root");
   c->SaveAs(name + ".eps");
   gSystem->Exec("convert -density 300 -quality 100 " + name + ".eps " + name + ".png");
}
//=======================================





