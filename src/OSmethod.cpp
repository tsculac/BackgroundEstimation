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
   
   _s_flavour.push_back("ele");
   _s_flavour.push_back("mu");
   
   DeclareHistos();
}
//============================================================



// Destructor
//====================
OSmethod::~OSmethod()
{
}
//====================


//===============================================================================
void OSmethod::FillHistos( TString input_file_data_name )
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
      if( fabs(Z1Mass - 91.2) < 7)
      {
         // Final event weight
         _event_weight = (_lumi * 1000 * xsec * overallEventWeight) / gen_sum_weights;

         if(LepisID->at(2) && LepCombRelIsoPF->at(2) < 0.35)
         {
            if(fabs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.447) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.447) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
         {
            if(fabs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.447) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), (abs(LepEta->at(2)) < 1.447) ? 0.5 : 1.5 , (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
      }
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
      
      if ( !CRflag ) continue;
      if ( !test_bit(CRflag, CRZLLss) ) continue;
      
      _current_final_state = FindFinalStateZX();
      
      _current_category = categoryMor17(nExtraLep, nExtraZ, nCleanedJetsPt30, nCleanedJetsPt30BTagged_bTagSF, jetQGL,
                                        p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,
                                        p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, p_HadWH_SIG_ghw1_1_JHUGen_JECNominal,
                                        p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, jetPhi, ZZMass, PFMET, true, false);
      
   }
   
}
//===============================================================================



//===============================================================
void OSmethod::DeclareHistos()
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
      
      _histo_name = "Passing_WZsubtracted_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);
      _histo_name = "Failing_WZsubtracted_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = new TH2F(_histo_name,"", 80, 0, 80, 2, 0, 2);

   }

}
//===============================================================

//===============================================================
void OSmethod::SaveHistos( TString file_name)
{
   TFile* fOutHistos = new TFile(file_name, "recreate");
   fOutHistos->cd();
   
   // Copy data histos to total histos, if there is no WZ subtraction this is the final histo for fake rate calculation
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing[Settings::Total][i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::Data][i_flav], 1.);
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
   
   cout << "[INFO] All histograms saved." << endl;
}
//===============================================================

//===============================================================
void OSmethod::GetHistos( TString file_name)
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
      
      _histo_name = "Passing_WZsubtracted_" + _s_flavour.at(i_flav);
      passing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      _histo_name = "Failing_WZsubtracted_" + _s_flavour.at(i_flav);
      failing[Settings::Total][i_flav] = (TH2F*)histo_file->Get(_histo_name);
      
   }
   
   cout << "[INFO] All histograms retrieved from file." << endl;
}
//===============================================================

//===============================================================
void OSmethod::SubtractWZ( bool remove_negative_bins)
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
//      passing[Settings::Total][i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      passing[Settings::Total][i_flav]->Add(passing[Settings::WZ][i_flav], -1.);
      
//      failing[Settings::Total][i_flav]->Add(failing[Settings::Data][i_flav], 1.);
      failing[Settings::Total][i_flav]->Add(failing[Settings::WZ][i_flav], -1.);
   }

   cout << "[INFO] WZ contribution subtracted." << endl;
   
   if ( remove_negative_bins )
   {
      for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
      {
         RemoveNegativeBins( passing[Settings::Total][i_flav] );
         RemoveNegativeBins( passing[Settings::Total][i_flav] );
         
         RemoveNegativeBins( failing[Settings::Total][i_flav] );
         RemoveNegativeBins( failing[Settings::Total][i_flav] );
      }
      cout << "[INFO] Negative bins removed." << endl;
   }
   
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
      //cout << "i_pT_bin: " << i_pT_bin << endl;
      
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
//         cout << "error Y = " << sqrt((temp_error_NP/temp_NP)*(temp_error_NP/temp_NP) + (temp_error_NF/temp_NF)*(temp_error_NF/temp_NF)) << endl;
         
         vector_X[Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::EB][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::EB][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::EB][i_flav].push_back(sqrt((temp_error_NP/temp_NP)*(temp_error_NP/temp_NP) + (temp_error_NF/temp_NF)*(temp_error_NF/temp_NF))); // simple error, consider updating
         
         temp_NP = passing[Settings::Total][i_flav]->IntegralAndError(passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),passing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NP);
         temp_NF = failing[Settings::Total][i_flav]->IntegralAndError(failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin]),failing[Settings::Total][i_flav]->GetXaxis()->FindBin(_pT_bins[i_pT_bin+1]) - 1, 2, 2, temp_error_NF);
         
         vector_X[Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin] + _pT_bins[i_pT_bin + 1])/2);
         vector_Y[Settings::EE][i_flav].push_back(temp_NP/(temp_NP+temp_NF));
         
         vector_EX[Settings::EE][i_flav].push_back((_pT_bins[i_pT_bin + 1] - _pT_bins[i_pT_bin])/2);
         vector_EY[Settings::EE][i_flav].push_back(sqrt((temp_error_NP/temp_NP)*(temp_error_NP/temp_NP) + (temp_error_NF/temp_NF)*(temp_error_NF/temp_NF))); // simple error, consider updating
      }
   }
   
   FR_OS_electron_EB = new TGraphErrors (vector_X[Settings::EB][Settings::ele].size(),
                                         &(vector_X[Settings::EB][Settings::ele][0]),
                                         &(vector_Y[Settings::EB][Settings::ele][0]),
                                         &(vector_EX[Settings::EB][Settings::ele][0]),
                                         &(vector_EY[Settings::EB][Settings::ele][0]));
   FR_OS_electron_EB->SetName("FR_OS_electron_EB");
   
   FR_OS_electron_EE = new TGraphErrors (vector_X[Settings::EE][Settings::ele].size(),
                                         &(vector_X[Settings::EE][Settings::ele][0]),
                                         &(vector_Y[Settings::EE][Settings::ele][0]),
                                         &(vector_EX[Settings::EE][Settings::ele][0]),
                                         &(vector_EY[Settings::EE][Settings::ele][0]));
   FR_OS_electron_EE->SetName("FR_OS_electron_EE");
   
   FR_OS_muon_EB = new TGraphErrors (vector_X[Settings::EB][Settings::mu].size(),
                                     &(vector_X[Settings::EB][Settings::mu][0]),
                                     &(vector_Y[Settings::EB][Settings::mu][0]),
                                     &(vector_EX[Settings::EB][Settings::mu][0]),
                                     &(vector_EY[Settings::EB][Settings::mu][0]));
   FR_OS_muon_EB->SetName("FR_OS_muon_EB");
   
   FR_OS_muon_EE = new TGraphErrors (vector_X[Settings::EE][Settings::mu].size(),
                                     &(vector_X[Settings::EE][Settings::mu][0]),
                                     &(vector_Y[Settings::EE][Settings::mu][0]),
                                     &(vector_EX[Settings::EE][Settings::mu][0]),
                                     &(vector_EY[Settings::EE][Settings::mu][0]));
   FR_OS_muon_EE->SetName("FR_OS_muon_EE");
   
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
void OSmethod::RemoveNegativeBins(TH2F *h) // too slow atm, can not be used like this
{
   for (int i_bin_x = 1; i_bin_x <= 80; i_bin_x++)
   {
      for (int i_bin_y = 1; i_bin_y <= 2; i_bin_x++)
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
   
   return current_process;
}
//==========================================================


//=============================
int OSmethod::FindFinalState()
{
   int final_state = -999;

   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs4e;
      else if ( Z2Flav == +169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
      }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs2mu2e;
      else if ( Z2Flav == +169 )
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

//=============================
int OSmethod::FindFinalStateZX()
{
   int final_state = -999;
   
   if ( Z1Flav == -121 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs4e;
      else if ( Z2Flav == +169 )
         final_state = Settings::fs2e2mu;
      else
         cerr << "[ERROR] in event " << RunNumber << ":" << LumiNumber << ":" << EventNumber << ", Z2Flav = " << Z2Flav << endl;
   }
   else if ( Z1Flav == -169 )
   {
      if ( Z2Flav == +121 )
         final_state = Settings::fs2mu2e;
      else if ( Z2Flav == +169 )
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




