// Include classes
#include "OSmethod.h"

// Constructor
//============================================================
OSmethod::OSmethod():Tree()
{
   _lumi = 36.8;
   _current_process = -999;
   _current_final_state = -999;
   _current_category = -999;
   
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
   Init( input_tree_data, input_file_data_name );
   
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
            if(fabs(LepLepId->at(2)) == 11 ) passing[_current_process][Settings::ele]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) passing[_current_process][Settings::mu]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
         else
         {
            if(fabs(LepLepId->at(2)) == 11 ) failing[_current_process][Settings::ele]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
            else if(fabs(LepLepId->at(2)) == 13 ) failing[_current_process][Settings::mu]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
         }
      }
   } // END events loop

}
//===============================================================================

//===============================================================
void OSmethod::DeclareHistos()
{
   _s_process.push_back("Data");
   _s_process.push_back("WZ");
   
   _s_flavour.push_back("ele");
   _s_flavour.push_back("mu");
   
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < num_of_processes; i_proc++)
      {
         _histo_name = "Passing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         passing[i_proc][i_flav] = new TH2F(_histo_name.c_str(),"", 200, 0, 200, 500, -2.5, 2.5);
         
         _histo_name = "Failing_" + _s_process.at(i_proc) + "_" + _s_flavour.at(i_flav);
         failing[i_proc][i_flav] = new TH2F(_histo_name.c_str(),"", 200, 0, 200, 500, -2.5, 2.5);

      }
      
      _histo_name = "Passing_WZsubtracted_" + _s_flavour.at(i_flav);
      passing_WZ_sub[i_flav] = new TH2F(_histo_name.c_str(),"", 200, 0, 200, 500, -2.5, 2.5);
      _histo_name = "Failing_WZsubtracted_" + _s_flavour.at(i_flav);
      failing_WZ_sub[i_flav] = new TH2F(_histo_name.c_str(),"", 200, 0, 200, 500, -2.5, 2.5);

   }
   
   
}
//===============================================================

//===============================================================
void OSmethod::SaveHistos()
{
   TFile* fOutHistos = new TFile("Histos.root", "recreate");
   fOutHistos->cd();
   
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      for (int i_proc = 0; i_proc < num_of_processes; i_proc++)
      {
         passing[i_proc][i_flav]->Write();
         failing[i_proc][i_flav]->Write();
      }
      
   passing_WZ_sub[i_flav]->Write();
   failing_WZ_sub[i_flav]->Write();
      
   }
   
   fOutHistos->Close();
   delete fOutHistos;
}
//===============================================================

//===============================================================
void OSmethod::SubtractWZ()
{
   for (int i_flav = 0; i_flav < num_of_flavours; i_flav++)
   {
      passing_WZ_sub[i_flav]->Add(passing[Settings::Data][i_flav], 1.);
      passing_WZ_sub[i_flav]->Add(passing[Settings::WZ][i_flav], -1.);
      
      failing_WZ_sub[i_flav]->Add(failing[Settings::Data][i_flav], 1.);
      failing_WZ_sub[i_flav]->Add(failing[Settings::WZ][i_flav], -1.);
   }

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



