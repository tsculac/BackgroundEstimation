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
   
   passing[0] = new TH2F("passing_Data","passing_Data", 1000, 0, 1000, 500, -2.5, 2.5);
   failing[0] = new TH2F("failing_Data","failing_Data", 1000, 0, 1000, 500, -2.5, 2.5);
   passing[1] = new TH2F("passing_WZ","passing_WZ", 1000, 0, 1000, 500, -2.5, 2.5);
   failing[1] = new TH2F("failing_WZ","failing_WZ", 1000, 0, 1000, 500, -2.5, 2.5);
   
}
//============================================================



// Destructor
//====================
OSmethod::~OSmethod()
{
}
//====================


//===============================================================================
void OSmethod::FillHist( TString input_file_data_name )
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

         if(LepisID->at(2) && LepCombRelIsoPF->at(2) < 0.35) passing[_current_process]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
         else failing[_current_process]->Fill(LepPt->at(2), LepEta->at(2), (_current_process == Settings::Data) ? 1 :  _event_weight);
      }
   } // END events loop

}
//===============================================================================

void OSmethod::SaveHistos()
{
   TFile* fOutHistos = new TFile("TEST.root", "recreate");
   fOutHistos->cd();
   passing[0]->Write();
   failing[0]->Write();
   passing[1]->Write();
   failing[1]->Write();
   fOutHistos->Close();
   delete fOutHistos;
}

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



