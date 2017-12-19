// C++
#include <iostream>
#include <fstream>
#include <string>

// ROOT
#include "TApplication.h"
#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

// My own files
#include "SSmethod.h"

using namespace std;

int main( int argc, char *argv[] )
{
   gROOT->ProcessLine(".L ./ext/setTDRStyle_cpp.so");
   gROOT->ProcessLine("setTDRStyle();");
   
   TString path = "Moriond_2017/";
   TString file_name = "/ZZ4lAnalysis.root";
   TString file_name_FR = "/FakeRate_SS_Moriond368.root";
   
   TString Data    = path + "Data"       + file_name;
   TString WZ      = path + "WZTo3LNu"   + file_name;
   TString ZZ      = path + "ZZTo4l"     + file_name;
   TString ttbar   = path + "TTJets"     + file_name;
   TString DY      = path + "DYJetsToLL" + file_name;
   TString FakeRates   = path + "FakeRates" + file_name_FR;
	
   bool SubtractWZ = true;
   bool Remove_NegBins_FR = true;
   bool SubtractMCContribution = true;
	
   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};

   SSmethod *ss = new SSmethod();
   ss->SetLumi(35.9);

   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   ss->FillDataMCPlots(Data);
   ss->FillDataMCPlots(WZ);
   ss->FillDataMCPlots(ZZ);
   ss->FillDataMCPlots(ttbar);
   ss->FillDataMCPlots(DY);
   ss->SaveDataMCHistos("DataMC_SS.root");

   ///////////////////////////////////
   // Fill passing/failling histos  //
   ///////////////////////////////////
   ss->FillFRHistos(Data);
   ss->FillFRHistos(WZ);
   ss->SaveFRHistos("Histos_SS.root", SubtractWZ, Remove_NegBins_FR);

   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   ss->GetFRHistos("Histos_SS.root");
   ss->Set_pT_binning(8, pT_bins);
   ss->ProduceFakeRates("FakeRates_SS_Moriond17.root", Data);

   ///////////////////////////////////
   // Calculate OS/SS ratios        //
   ///////////////////////////////////
   ss->Calculate_SSOS_Ratio( Data, ZZ, SubtractMCContribution);

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   ss->MakeHistogramsZX(Data, FakeRates);
   ss->SaveZXHistos("ZXHistos_SS.root");

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   ss->GetDataMCHistos("DataMC_SS.root");
   ss->PlotDataMC( "M4l", "Plots" );

   ///////////////////////////////////
   // Plot Z+X plots                //
   ///////////////////////////////////
   ss->GetZXHistos("ZXHistos_SS.root");
   ss->PlotZX("M4l", "Plots");
	
   delete ss;
}
