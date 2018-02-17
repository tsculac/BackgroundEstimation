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
#include "OSmethod.h"

using namespace std;

int main( int argc, char *argv[] )
{
   gROOT->ProcessLine(".L ./ext/setTDRStyle_cpp.so");
   gROOT->ProcessLine("setTDRStyle();");
   
   TString path = "/eos/user/t/tsculac/BigStuff/";
   TString file_name = "/ZZ4lAnalysis.root";
   
   TString Data    = path + "2017Production_Data_FinalJSON/AllData"       + file_name;
   TString WZ      = "Moriond_2017/WZTo3LNu"   + file_name;
   TString ZZ      = path + "2017Production_MC/ZZTo4l"     + file_name;
   TString ttbar   = "Moriond_2017/TTJets"     + file_name;
   TString DY      = path + "2017Production_MC/DYJetsToLL_M50" + file_name;
	
   bool SubtractWZ = true;
   bool Remove_NegBins_FR = true;
   bool Remove_NegBins_ZX = true;
	
   
   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};
   
   OSmethod *os = new OSmethod();

   os->SetLumi(41.37);

   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   os->FillDataMCPlots(Data);
   os->FillDataMCPlots(WZ);
   os->FillDataMCPlots(ZZ);
   os->FillDataMCPlots(ttbar);
   os->FillDataMCPlots(DY);
   os->SaveDataMCHistos("DataMC_OS_Moriond18.root");

   ///////////////////////////////////
   // Fill passing/failling histos  //
   ///////////////////////////////////
   os->FillFRHistos(Data);
   os->FillFRHistos(WZ);
   os->SaveFRHistos("Histos_OS_Moriond18.root", SubtractWZ, Remove_NegBins_FR);

   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   os->GetFRHistos("Histos_OS_Moriond18.root");
   os->Set_pT_binning(8, pT_bins);
   os->ProduceFakeRates("FakeRates_OS_Moriond18.root");

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   os->MakeHistogramsZX(Data, "FakeRates_OS_Moriond18.root");
   os->MakeZXMCContribution(ZZ, "FakeRates_OS_Moriond18.root");
   os->SaveZXHistos("ZXHistos_OS_Moriond18.root", Remove_NegBins_ZX);

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos_OS_Moriond18.root");
   os->GetDataMCHistos("DataMC_OS_Moriond18.root");
   os->PlotDataMC("M4l", "Plots");
   os->PlotDataMC_2P2F( "M4l", "Plots" );
   os->PlotDataMC_3P1F( "M4l", "Plots" );

   ///////////////////////////////////
   // Plot Z+X plots                //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos_OS_Moriond18.root");
   os->PlotZXContributions("Plots");
   os->FitZX("Plots");
	
   delete os;
}
