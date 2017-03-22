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
   
   TString path = "Moriond_2017/";
   TString file_name = "/ZZ4lAnalysis.root";
   TString file_name_FR = "/FakeRate_SS_Moriond368.root";
   
   TString Data    = path + "Data"       + file_name;
   TString WZ      = path + "WZTo3LNu"   + file_name;
   TString ZZ      = path + "ZZTo4l"     + file_name;
   TString ttbar   = path + "TTJets"     + file_name;
   TString DY      = path + "DYJetsToLL" + file_name;
   
   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};

   OSmethod *os = new OSmethod();

   os->SetLumi(35.9);
   
   ///////////////////////////////////
   // Fill control histos           //
   ///////////////////////////////////
   os->FillDataMCPlots(Data);
   os->FillDataMCPlots(WZ);
   os->FillDataMCPlots(ZZ);
   os->FillDataMCPlots(ttbar);
   os->FillDataMCPlots(DY);
   os->SaveDataMCHistos("DataMC.root");

   ///////////////////////////////////
   // Fill passing/failling histos  //
   ///////////////////////////////////
   os->FillFRHistos(Data);
   os->FillFRHistos(WZ);
   os->SubtractWZ();
   os->SaveFRHistos("Histos.root", true);
   
   ///////////////////////////////////
   // Calculate fake rates          //
   ///////////////////////////////////
   os->GetFRHistos("Histos.root");
   os->Set_pT_binning(8, pT_bins);
   os->ProduceFakeRates("FakeRates_OS_Moriond17.root");

   ///////////////////////////////////
   // Fill ZX contributions histos  //
   ///////////////////////////////////
   os->MakeHistogramsZX(Data, "FakeRates_OS_Moriond17.root");
   os->MakeZXMCContribution(ZZ, "FakeRates_OS_Moriond17.root");
   os->SaveZXHistos("ZXHistos.root");

   ///////////////////////////////////
   // Plot control plots            //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos.root");
   os->GetDataMCHistos("DataMC.root");
   os->PlotDataMC_2P2F( "M4l", "Plots" );
   os->PlotDataMC_3P1F( "M4l", "Plots" );
   
   ///////////////////////////////////
   // Plot Z+X plots                //
   ///////////////////////////////////
   os->GetZXHistos("ZXHistos.root");
   os->PlotZXContributions("Plots");
   
   delete os;
}
