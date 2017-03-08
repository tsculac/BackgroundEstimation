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
   TString path = "Moriond_2017/";
   TString file_name = "/ZZ4lAnalysis.root";
   TString file_name_FR = "/FakeRate_SS_Moriond368.root";
   
   TString Data    = path + "Data" + file_name;
   TString WZ      = path + "WZTo3LNu" + file_name;
   
   float pT_bins[] = {5, 7, 10, 20, 30, 40, 50, 80};

   OSmethod *os = new OSmethod();

   os->SetLumi(36.8);
   
   os->FillHistos(Data);
   os->FillHistos(WZ);
   os->SubtractWZ(true);
   os->SaveHistos("Histos.root");
   
   os->GetHistos("Histos.root");
   os->Set_pT_binning(8, pT_bins);
   os->ProduceFakeRates("FakeRates_OS_Moriond17.root");
   
   os->MakeHistogramsZX(Data, "FakeRates_OS_Moriond17.root");
   
   delete os;
}
