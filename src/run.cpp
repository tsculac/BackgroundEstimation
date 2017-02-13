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

   OSmethod *os = new OSmethod();

   os->FillHist(Data);
   os->FillHist(WZ);
   os->SaveHistos();
   
   delete os;
}
