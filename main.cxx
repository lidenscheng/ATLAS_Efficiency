#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include "EffStudy.h"

using namespace std;

int main( int argc, char** argv ) {
   TString inputname = "/export/share/data2/you/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root";
   TFile *inputfile = TFile::Open( inputname );
   TTree *tree = 0;
   inputfile->GetObject( "muon", tree );

   EffStudy effstudy(tree);
   effstudy.Loop();

   inputfile->Close();
}
