#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include "MuonTreeAnalysis.h"

using namespace std;

int main( int argc, char** argv ) {

   TString inputname = "/eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root";

   TFile *inputfile = TFile::Open( inputname );
   TTree *tree = 0;
   inputfile->GetObject( "muon", tree );

   MuonTreeAnalysis muontreeanalysis(tree);
   muontreeanalysis.Loop();

   inputfile->Close();
}
