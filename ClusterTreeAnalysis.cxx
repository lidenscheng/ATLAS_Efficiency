#define ClusterTreeAnalysis_cxx
#include "ClusterTreeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TImage.h> 
#include <iostream>
#include <vector>
#include <fstream> 

using namespace std;

void ClusterTreeAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ClusterTreeAnalysis.C
//      root> ClusterTreeAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   TString outputname = "ClusterTreeAnalysis.root";
   TFile *outputfile;
   outputfile = TFile::Open( outputname, "RECREATE" );




//   TH1F * h_hit_channel_phi  = new TH1F("hit_channel_phi", "hit_channel_phi", 50, 0, 49);
   TH1F * h_hit_channel_phi_oneLayer  = new TH1F("hit_channel_phi_oneLayer", "hit_channel_phi_oneLayer", 50, -49, 0);
//   TH2F * h2_hit_channel_phi_byLayers  = new TH2F("hit_channel_phi_byLayers", "hit_channel_phi_byLayers", 50, 0, 49, 170, -17, 17);




   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (jentry%100000==0) cout<<"Event: "<<jentry<<endl;

      // if (Cut(ientry) < 0) continue;



      for(UInt_t i=0; i<sector->size(); ++i) 
	  { // loop over hits

            Float_t fillnumber = sector->at(i)+0.25*(wlay->at(i)-1); // a combination of sector # and layer #

/* 			if( measphi->at(i)==1 ) 
			{
				h_hit_channel_phi->Fill(pstrip->at(i));

			}
*/

 			if( measphi->at(i)==1 ) 
			{
				if(fillnumber == 3.25)
					h_hit_channel_phi_oneLayer->Fill(-1.0*pstrip->at(i));					

			}


/* 			if( measphi->at(i)==1 ) 
			{
				h2_hit_channel_phi_byLayers->Fill(pstrip->at(i), fillnumber, 1);

			}
*/

	  }

   }



   outputfile->Write();
   outputfile->Close();

}















