#define MuonTreeAnalysis_cxx
#include "MuonTreeAnalysis.h"
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

void MuonTreeAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MuonTreeAnalysis.C
//      root> MuonTreeAnalysis t
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



   TString outputname = "MuonTreeAnalysis.root";
   TFile *outputfile;
   outputfile = TFile::Open( outputname, "RECREATE" );





   Double_t pi = 3.14159265359;
   // cut flow stats
   Long_t n_cut1 = 0, n_cut2 = 0;
   vector<Long_t> n_trkpersec; // number of matched trks per sector, eff denominator

   vector<Long_t> n_trkpersec_phi_good; // number of matched trks per sector when we take out all the problematic areas 
									// problematic areas meaning stuck-bit channels here  

   vector<Long_t> n_trkpersec_eta_good;	//taking out problematic areas which are stuck-bit channels and eta channels 0-49, 181-192 

   vector<vector<vector<Long_t>>> n_trkperchl; // number of matched trks per channel, eff denominator


   n_trkpersec.resize(33, 0); // magic number 33=16+16+1, 16 sectors on each end, the extra 1 is for convenience filling histogram
   n_trkpersec_phi_good.resize(33, 0); 
   n_trkpersec_eta_good.resize(33, 0); 

   n_trkperchl.resize(33);

   for(UInt_t i=0; i<n_trkperchl.size(); ++i) {
      n_trkperchl[i].resize(4); // 4 layers per sector
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
//		 n_trkperchl[i][j].resize(241, 0); // magic number 242=48+192+1, 48 phi channels, 192 eta channels, 1 for convenience
		 n_trkperchl[i][j].resize(49, 0); // magic number 242=48+192+1, 48 phi channels, 192 eta channels, 1 for convenience
		}
   }




//   TH1F * h_hit_channel_phi  = new TH1F("hit_channel_phi", "hit_channel_phi", 50, 0, 49);
//   TH1F * h_hit_channel_phi_oneLayer  = new TH1F("hit_channel_phi_oneLayer", "hit_channel_phi_oneLayer", 50, 0, 49);
   TH1F * h_hit_channel_phi_oneLayer  = new TH1F("hit_channel_phi_oneLayer", "hit_channel_phi_oneLayer", 50, -49, 0);

//   TH2F * h2_hit_channel_phi_byLayers  = new TH2F("hit_channel_phi_byLayers", "hit_channel_phi_byLayers", 50, 0, 49, 170, -17, 17);

//   TH2F * h2_eff_channel  = new TH2F("eff_channel", "eff_channel", 49, 0, 48, 170, -17, 17);
//   TH2F * h2_eff_channel  = new TH2F("eff_channel", "eff_channel", 242, -49, 193, 170, -17, 17);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      if (jentry%100000==0) cout<<"Event: "<<jentry<<endl;

//_______________________________________________________________________________________________________________

      // CUT1: at least one CSC trk that has P>60 GeV, i.e. good trk
      if(trkPt->size()==0) continue; // no trk, pass
      vector<UInt_t>   goodtrk_index;
      vector<Double_t> goodtrk_phi;
      for(UInt_t i=0; i<trkP->size(); ++i) { // loop over trks
         if( TMath::Abs(trkEta->at(i))<2. || TMath::Abs(trkEta->at(i))>2.7 ) continue; // cut out non-CSC region
         if( (trkP->at(i)/1.e3)<=60 ) continue; // cut on pt of trks
         goodtrk_index.push_back(i); // if no good trk in this event, it will not be filled. so no worry
         Double_t phi = 0.;
         phi = TMath::ATan2(trkPy->at(i), trkPx->at(i));
         goodtrk_phi.push_back(phi);
      } // end loop of trks

      if(goodtrk_index.size()==0) continue; // no good trk, pass
      ++n_cut1;

//_____________________________________________________________________________________________________________________

      // CUT2: at least one matched good trk
      vector<UInt_t>   matchedtrk_index;
      vector<Double_t> matchedtrk_phi;
      vector<UInt_t>   morethan2trk_index; // tracks that are matched by more than 2 sectors, should be good ones but dunno how to remove ill matched sectors
      for(UInt_t i=0; i<goodtrk_index.size(); ++i) { // loop over good trks
         Int_t sector = 0;
         Int_t sector_2ndary = 0;

         bool morethan2 = false; // if matched by more than 2 sectors
         for(UInt_t j=0; j<hit_sector->size(); ++j) {  // loop over hits

            if( hitToMuon->at(j)==goodtrk_index[i] ) {
 
               Double_t badphi = 0.; // to exclude mismatching
               if( hit_sector->at(j)<0 ) badphi = -(pi/15)*(2*hit_sector->at(j)+17);
               else badphi = (pi/15)*(2*hit_sector->at(j)-17);
               if( TMath::Abs(badphi-goodtrk_phi[i]) > 0.7 && TMath::Abs(badphi+2*pi-goodtrk_phi[i]) > 0.7 ) {
                  sector = hit_sector->at(j);
                  break;
               }
            }
         } // end loop of hits


         for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits again to find overlaped sector
            if( hitToMuon->at(j)==goodtrk_index[i] && hit_sector->at(j)!=sector && sector!=0 ) {
               Double_t badphi = 0.;
               if( hit_sector->at(j)<0 ) badphi = -(pi/15)*(2*hit_sector->at(j)+17);
               else badphi = (pi/15)*(2*hit_sector->at(j)-17);
               if( TMath::Abs(badphi-goodtrk_phi[i]) > 0.7 && TMath::Abs(badphi+2*pi-goodtrk_phi[i]) > 0.7 ) {
                  sector_2ndary = hit_sector->at(j);
                  break;
               }
            }
         } // end loop of hits
         for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits once more
            if( hitToMuon->at(j)==goodtrk_index[i] && hit_sector->at(j)!=sector && hit_sector->at(j)!=sector_2ndary && sector!=0 && sector_2ndary!=0 ) {
               Double_t badphi = 0.;
               if( hit_sector->at(j)<0 ) badphi = -(pi/15)*(2*hit_sector->at(j)+17);
               else badphi = (pi/15)*(2*hit_sector->at(j)-17);
               if( TMath::Abs(badphi-goodtrk_phi[i]) > 0.7 && TMath::Abs(badphi+2*pi-goodtrk_phi[i]) > 0.7 ) {
                  morethan2trk_index.push_back(goodtrk_index[i]);
                  morethan2 = true;
//                  cout<<"ERROR: more than two matched sectors!"<<endl;
 //                cout<<"info: "<<goodtrk_index[i]<<"  "<<trkP->at(goodtrk_index[i])/1.e3<<"  "<<trkEta->at(goodtrk_index[i])<<"  "<<goodtrk_phi[i]<<endl;
 //                 cout<<"      "<<sector<<"  "<<sector_2ndary<<"  "<<hit_sector->at(j)<<endl;
                  break;
               }
            }
         } // end loop of hits


         if( morethan2 ) continue; // abandon morethan2 trk

         if( sector!=0 && sector_2ndary!=0 && ( TMath::Abs(sector-sector_2ndary)==2 || TMath::Abs(sector-sector_2ndary)==14 ) ) continue; // abandon those have ill overlap

         if( sector!=0 ) { // if this good trk is matched with hits
            matchedtrk_index.push_back(goodtrk_index[i]);
            matchedtrk_phi.push_back(goodtrk_phi[i]);


            ++n_trkpersec[sector+16]; // why 16? see the last part of this code for conversion rule; ex if sec -16, then it would be 0 here
 
	        ++n_trkpersec_phi_good[sector+16];
	        ++n_trkpersec_eta_good[sector+16];

            Int_t avgchlphi = 0, avgchleta = 0; // average phi channel, average eta channel
            Int_t n_chlphi = 0, n_chleta = 0; // number of phi hits and eta hits, used to calculate average
            vector<Int_t> layer_phi; layer_phi.resize(4, 0); // number of hits on each layer
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);


            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector ) continue;

	
//	           ++n_trkperchl[sector+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48]; //n_trkperchl[sector#][layer#][channel#], conversion rule down below 


				if( hit_measphi->at(j)==1 ) {
	           		++n_trkperchl[sector+16][int(hit_wlay->at(j))-1][int(hit_pstrip->at(j))]; 

				
//               if( hit_measphi->at(j)==1 ) {
               		avgchlphi += hit_pstrip->at(j);
                 	++n_chlphi;
                  	++layer_phi[int(hit_wlay->at(j))-1]; 

               } 
				
/*				else{
		              avgchleta += hit_pstrip->at(j);
		              ++n_chleta;
		              ++layer_eta[int(hit_wlay->at(j))-1];
               }

*/

            } // end loop of hits

		
            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi); 
//            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
			 

            for(UInt_t j=0; j<layer_phi.size(); ++j) { // loop over all phi layers, accumulate the denominator


               if( layer_phi[j]==0 ) 
				{	
//					++n_trkperchl[sector+16][j][-1*avgchlphi+48];
					++n_trkperchl[sector+16][j][avgchlphi];

				}

			}
				
/*            for(UInt_t j=0; j<layer_eta.size(); ++j) { // loop over all eta layers, accumulate the denominator

               if( layer_eta[j]==0 ) 
				{	
					++n_trkperchl[sector+16][j][avgchleta+48];

				}
					
         	}
*/		

		} 

		 else { // if not matched

         }


         if( sector_2ndary!=0 ) { // if has overlap, this good trk is counted twice

            ++n_trkpersec[sector_2ndary+16];
	    	++n_trkpersec_phi_good[sector_2ndary+16];
	    	++n_trkpersec_eta_good[sector_2ndary+16];

            Int_t avgchlphi = 0, avgchleta = 0;
            Int_t n_chlphi = 0, n_chleta = 0;
            vector<Int_t> layer_phi; layer_phi.resize(4, 0);
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);


            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector_2ndary ) continue;


//				++n_trkperchl[sector_2ndary+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
				
               if( hit_measphi->at(j)==1 ) {
				  ++n_trkperchl[sector_2ndary+16][int(hit_wlay->at(j))-1][int(hit_pstrip->at(j))];
                  avgchlphi += hit_pstrip->at(j);
                  ++n_chlphi;
                  ++layer_phi[int(hit_wlay->at(j))-1];

               } 
/*				else {
		              avgchleta += hit_pstrip->at(j);
		              ++n_chleta;
		              ++layer_eta[int(hit_wlay->at(j))-1];

               }
*/

            }


            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
//            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);

            for(UInt_t j=0; j<layer_phi.size(); ++j) {
               if( layer_phi[j]==0 ) 
//					++n_trkperchl[sector_2ndary+16][j][-1*avgchlphi+48];
					++n_trkperchl[sector_2ndary+16][j][avgchlphi];

			}

/*            for(UInt_t j=0; j<layer_eta.size(); ++j) {
               if( layer_eta[j]==0 ) 
					++n_trkperchl[sector_2ndary+16][j][avgchleta+48];

         	}
*/

		} 
		 

      } // end loop of good trks

      if(matchedtrk_index.size()==0) continue; // no matched trk, pass
      ++n_cut2;


//________________________________________________________________________________________________________________________


      for(UInt_t i=0; i<hit_sector->size(); ++i) { // loop over hits

         bool matched = false;
         Double_t badphi = 0., goodphi = 0.; // to get rid of mismatched hits

         for(UInt_t j=0; j<matchedtrk_index.size(); ++j) {// loop over matched trks
            if( hitToMuon->at(i)==matchedtrk_index[j] ) {
               matched = true;
               goodphi = matchedtrk_phi[j];

            }
			
         } // end loop of matched trks


         if( !matched ) continue;
         if( hit_sector->at(i)<0 ) badphi = -(pi/15)*(2*hit_sector->at(i)+17);
         else badphi = (pi/15)*(2*hit_sector->at(i)-17);
         if( TMath::Abs(badphi-goodphi) <= 0.7 || TMath::Abs(badphi+2*pi-goodphi) <= 0.7) matched = false;
         if( !matched ) continue;

						
         Float_t fillnumber = hit_sector->at(i)+0.25*(hit_wlay->at(i)-1); // a combination of sector # and layer #

         if( hit_measphi->at(i)==1 ) 
		  {
 
            h2_eff_channel->Fill(hit_pstrip->at(i), fillnumber, 1);

			if(fillnumber== 3.25)
			{

				h_hit_channel_phi_oneLayer->Fill(-1.0*hit_pstrip->at(i));					
			}

         } 

/*		else {

            h2_eff_channel->Fill(hit_pstrip->at(i), fillnumber, 1);
       
	    }
*/

     } // end loop of hits




   }
//end loop of events 

//_______________________________________________________________________________________________________



   for(UInt_t i=0; i<n_trkperchl.size(); ++i) { // normalization
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
         for(UInt_t k=0; k<n_trkperchl[i][j].size(); ++k) {



            if( h2_eff_channel->GetBinContent( k, (i+1)*5+j+1 )!=0 ) {
               Double_t bincontent = h2_eff_channel->GetBinContent( k, (i+1)*5+j+1 );
//			   cout << k << "	" << (i+1)*5+j+1 <<"	" << bincontent << "	" << n_trkperchl[i][j][k] << endl; 
               Double_t binerror = bincontent*(n_trkperchl[i][j][k]-bincontent)/TMath::Power(n_trkperchl[i][j][k], 3);
               h2_eff_channel->SetBinContent( k, (i+1)*5+j+1, bincontent/n_trkperchl[i][j][k] );
               h2_eff_channel->SetBinError( k, (i+1)*5+j+1, binerror );

	     	}
			

         }
      }
   }



   outputfile->Write();
   outputfile->Close();


}
