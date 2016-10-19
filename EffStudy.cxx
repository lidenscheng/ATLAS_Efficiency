#define EffStudy_cxx
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <fstream> 
#include "EffStudy.h"

using namespace std;

void EffStudy::Loop()
{
   TString outputname = "EffStudy.root";
   TFile *outputfile;
   outputfile = TFile::Open( outputname, "RECREATE" );

//   ifstream inputPhi ("exclude_phi.dat", ifstream::in);
//   ifstream inputEta ("exclude_eta.dat", ifstream::in);
   ifstream input("exclude.dat", ifstream::in);

   vector<UInt_t> secVector, layerVector, channelVector;
//   vector<UInt_t> layerPhi, layerEta;
//   vector<UInt_t> channelPhi, channelEta; 
   UInt_t sec, channel, layer; 

   while(input >> sec >> layer >> channel)
   	{

		secVector.push_back(sec);
		layerVector.push_back(layer);
		channelVector.push_back(channel);

   }

/*   while(inputEta >> sec >> layer >> channel)
   {
		secEta.push_back(sec);
		layerEta.push_back(layer);
		channelEta.push_back(channel);
   }
*/

   
   Double_t pi = 3.14159265359;
   // cut flow stats
   Long_t n_cut1 = 0, n_cut2 = 0;
   vector<Long_t> n_trkpersec; // number of matched trks per sector, eff denominator
   vector<Long_t> n_trkpersec_noStuck; // number of matched trks per sector when we take out stuck bits, eff denominator

   vector<vector<vector<Long_t>>> n_trkperchl; // number of matched trks per channel, eff denominator
   vector<vector<vector<Long_t>>> n_trkperchl_modified; // used to count for the case of ignoring ch 1-48 

   n_trkpersec.resize(33, 0); // magic number 33=16+16+1, 16 sectors on each end, the extra 1 is for convenience filling histogram
   n_trkpersec_noStuck.resize(33, 0); 

   n_trkperchl.resize(33);
   n_trkperchl_modified.resize(33);

   for(UInt_t i=0; i<n_trkperchl.size(); ++i) {
      n_trkperchl[i].resize(4); // 4 layers per sector
      n_trkperchl_modified[i].resize(4); 
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
		 n_trkperchl[i][j].resize(241, 0); // magic number 242=48+192+1, 48 phi channels, 192 eta channels, 1 for convenience
		 n_trkperchl_modified[i][j].resize(241, 0);
		}
   }

   // book histograms
   TH1F * h_cutflow      = new TH1F("cutflow", "cutflow", 3, 0, 3);
//   TH1F * h_trk_eta_org  = new TH1F("trk_eta_org", "trk_eta_org", 66, -3.3, 3.3); // before any cut
   TH1F * h_trk_eta      = new TH1F("trk_eta", "trk_eta", 66, -3.3, 3.3);
   TH1F * h_trk_P        = new TH1F("trk_P", "trk_P", 100, 0, 600);
   TH1F * h_trk_phi      = new TH1F("trk_phi", "trk_phi", 500, -3.3, 3.3);
   TH1F * h_trk_sec      = new TH1F("trk_sec", "trk_sec", 34, -17, 17);
//   TH1F * h_hit_channel_org = new TH1F("hit_channel_org", "hit_channel_org", 242, -49, 193); // before any cut
   TH1F * h_hit_channel  = new TH1F("hit_channel", "hit_channel", 242, -49, 193);
   TH2F * h2_trk_sec_phi = new TH2F("trk_sec_phi", "trk_sec_phi", 34, -17, 17, 1000, -3.3, 3.3); // use large number of bins
   TH2F * h2_trk_P_eta_0 = new TH2F("trk_P_eta_0", "trk_P_eta_0", 100, 0, 900, 70, 2., 2.7); // _0 sets are good matched tracks
   TH2F * h2_trk_P_phi_0 = new TH2F("trk_P_phi_0", "trk_P_phi_0", 100, 0, 900, 66, -3.3, 3.3);
   TH2F * h2_trk_eta_phi_0 = new TH2F("trk_eta_phi_0", "trk_eta_phi_0", 70, 2., 2.7, 66, -3.3, 3.3);
   TH2F * h2_trk_P_eta_1 = new TH2F("trk_P_eta_1", "trk_P_eta_1", 100, 0, 900, 70, 2., 2.7); // _1 sets are unmatched tracks
   TH2F * h2_trk_P_phi_1 = new TH2F("trk_P_phi_1", "trk_P_phi_1", 100, 0, 900, 66, -3.3, 3.3);
   TH2F * h2_trk_eta_phi_1 = new TH2F("trk_eta_phi_1", "trk_eta_phi_1", 70, 2., 2.7, 66, -3.3, 3.3);
   TH1F * h_overlap      = new TH1F("overlap", "overlap", 18, 0, 18); // sector number difference between overlaped sectors
   TH1F * h_eff_phi      = new TH1F("eff_phi", "eff_phi", 170, -17, 17);
   TH1F * h_eff_eta      = new TH1F("eff_eta", "eff_eta", 170, -17, 17);

   TH2F * h2_eff_channel  = new TH2F("eff_channel", "eff_channel", 242, -49, 193, 170, -17, 17);
//   TH2F * h2_eff_channel_modified  = new TH2F("eff_channel_modified", "eff_channel_modified", 242, -49, 193, 170, -17, 17); //take out ch 1-48

   TH1F * h_eff_phi_noStuck      = new TH1F("eff_phi_noStuck", "eff_phi_noStuck", 170, -17, 17);
   TH1F * h_eff_eta_noStuck      = new TH1F("eff_eta_noStuck", "eff_eta_noStuck", 170, -17, 17);

   TH2F * h2_eff_channel_onlyNumerator = new TH2F("eff_channel_numerator", "eff_channel_numerator", 242, -49, 193, 170, -17, 17);

//   TObjArray *HList = new TObjArray(0); 
//vector<TH1F*> HList; 
   TList* HList = new TList(); 

   h_trk_phi->Sumw2();
   h_trk_sec->Sumw2();

   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
//   cout << "nentries = " << nentries << endl; 
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

//      if (jentry>10000) break;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//      if (jentry%100000==0) cout<<"Event: "<<jentry<<endl;
      h_cutflow->Fill(0);

/*
      for(UInt_t i=0; i<hit_sector->size(); ++i) { // loop over hits
         if( hit_measphi->at(i)==1 ) {
            h_hit_channel_org->Fill( -1.*hit_pstrip->at(i) );
         } else {
            h_hit_channel_org->Fill( hit_pstrip->at(i) );
         }
      }
*/

      // ===============================================================================================
      // CUT1: at least one CSC trk that has P>60 GeV, i.e. good trk
      if(trkPt->size()==0) continue; // no trk, pass
      vector<UInt_t>   goodtrk_index;
      vector<Double_t> goodtrk_phi;
      for(UInt_t i=0; i<trkP->size(); ++i) { // loop over trks
//         h_trk_eta_org->Fill(trkEta->at(i));  
         if( TMath::Abs(trkEta->at(i))<2. || TMath::Abs(trkEta->at(i))>2.7 ) continue; // cut out non-CSC region
         if( (trkP->at(i)/1.e3)<=60 ) continue; // cut on pt of trks
         goodtrk_index.push_back(i); // if no good trk in this event, it will not be filled. so no worry
         h_trk_P->Fill(trkP->at(i)/1.e3);
         h_trk_eta->Fill(trkEta->at(i));
         Double_t phi = 0.;
         phi = TMath::ATan2(trkPy->at(i), trkPx->at(i));
         goodtrk_phi.push_back(phi);
         h_trk_phi->Fill(phi);
      } // end loop of trks
      if(goodtrk_index.size()==0) continue; // no good trk, pass
//	  cout << "cut1: " << goodtrk_index.size() << endl; 
      ++n_cut1;
      h_cutflow->Fill(1);
      // ===============================================================================================

      // ===============================================================================================
      // CUT2: at least one matched good trk
      vector<UInt_t>   matchedtrk_index;
      vector<Double_t> matchedtrk_phi;
      vector<UInt_t>   morethan2trk_index; // tracks that are matched by more than 2 sectors, should be good ones but dunno how to remove ill matched sectors
      for(UInt_t i=0; i<goodtrk_index.size(); ++i) { // loop over good trks
//		cout << "size of pstrip vector: " << hit_pstrip->size() << endl; 
         Int_t sector = 0;
         Int_t sector_2ndary = 0;
         bool morethan2 = false; // if matched by more than 2 sectors
         for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
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
//                 cout<<"info: "<<goodtrk_index[i]<<"  "<<trkP->at(goodtrk_index[i])/1.e3<<"  "<<trkEta->at(goodtrk_index[i])<<"  "<<goodtrk_phi[i]<<endl;
//                  cout<<"      "<<sector<<"  "<<sector_2ndary<<"  "<<hit_sector->at(j)<<endl;
                  break;
               }
            }
         } // end loop of hits

         if( morethan2 ) continue; // abandon morethan2 trk

         if( sector!=0 && sector_2ndary!=0 && ( TMath::Abs(sector-sector_2ndary)==2 || TMath::Abs(sector-sector_2ndary)==14 ) ) continue; // abandon those have ill overlap

         if( sector!=0 ) { // if this good trk is matched with hits
//			cout << "sector: " << sector << endl; 
            matchedtrk_index.push_back(goodtrk_index[i]);
            matchedtrk_phi.push_back(goodtrk_phi[i]);
            h_trk_sec->Fill(sector);
            h2_trk_sec_phi->Fill(sector, goodtrk_phi[i]);
            h2_trk_P_eta_0->Fill(trkP->at(goodtrk_index[i])/1.e3, TMath::Abs(trkEta->at(goodtrk_index[i])));
            h2_trk_P_phi_0->Fill(trkP->at(goodtrk_index[i])/1.e3, goodtrk_phi[i]);
            h2_trk_eta_phi_0->Fill(TMath::Abs(trkEta->at(goodtrk_index[i])), goodtrk_phi[i]);

            ++n_trkpersec[sector+16]; // why 16? see the last part of this code for conversion rule; ex if sec -16, then it would be 0 here 
            ++n_trkpersec_noStuck[sector+16];  
            Int_t avgchlphi = 0, avgchleta = 0; // average phi channel, average eta channel
//			Int_t avgchleta_modified = 0; 
            Int_t n_chlphi = 0, n_chleta = 0; // number of phi hits and eta hits, used to calculate average
//			Int_t n_chleta_modified = 0; 
            vector<Int_t> layer_phi; layer_phi.resize(4, 0); // number of hits on each layer
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector ) continue;
		
//			   if(hit_pstrip->at(j)>0 && hit_pstrip->at(j)<=48){					
	               ++n_trkperchl[sector+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48]; //n_trkperchl[sector#][layer#][channel#], conversion rule down below
//				}
				
/*				else{
					++n_trkperchl[sector+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
					++n_trkperchl_modified[sector+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
				}
*/

               if( hit_measphi->at(j)==1 ) {
                  avgchlphi += hit_pstrip->at(j);
                  ++n_chlphi;
                  ++layer_phi[int(hit_wlay->at(j))-1];

               } else {
//						if(hit_pstrip->at(j)>0 && hit_pstrip->at(j)<=48){
		              		avgchleta += hit_pstrip->at(j);
		              		++n_chleta;
		              		++layer_eta[int(hit_wlay->at(j))-1];
//							}

/*						else{
							avgchleta += hit_pstrip->at(j);
							avgchleta_modified += hit_pstrip->at(j);
		              		++n_chleta;
		              		++n_chleta_modified;
		              		++layer_eta[int(hit_wlay->at(j))-1];
		              		++layer_eta_modified[int(hit_wlay->at(j))-1];
							}
*/

               }
            } // end loop of hits

		
            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
//            if( avgchleta_modified!=0) avgchleta_modified = int(avgchleta_modified/n_chleta_modified);

            for(UInt_t j=0; j<layer_phi.size(); ++j) { // loop over all phi layers, accumulate the denominator

               if( layer_phi[j]==0 ) 
				{	++n_trkperchl[sector+16][j][-1*avgchlphi+48];

				}


		}
				
            for(UInt_t j=0; j<layer_eta.size(); ++j) { // loop over all eta layers, accumulate the denominator

               if( layer_eta[j]==0 ) 
				{	++n_trkperchl[sector+16][j][avgchleta+48];

				}
					
         	}
		

/*            for(UInt_t j=0; j<layer_eta_modified.size(); ++j) { // loop over all eta layers, accumulate the denominator

               if( layer_eta_modified[j]==0 ) 
				{	++n_trkperchl_modified[sector+16][j][avgchleta_modified+48];

				}
					
         	}
*/

		} 

		 else { // if not matched
         h2_trk_P_eta_1->Fill(trkP->at(goodtrk_index[i])/1.e3, TMath::Abs(trkEta->at(goodtrk_index[i])));
         h2_trk_P_phi_1->Fill(trkP->at(goodtrk_index[i])/1.e3, goodtrk_phi[i]);
         h2_trk_eta_phi_1->Fill(TMath::Abs(trkEta->at(goodtrk_index[i])), goodtrk_phi[i]);
         }


         if( sector_2ndary!=0 ) { // if has overlap, this good trk is counted twice
//			cout << "sector_2ndary: " << sector_2ndary << endl; 
            h_overlap->Fill(TMath::Abs(sector-sector_2ndary));
            h_trk_sec->Fill(sector_2ndary);
            h2_trk_sec_phi->Fill(sector_2ndary, goodtrk_phi[i]);
            // but trk_P_eta and trk_P_phi will not be filled twice

            ++n_trkpersec[sector_2ndary+16];
            ++n_trkpersec_noStuck[sector_2ndary+16];
            Int_t avgchlphi = 0, avgchleta = 0;
//			Int_t avgchleta_modified = 0; 
            Int_t n_chlphi = 0, n_chleta = 0;
//			Int_t n_chleta_modified = 0; 
            vector<Int_t> layer_phi; layer_phi.resize(4, 0);
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector_2ndary ) continue;

//			   if(hit_pstrip->at(j)>0 && hit_pstrip->at(j)<=48){					
               		++n_trkperchl[sector_2ndary+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
//				}

/*				else{
					++n_trkperchl[sector_2ndary+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
               		++n_trkperchl_modified[sector_2ndary+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
				}
*/
               if( hit_measphi->at(j)==1 ) {
                  avgchlphi += hit_pstrip->at(j);
                  ++n_chlphi;
                  ++layer_phi[int(hit_wlay->at(j))-1];

               } else {

//						if(hit_pstrip->at(j)>0 && hit_pstrip->at(j)<=48){
		              		avgchleta += hit_pstrip->at(j);
		              		++n_chleta;
		              		++layer_eta[int(hit_wlay->at(j))-1];
//						}

/*						else{
		              		avgchleta += hit_pstrip->at(j);
		              		avgchleta_modified += hit_pstrip->at(j);
		              		++n_chleta;
		              		++n_chleta_modified;
		              		++layer_eta[int(hit_wlay->at(j))-1];		
		              		++layer_eta_modified[int(hit_wlay->at(j))-1];	
						}							
*/
               }
            }

            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
//            if( avgchleta_modified!=0) avgchleta_modified = int(avgchleta_modified/n_chleta_modified);

            for(UInt_t j=0; j<layer_phi.size(); ++j) {
               if( layer_phi[j]==0 ) 
					++n_trkperchl[sector_2ndary+16][j][-1*avgchlphi+48];

			}

            for(UInt_t j=0; j<layer_eta.size(); ++j) {
               if( layer_eta[j]==0 ) 
//					++n_trkperchl[sector_2ndary+16][j][avgchlphi+48];
					++n_trkperchl[sector_2ndary+16][j][avgchleta+48];

         	}

/*            for(UInt_t j=0; j<layer_eta_modified.size(); ++j) {
               if( layer_eta_modified[j]==0 ) 
					++n_trkperchl_modified[sector_2ndary+16][j][avgchlphi_modified+48];

         	}
*/
		}

      } // end loop of good trks
      if(matchedtrk_index.size()==0) continue; // no matched trk, pass
//	  cout << "cut2: " << goodtrk_index.size() << endl; 
//	  cout << "hitToMuon: " << hitToMuon->size() << endl; 
      ++n_cut2;
      h_cutflow->Fill(2);
      // ===============================================================================================

      // ===============================================================================================
      for(UInt_t i=0; i<hit_sector->size(); ++i) { // loop over hits
         bool matched = false;
         Double_t badphi = 0., goodphi = 0.; // to get rid of mismatched hits
         for(UInt_t j=0; j<matchedtrk_index.size(); ++j) {// loop over matched trks
            if( hitToMuon->at(i)==matchedtrk_index[j] ) {
               matched = true;
               goodphi = matchedtrk_phi[j];
            }
         } // end loop of matched trks
         for(UInt_t j=0; j<morethan2trk_index.size(); ++j) {
            if( hitToMuon->at(i)==morethan2trk_index[j] ) continue;
         }
         if( !matched ) continue;
         if( hit_sector->at(i)<0 ) badphi = -(pi/15)*(2*hit_sector->at(i)+17);
         else badphi = (pi/15)*(2*hit_sector->at(i)-17);
         if( TMath::Abs(badphi-goodphi) <= 0.7 || TMath::Abs(badphi+2*pi-goodphi) <= 0.7) matched = false;
         if( !matched ) continue;
						
         Float_t fillnumber = hit_sector->at(i)+0.25*(hit_wlay->at(i)-1); // a combination of sector # and layer #

         if( hit_measphi->at(i)==1 ) {
            h_hit_channel->Fill(-1.*hit_pstrip->at(i));
            h_eff_phi->Fill(fillnumber);
			h_eff_phi_noStuck->Fill(fillnumber);
 
            h2_eff_channel->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);
			h2_eff_channel_onlyNumerator->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);
//			h2_eff_channel_modified->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);

         } else {

//		   		if(hit_pstrip->at(i)>0 && hit_pstrip->at(i)<=48){
		        	h_hit_channel->Fill(hit_pstrip->at(i));
		        	h_eff_eta->Fill(fillnumber);
		        	h_eff_eta_noStuck->Fill(fillnumber);

		        	h2_eff_channel->Fill(hit_pstrip->at(i), fillnumber, 1);
					h2_eff_channel_onlyNumerator->Fill(hit_pstrip->at(i), fillnumber, 1);
//				}

/*				else{
		        	h_hit_channel->Fill(hit_pstrip->at(i));
		        	h_eff_eta->Fill(fillnumber);
		        	h_eff_eta_noStuck->Fill(fillnumber);

		        	h2_eff_channel->Fill(hit_pstrip->at(i), fillnumber, 1);
					h2_eff_channel_onlyNumerator->Fill(hit_pstrip->at(i), fillnumber, 1);
		        	h2_eff_channel_modified->Fill(hit_pstrip->at(i), fillnumber, 1);

				}
*/

         }

  
	  	 TString histoname = TString::Format("qsum_sector_%d_layer_%d", hit_sector->at(i), hit_wlay->at(i));
      
         TH1F *myhist = ((TH1F *)(HList->FindObject(histoname)));

         if (!myhist)
         {
//          	myhist = new TH1F(histoname, histoname; "qsum";;, 100, 0., 10000000.);
          	myhist = new TH1F(histoname, histoname, 100, 0., 6000000.);
          	HList->Add(myhist);
         } // if (!myhist) ...
      
      	 myhist->Fill(hit_qsum->at(i));
//		 myhist->Fit("landau"); 
//		 HList.push_back(myhist); 

     } // end loop of hits




      // ===============================================================================================

   } // end loop of events


//for(Int_t i=0; i<HList->GetEntries();++i)
//for(Int_t i=0; i<HList->size();++i)
//{
//	(HList[i])->Fit("landau"); 
//}

TH1F *source = (TH1F*)HList->First();
while(source)
{
	source->Fit("landau"); 
	source = (TH1F*)(HList->After(source));
}


//cout << "total events: " << Long64_t(h_cutflow->GetBinContent(1)) << endl; 
//cout << "good trk index: " << Long64_t(h_cutflow->GetBinContent(2)) << endl;
//cout << "matched trk index: " << Long64_t(h_cutflow->GetBinContent(3)) << endl; 


/*		 for(UInt_t m=0; m<secPhi.size(); ++m)
		 {
			UInt_t sectorNumberPhi = secPhi[m];
			--n_trkpersec_noStuck[sectorNumberPhi];
		 }

		 for(UInt_t m=0; m<secEta.size(); ++m)
		 {
			UInt_t sectorNumberEta = secEta[m];
			--n_trkpersec_noStuck[sectorNumberEta];
		 }
*/

/*	for(UInt_t i=0; i<n_trkpersec.size(); ++i)
	{
		cout << i <<"	" << n_trkpersec[i] << "	" << n_trkpersec_noStuck[i] << endl;  
	} 
*/


/*
   // useless fit
   TF1 *f_trk_phi = new TF1("f_trk_phi", "[0]",-3.3,3.3);
   h_trk_phi->Fit(f_trk_phi, "R");
*/

   // sector # | vector i | bin # (big block, block i is bin (i-1)*5+1~i*5)
   //   1~16   |  17~32   | 19~34
   // -16~-1   |   0~15   |  2~17
   for(UInt_t i=0; i<n_trkpersec.size(); ++i) { // normalization
      for(UInt_t j=0; j<4; ++j) {

		for(UInt_t m=0; m<secVector.size(); ++m)
		{
			if(secVector[m]==i && layerVector[m]==j)
				--n_trkpersec_noStuck[i];
		}

/*		for(UInt_t q=0; q<secEta.size(); ++q)
		{
			if(secEta[q]==i && layerEta[q]==j)
				--n_trkpersec_noStuck[i];
		}
*/

         if( h_eff_phi->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_phi->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
            h_eff_phi->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_phi->SetBinError((i+1)*5+j+1, binerror);
         }

         if( h_eff_phi_noStuck->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_phi_noStuck->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec_noStuck[i]-bincontent)/TMath::Power(n_trkpersec_noStuck[i], 3);
            h_eff_phi_noStuck->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec_noStuck[i]);
            h_eff_phi_noStuck->SetBinError((i+1)*5+j+1, binerror);
         }

         if( h_eff_eta->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_eta->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
            h_eff_eta->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_eta->SetBinError((i+1)*5+j+1, binerror);
         }

         if( h_eff_eta_noStuck->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_eta_noStuck->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec_noStuck[i]-bincontent)/TMath::Power(n_trkpersec_noStuck[i], 3);
            h_eff_eta_noStuck->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec_noStuck[i]);
            h_eff_eta_noStuck->SetBinError((i+1)*5+j+1, binerror);
         }

      }
   }


   // channel # | vector i | bin # (big block, block i is bin (i-1)*5+1~i*5)
   //   1~192   |  49~240  | 51~242
   // -48~-1    |   0~47   |  2~49
   for(UInt_t i=0; i<n_trkperchl.size(); ++i) { // normalization
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
         for(UInt_t k=0; k<n_trkperchl[i][j].size(); ++k) {
//			cout << i << "	" << j << "	" << k << "	" << n_trkperchl[i][j][k] << endl; 
            if( h2_eff_channel->GetBinContent( k+2, (i+1)*5+j+1 )!=0 ) {
               Double_t bincontent = h2_eff_channel->GetBinContent( k+2, (i+1)*5+j+1 );
               Double_t binerror = bincontent*(n_trkperchl[i][j][k]-bincontent)/TMath::Power(n_trkperchl[i][j][k], 3);
               h2_eff_channel->SetBinContent( k+2, (i+1)*5+j+1, bincontent/n_trkperchl[i][j][k] );
               h2_eff_channel->SetBinError( k+2, (i+1)*5+j+1, binerror );
            }
//			else {
//				cout << "i= " << i << "	" << "j= " << j << "	"<< "k= " << k << "	" << "bin k+2, (i+1)*5+j+1= " << k+2 <<", " << (i+1)*5+j+1 << endl;
//				cout << i <<"	" << j << "	" << k << endl;
//			}
         }
      }
   }


/*   for(UInt_t i=0; i<n_trkperchl_modified.size(); ++i) { // normalization
      for(UInt_t j=0; j<n_trkperchl_modified[i].size(); ++j) {
         for(UInt_t k=0; k<n_trkperchl_modified[i][j].size(); ++k) {
            if( h2_eff_channel_modified->GetBinContent( k+2, (i+1)*5+j+1 )!=0 ) {
               Double_t bincontent = h2_eff_channel_modified->GetBinContent( k+2, (i+1)*5+j+1 );
               Double_t binerror = bincontent*(n_trkperchl_modified[i][j][k]-bincontent)/TMath::Power(n_trkperchl_modified[i][j][k], 3);
               h2_eff_channel_modified->SetBinContent( k+2, (i+1)*5+j+1, bincontent/n_trkperchl_modified[i][j][k] );
               h2_eff_channel_modified->SetBinError( k+2, (i+1)*5+j+1, binerror );
            }

         }
      }
   }
*/

   outputfile->Write();
   outputfile->Close();
}
