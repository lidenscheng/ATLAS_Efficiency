#define EffStudy_cxx
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TImage.h> 
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

/*   ifstream inputPhi ("PhiExcludeList.dat", ifstream::in);
   ifstream inputEta ("EtaExcludeList.dat", ifstream::in);

   ifstream inputEtaEffLayer("output.txt_eta", ifstream::in); 
*/

//   ifstream inputEtaEffNormalLayers("etaEffNormalLayers.txt", ifstream::in);
//   ifstream inputEtaLowerLayers("etaEff-60VLayer.txt", ifstream::in); 

/*   vector<UInt_t> secPhi, secEta; 
   vector<UInt_t> layerPhi, layerEta;
   vector<UInt_t> channelPhi, channelEta; 
   UInt_t sec, layer, channel;  

   vector<UInt_t> binVector;
   vector<Double_t> effVector; 
   UInt_t bin;
   Double_t eff; 
*/

//   vector<Double_t> normalLayersEff, loweredLayersEff;
//   Double_t normalEff, loweredEff; 

/*   while(inputPhi >> sec >> layer >> channel)
   	{

		secPhi.push_back(sec);
		layerPhi.push_back(layer);
		channelPhi.push_back(channel);

   }

   while(inputEta >> sec >> layer >> channel)
   {
		secEta.push_back(sec);
		layerEta.push_back(layer);
		channelEta.push_back(channel);
   }

  while(inputEtaEffLayer >> bin >> eff)
  {	
		binVector.push_back(bin);
		effVector.push_back(eff);
  }
*/

/*
  while(inputEtaEffNormalLayers >> normalEff)
	normalLayersEff.push_back(normalEff); 

  while(inputEtaLowerLayers >> loweredEff)
	loweredLayersEff.push_back(loweredEff); 
*/

   
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
		 n_trkperchl[i][j].resize(241, 0); // magic number 242=48+192+1, 48 phi channels, 192 eta channels, 1 for convenience
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
//   TH2F * h2_eff_channel  = new TH2F("eff_channel", "eff_channel", 49, -49, 0, 170, -17, 17);

//eta and phi eff by layer without stuck-bit groups 
//   TH1F * h_eff_phi_noStuck      = new TH1F("eff_phi_noStuck", "eff_phi_noStuck", 170, -17, 17);
//   TH1F * h_eff_eta_noStuck      = new TH1F("eff_eta_noStuck", "eff_eta_noStuck", 170, -17, 17);

//eta and phi eff by layer without problematic areas 
   TH1F * h_eff_phi_good      = new TH1F("eff_phi_good", "eff_phi_good", 170, -17, 17);
   TH1F * h_eff_eta_good      = new TH1F("eff_eta_good", "eff_eta_good", 170, -17, 17);

   TH2F * h2_eff_channel_onlyNumerator = new TH2F("eff_channel_numerator", "eff_channel_numerator", 242, -49, 193, 170, -17, 17);

//   TH1F * h_eta_eff_normal_layers = new TH1F("h_eta_eff_normal_layers", "h_eta_eff_normal_layers", 11, 0.9, 1.0); 
//   TH1F * h_eta_eff_-60V_layers = new TH1F("h_eta_eff_-60V_layers", "h_eta_eff_-60V_layers", 11, 0.9, 1.0); 

//   TH1F * h_trk_sfit_unspoiled = new TH1F("trk_sfit_unspoiled", "trk_sfit_unspoiled", 9, -0.5, 8.5);
//   TH1F * h_trk_sfit_spoiled = new TH1F("trk_sfit_spoiled", "trk_sfit_spoiled", 9, -0.5, 8.5);
//   TH2F * h2_trk_sfit_both = new TH2F("trk_sfit_both", "trk_sfit_both", 9, -0.5, 8.5, 9, -0.5, 8.5); 

//   TH1F * h_eff_channel_forOneLayer = new TH1F("eff_channel_forOneLayer", "eff_channel_forOneLayer", 242, -49, 193); 
//   TH1F * h_eff_channel_forOneLayer = new TH1F("eff_channel_forOneLayer", "eff_channel_forOneLayer", 48, 0, 47); 

//   TH2F * h2_eta_eff_layer_etaValue = new TH2F("eta_eff_layer_etaValue", "eta_eff_layer_etaValue", 170, -17, 17, 70, 2., 2.7); 
//   TH2F * h2_eta_eff_layer_eta_regions = new TH2F("eta_eff_layer_eta_regions", "eta_eff_layer_eta_regions", 193, 0, 193, 170, -17, 17); 

//   TH2F * h2_qpeak_channelNumber = new TH2F("qpeak_channelNumber", "qpeak_channelNumber", 242, -49, 193, 3200, 0., 3200000.);  
//	 TH2F * h2_qpeak_channelNumber = new TH2F("qpeak_channelNumber", "qpeak_channelNumber", 242, -49, 193, 100, 0., 100000.);  

//   TH1F * h_phi_channel_eff  = new TH1F("phi_channel_eff", "phi_channel_eff", 170, -17, 17);
//   TH1F * h_phi_channel_eff  = new TH1F("phi_channel_eff", "phi_channel_eff", 100, 0., 1.0);

//   TH2F * h2_all_phi_channels_eff = new TH2F("all_phi_channels_eff", "all_phi_channels_eff", 100, 0., 1.0, 48, 0, 47); 

//   TH1F * h_hit_channel_phi  = new TH1F("hit_channel_phi", "hit_channel_phi", 49, -49, 0);
//   TH2F * h2_eff_channel_phi  = new TH2F("eff_channel_phi", "eff_channel_phi", 49, -49, 0, 170, -17, 17);

//	 TH1F * h_channel_segment = new TH1F("channel_segment", "channel_segment", 242, -49, 193); 


//histograms to show 4 layers of a particular sector just to see what each layer generally looks like; take out of code after quick check  
//   TH1F * h_hit_sector10_layer1  = new TH1F("hit_sector10_layer1", "hit_sector10_layer1", 242, -49, 193);
//   TH1F * h_hit_sector10_layer2  = new TH1F("hit_sector10_layer2", "hit_sector10_layer2", 242, -49, 193);
//   TH1F * h_hit_sector10_layer3  = new TH1F("hit_sector10_layer3", "hit_sector10_layer3", 242, -49, 193);
//   TH1F * h_hit_sector10_layer4  = new TH1F("hit_sector10_layer4", "hit_sector10_layer4", 242, -49, 193);

//   TObjArray *HList = new TObjArray(0); 
//vector<TH1F*> HList; 
   TList* HList = new TList(); 

   h_trk_phi->Sumw2();
   h_trk_sec->Sumw2();

//   TH1F * h_eta_cut0 = new TH1F("eta_cut0", "eta_cut0", 66, -3.3, 3.3);
//   TH1F * h_eta_cut1 = new TH1F("eta_cut1", "eta_cut1", 66, -3.3, 3.3);
//   TH1F * h_eta_cut2 = new TH1F("eta_cut2", "eta_cut2", 66, -3.3, 3.3);

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

      if (jentry%100000==0) cout<<"Event: "<<jentry<<endl;
      h_cutflow->Fill(0);


//	  for(UInt_t i=0; i<trkP->size(); ++i) {
//		 h_eta_cut0->Fill(trkEta->at(i)); 
//	  }
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
//		 h_eta_cut1->Fill(trkEta->at(i)); 
         Double_t phi = 0.;
         phi = TMath::ATan2(trkPy->at(i), trkPx->at(i));
         goodtrk_phi.push_back(phi);
         h_trk_phi->Fill(phi);
      } // end loop of trks

      if(goodtrk_index.size()==0) continue; // no good trk, pass
//	  cout << "cut1: " << goodtrk_index.size() << endl; 
      ++n_cut1;
      h_cutflow->Fill(1);


//count up unspoiled and spoiled hits on each track 
/*        	for(UInt_t i=0; i<goodtrk_index.size(); ++i) {
			UInt_t sfitUnspoiled = 0; 
			UInt_t sfitSpoiled = 0; 

         	for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits

           		if(hitToMuon->at(j)==goodtrk_index[i] && hit_measphi->at(j)==0) {  //don't need phi hits, only need eta hits 

					if(int(hit_sfit->at(j))==0) sfitUnspoiled++;  
					if(int(hit_sfit->at(j))!=0) sfitSpoiled++;    
 					}
			}

			h_trk_sfit_unspoiled->Fill(sfitUnspoiled); 
			h_trk_sfit_spoiled->Fill(sfitSpoiled); 
			h2_trk_sfit_both->Fill(sfitUnspoiled, sfitSpoiled); 
//			cout << sfitUnspoiled << "	" << sfitSpoiled << endl; 	
		}
*/

      // ===============================================================================================

      // ===============================================================================================
      // CUT2: at least one matched good trk

      vector<UInt_t>   matchedtrk_index;
      vector<Double_t> matchedtrk_phi;
      vector<UInt_t>   morethan2trk_index; // tracks that are matched by more than 2 sectors, should be good ones but dunno how to remove ill matched sectors
      for(UInt_t i=0; i<goodtrk_index.size(); ++i) { // loop over good trks
//		 cout << i << "	" << goodtrk_index.size() << endl; 
         Int_t sector = 0;
         Int_t sector_2ndary = 0;

//         vector<double>::iterator lowPhiIndex;
//         vector<double>::iterator lowEtaIndex;

//		UInt_t sfitNum = 0;

         bool morethan2 = false; // if matched by more than 2 sectors
         for(UInt_t j=0; j<hit_sector->size(); ++j) {  // loop over hits

            if( hitToMuon->at(j)==goodtrk_index[i] ) {
//            if( hitToSeg->at(j)==goodtrk_index[i] ) {

//				if(int(hit_sfit->at(j))==0) sfitNum++;  
//				cout << int(hit_sfit->at(j)) << endl; 
 
               Double_t badphi = 0.; // to exclude mismatching
               if( hit_sector->at(j)<0 ) badphi = -(pi/15)*(2*hit_sector->at(j)+17);
               else badphi = (pi/15)*(2*hit_sector->at(j)-17);
               if( TMath::Abs(badphi-goodtrk_phi[i]) > 0.7 && TMath::Abs(badphi+2*pi-goodtrk_phi[i]) > 0.7 ) {
                  sector = hit_sector->at(j);
                  break;
               }
            }
         } // end loop of hits

//			h_trk_sfit->Fill(sfitNum); 
//			cout << sfitNum << endl; //number of unspoiled hits on each track 

         for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits again to find overlaped sector
            if( hitToMuon->at(j)==goodtrk_index[i] && hit_sector->at(j)!=sector && sector!=0 ) {
//            if( hitToSeg->at(j)==goodtrk_index[i] && hit_sector->at(j)!=sector && sector!=0 ) {
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
//            if( hitToSeg->at(j)==goodtrk_index[i] && hit_sector->at(j)!=sector && hit_sector->at(j)!=sector_2ndary && sector!=0 && sector_2ndary!=0 ) {
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
//		UInt_t sfitNum = 0;
            matchedtrk_index.push_back(goodtrk_index[i]);
            matchedtrk_phi.push_back(goodtrk_phi[i]);
            h_trk_sec->Fill(sector);
            h2_trk_sec_phi->Fill(sector, goodtrk_phi[i]);
            h2_trk_P_eta_0->Fill(trkP->at(goodtrk_index[i])/1.e3, TMath::Abs(trkEta->at(goodtrk_index[i])));
            h2_trk_P_phi_0->Fill(trkP->at(goodtrk_index[i])/1.e3, goodtrk_phi[i]);
            h2_trk_eta_phi_0->Fill(TMath::Abs(trkEta->at(goodtrk_index[i])), goodtrk_phi[i]);
//		    h_eta_cut2->Fill(trkEta->at(goodtrk_index[i])); 

            ++n_trkpersec[sector+16]; // why 16? see the last part of this code for conversion rule; ex if sec -16, then it would be 0 here
 
	    ++n_trkpersec_phi_good[sector+16];
	    ++n_trkpersec_eta_good[sector+16];

            Int_t avgchlphi = 0, avgchleta = 0; // average phi channel, average eta channel
            Int_t n_chlphi = 0, n_chleta = 0; // number of phi hits and eta hits, used to calculate average
            vector<Int_t> layer_phi; layer_phi.resize(4, 0); // number of hits on each layer
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector ) continue;
//               if( hitToSeg->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector ) continue;

//				if(int(hit_sfit->at(j))==0) sfitNum++;  
	
	           ++n_trkperchl[sector+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48]; //n_trkperchl[sector#][layer#][channel#], conversion rule down below 

				
               if( hit_measphi->at(j)==1 ) {
               		avgchlphi += hit_pstrip->at(j);
                 	++n_chlphi;
                  	++layer_phi[int(hit_wlay->at(j))-1]; 

               } else{
		              avgchleta += hit_pstrip->at(j);
		              ++n_chleta;
		              ++layer_eta[int(hit_wlay->at(j))-1];
               }
            } // end loop of hits

//			h_trk_sfit->Fill(sfitNum); 
//			cout << sfitNum << endl; //number of unspoiled hits on each track

		
            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi); 
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
			 

            for(UInt_t j=0; j<layer_phi.size(); ++j) { // loop over all phi layers, accumulate the denominator

//				cout << "Event: "<<jentry <<"	" << "track " << i << "	" << "hit " << j << "	" << "sfit= " << int(hit_sfit->at(j)) << endl; 	

               if( layer_phi[j]==0 ) 
				{	
					++n_trkperchl[sector+16][j][-1*avgchlphi+48];

				}

			}
				
            for(UInt_t j=0; j<layer_eta.size(); ++j) { // loop over all eta layers, accumulate the denominator

               if( layer_eta[j]==0 ) 
				{	
					++n_trkperchl[sector+16][j][avgchleta+48];

				}
//				cout << "Event: "<<jentry <<"	" << "track " << i << "	" << "hit " << j << "	" << "sfit= " << int(hit_sfit->at(j)) << endl; 	
					
         	}
		

		} 

		 else { // if not matched
         h2_trk_P_eta_1->Fill(trkP->at(goodtrk_index[i])/1.e3, TMath::Abs(trkEta->at(goodtrk_index[i])));
         h2_trk_P_phi_1->Fill(trkP->at(goodtrk_index[i])/1.e3, goodtrk_phi[i]);
         h2_trk_eta_phi_1->Fill(TMath::Abs(trkEta->at(goodtrk_index[i])), goodtrk_phi[i]);
//			cout << "I got here" << endl; 
         }


         if( sector_2ndary!=0 ) { // if has overlap, this good trk is counted twice
//		UInt_t sfitNum = 0;
            h_overlap->Fill(TMath::Abs(sector-sector_2ndary));
            h_trk_sec->Fill(sector_2ndary);
            h2_trk_sec_phi->Fill(sector_2ndary, goodtrk_phi[i]);
            // but trk_P_eta and trk_P_phi will not be filled twice

            ++n_trkpersec[sector_2ndary+16];
	    	++n_trkpersec_phi_good[sector_2ndary+16];
	    	++n_trkpersec_eta_good[sector_2ndary+16];

            Int_t avgchlphi = 0, avgchleta = 0;
            Int_t n_chlphi = 0, n_chleta = 0;
            vector<Int_t> layer_phi; layer_phi.resize(4, 0);
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector_2ndary ) continue;
//               if( hitToSeg->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector_2ndary ) continue;

//				if(int(hit_sfit->at(j))==0) sfitNum++;  

				++n_trkperchl[sector_2ndary+16][int(hit_wlay->at(j))-1][TMath::Power(-1, int(hit_measphi->at(j)))*int(hit_pstrip->at(j))+48];
				
               if( hit_measphi->at(j)==1 ) {
                  avgchlphi += hit_pstrip->at(j);
                  ++n_chlphi;
                  ++layer_phi[int(hit_wlay->at(j))-1];

               } else {
		              avgchleta += hit_pstrip->at(j);
		              ++n_chleta;
		              ++layer_eta[int(hit_wlay->at(j))-1];

               }
            }

//			h_trk_sfit->Fill(sfitNum); 
//			cout << sfitNum << endl; //number of unspoiled hits on each track

            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);

            for(UInt_t j=0; j<layer_phi.size(); ++j) {
               if( layer_phi[j]==0 ) 
					++n_trkperchl[sector_2ndary+16][j][-1*avgchlphi+48];

			}

            for(UInt_t j=0; j<layer_eta.size(); ++j) {
               if( layer_eta[j]==0 ) 
//					++n_trkperchl[sector_2ndary+16][j][avgchlphi+48]; //the bug is here; YZ had avgchlphi instead of avgchleta here 
					++n_trkperchl[sector_2ndary+16][j][avgchleta+48];

         	}

		} 
		 

      } // end loop of good trks

      if(matchedtrk_index.size()==0) continue; // no matched trk, pass
//	  cout << "cut2: " << goodtrk_index.size() << endl; 
//	  cout << "hitToMuon: " << hitToMuon->size() << endl; 
      ++n_cut2;
      h_cutflow->Fill(2);



//count up spoiled and unspoiled hits for each track 
/*         for(UInt_t i=0; i<matchedtrk_index.size(); ++i) {// loop over matched trks
			UInt_t sfitUnspoiled = 0; 
			UInt_t sfitSpoiled = 0; 

         	for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits

           		if(hitToMuon->at(j)==matchedtrk_index[i] && hit_measphi->at(j)==0) {  //don't need phi hits, only need eta hits 

					if(int(hit_sfit->at(j))==0) sfitUnspoiled++;  
					if(int(hit_sfit->at(j))!=0) sfitSpoiled++;    
 					}
			}

			h_trk_sfit_unspoiled->Fill(sfitUnspoiled); 
			h_trk_sfit_spoiled->Fill(sfitSpoiled); 
			h2_trk_sfit_both->Fill(sfitUnspoiled, sfitSpoiled); 
//			cout << sfitUnspoiled << "	" << sfitSpoiled << endl; 	
		}
*/
      // ===============================================================================================

      // ===============================================================================================
      for(UInt_t i=0; i<hit_sector->size(); ++i) { // loop over hits

         bool matched = false;
         Double_t badphi = 0., goodphi = 0.; // to get rid of mismatched hits

         for(UInt_t j=0; j<matchedtrk_index.size(); ++j) {// loop over matched trks
            if( hitToMuon->at(i)==matchedtrk_index[j] ) {
//            if( hitToSeg->at(i)==matchedtrk_index[j] ) {

//         for(UInt_t j=0; j<goodtrk_index.size(); ++j) {// loop over matched trks
//            if( hitToMuon->at(i)==goodtrk_index[j] ) {



               matched = true;
               goodphi = matchedtrk_phi[j];
//		 goodphi = goodtrk_phi[j]; 
            }
			
         } // end loop of matched trks

//         for(UInt_t j=0; j<morethan2trk_index.size(); ++j) {
//            if( hitToMuon->at(i)==morethan2trk_index[j] ) continue;
//         }
         if( !matched ) continue;
         if( hit_sector->at(i)<0 ) badphi = -(pi/15)*(2*hit_sector->at(i)+17);
         else badphi = (pi/15)*(2*hit_sector->at(i)-17);
         if( TMath::Abs(badphi-goodphi) <= 0.7 || TMath::Abs(badphi+2*pi-goodphi) <= 0.7) matched = false;
         if( !matched ) continue;

						
         Float_t fillnumber = hit_sector->at(i)+0.25*(hit_wlay->at(i)-1); // a combination of sector # and layer #

         if( hit_measphi->at(i)==1 ) {

            h_hit_channel->Fill(-1.*hit_pstrip->at(i));
//            h_hit_channel_phi->Fill(-1.*hit_pstrip->at(i));

            h_eff_phi->Fill(fillnumber);
	    	h_eff_phi_good->Fill(fillnumber);
 
            h2_eff_channel->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);
	    	h2_eff_channel_onlyNumerator->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);
//			h2_eff_channel_phi->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);

//			cout << "Event: "<<jentry << "	" << "hit " << i << "	" << "sfit= " << int(hit_sfit->at(i)) << endl; 	


/*			if(fillnumber==-12.25)
//				h_eff_channel_forOneLayer->Fill((-1.*hit_pstrip->at(i))+48);  
				h_eff_channel_forOneLayer->Fill(-1.*hit_pstrip->at(i));  
*/		




         } else {

	    h_hit_channel->Fill(hit_pstrip->at(i));	
			
	    h_eff_eta->Fill(fillnumber);
	    h_eff_eta_good->Fill(fillnumber);

            h2_eff_channel->

Fill(hit_pstrip->at(i), fillnumber, 1);
//	    h2_eff_channel_onlyNumerator->Fill(hit_pstrip->at(i), fillnumber, 1);
       
	    }


/*	    if(hit_sector->at(i)==-13 and hit_wlay->at(i)==4)
	    {
		if(hit_measphi->at(i)==1)
		{
//			cout << -1.*hit_pstrip->at(i) << "	" << hit_qpeak->at(i) << endl; 
			h2_qpeak_channelNumber-> Fill(-1.*hit_pstrip->at(i), hit_qpeak->at(i)); 
		}
		
		else
			h2_qpeak_channelNumber-> Fill(hit_pstrip->at(i), hit_qpeak->at(i)); 				
	    }

*/		


//	 TString histoname = TString::Format("qsum_sector_%d_layer_%d", hit_sector->at(i), hit_wlay->at(i));
	 TString histoname = TString::Format("qpeak_sector_%d_layer_%d", hit_sector->at(i), hit_wlay->at(i));
      
         TH1F *myhist = ((TH1F *)(HList->FindObject(histoname)));

         if (!myhist)
         {
//          	myhist = new TH1F(histoname, histoname, 100, 0., 4000000.);
          	myhist = new TH1F(histoname, histoname, 100, 0., 6000000.);
//          myhist = new TH1F(histoname, histoname, 100, 0., 3300000.);
          	HList->Add(myhist);
         } 
      
//      	 myhist->Fill(hit_qsum->at(i));
       myhist->Fill(hit_qpeak->at(i));


     } // end loop of hits



      // ===============================================================================================

} // end loop of events





   // sector # | vector i | bin # (big block, block i is bin (i-1)*5+1~i*5)
   //   1~16   |  17~32   | 19~34
   // -16~-1   |   0~15   |  2~17
   for(UInt_t i=0; i<n_trkpersec.size(); ++i) { // normalization
      for(UInt_t j=0; j<4; ++j) {

         if( h_eff_phi->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_phi->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
            h_eff_phi->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_phi->SetBinError((i+1)*5+j+1, binerror);
         }

         if( h_eff_phi_good->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_phi_good->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
	        h_eff_phi_good->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_phi_good->SetBinError((i+1)*5+j+1, binerror);
         }


         if( h_eff_eta->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_eta->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
            h_eff_eta->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_eta->SetBinError((i+1)*5+j+1, binerror);
         }

         if( h_eff_eta_good->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_eta_good->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec[i]-bincontent)/TMath::Power(n_trkpersec[i], 3);
	        h_eff_eta_good->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec[i]);
            h_eff_eta_good->SetBinError((i+1)*5+j+1, binerror);
         }

//		cout << i << "	" << j <<"	" << h_eff_eta->GetBinContent((i+1)*5+j+1) << "	" << h_eff_eta_good->GetBinContent((i+1)*5+j+1) << endl; 

      }
   }

//update phi eff per layer by taking out stuck-bit groups and only averaging over leftover channels 
	h_eff_phi_good->SetBinContent(24,0.6725168333);
	h_eff_phi_good->SetBinContent(96,0.9329755833);
	h_eff_phi_good->SetBinContent(102,0.9368261389);
	h_eff_phi_good->SetBinContent(117,0.9145758889);
	h_eff_phi_good->SetBinContent(127,0.948701);
	h_eff_phi_good->SetBinContent(139,0.9257730556);
	h_eff_phi_good->SetBinContent(147,0.9549523611);
	h_eff_phi_good->SetBinContent(157,0.9501166667);


    h_eff_eta_good->SetBinContent(67, 0.9551121167); 
    h_eff_eta_good->SetBinContent(68, 0.9426124222); 
    h_eff_eta_good->SetBinContent(57, 0.9613561278) ; 
    h_eff_eta_good->SetBinContent(48, 0.9582207667); 
    h_eff_eta_good->SetBinContent(29, 0.9488718556); 
    h_eff_eta_good->SetBinContent(6, 0.9508530722); 
    h_eff_eta_good->SetBinContent(9, 0.9484674444); 
    h_eff_eta_good->SetBinContent(97, 0.957621); 
    h_eff_eta_good->SetBinContent(106, 0.9624483333); 
    h_eff_eta_good->SetBinContent(119, 0.9527053167); 
    h_eff_eta_good->SetBinContent(127, 0.95554915); 
    h_eff_eta_good->SetBinContent(136, 0.9600568929); 
    h_eff_eta_good->SetBinContent(146, 0.9626569048); 
    h_eff_eta_good->SetBinContent(148, 0.9567193833); 
    h_eff_eta_good->SetBinContent(82, 0.8572789889); 
    h_eff_eta_good->SetBinContent(64, 0.9511558778); 
    h_eff_eta_good->SetBinContent(31, 0.951617122); 
    h_eff_eta_good->SetBinContent(32, 0.9655316071); 
    h_eff_eta_good->SetBinContent(22, 0.9696598556); 
    h_eff_eta_good->SetBinContent(93, 0.9637039444); 
    h_eff_eta_good->SetBinContent(103, 0.9386933944); 
    h_eff_eta_good->SetBinContent(134, 0.9713682444); 
    h_eff_eta_good->SetBinContent(142, 0.9606662444) ; 
    h_eff_eta_good->SetBinContent(152, 0.9561088944); 
    h_eff_eta_good->SetBinContent(153, 0.9589092056); 
    h_eff_eta_good->SetBinContent(162, 0.96344795); 



   for(UInt_t i=0; i<n_trkpersec.size(); ++i) 
	{ 
      for(UInt_t j=0; j<4; ++j) 
		{
			cout << (i+1)*5+j+1 <<"	" << h_eff_phi_good->GetBinContent((i+1)*5+j+1) << "	" << h_eff_eta_good->GetBinContent((i+1)*5+j+1) << endl;
		}
	}



//update eff per layer by taking out stuck-bit groups and only averaging over channels 50-180 
/*    for(UInt_t i=0; i<binVector.size(); ++i)
	{
		h_eff_eta_good->SetBinContent(binVector[i], effVector[i]);
	}
*/


/*	for(UInt_t i=0; i<normalLayersEff.size(); ++i)
	{
		h_eta_eff_normal_layers->Fill(normalLayersEff[i]);
	}

	for(UInt_t i=0; i<loweredLayersEff.size(); ++i)
	{
		h_eta_eff_-60V_layers->Fill(loweredLayersEff[i]);
	}
*/

   // channel # | vector i | bin # (big block, block i is bin (i-1)*5+1~i*5)
   //   1~192   |  49~240  | 51~242
   // -48~-1    |   0~47   |  2~49
   for(UInt_t i=0; i<n_trkperchl.size(); ++i) { // normalization
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
         for(UInt_t k=0; k<n_trkperchl[i][j].size(); ++k) {

/*		for(UInt_t m=0; m<secPhi.size(); ++m)
		{
			if(secPhi[m]==i && layerPhi[m]==j && channelPhi[m]==k)
				--n_trkperchl[i][j][k];
		}
*/

            if( h2_eff_channel->GetBinContent( k+2, (i+1)*5+j+1 )!=0 ) {
               Double_t bincontent = h2_eff_channel->GetBinContent( k+2, (i+1)*5+j+1 );
               Double_t binerror = bincontent*(n_trkperchl[i][j][k]-bincontent)/TMath::Power(n_trkperchl[i][j][k], 3);
               h2_eff_channel->SetBinContent( k+2, (i+1)*5+j+1, bincontent/n_trkperchl[i][j][k] );
               h2_eff_channel->SetBinError( k+2, (i+1)*5+j+1, binerror );

	     	}
			


//			else {
//				cout << "i= " << i << "	" << "j= " << j << "	"<< "k= " << k << ", 	" << "bin k+2, (i+1)*5+j+1= " << k+2 <<", " << (i+1)*5+j+1 << endl;
//				cout << i <<"	" << j << "	" << k << endl;
//			}

         }
      }
   }


/*   for(UInt_t i=0; i<n_trkperchl.size(); ++i) { // normalization
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
         for(UInt_t k=0; k<48; ++k) {

            if( h2_eff_channel_phi->GetBinContent( k+2, (i+1)*5+j+1 )!=0 ) {
               Double_t bincontent = h2_eff_channel_phi->GetBinContent( k+2, (i+1)*5+j+1 );
               Double_t binerror = bincontent*(n_trkperchl[i][j][k]-bincontent)/TMath::Power(n_trkperchl[i][j][k], 3);
               h2_eff_channel_phi->SetBinContent( k+2, (i+1)*5+j+1, bincontent/n_trkperchl[i][j][k] );
               h2_eff_channel_phi->SetBinError( k+2, (i+1)*5+j+1, binerror );

	     	}
		}
	}
}
*/


/*	UInt_t i=3;
	UInt_t j=3; 

         for(UInt_t k=0; k<=48; ++k) {
		cout << k <<"	" << h2_eff_channel->GetBinContent( k+2, (i+1)*5+j+1 ) << endl; 
		
	}
*/



/*for(UInt_t k=0; k<48; ++k) {
	if( h2_eff_channel->GetBinContent(k+2, 97)!=0) {
		sum += h2_eff_channel->GetBinContent( k+2, 97);
		count++;
	}
}
*/

/*
for(UInt_t i=6; i<=169; ++i)
{
	Double_t sum=0.0;
	UInt_t count= 0;
	for(UInt_t k=0; k<48; ++k) {
		if( h2_eff_channel->GetBinContent(k+2, i)!=0) {
			sum += h2_eff_channel->GetBinContent( k+2, i);
			count++;
		}
	}

	cout << i << "	" << sum/count << endl; 

}
*/

/*
    Double_t sum=0.0;
    UInt_t count= 0;  
    for(UInt_t k=0; k<n_trkperchl[3][3].size(); ++k) {
	if( h_eff_channel_forOneLayer->GetBinContent(k+2)!=0) {
	    Double_t bincontent = h_eff_channel_forOneLayer->GetBinContent(k+2);
            Double_t binerror = bincontent*((n_trkperchl[3][3][k])-bincontent)/TMath::Power(n_trkperchl[3][3][k], 3);
            h_eff_channel_forOneLayer->SetBinContent( k+2, bincontent/(n_trkperchl[3][3][k]));
            h_eff_channel_forOneLayer->SetBinError( k+2, binerror );

	    sum += h_eff_channel_forOneLayer->GetBinContent(k+2);
	    count++;  
	}
     }

cout << "Phi Efficiency for this layer= " << sum/count << endl; 
*/

   outputfile->Write();
   outputfile->Close();
}
