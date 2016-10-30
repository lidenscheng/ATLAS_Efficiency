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

//   ifstream inputPhi ("exclude_phi.dat", ifstream::in);
//   ifstream inputEta ("exclude_eta.dat", ifstream::in);
   ifstream inputStuck("exclude.dat", ifstream::in);
   ifstream inputAll("allExclude.dat", ifstream::in);

   vector<UInt_t> secVectorStuck, layerVectorStuck, channelVectorStuck;
   vector<UInt_t> secVectorAll, layerVectorAll, channelVectorAll;
//   vector<Int_t> secPhi, secEta; 
//   vector<UInt_t> layerPhi, layerEta;
//   vector<Int_t> channelPhi, channelEta; 
   UInt_t sec, channel, layer;  

   while(inputStuck >> sec >> layer >> channel)
   	{
		secVectorStuck.push_back(sec);
		layerVectorStuck.push_back(layer);
		channelVectorStuck.push_back(channel);
   }

   while(inputAll >> sec >> layer >> channel)
   	{
		secVectorAll.push_back(sec);
		layerVectorAll.push_back(layer);
		channelVectorAll.push_back(channel);
   }

/*
   while(inputPhi >> sec >> layer >> channel)
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
*/

   
   Double_t pi = 3.14159265359;
   // cut flow stats
   Long_t n_cut1 = 0, n_cut2 = 0;
   vector<Long_t> n_trkpersec; // number of matched trks per sector, eff denominator
   vector<Long_t> n_trkpersec_noStuck; // number of matched trks per sector when we take out stuck bits, eff denominator
   vector<Long_t> n_trkpersec_good; // number of matched trks per sector when we take out all the problematic areas 

   vector<vector<vector<Long_t>>> n_trkperchl; // number of matched trks per channel, eff denominator
//   vector<vector<vector<Long_t>>> n_trkperchl_modified; // used to count for the case of ignoring ch 1-48 

   n_trkpersec.resize(33, 0); // magic number 33=16+16+1, 16 sectors on each end, the extra 1 is for convenience filling histogram
   n_trkpersec_noStuck.resize(33, 0); 
   n_trkpersec_good.resize(33, 0); 

   n_trkperchl.resize(33);
//   n_trkperchl_modified.resize(33);

   for(UInt_t i=0; i<n_trkperchl.size(); ++i) {
      n_trkperchl[i].resize(4); // 4 layers per sector
//      n_trkperchl_modified[i].resize(4); 
      for(UInt_t j=0; j<n_trkperchl[i].size(); ++j) {
		 n_trkperchl[i][j].resize(241, 0); // magic number 242=48+192+1, 48 phi channels, 192 eta channels, 1 for convenience
//		 n_trkperchl_modified[i][j].resize(241, 0);
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

//eta and phi eff by layer without stuck-bit groups 
   TH1F * h_eff_phi_noStuck      = new TH1F("eff_phi_noStuck", "eff_phi_noStuck", 170, -17, 17);
   TH1F * h_eff_eta_noStuck      = new TH1F("eff_eta_noStuck", "eff_eta_noStuck", 170, -17, 17);

//eta and phi eff by layer without problematic areas 
   TH1F * h_eff_phi_good      = new TH1F("eff_phi_good", "eff_phi_good", 170, -17, 17);
   TH1F * h_eff_eta_good      = new TH1F("eff_eta_good", "eff_eta_good", 170, -17, 17);

   TH2F * h2_eff_channel_onlyNumerator = new TH2F("eff_channel_numerator", "eff_channel_numerator", 242, -49, 193, 170, -17, 17);

/*
   TH1D * h_mpv  = new TH1D("most_prob_hist", "most_prob_hist", 170, -17, 17);
   TH2D * h_mpv  = new TH2D("most_prob_hist", "most_prob_hist", 170, -17, 17, 100, 0., 1000000.);

   double layer_array[132]= {-16.0, -15.75, -15.5, -15.25, -15.0, -14.75, -14.5, -14.25, -14.0, -13.75, -13.5, -13.25, -13.0, -12.75, -12.5, -12.25, -12.0, -11.75, -11.5, -11.25, -11.0, -10.75, -10.5, -10.25, -10.0, -9.75, -9.5, -9.25, -9.0, -8.75, -8.5, -8.25, -8.0, -7.75, -7.5, -7.25, -7.0, -6.75, -6.5, -6.25, -6.0, -5.75, -5.5, -5.25, -5.0, -4.75, -4.5, -4.25, -4.0, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75}; 

   double mpv_array[132] = {303528.06430068915, 352312.1005885047, 350774.8768631402, 204390.7400442789, 332765.47051792423, 334289.9647395915, 348354.2306202321, 180851.40209413983, 329537.6425937077, 318710.4369112644, 342563.8213782609, 340329.5529845818, 173116.58391859286, 392913.47899459855, 354905.1989591159, 214904.10915724907, 377868.43931294366, 321009.95516722155, 214409.68203069476, 183728.4792106583, 328815.3625490051, 322569.1625439248, 362663.5982391932, 191441.9424071328, 341938.230147599, 362179.5637415457, 387517.46743221494, 372201.0867254998, 313575.29882395035, 380175.2864507449, 363300.72843936813, 191085.49463926585, 343844.6678254668, 371972.47279121453, 370432.99780511187, 401399.20450133807, 326583.7071226696, 364776.4213860212, 388100.9158108409, 167770.93437876328, 323536.56964583695, 392566.429281291, 342821.85691914574, 383555.40418755007, 303965.7436669281, 181069.4357997908, 341174.6907721191, 338610.9752395102, 335791.4714803609, 388271.5704068667, 354319.25518126704, 411295.3019850741, 92932.47255924661, 0.0, 285759.88491068356, 397822.97152320575, 379709.67859140475, 447221.72002702724, 457370.24823870347, 448918.24520109105, 0.0, 74286.54742576943, 173871.2362560137, 0.0, 0.0, 0.0, 0.0, 0.0, 307739.89480060287, 328619.5032972728, 376349.08997509937, 115015.53752787919, 324076.4862743475, 370181.49211164395, 210481.08416029543, 346462.04313298635, 312111.02919527143, 423065.07137281727, 285315.87024469406, 349564.50348059693, 316240.4885295946, 347859.54419401474, 360631.39310485835, 369236.30550170224, 205331.20568286377, 137380.33712485654, 209223.02702841736, 354994.54525545787, 356571.01070032665, 359903.99379437853, 359133.417178487, 394549.4746655582, 315386.75993599923, 345435.13307391666, 395806.34481468156, 375112.99727746326, 346583.4228484598, 334304.4833079832, 406169.08131788834, 347452.9822731214, 104365.11727833198, 126660.56976989412, 324005.84552865714, 403222.14475406735, 344646.28974536044, 342444.0441511258, 387342.37480819295, 354714.95443417, 313268.4248884245, 368929.46375732386, 337322.0635396888, 191843.27220877656, 351474.3358056345, 369232.9849002818, 385277.89482436253, 342966.5483308096, 309296.2351441557, 193497.2971971612, 361996.0814791713, 191164.77990730427, 358277.1764685993, 216610.68750623317, 186638.36258466306, 196036.9922129175, 367640.0292386261, 377394.2353071206, 270498.7533845858, 222940.1918283263, 352200.74397950823, 348253.1679893707, 414654.1809835943, 197442.10861584812}; 

	for(UInt_t i=0; i<132; ++i)
	{
		h_mpv->Fill(layer_array[i], mpv_array[i], 1);
	} 
*/

//histograms to show 4 layers of a particular sector just to see what each layer generally looks like; take out of code after quick check  
//   TH1F * h_eff_sector10_layer1  = new TH1F("eff_sector10_layer1", "eff_sector10_layer1", 242, -49, 193);
//   TH1F * h_eff_sector10_layer2  = new TH1F("eff_sector10_layer2", "eff_sector10_layer2", 242, -49, 193);
//   TH1F * h_eff_sector10_layer3  = new TH1F("eff_sector10_layer3", "eff_sector10_layer3", 242, -49, 193);
//   TH1F * h_eff_sector10_layer4  = new TH1F("eff_sector10_layer4", "eff_sector10_layer4", 242, -49, 193);

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

      if (jentry%100000==0) cout<<"Event: "<<jentry<<endl;
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
                  cout<<"ERROR: more than two matched sectors!"<<endl;
                 cout<<"info: "<<goodtrk_index[i]<<"  "<<trkP->at(goodtrk_index[i])/1.e3<<"  "<<trkEta->at(goodtrk_index[i])<<"  "<<goodtrk_phi[i]<<endl;
                  cout<<"      "<<sector<<"  "<<sector_2ndary<<"  "<<hit_sector->at(j)<<endl;
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
			++n_trkpersec_good[sector+16];

            Int_t avgchlphi = 0, avgchleta = 0; // average phi channel, average eta channel
//			Int_t avgchleta_modified = 0; 
            Int_t n_chlphi = 0, n_chleta = 0; // number of phi hits and eta hits, used to calculate average
//			Int_t n_chleta_modified = 0; 
            vector<Int_t> layer_phi; layer_phi.resize(4, 0); // number of hits on each layer
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector ) continue;
		
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

		
            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
//            if( avgchleta_modified!=0) avgchleta_modified = int(avgchleta_modified/n_chleta_modified);

            for(UInt_t j=0; j<layer_phi.size(); ++j) { // loop over all phi layers, accumulate the denominator

               if( layer_phi[j]==0 ) 
				{	
					++n_trkperchl[sector+16][j][-1*avgchlphi+48];

				}

/*			   for(UInt_t m=0; m<secPhi.size(); ++m)
			   {
					if(secPhi[m]==sector+16 && layerPhi[m]==j-1)
						--n_trkpersec_noStuck[sector+16]; 
			   }
*/
			}
				
            for(UInt_t j=0; j<layer_eta.size(); ++j) { // loop over all eta layers, accumulate the denominator

               if( layer_eta[j]==0 ) 
				{	
					++n_trkperchl[sector+16][j][avgchleta+48];

				}

/*			   for(UInt_t m=0; m<secEta.size(); ++m)
			   {
					if(secEta[m]==sector+16 && layerEta[m]==j-1)
						--n_trkpersec_noStuck[sector+16]; 
			   }
*/					
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
            h_overlap->Fill(TMath::Abs(sector-sector_2ndary));
            h_trk_sec->Fill(sector_2ndary);
            h2_trk_sec_phi->Fill(sector_2ndary, goodtrk_phi[i]);
            // but trk_P_eta and trk_P_phi will not be filled twice

            ++n_trkpersec[sector_2ndary+16];
            ++n_trkpersec_noStuck[sector_2ndary+16];
			++n_trkpersec_good[sector_2ndary+16];

            Int_t avgchlphi = 0, avgchleta = 0;
//			Int_t avgchleta_modified = 0; 
            Int_t n_chlphi = 0, n_chleta = 0;
//			Int_t n_chleta_modified = 0; 
            vector<Int_t> layer_phi; layer_phi.resize(4, 0);
            vector<Int_t> layer_eta; layer_eta.resize(4, 0);
//            vector<Int_t> layer_eta_modified; layer_eta_modified.resize(4, 0);

            for(UInt_t j=0; j<hit_sector->size(); ++j) { // loop over hits
               if( hitToMuon->at(j)!=goodtrk_index[i] || hit_sector->at(j)!=sector_2ndary ) continue;

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

            if( avgchlphi!=0) avgchlphi = int(avgchlphi/n_chlphi);
            if( avgchleta!=0) avgchleta = int(avgchleta/n_chleta);
//            if( avgchleta_modified!=0) avgchleta_modified = int(avgchleta_modified/n_chleta_modified);

            for(UInt_t j=0; j<layer_phi.size(); ++j) {
               if( layer_phi[j]==0 ) 
					++n_trkperchl[sector_2ndary+16][j][-1*avgchlphi+48];

/*			   for(UInt_t m=0; m<secPhi.size(); ++m)
			   {
					if(secPhi[m]==sector_2ndary+16 && layerPhi[m]==j-1)
						--n_trkpersec_noStuck[sector_2ndary+16]; 
			   }
*/
			}

            for(UInt_t j=0; j<layer_eta.size(); ++j) {
               if( layer_eta[j]==0 ) 
//					++n_trkperchl[sector_2ndary+16][j][avgchlphi+48]; //the bug is here; YZ had avgchlphi instead of avgchleta here 
					++n_trkperchl[sector_2ndary+16][j][avgchleta+48];

/*			   for(UInt_t m=0; m<secEta.size(); ++m)
			   {
					if(secEta[m]==sector_2ndary+16 && layerEta[m]==j-1)
						--n_trkpersec_noStuck[sector_2ndary+16]; 
			   }
*/

         	}





/*            for(UInt_t j=0; j<layer_eta_modified.size(); ++j) {
               if( layer_eta_modified[j]==0 ) 
					++n_trkperchl_modified[sector_2ndary+16][j][avgchleta_modified+48];

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
			h_eff_phi_good->Fill(fillnumber);
 
            h2_eff_channel->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);
			h2_eff_channel_onlyNumerator->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);

// filling the histograms for 4 layers of a particular sector; just a quick check, don't keep in code 
/*			if(hit_sector->at(i)==10 && hit_wlay->at(i)==1) 
				h_eff_sector10_layer1->Fill(-1.*hit_pstrip->at(i)); 

			if(hit_sector->at(i)==10 && hit_wlay->at(i)==2) 
				h_eff_sector10_layer2->Fill(-1.*hit_pstrip->at(i));

			if(hit_sector->at(i)==10 && hit_wlay->at(i)==3) 
				h_eff_sector10_layer3->Fill(-1.*hit_pstrip->at(i));

			if(hit_sector->at(i)==10 && hit_wlay->at(i)==4) 
				h_eff_sector10_layer4->Fill(-1.*hit_pstrip->at(i));
*/


//			h2_eff_channel_modified->Fill(-1.*hit_pstrip->at(i), fillnumber, 1);

         } else {

		        h_hit_channel->Fill(hit_pstrip->at(i));	
			
		        h_eff_eta->Fill(fillnumber);
		        h_eff_eta_noStuck->Fill(fillnumber);
		        h_eff_eta_good->Fill(fillnumber);

		        h2_eff_channel->Fill(hit_pstrip->at(i), fillnumber, 1);
				h2_eff_channel_onlyNumerator->Fill(hit_pstrip->at(i), fillnumber, 1);

// filling the histograms for 4 layers of a particular sector; just a quick check, don't keep in code 
/*					if(hit_sector->at(i)==10 && hit_wlay->at(i)==1) 
						h_eff_sector10_layer1->Fill(hit_pstrip->at(i));

					if(hit_sector->at(i)==10 && hit_wlay->at(i)==2) 
						h_eff_sector10_layer2->Fill(hit_pstrip->at(i));

					if(hit_sector->at(i)==10 && hit_wlay->at(i)==3) 
						h_eff_sector10_layer3->Fill(hit_pstrip->at(i));

					if(hit_sector->at(i)==10 && hit_wlay->at(i)==4) 
						h_eff_sector10_layer4->Fill(hit_pstrip->at(i));
*/
         
		}

  
	  	 TString histoname = TString::Format("qsum_sector_%d_layer_%d", hit_sector->at(i), hit_wlay->at(i));
//	  	 TString histoname = TString::Format("qsum_layer_%f", hit_sector->at(i)+0.25*(hit_wlay->at(i)-1));
      
         TH1F *myhist = ((TH1F *)(HList->FindObject(histoname)));

         if (!myhist)
         {
//          	myhist = new TH1F(histoname, histoname; "qsum";;, 100, 0., 10000000.);
          	myhist = new TH1F(histoname, histoname, 100, 0., 6000000.);
          	HList->Add(myhist);
         } 
      
      	 myhist->Fill(hit_qsum->at(i));

     } // end loop of hits




      // ===============================================================================================

} // end loop of events


//for(Int_t i=0; i<HList->GetEntries();++i)
//for(Int_t i=0; i<HList->size();++i)
//{
//	(HList[i])->Fit("landau"); 
//}

/*
TH1F *source = (TH1F*)HList->First();
while(source)
{
	source->Draw(); 

	source->Fit("landau");

	gStyle->SetOptFit(111);
	source = (TH1F*)(HList->After(source));
}
*/

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

		for(UInt_t m=0; m<secVectorStuck.size(); ++m)
		{
			if(secVectorStuck[m]==i && layerVectorStuck[m]==j)
				--n_trkpersec_noStuck[i];
		}

		for(UInt_t m=0; m<secVectorAll.size(); ++m)
		{
			if(secVectorAll[m]==i && layerVectorAll[m]==j)
				--n_trkpersec_good[i];
		}


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

         if( h_eff_phi_good->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_phi_good->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec_good[i]-bincontent)/TMath::Power(n_trkpersec_good[i], 3);
            h_eff_phi_good->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec_good[i]);
            h_eff_phi_good->SetBinError((i+1)*5+j+1, binerror);
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

         if( h_eff_eta_good->GetBinContent((i+1)*5+j+1)!=0 ) {
            Double_t bincontent = h_eff_eta_good->GetBinContent((i+1)*5+j+1);
            Double_t binerror = bincontent*(n_trkpersec_good[i]-bincontent)/TMath::Power(n_trkpersec_good[i], 3);
            h_eff_eta_good->SetBinContent((i+1)*5+j+1, bincontent/n_trkpersec_good[i]);
            h_eff_eta_good->SetBinError((i+1)*5+j+1, binerror);
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
