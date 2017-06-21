//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 06:22:04 2017 by ROOT version 6.04/14
// from TTree muon/Muon Track Information
// found on file: /eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root
//////////////////////////////////////////////////////////

#ifndef MuonTreeAnalysis_h
#define MuonTreeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

using namespace std;

class MuonTreeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumiBlockNum;
   vector<int>     *nPV;
   Int_t           bcid;
   Char_t          L1RD0;
   Char_t          L1Mu;
   Int_t           L2Mu;
   Int_t           EFMu;
   vector<unsigned short> *muAuthors;
   vector<float>   *trkPx;
   vector<float>   *trk_me_p;
   vector<float>   *trk_me_eta;
   vector<float>   *trk_me_phi;
   vector<float>   *trk_id_p;
   vector<float>   *trk_id_eta;
   vector<float>   *trk_id_phi;
   vector<float>   *trkP_onTrack;
   vector<float>   *trkPy;
   vector<float>   *trkPz;
   vector<float>   *trkP;
   vector<float>   *trkPt;
   vector<float>   *trkEta;
   vector<float>   *trkPVXQoverP;
   vector<float>   *trkPVXd0;
   vector<float>   *trkPVXZ0;
   vector<float>   *trkPVXTheta;
   vector<float>   *trkPVXPhi0;
   vector<float>   *trkChisq;
   vector<unsigned char> *trkDoF;
   vector<unsigned char> *nBLYHit;
   vector<unsigned char> *nPixelHit;
   vector<unsigned char> *nSCTHit;
   vector<unsigned char> *nTRTHit;
   vector<unsigned char> *nTHTHit;
   vector<unsigned char> *nMDTHit;
   vector<unsigned char> *nTGCPHit;
   vector<unsigned char> *nTGCEHit;
   vector<unsigned char> *nRPCPHit;
   vector<unsigned char> *nRPCEHit;
   vector<unsigned char> *nCSCPHit;
   vector<unsigned char> *nCSCEHit;
   vector<unsigned char> *nMuSeg;
   vector<unsigned char> *hitToMuon;
   vector<unsigned char> *hitToSeg;
   vector<short>   *hit_sector;
   vector<unsigned char> *hit_wlay;
   vector<bool>    *hit_measphi;
   vector<unsigned char> *hit_pstrip;
   vector<float>   *hitRecoTime;
   vector<float>   *hitbeforeT0CorrTime;
   vector<float>   *hitbeforeBPCorrTime;
   vector<float>   *hit_time;
   vector<unsigned char> *hit_sfit;
   vector<unsigned char> *hit_tfit;
   vector<float>   *hit_pos;
   vector<float>   *hit_PhiPosReFit;
   vector<float>   *hit_dpos;
   vector<unsigned char> *hit_nstr;
   vector<unsigned char> *hit_str0;
   vector<float>   *hit_qsum;
   vector<float>   *hit_qpeak;
   vector<float>   *hit_qleft;
   vector<float>   *hit_qright;
   vector<float>   *hit_dqpeak;
   vector<float>   *hit_dqleft;
   vector<float>   *hit_dqright;
   vector<bool>    *hit_phase;
   vector<float>   *hit_posrefit;
   vector<float>   *hit_phiposrefit;
   vector<float>   *hit_dposrefit;
   vector<float>   *hit_qfitsig;
   vector<float>   *hit_qfitdiff;
   vector<unsigned char> *hit_srefit;
   vector<unsigned char> *segToMuon;
   vector<float>   *seg_time;
   vector<float>   *seg_chsq;
   vector<float>   *seg_ndof;
   vector<float>   *seg_gposx;
   vector<float>   *seg_gposy;
   vector<float>   *seg_gposz;
   vector<float>   *seg_gdirx;
   vector<float>   *seg_gdiry;
   vector<float>   *seg_gdirz;
   vector<float>   *seg_posTheta;
   vector<float>   *seg_posEta;
   vector<float>   *seg_posPhi;
   vector<float>   *seg_dirTheta;
   vector<float>   *seg_dirEta;
   vector<float>   *seg_dirPhi;
   vector<float>   *seg_etalocpos;
   vector<float>   *seg_philocpos;
   vector<float>   *seg_etalocdir;
   vector<float>   *seg_philocdir;
   vector<float>   *seg_dy;
   vector<float>   *seg_dz;
   vector<float>   *seg_day;
   vector<float>   *seg_daz;
   vector<float>   *seg_eyz;
   vector<float>   *seg_eyay;
   vector<float>   *seg_eyaz;
   vector<float>   *seg_ezay;
   vector<float>   *seg_ezaz;
   vector<float>   *seg_eayaz;
   vector<short>   *seg_sector;
   vector<vector<float> > *seg_posesX;
   vector<vector<float> > *seg_poses;
   vector<vector<float> > *seg_posesZ;
   vector<vector<float> > *seg_dposes;
   vector<vector<int> > *seg_pstrips;
   vector<vector<float> > *seg_charges;
   vector<vector<float> > *seg_times;
   vector<vector<int> > *seg_sfits;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumiBlockNum;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_L1RD0;   //!
   TBranch        *b_L1Mu;   //!
   TBranch        *b_L2Mu;   //!
   TBranch        *b_EFMu;   //!
   TBranch        *b_muAuthors;   //!
   TBranch        *b_trkPx;   //!
   TBranch        *b_trk_me_p;   //!
   TBranch        *b_trk_me_eta;   //!
   TBranch        *b_trk_me_phi;   //!
   TBranch        *b_trk_id_p;   //!
   TBranch        *b_trk_id_eta;   //!
   TBranch        *b_trk_id_phi;   //!
   TBranch        *b_trkP_onTrack;   //!
   TBranch        *b_trkPy;   //!
   TBranch        *b_trkPz;   //!
   TBranch        *b_trkP;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPVXQoverP;   //!
   TBranch        *b_trkPVXd0;   //!
   TBranch        *b_trkPVXZ0;   //!
   TBranch        *b_trkPVXTheta;   //!
   TBranch        *b_trkPVXPhi0;   //!
   TBranch        *b_trkChisq;   //!
   TBranch        *b_trkDoF;   //!
   TBranch        *b_nBLYHit;   //!
   TBranch        *b_nPixelHit;   //!
   TBranch        *b_nSCTHit;   //!
   TBranch        *b_nTRTHit;   //!
   TBranch        *b_nTHTHit;   //!
   TBranch        *b_nMDTHit;   //!
   TBranch        *b_nTGCPHit;   //!
   TBranch        *b_nTGCEHit;   //!
   TBranch        *b_nRPCPHit;   //!
   TBranch        *b_nRPCEHit;   //!
   TBranch        *b_nCSCPHit;   //!
   TBranch        *b_nCSCEHit;   //!
   TBranch        *b_nMuSeg;   //!
   TBranch        *b_hitToMuon;   //!
   TBranch        *b_hitToSeg;   //!
   TBranch        *b_hit_sector;   //!
   TBranch        *b_hit_wlay;   //!
   TBranch        *b_hit_measphi;   //!
   TBranch        *b_hit_pstrip;   //!
   TBranch        *b_hitRecoTime;   //!
   TBranch        *b_hitbeforeT0CorrTime;   //!
   TBranch        *b_hitbeforeBPCorrTime;   //!
   TBranch        *b_hit_time;   //!
   TBranch        *b_hit_sfit;   //!
   TBranch        *b_hit_tfit;   //!
   TBranch        *b_hit_pos;   //!
   TBranch        *b_hit_PhiPosReFit;   //!
   TBranch        *b_hit_dpos;   //!
   TBranch        *b_hit_nstr;   //!
   TBranch        *b_hit_str0;   //!
   TBranch        *b_hit_qsum;   //!
   TBranch        *b_hit_qpeak;   //!
   TBranch        *b_hit_qleft;   //!
   TBranch        *b_hit_qright;   //!
   TBranch        *b_hit_dqpeak;   //!
   TBranch        *b_hit_dqleft;   //!
   TBranch        *b_hit_dqright;   //!
   TBranch        *b_hit_phase;   //!
   TBranch        *b_hit_posrefit;   //!
   TBranch        *b_hit_phiposrefit;   //!
   TBranch        *b_hit_dposrefit;   //!
   TBranch        *b_hit_qfitsig;   //!
   TBranch        *b_hit_qfitdiff;   //!
   TBranch        *b_hit_srefit;   //!
   TBranch        *b_segToMuon;   //!
   TBranch        *b_seg_time;   //!
   TBranch        *b_seg_chsq;   //!
   TBranch        *b_seg_ndof;   //!
   TBranch        *b_seg_gposx;   //!
   TBranch        *b_seg_gposy;   //!
   TBranch        *b_seg_gposz;   //!
   TBranch        *b_seg_gdirx;   //!
   TBranch        *b_seg_gdiry;   //!
   TBranch        *b_seg_gdirz;   //!
   TBranch        *b_seg_posTheta;   //!
   TBranch        *b_seg_posEta;   //!
   TBranch        *b_seg_posPhi;   //!
   TBranch        *b_seg_dirTheta;   //!
   TBranch        *b_seg_dirEta;   //!
   TBranch        *b_seg_dirPhi;   //!
   TBranch        *b_seg_etalocpos;   //!
   TBranch        *b_seg_philocpos;   //!
   TBranch        *b_seg_etalocdir;   //!
   TBranch        *b_seg_philocdir;   //!
   TBranch        *b_seg_dy;   //!
   TBranch        *b_seg_dz;   //!
   TBranch        *b_seg_day;   //!
   TBranch        *b_seg_daz;   //!
   TBranch        *b_seg_eyz;   //!
   TBranch        *b_seg_eyay;   //!
   TBranch        *b_seg_eyaz;   //!
   TBranch        *b_seg_ezay;   //!
   TBranch        *b_seg_ezaz;   //!
   TBranch        *b_seg_eayaz;   //!
   TBranch        *b_seg_sector;   //!
   TBranch        *b_seg_posesX;   //!
   TBranch        *b_seg_poses;   //!
   TBranch        *b_seg_posesZ;   //!
   TBranch        *b_seg_dposes;   //!
   TBranch        *b_seg_pstrips;   //!
   TBranch        *b_seg_charges;   //!
   TBranch        *b_seg_times;   //!
   TBranch        *b_seg_sfits;   //!

   MuonTreeAnalysis(TTree *tree=0);
   virtual ~MuonTreeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MuonTreeAnalysis_cxx
MuonTreeAnalysis::MuonTreeAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root");
      }
      f->GetObject("muon",tree);

   }
   Init(tree);
}

MuonTreeAnalysis::~MuonTreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MuonTreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MuonTreeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MuonTreeAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   nPV = 0;
   muAuthors = 0;
   trkPx = 0;
   trk_me_p = 0;
   trk_me_eta = 0;
   trk_me_phi = 0;
   trk_id_p = 0;
   trk_id_eta = 0;
   trk_id_phi = 0;
   trkP_onTrack = 0;
   trkPy = 0;
   trkPz = 0;
   trkP = 0;
   trkPt = 0;
   trkEta = 0;
   trkPVXQoverP = 0;
   trkPVXd0 = 0;
   trkPVXZ0 = 0;
   trkPVXTheta = 0;
   trkPVXPhi0 = 0;
   trkChisq = 0;
   trkDoF = 0;
   nBLYHit = 0;
   nPixelHit = 0;
   nSCTHit = 0;
   nTRTHit = 0;
   nTHTHit = 0;
   nMDTHit = 0;
   nTGCPHit = 0;
   nTGCEHit = 0;
   nRPCPHit = 0;
   nRPCEHit = 0;
   nCSCPHit = 0;
   nCSCEHit = 0;
   nMuSeg = 0;
   hitToMuon = 0;
   hitToSeg = 0;
   hit_sector = 0;
   hit_wlay = 0;
   hit_measphi = 0;
   hit_pstrip = 0;
   hitRecoTime = 0;
   hitbeforeT0CorrTime = 0;
   hitbeforeBPCorrTime = 0;
   hit_time = 0;
   hit_sfit = 0;
   hit_tfit = 0;
   hit_pos = 0;
   hit_PhiPosReFit = 0;
   hit_dpos = 0;
   hit_nstr = 0;
   hit_str0 = 0;
   hit_qsum = 0;
   hit_qpeak = 0;
   hit_qleft = 0;
   hit_qright = 0;
   hit_dqpeak = 0;
   hit_dqleft = 0;
   hit_dqright = 0;
   hit_phase = 0;
   hit_posrefit = 0;
   hit_phiposrefit = 0;
   hit_dposrefit = 0;
   hit_qfitsig = 0;
   hit_qfitdiff = 0;
   hit_srefit = 0;
   segToMuon = 0;
   seg_time = 0;
   seg_chsq = 0;
   seg_ndof = 0;
   seg_gposx = 0;
   seg_gposy = 0;
   seg_gposz = 0;
   seg_gdirx = 0;
   seg_gdiry = 0;
   seg_gdirz = 0;
   seg_posTheta = 0;
   seg_posEta = 0;
   seg_posPhi = 0;
   seg_dirTheta = 0;
   seg_dirEta = 0;
   seg_dirPhi = 0;
   seg_etalocpos = 0;
   seg_philocpos = 0;
   seg_etalocdir = 0;
   seg_philocdir = 0;
   seg_dy = 0;
   seg_dz = 0;
   seg_day = 0;
   seg_daz = 0;
   seg_eyz = 0;
   seg_eyay = 0;
   seg_eyaz = 0;
   seg_ezay = 0;
   seg_ezaz = 0;
   seg_eayaz = 0;
   seg_sector = 0;
   seg_posesX = 0;
   seg_poses = 0;
   seg_posesZ = 0;
   seg_dposes = 0;
   seg_pstrips = 0;
   seg_charges = 0;
   seg_times = 0;
   seg_sfits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumiBlockNum", &lumiBlockNum, &b_lumiBlockNum);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("L1RD0", &L1RD0, &b_L1RD0);
   fChain->SetBranchAddress("L1Mu", &L1Mu, &b_L1Mu);
   fChain->SetBranchAddress("L2Mu", &L2Mu, &b_L2Mu);
   fChain->SetBranchAddress("EFMu", &EFMu, &b_EFMu);
   fChain->SetBranchAddress("muAuthors", &muAuthors, &b_muAuthors);
   fChain->SetBranchAddress("trkPx", &trkPx, &b_trkPx);
   fChain->SetBranchAddress("trk_me_p", &trk_me_p, &b_trk_me_p);
   fChain->SetBranchAddress("trk_me_eta", &trk_me_eta, &b_trk_me_eta);
   fChain->SetBranchAddress("trk_me_phi", &trk_me_phi, &b_trk_me_phi);
   fChain->SetBranchAddress("trk_id_p", &trk_id_p, &b_trk_id_p);
   fChain->SetBranchAddress("trk_id_eta", &trk_id_eta, &b_trk_id_eta);
   fChain->SetBranchAddress("trk_id_phi", &trk_id_phi, &b_trk_id_phi);
   fChain->SetBranchAddress("trkP_onTrack", &trkP_onTrack, &b_trkP_onTrack);
   fChain->SetBranchAddress("trkPy", &trkPy, &b_trkPy);
   fChain->SetBranchAddress("trkPz", &trkPz, &b_trkPz);
   fChain->SetBranchAddress("trkP", &trkP, &b_trkP);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPVXQoverP", &trkPVXQoverP, &b_trkPVXQoverP);
   fChain->SetBranchAddress("trkPVXd0", &trkPVXd0, &b_trkPVXd0);
   fChain->SetBranchAddress("trkPVXZ0", &trkPVXZ0, &b_trkPVXZ0);
   fChain->SetBranchAddress("trkPVXTheta", &trkPVXTheta, &b_trkPVXTheta);
   fChain->SetBranchAddress("trkPVXPhi0", &trkPVXPhi0, &b_trkPVXPhi0);
   fChain->SetBranchAddress("trkChisq", &trkChisq, &b_trkChisq);
   fChain->SetBranchAddress("trkDoF", &trkDoF, &b_trkDoF);
   fChain->SetBranchAddress("nBLYHit", &nBLYHit, &b_nBLYHit);
   fChain->SetBranchAddress("nPixelHit", &nPixelHit, &b_nPixelHit);
   fChain->SetBranchAddress("nSCTHit", &nSCTHit, &b_nSCTHit);
   fChain->SetBranchAddress("nTRTHit", &nTRTHit, &b_nTRTHit);
   fChain->SetBranchAddress("nTHTHit", &nTHTHit, &b_nTHTHit);
   fChain->SetBranchAddress("nMDTHit", &nMDTHit, &b_nMDTHit);
   fChain->SetBranchAddress("nTGCPHit", &nTGCPHit, &b_nTGCPHit);
   fChain->SetBranchAddress("nTGCEHit", &nTGCEHit, &b_nTGCEHit);
   fChain->SetBranchAddress("nRPCPHit", &nRPCPHit, &b_nRPCPHit);
   fChain->SetBranchAddress("nRPCEHit", &nRPCEHit, &b_nRPCEHit);
   fChain->SetBranchAddress("nCSCPHit", &nCSCPHit, &b_nCSCPHit);
   fChain->SetBranchAddress("nCSCEHit", &nCSCEHit, &b_nCSCEHit);
   fChain->SetBranchAddress("nMuSeg", &nMuSeg, &b_nMuSeg);
   fChain->SetBranchAddress("hitToMuon", &hitToMuon, &b_hitToMuon);
   fChain->SetBranchAddress("hitToSeg", &hitToSeg, &b_hitToSeg);
   fChain->SetBranchAddress("hit_sector", &hit_sector, &b_hit_sector);
   fChain->SetBranchAddress("hit_wlay", &hit_wlay, &b_hit_wlay);
   fChain->SetBranchAddress("hit_measphi", &hit_measphi, &b_hit_measphi);
   fChain->SetBranchAddress("hit_pstrip", &hit_pstrip, &b_hit_pstrip);
   fChain->SetBranchAddress("hitRecoTime", &hitRecoTime, &b_hitRecoTime);
   fChain->SetBranchAddress("hitbeforeT0CorrTime", &hitbeforeT0CorrTime, &b_hitbeforeT0CorrTime);
   fChain->SetBranchAddress("hitbeforeBPCorrTime", &hitbeforeBPCorrTime, &b_hitbeforeBPCorrTime);
   fChain->SetBranchAddress("hit_time", &hit_time, &b_hit_time);
   fChain->SetBranchAddress("hit_sfit", &hit_sfit, &b_hit_sfit);
   fChain->SetBranchAddress("hit_tfit", &hit_tfit, &b_hit_tfit);
   fChain->SetBranchAddress("hit_pos", &hit_pos, &b_hit_pos);
   fChain->SetBranchAddress("hit_PhiPosReFit", &hit_PhiPosReFit, &b_hit_PhiPosReFit);
   fChain->SetBranchAddress("hit_dpos", &hit_dpos, &b_hit_dpos);
   fChain->SetBranchAddress("hit_nstr", &hit_nstr, &b_hit_nstr);
   fChain->SetBranchAddress("hit_str0", &hit_str0, &b_hit_str0);
   fChain->SetBranchAddress("hit_qsum", &hit_qsum, &b_hit_qsum);
   fChain->SetBranchAddress("hit_qpeak", &hit_qpeak, &b_hit_qpeak);
   fChain->SetBranchAddress("hit_qleft", &hit_qleft, &b_hit_qleft);
   fChain->SetBranchAddress("hit_qright", &hit_qright, &b_hit_qright);
   fChain->SetBranchAddress("hit_dqpeak", &hit_dqpeak, &b_hit_dqpeak);
   fChain->SetBranchAddress("hit_dqleft", &hit_dqleft, &b_hit_dqleft);
   fChain->SetBranchAddress("hit_dqright", &hit_dqright, &b_hit_dqright);
   fChain->SetBranchAddress("hit_phase", &hit_phase, &b_hit_phase);
   fChain->SetBranchAddress("hit_posrefit", &hit_posrefit, &b_hit_posrefit);
   fChain->SetBranchAddress("hit_phiposrefit", &hit_phiposrefit, &b_hit_phiposrefit);
   fChain->SetBranchAddress("hit_dposrefit", &hit_dposrefit, &b_hit_dposrefit);
   fChain->SetBranchAddress("hit_qfitsig", &hit_qfitsig, &b_hit_qfitsig);
   fChain->SetBranchAddress("hit_qfitdiff", &hit_qfitdiff, &b_hit_qfitdiff);
   fChain->SetBranchAddress("hit_srefit", &hit_srefit, &b_hit_srefit);
   fChain->SetBranchAddress("segToMuon", &segToMuon, &b_segToMuon);
   fChain->SetBranchAddress("seg_time", &seg_time, &b_seg_time);
   fChain->SetBranchAddress("seg_chsq", &seg_chsq, &b_seg_chsq);
   fChain->SetBranchAddress("seg_ndof", &seg_ndof, &b_seg_ndof);
   fChain->SetBranchAddress("seg_gposx", &seg_gposx, &b_seg_gposx);
   fChain->SetBranchAddress("seg_gposy", &seg_gposy, &b_seg_gposy);
   fChain->SetBranchAddress("seg_gposz", &seg_gposz, &b_seg_gposz);
   fChain->SetBranchAddress("seg_gdirx", &seg_gdirx, &b_seg_gdirx);
   fChain->SetBranchAddress("seg_gdiry", &seg_gdiry, &b_seg_gdiry);
   fChain->SetBranchAddress("seg_gdirz", &seg_gdirz, &b_seg_gdirz);
   fChain->SetBranchAddress("seg_posTheta", &seg_posTheta, &b_seg_posTheta);
   fChain->SetBranchAddress("seg_posEta", &seg_posEta, &b_seg_posEta);
   fChain->SetBranchAddress("seg_posPhi", &seg_posPhi, &b_seg_posPhi);
   fChain->SetBranchAddress("seg_dirTheta", &seg_dirTheta, &b_seg_dirTheta);
   fChain->SetBranchAddress("seg_dirEta", &seg_dirEta, &b_seg_dirEta);
   fChain->SetBranchAddress("seg_dirPhi", &seg_dirPhi, &b_seg_dirPhi);
   fChain->SetBranchAddress("seg_etalocpos", &seg_etalocpos, &b_seg_etalocpos);
   fChain->SetBranchAddress("seg_philocpos", &seg_philocpos, &b_seg_philocpos);
   fChain->SetBranchAddress("seg_etalocdir", &seg_etalocdir, &b_seg_etalocdir);
   fChain->SetBranchAddress("seg_philocdir", &seg_philocdir, &b_seg_philocdir);
   fChain->SetBranchAddress("seg_dy", &seg_dy, &b_seg_dy);
   fChain->SetBranchAddress("seg_dz", &seg_dz, &b_seg_dz);
   fChain->SetBranchAddress("seg_day", &seg_day, &b_seg_day);
   fChain->SetBranchAddress("seg_daz", &seg_daz, &b_seg_daz);
   fChain->SetBranchAddress("seg_eyz", &seg_eyz, &b_seg_eyz);
   fChain->SetBranchAddress("seg_eyay", &seg_eyay, &b_seg_eyay);
   fChain->SetBranchAddress("seg_eyaz", &seg_eyaz, &b_seg_eyaz);
   fChain->SetBranchAddress("seg_ezay", &seg_ezay, &b_seg_ezay);
   fChain->SetBranchAddress("seg_ezaz", &seg_ezaz, &b_seg_ezaz);
   fChain->SetBranchAddress("seg_eayaz", &seg_eayaz, &b_seg_eayaz);
   fChain->SetBranchAddress("seg_sector", &seg_sector, &b_seg_sector);
   fChain->SetBranchAddress("seg_posesX", &seg_posesX, &b_seg_posesX);
   fChain->SetBranchAddress("seg_poses", &seg_poses, &b_seg_poses);
   fChain->SetBranchAddress("seg_posesZ", &seg_posesZ, &b_seg_posesZ);
   fChain->SetBranchAddress("seg_dposes", &seg_dposes, &b_seg_dposes);
   fChain->SetBranchAddress("seg_pstrips", &seg_pstrips, &b_seg_pstrips);
   fChain->SetBranchAddress("seg_charges", &seg_charges, &b_seg_charges);
   fChain->SetBranchAddress("seg_times", &seg_times, &b_seg_times);
   fChain->SetBranchAddress("seg_sfits", &seg_sfits, &b_seg_sfits);
   Notify();
}

Bool_t MuonTreeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MuonTreeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MuonTreeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MuonTreeAnalysis_cxx
