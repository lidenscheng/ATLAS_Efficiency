//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 18 06:25:01 2017 by ROOT version 6.04/14
// from TTree csc_cluster/list of any cluster
// found on file: /eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root
//////////////////////////////////////////////////////////

#ifndef ClusterTreeAnalysis_h
#define ClusterTreeAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

using namespace std;

class ClusterTreeAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumiBlockNum;
   Int_t           bcid;
   vector<short>   *cluToMuon;
   vector<float>   *x;
   vector<float>   *y;
   vector<float>   *z;
   vector<float>   *r;
   vector<bool>    *measphi;
   vector<unsigned char> *wlay;
   vector<short>   *sector;
   vector<float>   *pos;
   vector<float>   *error;
   vector<unsigned char> *sfit;
   vector<unsigned char> *tfit;
   vector<float>   *time;
   vector<float>   *qsum;
   vector<unsigned char> *pstrip;
   vector<unsigned char> *nstrip;
   vector<unsigned char> *strip0;
   vector<float>   *qpeak;
   vector<float>   *qleft;
   vector<float>   *qright;
   vector<float>   *dqpeak;
   vector<float>   *dqlef;
   vector<float>   *dqright;
   vector<float>   *posrefit;
   vector<float>   *errrefit;
   vector<float>   *qfitdiff;
   vector<float>   *qfitsig;
   vector<unsigned char> *srefit;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumiBlockNum;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_cluToMuon;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_r;   //!
   TBranch        *b_measphi;   //!
   TBranch        *b_wlay;   //!
   TBranch        *b_sector;   //!
   TBranch        *b_pos;   //!
   TBranch        *b_error;   //!
   TBranch        *b_sfit;   //!
   TBranch        *b_tfit;   //!
   TBranch        *b_time;   //!
   TBranch        *b_qsum;   //!
   TBranch        *b_pstrip;   //!
   TBranch        *b_nstrip;   //!
   TBranch        *b_strip0;   //!
   TBranch        *b_qpeak;   //!
   TBranch        *b_qleft;   //!
   TBranch        *b_qright;   //!
   TBranch        *b_dqpeak;   //!
   TBranch        *b_dqlef;   //!
   TBranch        *b_dqright;   //!
   TBranch        *b_posrefit;   //!
   TBranch        *b_errrefit;   //!
   TBranch        *b_qfitdiff;   //!
   TBranch        *b_qfitsig;   //!
   TBranch        *b_srefit;   //!

   ClusterTreeAnalysis(TTree *tree=0);
   virtual ~ClusterTreeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ClusterTreeAnalysis_cxx
ClusterTreeAnalysis::ClusterTreeAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/atlas/user/e/ekarentz/ATLAS_RUN_CSC/run00302393_DESDM_MCP_stableBeam_merged_extValidNtuple.root");
      }
      f->GetObject("csc_cluster",tree);

   }
   Init(tree);
}

ClusterTreeAnalysis::~ClusterTreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ClusterTreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ClusterTreeAnalysis::LoadTree(Long64_t entry)
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

void ClusterTreeAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   cluToMuon = 0;
   x = 0;
   y = 0;
   z = 0;
   r = 0;
   measphi = 0;
   wlay = 0;
   sector = 0;
   pos = 0;
   error = 0;
   sfit = 0;
   tfit = 0;
   time = 0;
   qsum = 0;
   pstrip = 0;
   nstrip = 0;
   strip0 = 0;
   qpeak = 0;
   qleft = 0;
   qright = 0;
   dqpeak = 0;
   dqlef = 0;
   dqright = 0;
   posrefit = 0;
   errrefit = 0;
   qfitdiff = 0;
   qfitsig = 0;
   srefit = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumiBlockNum", &lumiBlockNum, &b_lumiBlockNum);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("cluToMuon", &cluToMuon, &b_cluToMuon);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("r", &r, &b_r);
   fChain->SetBranchAddress("measphi", &measphi, &b_measphi);
   fChain->SetBranchAddress("wlay", &wlay, &b_wlay);
   fChain->SetBranchAddress("sector", &sector, &b_sector);
   fChain->SetBranchAddress("pos", &pos, &b_pos);
   fChain->SetBranchAddress("error", &error, &b_error);
   fChain->SetBranchAddress("sfit", &sfit, &b_sfit);
   fChain->SetBranchAddress("tfit", &tfit, &b_tfit);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("qsum", &qsum, &b_qsum);
   fChain->SetBranchAddress("pstrip", &pstrip, &b_pstrip);
   fChain->SetBranchAddress("nstrip", &nstrip, &b_nstrip);
   fChain->SetBranchAddress("strip0", &strip0, &b_strip0);
   fChain->SetBranchAddress("qpeak", &qpeak, &b_qpeak);
   fChain->SetBranchAddress("qleft", &qleft, &b_qleft);
   fChain->SetBranchAddress("qright", &qright, &b_qright);
   fChain->SetBranchAddress("dqpeak", &dqpeak, &b_dqpeak);
   fChain->SetBranchAddress("dqlef", &dqlef, &b_dqlef);
   fChain->SetBranchAddress("dqright", &dqright, &b_dqright);
   fChain->SetBranchAddress("posrefit", &posrefit, &b_posrefit);
   fChain->SetBranchAddress("errrefit", &errrefit, &b_errrefit);
   fChain->SetBranchAddress("qfitdiff", &qfitdiff, &b_qfitdiff);
   fChain->SetBranchAddress("qfitsig", &qfitsig, &b_qfitsig);
   fChain->SetBranchAddress("srefit", &srefit, &b_srefit);
   Notify();
}

Bool_t ClusterTreeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ClusterTreeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ClusterTreeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ClusterTreeAnalysis_cxx
