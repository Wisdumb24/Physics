//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 26 13:48:01 2025 by ROOT version 6.32.08
// from TTree trackTree/v1
// found on file: pp_parton_cascade_1.root
//////////////////////////////////////////////////////////

#ifndef trackTree_h
#define trackTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class trackTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *par_pdgid;
   vector<float>   *par_px;
   vector<float>   *par_py;
   vector<float>   *par_pz;
   vector<float>   *par_e;
   vector<float>   *par_x;
   vector<float>   *par_y;
   vector<float>   *par_z;
   vector<float>   *par_t;
   vector<float>   *px;
   vector<float>   *py;
   vector<float>   *pz;
   vector<float>   *m;
   vector<int>     *pid;
   vector<int>     *chg;
   vector<float>   *genJetEta;
   vector<float>   *genJetPt;
   vector<float>   *genJetPhi;
   vector<int>     *genJetChargedMultiplicity;
   vector<vector<int> > *genDau_chg;
   vector<vector<int> > *genDau_pid;
   vector<vector<float> > *genDau_pt;
   vector<vector<float> > *genDau_eta;
   vector<vector<float> > *genDau_phi;

   // List of branches
   TBranch        *b_par_pdgid;   //!
   TBranch        *b_par_px;   //!
   TBranch        *b_par_py;   //!
   TBranch        *b_par_pz;   //!
   TBranch        *b_par_e;   //!
   TBranch        *b_par_x;   //!
   TBranch        *b_par_y;   //!
   TBranch        *b_par_z;   //!
   TBranch        *b_par_t;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_m;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_chg;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetChargedMultiplicity;   //!
   TBranch        *b_genDau_chg;   //!
   TBranch        *b_genDau_pid;   //!
   TBranch        *b_genDau_pt;   //!
   TBranch        *b_genDau_eta;   //!
   TBranch        *b_genDau_phi;   //!

   trackTree(TTree *tree=0);
   virtual ~trackTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef trackTree_cxx
trackTree::trackTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pp_parton_cascade_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pp_parton_cascade_1.root");
      }
      f->GetObject("trackTree",tree);

   }
   Init(tree);
}

trackTree::~trackTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trackTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trackTree::LoadTree(Long64_t entry)
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

void trackTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   par_pdgid = 0;
   par_px = 0;
   par_py = 0;
   par_pz = 0;
   par_e = 0;
   par_x = 0;
   par_y = 0;
   par_z = 0;
   par_t = 0;
   px = 0;
   py = 0;
   pz = 0;
   m = 0;
   pid = 0;
   chg = 0;
   genJetEta = 0;
   genJetPt = 0;
   genJetPhi = 0;
   genJetChargedMultiplicity = 0;
   genDau_chg = 0;
   genDau_pid = 0;
   genDau_pt = 0;
   genDau_eta = 0;
   genDau_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("par_pdgid", &par_pdgid, &b_par_pdgid);
   fChain->SetBranchAddress("par_px", &par_px, &b_par_px);
   fChain->SetBranchAddress("par_py", &par_py, &b_par_py);
   fChain->SetBranchAddress("par_pz", &par_pz, &b_par_pz);
   fChain->SetBranchAddress("par_e", &par_e, &b_par_e);
   fChain->SetBranchAddress("par_x", &par_x, &b_par_x);
   fChain->SetBranchAddress("par_y", &par_y, &b_par_y);
   fChain->SetBranchAddress("par_z", &par_z, &b_par_z);
   fChain->SetBranchAddress("par_t", &par_t, &b_par_t);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("m", &m, &b_m);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("chg", &chg, &b_chg);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetChargedMultiplicity", &genJetChargedMultiplicity, &b_genJetChargedMultiplicity);
   fChain->SetBranchAddress("genDau_chg", &genDau_chg, &b_genDau_chg);
   fChain->SetBranchAddress("genDau_pid", &genDau_pid, &b_genDau_pid);
   fChain->SetBranchAddress("genDau_pt", &genDau_pt, &b_genDau_pt);
   fChain->SetBranchAddress("genDau_eta", &genDau_eta, &b_genDau_eta);
   fChain->SetBranchAddress("genDau_phi", &genDau_phi, &b_genDau_phi);
   Notify();
}

bool trackTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void trackTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trackTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trackTree_cxx
