#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector2.h"
#include "TVector3.h"
#include "coordinateTools.h"
#include <sstream>

using namespace std;

//----------------------------------------------------------------------
// 数据结构定义
//----------------------------------------------------------------------

struct Parton {
    int   pdgid;
    float px, py, pz, e;
    float x, y, z, t;
    float pt, eta, phi;
    float eta_s, tau;
    int   jetID;
    float jet_par_x, jet_par_y, jet_par_z;
};

struct Jet {
    float Pt, Eta, Phi;
    int   genJetChargedMultiplicity;
};

//----------------------------------------------------------------------
// 主程序
//----------------------------------------------------------------------

int pp() {
    TFile* file = TFile::Open("pp_parton_cascade_0.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: cannot open file pp_parton_cascade_0.root" << endl;
        return 1;
    }
    TTree* tree = dynamic_cast<TTree*>(file->Get("trackTree"));
    if (!tree) {
        cerr << "Error: cannot find TTree 'trackTree'" << endl;
        file->Close();
        return 1;
    }

    // Branches
    vector<int>*    b_par_pdgid = nullptr;
    vector<double>* b_par_px    = nullptr;
    vector<double>* b_par_py    = nullptr;
    vector<double>* b_par_pz    = nullptr;
    vector<double>* b_par_e     = nullptr;
    vector<double>* b_par_x     = nullptr;
    vector<double>* b_par_y     = nullptr;
    vector<double>* b_par_z     = nullptr;
    vector<double>* b_par_t     = nullptr;
    tree->SetBranchAddress("par_pdgid", &b_par_pdgid);
    tree->SetBranchAddress("par_px",    &b_par_px);
    tree->SetBranchAddress("par_py",    &b_par_py);
    tree->SetBranchAddress("par_pz",    &b_par_pz);
    tree->SetBranchAddress("par_e",     &b_par_e);
    tree->SetBranchAddress("par_x",     &b_par_x);
    tree->SetBranchAddress("par_y",     &b_par_y);
    tree->SetBranchAddress("par_z",     &b_par_z);
    tree->SetBranchAddress("par_t",     &b_par_t);

    vector<float>* b_genJetPt  = nullptr;
    vector<float>* b_genJetEta = nullptr;
    vector<float>* b_genJetPhi = nullptr;
    vector<int>*   b_genJetChargedMult = nullptr;
    tree->SetBranchAddress("genJetPt",  &b_genJetPt);
    tree->SetBranchAddress("genJetEta", &b_genJetEta);
    tree->SetBranchAddress("genJetPhi", &b_genJetPhi);
    tree->SetBranchAddress("genJetChargedMultiplicity", &b_genJetChargedMult);

    // Histograms
    const float maxEtas = 8.0;
    const int   Nbin    = 8;
    TH1F*  hEtaS = new TH1F("hEtaS","Space-time rapidity #eta_s;#eta_s;dN/d#eta_s",Nbin,0,maxEtas);
    TH2F*  hEcc  = new TH2F("hEcc","Eccentricity vs #eta_s",Nbin,0,maxEtas,100,0,1);
    TH2F*  hTauEta = new TH2F("hTauEta","Proper time #tau vs #eta_s",Nbin,0,maxEtas,100,0,10);
    TH2F*  hMult = new TH2F("hMult","Parton multiplicity vs #eta_s",Nbin,0,maxEtas,100,0,100);
    TH2F*  hJetMultParton = new TH2F("hJetMultParton",";N_{ch}^j;N_p^j",50,0,100,25,0,50);

    vector<int> selectedEvents = {52,138,515,815};
    unordered_map<int,vector<TH1F*>> phiDists;

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    vector<Parton> partons;
    vector<Jet>    jets;
    float xlow = hEcc->GetXaxis()->GetXmin();
    float xup  = hEcc->GetXaxis()->GetXmax();
    float binw = (xup-xlow)/Nbin;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t ie=0; ie<nentries; ++ie) {
        bool isSelected = find(selectedEvents.begin(),selectedEvents.end(),ie)!=selectedEvents.end();
        tree->GetEntry(ie);

        // Load jets
        jets.clear();
        for (size_t j=0;j<b_genJetPt->size();++j)
            jets.push_back({b_genJetPt->at(j),b_genJetEta->at(j),b_genJetPhi->at(j),b_genJetChargedMult->at(j)});

        // Load partons
        partons.clear();
        for (size_t j=0;j<b_par_pdgid->size();++j) {
            Parton p;
            p.pdgid = b_par_pdgid->at(j);
            p.px = b_par_px->at(j); p.py = b_par_py->at(j);
            p.pz = b_par_pz->at(j); p.e  = b_par_e->at(j);
            p.x  = b_par_x->at(j);  p.y  = b_par_y->at(j);
            p.z  = b_par_z->at(j);  p.t  = b_par_t->at(j);
            p.pt  = sqrt(p.px*p.px + p.py*p.py);
            p.phi = atan2(p.py,p.px);
            p.eta = -log(tan(atan2(p.pt,p.pz)/2.0f));
            p.jetID= -1;
            p.jet_par_x=p.jet_par_y=p.jet_par_z=0;
            partons.push_back(p);
        }

        // Associate
        for (auto &p:partons) {
            for (size_t j=0;j<jets.size();++j) {
                float dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
                float dEta = p.eta - jets[j].Eta;
                if (sqrt(dphi*dphi + dEta*dEta)<0.8f) {p.jetID=j; break;}
            }
        }

        // Compute parton multiplicities per jet
        vector<int> partonMult(jets.size(),0);
        for (auto &p:partons)
            if (p.jetID>=0)
                partonMult[p.jetID]++;

        // Output counts using parton multiplicity
        cout<<"Event "<<ie
            <<": partons="<<partons.size()
            <<", jets="<<jets.size()
            <<", first jet parton mult="<<(partonMult.size()>0?partonMult[0]:0)
            <<", second jet parton mult="<<(partonMult.size()>1?partonMult[1]:0)
            <<"\n";

        // Transform coordinates and fill histograms...
        // (rest of code unchanged)
    }

    // Draw/save canvases...
    return 0;
}
