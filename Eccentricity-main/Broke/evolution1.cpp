#include <vector>
#include <unordered_map>
#include <memory>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector2.h"
#include "coordinateTools.h"

struct Parton {
    int    pdgid;
    float  px, py, pz, e;
    float  x, y, z, t;
    float  pt, eta, phi;
    int    jetID = -1;
    float  jet_par_x = 0, jet_par_y = 0, jet_par_z = 0;
    float  jet_par_px = 0, jet_par_py = 0, jet_par_pz = 0;
};

struct Jet {
    float Pt, Eta, Phi;
    int   ChargedMult;
};

bool OpenTree(const std::string& path,
              std::unique_ptr<TFile>& file,
              TTree*& tree)
{
    file.reset(TFile::Open(path.c_str(), "READ"));
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open file " << path << "\n";
        return false;
    }
    tree = dynamic_cast<TTree*>(file->Get("trackTree"));
    if (!tree) {
        std::cerr << "Error: cannot find TTree 'trackTree'\n";
        return false;
    }
    return true;
}

void SetBranchAddresses(TTree* tree,
    std::vector<int>*&    b_pdgid,
    std::vector<double>*& b_px,
    std::vector<double>*& b_py,
    std::vector<double>*& b_pz,
    std::vector<double>*& b_e,
    std::vector<double>*& b_x,
    std::vector<double>*& b_y,
    std::vector<double>*& b_z,
    std::vector<double>*& b_t,
    std::vector<float>*&  b_gPt,
    std::vector<float>*&  b_gEta,
    std::vector<float>*&  b_gPhi,
    std::vector<int>*&    b_gMult)
{
    tree->SetBranchAddress("par_pdgid", &b_pdgid);
    tree->SetBranchAddress("par_px",    &b_px);
    tree->SetBranchAddress("par_py",    &b_py);
    tree->SetBranchAddress("par_pz",    &b_pz);
    tree->SetBranchAddress("par_e",     &b_e);
    tree->SetBranchAddress("par_x",     &b_x);
    tree->SetBranchAddress("par_y",     &b_y);
    tree->SetBranchAddress("par_z",     &b_z);
    tree->SetBranchAddress("par_t",     &b_t);

    tree->SetBranchAddress("genJetPt",                  &b_gPt);
    tree->SetBranchAddress("genJetEta",                 &b_gEta);
    tree->SetBranchAddress("genJetPhi",                 &b_gPhi);
    tree->SetBranchAddress("genJetChargedMultiplicity", &b_gMult);
}

std::vector<Jet> ExtractJets(const std::vector<float>* bPt,
                             const std::vector<float>* bEta,
                             const std::vector<float>* bPhi,
                             const std::vector<int>*   bMult)
{
    std::vector<Jet> jets;
    jets.reserve(bPt->size());
    for (size_t i = 0; i < bPt->size(); ++i) {
        jets.push_back(Jet{ (*bPt)[i], (*bEta)[i], (*bPhi)[i], (*bMult)[i] });
    }
    return jets;
}

std::vector<Parton> ExtractPartons(
    const std::vector<int>*    b_pdgid,
    const std::vector<double>* b_px,
    const std::vector<double>* b_py,
    const std::vector<double>* b_pz,
    const std::vector<double>* b_e,
    const std::vector<double>* b_x,
    const std::vector<double>* b_y,
    const std::vector<double>* b_z,
    const std::vector<double>* b_t)
{
    size_t n = b_pdgid->size();
    std::vector<Parton> parts;
    parts.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        Parton p;
        p.pdgid = (*b_pdgid)[i]; p.px = (*b_px)[i]; p.py = (*b_py)[i]; p.pz = (*b_pz)[i]; p.e = (*b_e)[i];
        p.x = (*b_x)[i]; p.y = (*b_y)[i]; p.z = (*b_z)[i]; p.t = (*b_t)[i];
        p.pt = std::hypot(p.px, p.py);
        p.phi = std::atan2(p.py, p.px);
        float theta = std::atan2(p.pt, p.pz);
        p.eta = -std::log(std::tan(theta/2.0f));
        parts.push_back(p);
    }
    return parts;
}

void MatchPartonsToJets(std::vector<Parton>& parts,
                        const std::vector<Jet>& jets,
                        float dRcut = 0.8f)
{
    for (auto& p : parts) {
        for (size_t j = 0; j < jets.size(); ++j) {
            float dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
            float dEta = p.eta - jets[j].Eta;
            if (std::hypot(dEta, dphi) < dRcut) {
                p.jetID = int(j);
                break;
            }
        }
    }
}

void TransformToJetFrame(std::vector<Parton>& parts,
                         const std::vector<Jet>& jets)
{
    for (auto& p : parts) {
        if (p.jetID < 0) continue;
        TVector3 jv; jv.SetPtEtaPhi(jets[p.jetID].Pt, jets[p.jetID].Eta, jets[p.jetID].Phi);
        TVector3 pos(p.x, p.y, p.z);
        float ptJ = ptWRTJet(jv, pos);
        float thJ = thetaWRTJet(jv, pos);
        float phJ = phiWRTJet(jv, pos);
        p.jet_par_x  = ptJ * std::cos(phJ);
        p.jet_par_y  = ptJ * std::sin(phJ);
        p.jet_par_z  = ptJ / std::tan(thJ);
        TVector3 mom(p.px, p.py, p.pz);
        float ptJm = ptWRTJet(jv, mom);
        float thJm = thetaWRTJet(jv, mom);
        float phJm = phiWRTJet(jv, mom);
        p.jet_par_px = ptJm * std::cos(phJm);
        p.jet_par_py = ptJm * std::sin(phJm);
        p.jet_par_pz = ptJm / std::tan(phJm);
    }
}

// FillPsi2AtTime: uses position-based phi (phiLoc) to compute ψ2,
// and momentum-based phi (phiMom) for φ–ψ₂ distributions
void FillPsi2AtTime(const std::vector<Parton>& parts,
                    double tTarget,
                    TH1F* hPsi2,
                    TH1F* hPhiRel,
                    TH2F* hPsi2VsMult,
                    TH2F* hPhiRelVsMult)
{
    std::unordered_map<int, std::vector<int>> byJet;
    byJet.reserve(16);
    for (size_t i = 0; i < parts.size(); ++i) {
        if (parts[i].jetID >= 0) byJet[parts[i].jetID].push_back(int(i));
    }
    for (auto const& kv : byJet) {
        auto const& idxs = kv.second;
        if (idxs.size() < 2) continue;
        double sumE=0, sumx=0, sumy=0;
        for (int i : idxs) {
            auto const& P = parts[i];
            if (P.t > tTarget) continue;
            double dt = tTarget - P.t;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            sumE += P.e; sumx += P.e * x; sumy += P.e * y;
        }
        if (sumE == 0) continue;
        double xm = sumx / sumE, ym = sumy / sumE;
        double Re=0, Im=0;
        // position-based phi for Q-vector
        for (int i : idxs) {
            auto const& P = parts[i];
            if (P.t > tTarget) continue;
            double dt = tTarget - P.t;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            double r2 = (x-xm)*(x-xm) + (y-ym)*(y-ym);
            double phiLoc = std::atan2(y-ym, x-xm);
            double w = P.e * r2;
            Re += w * std::cos(2*phiLoc);
            Im += w * std::sin(2*phiLoc);
        }
        if (Re==0 && Im==0) continue;
        double psi2 = 0.5 * std::atan2(Im, Re);
        int mult = int(idxs.size());
        hPsi2->Fill(psi2);
        hPsi2VsMult->Fill(psi2, mult);
        // momentum-based phi for phi-psi
        for (int i : idxs) {
            auto const& P = parts[i];
            if (P.t > tTarget) continue;
            double phiMom = std::atan2(P.jet_par_py, P.jet_par_px);
            double dphi   = TVector2::Phi_mpi_pi(phiMom - psi2);
            hPhiRel->Fill(dphi);
            hPhiRelVsMult->Fill(dphi, mult);
        }
    }
}

int evolution1() {
    constexpr double tTarget = 2.0;  // fm/c
    std::unique_ptr<TFile> file;
    TTree* tree = nullptr;
    if (!OpenTree("pp_parton_cascade_1.root", file, tree)) return 1;
    std::vector<int>*    b_pdgid = nullptr;
    std::vector<double>* b_px    = nullptr;
    std::vector<double>* b_py    = nullptr;
    std::vector<double>* b_pz    = nullptr;
    std::vector<double>* b_e     = nullptr;
    std::vector<double>* b_x     = nullptr;
    std::vector<double>* b_y     = nullptr;
    std::vector<double>* b_z     = nullptr;
    std::vector<double>* b_t     = nullptr;
    std::vector<float>*  b_gPt   = nullptr;
    std::vector<float>*  b_gEta  = nullptr;
    std::vector<float>*  b_gPhi  = nullptr;
    std::vector<int>*    b_gMult = nullptr;
    SetBranchAddresses(tree,
        b_pdgid, b_px, b_py, b_pz, b_e,
        b_x, b_y, b_z, b_t,
        b_gPt, b_gEta, b_gPhi, b_gMult);
    TH1F* hPsi2         = new TH1F("hPsi2",   "Psi2 @ 2fm/c;#Psi_{2};Counts", 50, -M_PI/2, M_PI/2);
    TH1F* hPhiRel       = new TH1F("hPhiRel", "phi-Psi2 @ 2fm/c;#Delta#phi;Counts", 100, -M_PI, M_PI);
    TH2F* hPsi2VsMult   = new TH2F("hPsi2VsMult",   "#Psi_{2} vs N @ 2fm/c;#Psi_{2};N_{partons}", 50, -M_PI/2, M_PI/2, 50, 0, 50);
    TH2F* hPhiRelVsMult = new TH2F("hPhiRelVsMult", "phi-Psi2 vs N @ 2fm/c;#Delta#phi;N_{partons}", 100, -M_PI, M_PI, 50, 0, 50);
    gStyle->SetOptStat(0);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);
        auto jets  = ExtractJets(b_gPt, b_gEta, b_gPhi, b_gMult);
        auto parts = ExtractPartons(b_pdgid, b_px, b_py, b_pz, b_e, b_x, b_y, b_z, b_t);
        MatchPartonsToJets(parts, jets);
        TransformToJetFrame(parts, jets);
        FillPsi2AtTime(parts, tTarget, hPsi2, hPhiRel, hPsi2VsMult, hPhiRelVsMult);
    }
    {
        TCanvas c1("c1","Psi2",800,600); hPsi2->Draw("HIST"); c1.SaveAs("Psi2_2fm.png");
    }
    {
        TCanvas c2("c2","PhiRel",800,600); hPhiRel->Draw("HIST"); c2.SaveAs("PhiRel_2fm.png");
    }
    {
        TCanvas c3("c3","Psi2VsMult",800,600); hPsi2VsMult->Draw("COLZ"); c3.SaveAs("Psi2VsMult_2fm.png");
    }
    {
        TCanvas c4("c4","PhiRelVsMult",800,600); hPhiRelVsMult->Draw("COLZ"); c4.SaveAs("PhiRelVsMult_2fm.png");
    }
    file->Close();
    return 0;
}
