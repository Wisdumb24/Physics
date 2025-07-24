#include <vector>
#include <unordered_map>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector2.h"
#include "coordinateTools.h"

using namespace std;

//C:\Users\david\OneDrive\Desktop\Rice-Physics\Eccentricity-main

//----------------------------------------------------------------------
// Defining Data Structures
//----------------------------------------------------------------------

struct Parton {
    int   pdgid;

    float px, py, pz, e;
    float x,  y,  z,  t;
    float pt, eta, phi;
    float eta_s, tau; // defined in jet frame

    int   jetID;

    float jet_par_x;
    float jet_par_y;
    float jet_par_z;
    float jet_par_px;
    float jet_par_py;
    float jet_par_pz;
};

struct Jet {
    float Pt;
    float Eta;
    float Phi;
    int   genJetChargedMultiplicity;
};

//----------------------------------------------------------------------
// Main Function
//----------------------------------------------------------------------

int evolution() {

    // Open ROOT File
    TFile* file = TFile::Open("pp_parton_cascade_1.root", "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: cannot open file pp_parton_cascade.root" << endl;
        return 1;
    }

    // Get TTree
    TTree* tree = dynamic_cast<TTree*>(file->Get("trackTree"));
    if (!tree) {
        cerr << "Error: cannot find TTree 'trackTree' in file" << endl;
        file->Close();
        return 1;
    }

    //----------------------------------------------------------------------
    // Set Branch Addresses
    //----------------------------------------------------------------------
    
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
    vector<float> *b_genJetPhi = nullptr;
    vector<int> *b_genJetChargedMult = nullptr;

    tree->SetBranchAddress("genJetPt",  &b_genJetPt);
    tree->SetBranchAddress("genJetEta", &b_genJetEta);
    tree->SetBranchAddress("genJetPhi", &b_genJetPhi);
    tree->SetBranchAddress("genJetChargedMultiplicity", &b_genJetChargedMult);

    //----------------------------------------------------------------------
    // Initialize Histograms
    //----------------------------------------------------------------------

    double maxTime = 4;

    TH1F *hTime = new TH1F("hTime", "Generation time distribution; t (fm/c); Entries", 100, 0, 100);

    TH2F *hMultTime = new TH2F("hMultTime", "Jet parton multiplicity vs time; t (fm/c); N_{parton}", 40, 0, maxTime, 50, 0, 50);
    TH2F *hEccTime = new TH2F("hEccTime", "Eccentricity vs time; t (fm/c); #epsilon_{2}", 40, 0, maxTime, 50, 0, 10);
    TH2F *hDimTime = new TH2F("hDimTime", "Dimension vs time; t (fm/c); r (fm)", 40, 0, maxTime, 50, 0, 2);
    TH2F *hDenTime = new TH2F("hDenTime", "Density vs time; t (fm/c); density (fm^{-2})", 40, 0, maxTime, 200, 0, 200);
    
    // Ψ2 / Δφ histograms  (early = bin-1 center, late = last-bin center)
    TH1F *hPsi2Early = new TH1F("hPsi2Early", "Early #Psi_{2};#Psi_{2};counts", 25, -TMath::Pi()/2, TMath::Pi()/2);
    TH1F *hPsi2Late  = new TH1F("hPsi2Late", "Late #Psi_{2};#Psi_{2};counts", 50, -TMath::Pi()/2, TMath::Pi()/2);

    TH1F *hPhiRelEarly = new TH1F("hPhiRelEarly", "#phi-#Psi_{2} (early);#Delta#phi;counts", 50, -TMath::Pi(), TMath::Pi());
    TH1F *hPhiRelLate  = new TH1F("hPhiRelLate", "#phi-#Psi_{2} (late);#Delta#phi;counts", 100, -TMath::Pi(), TMath::Pi());

    TH2F *hPsi2VsMultEarly = new TH2F("hPsi2VsMultEarly", "#Psi_{2} vs N (early);#Psi_{2};N_{partons}", 50, -TMath::Pi()/2, TMath::Pi()/2, 50, 0, 50);
    TH2F *hPsi2VsMultLate  = new TH2F("hPsi2VsMultLate", "#Psi_{2} vs N (late);#Psi_{2};N_{partons}", 50, -TMath::Pi()/2, TMath::Pi()/2, 50, 0, 50);

    TH2F *hPhiRelVsMultEarly = new TH2F("hPhiRelVsMultEarly", "#phi-#Psi_{2} vs N (early);#Delta#phi;N_{partons}", 100, -TMath::Pi(), TMath::Pi(), 50, 0, 50);
    TH2F *hPhiRelVsMultLate  = new TH2F("hPhiRelVsMultLate", "#phi-#Psi_{2} vs N (late);#Delta#phi;N_{partons}",100, -TMath::Pi(), TMath::Pi(), 50, 0, 50);

    //make sure histograms are empty
    // Reset histograms in case they persist from a previous run
    hTime->Reset();
    hMultTime->Reset();
    hEccTime->Reset();
    hDimTime->Reset();
    hDenTime->Reset();

    hPsi2Early->Reset();
    hPsi2Late->Reset();
    hPhiRelEarly->Reset();
    hPhiRelLate->Reset();
    hPsi2VsMultEarly->Reset();
    hPsi2VsMultLate->Reset();
    hPhiRelVsMultEarly->Reset();
    hPhiRelVsMultLate->Reset();

    // Style Settings
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    //----------------------------------------------------------------------
    // Prepare Containers
    //----------------------------------------------------------------------

    vector<Parton>           partons;
    vector<Jet>              jets;

    //----------------------------------------------------------------------
    // Event Loop
    //----------------------------------------------------------------------

    Long64_t nentries = tree->GetEntries();

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        // Read and store jets
        jets.clear();
        for (int j = 0; j < (int)b_genJetPt->size(); ++j) {
            Jet jet;
            jet.Pt  = b_genJetPt->at(j);
            jet.Eta = b_genJetEta->at(j);
            jet.Phi = b_genJetPhi->at(j);
            jet.genJetChargedMultiplicity = b_genJetChargedMult->at(j);
            jets.push_back(jet);
        }

        // Read and store partons
        partons.clear();

        for (int j = 0; j < (int)b_par_pdgid->size(); ++j) {
            Parton p;

            p.pdgid = b_par_pdgid->at(j);
            p.px    = b_par_px->at(j);
            p.py    = b_par_py->at(j);
            p.pz    = b_par_pz->at(j);
            p.e     = b_par_e->at(j);

            p.x     = b_par_x->at(j);
            p.y     = b_par_y->at(j);
            p.z     = b_par_z->at(j);
            p.t     = b_par_t->at(j);

            // Calculate parton kinematics
            p.pt   = sqrt(p.px*p.px + p.py*p.py);
            p.phi  = atan2(p.py, p.px);
            float theta = atan2(p.pt, p.pz);
            p.eta  = -log(tan(theta/2.0f));

            p.jetID       = -1;
            p.jet_par_x   = 0.0f;
            p.jet_par_y   = 0.0f;
            p.jet_par_z   = 0.0f;

            partons.push_back(p);
        }

        //----------------------------------------------------------------------
        // Match Partons to Jets
        //----------------------------------------------------------------------

        for (auto& p : partons) {
            for (int j = 0; j < (int)jets.size(); ++j) {
                float dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
                float dEta = p.eta - jets[j].Eta;
                float dR   = sqrt(dEta*dEta + dphi*dphi);

                if (dR < 0.8f) {
                    p.jetID = j;
                    break;
                }
            }
        }
        
        //Print event information
        /*
        std::cout << "Event " << ie
        << ":  partons=" << partons.size()
        << ", jets="    << jets.size() <<"\n";
        << ", first jet mult="    << jets[0].genJetChargedMultiplicity
        <<", second jet mult="    << jets[1].genJetChargedMultiplicity<<"\n";
        */

        //----------------------------------------------------------------------
        // Coordinate Transformation: Lab -> Jet Frame
        //----------------------------------------------------------------------

        for (auto& p : partons) {
            if (p.jetID < 0) continue;

            if ((p.x == 0) && (p.y == 0) && (p.z == 0))
            {
                p.jet_par_x = 0;
                p.jet_par_y = 0;
                p.jet_par_z = 0;
                continue;
            }

            TVector3 jv;
            jv.SetPtEtaPhi(
                jets[p.jetID].Pt,
                jets[p.jetID].Eta,
                jets[p.jetID].Phi
            );

            TVector3 pos(p.x, p.y, p.z);

            float ptJ = ptWRTJet  (jv, pos);
            float thJ = thetaWRTJet(jv, pos);
            float phJ = phiWRTJet  (jv, pos);

            p.jet_par_x = ptJ * cos(phJ);
            p.jet_par_y = ptJ * sin(phJ);
            p.jet_par_z = ptJ / tan(thJ);

            // rotate momentum into the jet frame
            TVector3 mom(p.px, p.py, p.pz);

            float ptJ_m = ptWRTJet  (jv, mom);
            float thJ_m = thetaWRTJet(jv, mom);
            float phJ_m = phiWRTJet  (jv, mom);

            p.jet_par_px = ptJ_m * cos(phJ_m);
            p.jet_par_py = ptJ_m * sin(phJ_m);
            p.jet_par_pz = ptJ_m / tan(thJ_m);


            p.tau = sqrt(p.t * p.t - p.jet_par_z * p.jet_par_z);

            // Calculate eta_s
            if (p.t > p.jet_par_z) {
                p.eta_s = 0.5f * log((p.t + p.jet_par_z) / (p.t - p.jet_par_z));
            }
            else {
                p.eta_s = -999.9f;
            //    cout<<"Warning: t <= z !!!"<<endl;
            //    cout<<"t="<<p.t<<"\tz="<<p.jet_par_z<<endl;
            }
        }

        //----------------------------------------------------------------------
        // Group by JetID and compute eccentricity
        //----------------------------------------------------------------------

 unordered_map<int, vector<int>> partonsByJet;

// -------------------------------------------------------------
// Helper: compute Ψ2 and fill requested histograms at tTarget
// -------------------------------------------------------------
auto fillPsi2Set = [&](double tTarget,
                       TH1F* hPsi2, TH1F* hPhiRel,
                       TH2F* hPsi2VsMult, TH2F* hPhiRelVsMult)
{
    for (auto &kv : partonsByJet)
    {
        const std::vector<int> &idxs = kv.second;
        if (idxs.size() <= 2) continue;        // skip jets with ≤1 partons

        // ------------ propagate & centroid ------------------
        double sumE=0, sumx=0, sumy=0;
        for (int idx : idxs) {
            const auto &P = partons[idx];
            if (P.t > tTarget) continue;
            double dt = tTarget - P.t;
            double x  = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y  = P.jet_par_y + dt * P.jet_par_py / P.e;
            sumE += P.e;
            sumx += P.e * x;
            sumy += P.e * y;
        }
        if (sumE==0) continue;                // no contributors
        double xm = sumx / sumE,  ym = sumy / sumE;

        // ------------ Q-vector ------------------------------
        double Re=0, Im=0;
        for (int idx : idxs) {
            const auto &P = partons[idx];
            if (P.t > tTarget) continue;
            double dt = tTarget - P.t;
            double x  = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y  = P.jet_par_y + dt * P.jet_par_py / P.e;
            double r2 = (x-xm)*(x-xm) + (y-ym)*(y-ym);
            double phi = atan2(y-ym, x-xm);
            double w   = P.e * r2;            // weight
            Re += w * cos(2*phi);
            Im += w * sin(2*phi);
        }
        if (Re==0 && Im==0) continue;
        double psi2 = 0.5 * atan2(Im, Re);
        int    mult = idxs.size();

        //if (fabs(psi2)>=0.075) continue;
        //cout << "psi2 = " << psi2 << endl;
        
        hPsi2->Fill(psi2);
        hPsi2VsMult->Fill(psi2, mult);

        for (int idx : idxs) {
            const auto &P = partons[idx];
            if (P.t > tTarget) continue;
            //double dt = tTarget - P.t;
            //double x  = P.jet_par_x + dt * P.jet_par_px / P.e;
            //double y  = P.jet_par_y + dt * P.jet_par_py / P.e;
            //double phi = atan2(y-ym, x-xm);
            //double dphi = TVector2::Phi_mpi_pi(phi - psi2);
            double phi_mom = atan2(P.jet_par_py, P.jet_par_px);  // momentum φ in jet frame
            double dphi    = TVector2::Phi_mpi_pi(phi_mom - psi2);  // use phi_mom here
            hPhiRel->Fill(dphi);
            hPhiRelVsMult->Fill(dphi, mult);
        }
    }
};

       

        for (int ip = 0; ip < (int)partons.size(); ++ip) {
            auto& p = partons[ip];

            if (p.jetID < 0) continue;

            partonsByJet[p.jetID].push_back(ip);
            hTime->Fill(p.t);
        }

        const int  earlyBin = 21;
        const int  lateBin  = hMultTime->GetNbinsX();
        const double tEarly = hMultTime->GetXaxis()->GetBinCenter(earlyBin);   
        const double tLate  = hMultTime->GetXaxis()->GetBinCenter(lateBin);   

        fillPsi2Set(tEarly, hPsi2Early, hPhiRelEarly, hPsi2VsMultEarly, hPhiRelVsMultEarly);
        fillPsi2Set(tLate,  hPsi2Late,  hPhiRelLate,  hPsi2VsMultLate,  hPhiRelVsMultLate);


        for (int it = 1; it <= hMultTime->GetNbinsX(); it++)
        {
            double t = hMultTime->GetXaxis()->GetBinCenter(it);

            for (auto &kv : partonsByJet)
            {
                auto &idxs = kv.second;
                if (idxs.empty())
                    continue;

                int nPartons = 0;
                double sumE = 0.0;
                double sumx = 0.0;
                double sumy = 0.0;

                for (int idx : idxs)
                {
                    auto &P = partons[idx];
                    if (P.t < t) 
                    {
                        nPartons++;

                        double x_jetzt = P.jet_par_x + (t - P.t) * P.jet_par_px / P.e;
                        double y_jetzt = P.jet_par_y + (t - P.t) * P.jet_par_py / P.e;

                        double w = P.e;
                        sumE += w;
                        sumx += w * x_jetzt;
                        sumy += w * y_jetzt;
                    }
                }
                
                hMultTime->Fill(t, nPartons);

                // eccentricity calculation, only proceed if it is well-defined
                if (nPartons < 2) continue;

                double xm = sumx / sumE;
                double ym = sumy / sumE;

                double real = 0.0;
                double imag = 0.0;
                double norm = 0.0;

                for (int idx : idxs) 
                {
                    auto &P = partons[idx];
                    if (P.t < t) 
                    {
                        double x_jetzt = P.jet_par_x + (t - P.t) * P.jet_par_px / P.e;
                        double y_jetzt = P.jet_par_y + (t - P.t) * P.jet_par_py / P.e;

                        double dx = x_jetzt - xm;
                        double dy = y_jetzt - ym;
                        double r2 = dx * dx + dy * dy;
                        double phi_loc = atan2(dy, dx);
                        double w = P.e * r2;

                        real += w * cos(2 * phi_loc);
                        imag += w * sin(2 * phi_loc);
                        norm += w;
                    }
                }

                double ecc = (norm > 0.0) ? sqrt(real * real + imag * imag) / norm : 0.0;
                double dimen = sqrt(norm / sumE);
                double dens = nPartons / pow(dimen, 2);

                hEccTime->Fill(t, ecc);
                hDimTime->Fill(t, dimen);
                hDenTime->Fill(t, dens);
                if (ecc > 1.0000001)
                    cout << "WARNING: "<<ecc << endl;
            }
        }
    }

    //----------------------------------------------------------------------
    // Draw and Save Plots
    //----------------------------------------------------------------------

    /*
    TCanvas *cTime = new TCanvas("cTime", "", 800, 600);
    //cTime->SetLogy();
    hTime->Draw();
    cTime->SaveAs("Zeit/TimeDist.png");

    TCanvas *cMultTime = new TCanvas("cMultTime", "Multiplicity vs Time", 800, 600);
    gStyle->SetOptStat(0);
    hMultTime->Draw("COLZ");
    cMultTime->SaveAs("Zeit/JetPartonMultiplicityTime.png");

    TCanvas* c1 = new TCanvas("c1", "", 800, 600);
    TProfile* pMultTime = hMultTime->ProfileX("pMultTime");
    pMultTime->SetMarkerStyle(20);
    pMultTime->SetLineWidth(2);
    pMultTime->SetLineColor(1);
    pMultTime->GetXaxis()->SetTitle("t (fm/c)");
    pMultTime->GetYaxis()->SetTitle("<N_{parton}>");
    pMultTime->Draw("E");
    c1->SaveAs("Zeit/JetPartonMultiplicityTime_Profile.png");

    */

    TCanvas *cEccTime = new TCanvas("cEccTime", "Eccentricity vs Time", 800, 600);
    gStyle->SetOptStat(0);
    hEccTime->Draw("COLZ");
    cEccTime->SaveAs("Zeit/EccentricityTime.png");

    TCanvas* c2 = new TCanvas("c2", "", 800, 600);
    TProfile* pEccTime = hEccTime->ProfileX("pEccTime");
    pEccTime->SetMarkerStyle(20);
    pEccTime->SetLineWidth(2);
    pEccTime->SetLineColor(1);
    pEccTime->GetXaxis()->SetTitle("t (fm/c)");
    pEccTime->GetYaxis()->SetTitle("<#epsilon_{2}>");
    pEccTime->Draw("E");
    c2->SaveAs("Zeit/EccentricityTime_Profile.png");

    /*

    TCanvas *cDimTime = new TCanvas("cDimTime", "Size vs Time", 800, 600);
    gStyle->SetOptStat(0);
    cDimTime->SetLogz();
    hDimTime->Draw("COLZ");
    cDimTime->SaveAs("Zeit/DimensionTime.png");

    TCanvas* c3 = new TCanvas("c3", "", 800, 600);
    TProfile* pDimTime = hDimTime->ProfileX("pDimTime");
    pDimTime->SetMarkerStyle(20);
    pDimTime->SetLineWidth(2);
    pDimTime->SetLineColor(1);
    pDimTime->GetXaxis()->SetTitle("t (fm/c)");
    pDimTime->GetYaxis()->SetTitle("<r> (fm)");
    pDimTime->Draw("E");
    c3->SaveAs("Zeit/DimensionTime_Profile.png");

    TCanvas *cDenTime = new TCanvas("cDenTime", "Density vs Time", 800, 600);
    gStyle->SetOptStat(0);
    cDenTime->SetLogz();
    hDenTime->Draw("COLZ");
    cDenTime->SaveAs("Zeit/DensityTime.png");

    TCanvas* c4 = new TCanvas("c4", "", 800, 600);
    TProfile* pDenTime = hDenTime->ProfileX("pDenTime");
    pDenTime->SetMarkerStyle(20);
    pDenTime->SetLineWidth(2);
    pDenTime->SetLineColor(1);
    pDenTime->GetXaxis()->SetTitle("t (fm/c)");
    pDenTime->GetYaxis()->SetTitle("<density> (fm^{-2})");
    pDenTime->Draw("E");
    c4->SaveAs("Zeit/DensityTime_Profile.png");

    */

    // -------------------------------------------------------------
    //  Ψ2 + Δφ  plots  (EARLY)
    // -------------------------------------------------------------
    TCanvas *cE1 = new TCanvas("cE1", "Psi2 distribution (early)", 800, 600);
    //hPsi2Early->SetLineWidth(2);
    //hPsi2Early->GetYaxis()->SetRangeUser(200, 1200);
    int mid   = hPsi2Early->FindBin(0.00001);          // bin that contains Ψ₂ ≈ 0
    int count = hPsi2Early->GetBinContent(mid);
    std::cout << mid << " " << "Central bin content = "<< count  << '\n';
    hPsi2Early->Draw("HIST");
    cE1->SaveAs("Zeit/E_Psi2_distribution_13.png");

    TCanvas *cE2 = new TCanvas("cE2", "Phi - Psi2 (early)", 800, 600);
    hPhiRelEarly->SetLineWidth(2);
    hPhiRelEarly->Draw("HIST");
    cE2->SaveAs("Zeit/E_PhiMinusPsi2_distribution.png");

    TCanvas *cE3 = new TCanvas("cE3", "Psi2 vs Multiplicity (early)", 800, 600);
    hPsi2VsMultEarly->Draw("COLZ");
    cE3->SaveAs("Zeit/E_Psi2_vs_Multiplicity.png");

    TCanvas *cE4 = new TCanvas("cE4", "Phi - Psi2 vs Multiplicity (early)", 800, 600);
    hPhiRelVsMultEarly->Draw("COLZ");
    cE4->SaveAs("Zeit/E_PhiMinusPsi2_vs_Multiplicity.png");


    // -------------------------------------------------------------
    //  Ψ2 + Δφ  plots  (LATE)
    // -------------------------------------------------------------
    
    /*
    TCanvas *cL3 = new TCanvas("cL3", "Ψ₂ (late)", 800, 600);
    hPsi2Late->SetLineWidth(2);
    hPsi2Late->Draw("HIST");
    cL3->SaveAs("Zeit/L_Psi2_distribution.png");

    TCanvas *cL4 = new TCanvas("cL4", "Δϕ (late)", 800, 600);
    hPhiRelLate->SetLineWidth(2);
    hPhiRelLate->Draw("HIST");
    cL4->SaveAs("Zeit/L_PhiMinusPsi2_distribution.png");


    TCanvas *cL1=new TCanvas("cL1","",800,600);
    hPsi2VsMultLate->Draw("COLZ");
    cL1->SaveAs("Zeit/L_Psi2_vs_Multiplicity.png");

    TCanvas *cL2=new TCanvas("cL2","",800,600);
    hPhiRelVsMultLate->Draw("COLZ");
    cL2->SaveAs("Zeit/L_PhiMinusPsi2_vs_Multiplicity.png");
    */

    file->Close();
    return 0;
}
