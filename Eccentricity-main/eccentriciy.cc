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
#include <sstream>

using namespace std;

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

int eccentriciy() {

    // Open ROOT file
    TFile* file = TFile::Open("pp_parton_cascade_0.root", "READ");
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
    // Set Branch Address
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
    
    const float maxEtas = 8.0;
    const float Nbin = 8;
    
    TH1F*  hEtaS = new TH1F("hEtaS", "Space-time rapidity #eta_{s}", 8, 0, 8.0);
    hEtaS->GetXaxis()->SetTitle("#eta_{s}*");
    hEtaS->GetYaxis()->SetTitle("dN/d#eta_{s}*");

    TH2F*  hEcc  = new TH2F("hEcc", "Eccentricity vs #eta_{s}", Nbin, 0, maxEtas, 100, 0.0, 1.0);
    hEcc->GetXaxis()->SetTitle("#eta_{s}*");
    hEcc->GetYaxis()->SetTitle("Eccentricity");
    
    TH2F* hTauEta   = new TH2F("hTauEta", "Proper time #tau vs #eta_{s}", Nbin, 0, maxEtas, 100, 0.0, 10.0);
    hTauEta->GetXaxis()->SetTitle("#eta_{s}*");
    hTauEta->GetYaxis()->SetTitle("#tau* [fm/c]");

    TH2F *hMult = new TH2F("hMult", "Parton multiplicity vs #eta_{s}", Nbin, 0, maxEtas, 100, 0.0, 100.0);
    hMult->GetXaxis()->SetTitle("#eta_{s}*");
    hMult->GetYaxis()->SetTitle("Parton multiplicity");

    TH2F* hPsi2VsMult = new TH2F("hPsi2VsMult", "#Psi_{2} vs Parton Multiplicity;#Psi_{2};N_{partons}", 100, -TMath::Pi()/2, TMath::Pi()/2,50, 0, 50);
    TH2F* hphiminusPsi2VsMult = new TH2F("hPhiminusPsi2VsMult", "#Phi-Psi_{2} vs Parton Multiplicity;#phi-Psi_{2};N_{partons}", 100, -TMath::Pi()/2, TMath::Pi()/2,50, 0, 50);

    TH2F *hJetMultParton = new TH2F("hJetMultParton", ";N_{ch}^{j};N_{p}^{j}", 50, 0, 100, 25, 0, 50);

    // Event‐plane angle distribution
    TH1F* hPsi2 = new TH1F("hPsi2", "Event plane #Psi_{2};#Psi_{2};counts", 50, -TMath::Pi()/2, TMath::Pi()/2);

    // Distribution of φ – Ψ₂
    TH1F* hPhiRel = new TH1F("hPhiRel", "#phi - #Psi_{2};#Delta#phi;counts", 100, -TMath::Pi(), TMath::Pi());


    const vector<int> selectedEvents = {284}; // Events for phi dist plotting
    unordered_map<int, vector<TH1F*>> phiDists;   // event_id -> histograms

    // Style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    //----------------------------------------------------------------------
    // Containers
    //----------------------------------------------------------------------

    vector<Parton>           partons;
    vector<Jet>              jets;

    const float              xlow = hEcc->GetXaxis()->GetXmin();
    const float              xup  = hEcc->GetXaxis()->GetXmax();
    const float              binw = (xup - xlow) / Nbin;

    //----------------------------------------------------------------------
    // Event Loop
    //----------------------------------------------------------------------

    Long64_t nentries = tree->GetEntries();

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        bool isSelectedEvent = find(selectedEvents.begin(), selectedEvents.end(), ie) != selectedEvents.end();

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

            // Calculate Momentum-related quantities
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
        
        /*std::cout << "Event " << ie
                  << ":  partons=" << partons.size()
                  << ", jets="    << jets.size() 
                  << ", first jet mult="    << jets[0].genJetChargedMultiplicity
                  <<", second jet mult="    << jets[1].genJetChargedMultiplicity<<"\n";
        */         

        //----------------------------------------------------------------------
        // Coordinate transformation: Lab -> Jet Frame
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

            if (isSelectedEvent) {
                std::stringstream name;
                name << "phiDist_evt" << ie << "_jet" << p.jetID;
                if (phiDists[ie].size() <= p.jetID) phiDists[ie].resize(p.jetID + 1, nullptr);
                    if (!phiDists[ie][p.jetID]) {
                    phiDists[ie][p.jetID] = new TH1F(name.str().c_str(), ("Jet-frame #phi for event " + std::to_string(ie) + ", jet " + std::to_string(p.jetID) + ";#phi (jet frame);dN/d#phi").c_str(), 24, -TMath::Pi(), TMath::Pi());
                }
                phiDists[ie][p.jetID]->Fill(phJ);
            }


            p.jet_par_x = ptJ * cos(phJ);
            p.jet_par_y = ptJ * sin(phJ);
            p.jet_par_z = ptJ / tan(thJ);

            p.tau = sqrt(p.t * p.t - p.jet_par_z * p.jet_par_z);

            // Calculate space-time rapidity
            if (p.t > p.jet_par_z) {
                p.eta_s = 0.5f * log((p.t + p.jet_par_z) / (p.t - p.jet_par_z));
                hEtaS->Fill(p.eta_s);
                hTauEta->Fill(p.eta_s, p.tau);
            }
            else {
                p.eta_s = -999.9f;
            //    cout<<"Warning: t <= z !!!"<<endl;
            //    cout<<"t="<<p.t<<"\tz="<<p.jet_par_z<<endl;
            }
        }

        int countUsedPartons = 0;
        for (auto& p : partons) {
            if (p.jetID < 0) continue;
                countUsedPartons++;
        }

        // Throw out low parton multiplicity events
        if (countUsedPartons < 3) {
            std::cerr << "Skipping Ψ₂ for event " << ie << " due to only " << countUsedPartons << " partons\n";
            continue;
        }


        // compute Ψ₂ for this event
        double sumCos2 = 0.0, sumSin2 = 0.0;
        for (auto& p : partons) {
            if (p.jetID < 0) continue;
                TVector3 jv; 
                jv.SetPtEtaPhi(jets[p.jetID].Pt,
                    jets[p.jetID].Eta,
                    jets[p.jetID].Phi);
                TVector3 pos(p.x, p.y, p.z);
                float phJ = phiWRTJet(jv, pos);
                float r2  = p.jet_par_x*p.jet_par_x + p.jet_par_y*p.jet_par_y;
                sumCos2 += r2 * cos(2*phJ);
                sumSin2 += r2 * sin(2*phJ);
        }

        double psi2 = 0.5 * atan2(sumSin2, sumCos2);
        if (!std::isnan(psi2)) {
            hPsi2->Fill(psi2);
            hPsi2VsMult->Fill(psi2, countUsedPartons);
        }


        // fill φ−Ψ₂ histogram
        for (auto& p : partons) {
            if (p.jetID < 0) continue;
                TVector3 jv; 
                jv.SetPtEtaPhi(jets[p.jetID].Pt,
                    jets[p.jetID].Eta,
                    jets[p.jetID].Phi);
                TVector3 pos(p.x, p.y, p.z);
                float phJ = phiWRTJet(jv, pos);
                double dphi = TVector2::Phi_mpi_pi(phJ - psi2);
                hPhiRel->Fill(dphi);
                hphiminusPsi2VsMult->Fill(dphi, countUsedPartons);
        }


        //----------------------------------------------------------------------
        // Group by eta_s and jetID and calculate eccentricity
        //----------------------------------------------------------------------

        vector< unordered_map<int, vector<int>> > bins(Nbin);

        for (int ip = 0; ip < (int)partons.size(); ++ip) {
            auto& p = partons[ip];

            if (p.jetID < 0) continue;

            int ib = int((p.eta_s - xlow) / binw);
            if (ib < 0 || ib >= Nbin) continue;

            bins[ib][p.jetID].push_back(ip);
        }

        for (int ij = 0; ij < jets.size(); ij++)
        {
            int nPartons = 0;
            for (int ib = 0; ib < Nbin; ib++)
            {
                nPartons += bins[ib][ij].size();
            }

            hJetMultParton->Fill(jets[ij].genJetChargedMultiplicity, nPartons);
        }

        for (int ib = 0; ib < Nbin; ++ib) {
            float xcenter = xlow + (ib + 0.5f) * binw;

            for (auto& kv : bins[ib]) {
                auto& idxs = kv.second;
                if (idxs.empty()) continue;

                double sumE = 0.0;
                double sumx = 0.0;
                double sumy = 0.0;

                for (int idx : idxs) {
                    auto& P = partons[idx];
                    double w = P.e;
                    sumE += w;
                    sumx += w * P.jet_par_x;
                    sumy += w * P.jet_par_y;
                }

                double xm = sumx / sumE;
                double ym = sumy / sumE;

                double real = 0.0;
                double imag = 0.0;
                double norm = 0.0;
                double eventplane = 0.0;
               
                for (int idx : idxs) {
                    auto& P = partons[idx];

                    double dx = P.jet_par_x - xm;
                    double dy = P.jet_par_y - ym;
                    double r2 = dx*dx + dy*dy;
                    double phi_loc = atan2(dy, dx);
                    double w = P.e * r2;

                    real += w * cos(2*phi_loc);
                    imag += w * sin(2*phi_loc);
                    //double eventplane = 1/2 * atan(imag/real);
                    //cout<<"eventplane: " <<eventplane<<endl;
                    norm += w;
                }

                double ecc = (norm > 0.0)
                           ? sqrt(real*real + imag*imag) / norm
                           : 0.0;

                hEcc->Fill(xcenter, ecc);
                
                // fill parton multiplicity
                int mult = idxs.size();       //number of partons in eta_s bin for the given jet
                hMult->Fill(xcenter, mult);


                // compute Ψ₂ for each event for (n=2)
                double psi2 = 0.5 * atan(imag/real);

                // print event information
                cout << "Event " << ie
                << ", Psi2=" << psi2
                << ", imag=" << imag
                << ", real=" << real
                << ", partons=" << partons.size()
                << ", jets=" << jets.size()
                << ", first jet parton number=" << jets[0].genJetChargedMultiplicity
                << ", second jet parton number=" << jets[1].genJetChargedMultiplicity
                << "\n";

            }
        }

    
    }

    //----------------------------------------------------------------------
    // Draw and Save Plots
    //----------------------------------------------------------------------

    /*
    for (auto& evt : phiDists) {
    int ie = evt.first;
    auto& jetHists = evt.second;

        for (size_t j = 0; j < jetHists.size(); ++j) {
        TCanvas* cPhi = new TCanvas();
            if (jetHists[j]) {
            jetHists[j]->SetLineWidth(2);
            jetHists[j]->Draw("E");
            //jetHists[j]->Draw("HIST");
            //jetHists[j]->Draw("E SAME");

            }
        std::stringstream fname;
        fname << "Abb/phiDist_evt" << ie << "_jet" << j << ".png";
        cPhi->SaveAs(fname.str().c_str());
        //delete cPhi;
        }
    }


    // eta_s distribution
    TCanvas* c1 = new TCanvas("c1", "eta_s distribution", 800, 600);
    hEtaS->Scale(1. / hEtaS->GetXaxis()->GetBinWidth(1));
    hEtaS->SetMarkerStyle(20);
    hEtaS->SetLineWidth(2);
    hEtaS->Draw("E");
    c1->SaveAs("Abb/eta_s_distribution.png");

    // Eccentricity 2D heatmap (log Z)
    TCanvas* c2 = new TCanvas("c2", "eccentricity vs eta_s", 800, 600);
    c2->SetLogz();
    hEcc->Draw("COLZ");
    c2->SaveAs("Abb/eccentricity_vs_eta_s_COLZ.png");

    // Average eccentricity curve
    TCanvas* c3 = new TCanvas("c3", "Mean eccentricity vs eta_s", 800, 600);
    TProfile* pEcc = hEcc->ProfileX("pEcc");
    pEcc->SetMarkerStyle(20);
    pEcc->SetLineWidth(2);
    pEcc->GetXaxis()->SetTitle("#eta_{s}");
    pEcc->GetYaxis()->SetTitle("<eccentricity>");
    pEcc->Draw("E");
    c3->SaveAs("Abb/eccentricity_vs_eta_s_Profile.png");
    
    // plot tau vs eta_s
    TCanvas* c4 = new TCanvas("c4", "tau vs eta_s", 800, 600);
    hTauEta->Draw("COLZ");
    c4->SaveAs("Abb/tau_vs_eta_s_COLZ.png");

    // Average tau curve
    TCanvas* c5 = new TCanvas("c5", "Mean tau vs eta_s", 800, 600);
    TProfile* pTau = hTauEta->ProfileX("pTau");
    pTau->SetMarkerStyle(20);
    pTau->SetLineWidth(2);
    pTau->GetXaxis()->SetTitle("#eta_{s}");
    pTau->GetYaxis()->SetTitle("<#tau>");
    pTau->Draw("E");
    c5->SaveAs("Abb/tau_vs_eta_s_Profile.png");
    
    // --- draw parton multiplicity 2D heatmap (log Z) ---
    TCanvas* cMult = new TCanvas("cMult", "Parton multiplicity vs eta_s", 800, 600);
    cMult->SetLogz();
    hMult->Draw("COLZ");
    cMult->SaveAs("Abb/parton_multiplicity_vs_eta_s_COLZ.png");

    // Draw average parton multiplicity curve
    TCanvas* cMultProf = new TCanvas("cMultProf", "Mean parton multiplicity vs eta_s", 800, 600);
    TProfile* pMult = hMult->ProfileX("pMult");
    pMult->SetMarkerStyle(20);
    pMult->SetLineWidth(2);
    pMult->GetXaxis()->SetTitle("#eta_{s}");
    pMult->GetYaxis()->SetTitle("<# partons>");
    pMult->Draw("E");
    cMultProf->SaveAs("Abb/parton_multiplicity_vs_eta_s_Profile.png");
*/
    //Draw jet charged multiplicity vs parton multiplicity heatmap
    TCanvas* c6 = new TCanvas("c6", "", 800, 600);
    hJetMultParton->Draw("COLZ");
    c6->SaveAs("Abb/JetMult_vs_nPartons_COLZ.png");

    // Draw average parton count vs jet charged multiplicity
    TCanvas* c7 = new TCanvas("c7", "", 800, 600);
    TProfile* pJetMultParton = hJetMultParton->ProfileX("pJetMultParton");
    pJetMultParton->SetMarkerStyle(20);
    pJetMultParton->SetLineWidth(2);
    pJetMultParton->GetXaxis()->SetTitle("N_{ch}^{j}");
    pJetMultParton->GetYaxis()->SetTitle("<N_{p}^{j}>");
    pJetMultParton->Draw("E");
    c7->SaveAs("Abb/JetMult_vs_nPartons_Profile.png");


    // --- draw Psi2 distribution ---
    TCanvas* cPsi2 = new TCanvas("cPsi2","Psi2 distribution",800,600);
    hPsi2->Draw("HIST");
    cPsi2->SaveAs("Abb/psi2_distribution.png");

    // --- draw Phi minus Psi2 distribution ---
    TCanvas* cPhiRel = new TCanvas("cPhiRel","Phi minus Psi2",800,600);
    hPhiRel->Draw("HIST");
    cPhiRel->SaveAs("Abb/phi_minus_psi2_distribution.png");

    // --- draw Psi2 vs Parton Multiplicity 2D heatmap (log Z) ---
    TCanvas* cPsi2VsMult = new TCanvas("cPsi2VsMult", "#Psi_{2} vs Multiplicity", 800, 600);
    cPsi2VsMult->SetLogz(); // optional for contrast
    hPsi2VsMult->Draw("COLZ");
    cPsi2VsMult->SaveAs("Abb/psi2_vs_partonMultiplicity_COLZ.png");

    TCanvas* cphiminusPsi2VsMult = new TCanvas("phiminuscPsi2VsMult", "#phiminusPsi_{2} vs Multiplicity", 800, 600);
    cphiminusPsi2VsMult->SetLogz(); // optional for contrast
    hphiminusPsi2VsMult->Draw("COLZ");
    cphiminusPsi2VsMult->SaveAs("Abb/phiminuspsi2_vs_partonMultiplicity_COLZ.png");

    //file->Close();
    return 0;
}

