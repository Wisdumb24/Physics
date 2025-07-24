//C:\Users\david\OneDrive\Desktop\Rice-Physics\Eccentricity-main

#include "trackTree.C"
#include "coordinateTools.h"


const double maxTime = 10;

//initialize parton histograms
TH2D *hMultTime = new TH2D("hMultTime", "Jet parton multiplicity vs time; t (fm/c); N_{parton}", 40, 0, maxTime, 50, 0, 50);
TH1D *hPsi2Early = new TH1D("hPsi2Early123", "Early #Psi_{2};#Psi_{2};counts", 25, -TMath::Pi()/2, TMath::Pi()/2);
TH2D *hPsi2VsMultEarly = new TH2D("hPsi2VsMultEarly", "#Psi_{2} vs N (early);#Psi_{2};N_{partons}", 50, -TMath::Pi()/2, TMath::Pi()/2, 50, 0, 50);
TH2D *hPhiRelVsMultEarly = new TH2D("hPhiRelVsMultEarly", "#phi-#Psi_{2} vs N (early);#Delta#phi;N_{partons}", 100, -TMath::Pi(), TMath::Pi(), 50, 0, 50);
TH1D* hPsi2;
TH1D *hPhiRelEarly = new TH1D("hPhiRelEarly", "#phi-#Psi_{2} (early);#Delta#phi;counts", 50, -TMath::Pi(), TMath::Pi());

//initialize v2 and eccentricity histograms
TH2D* hEccVsTau = new TH2D("hEccVsTau", "#varepsilon_{2} vs proper time; #tau (fm/c); #varepsilon_{2}", 40, 0, maxTime, 50, 0, 1);
TH1D* hEccRMSvsTau = new TH1D("hEccRMSvsTau", "RMS #varepsilon_{2} vs #tau;#tau (fm/c); #varepsilon_{2}^{RMS}", 40, 0, maxTime);  // same binning as hEccVsTau
TH1D* hV2Early = new TH1D("hV2Early", "v_{2} at #tau = 2 fm/c; v_{2}; counts", 50, -1, 1);
TH2D* hV2VsTau = new TH2D("hV2VsTau", "v_{2}^{2} vs proper time; #tau (fm/c); v_{2}^{2}", 40, 0, maxTime, 50, 0, 1);
TH1D* hV2RMSvsTau = new TH1D("hV2RMSvsTau", "RMS v_{2} vs #tau;#tau (fm/c); v_{2}^{RMS}", 40, 0, maxTime);
TH2D* hV2vsEcc = new TH2D("hV2vsEcc", "v_{2} vs #varepsilon_{2};#varepsilon_{2};v_{2}", 50, 0, 1.2, 50, -1, 1.2);

struct Parton {
    int   pdgid;

    double px, py, pz, e;
    double x,  y,  z,  t;
    double pt, eta, phi;
    double eta_s, tau; // defined in jet frame

    int   jetID;

    double jet_par_x;
    double jet_par_y;
    double jet_par_z;
    double jet_par_px;
    double jet_par_py;
    double jet_par_pz;
};

struct Jet {
    double Pt;
    double Eta;
    double Phi;
    int   genJetChargedMultiplicity;
};


void fillPsi2Set(
    double tTarget,
    TH1D* hPsi2, TH1D* hPhiRel,
    TH2D* hPsi2VsMult, TH2D* hPhiRelVsMult,
    std::map<int, std::vector<int>>& partonsByJet,
    std::vector<Parton>& partons
)
{
    for (const auto& [jetID, idxs] : partonsByJet) {
        if (idxs.size() < 2) continue; // Skip jets with <2 partons

        // --- Compute centroid in jet frame at tTarget ---
        double sumE = 0, sumx = 0, sumy = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tTarget) continue;
            double dt = tTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            sumE += P.e;
            sumx += P.e * x;
            sumy += P.e * y;
        }
        if (sumE == 0) continue; // No contributors

        double xm = sumx / sumE;
        double ym = sumy / sumE;

        // --- Q-vector calculation ---
        double Re = 0, Im = 0;  // total weight for eccentricity denominator

        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tTarget) continue;
            double dt = tTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            double dx = x - xm;
            double dy = y - ym;
            double r2 = dx * dx + dy * dy;
            double phi = atan2(dy, dx);
            double w = P.e * r2;
            Re += w * cos(2 * phi);
            Im += w * sin(2 * phi);
        }

        if (Re == 0) continue;

        double psi2 = 0.5 * atan2(Im, Re);
        int mult = static_cast<int>(idxs.size());
		
        hPsi2->Fill(psi2);
        hPsi2VsMult->Fill(psi2, mult);

        // --- Fill Δφ histograms for each parton ---
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tTarget) continue;
            double phi_mom = atan2(P.jet_par_py, P.jet_par_px); // momentum φ in jet frame
            double dphi = TVector2::Phi_mpi_pi(phi_mom-psi2);
            hPhiRel->Fill(dphi);
            hPhiRelVsMult->Fill(dphi, mult);
        }
    }
}

// Read and store jets
void readJets(
	trackTree* tree,
	vector<Jet>& jets
)
{
	jets.clear();
	for (int j = 0; j < (int)tree->genJetPt->size(); ++j) {
		Jet jet;
		jet.Pt = tree->genJetPt->at(j);
		jet.Eta = tree->genJetEta->at(j);
		jet.Phi = tree->genJetPhi->at(j);
		jet.genJetChargedMultiplicity = tree->genJetChargedMultiplicity->at(j);
		jets.push_back(jet);
	}
}

// Read and store partons
void readPartons(
	trackTree* tree,
	vector<Parton>& partons
)
{
	partons.clear();

	for (int j = 0; j < (int)tree->par_pdgid->size(); ++j) {
		Parton p;

		p.pdgid = (double)tree->par_pdgid->at(j);
		p.px    = (double)tree->par_px->at(j);
		p.py    = (double)tree->par_py->at(j);
		p.pz    = (double)tree->par_pz->at(j);
		p.e     = (double)tree->par_e->at(j);

		p.x     = (double)tree->par_x->at(j);
		p.y     = (double)tree->par_y->at(j);
		p.z     = (double)tree->par_z->at(j);
		p.t     = (double)tree->par_t->at(j);

		// Calculate parton kinematics
		p.pt   = sqrt(p.px * p.px + p.py * p.py);
		p.phi  = atan2(p.py, p.px);
		double theta = atan2(p.pt, p.pz);
		p.eta  = -log(tan(theta / 2.0));

		p.jetID       = -1;
		p.jet_par_x   = 0.0;
		p.jet_par_y   = 0.0;
		p.jet_par_z   = 0.0;

		partons.push_back(p);
	}
}

// Match Partons to Jets
void matchPartonsToJets(
	vector<Parton>& partons,
	const vector<Jet>& jets
)
{
	for (auto& p : partons) {
		for (int j = 0; j < (int)jets.size(); ++j) {
			double dphi = TVector2::Phi_mpi_pi(p.phi - jets[j].Phi);
			double dEta = p.eta - jets[j].Eta;
			double dR   = sqrt(dEta * dEta + dphi * dphi);

			if (dR < 0.8) {
				p.jetID = j;
				break;
			}
		}
	}
}

// Coordinate Transformation from Lab to Jet Frame
void transformToJetFrame(
	vector<Parton>& partons,
	const vector<Jet>& jets
)
{
	for (auto& p : partons) {
		if (p.jetID < 0) continue;

		TVector3 jv;
		jv.SetPtEtaPhi(
			jets[p.jetID].Pt,
			jets[p.jetID].Eta,
			jets[p.jetID].Phi
		);

		// rotate momentum into the jet frame
		TVector3 mom(p.px, p.py, p.pz);

		double ptJ_m = ptWRTJet(jv, mom);
		double thJ_m = thetaWRTJet(jv, mom);
		double phJ_m = phiWRTJet(jv, mom);

		p.jet_par_px = ptJ_m * cos(phJ_m);
		p.jet_par_py = ptJ_m * sin(phJ_m);
		p.jet_par_pz = ptJ_m / tan(thJ_m);

		if ((p.x == 0) && (p.y == 0) && (p.z == 0)) {
			p.jet_par_x = 0;
			p.jet_par_y = 0;
			p.jet_par_z = 0;
			continue;
		}

		TVector3 pos(p.x, p.y, p.z);

		double ptJ = ptWRTJet(jv, pos);
		double thJ = thetaWRTJet(jv, pos);
		double phJ = phiWRTJet(jv, pos);

		p.jet_par_x = ptJ * cos(phJ);
		p.jet_par_y = ptJ * sin(phJ);
		p.jet_par_z = ptJ / tan(thJ);

	}
}

// Group by JetID and compute eccentricity
void groupByJet(
	vector<Parton>& partons,
	std::map<int, std::vector<int>>& partonsByJet
)
{
	partonsByJet.clear();
	int npartons = partons.size();
	for (int ip = 0; ip < npartons; ++ip) {
		auto& p = partons[ip];

		if (p.jetID < 0) continue;

		partonsByJet[p.jetID].push_back(ip);
	}
}

// compute eccentricity vs tau distribution
void fillEccVsTau(
    TH2D* hEccVsTau,
    std::map<int, std::vector<int>>& partonsByJet,
    std::vector<Parton>& partons
)
{
    const int nBins = hEccVsTau->GetNbinsX();

    for (int it = 1; it <= nBins; ++it) {
        double tauTarget = hEccVsTau->GetXaxis()->GetBinCenter(it);

        for (const auto& [jetID, idxs] : partonsByJet) {
            if (idxs.size() <= 2) continue;

            // Compute centroid
            double sumE = 0, sumx = 0, sumy = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau > tauTarget) continue;
                double dt = tauTarget - P.tau;
                double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                sumE += P.e;
                sumx += P.e * x;
                sumy += P.e * y;
            }
            if (sumE == 0) continue;

            double xm = sumx / sumE;
            double ym = sumy / sumE;

            // Q-vector
            double Re = 0, Im = 0, Wtot = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau > tauTarget) continue;
                double dt = tauTarget - P.tau;
                double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                double dx = x - xm;
                double dy = y - ym;
                double r2 = dx * dx + dy * dy;
                double phi = atan2(dy, dx);
                double w = P.e * r2;
                Re += w * cos(2 * phi);
                Im += w * sin(2 * phi);
                Wtot += w;
            }
            if (Wtot == 0) continue;

            double eps2 = sqrt(Re * Re + Im * Im) / Wtot;
            hEccVsTau->Fill(tauTarget, eps2 * eps2);
        }
    }
}

// Compute RMS of eccentricity vs tau
void computeEccentricityRMSvsTau(
    TH2D* hEccVsTau,
    TH1D* hEccRMSvsTau
)
{
    int nBinsX = hEccVsTau->GetNbinsX();

    for (int i = 1; i <= nBinsX; ++i) {
        double sum = 0;
        double count = 0;

        for (int j = 1; j <= hEccVsTau->GetNbinsY(); ++j) {
            double e2 = hEccVsTau->GetYaxis()->GetBinCenter(j);
            double weight = hEccVsTau->GetBinContent(i, j);
            sum += weight * e2;
            count += weight;
        }

        if (count > 0) {
            double meanSq = sum / count;
            double rms = sqrt(meanSq);
            hEccRMSvsTau->SetBinContent(i, rms);
        }
    }
}

// compute v2 distribution at tau tauTarget
void fillV2AtTau(
    double tauTarget,
    TH1D* hV2,
    std::map<int, std::vector<int>>& partonsByJet,
    std::vector<Parton>& partons
)
{
    for (const auto& [jetID, idxs] : partonsByJet) {
        if (idxs.size() < 2) continue;

        // Compute event plane Ψ₂ using spatial eccentricity
        double sumE = 0, sumx = 0, sumy = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double dt = tauTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            sumE += P.e;
            sumx += P.e * x;
            sumy += P.e * y;
        }
        if (sumE == 0) continue;

        double xm = sumx / sumE;
        double ym = sumy / sumE;

        double Re = 0, Im = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double dt = tauTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            double dx = x - xm;
            double dy = y - ym;
            double r2 = dx * dx + dy * dy;
            double phi = atan2(dy, dx);
            double w = P.e * r2;
            Re += w * cos(2 * phi);
            Im += w * sin(2 * phi);
        }

        if (Re == 0 && Im == 0) continue;

        double psi2 = 0.5 * atan2(Im, Re);

        //Compute v₂ using momenta relative to Ψ₂
        double sumV2 = 0;
        int n = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double phi_mom = atan2(P.jet_par_py, P.jet_par_px);
            sumV2 += cos(2 * (phi_mom - psi2));
            ++n;
        }

        if (n > 0) {
            double v2 = sumV2 / n;
            hV2->Fill(v2);
        }
    }
}

//V2 vs tau functions
void fillV2VsTau(
    TH2D* hV2VsTau,
    std::map<int, std::vector<int>>& partonsByJet,
    std::vector<Parton>& partons
)
{
    const int nBins = hV2VsTau->GetNbinsX();

    for (int it = 1; it <= nBins; ++it) {
        double tauTarget = hV2VsTau->GetXaxis()->GetBinCenter(it);

        for (const auto& [jetID, idxs] : partonsByJet) {
            if (idxs.size() < 2) continue;

            // Compute Ψ₂
            double sumE = 0, sumx = 0, sumy = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau > tauTarget) continue;
                double dt = tauTarget - P.tau;
                double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                sumE += P.e;
                sumx += P.e * x;
                sumy += P.e * y;
            }
            if (sumE == 0) continue;

            double xm = sumx / sumE;
            double ym = sumy / sumE;

            double Re = 0, Im = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau > tauTarget) continue;
                double dt = tauTarget - P.tau;
                double x = P.jet_par_x + dt * P.jet_par_px / P.e;
                double y = P.jet_par_y + dt * P.jet_par_py / P.e;
                double dx = x - xm;
                double dy = y - ym;
                double phi = atan2(dy, dx);
                double r2 = dx * dx + dy * dy;
                double w = P.e * r2;
                Re += w * cos(2 * phi);
                Im += w * sin(2 * phi);
            }

            if (Re == 0 && Im == 0) continue;

            double psi2 = 0.5 * atan2(Im, Re);

            double sumV2 = 0;
            int n = 0;
            for (int idx : idxs) {
                const auto& P = partons[idx];
                if (P.tau > tauTarget) continue;
                double phi_mom = atan2(P.jet_par_py, P.jet_par_px);
                sumV2 += cos(2 * (phi_mom - psi2));
                ++n;
            }

            if (n > 0) {
                double v2 = sumV2 / n;
                hV2VsTau->Fill(tauTarget, v2 * v2);
            }
        }
    }
}

//compute RMS of v2 vs tau
void computeV2RMSvsTau(
    TH2D* hV2VsTau,
    TH1D* hV2RMSvsTau
)
{
    int nBinsX = hV2VsTau->GetNbinsX();

    for (int i = 1; i <= nBinsX; ++i) {
        double sum = 0;
        double count = 0;

        for (int j = 1; j <= hV2VsTau->GetNbinsY(); ++j) {
            double v2sq = hV2VsTau->GetYaxis()->GetBinCenter(j);
            double weight = hV2VsTau->GetBinContent(i, j);
            sum += weight * v2sq;
            count += weight;
        }

        if (count > 0) {
            double meanSq = sum / count;
            double rms = sqrt(meanSq);
            hV2RMSvsTau->SetBinContent(i, rms);
        }
    }
}

void fillV2VsEccAtTau(
    double tauTarget,
    TH2D* hV2vsEcc,
    std::map<int, std::vector<int>>& partonsByJet,
    std::vector<Parton>& partons
)
{
    for (const auto& [jetID, idxs] : partonsByJet) {
        if (idxs.size() < 2) continue;

        // --- Compute spatial centroid ---
        double sumE = 0, sumx = 0, sumy = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double dt = tauTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            sumE += P.e;
            sumx += P.e * x;
            sumy += P.e * y;
        }
        if (sumE == 0) continue;

        double xm = sumx / sumE;
        double ym = sumy / sumE;

        // --- Eccentricity Q-vector ---
        double Re = 0, Im = 0, Wtot = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double dt = tauTarget - P.tau;
            double x = P.jet_par_x + dt * P.jet_par_px / P.e;
            double y = P.jet_par_y + dt * P.jet_par_py / P.e;
            double dx = x - xm;
            double dy = y - ym;
            double r2 = dx * dx + dy * dy;
            double phi = atan2(dy, dx);
            double w = P.e * r2;
            Re += w * cos(2 * phi);
            Im += w * sin(2 * phi);
            Wtot += w;
        }

        if (Wtot == 0) continue;

        double eps2 = sqrt(Re * Re + Im * Im) / Wtot;
        double psi2 = 0.5 * atan2(Im, Re);

        // --- Compute v₂ ---
        double sumV2 = 0;
        int n = 0;
        for (int idx : idxs) {
            const auto& P = partons[idx];
            if (P.tau > tauTarget) continue;
            double phi_mom = atan2(P.jet_par_py, P.jet_par_px);
            sumV2 += cos(2 * (phi_mom - psi2));
            ++n;
        }

        if (n == 0) continue;

        double v2 = sumV2 / n;
        hV2vsEcc->Fill(eps2, v2);
    }
}

void test()
{
	const int  earlyBin = 21;
	const int  lateBin  = hMultTime->GetNbinsX();
	//const double tEarly = hMultTime->GetXaxis()->GetBinCenter(earlyBin);  
	const double tauEarly = 2.0; // Proper time in fm/c 
	const double tLate  = hMultTime->GetXaxis()->GetBinCenter(lateBin);

	TFile* inFile = TFile::Open("pp_parton_cascade_1.root", "READ");

	TTree *inTree = (TTree*)inFile->Get("trackTree");

	trackTree *tree = new trackTree(inTree);

	int nentries = inTree->GetEntries();
	int max_entries = 100000;
	//cout << "Total Entries: " << inTree->GetEntries() << endl;

	//----------------------------------------------------------------------
    // Prepare Containers
    //----------------------------------------------------------------------

    vector<Parton>           partons;
    vector<Jet>              jets;

    //----------------------------------------------------------------------
    // Event Loop
    //----------------------------------------------------------------------
    std::map<int, vector<int>> partonsByJet;


	for (int ientry = 0; ientry < nentries; ientry++)
	{
		tree->GetEntry(ientry);
		double pT = tree->par_t->size();

		readJets(tree, jets);
		readPartons(tree, partons);

		matchPartonsToJets(partons, jets);
		transformToJetFrame(partons, jets);

		for (auto& p : partons) {
    		double z2 = p.jet_par_z * p.jet_par_z;
    		p.tau = (p.t * p.t > z2) ? sqrt(p.t * p.t - z2) : 0.0;
		}

		groupByJet(partons, partonsByJet);
		fillEccVsTau(hEccVsTau, partonsByJet, partons);

		// fill histograms
		fillPsi2Set(
			tauEarly, hPsi2Early, hPhiRelEarly, hPsi2VsMultEarly, hPhiRelVsMultEarly,
			partonsByJet, partons
		);

		fillV2AtTau(tauEarly, hV2Early, partonsByJet, partons);
        fillV2VsTau(hV2VsTau, partonsByJet, partons);
        fillV2VsEccAtTau(tauEarly, hV2vsEcc, partonsByJet, partons);

	}

	computeEccentricityRMSvsTau(hEccVsTau, hEccRMSvsTau);
    computeV2RMSvsTau(hV2VsTau, hV2RMSvsTau);

	/*
    // Evenet plane distribution
	TCanvas *cE1 = new TCanvas();
	hPsi2Early->SetLineWidth(2);
	hPsi2Early->GetYaxis()->SetRangeUser(200, 1200);
	hPsi2Early->Draw();
	int mid   = hPsi2Early->FindBin(0.00001);          
	int count = hPsi2Early->GetBinContent(mid);
	// get number of entries
	//cout << "Number of entries in Ψ₂ distribution (early): " << hPsi2Early->GetEntries() << std::endl;
    //std::cout << mid << " " << "Central bin content = "<< count  << std::endl;
	cE1->SaveAs("David_Plots/Psi2Early.png");

	// Δφ distribution
	TCanvas *cPhiRel = new TCanvas();
	hPhiRelEarly->SetLineWidth(2);
	hPhiRelEarly->SetLineColor(kBlue);
	hPhiRelEarly->Draw();
	cPhiRel->SaveAs("David_Plots/PhiRelEarly.png");

	// Δφ vs multiplicity
	TCanvas *cPhiRel2D = new TCanvas();
	hPhiRelVsMultEarly->Draw("COLZ");
	cPhiRel2D->SaveAs("David_Plots/PhiRelVsMultEarly.png");

	// Ψ₂ vs multiplicity
	TCanvas* cPsi2VsN = new TCanvas();
	hPsi2VsMultEarly->Draw("COLZ");
	cPsi2VsN->SaveAs("David_Plots/Psi2VsMultEarly.png");
	*/

	// Eccentricity vs tau distribution
	//TCanvas* cEccT = new TCanvas();
	//hEccVsTau->GetYaxis()->SetRangeUser(0, 10);  // Set Y-axis range from 0 to 0.5
	//hEccVsTau->Draw("COLZ");
	//cEccT->SaveAs("David_Plots/EccentricityVsTau.png");

    // RMS Eccentricity vs tau
	TCanvas* cEccRMS = new TCanvas("cEccRMS", "RMS Eccentricity vs Tau");
	hEccRMSvsTau->SetLineWidth(2);
	hEccRMSvsTau->Draw("HIST");
	cEccRMS->SaveAs("David_Plots/EccentricityRMSvsTau.png");

    // v2 distribution at tau tauTarget
	TCanvas* cV2 = new TCanvas("cV2", "v2 at tau = 2 fm/c");
	hV2Early->SetLineWidth(2);
	hV2Early->SetLineColor(kBlue);
	hV2Early->Draw();
	cV2->SaveAs("David_Plots/V2Early.png");

    // RMS v2 vs tau
    TCanvas* cV2RMS = new TCanvas("cV2RMS", "RMS v2 vs Tau");
    hV2RMSvsTau->SetLineWidth(2);
    hV2RMSvsTau->Draw("HIST");
    cV2RMS->SaveAs("David_Plots/V2RMSvsTau.png");

    // v2 vs eccentricity at tau tauTarget scatterplot
    TCanvas* cV2vsEcc = new TCanvas("cV2vsEcc", "v2 vs Eccentricity");
    hV2vsEcc->Draw("COLZ");
    cV2vsEcc->SaveAs("David_Plots/V2vsEcc.png");

	inFile->Close();
	delete tree;
}