#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TH3F.h"
#include "Cluster/binning.h"

int drawPartonJetBin() {
    // 打开 ROOT 文件
    TFile *file = TFile::Open("parton_jet_bins.root", "READ");
    if (!file || file->IsZombie()) {
        printf("Error: 无法打开文件 parton_jet.root\n");
        return 1;
    }

    // 读取 TH2F 对象
    TH2F *hMultTime[trackbin];
    TH2F *hEccTime[trackbin];
    TH2F *hDimTime[trackbin];
    TH3F *hDenTime[trackbin];
    TH3F *hNpTime[trackbin];
    TH3F *hAreaTime[trackbin];

    for (int ibin = 0; ibin < trackbin; ibin++)
    {
        hMultTime[ibin] = (TH2F *)file->Get(Form("hMultTime_%d", ibin));
        hEccTime[ibin] = (TH2F *)file->Get(Form("hEccTime_%d", ibin));
        hDimTime[ibin] = (TH2F *)file->Get(Form("hDimTime_%d", ibin));
        hDenTime[ibin] = (TH3F *)file->Get(Form("hDenTime_%d", ibin));
        hNpTime[ibin] = (TH3F *)file->Get(Form("hNpTime_%d", ibin));
        hAreaTime[ibin] = (TH3F *)file->Get(Form("hAreaTime_%d", ibin));

        //----------------------------------------------------------------------
        // 绘制 hMultTime 并保存
        //----------------------------------------------------------------------
        if (hMultTime[ibin])
        {
            TCanvas *cMultTime = new TCanvas("cMultTime", "Multiplicity vs Time", 800, 600);
            gStyle->SetOptStat(0);
            hMultTime[ibin]->Draw("COLZ");
            cMultTime->SaveAs(Form("Cluster/Zeit/jetBins/JetPartonMultiplicityTime_%d.png", ibin));
            delete cMultTime;

            TProfile *pMultTime = hMultTime[ibin]->ProfileX("pMultTime", 1, hMultTime[ibin]->GetNbinsY());
            if (pMultTime)
            {
                pMultTime->SetMarkerStyle(20);
                pMultTime->SetLineWidth(2);
                pMultTime->SetLineColor(1);
                pMultTime->GetXaxis()->SetTitle("t (fm/c)");
                pMultTime->GetYaxis()->SetTitle("<N_{ch}^{j}>");
                TCanvas *c1 = new TCanvas("c1", "Multiplicity vs Time", 800, 600);
                pMultTime->Draw("E");
                c1->SaveAs(Form("Cluster/Zeit/jetBins/JetPartonMultiplicityTime_Profile_%d.png", ibin));
                delete c1;
                delete pMultTime;
            }
        }

        //----------------------------------------------------------------------
        // 绘制 hEccTime 并保存
        //----------------------------------------------------------------------     
        if (hEccTime[ibin])
        {
            TCanvas *cEccTime = new TCanvas("cEccTime", "Eccentricity vs Time", 800, 600);
            gStyle->SetOptStat(0);
            hEccTime[ibin]->Draw("COLZ");
            cEccTime->SaveAs(Form("Cluster/Zeit/jetBins/EccentricityTime_%d.png", ibin));
            delete cEccTime;
        }

        //----------------------------------------------------------------------
        // 从 hEccTime 创建 TProfile 并保存
        //----------------------------------------------------------------------
        if (hEccTime[ibin])
        {
            TProfile *pEccTime = hEccTime[ibin]->ProfileX("pEccTime", 1, hEccTime[ibin]->GetNbinsY());
            if (pEccTime)
            {
                pEccTime->SetMarkerStyle(20);
                pEccTime->SetLineWidth(2);
                pEccTime->SetLineColor(1);
                pEccTime->GetXaxis()->SetTitle("t (fm/c)");
                pEccTime->GetYaxis()->SetTitle("<#epsilon_{2}>");
                TCanvas *c2 = new TCanvas("c2", "Eccentricity Profile vs Time", 800, 600);
                pEccTime->Draw("E");
                c2->SaveAs(Form("Cluster/Zeit/jetBins/EccentricityTime_Profile_%d.png", ibin));
                delete c2;
                delete pEccTime;
            }
        }

        //----------------------------------------------------------------------
        // 绘制 hDimTime 并保存
        //----------------------------------------------------------------------
        if (hDimTime[ibin])
        {
            TCanvas *cDimTime = new TCanvas("cDimTime", "Dimensionality vs Time", 800, 600);
            gStyle->SetOptStat(0);
            hDimTime[ibin]->Draw("COLZ");
            cDimTime->SaveAs(Form("Cluster/Zeit/jetBins/DimensionTime_%d.png", ibin));
            delete cDimTime;
        }

        //----------------------------------------------------------------------
        // 从 hDimTime 创建 TProfile 并保存
        //----------------------------------------------------------------------
        if (hDimTime[ibin])
        {
            TProfile *pDimTime = hDimTime[ibin]->ProfileX("pDimTime", 1, hDimTime[ibin]->GetNbinsY());
            if (pDimTime)
            {
                pDimTime->SetMarkerStyle(20);
                pDimTime->SetLineWidth(2);
                pDimTime->SetLineColor(1);
                pDimTime->GetXaxis()->SetTitle("t (fm/c)");
                pDimTime->GetYaxis()->SetTitle("<r> (fm)");
                TCanvas *c3 = new TCanvas("c3", "Transverse size vs Time", 800, 600);
                pDimTime->Draw("E");
                c3->SaveAs(Form("Cluster/Zeit/jetBins/DimensionTime_Profile_%d.png", ibin));
                delete c3;
                delete pDimTime;
            }
        }

        if (hDenTime[ibin])
        {
            /*
            TCanvas *cDenTime = new TCanvas("cDenTime", "Density vs Time vs #eta_{s}", 800, 600);
            gStyle->SetOptStat(0);
            hDenTime[ibin]->Draw("box");
            cDenTime->SaveAs(Form("Cluster/Zeit/jetBins/DensityTime_%d.png", ibin));
            delete cDenTime;
            */

            TProfile2D *pDenTime = hDenTime[ibin]->Project3DProfile("yx");
            if (pDenTime)
            {
                pDenTime->SetTitle(Form("Density vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; dN/#eta_{s}/S (fm^{-2})", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]));
                TCanvas *c4 = new TCanvas("c4", "Density vs Time vs #eta_{s}", 800, 600);
                pDenTime->Draw("colz");
                c4->SaveAs(Form("Cluster/Zeit/jetBins/DensityTime_Profile_%d.png", ibin));
                delete c4;
                delete pDenTime;
            }
        }

        TProfile2D *pNpTime = hNpTime[ibin]->Project3DProfile("yx");
        pNpTime->SetTitle(Form("Parton number vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; # parton", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]));
        TCanvas *c5 = new TCanvas("c5", "", 800, 600);
        pNpTime->Draw("colz");
        c5->SaveAs(Form("Cluster/Zeit/jetBins/NPartonTime_Profile_%d.png", ibin));
        delete c5;
        delete pNpTime;

        TProfile2D *pAreaTime = hAreaTime[ibin]->Project3DProfile("yx");
        pAreaTime->SetTitle(Form("Area vs time vs #eta_{s}, %d < N_{ch}^{j} < %d; t (fm/c); #eta_{s}; S (fm^{2})", trackbinbounds_MC[ibin], trackbinboundsUpper_MC[ibin]));
        TCanvas *c6 = new TCanvas("c6", "", 800, 600);
        pAreaTime->Draw("colz");
        c6->SaveAs(Form("Cluster/Zeit/jetBins/AreaTime_Profile_%d.png", ibin));
        delete c6;
        delete pAreaTime;
    }

    // 关闭文件
    file->Close();
    delete file;

    return 0;
}
