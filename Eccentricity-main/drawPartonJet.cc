#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TStyle.h>

int drawPartonJet() {
    // 打开 ROOT 文件
    TFile *file = TFile::Open("parton_jet.root", "READ");
    if (!file || file->IsZombie()) {
        printf("Error: 无法打开文件 parton_jet.root\n");
        return 1;
    }

    // 读取 TH2F 对象
    TH2F *hMultTime = (TH2F*)file->Get("hMultTime");
    TH2F *hEccTime  = (TH2F*)file->Get("hEccTime");
    TH2F *hDimTime  = (TH2F*)file->Get("hDimTime");

    if (!hMultTime) {
        printf("Error: 无法找到 hMultTime\n");
    }
    if (!hEccTime) {
        printf("Error: 无法找到 hEccTime\n");
    }
    if (!hDimTime) {
        printf("Error: 无法找到 hDimTime\n");
    }

    //----------------------------------------------------------------------
    // 绘制 hMultTime 并保存
    //----------------------------------------------------------------------
    if (hMultTime) {
        TCanvas *cMultTime = new TCanvas("cMultTime", "Multiplicity vs Time", 800, 600);
        gStyle->SetOptStat(0);
        hMultTime->Draw("COLZ");
        cMultTime->SaveAs("Cluster/Zeit/JetPartonMultiplicityTime.png");
        delete cMultTime;
    }

    //----------------------------------------------------------------------
    // 从 hMultTime 创建 TProfile 并保存
    //----------------------------------------------------------------------
    if (hMultTime) {
        TProfile *pMultTime = hMultTime->ProfileX("pMultTime", 1, hMultTime->GetNbinsY());
        if (pMultTime) {
            pMultTime->GetXaxis()->SetTitle("t (fm/c)");
            pMultTime->GetYaxis()->SetTitle("<N_{ch}^{j}>");
            TCanvas *c1 = new TCanvas("c1", "Multiplicity vs Time", 800, 600);
            pMultTime->Draw("E");
            c1->SaveAs("Cluster/Zeit/JetPartonMultiplicityTime_Profile.png");
            delete c1;
            delete pMultTime;
        }
    }

    //----------------------------------------------------------------------
    // 绘制 hEccTime 并保存
    //----------------------------------------------------------------------
    if (hEccTime) {
        TCanvas *cEccTime = new TCanvas("cEccTime", "Eccentricity vs Time", 800, 600);
        gStyle->SetOptStat(0);
        hEccTime->Draw("COLZ");
        cEccTime->SaveAs("Cluster/Zeit/EccentricityTime.png");
        delete cEccTime;
    }

    //----------------------------------------------------------------------
    // 从 hEccTime 创建 TProfile 并保存
    //----------------------------------------------------------------------
    if (hEccTime) {
        TProfile *pEccTime = hEccTime->ProfileX("pEccTime", 1, hEccTime->GetNbinsY());
        if (pEccTime) {
            pEccTime->GetXaxis()->SetTitle("t (fm/c)");
            pEccTime->GetYaxis()->SetTitle("<#epsilon_{2}>");
            TCanvas *c2 = new TCanvas("c2", "Eccentricity vs Time", 800, 600);
            pEccTime->Draw("E");
            c2->SaveAs("Cluster/Zeit/EccentricityTime_Profile.png");
            delete c2;
            delete pEccTime;
        }
    }

    //----------------------------------------------------------------------
    // 绘制 hDimTime 并保存
    //----------------------------------------------------------------------
    if (hDimTime) {
        TCanvas *cDimTime = new TCanvas("cDimTime", "Dimensionality vs Time", 800, 600);
        gStyle->SetOptStat(0);
        hDimTime->Draw("COLZ");
        cDimTime->SaveAs("Cluster/Zeit/DimensionTime.png");
        delete cDimTime;
    }

    //----------------------------------------------------------------------
    // 从 hDimTime 创建 TProfile 并保存
    //----------------------------------------------------------------------
    if (hDimTime) {
        TProfile *pDimTime = hDimTime->ProfileX("pDimTime", 1, hDimTime->GetNbinsY());
        if (pDimTime) {
            pDimTime->GetXaxis()->SetTitle("t (fm/c)");
            pDimTime->GetYaxis()->SetTitle("<r> (fm)");
            TCanvas *c3 = new TCanvas("c3", "Transverse size vs Time", 800, 600);
            pDimTime->Draw("E");
            c3->SaveAs("Cluster/Zeit/DimensionTime_Profile.png");
            delete c3;
            delete pDimTime;
        }
    }

    // 关闭文件
    file->Close();
    delete file;

    return 0;
}
