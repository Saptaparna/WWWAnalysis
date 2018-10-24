#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TCanvas.h>
using std::string;

void makeGraphs_pValue_All()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double p_value_noST[11] = {0.144235, 0.240737, 0.342305, 0.427696, 0.482399, 0.5, 0.478365, 0.420009, 0.331509, 0.229094, 0.13437};

    double p_value_ST1000[11] = {6.04863e-07, 0.000176481, 0.00928318, 0.105986, 0.35638, 0.5, 0.342286, 0.096894, 0.00777221, 0.000133828, 4.35856e-07};

    double p_value_ST1500[11] = {5.34307e-11, 3.71908e-07, 0.000256515, 0.0188431, 0.210674, 0.5, 0.203202, 0.0175715, 0.000220472, 3.16381e-07, 4.26177e-11};

    double p_value_ST2000[11] = {4.5312e-09, 4.68642e-06, 0.000814214, 0.0273138, 0.214777, 0.5, 0.209442, 0.0260321, 0.000764629, 4.10714e-06, 3.94965e-09}; 

    double p_value_ST2500[11] = {1.08031e-05, 0.000494848, 0.00905416, 0.0710233, 0.264023, 0.499999, 0.259189, 0.0687127, 0.00862136, 0.000462987, 9.92006e-06}; 

    double p_value_disc[11] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(11, param, p_value_noST);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param, p_value_ST1000);
    gr_1->SetMarkerStyle(20);
    gr_1->SetMarkerColor(kCyan);
    gr_1->SetLineColor(kCyan);
    gr_1->SetMarkerSize(0.9);
    gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(11, param, p_value_ST1500);
    gr_2->SetMarkerStyle(20);
    gr_2->SetMarkerColor(kBlack);
    gr_2->SetLineColor(kBlack);
    gr_2->SetMarkerSize(0.9);
    gr_2->Draw("LP*");

    TGraph *gr_3 = new TGraph(11, param, p_value_ST2000);
    gr_3->SetMarkerStyle(20);
    gr_3->SetMarkerColor(kGreen+3);
    gr_3->SetLineColor(kGreen+3);
    gr_3->SetMarkerSize(0.9);
    gr_3->Draw("LP*");

    TGraph *gr_4 = new TGraph(11, param, p_value_ST2500);
    gr_4->SetMarkerStyle(20);
    gr_4->SetMarkerColor(kMagenta);
    gr_4->SetLineColor(kMagenta);
    gr_4->SetMarkerSize(0.9);
    gr_4->Draw("LP*");

    TGraph *gr_5 = new TGraph(11, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.75,0.70,0.89,0.89,NULL,"brNDC");
    leg->AddEntry(gr, "no ST cut", "pl");
    leg->AddEntry(gr_1, "ST 1000", "pl");
    leg->AddEntry(gr_2, "ST 1500", "pl");
    leg->AddEntry(gr_3, "ST 2000", "pl");
    leg->AddEntry(gr_4, "ST 2500", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr_1->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit ST 1000 = " << m << endl;
    
    while (m<=2.0 and gr_1->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit ST 1000 = " << m << endl;
   
    double m1=-2.0; 
    while (m1<=2.0 and gr_2->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit ST 1500 = " << m1 << endl;
    
    while (m1<=2.0 and gr_2->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit ST 1500 = " << m1 << endl;

    c1.SaveAs("Parabola_ft0_p_value_All.pdf");
    c1.SaveAs("Parabola_ft0_p_value_All.png");
    
}

void makeGraphs_pValue_combined()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
   
    double p_value_ST1500[11] = {3.09056e-14, 5.27986e-09, 3.67597e-05, 0.0102541, 0.194951, 0.5, 0.18665, 0.00931827, 3.00175e-05, 4.14089e-09, 2.16566e-14};
    //double p_value_ST1500[11] = {1.43907e-08, 2.1064e-05, 0.00360288, 0.0794799, 0.341577, 0.5, 0.326363, 0.0711897, 0.00287657, 1.50716e-05, 9.50438e-09};
    
    double p_value_disc[11] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(11, param, p_value_ST1500);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000000000000000000001);
    //gr->GetXaxis()->SetLimits(-1.2, 1.2);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_5 = new TGraph(11, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.75,0.70,0.89,0.89,NULL,"brNDC");
    leg->AddEntry(gr, "ST 1500", "pl");
    //leg->AddEntry(gr_1, "ST 1000", "pl");
    //leg->AddEntry(gr_2, "ST 1500", "pl");
    //leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit ST 1500 = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit ST 1500 = " << m << endl;
    
    
    c1.SaveAs("Parabola_ft0_p_value_All.pdf");
    c1.SaveAs("Parabola_ft0_p_value_All.png");
    
}

