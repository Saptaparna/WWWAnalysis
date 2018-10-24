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


void makeGraphs_noST()
{
  gROOT->SetStyle("Plain");
  TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

  double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

  /*double yield1 = 8.98;
  double yield2 = 6.44;
  double yield3 = 4.47;
  double yield4 = 3.07;
  double yield5 = 2.24;
  double yield6 = 1.98;
  double yield7 = 2.3;
  double yield8 = 3.19;
  double yield9 = 4.66;
  double yield10 = 6.7;
  double yield11 = 9.31;*/

  double yield1 = 8.98-1.98;
  double yield2 = 6.44-1.98;
  double yield3 = 4.47-1.98;
  double yield4 = 3.07-1.98;
  double yield5 = 2.24-1.98;
  double yield6 = 1.98-1.98;
  double yield7 = 2.3-1.98;
  double yield8 = 3.19-1.98;
  double yield9 = 4.66-1.98;
  double yield10 = 6.7-1.98;
  double yield11 = 9.31-1.98;

  double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

  TGraph *gr = new TGraph(11, param, rawYield);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(0.9);
  gr->SetMaximum(10);
  gr->SetMinimum(0.0);
  gr->GetXaxis()->SetLimits(-3.0, 3.0);
  gr->SetTitle("");
  gr->SetTitle("Param (ft0) Vs Yield");
  gr->GetXaxis()->SetTitle("Param (ft0)");
  gr->GetYaxis()->SetTitle("Signal yield");
  gr->Draw("AP*");

  TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
  fit_yield->SetLineColor(kRed);
  fit_yield->SetLineWidth(2);
  fit_yield->SetLineStyle(9);
  gr->Fit("fit_yield", "", "", -2.5, 2.5);

  c1.SaveAs("Parabola_ft0_noST.pdf");
  c1.SaveAs("Parabola_ft0_noST.png");

  TFile *outputFile = new TFile("output_noST.root", "RECREATE");
  fit_yield->Write("fit_yield");
  outputFile->Close();

}

void makeGraphs_ST1000()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    double yield1 = 6.00;
    double yield2 = 3.84;
    double yield3 = 2.17;
    double yield4 = 0.982;
    double yield5 = 0.273;
    double yield6 = 0.0455;
    double yield7 = 0.299;
    double yield8 = 1.03;
    double yield9 = 2.25;
    double yield10 = 3.95;
    double yield11 = 6.12;
    
    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};
    
    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(10);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");
    
    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);
    
    c1.SaveAs("Parabola_ft0_ST1000.pdf");
    c1.SaveAs("Parabola_ft0_ST1000.png");
    
    TFile *outputFile = new TFile("output_ST1000.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
    
}

void makeGraphs_ST1500()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    double yield1 = 4.45;
    double yield2 = 2.85;
    double yield3 = 1.6;
    double yield4 = 0.716;
    double yield5 = 0.184;
    double yield6 = 0.00985;
    double yield7 = 0.192;
    double yield8 = 0.731;
    double yield9 = 1.63;
    double yield10 = 2.88;
    double yield11 = 4.49;
    
    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};
    
    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(10);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");
    
    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);
    
    c1.SaveAs("Parabola_ft0_ST1500.pdf");
    c1.SaveAs("Parabola_ft0_ST1500.png");
    
    TFile *outputFile = new TFile("output_ST1500.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
    
}

void makeGraphs_ST2000()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);
    
    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double yield1 = 2.86;
    double yield2 = 1.83;
    double yield3 = 1.03;
    double yield4 = 0.458;
    double yield5 = 0.116;
    double yield6 = 0.00296;
    double yield7 = 0.12;
    double yield8 = 0.466;
    double yield9 = 1.04;
    double yield10 = 1.85;
    double yield11 = 2.88;
    
    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};
    
    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(10);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");
    
    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);
    
    c1.SaveAs("Parabola_ft0_ST2000.pdf");
    c1.SaveAs("Parabola_ft0_ST2000.png");
    
    TFile *outputFile = new TFile("output_ST2000.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
    
}

void makeGraphs_ST2500()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double yield1 = 1.26; 
    double yield2 = 0.806;
    double yield3 = 0.453;
    double yield4 = 0.201;
    double yield5 = 0.05; 
    double yield6 = 0.000344;
    double yield7 = 0.0519;
    double yield8 = 0.205;
    double yield9 = 0.459;
    double yield10 = 0.814;
    double yield11 = 1.27;

    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(10);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");

    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);

    c1.SaveAs("Parabola_ft0_ST2500.pdf");
    c1.SaveAs("Parabola_ft0_ST2500.png");

    TFile *outputFile = new TFile("output_ST2500.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();

}

void makeGraphs_ST1500_0SFOS()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double yield1 = 1.15;
    double yield2 = 0.739;
    double yield3 = 0.417;
    double yield4 = 0.186;
    double yield5 = 0.0483;
    double yield6 = 0.00271;
    double yield7 = 0.0495;
    double yield8 = 0.189;
    double yield9 = 0.42;
    double yield10 = 0.744;
    double yield11 = 1.16;

    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(3.0);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");

    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);

    c1.SaveAs("Parabola_ft0_ST1500_0SFOS.pdf");
    c1.SaveAs("Parabola_ft0_ST1500_0SFOS.png");

    TFile *outputFile = new TFile("output_ST1500_0SFOS.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
}

void makeGraphs_ST1500_1SFOS()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double yield1 = 1.82;
    double yield2 = 1.16;
    double yield3 = 0.655;
    double yield4 = 0.291;
    double yield5 = 0.0744;
    double yield6 = 0.00329;
    double yield7 = 0.0782;
    double yield8 = 0.299;
    double yield9 = 0.666;
    double yield10 = 1.18;
    double yield11 = 1.84;

    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(3.0);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");

    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);

    c1.SaveAs("Parabola_ft0_ST1500_1SFOS.pdf");
    c1.SaveAs("Parabola_ft0_ST1500_1SFOS.png");

    TFile *outputFile = new TFile("output_ST1500_1SFOS.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
}

void makeGraphs_ST1500_2SFOS()
{
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

    double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

    double yield1 = 1.48;
    double yield2 = 0.946;
    double yield3 = 0.533;
    double yield4 = 0.238;
    double yield5 = 0.0618;
    double yield6 = 0.00385;
    double yield7 = 0.0643;
    double yield8 = 0.243;
    double yield9 = 0.54;
    double yield10 = 0.956;
    double yield11 = 1.49;

    double rawYield[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

    TGraph *gr = new TGraph(11, param, rawYield);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(3.0);
    gr->SetMinimum(0.0);
    gr->GetXaxis()->SetLimits(-3.0, 3.0);
    gr->SetTitle("");
    gr->SetTitle("Param (ft0) Vs Yield");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("Signal yield");
    gr->Draw("AP*");

    TF1 *fit_yield = new TF1("fit_yield","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
    fit_yield->SetLineColor(kRed);
    fit_yield->SetLineWidth(2);
    fit_yield->SetLineStyle(9);
    gr->Fit("fit_yield", "", "", -2.5, 2.5);

    c1.SaveAs("Parabola_ft0_ST1500_2SFOS.pdf");
    c1.SaveAs("Parabola_ft0_ST1500_2SFOS.png");

    TFile *outputFile = new TFile("output_ST1500_2SFOS.root", "RECREATE");
    fit_yield->Write("fit_yield");
    outputFile->Close();
}

