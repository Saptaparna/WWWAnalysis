#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
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
#include <TMath.h>
#include <TH2D.h>

using std::string;

void makeComparisonPlot()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("higgsCombineTest.MultiDimFit.mH120_onlySS.root");
  TTree* limit = (TTree*) file1->Get("limit"); 
  long int n = limit->Draw("2*deltaNLL:ft0", "", "goff", 100000, 0);
  std::cout << "n = " << n << std::endl;
  TGraph *g = new TGraph(n,limit->GetV2(),limit->GetV1());
  g->SetLineColor(kBlue);
  g->SetMarkerColor(kBlue);
  g->SetLineWidth(2);
  g->SetMaximum(18);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("f_{T0}");
  g->GetYaxis()->SetTitle("Test statistic (-2ln#lambda )");
  g->Draw("ap");  

  TFile* file3 = TFile::Open("higgsCombineTest.MultiDimFit.mH120_3l_ft0.root");
  TTree* limit3 = (TTree*) file3->Get("limit");  
  long int n3 = limit3->Draw("2*deltaNLL:ft0", "", "goff", 100000, 0);
  std::cout << "n3 = " << n3 << std::endl;
  TGraph *g3 = new TGraph(n3,limit3->GetV2(),limit3->GetV1());
  g3->SetLineColor(kBlack);
  g3->SetMarkerColor(kBlack);
  g3->SetLineWidth(2);
  g3->Draw("p same");

  TFile* file2 = TFile::Open("higgsCombineTest.MultiDimFit.mH120_SS_3l.root");
  TTree* limit2 = (TTree*) file2->Get("limit");
  long int n2 = limit2->Draw("2*deltaNLL:ft0", "", "goff", 100000, 0);
  std::cout << "n2 = " << n2 << std::endl;
  TGraph *g2 = new TGraph(n2,limit2->GetV2(),limit2->GetV1());
  g2->SetLineColor(kRed);
  g2->SetMarkerColor(kRed);
  g2->SetLineWidth(2);
  g2->Draw("p same");

  TFile* file4 = TFile::Open("ft0_data/higgsCombineTest.MultiDimFit.mH120.root");
  TTree* limit4 = (TTree*) file4->Get("limit");
  long int n4 = limit4->Draw("2*deltaNLL:ft0", "", "goff", 100000, 0);
  std::cout << "n4 = " << n4 << std::endl;
  TGraph *g4 = new TGraph(n4, limit4->GetV2(), limit4->GetV1());
  g4->SetLineColor(kCyan);
  g4->SetMarkerColor(kCyan);
  g4->SetLineWidth(2);
  g4->Draw("p same");

  TLegend *leg = new TLegend(0.45,0.70,0.65,0.85,NULL,"brNDC");
  leg->AddEntry(g, "Dilepton only", "pl");
  leg->AddEntry(g3, "Trilepton only", "pl");
  leg->AddEntry(g2, "Combined", "pl");
  leg->AddEntry(g4, "Observed", "pl");
  leg->SetBorderSize(0);
  leg->Draw();

  c1.SaveAs("ft0_comparison_SS_3l_data.png");
  c1.SaveAs("ft0_comparison_SS_3l_data.pdf");


} 
