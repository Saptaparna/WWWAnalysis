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

void makeCombined()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("higgsCombineTest.MultiDimFit.mH120_SS_3l.root");
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

  TLegend *leg = new TLegend(0.45,0.70,0.65,0.85,NULL,"brNDC");
  leg->AddEntry(g, "Combined", "pl");
  leg->SetBorderSize(0);
  leg->Draw();

  c1.SaveAs("ft0_combined.png");
  c1.SaveAs("ft0_combined.pdf");


} 
