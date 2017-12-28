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
#include <TCanvas.h>

using std::string;

void sumHist()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("test_LowPtSUSY_Tree_Triboson_WWW.root");
  TH1F *h1 = (TH1F*)file1->Get("h_WJ_PT_ElEl");
  TH1F *h2 = (TH1F*)file1->Get("h_WJ_PT_MuMu");
  TH1F *h3 = (TH1F*)file1->Get("h_WJ_PT_ElMu");

  h1->Add(h2);
  h1->Add(h3);

  h1->Rebin(20.0);
  h1->Scale(0.1124*0.67*0.22*0.22*2*100*1000.0/10000.0);
  h1->SetMinimum(10.0);
 
  h1->Draw("HIST");
  h1->GetYaxis()->SetTitle("Events/20 GeV");
  h1->GetXaxis()->SetTitle("W_{Had} p_{T} [GeV]");
  h1->SetTitle("");
  //h1->SetStats(kFALSE);
  h1->SetMinimum(0.1);

  std::cout << "h1->Integral() = " << h1->Integral() << std::endl;
  c1.SetLogy();
  c1.SaveAs("WHad_pT_AllFlavor_Triboson.pdf");
  c1.SaveAs("WHad_pT_AllFlavor_Triboson.png");

} 
