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

void sumVHadpT()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("test_LowPtSUSY_Tree_Diboson_MJCut.root"); 
     
  TH1F *h1 = (TH1F*)file1->Get("h_WJ_PT");
  h1->Rebin(40.0);
  h1->GetXaxis()->SetRangeUser(200, 10000);
  h1->Scale(190.94);    
  
  TFile* file2 = TFile::Open("test_LowPtSUSY_Tree_Diboson_WZ_MJCut.root");

  TH1F *h2 = (TH1F*)file2->Get("h_WJ_PT");
  h2->Rebin(40.0);
  h2->GetXaxis()->SetRangeUser(200, 10000);
  h2->Scale(72.99);
  
  h2->Add(h1);
  
  h2->Draw("");
  h2->GetYaxis()->SetTitle("Events/40 GeV");
  h2->GetXaxis()->SetTitle("V_{Had} p_{T} GeV");
  h2->SetTitle("");
  h2->SetMinimum(0.1); 

  c1.SetLogy(); 
  c1.SaveAs("Vhad_pT_MJCut.pdf");
  c1.SaveAs("Vhad_pT_MJCut.png");
}

void WpT()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("test_LowPtSUSY_Tree_Diboson_MJCut.root");

  TH1F *h1 = (TH1F*)file1->Get("h_WJ_PT");
  h1->Rebin(40.0);
  h1->GetXaxis()->SetRangeUser(200, 10000);
  h1->Scale(190.94);
  h1->SetMinimum(0.1);

  h1->Draw("");
  h1->GetYaxis()->SetTitle("Events/40 GeV");
  h1->GetXaxis()->SetTitle("W_{Had} p_{T} [GeV]");
  h1->SetTitle("");
  h1->SetMinimum(0.1);

  c1.SetLogy();
  c1.SaveAs("WHad_pT_MJCut.pdf");
  c1.SaveAs("WHad_pT_MJCut.png");

} 

void ZpT()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file2 = TFile::Open("test_LowPtSUSY_Tree_Diboson_WZ_MJCut.root");

  TH1F *h2 = (TH1F*)file2->Get("h_WJ_PT");
  h2->Rebin(40.0);
  h2->GetXaxis()->SetRangeUser(200, 10000);
  h2->Scale(72.99);
  h2->SetMinimum(0.1);

  h2->Draw("");
  h2->GetYaxis()->SetTitle("Events/40 GeV");
  h2->GetXaxis()->SetTitle("Z_{Had} p_{T} [GeV]");
  h2->SetTitle("");
  h2->SetMinimum(0.1);

  c1.SetLogy();
  c1.SaveAs("ZHad_pT_MJCut.pdf");
  c1.SaveAs("ZHad_pT_MJCut.png");

}

void ZpTLeptonic()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file3 = TFile::Open("test_LowPtSUSY_Tree_Diboson_WZ_Z_leptonic.root");

  TH1F *h3 = (TH1F*)file3->Get("h_WJ_PT");
  h3->Rebin(40.0);
  h3->GetXaxis()->SetRangeUser(200, 10000);
  h3->Scale(19.90);
  h3->SetMinimum(0.1);

  h3->Draw("");
  h3->GetYaxis()->SetTitle("Events/40 GeV");
  h3->GetXaxis()->SetTitle("W_{Had} p_{T} [GeV]");
  h3->SetTitle("W_{Had} p_{T}, Z #rightarrow l+ l-");
  h3->SetMinimum(0.1);

  c1.SetLogy();
  c1.SaveAs("WHad_Zleptonic_pT.pdf");
  c1.SaveAs("WHad_Zleptonic_pT.png");

}
 
void sumVHadpTWith8TeVComparison()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  TFile* file1 = TFile::Open("test_LowPtSUSY_Tree_Diboson_MJCut.root");

  TH1F *h1 = (TH1F*)file1->Get("h_WJ_PT");
  h1->Rebin(40.0);
  h1->GetXaxis()->SetRangeUser(200, 10000);
  h1->Scale(190.94);

  TFile* file2 = TFile::Open("test_LowPtSUSY_Tree_Diboson_WZ_MJCut.root");

  TH1F *h2 = (TH1F*)file2->Get("h_WJ_PT");
  h2->Rebin(40.0);
  h2->GetXaxis()->SetRangeUser(200, 10000);
  h2->Scale(72.99);

  h2->Add(h1);

  h2->Draw("");
  h2->GetYaxis()->SetTitle("Events/40 GeV");
  h2->GetXaxis()->SetTitle("V_{Had} p_{T} GeV");
  h2->SetTitle("");
  h2->SetStats(kFALSE);
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);
  h2->SetMinimum(0.1);

  TH1F *h3 = new TH1F("h3", "h3", 15, 200.0, 800.0);
  h3->SetBinContent(1, 1100.0);
  h3->SetBinContent(3, 500);
  h3->SetBinContent(6, 60);
  h3->SetBinContent(8, 30);
  h3->SetBinContent(11, 8);
  h3->SetMinimum(0.1);
  h3->SetLineColor(kBlack);
  h3->SetLineWidth(2);
  h3->Draw("pe same");

  c1.SetLogy();
  c1.SaveAs("Vhad_pT_MJCut.pdf");
  c1.SaveAs("Vhad_pT_MJCut.png");
}

void overSigBkg()
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
  h1->SetLineColor(kBlue);
  h1->SetLineWidth(2);

  h1->Draw("HIST");
  h1->GetYaxis()->SetTitle("Events/20 GeV");
  h1->GetXaxis()->SetTitle("W_{Had} p_{T} [GeV]");
  h1->SetTitle("");
  h1->SetStats(kFALSE);
  h1->SetMinimum(0.1);

  TFile* file2 = TFile::Open("test_LowPtSUSY_Tree_WW_JJ.root");
  TH1F *h4 = (TH1F*)file2->Get("h_jet_pt_leading");
  h4->SetLineColor(kBlack);
  h4->SetLineWidth(2);
  h4->Rebin(20.0);
  h4->Scale(0.01289*1000*100/10000);
  h4->SetStats(kFALSE);
  h4->Draw("HIST SAME");


  std::cout << "h1->Integral() = " << h1->Integral() << std::endl;
  std::cout << "h4->Integral() = " << h4->Integral() << std::endl; 
  
  std::cout << "h1->GetMean() = " << h1->GetMean() << std::endl;
  std::cout << "h4->GetMean() = " << h4->GetMean() << std::endl;
  c1.SetLogy();
  c1.SaveAs("WHad_pT_AllFlavor_Triboson_OverLay_Bkg.pdf");
  c1.SaveAs("WHad_pT_AllFlavor_Triboson_OverLay_Bkg.png");



}

