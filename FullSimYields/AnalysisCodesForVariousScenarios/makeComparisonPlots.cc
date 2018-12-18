#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLegend.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double LUMI = 35.9;
double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

void makeComparisonPlots_signal_ft()
{

  TCanvas c1("c1","Stacked Histogram",10,10,700,800);
  //c1.SetLogx();
  c1.SetLogy();

  TFile* fileWWW_run_01_ft_neg2_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_01_ft_neg2_5 = (TH1F*)fileWWW_run_01_ft_neg2_5->Get("h_ST_rwgt67");
  hWWW_run_01_ft_neg2_5->SetLineColor(kBlack);
  hWWW_run_01_ft_neg2_5->SetLineWidth(2);
  hWWW_run_01_ft_neg2_5->Rebin(1);
  hWWW_run_01_ft_neg2_5->SetTitle("");
  hWWW_run_01_ft_neg2_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_01_ft_neg2_5->Scale(LUMI*BR);
  hWWW_run_01_ft_neg2_5->SetMaximum(10);
  hWWW_run_01_ft_neg2_5->SetStats(kFALSE);
  hWWW_run_01_ft_neg2_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_01_ft_neg2_5->GetXaxis()->SetTitle("ST [GeV]");
  hWWW_run_01_ft_neg2_5->GetYaxis()->SetTitle("Events");
  hWWW_run_01_ft_neg2_5->Draw("hist");
  
  TFile* fileWWW_run_02_ft_neg2_0 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_02_ft_neg2_0 = (TH1F*)fileWWW_run_02_ft_neg2_0->Get("h_ST_rwgt68");
  hWWW_run_02_ft_neg2_0->SetLineColor(kBlue);
  hWWW_run_02_ft_neg2_0->SetLineWidth(2);
  hWWW_run_02_ft_neg2_0->Rebin(1);
  hWWW_run_02_ft_neg2_0->SetStats(kFALSE);
  hWWW_run_02_ft_neg2_0->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_02_ft_neg2_0->Scale(LUMI*BR);
  hWWW_run_02_ft_neg2_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_02_ft_neg2_0->Draw("hist same");

  TFile* fileWWW_run_03_ft_neg1_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_03_ft_neg1_5 = (TH1F*)fileWWW_run_03_ft_neg1_5->Get("h_ST_rwgt69");
  hWWW_run_03_ft_neg1_5->SetLineColor(kGreen+3);
  hWWW_run_03_ft_neg1_5->SetLineWidth(2);
  hWWW_run_03_ft_neg1_5->Rebin(1);
  hWWW_run_03_ft_neg1_5->SetStats(kFALSE);
  hWWW_run_03_ft_neg1_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_03_ft_neg1_5->Scale(LUMI*BR);
  hWWW_run_03_ft_neg1_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_03_ft_neg1_5->Draw("hist same");
 
  TFile* fileWWW_run_04_ft_neg1_0 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_04_ft_neg1_0 = (TH1F*)fileWWW_run_04_ft_neg1_0->Get("h_ST_rwgt70");
  hWWW_run_04_ft_neg1_0->SetLineColor(kCyan);
  hWWW_run_04_ft_neg1_0->SetLineWidth(2);
  hWWW_run_04_ft_neg1_0->Rebin(1);
  hWWW_run_04_ft_neg1_0->SetStats(kFALSE);
  hWWW_run_04_ft_neg1_0->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_04_ft_neg1_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_04_ft_neg1_0->Scale(LUMI*BR);
  hWWW_run_04_ft_neg1_0->Draw("hist same");
 
  TFile* fileWWW_run_05_ft_neg0_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_05_ft_neg0_5 = (TH1F*)fileWWW_run_05_ft_neg0_5->Get("h_ST_rwgt71");
  hWWW_run_05_ft_neg0_5->SetLineColor(kMagenta);
  hWWW_run_05_ft_neg0_5->SetLineWidth(2);
  hWWW_run_05_ft_neg0_5->Rebin(1);
  hWWW_run_05_ft_neg0_5->SetStats(kFALSE);
  hWWW_run_05_ft_neg0_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_05_ft_neg0_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_05_ft_neg0_5->Scale(LUMI*BR);
  hWWW_run_05_ft_neg0_5->Draw("hist same");
 
  TFile* fileWWW_run_11_ft_0_0 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_11_ft_0_0 = (TH1F*)fileWWW_run_11_ft_0_0->Get("h_ST_rwgt72");
  hWWW_run_11_ft_0_0->SetLineColor(kViolet+1);
  hWWW_run_11_ft_0_0->SetLineWidth(2);
  hWWW_run_11_ft_0_0->Rebin(1);
  hWWW_run_11_ft_0_0->SetLineStyle(2);
  hWWW_run_11_ft_0_0->SetStats(kFALSE);
  hWWW_run_11_ft_0_0->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_11_ft_0_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_11_ft_0_0->Scale(LUMI*BR);
  hWWW_run_11_ft_0_0->Draw("hist same");
 
  TFile* fileWWW_run_06_ft_pos0_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_06_ft_pos0_5 = (TH1F*)fileWWW_run_06_ft_pos0_5->Get("h_ST_rwgt73");
  hWWW_run_06_ft_pos0_5->SetLineColor(kOrange+3);
  hWWW_run_06_ft_pos0_5->SetLineWidth(2);
  hWWW_run_06_ft_pos0_5->Rebin(1);
  hWWW_run_06_ft_pos0_5->SetStats(kFALSE);
  hWWW_run_06_ft_pos0_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_06_ft_pos0_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_06_ft_pos0_5->Scale(LUMI*BR);
  hWWW_run_06_ft_pos0_5->Draw("hist same");

  TFile* fileWWW_run_07_ft_pos1_0 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_07_ft_pos1_0 = (TH1F*)fileWWW_run_07_ft_pos1_0->Get("h_ST_rwgt74");
  hWWW_run_07_ft_pos1_0->SetLineColor(kRed);
  hWWW_run_07_ft_pos1_0->SetLineWidth(2);
  hWWW_run_07_ft_pos1_0->Rebin(1);
  hWWW_run_07_ft_pos1_0->SetStats(kFALSE);
  hWWW_run_07_ft_pos1_0->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_07_ft_pos1_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_07_ft_pos1_0->Scale(LUMI*BR);
  hWWW_run_07_ft_pos1_0->Draw("hist same");

  TFile* fileWWW_run_08_ft_pos1_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_08_ft_pos1_5 = (TH1F*)fileWWW_run_08_ft_pos1_5->Get("h_ST_rwgt75");
  hWWW_run_08_ft_pos1_5->SetLineColor(kYellow);
  hWWW_run_08_ft_pos1_5->SetLineWidth(2);
  hWWW_run_08_ft_pos1_5->Rebin(1);
  hWWW_run_08_ft_pos1_5->SetStats(kFALSE);
  hWWW_run_08_ft_pos1_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_08_ft_pos1_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_08_ft_pos1_5->Scale(LUMI*BR);
  hWWW_run_08_ft_pos1_5->Draw("hist same"); 

  TFile* fileWWW_run_09_ft_pos2_0 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_09_ft_pos2_0 = (TH1F*)fileWWW_run_09_ft_pos2_0->Get("h_ST_rwgt76");
  hWWW_run_09_ft_pos2_0->SetLineColor(kGreen-3);
  hWWW_run_09_ft_pos2_0->SetLineWidth(2);
  hWWW_run_09_ft_pos2_0->Rebin(1);
  hWWW_run_09_ft_pos2_0->SetStats(kFALSE);
  hWWW_run_09_ft_pos2_0->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_09_ft_pos2_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_09_ft_pos2_0->Scale(LUMI*BR);
  hWWW_run_09_ft_pos2_0->Draw("hist same");

  TFile* fileWWW_run_10_ft_pos2_5 = new TFile("output_www_aqgc_sapta_skim_1_1.root");
  TH1F *hWWW_run_10_ft_pos2_5 = (TH1F*)fileWWW_run_10_ft_pos2_5->Get("h_ST_rwgt77");
  hWWW_run_10_ft_pos2_5->SetLineColor(kBlue-7);
  hWWW_run_10_ft_pos2_5->SetLineWidth(2);
  hWWW_run_10_ft_pos2_5->Rebin(1);
  hWWW_run_10_ft_pos2_5->SetStats(kFALSE);
  hWWW_run_10_ft_pos2_5->GetXaxis()->SetRangeUser(50, 3000.0);
  hWWW_run_10_ft_pos2_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_10_ft_pos2_5->Scale(LUMI*BR);
  hWWW_run_10_ft_pos2_5->Draw("hist same");

  TLegend *leg1 = new TLegend(0.55,0.58,0.85,0.88,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->SetHeader("ft [TeV^{-4}]");
  leg1->AddEntry(hWWW_run_01_ft_neg2_5, "-2_5", "l");
  leg1->AddEntry(hWWW_run_02_ft_neg2_0, "-2_0", "l");
  leg1->AddEntry(hWWW_run_03_ft_neg1_5, "-1_5", "l");
  leg1->AddEntry(hWWW_run_04_ft_neg1_0, "-1_0", "l");
  leg1->AddEntry(hWWW_run_05_ft_neg0_5, "-0_5", "l");
  leg1->AddEntry(hWWW_run_06_ft_pos0_5, " 0_5", "l");
  leg1->AddEntry(hWWW_run_07_ft_pos1_0, " 1_0", "l");
  leg1->AddEntry(hWWW_run_08_ft_pos1_5, " 1_5", "l");
  leg1->AddEntry(hWWW_run_09_ft_pos2_0, " 2_0", "l");
  leg1->AddEntry(hWWW_run_10_ft_pos2_5, " 2_5", "l");
  leg1->AddEntry(hWWW_run_11_ft_0_0, "0 (SM)", "l");
  leg1->Draw();

  c1.SaveAs("AC_ST.png");
  c1.SaveAs("AC_ST.pdf");

}
