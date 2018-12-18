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

void makeComparisonPlots_signal_fm0()
{

  TCanvas c1("c1","Stacked Histogram",10,10,700,800);
  //c1.SetLogx();
  c1.SetLogy();

  TFile* fileWWW_run_01_fm0_neg2_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_01_fm0_neg2_5 = (TH1F*)fileWWW_run_01_fm0_neg2_5->Get("h_ST_rwgt23");
  hWWW_run_01_fm0_neg2_5->SetLineColor(kBlack);
  hWWW_run_01_fm0_neg2_5->SetLineWidth(2);
  hWWW_run_01_fm0_neg2_5->Rebin(2);
  hWWW_run_01_fm0_neg2_5->SetTitle("");
  hWWW_run_01_fm0_neg2_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_01_fm0_neg2_5->Scale(LUMI*BR);
  hWWW_run_01_fm0_neg2_5->SetMaximum(50);
  hWWW_run_01_fm0_neg2_5->SetMinimum(0.001);
  hWWW_run_01_fm0_neg2_5->SetStats(kFALSE);
  hWWW_run_01_fm0_neg2_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_01_fm0_neg2_5->GetXaxis()->SetTitle("ST [GeV]");
  hWWW_run_01_fm0_neg2_5->GetYaxis()->SetTitle("Events");
  hWWW_run_01_fm0_neg2_5->Draw("hist");
  
  TFile* fileWWW_run_02_fm0_neg2_0 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_02_fm0_neg2_0 = (TH1F*)fileWWW_run_02_fm0_neg2_0->Get("h_ST_rwgt24");
  hWWW_run_02_fm0_neg2_0->SetLineColor(kBlue);
  hWWW_run_02_fm0_neg2_0->SetLineWidth(2);
  hWWW_run_02_fm0_neg2_0->Rebin(2);
  hWWW_run_02_fm0_neg2_0->SetStats(kFALSE);
  hWWW_run_02_fm0_neg2_0->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_02_fm0_neg2_0->Scale(LUMI*BR);
  hWWW_run_02_fm0_neg2_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_02_fm0_neg2_0->Draw("hist same");

  TFile* fileWWW_run_03_fm0_neg1_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_03_fm0_neg1_5 = (TH1F*)fileWWW_run_03_fm0_neg1_5->Get("h_ST_rwgt25");
  hWWW_run_03_fm0_neg1_5->SetLineColor(kGreen+3);
  hWWW_run_03_fm0_neg1_5->SetLineWidth(2);
  hWWW_run_03_fm0_neg1_5->Rebin(2);
  hWWW_run_03_fm0_neg1_5->SetStats(kFALSE);
  hWWW_run_03_fm0_neg1_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_03_fm0_neg1_5->Scale(LUMI*BR);
  hWWW_run_03_fm0_neg1_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_03_fm0_neg1_5->Draw("hist same");
 
  TFile* fileWWW_run_04_fm0_neg1_0 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_04_fm0_neg1_0 = (TH1F*)fileWWW_run_04_fm0_neg1_0->Get("h_ST_rwgt26");
  hWWW_run_04_fm0_neg1_0->SetLineColor(kCyan);
  hWWW_run_04_fm0_neg1_0->SetLineWidth(2);
  hWWW_run_04_fm0_neg1_0->Rebin(2);
  hWWW_run_04_fm0_neg1_0->SetStats(kFALSE);
  hWWW_run_04_fm0_neg1_0->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_04_fm0_neg1_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_04_fm0_neg1_0->Scale(LUMI*BR);
  hWWW_run_04_fm0_neg1_0->Draw("hist same");
 
  TFile* fileWWW_run_05_fm0_neg0_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_05_fm0_neg0_5 = (TH1F*)fileWWW_run_05_fm0_neg0_5->Get("h_ST_rwgt27");
  hWWW_run_05_fm0_neg0_5->SetLineColor(kMagenta);
  hWWW_run_05_fm0_neg0_5->SetLineWidth(2);
  hWWW_run_05_fm0_neg0_5->Rebin(2);
  hWWW_run_05_fm0_neg0_5->SetStats(kFALSE);
  hWWW_run_05_fm0_neg0_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_05_fm0_neg0_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_05_fm0_neg0_5->Scale(LUMI*BR);
  hWWW_run_05_fm0_neg0_5->Draw("hist same");
 
  TFile* fileWWW_run_11_fm0_0_0 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_11_fm0_0_0 = (TH1F*)fileWWW_run_11_fm0_0_0->Get("h_ST_rwgt28");
  hWWW_run_11_fm0_0_0->SetLineColor(kViolet+1);
  hWWW_run_11_fm0_0_0->SetLineWidth(2);
  hWWW_run_11_fm0_0_0->Rebin(2);
  hWWW_run_11_fm0_0_0->SetLineStyle(2);
  hWWW_run_11_fm0_0_0->SetStats(kFALSE);
  hWWW_run_11_fm0_0_0->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_11_fm0_0_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_11_fm0_0_0->Scale(LUMI*BR);
  hWWW_run_11_fm0_0_0->Draw("hist same");
 
  TFile* fileWWW_run_06_fm0_pos0_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_06_fm0_pos0_5 = (TH1F*)fileWWW_run_06_fm0_pos0_5->Get("h_ST_rwgt29");
  hWWW_run_06_fm0_pos0_5->SetLineColor(kOrange+3);
  hWWW_run_06_fm0_pos0_5->SetLineWidth(2);
  hWWW_run_06_fm0_pos0_5->Rebin(2);
  hWWW_run_06_fm0_pos0_5->SetStats(kFALSE);
  hWWW_run_06_fm0_pos0_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_06_fm0_pos0_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_06_fm0_pos0_5->Scale(LUMI*BR);
  hWWW_run_06_fm0_pos0_5->Draw("hist same");

  TFile* fileWWW_run_07_fm0_pos1_0 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_07_fm0_pos1_0 = (TH1F*)fileWWW_run_07_fm0_pos1_0->Get("h_ST_rwgt30");
  hWWW_run_07_fm0_pos1_0->SetLineColor(kRed);
  hWWW_run_07_fm0_pos1_0->SetLineWidth(2);
  hWWW_run_07_fm0_pos1_0->Rebin(2);
  hWWW_run_07_fm0_pos1_0->SetStats(kFALSE);
  hWWW_run_07_fm0_pos1_0->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_07_fm0_pos1_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_07_fm0_pos1_0->Scale(LUMI*BR);
  hWWW_run_07_fm0_pos1_0->Draw("hist same");

  TFile* fileWWW_run_08_fm0_pos1_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_08_fm0_pos1_5 = (TH1F*)fileWWW_run_08_fm0_pos1_5->Get("h_ST_rwgt31");
  hWWW_run_08_fm0_pos1_5->SetLineColor(kYellow);
  hWWW_run_08_fm0_pos1_5->SetLineWidth(2);
  hWWW_run_08_fm0_pos1_5->Rebin(2);
  hWWW_run_08_fm0_pos1_5->SetStats(kFALSE);
  hWWW_run_08_fm0_pos1_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_08_fm0_pos1_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_08_fm0_pos1_5->Scale(LUMI*BR);
  hWWW_run_08_fm0_pos1_5->Draw("hist same"); 

  TFile* fileWWW_run_09_fm0_pos2_0 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_09_fm0_pos2_0 = (TH1F*)fileWWW_run_09_fm0_pos2_0->Get("h_ST_rwgt32");
  hWWW_run_09_fm0_pos2_0->SetLineColor(kGreen-3);
  hWWW_run_09_fm0_pos2_0->SetLineWidth(2);
  hWWW_run_09_fm0_pos2_0->Rebin(2);
  hWWW_run_09_fm0_pos2_0->SetStats(kFALSE);
  hWWW_run_09_fm0_pos2_0->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_09_fm0_pos2_0->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_09_fm0_pos2_0->Scale(LUMI*BR);
  hWWW_run_09_fm0_pos2_0->Draw("hist same");

  TFile* fileWWW_run_10_fm0_pos2_5 = new TFile("output_www_aqgc_sapta_skim_1_1_fm0.root");
  TH1F *hWWW_run_10_fm0_pos2_5 = (TH1F*)fileWWW_run_10_fm0_pos2_5->Get("h_ST_rwgt33");
  hWWW_run_10_fm0_pos2_5->SetLineColor(kBlue-7);
  hWWW_run_10_fm0_pos2_5->SetLineWidth(2);
  hWWW_run_10_fm0_pos2_5->Rebin(2);
  hWWW_run_10_fm0_pos2_5->SetStats(kFALSE);
  hWWW_run_10_fm0_pos2_5->GetXaxis()->SetRangeUser(50, 5000.0);
  hWWW_run_10_fm0_pos2_5->GetYaxis()->SetTitleOffset(1.3);
  hWWW_run_10_fm0_pos2_5->Scale(LUMI*BR);
  hWWW_run_10_fm0_pos2_5->Draw("hist same");

  TLegend *leg1 = new TLegend(0.42,0.60,0.80,0.85,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hWWW_run_01_fm0_neg2_5, "fm0_neg2_5", "l");
  leg1->AddEntry(hWWW_run_02_fm0_neg2_0, "fm0_neg2_0", "l");
  leg1->AddEntry(hWWW_run_03_fm0_neg1_5, "fm0_neg1_5", "l");
  leg1->AddEntry(hWWW_run_04_fm0_neg1_0, "fm0_neg1_0", "l");
  leg1->AddEntry(hWWW_run_05_fm0_neg0_5, "fm0_neg0_5", "l");
  leg1->AddEntry(hWWW_run_06_fm0_pos0_5, "fm0_pos0_5", "l");
  leg1->AddEntry(hWWW_run_07_fm0_pos1_0, "fm0_pos1_0", "l");
  leg1->AddEntry(hWWW_run_08_fm0_pos1_5, "fm0_pos1_5", "l");
  leg1->AddEntry(hWWW_run_09_fm0_pos2_0, "fm0_pos2_0", "l");
  leg1->AddEntry(hWWW_run_10_fm0_pos2_5, "fm0_pos2_5", "l");
  leg1->AddEntry(hWWW_run_11_fm0_0_0, "fm0_0_0", "l");
  leg1->Draw();

  c1.SaveAs("AC_ST_fm0.png");
  c1.SaveAs("AC_ST_fm0.pdf");

}
