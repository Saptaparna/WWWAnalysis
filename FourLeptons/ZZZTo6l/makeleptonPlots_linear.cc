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

void makeleptonPlots_linear()
{

  TCanvas c1("c1","Stacked Histogram",10,12,700,850);
  //c1.SetLogx();
  //c1.SetLogy();

  TFile* fileLeptonPtOne = new TFile("test.root");
  TH1F *hLeptonPtOne = (TH1F*)fileLeptonPtOne->Get("h_Lepton_1_Pt");
  hLeptonPtOne->SetLineColor(kBlack);
  hLeptonPtOne->SetMaximum(10000);
  hLeptonPtOne->SetTitle("");
  hLeptonPtOne->SetLineWidth(2);
  hLeptonPtOne->Rebin(2);
  hLeptonPtOne->GetXaxis()->SetRangeUser(0, 60.0);
  hLeptonPtOne->SetMinimum(5.0);
  hLeptonPtOne->SetStats(kFALSE);
  //hLeptonPtOne->GetYaxis()->SetTitleOffset(1.8);
  hLeptonPtOne->GetXaxis()->SetTitle("Transverse momentum of leptons in a ZZZ sample");
  //hLeptonPtOne->GetYaxis()->SetTitle("Events");
  hLeptonPtOne->Draw("HIST");

  TFile* fileLeptonPtTwo = new TFile("test.root");
  TH1F *hLeptonPtTwo = (TH1F*)fileLeptonPtTwo->Get("h_Lepton_2_Pt");
  hLeptonPtTwo->SetLineColor(kRed);
  hLeptonPtTwo->SetLineWidth(2);
  hLeptonPtTwo->Rebin(2);
  hLeptonPtTwo->GetXaxis()->SetRangeUser(0, 750.0);
  hLeptonPtTwo->SetStats(kFALSE);
  hLeptonPtTwo->GetYaxis()->SetTitleOffset(1.3);
  hLeptonPtTwo->Draw("HIST SAME");

  TFile* fileLeptonPtThree = new TFile("test.root");
  TH1F *LeptonPtThree = (TH1F*)fileLeptonPtThree->Get("h_Lepton_3_Pt");
  LeptonPtThree->SetLineColor(kBlue);
  LeptonPtThree->SetLineWidth(2);
  LeptonPtThree->Rebin(2);
  LeptonPtThree->GetXaxis()->SetRangeUser(0, 750.0);
  LeptonPtThree->SetStats(kFALSE);
  LeptonPtThree->GetYaxis()->SetTitleOffset(1.3);
  LeptonPtThree->Draw("HIST SAME");

  TFile* fileLeptonPtFour = new TFile("test.root");
  TH1F *LeptonPtFour = (TH1F*)fileLeptonPtFour->Get("h_Lepton_4_Pt");
  LeptonPtFour->SetLineColor(kMagenta);
  LeptonPtFour->SetLineWidth(2);
  LeptonPtFour->Rebin(2);
  LeptonPtFour->GetXaxis()->SetRangeUser(0, 750.0);
  LeptonPtFour->SetStats(kFALSE);
  LeptonPtFour->GetYaxis()->SetTitleOffset(1.3);
  LeptonPtFour->Draw("HIST SAME");

  TFile* fileLeptonPtFive = new TFile("test.root");
  TH1F *LeptonPtFive = (TH1F*)fileLeptonPtFive->Get("h_Lepton_5_Pt");
  LeptonPtFive->SetLineColor(kCyan);
  LeptonPtFive->SetLineWidth(2);
  LeptonPtFive->Rebin(2);
  LeptonPtFive->GetXaxis()->SetRangeUser(0, 750.0);
  LeptonPtFive->SetStats(kFALSE);
  LeptonPtFive->GetYaxis()->SetTitleOffset(1.3);
  LeptonPtFive->Draw("HIST SAME");

  TFile* fileLeptonPtSix = new TFile("test.root");
  TH1F *LeptonPtSix = (TH1F*)fileLeptonPtSix->Get("h_Lepton_6_Pt");
  LeptonPtSix->SetLineColor(kOrange);
  LeptonPtSix->SetLineWidth(2);
  LeptonPtSix->Rebin(2);
  LeptonPtSix->GetXaxis()->SetRangeUser(0, 750.0);
  LeptonPtSix->SetStats(kFALSE);
  LeptonPtSix->GetYaxis()->SetTitleOffset(1.3);
  LeptonPtSix->Draw("HIST SAME");

  TLegend *leg1 = new TLegend(0.28,0.55,0.48,0.82,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hLeptonPtOne,"Leading lepton p_{T} Mean: 127 GeV","l");
  leg1->AddEntry(hLeptonPtTwo,"Second lepton p_{T} Mean: 98 GeV ","l");
  leg1->AddEntry(LeptonPtThree,"Third lepton p_{T} Mean: 68 GeV","l");
  leg1->AddEntry(LeptonPtFour,"Fourth lepton p_{T} Mean: 50 GeV","l");
  leg1->AddEntry(LeptonPtFive,"Fifth lepton p_{T} Mean: 35 GeV","l");
  leg1->AddEntry(LeptonPtSix,"Sixth lepton p_{T} Mean: 23 GeV","l");
  //leg1->Draw();

  c1.SaveAs("lepton_ZZZ_linear.png");
  c1.SaveAs("lepton_ZZZ_linear.pdf");
}
