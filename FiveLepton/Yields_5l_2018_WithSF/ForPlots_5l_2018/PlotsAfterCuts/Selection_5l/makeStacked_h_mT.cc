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
#include "TColor.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

void makeStacked_h_mT(){

  gROOT->SetStyle("Plain");
  gStyle->SetErrorX(0);

  TColor *lightblue  = new TColor(2001,91/255.,187/255.,241/255.);
  TColor *blue       = new TColor(2002,60/255.,144/255.,196/255.);
  TColor *orange     = new TColor(2003,230/255.,159/255.,0/255.);
  TColor *brown      = new TColor(2004,180/255.,117/255.,0/255.);
  TColor *yellow     = new TColor(2005,245/255.,236/255.,69/255.);
  TColor *darkyellow = new TColor(2006,215/255.,200/255.,0/255.);
  TColor *blueviolet = new TColor(2007,70/255.,109/255.,171/255.);
  TColor *violet     = new TColor(2008,70/255.,90/255.,134/255.);
  TColor *darkviolet = new TColor(2009,55/255.,65/255.,100/255.);
  TColor *lightgreen = new TColor(2010,120/255.,160/255.,0/255.);
  TColor *green      = new TColor(2011,0/255.,158/255.,115/255.);
  TColor *pink       = new TColor(2012,204/255.,121/255.,167/255.);
  int ci = TColor::GetColor("#466dab");
  
  TCanvas c1("c1","Stacked Histogram",10,10,700,800);
  //TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.25);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.15, 1, 1);
  p_1->SetTopMargin(0.05);
  p_1->SetBottomMargin(0.15);
  p_1->SetFrameFillColor(0);
  //p_2->SetTopMargin(0.0);
  //p_2->SetBottomMargin(0.15);
  //p_2->SetFrameFillColor(0);
  //p_1->SetLogy();
  p_1->Draw();
  //p_2->Draw();
  p_1->cd();

  //double LUMI = 59.7;
  double LUMI = 137;

  THStack *hs = new THStack("hs","");

  TFile* file2 = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *h2 = (TH1F*)file2->Get("h_mT");
  h2->Rebin(5);
  h2->GetXaxis()->SetRangeUser(0.0, 200.0);
  h2->Scale(LUMI);
  h2->SetLineColor(2011);
  h2->SetFillColor(2011);
  h2->Draw();
  hs->Add(h2);

  TFile* file16 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *h16 = (TH1F*)file16->Get("h_mT");
  h16->Rebin(5);
  h16->Scale(LUMI);
  h16->GetXaxis()->SetRangeUser(0.0, 200.0);
  h16->SetLineColor(2003);
  h16->SetFillColor(2003);
  h16->Draw();
  hs->Add(h16);

  TFile* file17 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *h17 = (TH1F*)file17->Get("h_mT");
  h17->Rebin(5);
  h17->Scale(LUMI);
  h17->GetXaxis()->SetRangeUser(0.0, 200.0);
  h17->SetLineColor(2012);
  h17->SetFillColor(2012);
  h17->Draw();
  hs->Add(h17);
 
  TFile* file0 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *h0 = (TH1F*)file0->Get("h_mT");
  h0->Rebin(5);
  h0->GetXaxis()->SetRangeUser(0.0, 200.0);
  h0->Scale(LUMI); 
  h0->SetLineColor(2001);
  h0->SetFillColor(2001);
  h0->Draw();
  hs->Add(h0);

  TFile* file1 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *h1 = (TH1F*)file1->Get("h_mT");
  h1->Rebin(5);
  h1->GetXaxis()->SetRangeUser(0.0, 200.0);
  h1->Scale(LUMI);
  h1->SetLineColor(2001);
  h1->SetFillColor(2001);
  h1->Draw();
  hs->Add(h1);

  TFile* file3 = new TFile("output_twz_ll_madgraph_1.root");
  TH1F *h3 = (TH1F*)file3->Get("h_mT");
  h3->Rebin(5);
  h3->GetXaxis()->SetRangeUser(0.0, 200.0);
  h3->Scale(LUMI);
  h3->SetLineColor(2001);
  h3->SetFillColor(2001);
  h3->Draw();
  hs->Add(h3); 

  TFile* file4 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *h4 = (TH1F*)file4->Get("h_mT");
  h4->Rebin(5);
  h4->Scale(LUMI);
  h4->GetXaxis()->SetRangeUser(0.0, 200.0);
  h4->SetLineColor(2002);
  h4->SetFillColor(2002);
  h4->Draw();
  hs->Add(h4);

  TFile* file5 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *h5 = (TH1F*)file5->Get("h_mT");
  h5->Rebin(5);
  h5->Scale(LUMI);
  h5->GetXaxis()->SetRangeUser(0.0, 200.0);
  h5->SetLineColor(2002);
  h5->SetFillColor(2002);
  h5->Draw();
  hs->Add(h5);

  TFile* file6 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *h6 = (TH1F*)file6->Get("h_mT");
  h6->Rebin(5);
  h6->Scale(LUMI);
  h6->GetXaxis()->SetRangeUser(0.0, 200.0);
  h6->SetLineColor(2002);
  h6->SetFillColor(2002);
  h6->Draw();
  hs->Add(h6);
  
  TFile* file7 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *h7 = (TH1F*)file7->Get("h_mT");
  h7->Rebin(5);
  h7->Scale(LUMI);
  h7->GetXaxis()->SetRangeUser(0.0, 200.0);
  h7->SetLineColor(2002);
  h7->SetFillColor(2002);
  h7->Draw();
  hs->Add(h7);

  TFile* file8 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *h8 = (TH1F*)file8->Get("h_mT");
  h8->Rebin(5);
  h8->Scale(LUMI);
  h8->GetXaxis()->SetRangeUser(0.0, 200.0);
  h8->SetLineColor(2002);
  h8->SetFillColor(2002);
  h8->Draw();
  hs->Add(h8);

  TFile* file9 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *h9 = (TH1F*)file9->Get("h_mT");
  h9->Rebin(5);
  h9->Scale(LUMI);
  h9->GetXaxis()->SetRangeUser(0.0, 200.0);
  h9->SetLineColor(2002);
  h9->SetFillColor(2002);
  h9->Draw();
  hs->Add(h9);
  
  TFile* file10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *h10 = (TH1F*)file10->Get("h_mT");
  h10->Rebin(5);
  h10->Scale(LUMI);
  h10->GetXaxis()->SetRangeUser(0.0, 200.0);
  h10->SetLineColor(2002);
  h10->SetFillColor(2002);
  h10->Draw();
  hs->Add(h10);

  TFile* file11 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *h11 = (TH1F*)file11->Get("h_mT");
  h11->Rebin(5);
  h11->Scale(LUMI);
  h11->GetXaxis()->SetRangeUser(0.0, 200.0);
  h11->SetLineColor(2002);
  h11->SetFillColor(2002);
  h11->Draw();
  hs->Add(h11);

  TFile* file12 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *h12 = (TH1F*)file12->Get("h_mT");
  h12->Rebin(5);
  h12->Scale(LUMI);
  h12->GetXaxis()->SetRangeUser(0.0, 200.0);
  h12->SetLineColor(2002);
  h12->SetFillColor(2002);
  h12->Draw();
  hs->Add(h12);

  TFile* file13 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *h13 = (TH1F*)file13->Get("h_mT");
  h13->Rebin(5);
  h13->Scale(LUMI);
  h13->GetXaxis()->SetRangeUser(0.0, 200.0);
  h13->SetLineColor(2002);
  h13->SetFillColor(2002);
  h13->Draw();
  hs->Add(h13);

  TFile* file14 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *h14 = (TH1F*)file14->Get("h_mT");
  h14->Rebin(5);
  h14->Scale(LUMI);
  h14->GetXaxis()->SetRangeUser(0.0, 200.0);
  h14->SetLineColor(2005);
  h14->SetFillColor(2005);
  h14->Draw();
  hs->Add(h14);

  TFile* file15 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *h15 = (TH1F*)file15->Get("h_mT");
  h15->Rebin(5);
  h15->Scale(LUMI);
  h15->GetXaxis()->SetRangeUser(0.0, 200.0);
  h15->SetLineColor(2005);
  h15->SetFillColor(2005);
  h15->Draw();
  hs->Add(h15);

  hs->Draw("HIST");
  //hs->SetMaximum(500000.0);
  hs->SetMaximum(1.0);
  hs->GetXaxis()->SetRangeUser(0.0, 200.0);
  hs->SetMinimum(0.1);
  hs->GetXaxis()->SetTitle("M_{T} constructed with remaining lepton");
  hs->GetYaxis()->SetLabelSize(0.03);
  hs->GetYaxis()->SetTitleSize(0.03);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetXaxis()->SetTitleOffset(1.3);
  hs->GetYaxis()->SetTitleOffset(1.3);
  hs->GetXaxis()->SetLabelSize(0.03);
  hs->GetXaxis()->SetTitleSize(0.03);

  TFile* file20 = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *h20 = (TH1F*)file20->Get("h_mT");
  h20->Rebin(5);
  h20->Scale(LUMI);
  h20->GetXaxis()->SetRangeUser(0.0, 200.0);
  h20->SetLineColor(kRed);
  h20->SetLineWidth(3);
  h20->Draw("SAME HIST");

  TH1F *h_AllMC=(TH1F*)h20->Clone("h_AllMC");
  h_AllMC->Reset();
  h_AllMC->Add(h0);
  h_AllMC->Add(h1);
  h_AllMC->Add(h2);
  h_AllMC->Add(h3);
  h_AllMC->Add(h4);
  h_AllMC->Add(h5);
  h_AllMC->Add(h6);
  h_AllMC->Add(h7);
  h_AllMC->Add(h8);
  h_AllMC->Add(h9);
  h_AllMC->Add(h10);
  h_AllMC->Add(h11);
  h_AllMC->Add(h12);
  h_AllMC->Add(h13);
  h_AllMC->Add(h14);
  h_AllMC->Add(h15);
  h_AllMC->Add(h16);
  h_AllMC->Add(h17);
  /*h_AllMC->Add(h18);
  h_AllMC->Add(h19);
  h_AllMC->Add(h20);
  h_AllMC->Add(h21);
  h_AllMC->Add(h22);
  h_AllMC->Add(h23);
  h_AllMC->Add(h24);
  h_AllMC->Add(h25);
  h_AllMC->Add(h26);
  h_AllMC->Add(h27);
  h_AllMC->Add(h28);
  */

  TH1F *h_AllMC_Unc=(TH1F*)h_AllMC->Clone("h_AllMC_Unc");
  for (int ibin = 1; ibin < h_AllMC_Unc->GetNbinsX()+1; ibin++){

    double uncStat = 0;
    double uncTot  = 0;
    uncStat += h0->GetBinError(ibin)*h0->GetBinError(ibin) +
               h1->GetBinError(ibin)*h1->GetBinError(ibin) + 
               h2->GetBinError(ibin)*h2->GetBinError(ibin) + 
               h3->GetBinError(ibin)*h3->GetBinError(ibin) +
               h4->GetBinError(ibin)*h4->GetBinError(ibin) +
               h5->GetBinError(ibin)*h5->GetBinError(ibin) +
               h6->GetBinError(ibin)*h6->GetBinError(ibin) +
               h7->GetBinError(ibin)*h7->GetBinError(ibin) +
               h8->GetBinError(ibin)*h8->GetBinError(ibin) + 
               h9->GetBinError(ibin)*h9->GetBinError(ibin) + 
               h10->GetBinError(ibin)*h10->GetBinError(ibin) + 
               h11->GetBinError(ibin)*h11->GetBinError(ibin) + 
               h12->GetBinError(ibin)*h12->GetBinError(ibin) +
               h13->GetBinError(ibin)*h13->GetBinError(ibin) + 
               h14->GetBinError(ibin)*h14->GetBinError(ibin) +
               h15->GetBinError(ibin)*h15->GetBinError(ibin) +
               h16->GetBinError(ibin)*h16->GetBinError(ibin) +
               h17->GetBinError(ibin)*h17->GetBinError(ibin);
               /*h18->GetBinError(ibin)*h18->GetBinError(ibin) +
               h19->GetBinError(ibin)*h19->GetBinError(ibin) +
               h20->GetBinError(ibin)*h20->GetBinError(ibin) +
               h21->GetBinError(ibin)*h21->GetBinError(ibin) +
               h22->GetBinError(ibin)*h22->GetBinError(ibin) +
               h23->GetBinError(ibin)*h23->GetBinError(ibin) +
               h24->GetBinError(ibin)*h24->GetBinError(ibin) +
               h25->GetBinError(ibin)*h25->GetBinError(ibin) +
               h26->GetBinError(ibin)*h26->GetBinError(ibin) +
               h27->GetBinError(ibin)*h27->GetBinError(ibin) +
               h28->GetBinError(ibin)*h28->GetBinError(ibin) ;*/
    uncTot = sqrt(uncStat);
    h_AllMC_Unc->SetBinError(ibin, uncTot);

  }

  gStyle->SetHatchesLineWidth(1);
  gStyle->SetErrorX(0.5);
  h_AllMC_Unc->SetFillStyle(3544);
  h_AllMC_Unc->SetLineColor(1);
  h_AllMC_Unc->SetLineStyle(1);
  h_AllMC_Unc->SetLineWidth(1);
  h_AllMC_Unc->SetFillColor(1);
  h_AllMC_Unc->SetMarkerStyle(1);
  h_AllMC_Unc->Draw("SAME E2");
   
 
  TLatex *tLumi = new TLatex(0.90,0.955,"137 fb^{-1} (13 TeV)");
  tLumi->SetNDC();
  tLumi->SetTextAlign(31);
  tLumi->SetTextFont(42);
  tLumi->SetTextSize(0.042);
  tLumi->SetLineWidth(2);
  TLatex *tECM = new TLatex(0.90,0.955,"(13 TeV)");
  tECM->SetNDC();
  tECM->SetTextAlign(31);
  tECM->SetTextFont(42);
  tECM->SetTextSize(0.042);
  tECM->SetLineWidth(2);
  tLumi->Draw();
  TLatex *tCMS = new TLatex(0.125,0.955,"CMS");
  tCMS->SetNDC();
  tCMS->SetTextAlign(11);
  tCMS->SetTextFont(61);
  tCMS->SetTextSize(0.0525);
  tCMS->SetLineWidth(2);
  tCMS->Draw();
  TLatex *tPrel = new TLatex(0.295,0.955,"Preliminary");
  tPrel->SetNDC();
  tPrel->SetTextAlign(11);
  tPrel->SetTextFont(52);
  tPrel->SetTextSize(0.042);
  tPrel->SetLineWidth(2);
  tPrel->Draw();
 
  TLegend *leg1 = new TLegend(0.2,0.73,0.5,0.9075,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.033);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(1);
  leg1->SetLineWidth(2);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1001);
  leg1->AddEntry(h0, "t#bar{t}Z", "f");
  leg1->AddEntry(h2, "t#bar{t}W", "f");
  leg1->AddEntry(h4, "ZZ + Higgs#rightarrow 4l", "f");
  //leg1->AddEntry(h14,"DY", "f");
  //leg1->AddEntry(h10,"Higgs#rightarrow 4l", "f");
  leg1->Draw();
  TLegend *leg1SS = new TLegend(0.5,0.73,0.8,0.9075,NULL,"brNDC");
  leg1SS->SetBorderSize(0);
  leg1SS->SetTextSize(0.035);
  leg1SS->SetLineColor(1);
  leg1SS->SetLineStyle(1);
  leg1SS->SetLineWidth(2);
  leg1SS->SetFillColor(0);
  leg1SS->SetFillStyle(1001);
  leg1SS->AddEntry(h14,"DY", "f");
  leg1SS->AddEntry(h16,"t#bar{t}", "f");
  leg1SS->AddEntry(h17, "WZ", "f");
  leg1SS->AddEntry(h20, "wzz signal", "l");
  leg1SS->Draw();
  TLegend *leg13l = new TLegend(0.1667,0.67,0.5,0.89,NULL,"brNDC");
  leg13l->SetBorderSize(0);
  leg13l->SetTextSize(0.035);
  leg13l->SetLineColor(1);
  leg13l->SetLineStyle(1);
  leg13l->SetLineWidth(2);
  leg13l->SetFillColor(0);
  leg13l->SetFillStyle(1001);
  /* 
  p_2->cd();
  p_2->SetGridy();

  TH1F *h_ratio=(TH1F*)h_data->Clone("h_ratio");
  h_ratio->SetLabelSize(0.05);
  h_ratio->SetTitleSize(0.10);
  h_ratio->SetTitle("; M_T constructed with remaining lepton");
  h_ratio->GetXaxis()->SetTitleOffset(1.1);
  h_ratio->GetXaxis()->SetTitleSize(0.10);
  h_ratio->GetYaxis()->SetTitleOffset(0.4);
  h_ratio->GetYaxis()->SetTitleSize(0.10);
  h_ratio->GetYaxis()->SetLabelSize(0.07);
  h_ratio->GetYaxis()->SetTitle("Data/MC Ratio  ");
  h_ratio->SetStats(kFALSE);
  h_ratio->Divide(h_AllMC);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->SetLabelSize(0.10);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->SetMinimum(-0.5);
  h_ratio->SetMaximum(2.5);
  h_ratio->Draw();

  TF1 *fit_ratio = new TF1("fit_ratio","[0]*x + [1]", 0.0, 1000.0);
  fit_ratio->SetParLimits(0,0.0,0.00001);
  fit_ratio->SetParLimits(1,0.0, 1.2);
  fit_ratio->SetLineColor(kRed);
  fit_ratio->SetLineWidth(3);
  h_ratio->Fit("fit_ratio", "", "", 0.0, 1000.0);
  
  
  
  TH1F *h_ratio_Unc=(TH1F*)h_ratio->Clone("h_ratio_Unc");
  for (int ibin = 1; ibin < h_ratio_Unc->GetNbinsX()+1; ibin++)
  {
    h_ratio_Unc->SetBinError(ibin, (h_ratio->GetBinContent(ibin)*h_AllMC_Unc->GetBinError(ibin))/h_AllMC->GetBinContent(ibin));
  }

  gStyle->SetHatchesLineWidth(1);
  gStyle->SetErrorX(0.5);
  h_ratio_Unc->SetFillStyle(3544);
  h_ratio_Unc->SetFillColor(1);
  h_ratio_Unc->SetMarkerStyle(1);
  h_ratio_Unc->Draw("SAME E2");
  */

  c1.SaveAs("h_mT.pdf");
  c1.SaveAs("h_mT.png");

}
