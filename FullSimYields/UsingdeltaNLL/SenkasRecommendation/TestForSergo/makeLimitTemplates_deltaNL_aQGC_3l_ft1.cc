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
#include "TChain.h"
#include "TGraph.h"
#include "TFile.h"
#include "TF1.h"
#include <TROOT.h>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double LUMI = 35.9;
double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

int makeLimitTemplates_deltaNL_aQGC_3l_ft1(std::string outfile)
{

  TH1D *h_signal = new TH1D("h_signal", "h_signal", 1.0, 0.0, 1.0); h_signal->Sumw2();  
  TH1D *h_bkg = new TH1D("h_bkg", "h_bkg", 1.0, 0.0, 1.0); h_bkg->Sumw2();

  //fill histos
  h_signal->SetBinContent(1.0, 1.0);
  h_bkg->SetBinContent(1.0, 0.013060);

  std::cout << "h_signal->GetNbinsX() = " << h_signal->GetNbinsX() << std::endl;
  std::cout << "h_bkg->GetNbinsX() = " << h_bkg->GetNbinsX() << std::endl;

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_signal->Write();
  h_bkg->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
}

int makeLimitTemplates_deltaNL_aQGC_3l_signal_ft1(std::string outfile)
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

  double param[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 0.564;
  double yield2 = 0.363;
  double yield3 = 0.207;
  double yield4 = 0.0964;
  double yield5 = 0.0306;
  double yield6 = 0.00985;
  double yield7 = 0.0342;
  double yield8 = 0.104;
  double yield9 = 0.218;
  double yield10 = 0.377;
  double yield11 = 0.582;

  double signal_beforeSub[11] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11};

  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  
  TGraph *gr = new TGraph(11, param, signal);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(0.9);
  gr->SetMaximum(15);
  gr->SetMinimum(0.0);
  gr->GetXaxis()->SetLimits(-6.0, 6.0);
  gr->SetTitle("");
  gr->SetTitle("Param (ft1) Vs Yield");
  gr->GetXaxis()->SetTitle("Param (ft1)");
  gr->GetYaxis()->SetTitle("Signal yield");
  gr->Draw("AP*");

  TF1 *bin_content_par1_1 = new TF1("bin_content_par1_1","[0]+[1]*x+[2]*x*x", -6.0, 6.0);
  bin_content_par1_1->SetLineColor(kRed);
  bin_content_par1_1->SetLineWidth(2);
  bin_content_par1_1->SetLineStyle(9);
  gr->Fit("bin_content_par1_1", "", "", -6.0, 6.0);

  c1.SaveAs("Parabola_ft1_ST1500.pdf");
  c1.SaveAs("Parabola_ft1_ST1500.png");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  bin_content_par1_1->Write("bin_content_par1_1");
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;

}
