#include <TROOT.h>
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

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double LUMI = 35.9;
double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

int makeLimitTemplates_deltaNL_aQGC_SS(std::string outfile)
{

  TH1F *h_signal = new TH1F("h_signal", "h_signal", 1.0, 0.0, 1.0); h_signal->Sumw2();  
  TH1F *h_bkg = new TH1F("h_bkg", "h_bkg", 1.0, 0.0, 1.0); h_bkg->Sumw2();

  //fill histos
  h_signal->SetBinContent(1.0, 1.0);
  h_bkg->SetBinContent(1.0, 0.2464*10);

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

int makeLimitTemplates_deltaNL_aQGC_SS_signal(std::string outfile)
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","Cleaning Plot", 10, 10, 600, 400);

  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"}; 

  double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5}; 

  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 4.48827*10;
  double yield2 = 4.13778*10;
  double yield3 = 3.80155*10;
  double yield4 = 3.47958*10;
  double yield5 = 3.17188*10;
  double yield6 = 2.87843*10;
  double yield7 = 2.59925*10;
  double yield8 = 2.33433*10;
  double yield9 = 2.08366*10;
  double yield10 = 1.84726*10;
  double yield11 = 1.62512*10;
  double yield12 = 1.41725*10;
  double yield13 = 1.22363*10;
  double yield14 = 1.04427*10;
  double yield15 = 0.879178*10;
  double yield16 = 0.728344*10;
  double yield17 = 0.591771*10;
  double yield18 = 0.469459*10;
  double yield19 = 0.361408*10;
  double yield20 = 0.267619*10;
  double yield21 = 0.18809*10;
  double yield22 = 0.122823*10;
  double yield23 = 0.0718162*10;
  double yield24 = 0.0350709*10;
  double yield25 = 0.0125866*10;
  double yield26 = 0.00436339*10;
  double yield27 = 0.0104013*10;
  double yield28 = 0.0307002*10;
  double yield29 = 0.0652603*10;
  double yield30 = 0.114081*10;
  double yield31 = 0.177164*10;
  double yield32 = 0.254507*10;
  double yield33 = 0.346111*10;
  double yield34 = 0.451977*10;
  double yield35 = 0.572103*10;
  double yield36 = 0.706491*10;
  double yield37 = 0.85514*10;
  double yield38 = 1.01805*10;
  double yield39 = 1.19522*10;
  double yield40 = 1.38665*10;
  double yield41 = 1.59235*10;
  double yield42 = 1.8123*10;
  double yield43 = 2.04651*10;
  double yield44 = 2.29499*10;
  double yield45 = 2.55773*10;
  double yield46 = 2.83473*10;
  double yield47 = 3.12599*10;
  double yield48 = 3.43151*10;
  double yield49 = 3.75129*10;
  double yield50 = 4.08533*10;
  double yield51 = 4.43363*10;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  TGraph *gr = new TGraph(51, param, signal); 
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(0.9);
  gr->SetMaximum(15);
  gr->SetMinimum(0.0);
  gr->GetXaxis()->SetLimits(-3.0, 3.0);
  gr->SetTitle("");
  gr->SetTitle("Param (ft0) Vs Yield");
  gr->GetXaxis()->SetTitle("Param (ft0)");
  gr->GetYaxis()->SetTitle("Signal yield");
  gr->Draw("AP*");

  
  TF1 *bin_content_par1_1 = new TF1("bin_content_par1_1","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
  bin_content_par1_1->SetLineColor(kRed);
  bin_content_par1_1->SetLineWidth(2);
  bin_content_par1_1->SetLineStyle(9);
  gr->Fit("bin_content_par1_1", "", "", -2.5, 2.5);
  
  /* 
  TF1 *aC_test_SS = new TF1("aC_test_SS","[0]+[1]*x+[2]*x*x", -3.0, 3.0);
  aC_test_SS->SetLineColor(kRed);
  aC_test_SS->SetLineWidth(2);
  aC_test_SS->SetLineStyle(9);
  gr->Fit("aC_test_SS", "", "", -2.5, 2.5);
  */
  c1.SaveAs("Parabola_ft0_ST2000.pdf");
  c1.SaveAs("Parabola_ft0_ST2000.png");

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  bin_content_par1_1->Write("bin_content_par1_1");
  //aC_test_SS->Write("aC_test_SS"); 
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
}
