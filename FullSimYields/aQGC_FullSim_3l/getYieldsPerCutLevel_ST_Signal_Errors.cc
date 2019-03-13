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

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double LUMI = 35.9;
double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

double total_noST_rwgt67, total_noST_rwgt68, total_noST_rwgt69, total_noST_rwgt70, total_noST_rwgt71, total_noST_rwgt72, total_noST_rwgt73, total_noST_rwgt74, total_noST_rwgt75, total_noST_rwgt76, total_noST_rwgt77;
double total_ST250_rwgt67, total_ST250_rwgt68, total_ST250_rwgt69, total_ST250_rwgt70, total_ST250_rwgt71, total_ST250_rwgt72, total_ST250_rwgt73, total_ST250_rwgt74, total_ST250_rwgt75, total_ST250_rwgt76, total_ST250_rwgt77;
double total_ST500_rwgt67, total_ST500_rwgt68, total_ST500_rwgt69, total_ST500_rwgt70, total_ST500_rwgt71, total_ST500_rwgt72, total_ST500_rwgt73, total_ST500_rwgt74, total_ST500_rwgt75, total_ST500_rwgt76, total_ST500_rwgt77;
double total_ST750_rwgt67, total_ST750_rwgt68, total_ST750_rwgt69, total_ST750_rwgt70, total_ST750_rwgt71, total_ST750_rwgt72, total_ST750_rwgt73, total_ST750_rwgt74, total_ST750_rwgt75, total_ST750_rwgt76, total_ST750_rwgt77;
double total_ST1000_rwgt67, total_ST1000_rwgt68, total_ST1000_rwgt69, total_ST1000_rwgt70, total_ST1000_rwgt71, total_ST1000_rwgt72, total_ST1000_rwgt73, total_ST1000_rwgt74, total_ST1000_rwgt75, total_ST1000_rwgt76, total_ST1000_rwgt77;
double total_ST1500_rwgt67, total_ST1500_rwgt68, total_ST1500_rwgt69, total_ST1500_rwgt70, total_ST1500_rwgt71, total_ST1500_rwgt72, total_ST1500_rwgt73, total_ST1500_rwgt74, total_ST1500_rwgt75, total_ST1500_rwgt76, total_ST1500_rwgt77;
double total_ST2000_rwgt67, total_ST2000_rwgt68, total_ST2000_rwgt69, total_ST2000_rwgt70, total_ST2000_rwgt71, total_ST2000_rwgt72, total_ST2000_rwgt73, total_ST2000_rwgt74, total_ST2000_rwgt75, total_ST2000_rwgt76, total_ST2000_rwgt77;
double total_ST2500_rwgt67, total_ST2500_rwgt68, total_ST2500_rwgt69, total_ST2500_rwgt70, total_ST2500_rwgt71, total_ST2500_rwgt72, total_ST2500_rwgt73, total_ST2500_rwgt74, total_ST2500_rwgt75, total_ST2500_rwgt76, total_ST2500_rwgt77;

double error_noST_rwgt67_0SFOS, error_noST_rwgt68_0SFOS, error_noST_rwgt69_0SFOS, error_noST_rwgt70_0SFOS, error_noST_rwgt71_0SFOS, error_noST_rwgt72_0SFOS, error_noST_rwgt73_0SFOS, error_noST_rwgt74_0SFOS, error_noST_rwgt75_0SFOS, error_noST_rwgt76_0SFOS, error_noST_rwgt77_0SFOS;
double error_ST250_rwgt67_0SFOS, error_ST250_rwgt68_0SFOS, error_ST250_rwgt69_0SFOS, error_ST250_rwgt70_0SFOS, error_ST250_rwgt71_0SFOS, error_ST250_rwgt72_0SFOS, error_ST250_rwgt73_0SFOS, error_ST250_rwgt74_0SFOS, error_ST250_rwgt75_0SFOS, error_ST250_rwgt76_0SFOS, error_ST250_rwgt77_0SFOS;
double error_ST500_rwgt67_0SFOS, error_ST500_rwgt68_0SFOS, error_ST500_rwgt69_0SFOS, error_ST500_rwgt70_0SFOS, error_ST500_rwgt71_0SFOS, error_ST500_rwgt72_0SFOS, error_ST500_rwgt73_0SFOS, error_ST500_rwgt74_0SFOS, error_ST500_rwgt75_0SFOS, error_ST500_rwgt76_0SFOS, error_ST500_rwgt77_0SFOS;
double error_ST750_rwgt67_0SFOS, error_ST750_rwgt68_0SFOS, error_ST750_rwgt69_0SFOS, error_ST750_rwgt70_0SFOS, error_ST750_rwgt71_0SFOS, error_ST750_rwgt72_0SFOS, error_ST750_rwgt73_0SFOS, error_ST750_rwgt74_0SFOS, error_ST750_rwgt75_0SFOS, error_ST750_rwgt76_0SFOS, error_ST750_rwgt77_0SFOS;
double error_ST1000_rwgt67_0SFOS, error_ST1000_rwgt68_0SFOS, error_ST1000_rwgt69_0SFOS, error_ST1000_rwgt70_0SFOS, error_ST1000_rwgt71_0SFOS, error_ST1000_rwgt72_0SFOS, error_ST1000_rwgt73_0SFOS, error_ST1000_rwgt74_0SFOS, error_ST1000_rwgt75_0SFOS, error_ST1000_rwgt76_0SFOS, error_ST1000_rwgt77_0SFOS;
double error_ST1500_rwgt67_0SFOS, error_ST1500_rwgt68_0SFOS, error_ST1500_rwgt69_0SFOS, error_ST1500_rwgt70_0SFOS, error_ST1500_rwgt71_0SFOS, error_ST1500_rwgt72_0SFOS, error_ST1500_rwgt73_0SFOS, error_ST1500_rwgt74_0SFOS, error_ST1500_rwgt75_0SFOS, error_ST1500_rwgt76_0SFOS, error_ST1500_rwgt77_0SFOS;
double error_ST2000_rwgt67_0SFOS, error_ST2000_rwgt68_0SFOS, error_ST2000_rwgt69_0SFOS, error_ST2000_rwgt70_0SFOS, error_ST2000_rwgt71_0SFOS, error_ST2000_rwgt72_0SFOS, error_ST2000_rwgt73_0SFOS, error_ST2000_rwgt74_0SFOS, error_ST2000_rwgt75_0SFOS, error_ST2000_rwgt76_0SFOS, error_ST2000_rwgt77_0SFOS;
double error_ST2500_rwgt67_0SFOS, error_ST2500_rwgt68_0SFOS, error_ST2500_rwgt69_0SFOS, error_ST2500_rwgt70_0SFOS, error_ST2500_rwgt71_0SFOS, error_ST2500_rwgt72_0SFOS, error_ST2500_rwgt73_0SFOS, error_ST2500_rwgt74_0SFOS, error_ST2500_rwgt75_0SFOS, error_ST2500_rwgt76_0SFOS, error_ST2500_rwgt77_0SFOS;

double error_noST_rwgt67_1SFOS, error_noST_rwgt68_1SFOS, error_noST_rwgt69_1SFOS, error_noST_rwgt70_1SFOS, error_noST_rwgt71_1SFOS, error_noST_rwgt72_1SFOS, error_noST_rwgt73_1SFOS, error_noST_rwgt74_1SFOS, error_noST_rwgt75_1SFOS, error_noST_rwgt76_1SFOS, error_noST_rwgt77_1SFOS;
double error_ST250_rwgt67_1SFOS, error_ST250_rwgt68_1SFOS, error_ST250_rwgt69_1SFOS, error_ST250_rwgt70_1SFOS, error_ST250_rwgt71_1SFOS, error_ST250_rwgt72_1SFOS, error_ST250_rwgt73_1SFOS, error_ST250_rwgt74_1SFOS, error_ST250_rwgt75_1SFOS, error_ST250_rwgt76_1SFOS, error_ST250_rwgt77_1SFOS;
double error_ST500_rwgt67_1SFOS, error_ST500_rwgt68_1SFOS, error_ST500_rwgt69_1SFOS, error_ST500_rwgt70_1SFOS, error_ST500_rwgt71_1SFOS, error_ST500_rwgt72_1SFOS, error_ST500_rwgt73_1SFOS, error_ST500_rwgt74_1SFOS, error_ST500_rwgt75_1SFOS, error_ST500_rwgt76_1SFOS, error_ST500_rwgt77_1SFOS;
double error_ST750_rwgt67_1SFOS, error_ST750_rwgt68_1SFOS, error_ST750_rwgt69_1SFOS, error_ST750_rwgt70_1SFOS, error_ST750_rwgt71_1SFOS, error_ST750_rwgt72_1SFOS, error_ST750_rwgt73_1SFOS, error_ST750_rwgt74_1SFOS, error_ST750_rwgt75_1SFOS, error_ST750_rwgt76_1SFOS, error_ST750_rwgt77_1SFOS;
double error_ST1000_rwgt67_1SFOS, error_ST1000_rwgt68_1SFOS, error_ST1000_rwgt69_1SFOS, error_ST1000_rwgt70_1SFOS, error_ST1000_rwgt71_1SFOS, error_ST1000_rwgt72_1SFOS, error_ST1000_rwgt73_1SFOS, error_ST1000_rwgt74_1SFOS, error_ST1000_rwgt75_1SFOS, error_ST1000_rwgt76_1SFOS, error_ST1000_rwgt77_1SFOS;
double error_ST1500_rwgt67_1SFOS, error_ST1500_rwgt68_1SFOS, error_ST1500_rwgt69_1SFOS, error_ST1500_rwgt70_1SFOS, error_ST1500_rwgt71_1SFOS, error_ST1500_rwgt72_1SFOS, error_ST1500_rwgt73_1SFOS, error_ST1500_rwgt74_1SFOS, error_ST1500_rwgt75_1SFOS, error_ST1500_rwgt76_1SFOS, error_ST1500_rwgt77_1SFOS;
double error_ST2000_rwgt67_1SFOS, error_ST2000_rwgt68_1SFOS, error_ST2000_rwgt69_1SFOS, error_ST2000_rwgt70_1SFOS, error_ST2000_rwgt71_1SFOS, error_ST2000_rwgt72_1SFOS, error_ST2000_rwgt73_1SFOS, error_ST2000_rwgt74_1SFOS, error_ST2000_rwgt75_1SFOS, error_ST2000_rwgt76_1SFOS, error_ST2000_rwgt77_1SFOS;
double error_ST2500_rwgt67_1SFOS, error_ST2500_rwgt68_1SFOS, error_ST2500_rwgt69_1SFOS, error_ST2500_rwgt70_1SFOS, error_ST2500_rwgt71_1SFOS, error_ST2500_rwgt72_1SFOS, error_ST2500_rwgt73_1SFOS, error_ST2500_rwgt74_1SFOS, error_ST2500_rwgt75_1SFOS, error_ST2500_rwgt76_1SFOS, error_ST2500_rwgt77_1SFOS;

double error_noST_rwgt67_2SFOS, error_noST_rwgt68_2SFOS, error_noST_rwgt69_2SFOS, error_noST_rwgt70_2SFOS, error_noST_rwgt71_2SFOS, error_noST_rwgt72_2SFOS, error_noST_rwgt73_2SFOS, error_noST_rwgt74_2SFOS, error_noST_rwgt75_2SFOS, error_noST_rwgt76_2SFOS, error_noST_rwgt77_2SFOS;
double error_ST250_rwgt67_2SFOS, error_ST250_rwgt68_2SFOS, error_ST250_rwgt69_2SFOS, error_ST250_rwgt70_2SFOS, error_ST250_rwgt71_2SFOS, error_ST250_rwgt72_2SFOS, error_ST250_rwgt73_2SFOS, error_ST250_rwgt74_2SFOS, error_ST250_rwgt75_2SFOS, error_ST250_rwgt76_2SFOS, error_ST250_rwgt77_2SFOS;
double error_ST500_rwgt67_2SFOS, error_ST500_rwgt68_2SFOS, error_ST500_rwgt69_2SFOS, error_ST500_rwgt70_2SFOS, error_ST500_rwgt71_2SFOS, error_ST500_rwgt72_2SFOS, error_ST500_rwgt73_2SFOS, error_ST500_rwgt74_2SFOS, error_ST500_rwgt75_2SFOS, error_ST500_rwgt76_2SFOS, error_ST500_rwgt77_2SFOS;
double error_ST750_rwgt67_2SFOS, error_ST750_rwgt68_2SFOS, error_ST750_rwgt69_2SFOS, error_ST750_rwgt70_2SFOS, error_ST750_rwgt71_2SFOS, error_ST750_rwgt72_2SFOS, error_ST750_rwgt73_2SFOS, error_ST750_rwgt74_2SFOS, error_ST750_rwgt75_2SFOS, error_ST750_rwgt76_2SFOS, error_ST750_rwgt77_2SFOS;
double error_ST1000_rwgt67_2SFOS, error_ST1000_rwgt68_2SFOS, error_ST1000_rwgt69_2SFOS, error_ST1000_rwgt70_2SFOS, error_ST1000_rwgt71_2SFOS, error_ST1000_rwgt72_2SFOS, error_ST1000_rwgt73_2SFOS, error_ST1000_rwgt74_2SFOS, error_ST1000_rwgt75_2SFOS, error_ST1000_rwgt76_2SFOS, error_ST1000_rwgt77_2SFOS;
double error_ST1500_rwgt67_2SFOS, error_ST1500_rwgt68_2SFOS, error_ST1500_rwgt69_2SFOS, error_ST1500_rwgt70_2SFOS, error_ST1500_rwgt71_2SFOS, error_ST1500_rwgt72_2SFOS, error_ST1500_rwgt73_2SFOS, error_ST1500_rwgt74_2SFOS, error_ST1500_rwgt75_2SFOS, error_ST1500_rwgt76_2SFOS, error_ST1500_rwgt77_2SFOS;
double error_ST2000_rwgt67_2SFOS, error_ST2000_rwgt68_2SFOS, error_ST2000_rwgt69_2SFOS, error_ST2000_rwgt70_2SFOS, error_ST2000_rwgt71_2SFOS, error_ST2000_rwgt72_2SFOS, error_ST2000_rwgt73_2SFOS, error_ST2000_rwgt74_2SFOS, error_ST2000_rwgt75_2SFOS, error_ST2000_rwgt76_2SFOS, error_ST2000_rwgt77_2SFOS;
double error_ST2500_rwgt67_2SFOS, error_ST2500_rwgt68_2SFOS, error_ST2500_rwgt69_2SFOS, error_ST2500_rwgt70_2SFOS, error_ST2500_rwgt71_2SFOS, error_ST2500_rwgt72_2SFOS, error_ST2500_rwgt73_2SFOS, error_ST2500_rwgt74_2SFOS, error_ST2500_rwgt75_2SFOS, error_ST2500_rwgt76_2SFOS, error_ST2500_rwgt77_2SFOS;


void getYieldsPerCutLevel_0SFOS_noST()
{

  total_noST_rwgt67 = total_noST_rwgt68 = total_noST_rwgt69 = total_noST_rwgt70 = total_noST_rwgt71 = total_noST_rwgt72 = total_noST_rwgt73 = total_noST_rwgt74 = total_noST_rwgt75 = total_noST_rwgt76 = total_noST_rwgt77 = 0.0;
  total_ST250_rwgt67 = total_ST250_rwgt68 = total_ST250_rwgt69 = total_ST250_rwgt70 = total_ST250_rwgt71 = total_ST250_rwgt72 = total_ST250_rwgt73 = total_ST250_rwgt74 = total_ST250_rwgt75 = total_ST250_rwgt76 = total_ST250_rwgt77 = 0.0; 
  total_ST500_rwgt67 = total_ST500_rwgt68 = total_ST500_rwgt69 = total_ST500_rwgt70 = total_ST500_rwgt71 = total_ST500_rwgt72 = total_ST500_rwgt73 = total_ST500_rwgt74 = total_ST500_rwgt75 = total_ST500_rwgt76 = total_ST500_rwgt77 = 0.0;
  total_ST750_rwgt67 = total_ST750_rwgt68 = total_ST750_rwgt69 = total_ST750_rwgt70 = total_ST750_rwgt71 = total_ST750_rwgt72 = total_ST750_rwgt73 = total_ST750_rwgt74 = total_ST750_rwgt75 = total_ST750_rwgt76 = total_ST750_rwgt77 = 0.0; 
  total_ST1000_rwgt67 = total_ST1000_rwgt68 = total_ST1000_rwgt69 = total_ST1000_rwgt70 = total_ST1000_rwgt71 = total_ST1000_rwgt72 = total_ST1000_rwgt73 = total_ST1000_rwgt74 = total_ST1000_rwgt75 = total_ST1000_rwgt76 = total_ST1000_rwgt77 = 0.0;
  total_ST1500_rwgt67 = total_ST1500_rwgt68 = total_ST1500_rwgt69 = total_ST1500_rwgt70 = total_ST1500_rwgt71 = total_ST1500_rwgt72 = total_ST1500_rwgt73 = total_ST1500_rwgt74 = total_ST1500_rwgt75 = total_ST1500_rwgt76 = total_ST1500_rwgt77 = 0.0;
  total_ST2000_rwgt67 = total_ST2000_rwgt68 = total_ST2000_rwgt69 = total_ST2000_rwgt70 = total_ST2000_rwgt71 = total_ST2000_rwgt72 = total_ST2000_rwgt73 = total_ST2000_rwgt74 = total_ST2000_rwgt75 = total_ST2000_rwgt76 = total_ST2000_rwgt77 = 0.0;
  total_ST2500_rwgt67 = total_ST2500_rwgt68 = total_ST2500_rwgt69 = total_ST2500_rwgt70 = total_ST2500_rwgt71 = total_ST2500_rwgt72 = total_ST2500_rwgt73 = total_ST2500_rwgt74 = total_ST2500_rwgt75 = total_ST2500_rwgt76 = total_ST2500_rwgt77 = 0.0;

  double error_noST_rwgt67, error_noST_rwgt68, error_noST_rwgt69, error_noST_rwgt70, error_noST_rwgt71, error_noST_rwgt72, error_noST_rwgt73, error_noST_rwgt74, error_noST_rwgt75, error_noST_rwgt76, error_noST_rwgt77;
double error_ST250_rwgt67, error_ST250_rwgt68, error_ST250_rwgt69, error_ST250_rwgt70, error_ST250_rwgt71, error_ST250_rwgt72, error_ST250_rwgt73, error_ST250_rwgt74, error_ST250_rwgt75, error_ST250_rwgt76, error_ST250_rwgt77;
double error_ST500_rwgt67, error_ST500_rwgt68, error_ST500_rwgt69, error_ST500_rwgt70, error_ST500_rwgt71, error_ST500_rwgt72, error_ST500_rwgt73, error_ST500_rwgt74, error_ST500_rwgt75, error_ST500_rwgt76, error_ST500_rwgt77;
double error_ST750_rwgt67, error_ST750_rwgt68, error_ST750_rwgt69, error_ST750_rwgt70, error_ST750_rwgt71, error_ST750_rwgt72, error_ST750_rwgt73, error_ST750_rwgt74, error_ST750_rwgt75, error_ST750_rwgt76, error_ST750_rwgt77;
double error_ST1000_rwgt67, error_ST1000_rwgt68, error_ST1000_rwgt69, error_ST1000_rwgt70, error_ST1000_rwgt71, error_ST1000_rwgt72, error_ST1000_rwgt73, error_ST1000_rwgt74, error_ST1000_rwgt75, error_ST1000_rwgt76, error_ST1000_rwgt77;
double error_ST1500_rwgt67, error_ST1500_rwgt68, error_ST1500_rwgt69, error_ST1500_rwgt70, error_ST1500_rwgt71, error_ST1500_rwgt72, error_ST1500_rwgt73, error_ST1500_rwgt74, error_ST1500_rwgt75, error_ST1500_rwgt76, error_ST1500_rwgt77;
double error_ST2000_rwgt67, error_ST2000_rwgt68, error_ST2000_rwgt69, error_ST2000_rwgt70, error_ST2000_rwgt71, error_ST2000_rwgt72, error_ST2000_rwgt73, error_ST2000_rwgt74, error_ST2000_rwgt75, error_ST2000_rwgt76, error_ST2000_rwgt77;
double error_ST2500_rwgt67, error_ST2500_rwgt68, error_ST2500_rwgt69, error_ST2500_rwgt70, error_ST2500_rwgt71, error_ST2500_rwgt72, error_ST2500_rwgt73, error_ST2500_rwgt74, error_ST2500_rwgt75, error_ST2500_rwgt76, error_ST2500_rwgt77;

  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");    
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");
  
  std::cout << "0SFOS channel" << std::endl;
  
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl; 
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts except ST}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_noST_rwgt67 = h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt68 = h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt69 = h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt70 = h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt71 = h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt72 = h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt73 = h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt74 = h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt75 = h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt76 = h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt77 = h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR;

  error_noST_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(2)*LUMI*BR;
  error_noST_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(2)*LUMI*BR;
  error_noST_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(2)*LUMI*BR;
  error_noST_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(2)*LUMI*BR;
  error_noST_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(2)*LUMI*BR;
  error_noST_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(2)*LUMI*BR;
  error_noST_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(2)*LUMI*BR;
  error_noST_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(2)*LUMI*BR;
  error_noST_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(2)*LUMI*BR;
  error_noST_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(2)*LUMI*BR;
  error_noST_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(2)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST250()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;    
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST250_rwgt67 = h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt68 = h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt69 = h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt70 = h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt71 = h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt72 = h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt73 = h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt74 = h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt75 = h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt76 = h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt77 = h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR;

  error_ST250_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(3)*LUMI*BR;
  error_ST250_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(3)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;  
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST500_rwgt67 = h_aQGC_rwgt67->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt68 = h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt69 = h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt70 = h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt71 = h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt72 = h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt73 = h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt74 = h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt75 = h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt76 = h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt76 = h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR;

  error_ST500_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;
  error_ST500_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST750()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST750_rwgt67 = h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt68 = h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt69 = h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt70 = h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt71 = h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt72 = h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt73 = h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt74 = h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt75 = h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt76 = h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt77 = h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR;

  error_ST750_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(5)*LUMI*BR;
  error_ST750_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(5)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST1000()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

    
  total_ST1000_rwgt67 = h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt68 = h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt69 = h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt70 = h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt71 = h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt72 = h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt73 = h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt74 = h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt75 = h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt76 = h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt77 = h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR;

  error_ST1000_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(6)*LUMI*BR;
  error_ST1000_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(6)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST1500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST1500_rwgt67 = h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt68 = h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt69 = h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt70 = h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt71 = h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt72 = h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt73 = h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt74 = h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt75 = h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt76 = h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt77 = h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR;

  error_ST1500_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(7)*LUMI*BR;
  error_ST1500_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(7)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST2000()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST2000_rwgt67 = h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt68 = h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt69 = h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt70 = h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt71 = h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt72 = h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt73 = h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt74 = h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt75 = h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt76 = h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt77 = h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR;

  error_ST2000_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(8)*LUMI*BR;
  error_ST2000_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(8)*LUMI*BR;

}

void getYieldsPerCutLevel_0SFOS_ST2500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
  TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt67");
  TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt68");
  TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt69");
  TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt70");
  TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt71");
  TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt72");
  TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt73");
  TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt74");
  TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt75");
  TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt76");
  TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt77");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST2500_rwgt67 = h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt68 = h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt69 = h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt70 = h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt71 = h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt72 = h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt73 = h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt74 = h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt75 = h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt76 = h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt77 = h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR;

  error_ST2500_rwgt67_0SFOS = h_aQGC_rwgt67->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt68_0SFOS = h_aQGC_rwgt68->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt69_0SFOS = h_aQGC_rwgt69->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt70_0SFOS = h_aQGC_rwgt70->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt71_0SFOS = h_aQGC_rwgt71->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt72_0SFOS = h_aQGC_rwgt72->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt73_0SFOS = h_aQGC_rwgt73->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt74_0SFOS = h_aQGC_rwgt74->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt75_0SFOS = h_aQGC_rwgt75->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt76_0SFOS = h_aQGC_rwgt76->GetBinError(9)*LUMI*BR;
  error_ST2500_rwgt77_0SFOS = h_aQGC_rwgt77->GetBinError(9)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_noST()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST_rwgt67 += h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt68 += h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt69 += h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt70 += h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt71 += h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt72 += h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt73 += h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt74 += h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt75 += h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt76 += h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt77 += h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR;

    error_noST_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(2)*LUMI*BR;
    error_noST_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(2)*LUMI*BR;
    error_noST_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(2)*LUMI*BR;
    error_noST_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(2)*LUMI*BR;
    error_noST_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(2)*LUMI*BR;
    error_noST_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(2)*LUMI*BR;
    error_noST_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(2)*LUMI*BR;
    error_noST_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(2)*LUMI*BR;
    error_noST_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(2)*LUMI*BR;
    error_noST_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(2)*LUMI*BR;
    error_noST_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(2)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST250()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl; 
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

    total_ST250_rwgt67 += h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt68 += h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt69 += h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt70 += h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt71 += h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt72 += h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt73 += h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt74 += h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt75 += h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt76 += h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt77 += h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR;

    error_ST250_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(3)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST500_rwgt67 += h_aQGC_rwgt67->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt68 += h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt69 += h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt70 += h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt71 += h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt72 += h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt73 += h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt74 += h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt75 += h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt76 += h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt77 += h_aQGC_rwgt77->GetBinContent(4)*LUMI*BR;

    error_ST500_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST750()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST750_rwgt67 += h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt68 += h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt69 += h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt70 += h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt71 += h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt72 += h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt73 += h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt74 += h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt75 += h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt76 += h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt77 += h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR;

    error_ST750_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(5)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST1000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST1000_rwgt67 += h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt68 += h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt69 += h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt70 += h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt71 += h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt72 += h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt73 += h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt74 += h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt75 += h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt76 += h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt77 += h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR;

    error_ST1000_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(6)*LUMI*BR;


}

void getYieldsPerCutLevel_1SFOS_ST1500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1500_rwgt67 += h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt68 += h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt69 += h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt70 += h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt71 += h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt72 += h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt73 += h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt74 += h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt75 += h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt76 += h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt77 += h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR;

    
    error_ST1500_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(7)*LUMI*BR;


}

void getYieldsPerCutLevel_1SFOS_ST2000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2000_rwgt67 += h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt68 += h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt69 += h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt70 += h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt71 += h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt72 += h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt73 += h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt74 += h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt75 += h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt76 += h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt77 += h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR;

    error_ST2000_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(8)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST2500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt77");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2500_rwgt67 += h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt68 += h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt69 += h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt70 += h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt71 += h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt72 += h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt73 += h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt74 += h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt75 += h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt76 += h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt77 += h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR;

    error_ST2500_rwgt67_1SFOS = h_aQGC_rwgt67->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt68_1SFOS = h_aQGC_rwgt68->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt69_1SFOS = h_aQGC_rwgt69->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt70_1SFOS = h_aQGC_rwgt70->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt71_1SFOS = h_aQGC_rwgt71->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt72_1SFOS = h_aQGC_rwgt72->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt73_1SFOS = h_aQGC_rwgt73->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt74_1SFOS = h_aQGC_rwgt74->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt75_1SFOS = h_aQGC_rwgt75->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt76_1SFOS = h_aQGC_rwgt76->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt77_1SFOS = h_aQGC_rwgt77->GetBinError(9)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_noST()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(2)*LUMI*BR << " \\\\ " << std::endl; 
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST_rwgt67 += h_aQGC_rwgt67->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt68 += h_aQGC_rwgt68->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt69 += h_aQGC_rwgt69->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt70 += h_aQGC_rwgt70->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt71 += h_aQGC_rwgt71->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt72 += h_aQGC_rwgt72->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt73 += h_aQGC_rwgt73->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt74 += h_aQGC_rwgt74->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt75 += h_aQGC_rwgt75->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt76 += h_aQGC_rwgt76->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt77 += h_aQGC_rwgt77->GetBinContent(2)*LUMI*BR;

    error_noST_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(2)*LUMI*BR;
    error_noST_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(2)*LUMI*BR;
    error_noST_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(2)*LUMI*BR;
    error_noST_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(2)*LUMI*BR;
    error_noST_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(2)*LUMI*BR;
    error_noST_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(2)*LUMI*BR;
    error_noST_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(2)*LUMI*BR;
    error_noST_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(2)*LUMI*BR;
    error_noST_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(2)*LUMI*BR;
    error_noST_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(2)*LUMI*BR;
    error_noST_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(2)*LUMI*BR;


}

void getYieldsPerCutLevel_2SFOS_ST250()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(3)*LUMI*BR << " \\\\ " << std::endl; 
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST250_rwgt67 += h_aQGC_rwgt67->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt68 += h_aQGC_rwgt68->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt69 += h_aQGC_rwgt69->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt70 += h_aQGC_rwgt70->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt71 += h_aQGC_rwgt71->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt72 += h_aQGC_rwgt72->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt73 += h_aQGC_rwgt73->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt74 += h_aQGC_rwgt74->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt75 += h_aQGC_rwgt75->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt76 += h_aQGC_rwgt76->GetBinContent(3)*LUMI*BR; 
    total_ST250_rwgt77 += h_aQGC_rwgt77->GetBinContent(3)*LUMI*BR;

    error_ST250_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(3)*LUMI*BR;
    error_ST250_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(3)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_ST500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(4)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(4)*LUMI*BR << " \\\\ " << std::endl; 
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST500_rwgt67 += h_aQGC_rwgt67->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt68 += h_aQGC_rwgt68->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt69 += h_aQGC_rwgt69->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt70 += h_aQGC_rwgt70->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt71 += h_aQGC_rwgt71->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt72 += h_aQGC_rwgt72->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt73 += h_aQGC_rwgt73->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt74 += h_aQGC_rwgt74->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt75 += h_aQGC_rwgt75->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt76 += h_aQGC_rwgt76->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt77 += h_aQGC_rwgt77->GetBinContent(4)*LUMI*BR;

    error_ST500_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;
    error_ST500_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(4)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_ST750()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST750_rwgt67 += h_aQGC_rwgt67->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt68 += h_aQGC_rwgt68->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt69 += h_aQGC_rwgt69->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt70 += h_aQGC_rwgt70->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt71 += h_aQGC_rwgt71->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt72 += h_aQGC_rwgt72->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt73 += h_aQGC_rwgt73->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt74 += h_aQGC_rwgt74->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt75 += h_aQGC_rwgt75->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt76 += h_aQGC_rwgt76->GetBinContent(5)*LUMI*BR; 
    total_ST750_rwgt76 += h_aQGC_rwgt77->GetBinContent(5)*LUMI*BR;

    error_ST750_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(5)*LUMI*BR;
    error_ST750_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(5)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_ST1000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1000_rwgt67 += h_aQGC_rwgt67->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt68 += h_aQGC_rwgt68->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt69 += h_aQGC_rwgt69->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt70 += h_aQGC_rwgt70->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt71 += h_aQGC_rwgt71->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt72 += h_aQGC_rwgt72->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt73 += h_aQGC_rwgt73->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt74 += h_aQGC_rwgt74->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt75 += h_aQGC_rwgt75->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt76 += h_aQGC_rwgt76->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt77 += h_aQGC_rwgt77->GetBinContent(6)*LUMI*BR;

    error_ST1000_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(6)*LUMI*BR;
    error_ST1000_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(6)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_ST1500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST1500_rwgt67 += h_aQGC_rwgt67->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt68 += h_aQGC_rwgt68->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt69 += h_aQGC_rwgt69->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt70 += h_aQGC_rwgt70->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt71 += h_aQGC_rwgt71->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt72 += h_aQGC_rwgt72->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt73 += h_aQGC_rwgt73->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt74 += h_aQGC_rwgt74->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt75 += h_aQGC_rwgt75->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt76 += h_aQGC_rwgt76->GetBinContent(7)*LUMI*BR; 
    total_ST1500_rwgt77 += h_aQGC_rwgt77->GetBinContent(7)*LUMI*BR;

    error_ST1500_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(7)*LUMI*BR;
    error_ST1500_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(7)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST2000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2000_rwgt67 += h_aQGC_rwgt67->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt68 += h_aQGC_rwgt68->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt69 += h_aQGC_rwgt69->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt70 += h_aQGC_rwgt70->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt71 += h_aQGC_rwgt71->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt72 += h_aQGC_rwgt72->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt73 += h_aQGC_rwgt73->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt74 += h_aQGC_rwgt74->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt75 += h_aQGC_rwgt75->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt76 += h_aQGC_rwgt76->GetBinContent(8)*LUMI*BR; 
    total_ST2000_rwgt77 += h_aQGC_rwgt77->GetBinContent(8)*LUMI*BR;
    
    error_ST2000_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(8)*LUMI*BR;
    error_ST2000_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(8)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_ST2500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l.root");
    TH1F *h_aQGC_rwgt67 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt67");
    TH1F *h_aQGC_rwgt68 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt68");
    TH1F *h_aQGC_rwgt69 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt69");
    TH1F *h_aQGC_rwgt70 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt70");
    TH1F *h_aQGC_rwgt71 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt71");
    TH1F *h_aQGC_rwgt72 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt72");
    TH1F *h_aQGC_rwgt73 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt73");
    TH1F *h_aQGC_rwgt74 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt74");
    TH1F *h_aQGC_rwgt75 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt75");
    TH1F *h_aQGC_rwgt76 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt76");
    TH1F *h_aQGC_rwgt77 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt77");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt67->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt68->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt69->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt70->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt71->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt72->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt73->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt74->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt75->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt76->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t0}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt77->GetBinError(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2500_rwgt67 += h_aQGC_rwgt67->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt68 += h_aQGC_rwgt68->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt69 += h_aQGC_rwgt69->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt70 += h_aQGC_rwgt70->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt71 += h_aQGC_rwgt71->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt72 += h_aQGC_rwgt72->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt73 += h_aQGC_rwgt73->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt74 += h_aQGC_rwgt74->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt75 += h_aQGC_rwgt75->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt76 += h_aQGC_rwgt76->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt77 += h_aQGC_rwgt77->GetBinContent(9)*LUMI*BR; 

    error_ST2500_rwgt67_2SFOS = h_aQGC_rwgt67->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt68_2SFOS = h_aQGC_rwgt68->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt69_2SFOS = h_aQGC_rwgt69->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt70_2SFOS = h_aQGC_rwgt70->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt71_2SFOS = h_aQGC_rwgt71->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt72_2SFOS = h_aQGC_rwgt72->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt73_2SFOS = h_aQGC_rwgt73->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt74_2SFOS = h_aQGC_rwgt74->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt75_2SFOS = h_aQGC_rwgt75->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt76_2SFOS = h_aQGC_rwgt76->GetBinError(9)*LUMI*BR;
    error_ST2500_rwgt77_2SFOS = h_aQGC_rwgt77->GetBinError(9)*LUMI*BR;

}

void total()
{

  std::cout << "total_noST_rwgt67 = " << total_noST_rwgt67 << " $\\pm$ " << sqrt(pow(error_noST_rwgt67_0SFOS,2) + pow(error_noST_rwgt67_1SFOS,2) + pow(error_noST_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt68 = " << total_noST_rwgt68 << " $\\pm$ " << sqrt(pow(error_noST_rwgt68_0SFOS,2) + pow(error_noST_rwgt68_1SFOS,2) + pow(error_noST_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt69 = " << total_noST_rwgt69 << " $\\pm$ " << sqrt(pow(error_noST_rwgt69_0SFOS,2) + pow(error_noST_rwgt69_1SFOS,2) + pow(error_noST_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt70 = " << total_noST_rwgt70 << " $\\pm$ " << sqrt(pow(error_noST_rwgt70_0SFOS,2) + pow(error_noST_rwgt70_1SFOS,2) + pow(error_noST_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt71 = " << total_noST_rwgt71 << " $\\pm$ " << sqrt(pow(error_noST_rwgt71_0SFOS,2) + pow(error_noST_rwgt71_1SFOS,2) + pow(error_noST_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt72 = " << total_noST_rwgt72 << " $\\pm$ " << sqrt(pow(error_noST_rwgt72_0SFOS,2) + pow(error_noST_rwgt72_1SFOS,2) + pow(error_noST_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt73 = " << total_noST_rwgt73 << " $\\pm$ " << sqrt(pow(error_noST_rwgt73_0SFOS,2) + pow(error_noST_rwgt73_1SFOS,2) + pow(error_noST_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt74 = " << total_noST_rwgt74 << " $\\pm$ " << sqrt(pow(error_noST_rwgt74_0SFOS,2) + pow(error_noST_rwgt74_1SFOS,2) + pow(error_noST_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt75 = " << total_noST_rwgt75 << " $\\pm$ " << sqrt(pow(error_noST_rwgt75_0SFOS,2) + pow(error_noST_rwgt75_1SFOS,2) + pow(error_noST_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt76 = " << total_noST_rwgt76 << " $\\pm$ " << sqrt(pow(error_noST_rwgt76_0SFOS,2) + pow(error_noST_rwgt76_1SFOS,2) + pow(error_noST_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_noST_rwgt77 = " << total_noST_rwgt77 << " $\\pm$ " << sqrt(pow(error_noST_rwgt77_0SFOS,2) + pow(error_noST_rwgt77_1SFOS,2) + pow(error_noST_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST250_rwgt67 = " << total_ST250_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt67_0SFOS,2) + pow(error_ST250_rwgt67_1SFOS,2) + pow(error_ST250_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt68 = " << total_ST250_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt68_0SFOS,2) + pow(error_ST250_rwgt68_1SFOS,2) + pow(error_ST250_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt69 = " << total_ST250_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt69_0SFOS,2) + pow(error_ST250_rwgt69_1SFOS,2) + pow(error_ST250_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt70 = " << total_ST250_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt70_0SFOS,2) + pow(error_ST250_rwgt70_1SFOS,2) + pow(error_ST250_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt71 = " << total_ST250_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt71_0SFOS,2) + pow(error_ST250_rwgt71_1SFOS,2) + pow(error_ST250_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt72 = " << total_ST250_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt72_0SFOS,2) + pow(error_ST250_rwgt72_1SFOS,2) + pow(error_ST250_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt73 = " << total_ST250_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt73_0SFOS,2) + pow(error_ST250_rwgt73_1SFOS,2) + pow(error_ST250_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt74 = " << total_ST250_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt74_0SFOS,2) + pow(error_ST250_rwgt74_1SFOS,2) + pow(error_ST250_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt75 = " << total_ST250_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt75_0SFOS,2) + pow(error_ST250_rwgt75_1SFOS,2) + pow(error_ST250_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt76 = " << total_ST250_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt76_0SFOS,2) + pow(error_ST250_rwgt76_1SFOS,2) + pow(error_ST250_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST250_rwgt77 = " << total_ST250_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST250_rwgt77_0SFOS,2) + pow(error_ST250_rwgt77_1SFOS,2) + pow(error_ST250_rwgt77_2SFOS,2)) << std::endl;


  std::cout << "total_ST500_rwgt67 = " << total_ST500_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt67_0SFOS,2) + pow(error_ST500_rwgt67_1SFOS,2) + pow(error_ST500_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt68 = " << total_ST500_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt68_0SFOS,2) + pow(error_ST500_rwgt68_1SFOS,2) + pow(error_ST500_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt69 = " << total_ST500_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt69_0SFOS,2) + pow(error_ST500_rwgt69_1SFOS,2) + pow(error_ST500_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt70 = " << total_ST500_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt70_0SFOS,2) + pow(error_ST500_rwgt70_1SFOS,2) + pow(error_ST500_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt71 = " << total_ST500_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt71_0SFOS,2) + pow(error_ST500_rwgt71_1SFOS,2) + pow(error_ST500_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt72 = " << total_ST500_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt72_0SFOS,2) + pow(error_ST500_rwgt72_1SFOS,2) + pow(error_ST500_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt73 = " << total_ST500_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt73_0SFOS,2) + pow(error_ST500_rwgt73_1SFOS,2) + pow(error_ST500_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt74 = " << total_ST500_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt74_0SFOS,2) + pow(error_ST500_rwgt74_1SFOS,2) + pow(error_ST500_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt75 = " << total_ST500_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt75_0SFOS,2) + pow(error_ST500_rwgt75_1SFOS,2) + pow(error_ST500_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt76 = " << total_ST500_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt76_0SFOS,2) + pow(error_ST500_rwgt76_1SFOS,2) + pow(error_ST500_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST500_rwgt77 = " << total_ST500_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST500_rwgt77_0SFOS,2) + pow(error_ST500_rwgt77_1SFOS,2) + pow(error_ST500_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST750_rwgt67 = " << total_ST750_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt67_0SFOS,2) + pow(error_ST750_rwgt67_1SFOS,2) + pow(error_ST750_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt68 = " << total_ST750_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt68_0SFOS,2) + pow(error_ST750_rwgt68_1SFOS,2) + pow(error_ST750_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt69 = " << total_ST750_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt69_0SFOS,2) + pow(error_ST750_rwgt69_1SFOS,2) + pow(error_ST750_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt70 = " << total_ST750_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt70_0SFOS,2) + pow(error_ST750_rwgt70_1SFOS,2) + pow(error_ST750_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt71 = " << total_ST750_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt71_0SFOS,2) + pow(error_ST750_rwgt71_1SFOS,2) + pow(error_ST750_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt72 = " << total_ST750_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt72_0SFOS,2) + pow(error_ST750_rwgt72_1SFOS,2) + pow(error_ST750_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt73 = " << total_ST750_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt73_0SFOS,2) + pow(error_ST750_rwgt73_1SFOS,2) + pow(error_ST750_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt74 = " << total_ST750_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt74_0SFOS,2) + pow(error_ST750_rwgt74_1SFOS,2) + pow(error_ST750_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt75 = " << total_ST750_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt75_0SFOS,2) + pow(error_ST750_rwgt75_1SFOS,2) + pow(error_ST750_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt76 = " << total_ST750_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt76_0SFOS,2) + pow(error_ST750_rwgt76_1SFOS,2) + pow(error_ST750_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST750_rwgt77 = " << total_ST750_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST750_rwgt77_0SFOS,2) + pow(error_ST750_rwgt77_1SFOS,2) + pow(error_ST750_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST1000_rwgt67 = " << total_ST1000_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt67_0SFOS,2) + pow(error_ST1000_rwgt67_1SFOS,2) + pow(error_ST1000_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt68 = " << total_ST1000_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt68_0SFOS,2) + pow(error_ST1000_rwgt68_1SFOS,2) + pow(error_ST1000_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt69 = " << total_ST1000_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt69_0SFOS,2) + pow(error_ST1000_rwgt69_1SFOS,2) + pow(error_ST1000_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt70 = " << total_ST1000_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt70_0SFOS,2) + pow(error_ST1000_rwgt70_1SFOS,2) + pow(error_ST1000_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt71 = " << total_ST1000_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt71_0SFOS,2) + pow(error_ST1000_rwgt71_1SFOS,2) + pow(error_ST1000_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt72 = " << total_ST1000_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt72_0SFOS,2) + pow(error_ST1000_rwgt72_1SFOS,2) + pow(error_ST1000_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt73 = " << total_ST1000_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt73_0SFOS,2) + pow(error_ST1000_rwgt73_1SFOS,2) + pow(error_ST1000_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt74 = " << total_ST1000_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt74_0SFOS,2) + pow(error_ST1000_rwgt74_1SFOS,2) + pow(error_ST1000_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt75 = " << total_ST1000_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt75_0SFOS,2) + pow(error_ST1000_rwgt75_1SFOS,2) + pow(error_ST1000_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt76 = " << total_ST1000_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt76_0SFOS,2) + pow(error_ST1000_rwgt76_1SFOS,2) + pow(error_ST1000_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST1000_rwgt77 = " << total_ST1000_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST1000_rwgt77_0SFOS,2) + pow(error_ST1000_rwgt77_1SFOS,2) + pow(error_ST1000_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST1500_rwgt67 = " << total_ST1500_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt67_0SFOS,2) + pow(error_ST1500_rwgt67_1SFOS,2) + pow(error_ST1500_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt68 = " << total_ST1500_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt68_0SFOS,2) + pow(error_ST1500_rwgt68_1SFOS,2) + pow(error_ST1500_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt69 = " << total_ST1500_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt69_0SFOS,2) + pow(error_ST1500_rwgt69_1SFOS,2) + pow(error_ST1500_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt70 = " << total_ST1500_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt70_0SFOS,2) + pow(error_ST1500_rwgt70_1SFOS,2) + pow(error_ST1500_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt71 = " << total_ST1500_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt71_0SFOS,2) + pow(error_ST1500_rwgt71_1SFOS,2) + pow(error_ST1500_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt72 = " << total_ST1500_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt72_0SFOS,2) + pow(error_ST1500_rwgt72_1SFOS,2) + pow(error_ST1500_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt73 = " << total_ST1500_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt73_0SFOS,2) + pow(error_ST1500_rwgt73_1SFOS,2) + pow(error_ST1500_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt74 = " << total_ST1500_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt74_0SFOS,2) + pow(error_ST1500_rwgt74_1SFOS,2) + pow(error_ST1500_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt75 = " << total_ST1500_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt75_0SFOS,2) + pow(error_ST1500_rwgt75_1SFOS,2) + pow(error_ST1500_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt76 = " << total_ST1500_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt76_0SFOS,2) + pow(error_ST1500_rwgt76_1SFOS,2) + pow(error_ST1500_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST1500_rwgt77 = " << total_ST1500_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST1500_rwgt77_0SFOS,2) + pow(error_ST1500_rwgt77_1SFOS,2) + pow(error_ST1500_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST2000_rwgt67 = " << total_ST2000_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt67_0SFOS,2) + pow(error_ST2000_rwgt67_1SFOS,2) + pow(error_ST2000_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt68 = " << total_ST2000_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt68_0SFOS,2) + pow(error_ST2000_rwgt68_1SFOS,2) + pow(error_ST2000_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt69 = " << total_ST2000_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt69_0SFOS,2) + pow(error_ST2000_rwgt69_1SFOS,2) + pow(error_ST2000_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt70 = " << total_ST2000_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt70_0SFOS,2) + pow(error_ST2000_rwgt70_1SFOS,2) + pow(error_ST2000_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt71 = " << total_ST2000_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt71_0SFOS,2) + pow(error_ST2000_rwgt71_1SFOS,2) + pow(error_ST2000_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt72 = " << total_ST2000_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt72_0SFOS,2) + pow(error_ST2000_rwgt72_1SFOS,2) + pow(error_ST2000_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt73 = " << total_ST2000_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt73_0SFOS,2) + pow(error_ST2000_rwgt73_1SFOS,2) + pow(error_ST2000_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt74 = " << total_ST2000_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt74_0SFOS,2) + pow(error_ST2000_rwgt74_1SFOS,2) + pow(error_ST2000_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt75 = " << total_ST2000_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt75_0SFOS,2) + pow(error_ST2000_rwgt75_1SFOS,2) + pow(error_ST2000_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt76 = " << total_ST2000_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt76_0SFOS,2) + pow(error_ST2000_rwgt76_1SFOS,2) + pow(error_ST2000_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST2000_rwgt77 = " << total_ST2000_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST2000_rwgt77_0SFOS,2) + pow(error_ST2000_rwgt77_1SFOS,2) + pow(error_ST2000_rwgt77_2SFOS,2)) << std::endl;

  std::cout << "total_ST2500_rwgt67 = " << total_ST2500_rwgt67 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt67_0SFOS,2) + pow(error_ST2500_rwgt67_1SFOS,2) + pow(error_ST2500_rwgt67_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt68 = " << total_ST2500_rwgt68 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt68_0SFOS,2) + pow(error_ST2500_rwgt68_1SFOS,2) + pow(error_ST2500_rwgt68_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt69 = " << total_ST2500_rwgt69 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt69_0SFOS,2) + pow(error_ST2500_rwgt69_1SFOS,2) + pow(error_ST2500_rwgt69_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt70 = " << total_ST2500_rwgt70 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt70_0SFOS,2) + pow(error_ST2500_rwgt70_1SFOS,2) + pow(error_ST2500_rwgt70_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt71 = " << total_ST2500_rwgt71 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt71_0SFOS,2) + pow(error_ST2500_rwgt71_1SFOS,2) + pow(error_ST2500_rwgt71_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt72 = " << total_ST2500_rwgt72 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt72_0SFOS,2) + pow(error_ST2500_rwgt72_1SFOS,2) + pow(error_ST2500_rwgt72_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt73 = " << total_ST2500_rwgt73 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt73_0SFOS,2) + pow(error_ST2500_rwgt73_1SFOS,2) + pow(error_ST2500_rwgt73_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt74 = " << total_ST2500_rwgt74 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt74_0SFOS,2) + pow(error_ST2500_rwgt74_1SFOS,2) + pow(error_ST2500_rwgt74_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt75 = " << total_ST2500_rwgt75 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt75_0SFOS,2) + pow(error_ST2500_rwgt75_1SFOS,2) + pow(error_ST2500_rwgt75_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt76 = " << total_ST2500_rwgt76 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt76_0SFOS,2) + pow(error_ST2500_rwgt76_1SFOS,2) + pow(error_ST2500_rwgt76_2SFOS,2)) << std::endl;
  std::cout << "total_ST2500_rwgt77 = " << total_ST2500_rwgt77 << " $\\pm$ " << sqrt(pow(error_ST2500_rwgt77_0SFOS,2) + pow(error_ST2500_rwgt77_1SFOS,2) + pow(error_ST2500_rwgt77_2SFOS,2)) << std::endl;

/*
  std::cout << "total_ST250_rwgt67 = " << total_ST250_rwgt67 << std::endl;
  std::cout << "total_ST250_rwgt68 = " << total_ST250_rwgt68 << std::endl;
  std::cout << "total_ST250_rwgt69 = " << total_ST250_rwgt69 << std::endl;
  std::cout << "total_ST250_rwgt70 = " << total_ST250_rwgt70 << std::endl;
  std::cout << "total_ST250_rwgt71 = " << total_ST250_rwgt71 << std::endl;
  std::cout << "total_ST250_rwgt72 = " << total_ST250_rwgt72 << std::endl;
  std::cout << "total_ST250_rwgt73 = " << total_ST250_rwgt73 << std::endl;
  std::cout << "total_ST250_rwgt74 = " << total_ST250_rwgt74 << std::endl;
  std::cout << "total_ST250_rwgt75 = " << total_ST250_rwgt75 << std::endl;
  std::cout << "total_ST250_rwgt76 = " << total_ST250_rwgt76 << std::endl;
  std::cout << "total_ST250_rwgt77 = " << total_ST250_rwgt77 << std::endl;

  std::cout << "total_ST500_rwgt67 = " << total_ST500_rwgt67 << std::endl;
  std::cout << "total_ST500_rwgt68 = " << total_ST500_rwgt68 << std::endl;
  std::cout << "total_ST500_rwgt69 = " << total_ST500_rwgt69 << std::endl;
  std::cout << "total_ST500_rwgt70 = " << total_ST500_rwgt70 << std::endl;
  std::cout << "total_ST500_rwgt71 = " << total_ST500_rwgt71 << std::endl;
  std::cout << "total_ST500_rwgt72 = " << total_ST500_rwgt72 << std::endl;
  std::cout << "total_ST500_rwgt73 = " << total_ST500_rwgt73 << std::endl;
  std::cout << "total_ST500_rwgt74 = " << total_ST500_rwgt74 << std::endl;
  std::cout << "total_ST500_rwgt75 = " << total_ST500_rwgt75 << std::endl;
  std::cout << "total_ST500_rwgt76 = " << total_ST500_rwgt76 << std::endl;
  std::cout << "total_ST500_rwgt77 = " << total_ST500_rwgt77 << std::endl;

  std::cout << "total_ST750_rwgt67 = " << total_ST750_rwgt67 << std::endl;
  std::cout << "total_ST750_rwgt68 = " << total_ST750_rwgt68 << std::endl;
  std::cout << "total_ST750_rwgt69 = " << total_ST750_rwgt69 << std::endl;
  std::cout << "total_ST750_rwgt70 = " << total_ST750_rwgt70 << std::endl;
  std::cout << "total_ST750_rwgt71 = " << total_ST750_rwgt71 << std::endl;
  std::cout << "total_ST750_rwgt72 = " << total_ST750_rwgt72 << std::endl;
  std::cout << "total_ST750_rwgt73 = " << total_ST750_rwgt73 << std::endl;
  std::cout << "total_ST750_rwgt74 = " << total_ST750_rwgt74 << std::endl;
  std::cout << "total_ST750_rwgt75 = " << total_ST750_rwgt75 << std::endl;
  std::cout << "total_ST750_rwgt76 = " << total_ST750_rwgt76 << std::endl;
  std::cout << "total_ST750_rwgt77 = " << total_ST750_rwgt77 << std::endl;

  std::cout << "total_ST1000_rwgt67 = " << total_ST1000_rwgt67 << std::endl;
  std::cout << "total_ST1000_rwgt68 = " << total_ST1000_rwgt68 << std::endl;
  std::cout << "total_ST1000_rwgt69 = " << total_ST1000_rwgt69 << std::endl;
  std::cout << "total_ST1000_rwgt70 = " << total_ST1000_rwgt70 << std::endl;
  std::cout << "total_ST1000_rwgt71 = " << total_ST1000_rwgt71 << std::endl;
  std::cout << "total_ST1000_rwgt72 = " << total_ST1000_rwgt72 << std::endl;
  std::cout << "total_ST1000_rwgt73 = " << total_ST1000_rwgt73 << std::endl;
  std::cout << "total_ST1000_rwgt74 = " << total_ST1000_rwgt74 << std::endl;
  std::cout << "total_ST1000_rwgt75 = " << total_ST1000_rwgt75 << std::endl;
  std::cout << "total_ST1000_rwgt76 = " << total_ST1000_rwgt76 << std::endl;
  std::cout << "total_ST1000_rwgt77 = " << total_ST1000_rwgt77 << std::endl;

  std::cout << "total_ST1500_rwgt67 = " << total_ST1500_rwgt67 << std::endl;
  std::cout << "total_ST1500_rwgt68 = " << total_ST1500_rwgt68 << std::endl;
  std::cout << "total_ST1500_rwgt69 = " << total_ST1500_rwgt69 << std::endl;
  std::cout << "total_ST1500_rwgt70 = " << total_ST1500_rwgt70 << std::endl;
  std::cout << "total_ST1500_rwgt71 = " << total_ST1500_rwgt71 << std::endl;
  std::cout << "total_ST1500_rwgt72 = " << total_ST1500_rwgt72 << std::endl;
  std::cout << "total_ST1500_rwgt73 = " << total_ST1500_rwgt73 << std::endl;
  std::cout << "total_ST1500_rwgt74 = " << total_ST1500_rwgt74 << std::endl;
  std::cout << "total_ST1500_rwgt75 = " << total_ST1500_rwgt75 << std::endl;
  std::cout << "total_ST1500_rwgt76 = " << total_ST1500_rwgt76 << std::endl;
  std::cout << "total_ST1500_rwgt77 = " << total_ST1500_rwgt77 << std::endl;

  std::cout << "total_ST2000_rwgt67 = " << total_ST2000_rwgt67 << std::endl;
  std::cout << "total_ST2000_rwgt68 = " << total_ST2000_rwgt68 << std::endl;
  std::cout << "total_ST2000_rwgt69 = " << total_ST2000_rwgt69 << std::endl;
  std::cout << "total_ST2000_rwgt70 = " << total_ST2000_rwgt70 << std::endl;
  std::cout << "total_ST2000_rwgt71 = " << total_ST2000_rwgt71 << std::endl;
  std::cout << "total_ST2000_rwgt72 = " << total_ST2000_rwgt72 << std::endl;
  std::cout << "total_ST2000_rwgt73 = " << total_ST2000_rwgt73 << std::endl;
  std::cout << "total_ST2000_rwgt74 = " << total_ST2000_rwgt74 << std::endl;
  std::cout << "total_ST2000_rwgt75 = " << total_ST2000_rwgt75 << std::endl;
  std::cout << "total_ST2000_rwgt76 = " << total_ST2000_rwgt76 << std::endl;
  std::cout << "total_ST2000_rwgt77 = " << total_ST2000_rwgt77 << std::endl;

  std::cout << "total_ST2500_rwgt67 = " << total_ST2500_rwgt67 << std::endl;
  std::cout << "total_ST2500_rwgt68 = " << total_ST2500_rwgt68 << std::endl;
  std::cout << "total_ST2500_rwgt69 = " << total_ST2500_rwgt69 << std::endl;
  std::cout << "total_ST2500_rwgt70 = " << total_ST2500_rwgt70 << std::endl;
  std::cout << "total_ST2500_rwgt71 = " << total_ST2500_rwgt71 << std::endl;
  std::cout << "total_ST2500_rwgt72 = " << total_ST2500_rwgt72 << std::endl;
  std::cout << "total_ST2500_rwgt73 = " << total_ST2500_rwgt73 << std::endl;
  std::cout << "total_ST2500_rwgt74 = " << total_ST2500_rwgt74 << std::endl;
  std::cout << "total_ST2500_rwgt75 = " << total_ST2500_rwgt75 << std::endl;
  std::cout << "total_ST2500_rwgt76 = " << total_ST2500_rwgt76 << std::endl;
  std::cout << "total_ST2500_rwgt77 = " << total_ST2500_rwgt77 << std::endl;
*/

}
