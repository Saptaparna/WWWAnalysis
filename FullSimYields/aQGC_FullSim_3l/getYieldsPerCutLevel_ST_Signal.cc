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
//double BR = 1.0;

double total_noST_rwgt78, total_noST_rwgt79, total_noST_rwgt80, total_noST_rwgt81, total_noST_rwgt82, total_noST_rwgt83, total_noST_rwgt84, total_noST_rwgt85, total_noST_rwgt86, total_noST_rwgt87, total_noST_rwgt88;
double total_ST250_rwgt78, total_ST250_rwgt79, total_ST250_rwgt80, total_ST250_rwgt81, total_ST250_rwgt82, total_ST250_rwgt83, total_ST250_rwgt84, total_ST250_rwgt85, total_ST250_rwgt86, total_ST250_rwgt87, total_ST250_rwgt88;
double total_ST500_rwgt78, total_ST500_rwgt79, total_ST500_rwgt80, total_ST500_rwgt81, total_ST500_rwgt82, total_ST500_rwgt83, total_ST500_rwgt84, total_ST500_rwgt85, total_ST500_rwgt86, total_ST500_rwgt87, total_ST500_rwgt88;
double total_ST750_rwgt78, total_ST750_rwgt79, total_ST750_rwgt80, total_ST750_rwgt81, total_ST750_rwgt82, total_ST750_rwgt83, total_ST750_rwgt84, total_ST750_rwgt85, total_ST750_rwgt86, total_ST750_rwgt87, total_ST750_rwgt88;
double total_ST1000_rwgt78, total_ST1000_rwgt79, total_ST1000_rwgt80, total_ST1000_rwgt81, total_ST1000_rwgt82, total_ST1000_rwgt83, total_ST1000_rwgt84, total_ST1000_rwgt85, total_ST1000_rwgt86, total_ST1000_rwgt87, total_ST1000_rwgt88;
double total_ST1500_rwgt78, total_ST1500_rwgt79, total_ST1500_rwgt80, total_ST1500_rwgt81, total_ST1500_rwgt82, total_ST1500_rwgt83, total_ST1500_rwgt84, total_ST1500_rwgt85, total_ST1500_rwgt86, total_ST1500_rwgt87, total_ST1500_rwgt88;
double total_ST2000_rwgt78, total_ST2000_rwgt79, total_ST2000_rwgt80, total_ST2000_rwgt81, total_ST2000_rwgt82, total_ST2000_rwgt83, total_ST2000_rwgt84, total_ST2000_rwgt85, total_ST2000_rwgt86, total_ST2000_rwgt87, total_ST2000_rwgt88;
double total_ST2500_rwgt78, total_ST2500_rwgt79, total_ST2500_rwgt80, total_ST2500_rwgt81, total_ST2500_rwgt82, total_ST2500_rwgt83, total_ST2500_rwgt84, total_ST2500_rwgt85, total_ST2500_rwgt86, total_ST2500_rwgt87, total_ST2500_rwgt88;

void getYieldsPerCutLevel_0SFOS_noST()
{

  total_noST_rwgt78 = total_noST_rwgt79 = total_noST_rwgt80 = total_noST_rwgt81 = total_noST_rwgt82 = total_noST_rwgt83 = total_noST_rwgt84 = total_noST_rwgt85 = total_noST_rwgt86 = total_noST_rwgt87 = total_noST_rwgt88 = 0.0;
  total_ST250_rwgt78 = total_ST250_rwgt79 = total_ST250_rwgt80 = total_ST250_rwgt81 = total_ST250_rwgt82 = total_ST250_rwgt83 = total_ST250_rwgt84 = total_ST250_rwgt85 = total_ST250_rwgt86 = total_ST250_rwgt87 = total_ST250_rwgt88 = 0.0; 
  total_ST500_rwgt78 = total_ST500_rwgt79 = total_ST500_rwgt80 = total_ST500_rwgt81 = total_ST500_rwgt82 = total_ST500_rwgt83 = total_ST500_rwgt84 = total_ST500_rwgt85 = total_ST500_rwgt86 = total_ST500_rwgt87 = total_ST500_rwgt88 = 0.0;
  total_ST750_rwgt78 = total_ST750_rwgt79 = total_ST750_rwgt80 = total_ST750_rwgt81 = total_ST750_rwgt82 = total_ST750_rwgt83 = total_ST750_rwgt84 = total_ST750_rwgt85 = total_ST750_rwgt86 = total_ST750_rwgt87 = total_ST750_rwgt88 = 0.0; 
  total_ST1000_rwgt78 = total_ST1000_rwgt79 = total_ST1000_rwgt80 = total_ST1000_rwgt81 = total_ST1000_rwgt82 = total_ST1000_rwgt83 = total_ST1000_rwgt84 = total_ST1000_rwgt85 = total_ST1000_rwgt86 = total_ST1000_rwgt87 = total_ST1000_rwgt88 = 0.0;
  total_ST1500_rwgt78 = total_ST1500_rwgt79 = total_ST1500_rwgt80 = total_ST1500_rwgt81 = total_ST1500_rwgt82 = total_ST1500_rwgt83 = total_ST1500_rwgt84 = total_ST1500_rwgt85 = total_ST1500_rwgt86 = total_ST1500_rwgt87 = total_ST1500_rwgt88 = 0.0;
  total_ST2000_rwgt78 = total_ST2000_rwgt79 = total_ST2000_rwgt80 = total_ST2000_rwgt81 = total_ST2000_rwgt82 = total_ST2000_rwgt83 = total_ST2000_rwgt84 = total_ST2000_rwgt85 = total_ST2000_rwgt86 = total_ST2000_rwgt87 = total_ST2000_rwgt88 = 0.0;
  total_ST2500_rwgt78 = total_ST2500_rwgt79 = total_ST2500_rwgt80 = total_ST2500_rwgt81 = total_ST2500_rwgt82 = total_ST2500_rwgt83 = total_ST2500_rwgt84 = total_ST2500_rwgt85 = total_ST2500_rwgt86 = total_ST2500_rwgt87 = total_ST2500_rwgt88 = 0.0;

  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");    
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");
  
  std::cout << "0SFOS channel" << std::endl;
  
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl; 
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts except ST}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_noST_rwgt78 = h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt79 = h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt80 = h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt81 = h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt82 = h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt83 = h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt84 = h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt85 = h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt86 = h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt87 = h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR;
  total_noST_rwgt88 = h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST250()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST250_rwgt78 = h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt79 = h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt80 = h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt81 = h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt82 = h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt83 = h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt84 = h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt85 = h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt86 = h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt87 = h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR;
  total_ST250_rwgt88 = h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST500_rwgt78 = h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt79 = h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt80 = h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt81 = h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt82 = h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt83 = h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt84 = h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt85 = h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt86 = h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt87 = h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR;
  total_ST500_rwgt87 = h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST750()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST750_rwgt78 = h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt79 = h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt80 = h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt81 = h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt82 = h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt83 = h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt84 = h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt85 = h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt86 = h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt87 = h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR;
  total_ST750_rwgt88 = h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST1000()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST1000_rwgt78 = h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt79 = h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt80 = h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt81 = h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt82 = h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt83 = h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt84 = h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt85 = h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt86 = h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt87 = h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR;
  total_ST1000_rwgt88 = h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST1500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt79->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt81->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt82->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST1500_rwgt78 = h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt79 = h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt80 = h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt81 = h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt82 = h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt83 = h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt84 = h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt85 = h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt86 = h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt87 = h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR;
  total_ST1500_rwgt88 = h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST2000()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST2000_rwgt78 = h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt79 = h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt80 = h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt81 = h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt82 = h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt83 = h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt84 = h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt85 = h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt86 = h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt87 = h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR;
  total_ST2000_rwgt88 = h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR;
}

void getYieldsPerCutLevel_0SFOS_ST2500()
{
  TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
  TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt78");
  TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt79");
  TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt80");
  TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt81");
  TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt82");
  TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt83");
  TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt84");
  TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt85");
  TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt86");
  TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt87");
  TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_0SFOS_ST_rwgt88");

  std::cout << "0SFOS channel" << std::endl;

  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_ST2500_rwgt78 = h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt79 = h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt80 = h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt81 = h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt82 = h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt83 = h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt84 = h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt85 = h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt86 = h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt87 = h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR;
  total_ST2500_rwgt88 = h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_noST()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST_rwgt78 += h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt79 += h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt80 += h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt81 += h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt82 += h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt83 += h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt84 += h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt85 += h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt86 += h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt87 += h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt88 += h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_ST250()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
    std::cout << "\\end{table}" << std::endl;

    total_ST250_rwgt78 += h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt79 += h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt80 += h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt81 += h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt82 += h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt83 += h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt84 += h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt85 += h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt86 += h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt87 += h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt88 += h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST500_rwgt78 += h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt79 += h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt80 += h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt81 += h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt82 += h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt83 += h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt84 += h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt85 += h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt86 += h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt87 += h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt88 += h_aQGC_rwgt88->GetBinContent(4)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_ST750()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST750_rwgt78 += h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt79 += h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt80 += h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt81 += h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt82 += h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt83 += h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt84 += h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt85 += h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt86 += h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt87 += h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt88 += h_aQGC_rwgt88->GetBinContent(5)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_ST1000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST1000_rwgt78 += h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt79 += h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt80 += h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt81 += h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt82 += h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt83 += h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt84 += h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt85 += h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt86 += h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt87 += h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt88 += h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_ST1500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt79->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt81->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt82->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1500_rwgt78 += h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt79 += h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt80 += h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt81 += h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt82 += h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt83 += h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt84 += h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt85 += h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt86 += h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt87 += h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt88 += h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR;

}

void getYieldsPerCutLevel_1SFOS_ST2000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2000_rwgt78 += h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt79 += h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt80 += h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt81 += h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt82 += h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt83 += h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt84 += h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt85 += h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt86 += h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt87 += h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt88 += h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR;
}

void getYieldsPerCutLevel_1SFOS_ST2500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_1SFOS_ST_rwgt88");
    
    std::cout << "1SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2500_rwgt78 += h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt79 += h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt80 += h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt81 += h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt82 += h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt83 += h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt84 += h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt85 += h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt86 += h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt87 += h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt88 += h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR;

}

void getYieldsPerCutLevel_2SFOS_noST()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST_rwgt78 += h_aQGC_rwgt78->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt79 += h_aQGC_rwgt79->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt80 += h_aQGC_rwgt80->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt81 += h_aQGC_rwgt81->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt82 += h_aQGC_rwgt82->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt83 += h_aQGC_rwgt83->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt84 += h_aQGC_rwgt84->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt85 += h_aQGC_rwgt85->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt86 += h_aQGC_rwgt86->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt87 += h_aQGC_rwgt87->GetBinContent(2)*LUMI*BR;
    total_noST_rwgt88 += h_aQGC_rwgt88->GetBinContent(2)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST250()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 250}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST250_rwgt78 += h_aQGC_rwgt78->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt79 += h_aQGC_rwgt79->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt80 += h_aQGC_rwgt80->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt81 += h_aQGC_rwgt81->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt82 += h_aQGC_rwgt82->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt83 += h_aQGC_rwgt83->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt84 += h_aQGC_rwgt84->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt85 += h_aQGC_rwgt85->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt86 += h_aQGC_rwgt86->GetBinContent(3)*LUMI*BR;
    total_ST250_rwgt87 += h_aQGC_rwgt87->GetBinContent(3)*LUMI*BR; 
    total_ST500_rwgt88 += h_aQGC_rwgt88->GetBinContent(3)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(4)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST500_rwgt78 += h_aQGC_rwgt78->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt79 += h_aQGC_rwgt79->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt80 += h_aQGC_rwgt80->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt81 += h_aQGC_rwgt81->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt82 += h_aQGC_rwgt82->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt83 += h_aQGC_rwgt83->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt84 += h_aQGC_rwgt84->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt85 += h_aQGC_rwgt85->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt86 += h_aQGC_rwgt86->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt87 += h_aQGC_rwgt87->GetBinContent(4)*LUMI*BR;
    total_ST500_rwgt88 += h_aQGC_rwgt88->GetBinContent(4)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST750()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(5)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 750}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST750_rwgt78 += h_aQGC_rwgt78->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt79 += h_aQGC_rwgt79->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt80 += h_aQGC_rwgt80->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt81 += h_aQGC_rwgt81->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt82 += h_aQGC_rwgt82->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt83 += h_aQGC_rwgt83->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt84 += h_aQGC_rwgt84->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt85 += h_aQGC_rwgt85->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt86 += h_aQGC_rwgt86->GetBinContent(5)*LUMI*BR;
    total_ST750_rwgt87 += h_aQGC_rwgt87->GetBinContent(5)*LUMI*BR; 
    total_ST750_rwgt87 += h_aQGC_rwgt88->GetBinContent(5)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST1000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1000_rwgt78 += h_aQGC_rwgt78->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt79 += h_aQGC_rwgt79->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt80 += h_aQGC_rwgt80->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt81 += h_aQGC_rwgt81->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt82 += h_aQGC_rwgt82->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt83 += h_aQGC_rwgt83->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt84 += h_aQGC_rwgt84->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt85 += h_aQGC_rwgt85->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt86 += h_aQGC_rwgt86->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt87 += h_aQGC_rwgt87->GetBinContent(6)*LUMI*BR;
    total_ST1000_rwgt88 += h_aQGC_rwgt88->GetBinContent(6)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST1500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt79->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt81->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR << " $\\pm$ " << h_aQGC_rwgt82->GetBinError(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST1500_rwgt78 += h_aQGC_rwgt78->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt79 += h_aQGC_rwgt79->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt80 += h_aQGC_rwgt80->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt81 += h_aQGC_rwgt81->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt82 += h_aQGC_rwgt82->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt83 += h_aQGC_rwgt83->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt84 += h_aQGC_rwgt84->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt85 += h_aQGC_rwgt85->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt86 += h_aQGC_rwgt86->GetBinContent(7)*LUMI*BR;
    total_ST1500_rwgt87 += h_aQGC_rwgt87->GetBinContent(7)*LUMI*BR; 
    total_ST1500_rwgt88 += h_aQGC_rwgt88->GetBinContent(7)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST2000()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2000_rwgt78 += h_aQGC_rwgt78->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt79 += h_aQGC_rwgt79->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt80 += h_aQGC_rwgt80->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt81 += h_aQGC_rwgt81->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt82 += h_aQGC_rwgt82->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt83 += h_aQGC_rwgt83->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt84 += h_aQGC_rwgt84->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt85 += h_aQGC_rwgt85->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt86 += h_aQGC_rwgt86->GetBinContent(8)*LUMI*BR;
    total_ST2000_rwgt87 += h_aQGC_rwgt87->GetBinContent(8)*LUMI*BR; 
    total_ST2000_rwgt88 += h_aQGC_rwgt88->GetBinContent(8)*LUMI*BR;
}

void getYieldsPerCutLevel_2SFOS_ST2500()
{
    TFile* fileaQGC = new TFile("output_www_aqgc_sapta_skim_1_1_3l_ft1.root");
    TH1F *h_aQGC_rwgt78 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt78");
    TH1F *h_aQGC_rwgt79 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt79");
    TH1F *h_aQGC_rwgt80 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt80");
    TH1F *h_aQGC_rwgt81 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt81");
    TH1F *h_aQGC_rwgt82 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt82");
    TH1F *h_aQGC_rwgt83 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt83");
    TH1F *h_aQGC_rwgt84 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt84");
    TH1F *h_aQGC_rwgt85 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt85");
    TH1F *h_aQGC_rwgt86 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt86");
    TH1F *h_aQGC_rwgt87 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt87");
    TH1F *h_aQGC_rwgt88 = (TH1F*)fileaQGC->Get("h_TotalEvents_2SFOS_ST_rwgt88");
    
    std::cout << "2SFOS channel" << std::endl;
    
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$: -0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  0.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  1.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.0"              << " & "   << std::setprecision(3) << h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "$aQGC f_{t1}$:  2.5"              << " & "   << std::setprecision(3) << h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST2500_rwgt78 += h_aQGC_rwgt78->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt79 += h_aQGC_rwgt79->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt80 += h_aQGC_rwgt80->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt81 += h_aQGC_rwgt81->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt82 += h_aQGC_rwgt82->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt83 += h_aQGC_rwgt83->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt84 += h_aQGC_rwgt84->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt85 += h_aQGC_rwgt85->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt86 += h_aQGC_rwgt86->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt87 += h_aQGC_rwgt87->GetBinContent(9)*LUMI*BR;
    total_ST2500_rwgt88 += h_aQGC_rwgt88->GetBinContent(9)*LUMI*BR; 
}

void total()
{

  std::cout << "total_noST_rwgt78 = " << total_noST_rwgt78 << std::endl;
  std::cout << "total_noST_rwgt79 = " << total_noST_rwgt79 << std::endl;
  std::cout << "total_noST_rwgt80 = " << total_noST_rwgt80 << std::endl;
  std::cout << "total_noST_rwgt81 = " << total_noST_rwgt81 << std::endl;
  std::cout << "total_noST_rwgt82 = " << total_noST_rwgt82 << std::endl;
  std::cout << "total_noST_rwgt83 = " << total_noST_rwgt83 << std::endl;
  std::cout << "total_noST_rwgt84 = " << total_noST_rwgt84 << std::endl;
  std::cout << "total_noST_rwgt85 = " << total_noST_rwgt85 << std::endl;
  std::cout << "total_noST_rwgt86 = " << total_noST_rwgt86 << std::endl;
  std::cout << "total_noST_rwgt87 = " << total_noST_rwgt87 << std::endl;
  std::cout << "total_noST_rwgt88 = " << total_noST_rwgt88 << std::endl;

  std::cout << "total_ST250_rwgt78 = " << total_ST250_rwgt78 << std::endl;
  std::cout << "total_ST250_rwgt79 = " << total_ST250_rwgt79 << std::endl;
  std::cout << "total_ST250_rwgt80 = " << total_ST250_rwgt80 << std::endl;
  std::cout << "total_ST250_rwgt81 = " << total_ST250_rwgt81 << std::endl;
  std::cout << "total_ST250_rwgt82 = " << total_ST250_rwgt82 << std::endl;
  std::cout << "total_ST250_rwgt83 = " << total_ST250_rwgt83 << std::endl;
  std::cout << "total_ST250_rwgt84 = " << total_ST250_rwgt84 << std::endl;
  std::cout << "total_ST250_rwgt85 = " << total_ST250_rwgt85 << std::endl;
  std::cout << "total_ST250_rwgt86 = " << total_ST250_rwgt86 << std::endl;
  std::cout << "total_ST250_rwgt87 = " << total_ST250_rwgt87 << std::endl;
  std::cout << "total_ST250_rwgt88 = " << total_ST250_rwgt88 << std::endl;

  std::cout << "total_ST500_rwgt78 = " << total_ST500_rwgt78 << std::endl;
  std::cout << "total_ST500_rwgt79 = " << total_ST500_rwgt79 << std::endl;
  std::cout << "total_ST500_rwgt80 = " << total_ST500_rwgt80 << std::endl;
  std::cout << "total_ST500_rwgt81 = " << total_ST500_rwgt81 << std::endl;
  std::cout << "total_ST500_rwgt82 = " << total_ST500_rwgt82 << std::endl;
  std::cout << "total_ST500_rwgt83 = " << total_ST500_rwgt83 << std::endl;
  std::cout << "total_ST500_rwgt84 = " << total_ST500_rwgt84 << std::endl;
  std::cout << "total_ST500_rwgt85 = " << total_ST500_rwgt85 << std::endl;
  std::cout << "total_ST500_rwgt86 = " << total_ST500_rwgt86 << std::endl;
  std::cout << "total_ST500_rwgt87 = " << total_ST500_rwgt87 << std::endl;
  std::cout << "total_ST500_rwgt88 = " << total_ST500_rwgt88 << std::endl;

  std::cout << "total_ST750_rwgt78 = " << total_ST750_rwgt78 << std::endl;
  std::cout << "total_ST750_rwgt79 = " << total_ST750_rwgt79 << std::endl;
  std::cout << "total_ST750_rwgt80 = " << total_ST750_rwgt80 << std::endl;
  std::cout << "total_ST750_rwgt81 = " << total_ST750_rwgt81 << std::endl;
  std::cout << "total_ST750_rwgt82 = " << total_ST750_rwgt82 << std::endl;
  std::cout << "total_ST750_rwgt83 = " << total_ST750_rwgt83 << std::endl;
  std::cout << "total_ST750_rwgt84 = " << total_ST750_rwgt84 << std::endl;
  std::cout << "total_ST750_rwgt85 = " << total_ST750_rwgt85 << std::endl;
  std::cout << "total_ST750_rwgt86 = " << total_ST750_rwgt86 << std::endl;
  std::cout << "total_ST750_rwgt87 = " << total_ST750_rwgt87 << std::endl;
  std::cout << "total_ST750_rwgt88 = " << total_ST750_rwgt88 << std::endl;

  std::cout << "total_ST1000_rwgt78 = " << total_ST1000_rwgt78 << std::endl;
  std::cout << "total_ST1000_rwgt79 = " << total_ST1000_rwgt79 << std::endl;
  std::cout << "total_ST1000_rwgt80 = " << total_ST1000_rwgt80 << std::endl;
  std::cout << "total_ST1000_rwgt81 = " << total_ST1000_rwgt81 << std::endl;
  std::cout << "total_ST1000_rwgt82 = " << total_ST1000_rwgt82 << std::endl;
  std::cout << "total_ST1000_rwgt83 = " << total_ST1000_rwgt83 << std::endl;
  std::cout << "total_ST1000_rwgt84 = " << total_ST1000_rwgt84 << std::endl;
  std::cout << "total_ST1000_rwgt85 = " << total_ST1000_rwgt85 << std::endl;
  std::cout << "total_ST1000_rwgt86 = " << total_ST1000_rwgt86 << std::endl;
  std::cout << "total_ST1000_rwgt87 = " << total_ST1000_rwgt87 << std::endl;
  std::cout << "total_ST1000_rwgt88 = " << total_ST1000_rwgt88 << std::endl;

  std::cout << "total_ST1500_rwgt78 = " << total_ST1500_rwgt78 << std::endl;
  std::cout << "total_ST1500_rwgt79 = " << total_ST1500_rwgt79 << std::endl;
  std::cout << "total_ST1500_rwgt80 = " << total_ST1500_rwgt80 << std::endl;
  std::cout << "total_ST1500_rwgt81 = " << total_ST1500_rwgt81 << std::endl;
  std::cout << "total_ST1500_rwgt82 = " << total_ST1500_rwgt82 << std::endl;
  std::cout << "total_ST1500_rwgt83 = " << total_ST1500_rwgt83 << std::endl;
  std::cout << "total_ST1500_rwgt84 = " << total_ST1500_rwgt84 << std::endl;
  std::cout << "total_ST1500_rwgt85 = " << total_ST1500_rwgt85 << std::endl;
  std::cout << "total_ST1500_rwgt86 = " << total_ST1500_rwgt86 << std::endl;
  std::cout << "total_ST1500_rwgt87 = " << total_ST1500_rwgt87 << std::endl;
  std::cout << "total_ST1500_rwgt88 = " << total_ST1500_rwgt88 << std::endl;

  std::cout << "total_ST2000_rwgt78 = " << total_ST2000_rwgt78 << std::endl;
  std::cout << "total_ST2000_rwgt79 = " << total_ST2000_rwgt79 << std::endl;
  std::cout << "total_ST2000_rwgt80 = " << total_ST2000_rwgt80 << std::endl;
  std::cout << "total_ST2000_rwgt81 = " << total_ST2000_rwgt81 << std::endl;
  std::cout << "total_ST2000_rwgt82 = " << total_ST2000_rwgt82 << std::endl;
  std::cout << "total_ST2000_rwgt83 = " << total_ST2000_rwgt83 << std::endl;
  std::cout << "total_ST2000_rwgt84 = " << total_ST2000_rwgt84 << std::endl;
  std::cout << "total_ST2000_rwgt85 = " << total_ST2000_rwgt85 << std::endl;
  std::cout << "total_ST2000_rwgt86 = " << total_ST2000_rwgt86 << std::endl;
  std::cout << "total_ST2000_rwgt87 = " << total_ST2000_rwgt87 << std::endl;
  std::cout << "total_ST2000_rwgt88 = " << total_ST2000_rwgt88 << std::endl;

  std::cout << "total_ST2500_rwgt78 = " << total_ST2500_rwgt78 << std::endl;
  std::cout << "total_ST2500_rwgt79 = " << total_ST2500_rwgt79 << std::endl;
  std::cout << "total_ST2500_rwgt80 = " << total_ST2500_rwgt80 << std::endl;
  std::cout << "total_ST2500_rwgt81 = " << total_ST2500_rwgt81 << std::endl;
  std::cout << "total_ST2500_rwgt82 = " << total_ST2500_rwgt82 << std::endl;
  std::cout << "total_ST2500_rwgt83 = " << total_ST2500_rwgt83 << std::endl;
  std::cout << "total_ST2500_rwgt84 = " << total_ST2500_rwgt84 << std::endl;
  std::cout << "total_ST2500_rwgt85 = " << total_ST2500_rwgt85 << std::endl;
  std::cout << "total_ST2500_rwgt86 = " << total_ST2500_rwgt86 << std::endl;
  std::cout << "total_ST2500_rwgt87 = " << total_ST2500_rwgt87 << std::endl;
  std::cout << "total_ST2500_rwgt88 = " << total_ST2500_rwgt88 << std::endl;


}
