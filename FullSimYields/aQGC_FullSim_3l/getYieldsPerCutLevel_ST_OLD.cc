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

double total_noST;
double total_ST1000;
double total_ST1500;
double total_ST2000;
double total_ST2500;

void getYieldsPerCutLevel_0SFOS_noST()
{
  total_noST=total_ST1000=total_ST1500=total_ST2000=total_ST2500=0.0;
 
  TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
  TH1F *hTTW_0SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
  TH1F *hTTZ_0SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
  TH1F *hTTbar1_0SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
  TH1F *hTTbar2_0SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
  TH1F *hTTbar3_0SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
  TH1F *hDY1050_0SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_0SFOS_ST");
  
  TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
  TH1F *hDY50_0SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root"); 
  TH1F *hWJets1_0SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
  TH1F *hWJets2_0SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_0SFOS_ST");
  
  TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
  TH1F *hWJets3_0SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
  TH1F *hWJets4_0SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
  TH1F *hWJets5_0SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
  TH1F *hWJets6_0SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
  TH1F *hWJets7_0SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root"); 
  TH1F *hWWJJ_0SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
  TH1F *hZZ1_0SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
  TH1F *hZZ2_0SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_0SFOS_ST");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
  TH1F *hZZ3_0SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_0SFOS_ST");
 
  TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
  TH1F *hWZ_0SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_0SFOS_ST");

  TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
  TH1F *hSignal_0SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_0SFOS_ST");

  std::cout << "0SFOS channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_0SFOS->GetBinContent(2)*LUMI + hZZ2_0SFOS->GetBinContent(2)*LUMI + hZZ3_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_0SFOS->GetBinContent(2)*LUMI + hWJets2_0SFOS->GetBinContent(2)*LUMI + hWJets3_0SFOS->GetBinContent(2)*LUMI + hWJets4_0SFOS->GetBinContent(2)*LUMI + hWJets5_0SFOS->GetBinContent(2)*LUMI + hWJets6_0SFOS->GetBinContent(2)*LUMI + hWJets7_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "WW"                << " & "   << hWWJJ_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_0SFOS->GetBinContent(2)*LUMI + hTTbar2_0SFOS->GetBinContent(2)*LUMI + hTTbar3_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_0SFOS->GetBinContent(2)*LUMI + hDY50_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(2)*LUMI+hTTW_0SFOS->GetBinContent(2)*LUMI+hTTZ_0SFOS->GetBinContent(2)*LUMI+hZZ1_0SFOS->GetBinContent(2)*LUMI + hZZ2_0SFOS->GetBinContent(2)*LUMI + hZZ3_0SFOS->GetBinContent(2)*LUMI+hWJets1_0SFOS->GetBinContent(2)*LUMI + hWJets2_0SFOS->GetBinContent(2)*LUMI + hWJets3_0SFOS->GetBinContent(2)*LUMI + hWJets4_0SFOS->GetBinContent(2)*LUMI + hWJets5_0SFOS->GetBinContent(2)*LUMI + hWJets6_0SFOS->GetBinContent(2)*LUMI + hWJets7_0SFOS->GetBinContent(2)*LUMI+hWWJJ_0SFOS->GetBinContent(2)*LUMI+hTTbar1_0SFOS->GetBinContent(2)*LUMI + hTTbar2_0SFOS->GetBinContent(2)*LUMI + hTTbar3_0SFOS->GetBinContent(2)*LUMI +hDY1050_0SFOS->GetBinContent(2)*LUMI + hDY50_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl; 
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM WWW"    << " & "   << hSignal_0SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts except ST}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  total_noST = hWZ_0SFOS->GetBinContent(2)*LUMI+hTTW_0SFOS->GetBinContent(2)*LUMI+hTTZ_0SFOS->GetBinContent(2)*LUMI+hZZ1_0SFOS->GetBinContent(2)*LUMI + hZZ2_0SFOS->GetBinContent(2)*LUMI + hZZ3_0SFOS->GetBinContent(2)*LUMI+hWJets1_0SFOS->GetBinContent(2)*LUMI + hWJets2_0SFOS->GetBinContent(2)*LUMI + hWJets3_0SFOS->GetBinContent(2)*LUMI + hWJets4_0SFOS->GetBinContent(2)*LUMI + hWJets5_0SFOS->GetBinContent(2)*LUMI + hWJets6_0SFOS->GetBinContent(2)*LUMI + hWJets7_0SFOS->GetBinContent(2)*LUMI+hWWJJ_0SFOS->GetBinContent(2)*LUMI+hTTbar1_0SFOS->GetBinContent(2)*LUMI + hTTbar2_0SFOS->GetBinContent(2)*LUMI + hTTbar3_0SFOS->GetBinContent(2)*LUMI +hDY1050_0SFOS->GetBinContent(2)*LUMI + hDY50_0SFOS->GetBinContent(2)*LUMI; 

}

void getYieldsPerCutLevel_0SFOS_ST1000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_0SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_0SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_0SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_0SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_0SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_0SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_0SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_0SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_0SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_0SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_0SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_0SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_0SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_0SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_0SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_0SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_0SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_0SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_0SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_0SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_0SFOS_ST");
    
    std::cout << "0SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_0SFOS->GetBinContent(3)*LUMI + hZZ2_0SFOS->GetBinContent(3)*LUMI + hZZ3_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_0SFOS->GetBinContent(3)*LUMI + hWJets2_0SFOS->GetBinContent(3)*LUMI + hWJets3_0SFOS->GetBinContent(3)*LUMI + hWJets4_0SFOS->GetBinContent(3)*LUMI + hWJets5_0SFOS->GetBinContent(3)*LUMI + hWJets6_0SFOS->GetBinContent(3)*LUMI + hWJets7_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_0SFOS->GetBinContent(3)*LUMI + hTTbar2_0SFOS->GetBinContent(3)*LUMI + hTTbar3_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_0SFOS->GetBinContent(3)*LUMI + hDY50_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(3)*LUMI+hTTW_0SFOS->GetBinContent(3)*LUMI+hTTZ_0SFOS->GetBinContent(3)*LUMI+hZZ1_0SFOS->GetBinContent(3)*LUMI + hZZ2_0SFOS->GetBinContent(3)*LUMI + hZZ3_0SFOS->GetBinContent(3)*LUMI+hWJets1_0SFOS->GetBinContent(3)*LUMI + hWJets2_0SFOS->GetBinContent(3)*LUMI + hWJets3_0SFOS->GetBinContent(3)*LUMI + hWJets4_0SFOS->GetBinContent(3)*LUMI + hWJets5_0SFOS->GetBinContent(3)*LUMI + hWJets6_0SFOS->GetBinContent(3)*LUMI + hWJets7_0SFOS->GetBinContent(3)*LUMI+hWWJJ_0SFOS->GetBinContent(3)*LUMI+hTTbar1_0SFOS->GetBinContent(3)*LUMI + hTTbar2_0SFOS->GetBinContent(3)*LUMI + hTTbar3_0SFOS->GetBinContent(3)*LUMI +hDY1050_0SFOS->GetBinContent(3)*LUMI + hDY50_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_0SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1000 = hWZ_0SFOS->GetBinContent(3)*LUMI+hTTW_0SFOS->GetBinContent(3)*LUMI+hTTZ_0SFOS->GetBinContent(3)*LUMI+hZZ1_0SFOS->GetBinContent(3)*LUMI + hZZ2_0SFOS->GetBinContent(3)*LUMI + hZZ3_0SFOS->GetBinContent(3)*LUMI+hWJets1_0SFOS->GetBinContent(3)*LUMI + hWJets2_0SFOS->GetBinContent(3)*LUMI + hWJets3_0SFOS->GetBinContent(3)*LUMI + hWJets4_0SFOS->GetBinContent(3)*LUMI + hWJets5_0SFOS->GetBinContent(3)*LUMI + hWJets6_0SFOS->GetBinContent(3)*LUMI + hWJets7_0SFOS->GetBinContent(3)*LUMI+hWWJJ_0SFOS->GetBinContent(3)*LUMI+hTTbar1_0SFOS->GetBinContent(3)*LUMI + hTTbar2_0SFOS->GetBinContent(3)*LUMI + hTTbar3_0SFOS->GetBinContent(3)*LUMI +hDY1050_0SFOS->GetBinContent(3)*LUMI + hDY50_0SFOS->GetBinContent(3)*LUMI;
}


void getYieldsPerCutLevel_0SFOS_ST1500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_0SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_0SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_0SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_0SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_0SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_0SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_0SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_0SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_0SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_0SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_0SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_0SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_0SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_0SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_0SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_0SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_0SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_0SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_0SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_0SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_0SFOS_ST");
    
    std::cout << "0SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_0SFOS->GetBinContent(4)*LUMI + hZZ2_0SFOS->GetBinContent(4)*LUMI + hZZ3_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_0SFOS->GetBinContent(4)*LUMI + hWJets2_0SFOS->GetBinContent(4)*LUMI + hWJets3_0SFOS->GetBinContent(4)*LUMI + hWJets4_0SFOS->GetBinContent(4)*LUMI + hWJets5_0SFOS->GetBinContent(4)*LUMI + hWJets6_0SFOS->GetBinContent(4)*LUMI + hWJets7_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_0SFOS->GetBinContent(4)*LUMI + hTTbar2_0SFOS->GetBinContent(4)*LUMI + hTTbar3_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_0SFOS->GetBinContent(4)*LUMI + hDY50_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(4)*LUMI+hTTW_0SFOS->GetBinContent(4)*LUMI+hTTZ_0SFOS->GetBinContent(4)*LUMI+hZZ1_0SFOS->GetBinContent(4)*LUMI + hZZ2_0SFOS->GetBinContent(4)*LUMI + hZZ3_0SFOS->GetBinContent(4)*LUMI+hWJets1_0SFOS->GetBinContent(4)*LUMI + hWJets2_0SFOS->GetBinContent(4)*LUMI + hWJets3_0SFOS->GetBinContent(4)*LUMI + hWJets4_0SFOS->GetBinContent(4)*LUMI + hWJets5_0SFOS->GetBinContent(4)*LUMI + hWJets6_0SFOS->GetBinContent(4)*LUMI + hWJets7_0SFOS->GetBinContent(4)*LUMI+hWWJJ_0SFOS->GetBinContent(4)*LUMI+hTTbar1_0SFOS->GetBinContent(4)*LUMI + hTTbar2_0SFOS->GetBinContent(4)*LUMI + hTTbar3_0SFOS->GetBinContent(4)*LUMI +hDY1050_0SFOS->GetBinContent(4)*LUMI + hDY50_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_0SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
   
    total_ST1500 = hWZ_0SFOS->GetBinContent(4)*LUMI+hTTW_0SFOS->GetBinContent(4)*LUMI+hTTZ_0SFOS->GetBinContent(4)*LUMI+hZZ1_0SFOS->GetBinContent(4)*LUMI + hZZ2_0SFOS->GetBinContent(4)*LUMI + hZZ3_0SFOS->GetBinContent(4)*LUMI+hWJets1_0SFOS->GetBinContent(4)*LUMI + hWJets2_0SFOS->GetBinContent(4)*LUMI + hWJets3_0SFOS->GetBinContent(4)*LUMI + hWJets4_0SFOS->GetBinContent(4)*LUMI + hWJets5_0SFOS->GetBinContent(4)*LUMI + hWJets6_0SFOS->GetBinContent(4)*LUMI + hWJets7_0SFOS->GetBinContent(4)*LUMI+hWWJJ_0SFOS->GetBinContent(4)*LUMI+hTTbar1_0SFOS->GetBinContent(4)*LUMI + hTTbar2_0SFOS->GetBinContent(4)*LUMI + hTTbar3_0SFOS->GetBinContent(4)*LUMI +hDY1050_0SFOS->GetBinContent(4)*LUMI + hDY50_0SFOS->GetBinContent(4)*LUMI; 
}

void getYieldsPerCutLevel_0SFOS_ST2000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_0SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_0SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_0SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_0SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_0SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_0SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_0SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_0SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_0SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_0SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_0SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_0SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_0SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_0SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_0SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_0SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_0SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_0SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_0SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_0SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_0SFOS_ST");
    
    std::cout << "0SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_0SFOS->GetBinContent(5)*LUMI + hZZ2_0SFOS->GetBinContent(5)*LUMI + hZZ3_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_0SFOS->GetBinContent(5)*LUMI + hWJets2_0SFOS->GetBinContent(5)*LUMI + hWJets3_0SFOS->GetBinContent(5)*LUMI + hWJets4_0SFOS->GetBinContent(5)*LUMI + hWJets5_0SFOS->GetBinContent(5)*LUMI + hWJets6_0SFOS->GetBinContent(5)*LUMI + hWJets7_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_0SFOS->GetBinContent(5)*LUMI + hTTbar2_0SFOS->GetBinContent(5)*LUMI + hTTbar3_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_0SFOS->GetBinContent(5)*LUMI + hDY50_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(5)*LUMI+hTTW_0SFOS->GetBinContent(5)*LUMI+hTTZ_0SFOS->GetBinContent(5)*LUMI+hZZ1_0SFOS->GetBinContent(5)*LUMI + hZZ2_0SFOS->GetBinContent(5)*LUMI + hZZ3_0SFOS->GetBinContent(5)*LUMI+hWJets1_0SFOS->GetBinContent(5)*LUMI + hWJets2_0SFOS->GetBinContent(5)*LUMI + hWJets3_0SFOS->GetBinContent(5)*LUMI + hWJets4_0SFOS->GetBinContent(5)*LUMI + hWJets5_0SFOS->GetBinContent(5)*LUMI + hWJets6_0SFOS->GetBinContent(5)*LUMI + hWJets7_0SFOS->GetBinContent(5)*LUMI+hWWJJ_0SFOS->GetBinContent(5)*LUMI+hTTbar1_0SFOS->GetBinContent(5)*LUMI + hTTbar2_0SFOS->GetBinContent(5)*LUMI + hTTbar3_0SFOS->GetBinContent(5)*LUMI +hDY1050_0SFOS->GetBinContent(5)*LUMI + hDY50_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_0SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2000 = hWZ_0SFOS->GetBinContent(5)*LUMI+hTTW_0SFOS->GetBinContent(5)*LUMI+hTTZ_0SFOS->GetBinContent(5)*LUMI+hZZ1_0SFOS->GetBinContent(5)*LUMI + hZZ2_0SFOS->GetBinContent(5)*LUMI + hZZ3_0SFOS->GetBinContent(5)*LUMI+hWJets1_0SFOS->GetBinContent(5)*LUMI + hWJets2_0SFOS->GetBinContent(5)*LUMI + hWJets3_0SFOS->GetBinContent(5)*LUMI + hWJets4_0SFOS->GetBinContent(5)*LUMI + hWJets5_0SFOS->GetBinContent(5)*LUMI + hWJets6_0SFOS->GetBinContent(5)*LUMI + hWJets7_0SFOS->GetBinContent(5)*LUMI+hWWJJ_0SFOS->GetBinContent(5)*LUMI+hTTbar1_0SFOS->GetBinContent(5)*LUMI + hTTbar2_0SFOS->GetBinContent(5)*LUMI + hTTbar3_0SFOS->GetBinContent(5)*LUMI +hDY1050_0SFOS->GetBinContent(5)*LUMI + hDY50_0SFOS->GetBinContent(5)*LUMI;
}

void getYieldsPerCutLevel_0SFOS_ST2500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_0SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_0SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_0SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_0SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_0SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_0SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_0SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_0SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_0SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_0SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_0SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_0SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_0SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_0SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_0SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_0SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_0SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_0SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_0SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_0SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_0SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_0SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_0SFOS_ST");
    
    std::cout << "0SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_0SFOS->GetBinContent(6)*LUMI + hZZ2_0SFOS->GetBinContent(6)*LUMI + hZZ3_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_0SFOS->GetBinContent(6)*LUMI + hWJets2_0SFOS->GetBinContent(6)*LUMI + hWJets3_0SFOS->GetBinContent(6)*LUMI + hWJets4_0SFOS->GetBinContent(6)*LUMI + hWJets5_0SFOS->GetBinContent(6)*LUMI + hWJets6_0SFOS->GetBinContent(6)*LUMI + hWJets7_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_0SFOS->GetBinContent(6)*LUMI + hTTbar2_0SFOS->GetBinContent(6)*LUMI + hTTbar3_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_0SFOS->GetBinContent(6)*LUMI + hDY50_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_0SFOS->GetBinContent(6)*LUMI+hTTW_0SFOS->GetBinContent(6)*LUMI+hTTZ_0SFOS->GetBinContent(6)*LUMI+hZZ1_0SFOS->GetBinContent(6)*LUMI + hZZ2_0SFOS->GetBinContent(6)*LUMI + hZZ3_0SFOS->GetBinContent(6)*LUMI+hWJets1_0SFOS->GetBinContent(6)*LUMI + hWJets2_0SFOS->GetBinContent(6)*LUMI + hWJets3_0SFOS->GetBinContent(6)*LUMI + hWJets4_0SFOS->GetBinContent(6)*LUMI + hWJets5_0SFOS->GetBinContent(6)*LUMI + hWJets6_0SFOS->GetBinContent(6)*LUMI + hWJets7_0SFOS->GetBinContent(6)*LUMI+hWWJJ_0SFOS->GetBinContent(6)*LUMI+hTTbar1_0SFOS->GetBinContent(6)*LUMI + hTTbar2_0SFOS->GetBinContent(6)*LUMI + hTTbar3_0SFOS->GetBinContent(6)*LUMI +hDY1050_0SFOS->GetBinContent(6)*LUMI + hDY50_0SFOS->GetBinContent(6)*LUMI << " \\\\ "  << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_0SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 0SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2500 = hWZ_0SFOS->GetBinContent(6)*LUMI+hTTW_0SFOS->GetBinContent(6)*LUMI+hTTZ_0SFOS->GetBinContent(6)*LUMI+hZZ1_0SFOS->GetBinContent(6)*LUMI + hZZ2_0SFOS->GetBinContent(6)*LUMI + hZZ3_0SFOS->GetBinContent(6)*LUMI+hWJets1_0SFOS->GetBinContent(6)*LUMI + hWJets2_0SFOS->GetBinContent(6)*LUMI + hWJets3_0SFOS->GetBinContent(6)*LUMI + hWJets4_0SFOS->GetBinContent(6)*LUMI + hWJets5_0SFOS->GetBinContent(6)*LUMI + hWJets6_0SFOS->GetBinContent(6)*LUMI + hWJets7_0SFOS->GetBinContent(6)*LUMI+hWWJJ_0SFOS->GetBinContent(6)*LUMI+hTTbar1_0SFOS->GetBinContent(6)*LUMI + hTTbar2_0SFOS->GetBinContent(6)*LUMI + hTTbar3_0SFOS->GetBinContent(6)*LUMI +hDY1050_0SFOS->GetBinContent(6)*LUMI + hDY50_0SFOS->GetBinContent(6)*LUMI;
}


void getYieldsPerCutLevel_1SFOS_noST()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_1SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_1SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_1SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_1SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_1SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_1SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_1SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_1SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_1SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_1SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_1SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_1SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_1SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_1SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_1SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_1SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_1SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_1SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_1SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_1SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_1SFOS_ST");
    
    std::cout << "1SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_1SFOS->GetBinContent(2)*LUMI + hZZ2_1SFOS->GetBinContent(2)*LUMI + hZZ3_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_1SFOS->GetBinContent(2)*LUMI + hWJets2_1SFOS->GetBinContent(2)*LUMI + hWJets3_1SFOS->GetBinContent(2)*LUMI + hWJets4_1SFOS->GetBinContent(2)*LUMI + hWJets5_1SFOS->GetBinContent(2)*LUMI + hWJets6_1SFOS->GetBinContent(2)*LUMI + hWJets7_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_1SFOS->GetBinContent(2)*LUMI + hTTbar2_1SFOS->GetBinContent(2)*LUMI + hTTbar3_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_1SFOS->GetBinContent(2)*LUMI + hDY50_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(2)*LUMI+hTTW_1SFOS->GetBinContent(2)*LUMI+hTTZ_1SFOS->GetBinContent(2)*LUMI+hZZ1_1SFOS->GetBinContent(2)*LUMI + hZZ2_1SFOS->GetBinContent(2)*LUMI + hZZ3_1SFOS->GetBinContent(2)*LUMI+hWJets1_1SFOS->GetBinContent(2)*LUMI + hWJets2_1SFOS->GetBinContent(2)*LUMI + hWJets3_1SFOS->GetBinContent(2)*LUMI + hWJets4_1SFOS->GetBinContent(2)*LUMI + hWJets5_1SFOS->GetBinContent(2)*LUMI + hWJets6_1SFOS->GetBinContent(2)*LUMI + hWJets7_1SFOS->GetBinContent(2)*LUMI+hWWJJ_1SFOS->GetBinContent(2)*LUMI+hTTbar1_1SFOS->GetBinContent(2)*LUMI + hTTbar2_1SFOS->GetBinContent(2)*LUMI + hTTbar3_1SFOS->GetBinContent(2)*LUMI +hDY1050_1SFOS->GetBinContent(2)*LUMI + hDY50_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_1SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST += hWZ_1SFOS->GetBinContent(2)*LUMI+hTTW_1SFOS->GetBinContent(2)*LUMI+hTTZ_1SFOS->GetBinContent(2)*LUMI+hZZ1_1SFOS->GetBinContent(2)*LUMI + hZZ2_1SFOS->GetBinContent(2)*LUMI + hZZ3_1SFOS->GetBinContent(2)*LUMI+hWJets1_1SFOS->GetBinContent(2)*LUMI + hWJets2_1SFOS->GetBinContent(2)*LUMI + hWJets3_1SFOS->GetBinContent(2)*LUMI + hWJets4_1SFOS->GetBinContent(2)*LUMI + hWJets5_1SFOS->GetBinContent(2)*LUMI + hWJets6_1SFOS->GetBinContent(2)*LUMI + hWJets7_1SFOS->GetBinContent(2)*LUMI+hWWJJ_1SFOS->GetBinContent(2)*LUMI+hTTbar1_1SFOS->GetBinContent(2)*LUMI + hTTbar2_1SFOS->GetBinContent(2)*LUMI + hTTbar3_1SFOS->GetBinContent(2)*LUMI +hDY1050_1SFOS->GetBinContent(2)*LUMI + hDY50_1SFOS->GetBinContent(2)*LUMI;
    
}

void getYieldsPerCutLevel_1SFOS_ST1000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_1SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_1SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_1SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_1SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_1SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_1SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_1SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_1SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_1SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_1SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_1SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_1SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_1SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_1SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_1SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_1SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_1SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_1SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_1SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_1SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_1SFOS_ST");
    
    std::cout << "1SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_1SFOS->GetBinContent(3)*LUMI + hZZ2_1SFOS->GetBinContent(3)*LUMI + hZZ3_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_1SFOS->GetBinContent(3)*LUMI + hWJets2_1SFOS->GetBinContent(3)*LUMI + hWJets3_1SFOS->GetBinContent(3)*LUMI + hWJets4_1SFOS->GetBinContent(3)*LUMI + hWJets5_1SFOS->GetBinContent(3)*LUMI + hWJets6_1SFOS->GetBinContent(3)*LUMI + hWJets7_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_1SFOS->GetBinContent(3)*LUMI + hTTbar2_1SFOS->GetBinContent(3)*LUMI + hTTbar3_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_1SFOS->GetBinContent(3)*LUMI + hDY50_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(3)*LUMI+hTTW_1SFOS->GetBinContent(3)*LUMI+hTTZ_1SFOS->GetBinContent(3)*LUMI+hZZ1_1SFOS->GetBinContent(3)*LUMI + hZZ2_1SFOS->GetBinContent(3)*LUMI + hZZ3_1SFOS->GetBinContent(3)*LUMI+hWJets1_1SFOS->GetBinContent(3)*LUMI + hWJets2_1SFOS->GetBinContent(3)*LUMI + hWJets3_1SFOS->GetBinContent(3)*LUMI + hWJets4_1SFOS->GetBinContent(3)*LUMI + hWJets5_1SFOS->GetBinContent(3)*LUMI + hWJets6_1SFOS->GetBinContent(3)*LUMI + hWJets7_1SFOS->GetBinContent(3)*LUMI+hWWJJ_1SFOS->GetBinContent(3)*LUMI+hTTbar1_1SFOS->GetBinContent(3)*LUMI + hTTbar2_1SFOS->GetBinContent(3)*LUMI + hTTbar3_1SFOS->GetBinContent(3)*LUMI +hDY1050_1SFOS->GetBinContent(3)*LUMI + hDY50_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_1SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1000 += hWZ_1SFOS->GetBinContent(3)*LUMI+hTTW_1SFOS->GetBinContent(3)*LUMI+hTTZ_1SFOS->GetBinContent(3)*LUMI+hZZ1_1SFOS->GetBinContent(3)*LUMI + hZZ2_1SFOS->GetBinContent(3)*LUMI + hZZ3_1SFOS->GetBinContent(3)*LUMI+hWJets1_1SFOS->GetBinContent(3)*LUMI + hWJets2_1SFOS->GetBinContent(3)*LUMI + hWJets3_1SFOS->GetBinContent(3)*LUMI + hWJets4_1SFOS->GetBinContent(3)*LUMI + hWJets5_1SFOS->GetBinContent(3)*LUMI + hWJets6_1SFOS->GetBinContent(3)*LUMI + hWJets7_1SFOS->GetBinContent(3)*LUMI+hWWJJ_1SFOS->GetBinContent(3)*LUMI+hTTbar1_1SFOS->GetBinContent(3)*LUMI + hTTbar2_1SFOS->GetBinContent(3)*LUMI + hTTbar3_1SFOS->GetBinContent(3)*LUMI +hDY1050_1SFOS->GetBinContent(3)*LUMI + hDY50_1SFOS->GetBinContent(3)*LUMI;
}


void getYieldsPerCutLevel_1SFOS_ST1500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_1SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_1SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_1SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_1SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_1SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_1SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_1SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_1SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_1SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_1SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_1SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_1SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_1SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_1SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_1SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_1SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_1SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_1SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_1SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_1SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_1SFOS_ST");
    
    std::cout << "1SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_1SFOS->GetBinContent(4)*LUMI + hZZ2_1SFOS->GetBinContent(4)*LUMI + hZZ3_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_1SFOS->GetBinContent(4)*LUMI + hWJets2_1SFOS->GetBinContent(4)*LUMI + hWJets3_1SFOS->GetBinContent(4)*LUMI + hWJets4_1SFOS->GetBinContent(4)*LUMI + hWJets5_1SFOS->GetBinContent(4)*LUMI + hWJets6_1SFOS->GetBinContent(4)*LUMI + hWJets7_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_1SFOS->GetBinContent(4)*LUMI + hTTbar2_1SFOS->GetBinContent(4)*LUMI + hTTbar3_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_1SFOS->GetBinContent(4)*LUMI + hDY50_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(4)*LUMI+hTTW_1SFOS->GetBinContent(4)*LUMI+hTTZ_1SFOS->GetBinContent(4)*LUMI+hZZ1_1SFOS->GetBinContent(4)*LUMI + hZZ2_1SFOS->GetBinContent(4)*LUMI + hZZ3_1SFOS->GetBinContent(4)*LUMI+hWJets1_1SFOS->GetBinContent(4)*LUMI + hWJets2_1SFOS->GetBinContent(4)*LUMI + hWJets3_1SFOS->GetBinContent(4)*LUMI + hWJets4_1SFOS->GetBinContent(4)*LUMI + hWJets5_1SFOS->GetBinContent(4)*LUMI + hWJets6_1SFOS->GetBinContent(4)*LUMI + hWJets7_1SFOS->GetBinContent(4)*LUMI+hWWJJ_1SFOS->GetBinContent(4)*LUMI+hTTbar1_1SFOS->GetBinContent(4)*LUMI + hTTbar2_1SFOS->GetBinContent(4)*LUMI + hTTbar3_1SFOS->GetBinContent(4)*LUMI +hDY1050_1SFOS->GetBinContent(4)*LUMI + hDY50_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_1SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1500 += hWZ_1SFOS->GetBinContent(4)*LUMI+hTTW_1SFOS->GetBinContent(4)*LUMI+hTTZ_1SFOS->GetBinContent(4)*LUMI+hZZ1_1SFOS->GetBinContent(4)*LUMI + hZZ2_1SFOS->GetBinContent(4)*LUMI + hZZ3_1SFOS->GetBinContent(4)*LUMI+hWJets1_1SFOS->GetBinContent(4)*LUMI + hWJets2_1SFOS->GetBinContent(4)*LUMI + hWJets3_1SFOS->GetBinContent(4)*LUMI + hWJets4_1SFOS->GetBinContent(4)*LUMI + hWJets5_1SFOS->GetBinContent(4)*LUMI + hWJets6_1SFOS->GetBinContent(4)*LUMI + hWJets7_1SFOS->GetBinContent(4)*LUMI+hWWJJ_1SFOS->GetBinContent(4)*LUMI+hTTbar1_1SFOS->GetBinContent(4)*LUMI + hTTbar2_1SFOS->GetBinContent(4)*LUMI + hTTbar3_1SFOS->GetBinContent(4)*LUMI +hDY1050_1SFOS->GetBinContent(4)*LUMI + hDY50_1SFOS->GetBinContent(4)*LUMI;
}

void getYieldsPerCutLevel_1SFOS_ST2000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_1SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_1SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_1SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_1SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_1SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_1SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_1SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_1SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_1SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_1SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_1SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_1SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_1SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_1SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_1SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_1SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_1SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_1SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_1SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_1SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_1SFOS_ST");
    
    std::cout << "1SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_1SFOS->GetBinContent(5)*LUMI + hZZ2_1SFOS->GetBinContent(5)*LUMI + hZZ3_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_1SFOS->GetBinContent(5)*LUMI + hWJets2_1SFOS->GetBinContent(5)*LUMI + hWJets3_1SFOS->GetBinContent(5)*LUMI + hWJets4_1SFOS->GetBinContent(5)*LUMI + hWJets5_1SFOS->GetBinContent(5)*LUMI + hWJets6_1SFOS->GetBinContent(5)*LUMI + hWJets7_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_1SFOS->GetBinContent(5)*LUMI + hTTbar2_1SFOS->GetBinContent(5)*LUMI + hTTbar3_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_1SFOS->GetBinContent(5)*LUMI + hDY50_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(5)*LUMI+hTTW_1SFOS->GetBinContent(5)*LUMI+hTTZ_1SFOS->GetBinContent(5)*LUMI+hZZ1_1SFOS->GetBinContent(5)*LUMI + hZZ2_1SFOS->GetBinContent(5)*LUMI + hZZ3_1SFOS->GetBinContent(5)*LUMI+hWJets1_1SFOS->GetBinContent(5)*LUMI + hWJets2_1SFOS->GetBinContent(5)*LUMI + hWJets3_1SFOS->GetBinContent(5)*LUMI + hWJets4_1SFOS->GetBinContent(5)*LUMI + hWJets5_1SFOS->GetBinContent(5)*LUMI + hWJets6_1SFOS->GetBinContent(5)*LUMI + hWJets7_1SFOS->GetBinContent(5)*LUMI+hWWJJ_1SFOS->GetBinContent(5)*LUMI+hTTbar1_1SFOS->GetBinContent(5)*LUMI + hTTbar2_1SFOS->GetBinContent(5)*LUMI + hTTbar3_1SFOS->GetBinContent(5)*LUMI +hDY1050_1SFOS->GetBinContent(5)*LUMI + hDY50_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_1SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2000 += hWZ_1SFOS->GetBinContent(5)*LUMI+hTTW_1SFOS->GetBinContent(5)*LUMI+hTTZ_1SFOS->GetBinContent(5)*LUMI+hZZ1_1SFOS->GetBinContent(5)*LUMI + hZZ2_1SFOS->GetBinContent(5)*LUMI + hZZ3_1SFOS->GetBinContent(5)*LUMI+hWJets1_1SFOS->GetBinContent(5)*LUMI + hWJets2_1SFOS->GetBinContent(5)*LUMI + hWJets3_1SFOS->GetBinContent(5)*LUMI + hWJets4_1SFOS->GetBinContent(5)*LUMI + hWJets5_1SFOS->GetBinContent(5)*LUMI + hWJets6_1SFOS->GetBinContent(5)*LUMI + hWJets7_1SFOS->GetBinContent(5)*LUMI+hWWJJ_1SFOS->GetBinContent(5)*LUMI+hTTbar1_1SFOS->GetBinContent(5)*LUMI + hTTbar2_1SFOS->GetBinContent(5)*LUMI + hTTbar3_1SFOS->GetBinContent(5)*LUMI +hDY1050_1SFOS->GetBinContent(5)*LUMI + hDY50_1SFOS->GetBinContent(5)*LUMI;
}

void getYieldsPerCutLevel_1SFOS_ST2500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_1SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_1SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_1SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_1SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_1SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_1SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_1SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_1SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_1SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_1SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_1SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_1SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_1SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_1SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_1SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_1SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_1SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_1SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_1SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_1SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_1SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_1SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_1SFOS_ST");
    
    std::cout << "1SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_1SFOS->GetBinContent(6)*LUMI + hZZ2_1SFOS->GetBinContent(6)*LUMI + hZZ3_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_1SFOS->GetBinContent(6)*LUMI + hWJets2_1SFOS->GetBinContent(6)*LUMI + hWJets3_1SFOS->GetBinContent(6)*LUMI + hWJets4_1SFOS->GetBinContent(6)*LUMI + hWJets5_1SFOS->GetBinContent(6)*LUMI + hWJets6_1SFOS->GetBinContent(6)*LUMI + hWJets7_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_1SFOS->GetBinContent(6)*LUMI + hTTbar2_1SFOS->GetBinContent(6)*LUMI + hTTbar3_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_1SFOS->GetBinContent(6)*LUMI + hDY50_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_1SFOS->GetBinContent(6)*LUMI+hTTW_1SFOS->GetBinContent(6)*LUMI+hTTZ_1SFOS->GetBinContent(6)*LUMI+hZZ1_1SFOS->GetBinContent(6)*LUMI + hZZ2_1SFOS->GetBinContent(6)*LUMI + hZZ3_1SFOS->GetBinContent(6)*LUMI+hWJets1_1SFOS->GetBinContent(6)*LUMI + hWJets2_1SFOS->GetBinContent(6)*LUMI + hWJets3_1SFOS->GetBinContent(6)*LUMI + hWJets4_1SFOS->GetBinContent(6)*LUMI + hWJets5_1SFOS->GetBinContent(6)*LUMI + hWJets6_1SFOS->GetBinContent(6)*LUMI + hWJets7_1SFOS->GetBinContent(6)*LUMI+hWWJJ_1SFOS->GetBinContent(6)*LUMI+hTTbar1_1SFOS->GetBinContent(6)*LUMI + hTTbar2_1SFOS->GetBinContent(6)*LUMI + hTTbar3_1SFOS->GetBinContent(6)*LUMI +hDY1050_1SFOS->GetBinContent(6)*LUMI + hDY50_1SFOS->GetBinContent(6)*LUMI << " \\\\ "  << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_1SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 1SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2500 += hWZ_1SFOS->GetBinContent(6)*LUMI+hTTW_1SFOS->GetBinContent(6)*LUMI+hTTZ_1SFOS->GetBinContent(6)*LUMI+hZZ1_1SFOS->GetBinContent(6)*LUMI + hZZ2_1SFOS->GetBinContent(6)*LUMI + hZZ3_1SFOS->GetBinContent(6)*LUMI+hWJets1_1SFOS->GetBinContent(6)*LUMI + hWJets2_1SFOS->GetBinContent(6)*LUMI + hWJets3_1SFOS->GetBinContent(6)*LUMI + hWJets4_1SFOS->GetBinContent(6)*LUMI + hWJets5_1SFOS->GetBinContent(6)*LUMI + hWJets6_1SFOS->GetBinContent(6)*LUMI + hWJets7_1SFOS->GetBinContent(6)*LUMI+hWWJJ_1SFOS->GetBinContent(6)*LUMI+hTTbar1_1SFOS->GetBinContent(6)*LUMI + hTTbar2_1SFOS->GetBinContent(6)*LUMI + hTTbar3_1SFOS->GetBinContent(6)*LUMI +hDY1050_1SFOS->GetBinContent(6)*LUMI + hDY50_1SFOS->GetBinContent(6)*LUMI;
}


void getYieldsPerCutLevel_2SFOS_noST()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_2SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_2SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_2SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_2SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_2SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_2SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_2SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_2SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_2SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_2SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_2SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_2SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_2SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_2SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_2SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_2SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_2SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_2SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_2SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_2SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_2SFOS_ST");
    
    std::cout << "2SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_2SFOS->GetBinContent(2)*LUMI + hZZ2_2SFOS->GetBinContent(2)*LUMI + hZZ3_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_2SFOS->GetBinContent(2)*LUMI + hWJets2_2SFOS->GetBinContent(2)*LUMI + hWJets3_2SFOS->GetBinContent(2)*LUMI + hWJets4_2SFOS->GetBinContent(2)*LUMI + hWJets5_2SFOS->GetBinContent(2)*LUMI + hWJets6_2SFOS->GetBinContent(2)*LUMI + hWJets7_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_2SFOS->GetBinContent(2)*LUMI + hTTbar2_2SFOS->GetBinContent(2)*LUMI + hTTbar3_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_2SFOS->GetBinContent(2)*LUMI + hDY50_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(2)*LUMI+hTTW_2SFOS->GetBinContent(2)*LUMI+hTTZ_2SFOS->GetBinContent(2)*LUMI+hZZ1_2SFOS->GetBinContent(2)*LUMI + hZZ2_2SFOS->GetBinContent(2)*LUMI + hZZ3_2SFOS->GetBinContent(2)*LUMI+hWJets1_2SFOS->GetBinContent(2)*LUMI + hWJets2_2SFOS->GetBinContent(2)*LUMI + hWJets3_2SFOS->GetBinContent(2)*LUMI + hWJets4_2SFOS->GetBinContent(2)*LUMI + hWJets5_2SFOS->GetBinContent(2)*LUMI + hWJets6_2SFOS->GetBinContent(2)*LUMI + hWJets7_2SFOS->GetBinContent(2)*LUMI+hWWJJ_2SFOS->GetBinContent(2)*LUMI+hTTbar1_2SFOS->GetBinContent(2)*LUMI + hTTbar2_2SFOS->GetBinContent(2)*LUMI + hTTbar3_2SFOS->GetBinContent(2)*LUMI +hDY1050_2SFOS->GetBinContent(2)*LUMI + hDY50_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_2SFOS->GetBinContent(2)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts except ST}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_noST += hWZ_2SFOS->GetBinContent(2)*LUMI+hTTW_2SFOS->GetBinContent(2)*LUMI+hTTZ_2SFOS->GetBinContent(2)*LUMI+hZZ1_2SFOS->GetBinContent(2)*LUMI + hZZ2_2SFOS->GetBinContent(2)*LUMI + hZZ3_2SFOS->GetBinContent(2)*LUMI+hWJets1_2SFOS->GetBinContent(2)*LUMI + hWJets2_2SFOS->GetBinContent(2)*LUMI + hWJets3_2SFOS->GetBinContent(2)*LUMI + hWJets4_2SFOS->GetBinContent(2)*LUMI + hWJets5_2SFOS->GetBinContent(2)*LUMI + hWJets6_2SFOS->GetBinContent(2)*LUMI + hWJets7_2SFOS->GetBinContent(2)*LUMI+hWWJJ_2SFOS->GetBinContent(2)*LUMI+hTTbar1_2SFOS->GetBinContent(2)*LUMI + hTTbar2_2SFOS->GetBinContent(2)*LUMI + hTTbar3_2SFOS->GetBinContent(2)*LUMI +hDY1050_2SFOS->GetBinContent(2)*LUMI + hDY50_2SFOS->GetBinContent(2)*LUMI;
    
}

void getYieldsPerCutLevel_2SFOS_ST1000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_2SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_2SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_2SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_2SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_2SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_2SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_2SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_2SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_2SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_2SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_2SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_2SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_2SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_2SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_2SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_2SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_2SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_2SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_2SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_2SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_2SFOS_ST");
    
    std::cout << "2SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_2SFOS->GetBinContent(3)*LUMI + hZZ2_2SFOS->GetBinContent(3)*LUMI + hZZ3_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_2SFOS->GetBinContent(3)*LUMI + hWJets2_2SFOS->GetBinContent(3)*LUMI + hWJets3_2SFOS->GetBinContent(3)*LUMI + hWJets4_2SFOS->GetBinContent(3)*LUMI + hWJets5_2SFOS->GetBinContent(3)*LUMI + hWJets6_2SFOS->GetBinContent(3)*LUMI + hWJets7_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_2SFOS->GetBinContent(3)*LUMI + hTTbar2_2SFOS->GetBinContent(3)*LUMI + hTTbar3_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_2SFOS->GetBinContent(3)*LUMI + hDY50_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(3)*LUMI+hTTW_2SFOS->GetBinContent(3)*LUMI+hTTZ_2SFOS->GetBinContent(3)*LUMI+hZZ1_2SFOS->GetBinContent(3)*LUMI + hZZ2_2SFOS->GetBinContent(3)*LUMI + hZZ3_2SFOS->GetBinContent(3)*LUMI+hWJets1_2SFOS->GetBinContent(3)*LUMI + hWJets2_2SFOS->GetBinContent(3)*LUMI + hWJets3_2SFOS->GetBinContent(3)*LUMI + hWJets4_2SFOS->GetBinContent(3)*LUMI + hWJets5_2SFOS->GetBinContent(3)*LUMI + hWJets6_2SFOS->GetBinContent(3)*LUMI + hWJets7_2SFOS->GetBinContent(3)*LUMI+hWWJJ_2SFOS->GetBinContent(3)*LUMI+hTTbar1_2SFOS->GetBinContent(3)*LUMI + hTTbar2_2SFOS->GetBinContent(3)*LUMI + hTTbar3_2SFOS->GetBinContent(3)*LUMI +hDY1050_2SFOS->GetBinContent(3)*LUMI + hDY50_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_2SFOS->GetBinContent(3)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1000 += hWZ_2SFOS->GetBinContent(6)*LUMI+hTTW_2SFOS->GetBinContent(6)*LUMI+hTTZ_2SFOS->GetBinContent(6)*LUMI+hZZ1_2SFOS->GetBinContent(6)*LUMI + hZZ2_2SFOS->GetBinContent(6)*LUMI + hZZ3_2SFOS->GetBinContent(6)*LUMI+hWJets1_2SFOS->GetBinContent(6)*LUMI + hWJets2_2SFOS->GetBinContent(6)*LUMI + hWJets3_2SFOS->GetBinContent(6)*LUMI + hWJets4_2SFOS->GetBinContent(6)*LUMI + hWJets5_2SFOS->GetBinContent(6)*LUMI + hWJets6_2SFOS->GetBinContent(6)*LUMI + hWJets7_2SFOS->GetBinContent(6)*LUMI+hWWJJ_2SFOS->GetBinContent(6)*LUMI+hTTbar1_2SFOS->GetBinContent(6)*LUMI + hTTbar2_2SFOS->GetBinContent(6)*LUMI + hTTbar3_2SFOS->GetBinContent(6)*LUMI +hDY1050_2SFOS->GetBinContent(6)*LUMI + hDY50_2SFOS->GetBinContent(6)*LUMI;
}


void getYieldsPerCutLevel_2SFOS_ST1500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_2SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_2SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_2SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_2SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_2SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_2SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_2SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_2SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_2SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_2SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_2SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_2SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_2SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_2SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_2SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_2SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_2SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_2SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_2SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_2SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_2SFOS_ST");
    
    std::cout << "2SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_2SFOS->GetBinContent(4)*LUMI + hZZ2_2SFOS->GetBinContent(4)*LUMI + hZZ3_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_2SFOS->GetBinContent(4)*LUMI + hWJets2_2SFOS->GetBinContent(4)*LUMI + hWJets3_2SFOS->GetBinContent(4)*LUMI + hWJets4_2SFOS->GetBinContent(4)*LUMI + hWJets5_2SFOS->GetBinContent(4)*LUMI + hWJets6_2SFOS->GetBinContent(4)*LUMI + hWJets7_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"       << " & "   << hTTbar1_2SFOS->GetBinContent(4)*LUMI + hTTbar2_2SFOS->GetBinContent(4)*LUMI + hTTbar3_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_2SFOS->GetBinContent(4)*LUMI + hDY50_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(4)*LUMI+hTTW_2SFOS->GetBinContent(4)*LUMI+hTTZ_2SFOS->GetBinContent(4)*LUMI+hZZ1_2SFOS->GetBinContent(4)*LUMI + hZZ2_2SFOS->GetBinContent(4)*LUMI + hZZ3_2SFOS->GetBinContent(4)*LUMI+hWJets1_2SFOS->GetBinContent(4)*LUMI + hWJets2_2SFOS->GetBinContent(4)*LUMI + hWJets3_2SFOS->GetBinContent(4)*LUMI + hWJets4_2SFOS->GetBinContent(4)*LUMI + hWJets5_2SFOS->GetBinContent(4)*LUMI + hWJets6_2SFOS->GetBinContent(4)*LUMI + hWJets7_2SFOS->GetBinContent(4)*LUMI+hWWJJ_2SFOS->GetBinContent(4)*LUMI+hTTbar1_2SFOS->GetBinContent(4)*LUMI + hTTbar2_2SFOS->GetBinContent(4)*LUMI + hTTbar3_2SFOS->GetBinContent(4)*LUMI +hDY1050_2SFOS->GetBinContent(4)*LUMI + hDY50_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_2SFOS->GetBinContent(4)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 1500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST1500 += hWZ_2SFOS->GetBinContent(4)*LUMI+hTTW_2SFOS->GetBinContent(4)*LUMI+hTTZ_2SFOS->GetBinContent(4)*LUMI+hZZ1_2SFOS->GetBinContent(4)*LUMI + hZZ2_2SFOS->GetBinContent(4)*LUMI + hZZ3_2SFOS->GetBinContent(4)*LUMI+hWJets1_2SFOS->GetBinContent(4)*LUMI + hWJets2_2SFOS->GetBinContent(4)*LUMI + hWJets3_2SFOS->GetBinContent(4)*LUMI + hWJets4_2SFOS->GetBinContent(4)*LUMI + hWJets5_2SFOS->GetBinContent(4)*LUMI + hWJets6_2SFOS->GetBinContent(4)*LUMI + hWJets7_2SFOS->GetBinContent(4)*LUMI+hWWJJ_2SFOS->GetBinContent(4)*LUMI+hTTbar1_2SFOS->GetBinContent(4)*LUMI + hTTbar2_2SFOS->GetBinContent(4)*LUMI + hTTbar3_2SFOS->GetBinContent(4)*LUMI +hDY1050_2SFOS->GetBinContent(4)*LUMI + hDY50_2SFOS->GetBinContent(4)*LUMI;
}

void getYieldsPerCutLevel_2SFOS_ST2000()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_2SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_2SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_2SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_2SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_2SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_2SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_2SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_2SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_2SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_2SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_2SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_2SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_2SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_2SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_2SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_2SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_2SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_2SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_2SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_2SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_2SFOS_ST");
    
    std::cout << "2SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_2SFOS->GetBinContent(5)*LUMI + hZZ2_2SFOS->GetBinContent(5)*LUMI + hZZ3_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_2SFOS->GetBinContent(5)*LUMI + hWJets2_2SFOS->GetBinContent(5)*LUMI + hWJets3_2SFOS->GetBinContent(5)*LUMI + hWJets4_2SFOS->GetBinContent(5)*LUMI + hWJets5_2SFOS->GetBinContent(5)*LUMI + hWJets6_2SFOS->GetBinContent(5)*LUMI + hWJets7_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_2SFOS->GetBinContent(5)*LUMI + hTTbar2_2SFOS->GetBinContent(5)*LUMI + hTTbar3_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_2SFOS->GetBinContent(5)*LUMI + hDY50_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(5)*LUMI+hTTW_2SFOS->GetBinContent(5)*LUMI+hTTZ_2SFOS->GetBinContent(5)*LUMI+hZZ1_2SFOS->GetBinContent(5)*LUMI + hZZ2_2SFOS->GetBinContent(5)*LUMI + hZZ3_2SFOS->GetBinContent(5)*LUMI+hWJets1_2SFOS->GetBinContent(5)*LUMI + hWJets2_2SFOS->GetBinContent(5)*LUMI + hWJets3_2SFOS->GetBinContent(5)*LUMI + hWJets4_2SFOS->GetBinContent(5)*LUMI + hWJets5_2SFOS->GetBinContent(5)*LUMI + hWJets6_2SFOS->GetBinContent(5)*LUMI + hWJets7_2SFOS->GetBinContent(5)*LUMI+hWWJJ_2SFOS->GetBinContent(5)*LUMI+hTTbar1_2SFOS->GetBinContent(5)*LUMI + hTTbar2_2SFOS->GetBinContent(5)*LUMI + hTTbar3_2SFOS->GetBinContent(5)*LUMI +hDY1050_2SFOS->GetBinContent(5)*LUMI + hDY50_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_2SFOS->GetBinContent(5)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2000}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2000 += hWZ_2SFOS->GetBinContent(5)*LUMI+hTTW_2SFOS->GetBinContent(5)*LUMI+hTTZ_2SFOS->GetBinContent(5)*LUMI+hZZ1_2SFOS->GetBinContent(5)*LUMI + hZZ2_2SFOS->GetBinContent(5)*LUMI + hZZ3_2SFOS->GetBinContent(5)*LUMI+hWJets1_2SFOS->GetBinContent(5)*LUMI + hWJets2_2SFOS->GetBinContent(5)*LUMI + hWJets3_2SFOS->GetBinContent(5)*LUMI + hWJets4_2SFOS->GetBinContent(5)*LUMI + hWJets5_2SFOS->GetBinContent(5)*LUMI + hWJets6_2SFOS->GetBinContent(5)*LUMI + hWJets7_2SFOS->GetBinContent(5)*LUMI+hWWJJ_2SFOS->GetBinContent(5)*LUMI+hTTbar1_2SFOS->GetBinContent(5)*LUMI + hTTbar2_2SFOS->GetBinContent(5)*LUMI + hTTbar3_2SFOS->GetBinContent(5)*LUMI +hDY1050_2SFOS->GetBinContent(5)*LUMI + hDY50_2SFOS->GetBinContent(5)*LUMI;
}

void getYieldsPerCutLevel_2SFOS_ST2500()
{
    TFile* fileTTW = new TFile("output_ttw_incl_mgmlm_skim_1_1.root");
    TH1F *hTTW_2SFOS = (TH1F*)fileTTW->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTZ = new TFile("output_ttz_incl_mgmlm_skim_1_1.root");
    TH1F *hTTZ_2SFOS = (TH1F*)fileTTZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar1 = new TFile("output_ttbar_1ltbr_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar1_2SFOS = (TH1F*)fileTTbar1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar2 = new TFile("output_ttbar_1ltop_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar2_2SFOS = (TH1F*)fileTTbar2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileTTbar3 = new TFile("output_ttbar_dilep_mgmlm_ext1_skim_1_1.root");
    TH1F *hTTbar3_2SFOS = (TH1F*)fileTTbar3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY1050 = new TFile("output_dy_m1050_mgmlm_skim_1_1.root");
    TH1F *hDY1050_2SFOS = (TH1F*)fileDY1050->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileDY50 = new TFile("output_dy_m50_mgmlm_ext1_skim_1_1.root");
    TH1F *hDY50_2SFOS = (TH1F*)fileDY50->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets1 = new TFile("output_wjets_ht100_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets1_2SFOS = (TH1F*)fileWJets1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets2 = new TFile("output_wjets_ht200_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets2_2SFOS = (TH1F*)fileWJets2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets3 = new TFile("output_wjets_ht400_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets3_2SFOS = (TH1F*)fileWJets3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets4 = new TFile("output_wjets_ht600_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets4_2SFOS = (TH1F*)fileWJets4->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets5 = new TFile("output_wjets_ht800_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets5_2SFOS = (TH1F*)fileWJets5->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets6 = new TFile("output_wjets_ht1200_mgmlm_nonext_skim_1_1.root");
    TH1F *hWJets6_2SFOS = (TH1F*)fileWJets6->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWJets7 = new TFile("output_wjets_ht2500_mgmlm_ext1_skim_1_1.root");
    TH1F *hWJets7_2SFOS = (TH1F*)fileWJets7->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWWJJ = new TFile("output_wpwpjj_ewk-qcd_madgraph_skim_1_1.root");
    TH1F *hWWJJ_2SFOS = (TH1F*)fileWWJJ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ1 = new TFile("output_zz_2l2n_powheg_skim_1_1.root");
    TH1F *hZZ1_2SFOS = (TH1F*)fileZZ1->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_skim_1_1.root");
    TH1F *hZZ2_2SFOS = (TH1F*)fileZZ2->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileZZ3 = new TFile("output_zz_4l_powheg_skim_1_1.root");
    TH1F *hZZ3_2SFOS = (TH1F*)fileZZ3->Get("h_TotalEvents_2SFOS_ST");
    
    TFile* fileWZ = new TFile("output_wz_3lnu_powheg_skim_1_1.root");
    TH1F *hWZ_2SFOS = (TH1F*)fileWZ->Get("h_TotalEvents_2SFOS_ST");
    
    TFile *fileSignal = new TFile("output_www_2l_ext1_mia_skim_1_1.root");
    TH1F *hSignal_2SFOS = (TH1F*)fileSignal->Get("h_TotalEvents_2SFOS_ST");
    
    std::cout << "2SFOS channel" << std::endl;
    std::cout << "\\begin{table}[htb]" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(2) << hTTZ_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$ZZ$"              << " & "   << std::setprecision(2) << hZZ1_2SFOS->GetBinContent(6)*LUMI + hZZ2_2SFOS->GetBinContent(6)*LUMI + hZZ3_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "W+jets"            << " & "   << std::setprecision(2) << hWJets1_2SFOS->GetBinContent(6)*LUMI + hWJets2_2SFOS->GetBinContent(6)*LUMI + hWJets3_2SFOS->GetBinContent(6)*LUMI + hWJets4_2SFOS->GetBinContent(6)*LUMI + hWJets5_2SFOS->GetBinContent(6)*LUMI + hWJets6_2SFOS->GetBinContent(6)*LUMI + hWJets7_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "WW"                << " & "   << hWWJJ_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "$t\\bar{t}$"         << " & "   << hTTbar1_2SFOS->GetBinContent(6)*LUMI + hTTbar2_2SFOS->GetBinContent(6)*LUMI + hTTbar3_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "DY"                << " & "   << hDY1050_2SFOS->GetBinContent(6)*LUMI + hDY50_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "Total background"  << " & "   << std::setprecision(3) << hWZ_2SFOS->GetBinContent(6)*LUMI+hTTW_2SFOS->GetBinContent(6)*LUMI+hTTZ_2SFOS->GetBinContent(6)*LUMI+hZZ1_2SFOS->GetBinContent(6)*LUMI + hZZ2_2SFOS->GetBinContent(6)*LUMI + hZZ3_2SFOS->GetBinContent(6)*LUMI+hWJets1_2SFOS->GetBinContent(6)*LUMI + hWJets2_2SFOS->GetBinContent(6)*LUMI + hWJets3_2SFOS->GetBinContent(6)*LUMI + hWJets4_2SFOS->GetBinContent(6)*LUMI + hWJets5_2SFOS->GetBinContent(6)*LUMI + hWJets6_2SFOS->GetBinContent(6)*LUMI + hWJets7_2SFOS->GetBinContent(6)*LUMI+hWWJJ_2SFOS->GetBinContent(6)*LUMI+hTTbar1_2SFOS->GetBinContent(6)*LUMI + hTTbar2_2SFOS->GetBinContent(6)*LUMI + hTTbar3_2SFOS->GetBinContent(6)*LUMI +hDY1050_2SFOS->GetBinContent(6)*LUMI + hDY50_2SFOS->GetBinContent(6)*LUMI << " \\\\ "  << std::endl;
    std::cout << "Signal: SM WWW"    << " & "   << hSignal_2SFOS->GetBinContent(6)*LUMI << " \\\\ " << std::endl;
    std::cout << "\\hline\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\end{center}" << std::endl;
    std::cout << "\\caption{Yield table for the 2SFOS channel after the application of all cuts and ST $>$ 2500}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    
    total_ST2500 += hWZ_2SFOS->GetBinContent(6)*LUMI+hTTW_2SFOS->GetBinContent(6)*LUMI+hTTZ_2SFOS->GetBinContent(6)*LUMI+hZZ1_2SFOS->GetBinContent(6)*LUMI + hZZ2_2SFOS->GetBinContent(6)*LUMI + hZZ3_2SFOS->GetBinContent(6)*LUMI+hWJets1_2SFOS->GetBinContent(6)*LUMI + hWJets2_2SFOS->GetBinContent(6)*LUMI + hWJets3_2SFOS->GetBinContent(6)*LUMI + hWJets4_2SFOS->GetBinContent(6)*LUMI + hWJets5_2SFOS->GetBinContent(6)*LUMI + hWJets6_2SFOS->GetBinContent(6)*LUMI + hWJets7_2SFOS->GetBinContent(6)*LUMI+hWWJJ_2SFOS->GetBinContent(6)*LUMI+hTTbar1_2SFOS->GetBinContent(6)*LUMI + hTTbar2_2SFOS->GetBinContent(6)*LUMI + hTTbar3_2SFOS->GetBinContent(6)*LUMI +hDY1050_2SFOS->GetBinContent(6)*LUMI + hDY50_2SFOS->GetBinContent(6)*LUMI;
}


void total()
{

  std::cout << "total_noST = " << total_noST << std::endl;
  std::cout << "total_ST1000 = " << total_ST1000 << std::endl;
  std::cout << "total_ST1500 = " << total_ST1500 << std::endl;
  std::cout << "total_ST2000 = " << total_ST2000 << std::endl;
  std::cout << "total_ST2500 = " << total_ST2500 << std::endl;

}
