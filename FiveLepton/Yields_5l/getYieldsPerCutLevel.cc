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

double LUMI2016 = 35.9;
double LUMI2017 = 41.5;
double LUMI2018 = 59.7;

double total_noCut = 0.0;
double total_Ztag = 0.0;

void getYieldsPerCutLevel_4Mu()
{

  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_4Mu = (TH1F*)fileTTW->Get("h_TotalEvents_4Mu");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_4Mu1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_4Mu");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_4Mu2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_4Mu");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_4Mu = (TH1F*)fileTTbar1->Get("h_TotalEvents_4Mu");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_4Mu = (TH1F*)fileDY1050->Get("h_TotalEvents_4Mu");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_4Mu = (TH1F*)fileDY50->Get("h_TotalEvents_4Mu");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_4Mu = (TH1F*)fileZZ1->Get("h_TotalEvents_4Mu");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_4Mu = (TH1F*)fileZZ2->Get("h_TotalEvents_4Mu");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_4Mu = (TH1F*)fileZZ3->Get("h_TotalEvents_4Mu");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_4Mu = (TH1F*)fileZZ4->Get("h_TotalEvents_4Mu");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_4Mu = (TH1F*)fileZZ5->Get("h_TotalEvents_4Mu");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_4Mu = (TH1F*)fileZZ6->Get("h_TotalEvents_4Mu");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_4Mu = (TH1F*)fileZZ7->Get("h_TotalEvents_4Mu");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_4Mu = (TH1F*)fileZZ8->Get("h_TotalEvents_4Mu");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_4Mu = (TH1F*)fileZZ9->Get("h_TotalEvents_4Mu");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_4Mu = (TH1F*)fileZZ10->Get("h_TotalEvents_4Mu");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_4Mu1 = (TH1F*)fileWZ1->Get("h_TotalEvents_4Mu");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_4Mu2 = (TH1F*)fileWZ2->Get("h_TotalEvents_4Mu");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_4Mu = (TH1F*)fileWWG->Get("h_TotalEvents_4Mu");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_4Mu = (TH1F*)fileWZG->Get("h_TotalEvents_4Mu");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_4Mu = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_4Mu");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_4Mu = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_4Mu");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_4Mu = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_4Mu");

  std::cout << "4 Mu channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_4Mu1->GetBinContent(2)*LUMI2018 + hWZ_4Mu2->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_4Mu1->GetBinContent(2)*LUMI2018 + hTTZ_4Mu2->GetBinContent(2)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_4Mu->GetBinContent(2)*LUMI2018 + hZZ2_4Mu->GetBinContent(2)*LUMI2018 + hZZ3_4Mu->GetBinContent(2)*LUMI2018 + hZZ4_4Mu->GetBinContent(2)*LUMI2018 + hZZ5_4Mu->GetBinContent(2)*LUMI2018 + hZZ6_4Mu->GetBinContent(2)*LUMI2018 + hZZ7_4Mu->GetBinContent(2)*LUMI2018 + hZZ8_4Mu->GetBinContent(2)*LUMI2018 + hZZ9_4Mu->GetBinContent(2)*LUMI2018 + hZZ10_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_4Mu->GetBinContent(2)*LUMI2018 + hDY50_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_4Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the four muon channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}

void getYieldsPerCutLevel_4El()
{
  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_4El = (TH1F*)fileTTW->Get("h_TotalEvents_4El");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_4El1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_4El");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_4El2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_4El");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_4El = (TH1F*)fileTTbar1->Get("h_TotalEvents_4El");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_4El = (TH1F*)fileDY1050->Get("h_TotalEvents_4El");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_4El = (TH1F*)fileDY50->Get("h_TotalEvents_4El");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_4El = (TH1F*)fileZZ1->Get("h_TotalEvents_4El");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_4El = (TH1F*)fileZZ2->Get("h_TotalEvents_4El");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_4El = (TH1F*)fileZZ3->Get("h_TotalEvents_4El");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_4El = (TH1F*)fileZZ4->Get("h_TotalEvents_4El");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_4El = (TH1F*)fileZZ5->Get("h_TotalEvents_4El");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_4El = (TH1F*)fileZZ6->Get("h_TotalEvents_4El");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_4El = (TH1F*)fileZZ7->Get("h_TotalEvents_4El");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_4El = (TH1F*)fileZZ8->Get("h_TotalEvents_4El");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_4El = (TH1F*)fileZZ9->Get("h_TotalEvents_4El");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_4El = (TH1F*)fileZZ10->Get("h_TotalEvents_4El");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_4El1 = (TH1F*)fileWZ1->Get("h_TotalEvents_4El");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_4El2 = (TH1F*)fileWZ2->Get("h_TotalEvents_4El");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_4El = (TH1F*)fileWWG->Get("h_TotalEvents_4El");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_4El = (TH1F*)fileWZG->Get("h_TotalEvents_4El");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_4El = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_4El");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_4El = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_4El");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_4El = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_4El");

  std::cout << "4 El channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_4El1->GetBinContent(2)*LUMI2018 + hWZ_4El2->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_4El1->GetBinContent(2)*LUMI2018 + hTTZ_4El2->GetBinContent(2)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_4El->GetBinContent(2)*LUMI2018 + hZZ2_4El->GetBinContent(2)*LUMI2018 + hZZ3_4El->GetBinContent(2)*LUMI2018 + hZZ4_4El->GetBinContent(2)*LUMI2018 + hZZ5_4El->GetBinContent(2)*LUMI2018 + hZZ6_4El->GetBinContent(2)*LUMI2018 + hZZ7_4El->GetBinContent(2)*LUMI2018 + hZZ8_4El->GetBinContent(2)*LUMI2018 + hZZ9_4El->GetBinContent(2)*LUMI2018 + hZZ10_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_4El->GetBinContent(2)*LUMI2018 + hDY50_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_4El->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the four electron channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}

void getYieldsPerCutLevel_2El2Mu()
{
  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_2El2Mu = (TH1F*)fileTTW->Get("h_TotalEvents_2El2Mu");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_2El2Mu1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_2El2Mu");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_2El2Mu2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_2El2Mu");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_2El2Mu = (TH1F*)fileTTbar1->Get("h_TotalEvents_2El2Mu");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_2El2Mu = (TH1F*)fileDY1050->Get("h_TotalEvents_2El2Mu");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_2El2Mu = (TH1F*)fileDY50->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_2El2Mu = (TH1F*)fileZZ1->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_2El2Mu = (TH1F*)fileZZ2->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_2El2Mu = (TH1F*)fileZZ3->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_2El2Mu = (TH1F*)fileZZ4->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_2El2Mu = (TH1F*)fileZZ5->Get("h_TotalEvents_2El2Mu");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_2El2Mu = (TH1F*)fileZZ6->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_2El2Mu = (TH1F*)fileZZ7->Get("h_TotalEvents_2El2Mu");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_2El2Mu = (TH1F*)fileZZ8->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_2El2Mu = (TH1F*)fileZZ9->Get("h_TotalEvents_2El2Mu");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_2El2Mu = (TH1F*)fileZZ10->Get("h_TotalEvents_2El2Mu");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_2El2Mu1 = (TH1F*)fileWZ1->Get("h_TotalEvents_2El2Mu");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_2El2Mu2 = (TH1F*)fileWZ2->Get("h_TotalEvents_2El2Mu");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_2El2Mu = (TH1F*)fileWWG->Get("h_TotalEvents_2El2Mu");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_2El2Mu = (TH1F*)fileWZG->Get("h_TotalEvents_2El2Mu");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_2El2Mu = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_2El2Mu");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_2El2Mu = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_2El2Mu");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_2El2Mu = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_2El2Mu");

  std::cout << "2 El 2Mu channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_2El2Mu1->GetBinContent(2)*LUMI2018 + hWZ_2El2Mu2->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_2El2Mu1->GetBinContent(2)*LUMI2018 + hTTZ_2El2Mu2->GetBinContent(2)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ2_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ3_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ4_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ5_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ6_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ7_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ8_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ9_2El2Mu->GetBinContent(2)*LUMI2018 + hZZ10_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_2El2Mu->GetBinContent(2)*LUMI2018 + hDY50_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_2El2Mu->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the two electron and two muon channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}
