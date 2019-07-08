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
double LUMI2018 = 59.7+41.5+35.9;

double total_noCut = 0.0;
double total_Ztag = 0.0;

void getYieldsPerCutLevel_5l()
{

  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_5l = (TH1F*)fileTTW->Get("h_TotalEvents_5l");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_5l1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_5l");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_5l2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_5l");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_5l = (TH1F*)fileTTbar1->Get("h_TotalEvents_5l");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_5l = (TH1F*)fileDY1050->Get("h_TotalEvents_5l");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_5l = (TH1F*)fileDY50->Get("h_TotalEvents_5l");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_5l = (TH1F*)fileZZ1->Get("h_TotalEvents_5l");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_5l = (TH1F*)fileZZ2->Get("h_TotalEvents_5l");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_5l = (TH1F*)fileZZ3->Get("h_TotalEvents_5l");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_5l = (TH1F*)fileZZ4->Get("h_TotalEvents_5l");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_5l = (TH1F*)fileZZ5->Get("h_TotalEvents_5l");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_5l = (TH1F*)fileZZ6->Get("h_TotalEvents_5l");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_5l = (TH1F*)fileZZ7->Get("h_TotalEvents_5l");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_5l = (TH1F*)fileZZ8->Get("h_TotalEvents_5l");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_5l = (TH1F*)fileZZ9->Get("h_TotalEvents_5l");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_5l = (TH1F*)fileZZ10->Get("h_TotalEvents_5l");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_5l1 = (TH1F*)fileWZ1->Get("h_TotalEvents_5l");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_5l2 = (TH1F*)fileWZ2->Get("h_TotalEvents_5l");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_5l = (TH1F*)fileWWG->Get("h_TotalEvents_5l");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_5l = (TH1F*)fileWZG->Get("h_TotalEvents_5l");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_5l = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_5l = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_5l = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_5l");

  std::cout << "5 l channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_5l1->GetBinContent(2)*LUMI2018 + hWZ_5l2->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_5l1->GetBinContent(2)*LUMI2018 + hTTZ_5l2->GetBinContent(2)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_5l->GetBinContent(2)*LUMI2018 + hZZ2_5l->GetBinContent(2)*LUMI2018 + hZZ3_5l->GetBinContent(2)*LUMI2018 + hZZ4_5l->GetBinContent(2)*LUMI2018 + hZZ5_5l->GetBinContent(2)*LUMI2018 + hZZ6_5l->GetBinContent(2)*LUMI2018 + hZZ7_5l->GetBinContent(2)*LUMI2018 + hZZ8_5l->GetBinContent(2)*LUMI2018 + hZZ9_5l->GetBinContent(2)*LUMI2018 + hZZ10_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_5l->GetBinContent(2)*LUMI2018 + hDY50_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM wzz"    << " & "   << hSignal_WZZ_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl; 
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_5l->GetBinContent(2)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the five leptons channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}

void getYieldsPerCutLevel_5l_bVeto()
{

  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_5l = (TH1F*)fileTTW->Get("h_TotalEvents_5l");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_5l1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_5l");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_5l2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_5l");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_5l = (TH1F*)fileTTbar1->Get("h_TotalEvents_5l");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_5l = (TH1F*)fileDY1050->Get("h_TotalEvents_5l");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_5l = (TH1F*)fileDY50->Get("h_TotalEvents_5l");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_5l = (TH1F*)fileZZ1->Get("h_TotalEvents_5l");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_5l = (TH1F*)fileZZ2->Get("h_TotalEvents_5l");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_5l = (TH1F*)fileZZ3->Get("h_TotalEvents_5l");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_5l = (TH1F*)fileZZ4->Get("h_TotalEvents_5l");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_5l = (TH1F*)fileZZ5->Get("h_TotalEvents_5l");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_5l = (TH1F*)fileZZ6->Get("h_TotalEvents_5l");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_5l = (TH1F*)fileZZ7->Get("h_TotalEvents_5l");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_5l = (TH1F*)fileZZ8->Get("h_TotalEvents_5l");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_5l = (TH1F*)fileZZ9->Get("h_TotalEvents_5l");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_5l = (TH1F*)fileZZ10->Get("h_TotalEvents_5l");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_5l1 = (TH1F*)fileWZ1->Get("h_TotalEvents_5l");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_5l2 = (TH1F*)fileWZ2->Get("h_TotalEvents_5l");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_5l = (TH1F*)fileWWG->Get("h_TotalEvents_5l");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_5l = (TH1F*)fileWZG->Get("h_TotalEvents_5l");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_5l = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_5l = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_5l = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_5l");

  std::cout << "5 l channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_5l1->GetBinContent(3)*LUMI2018 + hWZ_5l2->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_5l1->GetBinContent(3)*LUMI2018 + hTTZ_5l2->GetBinContent(3)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_5l->GetBinContent(3)*LUMI2018 + hZZ2_5l->GetBinContent(3)*LUMI2018 + hZZ3_5l->GetBinContent(3)*LUMI2018 + hZZ4_5l->GetBinContent(3)*LUMI2018 + hZZ5_5l->GetBinContent(3)*LUMI2018 + hZZ6_5l->GetBinContent(3)*LUMI2018 + hZZ7_5l->GetBinContent(3)*LUMI2018 + hZZ8_5l->GetBinContent(3)*LUMI2018 + hZZ9_5l->GetBinContent(3)*LUMI2018 + hZZ10_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_5l->GetBinContent(3)*LUMI2018 + hDY50_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM wzz"    << " & "   << hSignal_WZZ_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl; 
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_5l->GetBinContent(3)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the five leptons channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}

void getYieldsPerCutLevel_5l_bVeto_MET()
{

  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_5l = (TH1F*)fileTTW->Get("h_TotalEvents_5l");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_5l1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_5l");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_5l2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_5l");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_5l = (TH1F*)fileTTbar1->Get("h_TotalEvents_5l");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_5l = (TH1F*)fileDY1050->Get("h_TotalEvents_5l");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_5l = (TH1F*)fileDY50->Get("h_TotalEvents_5l");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_5l = (TH1F*)fileZZ1->Get("h_TotalEvents_5l");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_5l = (TH1F*)fileZZ2->Get("h_TotalEvents_5l");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_5l = (TH1F*)fileZZ3->Get("h_TotalEvents_5l");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_5l = (TH1F*)fileZZ4->Get("h_TotalEvents_5l");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_5l = (TH1F*)fileZZ5->Get("h_TotalEvents_5l");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_5l = (TH1F*)fileZZ6->Get("h_TotalEvents_5l");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_5l = (TH1F*)fileZZ7->Get("h_TotalEvents_5l");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_5l = (TH1F*)fileZZ8->Get("h_TotalEvents_5l");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_5l = (TH1F*)fileZZ9->Get("h_TotalEvents_5l");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_5l = (TH1F*)fileZZ10->Get("h_TotalEvents_5l");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_5l1 = (TH1F*)fileWZ1->Get("h_TotalEvents_5l");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_5l2 = (TH1F*)fileWZ2->Get("h_TotalEvents_5l");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_5l = (TH1F*)fileWWG->Get("h_TotalEvents_5l");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_5l = (TH1F*)fileWZG->Get("h_TotalEvents_5l");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_5l = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_5l = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_5l = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_5l");

  std::cout << "5 l channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_5l1->GetBinContent(4)*LUMI2018 + hWZ_5l2->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_5l1->GetBinContent(4)*LUMI2018 + hTTZ_5l2->GetBinContent(4)*LUMI2018<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_5l->GetBinContent(4)*LUMI2018 + hZZ2_5l->GetBinContent(4)*LUMI2018 + hZZ3_5l->GetBinContent(4)*LUMI2018 + hZZ4_5l->GetBinContent(4)*LUMI2018 + hZZ5_5l->GetBinContent(4)*LUMI2018 + hZZ6_5l->GetBinContent(4)*LUMI2018 + hZZ7_5l->GetBinContent(4)*LUMI2018 + hZZ8_5l->GetBinContent(4)*LUMI2018 + hZZ9_5l->GetBinContent(4)*LUMI2018 + hZZ10_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_5l->GetBinContent(4)*LUMI2018 + hDY50_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM wzz"    << " & "   << hSignal_WZZ_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl; 
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_5l->GetBinContent(4)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the five leptons channel without additional cuts}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}

void getYieldsPerCutLevel_5l_ztag_wtag()
{

  total_noCut = total_Ztag = 0.0;
 
  TFile* fileTTW = new TFile("output_ttw_lnu_amcatnlo_1.root");
  TH1F *hTTW_5l = (TH1F*)fileTTW->Get("h_TotalEvents_5l");

  TFile* fileTTZ1 = new TFile("output_ttz_ll_mll1_amcatnlo_1.root");
  TH1F *hTTZ_5l1 = (TH1F*)fileTTZ1->Get("h_TotalEvents_5l");

  TFile* fileTTZ2 = new TFile("output_ttz_llvv_mll10_amcatnlo_1.root");
  TH1F *hTTZ_5l2 = (TH1F*)fileTTZ2->Get("h_TotalEvents_5l");

  TFile* fileTTbar1 = new TFile("output_ttbar_dilep_madgraph_1.root");
  TH1F *hTTbar1_5l = (TH1F*)fileTTbar1->Get("h_TotalEvents_5l");

  TFile* fileDY1050 = new TFile("output_dy_m1050_madgraph_1.root");
  TH1F *hDY1050_5l = (TH1F*)fileDY1050->Get("h_TotalEvents_5l");
  
  TFile* fileDY50 = new TFile("output_dy_m50_madgraph_1.root");
  TH1F *hDY50_5l = (TH1F*)fileDY50->Get("h_TotalEvents_5l");

  TFile* fileZZ1 = new TFile("output_zz_2l2v_powheg_1.root");
  TH1F *hZZ1_5l = (TH1F*)fileZZ1->Get("h_TotalEvents_5l");

  TFile* fileZZ2 = new TFile("output_zz_2l2q_powheg_1.root");
  TH1F *hZZ2_5l = (TH1F*)fileZZ2->Get("h_TotalEvents_5l");

  TFile* fileZZ3 = new TFile("output_zz_4l_powheg_1.root");
  TH1F *hZZ3_5l = (TH1F*)fileZZ3->Get("h_TotalEvents_5l");

  TFile* fileZZ4 = new TFile("output_ggzz_4e_mcfm_1.root");
  TH1F *hZZ4_5l = (TH1F*)fileZZ4->Get("h_TotalEvents_5l");

  TFile* fileZZ5 = new TFile("output_ggzz_4m_mcfm_1.root");
  TH1F *hZZ5_5l = (TH1F*)fileZZ5->Get("h_TotalEvents_5l");
 
  TFile* fileZZ6 = new TFile("output_ggzz_4t_mcfm_1.root");
  TH1F *hZZ6_5l = (TH1F*)fileZZ6->Get("h_TotalEvents_5l");

  TFile* fileZZ7 = new TFile("output_ggzz_2e2m_mcfm_1.root");
  TH1F *hZZ7_5l = (TH1F*)fileZZ7->Get("h_TotalEvents_5l");
  
  TFile* fileZZ8 = new TFile("output_ggzz_2e2t_mcfm_1.root");
  TH1F *hZZ8_5l = (TH1F*)fileZZ8->Get("h_TotalEvents_5l");

  TFile* fileZZ9 = new TFile("output_ggzz_2m2t_mcfm_1.root");
  TH1F *hZZ9_5l = (TH1F*)fileZZ9->Get("h_TotalEvents_5l");

  TFile* fileZZ10 = new TFile("output_ggh_hzz4l_powheg_1.root");
  TH1F *hZZ10_5l = (TH1F*)fileZZ10->Get("h_TotalEvents_5l");

  TFile* fileWZ1 = new TFile("output_wz_3lv_amcatnlo_1.root");
  TH1F *hWZ_5l1 = (TH1F*)fileWZ1->Get("h_TotalEvents_5l");

  TFile* fileWZ2 = new TFile("output_wz_2l2q_amcatnlo_1.root");
  TH1F *hWZ_5l2 = (TH1F*)fileWZ2->Get("h_TotalEvents_5l");

  TFile* fileWWG = new TFile("output_wwg_amcatnlo_1.root");
  TH1F *hWWG_5l = (TH1F*)fileWWG->Get("h_TotalEvents_5l");

  TFile* fileWZG = new TFile("output_wzg_amcatnlo_1.root");
  TH1F *hWZG_5l = (TH1F*)fileWZG->Get("h_TotalEvents_5l");

  TFile *fileSignal_WZZ = new TFile("output_wzz_amcatnlo_1.root");
  TH1F *hSignal_WZZ_5l = (TH1F*)fileSignal_WZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZZZ = new TFile("output_zzz_amcatnlo_1.root");
  TH1F *hSignal_ZZZ_5l = (TH1F*)fileSignal_ZZZ->Get("h_TotalEvents_5l");

  TFile *fileSignal_ZH_ZZ = new TFile("output_zh_zz_amcatnlo_1.root");
  TH1F *hSignal_ZH_ZZ_5l = (TH1F*)fileSignal_ZH_ZZ->Get("h_TotalEvents_5l");

  std::cout << "5 l channel" << std::endl;
  std::cout << "\\begin{table}[htb]" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\begin{tabular}{|| c | l ||}" << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Process  & Normalized yields  \\\\ [0.5ex]" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "$WZ$"              << " & "   << std::setprecision(3) << hWZ_5l1->GetBinContent(5)*LUMI2018 + hWZ_5l2->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hWZ_5l1->GetBinError(5)*LUMI2018*hWZ_5l1->GetBinError(5)*LUMI2018 + hWZ_5l2->GetBinError(5)*LUMI2018*hWZ_5l2->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}W$"      << " & "   << std::setprecision(2) << hTTW_5l->GetBinContent(5)*LUMI2018 <<  " $\\pm$ " << hTTW_5l->GetBinError(5)*LUMI2018 << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}Z$"      << " & "   << std::setprecision(3) << hTTZ_5l1->GetBinContent(5)*LUMI2018 + hTTZ_5l2->GetBinContent(5)*LUMI2018 <<  " $\\pm$ " << sqrt(hTTZ_5l1->GetBinError(5)*LUMI2018*hTTZ_5l1->GetBinError(5)*LUMI2018 + hTTZ_5l2->GetBinError(5)*LUMI2018*hTTZ_5l2->GetBinError(5)*LUMI2018)<< " \\\\ " << std::endl;
  std::cout << "$ZZ$"              << " & "   << std::setprecision(3) << hZZ1_5l->GetBinContent(5)*LUMI2018 + hZZ2_5l->GetBinContent(5)*LUMI2018 + hZZ3_5l->GetBinContent(5)*LUMI2018 + hZZ4_5l->GetBinContent(5)*LUMI2018 + hZZ5_5l->GetBinContent(5)*LUMI2018 + hZZ6_5l->GetBinContent(5)*LUMI2018 + hZZ7_5l->GetBinContent(5)*LUMI2018 + hZZ8_5l->GetBinContent(5)*LUMI2018 + hZZ9_5l->GetBinContent(5)*LUMI2018 + hZZ10_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hZZ1_5l->GetBinError(5)*LUMI2018*hZZ1_5l->GetBinError(5)*LUMI2018 + hZZ2_5l->GetBinError(5)*LUMI2018*hZZ2_5l->GetBinError(5)*LUMI2018 + hZZ3_5l->GetBinError(5)*LUMI2018*hZZ3_5l->GetBinError(5)*LUMI2018 + hZZ4_5l->GetBinError(5)*LUMI2018*hZZ4_5l->GetBinError(5)*LUMI2018 + hZZ5_5l->GetBinError(5)*LUMI2018*hZZ5_5l->GetBinError(5)*LUMI2018 + hZZ6_5l->GetBinError(5)*LUMI2018*hZZ6_5l->GetBinError(5)*LUMI2018 + hZZ7_5l->GetBinError(5)*LUMI2018*hZZ7_5l->GetBinError(5)*LUMI2018 + hZZ8_5l->GetBinError(5)*LUMI2018*hZZ8_5l->GetBinError(5)*LUMI2018 + hZZ9_5l->GetBinError(5)*LUMI2018*hZZ9_5l->GetBinError(5)*LUMI2018 + hZZ10_5l->GetBinError(5)*LUMI2018*hZZ10_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "$t\\bar{t}$"       << " & "   << std::setprecision(3) << hTTbar1_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hTTbar1_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "DY"                << " & "   << hDY1050_5l->GetBinContent(5)*LUMI2018 + hDY50_5l->GetBinContent(5)*LUMI2018 <<" $\\pm$ " << sqrt(hDY1050_5l->GetBinError(5)*LUMI2018*hDY1050_5l->GetBinError(5)*LUMI2018 + hDY50_5l->GetBinError(5)*LUMI2018*hDY50_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "ww$\\gamma$"                << " & "   << hWWG_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hWWG_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "wz$\\gamma$"                << " & "   << hWZG_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hWZG_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "Signal: SM wzz"    << " & "   << hSignal_WZZ_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hSignal_WZZ_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl; 
  std::cout << "Signal: SM zzz"    << " & "   << hSignal_ZZZ_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hSignal_ZZZ_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "Signal: zh $\\rightarrow$ zz"    << " & "   << hSignal_ZH_ZZ_5l->GetBinContent(5)*LUMI2018 << " $\\pm$ " << sqrt(hSignal_ZH_ZZ_5l->GetBinError(5)*LUMI2018) << " \\\\ " << std::endl;
  std::cout << "\\hline\\hline" << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\caption{Yield table for the five leptons channel with $ZZ$-tag and $W$-tag}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

}
