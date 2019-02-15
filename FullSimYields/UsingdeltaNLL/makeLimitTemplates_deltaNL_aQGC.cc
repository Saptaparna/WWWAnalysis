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

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

double LUMI = 35.9;
double BR = (0.33*0.33*0.67*3 + 0.33*0.33*0.33);

int makeLimitTemplates_deltaNL_aQGC(std::string outfile)
{

  TH1F *h_signal_neg2_5 = new TH1F("h_signal_neg2_5", "h_signal_neg2_5", 1.0, 0.0, 1.0); h_signal_neg2_5->Sumw2();
  TH1F *h_signal_neg2_4 = new TH1F("h_signal_neg2_4", "h_signal_neg2_4", 1.0, 0.0, 1.0); h_signal_neg2_4->Sumw2();
  TH1F *h_signal_neg2_3 = new TH1F("h_signal_neg2_3", "h_signal_neg2_3", 1.0, 0.0, 1.0); h_signal_neg2_3->Sumw2();
  TH1F *h_signal_neg2_2 = new TH1F("h_signal_neg2_2", "h_signal_neg2_2", 1.0, 0.0, 1.0); h_signal_neg2_2->Sumw2();
  TH1F *h_signal_neg2_1 = new TH1F("h_signal_neg2_1", "h_signal_neg2_1", 1.0, 0.0, 1.0); h_signal_neg2_1->Sumw2();
  TH1F *h_signal_neg2_0 = new TH1F("h_signal_neg2_0", "h_signal_neg2_0", 1.0, 0.0, 1.0); h_signal_neg2_0->Sumw2();
  TH1F *h_signal_neg1_9 = new TH1F("h_signal_neg1_9", "h_signal_neg1_9", 1.0, 0.0, 1.0); h_signal_neg1_9->Sumw2();
  TH1F *h_signal_neg1_8 = new TH1F("h_signal_neg1_8", "h_signal_neg1_8", 1.0, 0.0, 1.0); h_signal_neg1_8->Sumw2();
  TH1F *h_signal_neg1_7 = new TH1F("h_signal_neg1_7", "h_signal_neg1_7", 1.0, 0.0, 1.0); h_signal_neg1_7->Sumw2();
  TH1F *h_signal_neg1_6 = new TH1F("h_signal_neg1_6", "h_signal_neg1_6", 1.0, 0.0, 1.0); h_signal_neg1_6->Sumw2();
  TH1F *h_signal_neg1_5 = new TH1F("h_signal_neg1_5", "h_signal_neg1_5", 1.0, 0.0, 1.0); h_signal_neg1_5->Sumw2();
  TH1F *h_signal_neg1_4 = new TH1F("h_signal_neg1_4", "h_signal_neg1_4", 1.0, 0.0, 1.0); h_signal_neg1_4->Sumw2();
  TH1F *h_signal_neg1_3 = new TH1F("h_signal_neg1_3", "h_signal_neg1_3", 1.0, 0.0, 1.0); h_signal_neg1_3->Sumw2(); 
  TH1F *h_signal_neg1_2 = new TH1F("h_signal_neg1_2", "h_signal_neg1_2", 1.0, 0.0, 1.0); h_signal_neg1_2->Sumw2();
  TH1F *h_signal_neg1_1 = new TH1F("h_signal_neg1_1", "h_signal_neg1_1", 1.0, 0.0, 1.0); h_signal_neg1_1->Sumw2();
  TH1F *h_signal_neg1_0 = new TH1F("h_signal_neg1_0", "h_signal_neg1_0", 1.0, 0.0, 1.0); h_signal_neg1_0->Sumw2();
  TH1F *h_signal_neg0_9 = new TH1F("h_signal_neg0_9", "h_signal_neg0_9", 1.0, 0.0, 1.0); h_signal_neg0_9->Sumw2();
  TH1F *h_signal_neg0_8 = new TH1F("h_signal_neg0_8", "h_signal_neg0_8", 1.0, 0.0, 1.0); h_signal_neg0_8->Sumw2();
  TH1F *h_signal_neg0_7 = new TH1F("h_signal_neg0_7", "h_signal_neg0_7", 1.0, 0.0, 1.0); h_signal_neg0_7->Sumw2();
  TH1F *h_signal_neg0_6 = new TH1F("h_signal_neg0_6", "h_signal_neg0_6", 1.0, 0.0, 1.0); h_signal_neg0_6->Sumw2();
  TH1F *h_signal_neg0_5 = new TH1F("h_signal_neg0_5", "h_signal_neg0_5", 1.0, 0.0, 1.0); h_signal_neg0_5->Sumw2();
  TH1F *h_signal_neg0_4 = new TH1F("h_signal_neg0_4", "h_signal_neg0_4", 1.0, 0.0, 1.0); h_signal_neg0_4->Sumw2();
  TH1F *h_signal_neg0_3 = new TH1F("h_signal_neg0_3", "h_signal_neg0_3", 1.0, 0.0, 1.0); h_signal_neg0_3->Sumw2(); 
  TH1F *h_signal_neg0_2 = new TH1F("h_signal_neg0_2", "h_signal_neg0_2", 1.0, 0.0, 1.0); h_signal_neg0_2->Sumw2();
  TH1F *h_signal_neg0_1 = new TH1F("h_signal_neg0_1", "h_signal_neg0_1", 1.0, 0.0, 1.0); h_signal_neg0_1->Sumw2();
  TH1F *h_signal_0_0 = new TH1F("h_signal_0_0", "h_signal_0_0", 1.0, 0.0, 1.0); h_signal_0_0->Sumw2();  
  TH1F *h_signal_pos0_1 = new TH1F("h_signal_pos0_1", "h_signal_pos0_1", 1.0, 0.0, 1.0); h_signal_pos0_1->Sumw2();
  TH1F *h_signal_pos0_2 = new TH1F("h_signal_pos0_2", "h_signal_pos0_2", 1.0, 0.0, 1.0); h_signal_pos0_2->Sumw2();
  TH1F *h_signal_pos0_3 = new TH1F("h_signal_pos0_3", "h_signal_pos0_3", 1.0, 0.0, 1.0); h_signal_pos0_3->Sumw2();
  TH1F *h_signal_pos0_4 = new TH1F("h_signal_pos0_4", "h_signal_pos0_4", 1.0, 0.0, 1.0); h_signal_pos0_4->Sumw2();
  TH1F *h_signal_pos0_5 = new TH1F("h_signal_pos0_5", "h_signal_pos0_5", 1.0, 0.0, 1.0); h_signal_pos0_5->Sumw2();
  TH1F *h_signal_pos0_6 = new TH1F("h_signal_pos0_6", "h_signal_pos0_6", 1.0, 0.0, 1.0); h_signal_pos0_6->Sumw2();
  TH1F *h_signal_pos0_7 = new TH1F("h_signal_pos0_7", "h_signal_pos0_7", 1.0, 0.0, 1.0); h_signal_pos0_7->Sumw2();
  TH1F *h_signal_pos0_8 = new TH1F("h_signal_pos0_8", "h_signal_pos0_8", 1.0, 0.0, 1.0); h_signal_pos0_8->Sumw2();
  TH1F *h_signal_pos0_9 = new TH1F("h_signal_pos0_9", "h_signal_pos0_9", 1.0, 0.0, 1.0); h_signal_pos0_9->Sumw2();
  TH1F *h_signal_pos1_0 = new TH1F("h_signal_pos1_0", "h_signal_pos1_0", 1.0, 0.0, 1.0); h_signal_pos1_0->Sumw2();
  TH1F *h_signal_pos1_1 = new TH1F("h_signal_pos1_1", "h_signal_pos1_1", 1.0, 0.0, 1.0); h_signal_pos1_1->Sumw2();
  TH1F *h_signal_pos1_2 = new TH1F("h_signal_pos1_2", "h_signal_pos1_2", 1.0, 0.0, 1.0); h_signal_pos1_2->Sumw2();
  TH1F *h_signal_pos1_3 = new TH1F("h_signal_pos1_3", "h_signal_pos1_3", 1.0, 0.0, 1.0); h_signal_pos1_3->Sumw2();
  TH1F *h_signal_pos1_4 = new TH1F("h_signal_pos1_4", "h_signal_pos1_4", 1.0, 0.0, 1.0); h_signal_pos1_4->Sumw2();
  TH1F *h_signal_pos1_5 = new TH1F("h_signal_pos1_5", "h_signal_pos1_5", 1.0, 0.0, 1.0); h_signal_pos1_5->Sumw2();
  TH1F *h_signal_pos1_6 = new TH1F("h_signal_pos1_6", "h_signal_pos1_6", 1.0, 0.0, 1.0); h_signal_pos1_6->Sumw2();
  TH1F *h_signal_pos1_7 = new TH1F("h_signal_pos1_7", "h_signal_pos1_7", 1.0, 0.0, 1.0); h_signal_pos1_7->Sumw2();
  TH1F *h_signal_pos1_8 = new TH1F("h_signal_pos1_8", "h_signal_pos1_8", 1.0, 0.0, 1.0); h_signal_pos1_8->Sumw2();
  TH1F *h_signal_pos1_9 = new TH1F("h_signal_pos1_9", "h_signal_pos1_9", 1.0, 0.0, 1.0); h_signal_pos1_9->Sumw2();
  TH1F *h_signal_pos2_0 = new TH1F("h_signal_pos2_0", "h_signal_pos2_0", 1.0, 0.0, 1.0); h_signal_pos2_0->Sumw2();
  TH1F *h_signal_pos2_1 = new TH1F("h_signal_pos2_1", "h_signal_pos2_1", 1.0, 0.0, 1.0); h_signal_pos2_1->Sumw2();
  TH1F *h_signal_pos2_2 = new TH1F("h_signal_pos2_2", "h_signal_pos2_2", 1.0, 0.0, 1.0); h_signal_pos2_2->Sumw2();
  TH1F *h_signal_pos2_3 = new TH1F("h_signal_pos2_3", "h_signal_pos2_3", 1.0, 0.0, 1.0); h_signal_pos2_3->Sumw2();
  TH1F *h_signal_pos2_4 = new TH1F("h_signal_pos2_4", "h_signal_pos2_4", 1.0, 0.0, 1.0); h_signal_pos2_4->Sumw2();
  TH1F *h_signal_pos2_5 = new TH1F("h_signal_pos2_5", "h_signal_pos2_5", 1.0, 0.0, 1.0); h_signal_pos2_5->Sumw2();
  TH1F *h_bkg = new TH1F("h_bkg", "h_bkg", 1.0, 0.0, 1.0); h_bkg->Sumw2();

  //fill histos
  //h_signal_rgwt67->SetBinContent(1, );

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_signal_neg2_5->Write();
  h_signal_neg2_4->Write();
  h_signal_neg2_3->Write();
  h_signal_neg2_2->Write();
  h_signal_neg2_1->Write();
  h_signal_neg2_0->Write();
  h_signal_neg1_9->Write();
  h_signal_neg1_8->Write();
  h_signal_neg1_7->Write();
  h_signal_neg1_6->Write();
  h_signal_neg1_5->Write();
  h_signal_neg1_4->Write();
  h_signal_neg1_3->Write();
  h_signal_neg1_2->Write();
  h_signal_neg1_1->Write();
  h_signal_neg1_0->Write(); 
  h_signal_neg0_9->Write();
  h_signal_neg0_8->Write();
  h_signal_neg0_7->Write();
  h_signal_neg0_6->Write();
  h_signal_neg0_5->Write();
  h_signal_neg0_4->Write();
  h_signal_neg0_3->Write();
  h_signal_neg0_2->Write();
  h_signal_neg0_1->Write();
  h_signal_0_0->Write();
  h_signal_pos0_1->Write();
  h_signal_pos0_2->Write();
  h_signal_pos0_3->Write();
  h_signal_pos0_4->Write();
  h_signal_pos0_5->Write();
  h_signal_pos0_6->Write();
  h_signal_pos0_7->Write();
  h_signal_pos0_8->Write();
  h_signal_pos0_9->Write();
  h_signal_pos1_0->Write();
  h_signal_pos1_1->Write();
  h_signal_pos1_2->Write();
  h_signal_pos1_3->Write();
  h_signal_pos1_4->Write();
  h_signal_pos1_5->Write();
  h_signal_pos1_6->Write();
  h_signal_pos1_7->Write();
  h_signal_pos1_8->Write();
  h_signal_pos1_9->Write();
  h_signal_pos2_0->Write();
  h_signal_pos2_1->Write();
  h_signal_pos2_2->Write();
  h_signal_pos2_3->Write();
  h_signal_pos2_4->Write();
  h_signal_pos2_5->Write();
  h_bkg->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
}
