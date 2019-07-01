#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TRandom.h>

using std::string;

double MUON_MASS = 105.6583745*10e-03;
double ELECTRON_MASS = 511*10e-06;

bool sameVal(double a, double b)
{
   return fabs(a - b) < 20.0;
}

bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double mass)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, mass);
  return object_p4;
}
/*
typedef struct
{
  float pt; 
  float phi;
  float eta;
  float pfiso;
  float trkiso;
  int charge;
  float dz;
  float d0;
  float sip3d;
  bool id;
  float recoEW;
  float triggerEW;
  float trigger;
  bool isTightUCSD;
  bool isTightMIT;
  bool isveryLooseUCSD;
  bool isHLTsafeMIT;
} leptonInfo;
*/

typedef struct
{
  float pt; 
  float phi;
  float eta;
  float sip3d;
  float dz;
  float dxy;
  float mva;
  TLorentzVector lep_lv;
  int id; //charge information
} leptonInfo;

typedef struct
{
  float jetPt;
  float jetEta;
  float jetPhi;
  float jetMass;
  float jetCSV;
} jetInfo;

typedef struct
{
  float ak8JetPrunmass;
  float ak8JetTrimmass;
  float ak8Jetsd0;
  float ak8JetPt;
  float ak8JetEta;
  float ak8JetPhi;
  float ak8JetTau1;
  float ak8JetTau2;
} fatJetInfo;

typedef struct
{
    TLorentzVector JetLV;
    float BTag_CSV;
} AnalysisJetInfo;

bool sortLeptonsInDescendingpT(leptonInfo lep1, leptonInfo lep2)
{
  return (lep1.pt > lep2.pt);
}

bool sortJetVectorsInDescendingpT(AnalysisJetInfo jet1, AnalysisJetInfo jet2)
{
  return (jet1.JetLV.Pt() > jet2.JetLV.Pt());
}

typedef struct
{
  TH1F *h_mu_pt_leading;
  TH1F *h_mu_pt_second;
  TH1F *h_mu_pt_third;
  TH1F *h_mu_pt_fourth;
  TH1F *h_mu_pt_fifth;
  TH1F *h_l_pt_leading;
  TH1F *h_l_pt_second;
  TH1F *h_l_pt_third;
  TH1F *h_l_pt_fourth;
  TH1F *h_l_pt_fifth;
  TH1F *h_mu_eta_leading;
  TH1F *h_mu_eta_second;
  TH1F *h_mu_eta_third;
  TH1F *h_mu_eta_fourth;
  TH1F *h_mu_eta_fifth;
  TH1F *h_l_eta_leading;
  TH1F *h_l_eta_second;
  TH1F *h_l_eta_third;
  TH1F *h_l_eta_fourth;
  TH1F *h_l_eta_fifth;
  TH1F *h_mu_phi_leading;
  TH1F *h_mu_phi_second;
  TH1F *h_mu_phi_third;
  TH1F *h_mu_phi_fourth;
  TH1F *h_mu_phi_fifth;
  TH1F *h_l_phi_leading;
  TH1F *h_l_phi_second;
  TH1F *h_l_phi_third;
  TH1F *h_l_phi_fourth;
  TH1F *h_l_phi_fifth;
  TH1F *h_mu_pTRatio_leading;
  TH1F *h_mu_pTRatio_second;
  TH1F *h_mu_ip3d_leading;
  TH1F *h_mu_ip3d_second;
  TH1F *h_mu_ip3d_third;
  TH1F *h_mu_ip3d_fourth;
  TH1F *h_mu_ip3d_fifth;
  TH1F *h_l_ip3d_leading;
  TH1F *h_l_ip3d_second;
  TH1F *h_l_ip3d_third;
  TH1F *h_l_ip3d_fourth;
  TH1F *h_l_ip3d_fifth;
  TH1F *h_mu_d0_leading;
  TH1F *h_mu_d0_second;
  TH1F *h_mu_d0_third;
  TH1F *h_mu_d0_fourth;
  TH1F *h_mu_d0_fifth;
  TH1F *h_l_d0_leading;
  TH1F *h_l_d0_second;
  TH1F *h_l_d0_third;
  TH1F *h_l_d0_fourth;
  TH1F *h_l_d0_fifth;
  TH1F *h_mu_dz_leading;
  TH1F *h_mu_dz_second;
  TH1F *h_mu_dz_third;
  TH1F *h_mu_dz_fourth;
  TH1F *h_mu_dz_fifth;
  TH1F *h_l_dz_leading;
  TH1F *h_l_dz_second;
  TH1F *h_l_dz_third;
  TH1F *h_l_dz_fourth;
  TH1F *h_l_dz_fifth;
  TH1F *h_mu_pTrel_leading;
  TH1F *h_mu_pTrel_second;
  TH1F *h_mu_relIso_leading;
  TH1F *h_mu_relIso_second;
  TH1F *h_mu_deltaR_leading;
  TH1F *h_mu_deltaR_second;
  TH1F *h_InvariantMass;
  TH1F *h_MET;
  TH1F *h_nJets;
  TH1F *h_nbJets;
  TH1F *h_InvariantMassJJ;
  TH1F *h_jet_pt_leading;
  TH1F *h_jet_pt_second;
  TH1F *h_jet_eta_leading;
  TH1F *h_jet_eta_second;
  TH1F *h_jet_phi_leading;
  TH1F *h_jet_phi_second;
  TH1F *h_jet_csv_leading;
  TH1F *h_jet_csv_second;
  TH1F *h_bjet_pt_leading;
  TH1F *h_bjet_eta_leading;
  TH1F *h_bjet_phi_leading;
  TH1F *h_bjet_pt_second;
  TH1F *h_bjet_eta_second;
  TH1F *h_bjet_phi_second;
  TH1F *h_el_pt_leading;
  TH1F *h_el_pt_second;
  TH1F *h_el_pt_third;
  TH1F *h_el_pt_fourth;
  TH1F *h_el_pt_fifth;
  TH1F *h_el_eta_leading;
  TH1F *h_el_eta_second;
  TH1F *h_el_eta_third;
  TH1F *h_el_eta_fourth;
  TH1F *h_el_eta_fifth;
  TH1F *h_el_phi_leading;
  TH1F *h_el_phi_second;
  TH1F *h_el_phi_third;
  TH1F *h_el_phi_fourth;
  TH1F *h_el_phi_fifth;
  TH1F *h_el_pTRatio_leading;
  TH1F *h_el_pTRatio_second;
  TH1F *h_el_ip3d_leading;
  TH1F *h_el_ip3d_second;
  TH1F *h_el_ip3d_third;
  TH1F *h_el_ip3d_fourth;
  TH1F *h_el_ip3d_fifth;
  TH1F *h_el_d0_leading;
  TH1F *h_el_d0_second;
  TH1F *h_el_d0_third;
  TH1F *h_el_d0_fourth;
  TH1F *h_el_d0_fifth;
  TH1F *h_el_dz_leading;
  TH1F *h_el_dz_second;
  TH1F *h_el_dz_third;
  TH1F *h_el_dz_fourth;
  TH1F *h_el_dz_fifth;
  TH1F *h_el_pTrel_leading;
  TH1F *h_el_pTrel_second;
  TH1F *h_el_relIso_leading;
  TH1F *h_el_relIso_second;
  TH1F *h_el_deltaR_leading;
  TH1F *h_el_deltaR_second;
  TH1F *h_nfatJets;
  TH1F *h_fatjet_pt_leading;
  TH1F *h_fatjet_eta_leading;
  TH1F *h_fatjet_phi_leading;
  TH1F *h_fatjet_prunMass_leading;
  TH1F *h_fatjet_trimMass_leading;
  TH1F *h_fatjet_sd0Mass_leading;
  TH1F *h_fatjet_tau21_leading;
  TH1F *h_M4Mu;
  TH1F *h_M4El;
  TH1F *h_M12Mu;
  TH1F *h_M23Mu;
  TH1F *h_M34Mu;
  TH1F *h_M41Mu;
  TH1F *h_M12ElEl;;
  TH1F *h_M12MuMu;
} HistCollection;

void initializeHistCollection(HistCollection &histCol, std::string suffix)
{
   histCol.h_mu_pt_leading = new TH1F(("h_mu_pt_leading_"+suffix).c_str(), "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_leading->Sumw2();
   histCol.h_mu_pt_second = new TH1F(("h_mu_pt_second_"+suffix).c_str(), "Second muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_second->Sumw2();
   histCol.h_mu_pt_third = new TH1F(("h_mu_pt_third_"+suffix).c_str(), "Third muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_third->Sumw2();
   histCol.h_mu_pt_fourth = new TH1F(("h_mu_pt_fourth_"+suffix).c_str(), "Fourth muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_fourth->Sumw2();
   histCol.h_mu_pt_fifth = new TH1F(("h_mu_pt_fifth_"+suffix).c_str(), "Fifth muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_fifth->Sumw2();
   histCol.h_mu_eta_leading = new TH1F(("h_mu_eta_leading_"+suffix).c_str(), "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_leading->Sumw2();
   histCol.h_mu_eta_second = new TH1F(("h_mu_eta_second_"+suffix).c_str(), "Second muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_second->Sumw2();
   histCol.h_mu_eta_third = new TH1F(("h_mu_eta_third_"+suffix).c_str(), "Third muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_third->Sumw2();
   histCol.h_mu_eta_fourth = new TH1F(("h_mu_eta_fourth_"+suffix).c_str(), "Fourth muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_fourth->Sumw2();
   histCol.h_mu_eta_fifth = new TH1F(("h_mu_eta_fifth_"+suffix).c_str(), "Fifth muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_fifth->Sumw2();
   histCol.h_mu_phi_leading = new TH1F(("h_mu_phi_leading_"+suffix).c_str(), "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_leading->Sumw2();
   histCol.h_mu_phi_second = new TH1F(("h_mu_phi_second_"+suffix).c_str(), "Second muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_second->Sumw2();
   histCol.h_mu_phi_third = new TH1F(("h_mu_phi_third_"+suffix).c_str(), "Third muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_third->Sumw2();
   histCol.h_mu_phi_fourth = new TH1F(("h_mu_phi_fourth_"+suffix).c_str(), "Fourth muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_fourth->Sumw2();
   histCol.h_mu_phi_fifth = new TH1F(("h_mu_phi_fifth_"+suffix).c_str(), "Fifth muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_fifth->Sumw2();
   histCol.h_l_pt_leading = new TH1F(("h_l_pt_leading_"+suffix).c_str(), "Leading lepton pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_l_pt_leading->Sumw2();
   histCol.h_l_pt_second = new TH1F(("h_l_pt_second_"+suffix).c_str(), "Second lepton pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_l_pt_second->Sumw2();
   histCol.h_l_pt_third = new TH1F(("h_l_pt_third_"+suffix).c_str(), "Third lepton pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_l_pt_third->Sumw2();
   histCol.h_l_pt_fourth = new TH1F(("h_l_pt_fourth_"+suffix).c_str(), "Fourth lepton pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_l_pt_fourth->Sumw2();
   histCol.h_l_pt_fifth = new TH1F(("h_l_pt_fifth_"+suffix).c_str(), "Fifth lepton pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_l_pt_fifth->Sumw2();
   histCol.h_l_eta_leading = new TH1F(("h_l_eta_leading_"+suffix).c_str(), "Leading lepton #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_l_eta_leading->Sumw2();
   histCol.h_l_eta_second = new TH1F(("h_l_eta_second_"+suffix).c_str(), "Second lepton #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_l_eta_second->Sumw2();
   histCol.h_l_eta_third = new TH1F(("h_l_eta_third_"+suffix).c_str(), "Third lepton #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_l_eta_third->Sumw2();
   histCol.h_l_eta_fourth = new TH1F(("h_l_eta_fourth_"+suffix).c_str(), "Fourth lepton #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_l_eta_fourth->Sumw2();
   histCol.h_l_eta_fifth = new TH1F(("h_l_eta_fifth_"+suffix).c_str(), "Fifth lepton #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_l_eta_fifth->Sumw2();
   histCol.h_l_phi_leading = new TH1F(("h_l_phi_leading_"+suffix).c_str(), "Leading lepton #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_l_phi_leading->Sumw2();
   histCol.h_l_phi_second = new TH1F(("h_l_phi_second_"+suffix).c_str(), "Second lepton #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_l_phi_second->Sumw2();
   histCol.h_l_phi_third = new TH1F(("h_l_phi_third_"+suffix).c_str(), "Third lepton #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_l_phi_third->Sumw2();
   histCol.h_l_phi_fourth = new TH1F(("h_l_phi_fourth_"+suffix).c_str(), "Fourth lepton #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_l_phi_fourth->Sumw2();
   histCol.h_l_phi_fifth = new TH1F(("h_l_phi_fifth_"+suffix).c_str(), "Fifth lepton #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_l_phi_fifth->Sumw2();

   histCol.h_mu_pTRatio_leading = new TH1F(("h_mu_pTRatio_leading_"+suffix).c_str(), "Leading muon pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_mu_pTRatio_leading->Sumw2();
   histCol.h_mu_pTRatio_second = new TH1F(("h_mu_pTRatio_second_"+suffix).c_str(), "Trailing muon pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_mu_pTRatio_second->Sumw2();
   histCol.h_mu_ip3d_leading = new TH1F(("h_mu_ip3d_leading_"+suffix).c_str(), "Leading muon 3D impact parameter; Leading muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_leading->Sumw2();
   histCol.h_mu_ip3d_second = new TH1F(("h_mu_ip3d_second_"+suffix).c_str(), "Second muon 3D impact parameter; Trailing muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_second->Sumw2();
   histCol.h_mu_ip3d_third = new TH1F(("h_mu_ip3d_third_"+suffix).c_str(), "Third muon 3D impact parameter; Third muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_third->Sumw2();
   histCol.h_mu_ip3d_fourth = new TH1F(("h_mu_ip3d_fourth_"+suffix).c_str(), "Fourth muon 3D impact parameter; Fourth muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_fourth->Sumw2();
   histCol.h_mu_ip3d_fifth = new TH1F(("h_mu_ip3d_fifth_"+suffix).c_str(), "Fifth muon 3D impact parameter; Fifth muon 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_ip3d_fifth->Sumw2();

   histCol.h_l_ip3d_leading = new TH1F(("h_l_ip3d_leading_"+suffix).c_str(), "Leading lepton 3D impact parameter; Leading lepton 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_ip3d_leading->Sumw2();   
   histCol.h_l_ip3d_second = new TH1F(("h_l_ip3d_second_"+suffix).c_str(), "Second lepton 3D impact parameter; Trailing lepton 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_ip3d_second->Sumw2();
   histCol.h_l_ip3d_third = new TH1F(("h_l_ip3d_third_"+suffix).c_str(), "Third lepton 3D impact parameter; Third lepton 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_ip3d_third->Sumw2();
   histCol.h_l_ip3d_fourth = new TH1F(("h_l_ip3d_fourth_"+suffix).c_str(), "Fourth lepton 3D impact parameter; Fourth lepton 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_ip3d_fourth->Sumw2();
   histCol.h_l_ip3d_fifth = new TH1F(("h_l_ip3d_fifth_"+suffix).c_str(), "Fifth lepton 3D impact parameter; Fifth lepton 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_ip3d_fifth->Sumw2();

   histCol.h_mu_d0_leading = new TH1F(("h_mu_d0_leading_"+suffix).c_str(), "Leading muon d0 impact parameter; Leading muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_leading->Sumw2();
   histCol.h_mu_d0_second = new TH1F(("h_mu_d0_second_"+suffix).c_str(), "Second muon d0 impact parameter; Second muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_second->Sumw2();
   histCol.h_mu_d0_third = new TH1F(("h_mu_d0_third_"+suffix).c_str(), "Third muon d0 impact parameter; Third muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_third->Sumw2();
   histCol.h_mu_d0_fourth = new TH1F(("h_mu_d0_fourth_"+suffix).c_str(), "Fourth muon d0 impact parameter; Fourth muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_fourth->Sumw2();
   histCol.h_mu_d0_fifth = new TH1F(("h_mu_d0_fifth_"+suffix).c_str(), "Fifth muon d0 impact parameter; Fifth muon d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_d0_fifth->Sumw2();

   histCol.h_l_d0_leading = new TH1F(("h_l_d0_leading_"+suffix).c_str(), "Leading lepton d0 impact parameter; Leading lepton d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_d0_leading->Sumw2();
   histCol.h_l_d0_second = new TH1F(("h_l_d0_second_"+suffix).c_str(), "Second lepton d0 impact parameter; Second lepton d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_d0_second->Sumw2();
   histCol.h_l_d0_third = new TH1F(("h_l_d0_third_"+suffix).c_str(), "Third lepton d0 impact parameter; Third lepton d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_d0_third->Sumw2();
   histCol.h_l_d0_fourth = new TH1F(("h_l_d0_fourth_"+suffix).c_str(), "Fourth lepton d0 impact parameter; Fourth lepton d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_d0_fourth->Sumw2();
   histCol.h_l_d0_fifth = new TH1F(("h_l_d0_fifth_"+suffix).c_str(), "Fifth lepton d0 impact parameter; Fifth lepton d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_d0_fifth->Sumw2();
   
   histCol.h_mu_dz_leading = new TH1F(("h_mu_dz_leading_"+suffix).c_str(), "Leading muon dz impact parameter; Leading muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_leading->Sumw2();
   histCol.h_mu_dz_second = new TH1F(("h_mu_dz_second_"+suffix).c_str(), "Second muon dz impact parameter; Second muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_second->Sumw2();
   histCol.h_mu_dz_third = new TH1F(("h_mu_dz_third_"+suffix).c_str(), "Third muon dz impact parameter; Third muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_third->Sumw2();
   histCol.h_mu_dz_fourth = new TH1F(("h_mu_dz_fourth_"+suffix).c_str(), "Fourth muon dz impact parameter; Fourth muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_fourth->Sumw2();
   histCol.h_mu_dz_fifth = new TH1F(("h_mu_dz_fifth_"+suffix).c_str(), "Fifth muon dz impact parameter; Fifth muon dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_mu_dz_fifth->Sumw2();

   histCol.h_l_dz_leading = new TH1F(("h_l_dz_leading_"+suffix).c_str(), "Leading lepton dz impact parameter; Leading lepton dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_dz_leading->Sumw2();
   histCol.h_l_dz_second = new TH1F(("h_l_dz_second_"+suffix).c_str(), "Second lepton dz impact parameter; Second lepton dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_dz_second->Sumw2();
   histCol.h_l_dz_third = new TH1F(("h_l_dz_third_"+suffix).c_str(), "Third lepton dz impact parameter; Third lepton dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_dz_third->Sumw2();
   histCol.h_l_dz_fourth = new TH1F(("h_l_dz_fourth_"+suffix).c_str(), "Fourth lepton dz impact parameter; Fourth lepton dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_dz_fourth->Sumw2();
   histCol.h_l_dz_fifth = new TH1F(("h_l_dz_fifth_"+suffix).c_str(), "Fifth lepton dz impact parameter; Fifth lepton dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_l_dz_fifth->Sumw2();

   histCol.h_mu_pTrel_leading = new TH1F(("h_mu_pTrel_leading_"+suffix).c_str(), "Leading muon pT_{rel}; Leading muon pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_mu_pTrel_leading->Sumw2();
   histCol.h_mu_pTrel_second = new TH1F(("h_mu_pTrel_second_"+suffix).c_str(), "Trailing muon pT_{rel}; Trailing muon pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_mu_pTrel_second->Sumw2();
   histCol.h_mu_relIso_leading = new TH1F(("h_mu_relIso_leading_"+suffix).c_str(), "Leading relative isolation; Leading relative isolation; Events", 10000, 0.0, 1.0);histCol.h_mu_relIso_leading->Sumw2();
   histCol.h_mu_relIso_second = new TH1F(("h_mu_relIso_second_"+suffix).c_str(), "Trailing relative isolation; Trailing relative isolation; Events", 10000, 0.0, 1.0);histCol.h_mu_relIso_second->Sumw2();
   histCol.h_mu_deltaR_leading = new TH1F(("h_mu_deltaR_leading_"+suffix).c_str(), "Leading muon-jet #DeltaR; Leading muon-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_mu_deltaR_leading->Sumw2();
   histCol.h_mu_deltaR_second = new TH1F(("h_mu_deltaR_second_"+suffix).c_str(), "Trailing muon-jet #DeltaR; Trailing muon-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_mu_deltaR_second->Sumw2();
   histCol.h_InvariantMass=new TH1F(("h_InvariantMass_"+suffix).c_str(), "Di-lepton invariant mass; M_{ll} [GeV]; Events/GeV", 10000, 0.0, 1000.0); histCol.h_InvariantMass->Sumw2();
   histCol.h_MET=new TH1F(("h_MET_"+suffix).c_str(), "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); histCol.h_MET->Sumw2();
   histCol.h_nJets = new TH1F(("h_nJets_"+suffix).c_str(), "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);histCol.h_nJets->Sumw2();
   histCol.h_nbJets = new TH1F(("h_nbJets_"+suffix).c_str(), "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);histCol.h_nbJets->Sumw2();
   histCol.h_InvariantMassJJ = new TH1F(("h_InvariantMassJJ_"+suffix).c_str(), "Di-jet invariant mass M_{jj} [GeV]; Events/GeV", 10000.0, 0.0, 1000.0);histCol.h_InvariantMassJJ->Sumw2();
   histCol.h_jet_pt_leading=new TH1F(("h_jet_pt_leading_"+suffix).c_str(), "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_leading->Sumw2();
   histCol.h_jet_pt_second=new TH1F(("h_jet_pt_second_"+suffix).c_str(), "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_second->Sumw2();
   histCol.h_jet_eta_leading=new TH1F(("h_jet_eta_leading_"+suffix).c_str(), "Leading jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_leading->Sumw2();
   histCol.h_jet_eta_second=new TH1F(("h_jet_eta_second_"+suffix).c_str(), "Trailing jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_second->Sumw2();
   histCol.h_jet_phi_leading=new TH1F(("h_jet_phi_leading_"+suffix).c_str(), "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_leading->Sumw2();
   histCol.h_jet_phi_second=new TH1F(("h_jet_phi_second_"+suffix).c_str(), "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_second->Sumw2();
   histCol.h_jet_csv_leading=new TH1F(("h_jet_csv_leading_"+suffix).c_str(), "Leading jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_leading->Sumw2();
   histCol.h_jet_csv_second=new TH1F(("h_jet_csv_second_"+suffix).c_str(), "Trailing jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_second->Sumw2();
   histCol.h_bjet_pt_leading=new TH1F(("h_bjet_pt_leading_"+suffix).c_str(), "Leading b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_leading->Sumw2();
   histCol.h_bjet_eta_leading=new TH1F(("h_bjet_eta_leading_"+suffix).c_str(), "Leading b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_leading->Sumw2();
   histCol.h_bjet_phi_leading=new TH1F(("h_bjet_phi_leading_"+suffix).c_str(), "Leading b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_leading->Sumw2();
   histCol.h_bjet_pt_second=new TH1F(("h_bjet_pt_second_"+suffix).c_str(), "Trailing b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_second->Sumw2();
   histCol.h_bjet_eta_second=new TH1F(("h_bjet_eta_second_"+suffix).c_str(), "Trailing b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_second->Sumw2();
   histCol.h_bjet_phi_second=new TH1F(("h_bjet_phi_second_"+suffix).c_str(), "Trailing b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_second->Sumw2();
   histCol.h_el_pt_leading = new TH1F(("h_el_pt_leading_"+suffix).c_str(), "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_leading->Sumw2();     
   histCol.h_el_pt_second = new TH1F(("h_el_pt_second_"+suffix).c_str(), "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_second->Sumw2();
   histCol.h_el_pt_third = new TH1F(("h_el_pt_third_"+suffix).c_str(), "Third electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_third->Sumw2();
   histCol.h_el_pt_fourth = new TH1F(("h_el_pt_fourth_"+suffix).c_str(), "Fourth electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_fourth->Sumw2();
   histCol.h_el_pt_fifth = new TH1F(("h_el_pt_fifth_"+suffix).c_str(), "Fifth electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_fifth->Sumw2();
   histCol.h_el_eta_leading = new TH1F(("h_el_eta_leading_"+suffix).c_str(), "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_leading->Sumw2();
   histCol.h_el_eta_second = new TH1F(("h_el_eta_second_"+suffix).c_str(), "Trailing electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_second->Sumw2();
   histCol.h_el_eta_third = new TH1F(("h_el_eta_third_"+suffix).c_str(), "Third electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_third->Sumw2();
   histCol.h_el_eta_fourth = new TH1F(("h_el_eta_fourth_"+suffix).c_str(), "Fourth electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_fourth->Sumw2();
   histCol.h_el_eta_fifth = new TH1F(("h_el_eta_fifth_"+suffix).c_str(), "Fifth electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_fifth->Sumw2();
   histCol.h_el_phi_leading = new TH1F(("h_el_phi_leading_"+suffix).c_str(), "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_leading->Sumw2();
   histCol.h_el_phi_second = new TH1F(("h_el_phi_second_"+suffix).c_str(), "Trailing electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_second->Sumw2();
   histCol.h_el_phi_third = new TH1F(("h_el_phi_third_"+suffix).c_str(), "Third electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_third->Sumw2();
   histCol.h_el_phi_fourth = new TH1F(("h_el_phi_fourth_"+suffix).c_str(), "Fourth electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_fourth->Sumw2();
   histCol.h_el_phi_fifth = new TH1F(("h_el_phi_fifth_"+suffix).c_str(), "Fifth electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_fifth->Sumw2();
   histCol.h_el_pTRatio_leading = new TH1F(("h_el_pTRatio_leading_"+suffix).c_str(), "Leading electron pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_el_pTRatio_leading->Sumw2();
   histCol.h_el_pTRatio_second = new TH1F(("h_el_pTRatio_second_"+suffix).c_str(), "Trailing electron pT ratio ; pT ratio; Events", 150, 0.0, 1.5); histCol.h_el_pTRatio_second->Sumw2();
   histCol.h_el_ip3d_leading = new TH1F(("h_el_ip3d_leading_"+suffix).c_str(), "Leading electron 3D impact parameter; Leading electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_leading->Sumw2();
   histCol.h_el_ip3d_second = new TH1F(("h_el_ip3d_second_"+suffix).c_str(), "Trailing electron 3D impact parameter; Trailing electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_second->Sumw2();
   histCol.h_el_ip3d_third = new TH1F(("h_el_ip3d_third_"+suffix).c_str(), "Third electron 3D impact parameter; Third electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_third->Sumw2();
   histCol.h_el_ip3d_fourth = new TH1F(("h_el_ip3d_fourth_"+suffix).c_str(), "Fourth electron 3D impact parameter; Fourth electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_fourth->Sumw2();
   histCol.h_el_ip3d_fifth = new TH1F(("h_el_ip3d_fifth_"+suffix).c_str(), "Fifth electron 3D impact parameter; Fifth electron 3D impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_ip3d_fifth->Sumw2();
   histCol.h_el_d0_leading = new TH1F(("h_el_d0_leading_"+suffix).c_str(), "Leading electron d0 impact parameter; Leading electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_leading->Sumw2();
   histCol.h_el_d0_second = new TH1F(("h_el_d0_second_"+suffix).c_str(), "Trailing electron d0 impact parameter; Trailing electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_second->Sumw2();  
   histCol.h_el_d0_third = new TH1F(("h_el_d0_third_"+suffix).c_str(), "Third electron d0 impact parameter; Third electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_third->Sumw2();
   histCol.h_el_d0_fourth = new TH1F(("h_el_d0_fourth_"+suffix).c_str(), "Fourth electron d0 impact parameter; Fourth electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_fourth->Sumw2();
   histCol.h_el_d0_fifth = new TH1F(("h_el_d0_fifth_"+suffix).c_str(), "Fifth electron d0 impact parameter; Fifth electron d0 impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_d0_fifth->Sumw2();
   histCol.h_el_dz_leading = new TH1F(("h_el_dz_leading_"+suffix).c_str(), "Leading electron dz impact parameter; Leading electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_leading->Sumw2();
   histCol.h_el_dz_second = new TH1F(("h_el_dz_second_"+suffix).c_str(), "Trailing electron dz impact parameter; Trailing electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_second->Sumw2();
   histCol.h_el_dz_third = new TH1F(("h_el_dz_third_"+suffix).c_str(), "Third electron dz impact parameter; Third electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_third->Sumw2();
   histCol.h_el_dz_fourth = new TH1F(("h_el_dz_fourth_"+suffix).c_str(), "Fourth electron dz impact parameter; Fourth electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_fourth->Sumw2();
   histCol.h_el_dz_fifth = new TH1F(("h_el_dz_fifth_"+suffix).c_str(), "Fifth electron dz impact parameter; Fifth electron dz impact parameter; Events", 10000, 0.0, 1.0);histCol.h_el_dz_fifth->Sumw2();
   histCol.h_el_pTrel_leading = new TH1F(("h_el_pTrel_leading_"+suffix).c_str(), "Leading electron pT_{rel}; Leading electron pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_el_pTrel_leading->Sumw2();
   histCol.h_el_pTrel_second = new TH1F(("h_el_pTrel_second_"+suffix).c_str(), "Trailing electron pT_{rel}; Trailing electron pT_{rel} [GeV]; Events", 100000, 0.0, 1000.0);histCol.h_el_pTrel_second->Sumw2();
   histCol.h_el_relIso_leading = new TH1F(("h_el_relIso_leading_"+suffix).c_str(), "Leading relative isolation; Leading relative isolation; Events", 10000, 0.0, 1.0);histCol.h_el_relIso_leading->Sumw2();
   histCol.h_el_relIso_second = new TH1F(("h_el_relIso_second_"+suffix).c_str(), "Trailing relative isolation; Trailing relative isolation; Events", 10000, 0.0, 1.0);histCol.h_el_relIso_second->Sumw2();
   histCol.h_el_deltaR_leading = new TH1F(("h_el_deltaR_leading_"+suffix).c_str(), "Leading electron-jet #DeltaR; Leading electron-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_el_deltaR_leading->Sumw2();
   histCol.h_el_deltaR_second = new TH1F(("h_el_deltaR_second_"+suffix).c_str(), "Trailing electron-jet #DeltaR; Trailing electron-jet #DeltaR; Events", 10000, 0.0, 10.0);histCol.h_el_deltaR_second->Sumw2();
   histCol.h_fatjet_pt_leading=new TH1F(("h_fatjet_pt_leading_"+suffix).c_str(), "Leading fatjet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_fatjet_pt_leading->Sumw2();
   histCol.h_fatjet_phi_leading=new TH1F(("h_fatjet_phi_leading_"+suffix).c_str(), "Leading fatjet #phi; #phi; Events", 800, -4.0, 4.0); histCol.h_fatjet_phi_leading->Sumw2();
   histCol.h_fatjet_eta_leading=new TH1F(("h_fatjet_eta_leading_"+suffix).c_str(), "Leading fatjet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_fatjet_eta_leading->Sumw2();
   histCol.h_nfatJets = new TH1F(("h_nfatJets_"+suffix).c_str(), "Number of fatJets; Number of fatJets; Events", 20, -0.5, 19.5);histCol.h_nfatJets->Sumw2();
   histCol.h_fatjet_prunMass_leading = new TH1F(("h_fatjet_prunMass_leading_"+suffix).c_str(), "Pruned Mass [GeV]; Pruned Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_prunMass_leading->Sumw2();
   histCol.h_fatjet_trimMass_leading = new TH1F(("h_fatjet_trimMass_leading_"+suffix).c_str(), "Trimmed Mass [GeV]; Trimmed Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_trimMass_leading->Sumw2();
   histCol.h_fatjet_sd0Mass_leading = new TH1F(("h_fatjet_sd0Mass_leading_"+suffix).c_str(), "Soft drop Mass [GeV]; soft drop Mass [GeV]; Events", 10000, 0, 1000);histCol.h_fatjet_sd0Mass_leading->Sumw2();
   histCol.h_fatjet_tau21_leading = new TH1F(("h_fatjet_tau21_leading_"+suffix).c_str(), "N-subjettiness; #tau2/#tau1; Events", 1000, 0, 1.0);histCol.h_fatjet_tau21_leading->Sumw2();
   histCol.h_M4Mu = new TH1F(("h_M4Mu_"+suffix).c_str(), "4 muon invariant mass; mass_{4#mu}; Events", 1000, 0, 1000.0); histCol.h_M4Mu->Sumw2();
   histCol.h_M4El = new TH1F(("h_M4El_"+suffix).c_str(), "4 electron invariant mass; mass_{4e}; Events", 1000, 0, 1000.0); histCol.h_M4El->Sumw2(); 
   histCol.h_M12Mu = new TH1F(("h_M12Mu_"+suffix).c_str(), "2 muon invariant mass; mass_{12}; Events", 1000, 0, 1000.0); histCol.h_M12Mu->Sumw2();
   histCol.h_M23Mu = new TH1F(("h_M23Mu_"+suffix).c_str(), "2 muon invariant mass; mass_{23}; Events", 1000, 0, 1000.0); histCol.h_M23Mu->Sumw2();
   histCol.h_M34Mu = new TH1F(("h_M34Mu_"+suffix).c_str(), "2 muon invariant mass; mass_{34}; Events", 1000, 0, 1000.0); histCol.h_M34Mu->Sumw2();
   histCol.h_M41Mu = new TH1F(("h_M41Mu_"+suffix).c_str(), "2 muon invariant mass; mass_{41}; Events", 1000, 0, 1000.0); histCol.h_M41Mu->Sumw2();
   histCol.h_M12ElEl = new TH1F(("h_M12ElEl_"+suffix).c_str(), "2 electron invariant mass; mass_{ee}; Events", 1000, 0, 1000.0); histCol.h_M12ElEl->Sumw2();
   histCol.h_M12MuMu = new TH1F(("h_M12MuMu_"+suffix).c_str(), "2 muon invariant mass; mass_{#mu#mu}; Events", 1000, 0, 1000.0); histCol.h_M12MuMu->Sumw2();
}

void writeHistCollection(HistCollection &histCol)
{
   histCol.h_mu_pt_leading->Write();
   histCol.h_mu_pt_second->Write();
   histCol.h_mu_pt_third->Write();
   histCol.h_mu_pt_fourth->Write();
   histCol.h_mu_pt_fifth->Write();
   histCol.h_l_pt_leading->Write();
   histCol.h_l_pt_second->Write();
   histCol.h_l_pt_third->Write();
   histCol.h_l_pt_fourth->Write();
   histCol.h_l_pt_fifth->Write();
   histCol.h_mu_eta_leading->Write();
   histCol.h_mu_eta_second->Write();
   histCol.h_mu_eta_third->Write();
   histCol.h_mu_eta_fourth->Write();
   histCol.h_mu_eta_fifth->Write();
   histCol.h_l_eta_leading->Write();
   histCol.h_l_eta_second->Write();
   histCol.h_l_eta_third->Write();
   histCol.h_l_eta_fourth->Write();
   histCol.h_l_eta_fifth->Write();
   histCol.h_mu_phi_leading->Write();
   histCol.h_mu_phi_second->Write();
   histCol.h_mu_phi_third->Write();
   histCol.h_mu_phi_fourth->Write();
   histCol.h_mu_phi_fifth->Write();
   histCol.h_l_phi_leading->Write();
   histCol.h_l_phi_second->Write();
   histCol.h_l_phi_third->Write();
   histCol.h_l_phi_fourth->Write();
   histCol.h_l_phi_fifth->Write();
   histCol.h_mu_pTRatio_leading->Write();
   histCol.h_mu_pTRatio_second->Write();
   histCol.h_mu_ip3d_leading->Write();
   histCol.h_mu_ip3d_second->Write();
   histCol.h_mu_ip3d_third->Write();
   histCol.h_mu_ip3d_fourth->Write();
   histCol.h_mu_ip3d_fifth->Write();
   histCol.h_l_ip3d_leading->Write();
   histCol.h_l_ip3d_second->Write();
   histCol.h_l_ip3d_third->Write();
   histCol.h_l_ip3d_fourth->Write();
   histCol.h_l_ip3d_fifth->Write();
   histCol.h_mu_d0_leading->Write();
   histCol.h_mu_d0_second->Write();
   histCol.h_mu_d0_third->Write();
   histCol.h_mu_d0_fourth->Write();
   histCol.h_mu_d0_fifth->Write();
   histCol.h_l_d0_leading->Write();
   histCol.h_l_d0_second->Write();
   histCol.h_l_d0_third->Write();
   histCol.h_l_d0_fourth->Write();
   histCol.h_l_d0_fifth->Write();
   histCol.h_mu_dz_leading->Write();
   histCol.h_mu_dz_second->Write();
   histCol.h_mu_dz_third->Write();
   histCol.h_mu_dz_fourth->Write();
   histCol.h_mu_dz_fifth->Write();
   histCol.h_l_dz_leading->Write();
   histCol.h_l_dz_second->Write();
   histCol.h_l_dz_third->Write();
   histCol.h_l_dz_fourth->Write();
   histCol.h_l_dz_fifth->Write();
   histCol.h_mu_pTrel_leading->Write();
   histCol.h_mu_pTrel_second->Write();
   histCol.h_mu_relIso_leading->Write();
   histCol.h_mu_relIso_second->Write();
   histCol.h_mu_deltaR_leading->Write();
   histCol.h_mu_deltaR_second->Write();
   histCol.h_InvariantMass->Write();
   histCol.h_MET->Write();
   histCol.h_nJets->Write();
   histCol.h_nbJets->Write();
   histCol.h_InvariantMassJJ->Write();
   histCol.h_jet_pt_leading->Write();
   histCol.h_jet_pt_second->Write();
   histCol.h_jet_eta_leading->Write();
   histCol.h_jet_eta_second->Write();
   histCol.h_jet_phi_leading->Write();
   histCol.h_jet_csv_leading->Write();
   histCol.h_jet_csv_second->Write();
   histCol.h_bjet_pt_leading->Write();
   histCol.h_bjet_eta_leading->Write();
   histCol.h_bjet_phi_leading->Write();
   histCol.h_bjet_pt_second->Write();
   histCol.h_bjet_eta_second->Write();
   histCol.h_bjet_phi_second->Write();
   histCol.h_el_pt_leading->Write();
   histCol.h_el_pt_second->Write();
   histCol.h_el_pt_third->Write();
   histCol.h_el_pt_fourth->Write();
   histCol.h_el_pt_fifth->Write();
   histCol.h_el_eta_leading->Write();
   histCol.h_el_eta_second->Write();
   histCol.h_el_eta_third->Write();
   histCol.h_el_eta_fourth->Write();
   histCol.h_el_eta_fifth->Write();
   histCol.h_el_phi_leading->Write();
   histCol.h_el_phi_second->Write();
   histCol.h_el_phi_third->Write();
   histCol.h_el_phi_fourth->Write();
   histCol.h_el_phi_fifth->Write();
   histCol.h_el_pTRatio_leading->Write();
   histCol.h_el_pTRatio_second->Write();
   histCol.h_el_ip3d_leading->Write();
   histCol.h_el_ip3d_second->Write();
   histCol.h_el_ip3d_third->Write();
   histCol.h_el_ip3d_fourth->Write();
   histCol.h_el_ip3d_fifth->Write();
   histCol.h_el_d0_leading->Write();
   histCol.h_el_d0_second->Write();
   histCol.h_el_d0_third->Write();
   histCol.h_el_d0_fourth->Write();
   histCol.h_el_d0_fifth->Write();
   histCol.h_el_dz_leading->Write();
   histCol.h_el_dz_second->Write();
   histCol.h_el_dz_third->Write();
   histCol.h_el_dz_fourth->Write();
   histCol.h_el_dz_fifth->Write();
   histCol.h_el_pTrel_leading->Write();
   histCol.h_el_pTrel_second->Write();
   histCol.h_el_relIso_leading->Write();
   histCol.h_el_relIso_second->Write();
   histCol.h_el_deltaR_leading->Write();
   histCol.h_el_deltaR_second->Write();
   histCol.h_nfatJets->Write();
   histCol.h_fatjet_pt_leading->Write();
   histCol.h_fatjet_phi_leading->Write();
   histCol.h_fatjet_eta_leading->Write();
   histCol.h_fatjet_prunMass_leading->Write();
   histCol.h_fatjet_trimMass_leading->Write();
   histCol.h_fatjet_sd0Mass_leading->Write();
   histCol.h_fatjet_tau21_leading->Write();
   histCol.h_M4Mu->Write();
   histCol.h_M4El->Write();
   histCol.h_M12Mu->Write();
   histCol.h_M23Mu->Write();
   histCol.h_M34Mu->Write();
   histCol.h_M41Mu->Write();
   histCol.h_M12ElEl->Write();
   histCol.h_M12MuMu->Write();
}

void fillMuHistCollection(HistCollection &histCol, double mu1pt, double mu2pt, double mu1eta,  double mu2eta, double mu1phi, double mu2phi, double mu1pTRatio, double mu2pTRatio, double mu1ip3D, double mu2ip3D, double mu1d0, double mu2d0, double mu1dz, double mu2dz, double mu1pTRel, double mu2pTRel, double mu1relIso, double mu2relIso, double mu1DeltaR, double mu2DeltaR, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_second->Fill(mu2pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_second->Fill(mu2eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_second->Fill(mu2phi, weight);
  histCol.h_mu_pTRatio_leading->Fill(mu1pTRatio, weight);
  histCol.h_mu_pTRatio_second->Fill(mu2pTRatio, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_mu_ip3d_second->Fill(mu2ip3D, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_mu_d0_second->Fill(mu2d0, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_mu_dz_second->Fill(mu2dz, weight);
  histCol.h_mu_pTrel_leading->Fill(mu1pTRel, weight);
  histCol.h_mu_pTrel_second->Fill(mu2pTRel, weight);
  histCol.h_mu_relIso_leading->Fill(mu1relIso, weight);
  histCol.h_mu_relIso_second->Fill(mu2relIso, weight);
  histCol.h_mu_deltaR_leading->Fill(mu1DeltaR, weight);
  histCol.h_mu_deltaR_second->Fill(mu2DeltaR, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight);     
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight); 
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight);
}

void fillMuHistCollection4l(HistCollection &histCol, double mu1pt, double mu2pt, double mu3pt, double mu4pt, double mu1eta,  double mu2eta, double mu3eta, double mu4eta, double mu1phi,  double mu2phi, double mu3phi, double mu4phi, double mu1ip3D, double mu2ip3D, double mu3ip3D, double mu4ip3D, double mu1d0, double mu2d0, double mu3d0, double mu4d0, double mu1dz, double mu2dz, double mu3dz, double mu4dz, double MET, int njets, int nbjets, double m4l, double m12, double m23, double m34, double m41, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_second->Fill(mu2pt, weight);
  histCol.h_mu_pt_third->Fill(mu3pt, weight);
  histCol.h_mu_pt_fourth->Fill(mu4pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_second->Fill(mu2eta, weight);
  histCol.h_mu_eta_third->Fill(mu3eta, weight);
  histCol.h_mu_eta_fourth->Fill(mu4eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_second->Fill(mu2phi, weight);
  histCol.h_mu_phi_third->Fill(mu3phi, weight);
  histCol.h_mu_phi_fourth->Fill(mu4phi, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_mu_ip3d_second->Fill(mu2ip3D, weight);
  histCol.h_mu_ip3d_third->Fill(mu3ip3D, weight);
  histCol.h_mu_ip3d_fourth->Fill(mu4ip3D, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_mu_d0_second->Fill(mu2d0, weight);
  histCol.h_mu_d0_third->Fill(mu3d0, weight);
  histCol.h_mu_d0_fourth->Fill(mu4d0, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_mu_dz_second->Fill(mu2dz, weight); 
  histCol.h_mu_dz_third->Fill(mu3dz, weight);
  histCol.h_mu_dz_fourth->Fill(mu4dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_M4Mu->Fill(m4l, weight);
  histCol.h_M12Mu->Fill(m12, weight);
  histCol.h_M23Mu->Fill(m23, weight);
  histCol.h_M34Mu->Fill(m34, weight);
  histCol.h_M41Mu->Fill(m41, weight);
}

void fillElHistCollection4l(HistCollection &histCol, double el1pt, double el2pt, double el3pt, double el4pt, double el1eta,  double el2eta, double el3eta, double el4eta, double el1phi,  double el2phi, double el3phi, double el4phi, double el1ip3D, double el2ip3D, double el3ip3D, double el4ip3D, double el1d0, double el2d0, double el3d0, double el4d0, double el1dz, double el2dz, double el3dz, double el4dz, double MET, int njets, int nbjets, double m4l, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_second->Fill(el2pt, weight);
  histCol.h_el_pt_third->Fill(el3pt, weight);
  histCol.h_el_pt_fourth->Fill(el4pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_second->Fill(el2eta, weight);
  histCol.h_el_eta_third->Fill(el3eta, weight);
  histCol.h_el_eta_fourth->Fill(el4eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_second->Fill(el2phi, weight);
  histCol.h_el_phi_third->Fill(el3phi, weight);
  histCol.h_el_phi_fourth->Fill(el4phi, weight);
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_el_ip3d_second->Fill(el2ip3D, weight);
  histCol.h_el_ip3d_third->Fill(el3ip3D, weight);
  histCol.h_el_ip3d_fourth->Fill(el4ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_el_d0_second->Fill(el2d0, weight);
  histCol.h_el_d0_third->Fill(el3d0, weight);
  histCol.h_el_d0_fourth->Fill(el4d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_el_dz_second->Fill(el2dz, weight);
  histCol.h_el_dz_third->Fill(el3dz, weight);
  histCol.h_el_dz_fourth->Fill(el4dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_M4El->Fill(m4l, weight);
}

void fillElMuHistCollection4l(HistCollection &histCol, double el1pt, double el2pt, double mu1pt, double mu2pt, double el1eta,  double el2eta, double mu1eta, double mu2eta, double el1phi,  double el2phi, double mu1phi, double mu2phi, double el1ip3D, double el2ip3D, double mu1ip3D, double mu2ip3D, double el1d0, double el2d0, double mu1d0, double mu2d0, double el1dz, double el2dz, double mu1dz, double mu2dz, double MET, int njets, int nbjets, double m4l, double m12ee, double m12mm, double weight)
{ 
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_second->Fill(el2pt, weight);
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_second->Fill(mu2pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_second->Fill(el2eta, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_second->Fill(mu2eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_second->Fill(el2phi, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_second->Fill(mu2phi, weight);
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_el_ip3d_second->Fill(el2ip3D, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_mu_ip3d_second->Fill(mu2ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_el_d0_second->Fill(el2d0, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_mu_d0_second->Fill(mu2d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_el_dz_second->Fill(el2dz, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_mu_dz_second->Fill(mu2dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_M4El->Fill(m4l, weight);
  histCol.h_M12ElEl->Fill(m12ee, weight);
  histCol.h_M12MuMu->Fill(m12mm, weight);
}

void fillMuHistCollection5l(HistCollection &histCol, double mu1pt, double mu2pt, double mu3pt, double mu4pt, double mu5pt, double mu1eta,  double mu2eta, double mu3eta, double mu4eta, double mu5eta, double mu1phi,  double mu2phi, double mu3phi, double mu4phi, double mu5phi, double mu1ip3D, double mu2ip3D, double mu3ip3D, double mu4ip3D, double mu5ip3D, double mu1d0, double mu2d0, double mu3d0, double mu4d0, double mu5d0, double mu1dz, double mu2dz, double mu3dz, double mu4dz, double mu5dz, double MET, int njets, int nbjets, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_second->Fill(mu2pt, weight);
  histCol.h_mu_pt_third->Fill(mu3pt, weight);
  histCol.h_mu_pt_fourth->Fill(mu4pt, weight);
  histCol.h_mu_pt_fifth->Fill(mu5pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_second->Fill(mu2eta, weight);
  histCol.h_mu_eta_third->Fill(mu3eta, weight);
  histCol.h_mu_eta_fourth->Fill(mu4eta, weight);
  histCol.h_mu_eta_fifth->Fill(mu5eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_second->Fill(mu2phi, weight);
  histCol.h_mu_phi_third->Fill(mu3phi, weight);
  histCol.h_mu_phi_fourth->Fill(mu4phi, weight);
  histCol.h_mu_phi_fifth->Fill(mu5phi, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_mu_ip3d_second->Fill(mu2ip3D, weight);
  histCol.h_mu_ip3d_third->Fill(mu3ip3D, weight);
  histCol.h_mu_ip3d_fourth->Fill(mu4ip3D, weight);
  histCol.h_mu_ip3d_fifth->Fill(mu5ip3D, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_mu_d0_second->Fill(mu2d0, weight);
  histCol.h_mu_d0_third->Fill(mu3d0, weight);
  histCol.h_mu_d0_fourth->Fill(mu4d0, weight);
  histCol.h_mu_d0_fifth->Fill(mu5d0, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_mu_dz_second->Fill(mu2dz, weight);
  histCol.h_mu_dz_third->Fill(mu3dz, weight);
  histCol.h_mu_dz_fourth->Fill(mu4dz, weight);
  histCol.h_mu_dz_fifth->Fill(mu5dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
}

void fillLepHistCollection5l(HistCollection &histCol, double l1pt, double l2pt, double l3pt, double l4pt, double l5pt, double l1eta,  double l2eta, double l3eta, double l4eta, double l5eta, double l1phi,  double l2phi, double l3phi, double l4phi, double l5phi, double l1ip3D, double l2ip3D, double l3ip3D, double l4ip3D, double l5ip3D, double l1d0, double l2d0, double l3d0, double l4d0, double l5d0, double l1dz, double l2dz, double l3dz, double l4dz, double l5dz, double MET, int njets, int nbjets, double weight)
{
  histCol.h_l_pt_leading->Fill(l1pt, weight);
  histCol.h_l_pt_second->Fill(l2pt, weight);
  histCol.h_l_pt_third->Fill(l3pt, weight);
  histCol.h_l_pt_fourth->Fill(l4pt, weight);
  histCol.h_l_pt_fifth->Fill(l5pt, weight);
  histCol.h_l_eta_leading->Fill(l1eta, weight);
  histCol.h_l_eta_second->Fill(l2eta, weight);
  histCol.h_l_eta_third->Fill(l3eta, weight);
  histCol.h_l_eta_fourth->Fill(l4eta, weight);
  histCol.h_l_eta_fifth->Fill(l5eta, weight);
  histCol.h_l_phi_leading->Fill(l1phi, weight);
  histCol.h_l_phi_second->Fill(l2phi, weight);
  histCol.h_l_phi_third->Fill(l3phi, weight);
  histCol.h_l_phi_fourth->Fill(l4phi, weight);
  histCol.h_l_phi_fifth->Fill(l5phi, weight);
  histCol.h_l_ip3d_leading->Fill(l1ip3D, weight);
  histCol.h_l_ip3d_second->Fill(l2ip3D, weight);
  histCol.h_l_ip3d_third->Fill(l3ip3D, weight);
  histCol.h_l_ip3d_fourth->Fill(l4ip3D, weight);
  histCol.h_l_ip3d_fifth->Fill(l5ip3D, weight);
  histCol.h_l_d0_leading->Fill(l1d0, weight);
  histCol.h_l_d0_second->Fill(l2d0, weight);
  histCol.h_l_d0_third->Fill(l3d0, weight);
  histCol.h_l_d0_fourth->Fill(l4d0, weight);
  histCol.h_l_d0_fifth->Fill(l5d0, weight);
  histCol.h_l_dz_leading->Fill(l1dz, weight);
  histCol.h_l_dz_second->Fill(l2dz, weight);
  histCol.h_l_dz_third->Fill(l3dz, weight);
  histCol.h_l_dz_fourth->Fill(l4dz, weight);
  histCol.h_l_dz_fifth->Fill(l5dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
}

void fillElHistCollection5l(HistCollection &histCol, double el1pt, double el2pt, double el3pt, double el4pt, double el5pt, double el1eta,  double el2eta, double el3eta, double el4eta, double el5eta, double el1phi,  double el2phi, double el3phi, double el4phi, double el5phi, double el1ip3D, double el2ip3D, double el3ip3D, double el4ip3D, double el5ip3D, double el1d0, double el2d0, double el3d0, double el4d0, double el5d0, double el1dz, double el2dz, double el3dz, double el4dz, double el5dz, double MET, int njets, int nbjets, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_second->Fill(el2pt, weight);
  histCol.h_el_pt_third->Fill(el3pt, weight);
  histCol.h_el_pt_fourth->Fill(el4pt, weight);
  histCol.h_el_pt_fifth->Fill(el5pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_second->Fill(el2eta, weight);
  histCol.h_el_eta_third->Fill(el3eta, weight);
  histCol.h_el_eta_fourth->Fill(el4eta, weight);
  histCol.h_el_eta_fifth->Fill(el5eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_second->Fill(el2phi, weight);
  histCol.h_el_phi_third->Fill(el3phi, weight);
  histCol.h_el_phi_fourth->Fill(el4phi, weight);
  histCol.h_el_phi_fifth->Fill(el5phi, weight);
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_el_ip3d_second->Fill(el2ip3D, weight);
  histCol.h_el_ip3d_third->Fill(el3ip3D, weight);
  histCol.h_el_ip3d_fourth->Fill(el4ip3D, weight);
  histCol.h_el_ip3d_fifth->Fill(el5ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_el_d0_second->Fill(el2d0, weight);
  histCol.h_el_d0_third->Fill(el3d0, weight);
  histCol.h_el_d0_fourth->Fill(el4d0, weight);
  histCol.h_el_d0_fifth->Fill(el5d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_el_dz_second->Fill(el2dz, weight);
  histCol.h_el_dz_third->Fill(el3dz, weight);
  histCol.h_el_dz_fourth->Fill(el4dz, weight);
  histCol.h_el_dz_fifth->Fill(el5dz, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
}


void fillElHistCollection(HistCollection &histCol, double el1pt, double el2pt, double el1eta,  double el2eta, double el1phi, double el2phi,  double el1pTRatio, double el2pTRatio, double el1ip3D, double el2ip3D, double el1d0, double el2d0, double el1dz, double el2dz, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_second->Fill(el2pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_second->Fill(el2eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_second->Fill(el2phi, weight);
  histCol.h_el_pTRatio_leading->Fill(el1pTRatio, weight);
  histCol.h_el_pTRatio_second->Fill(el2pTRatio, weight); 
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_el_ip3d_second->Fill(el2ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_el_d0_second->Fill(el2d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_el_dz_second->Fill(el2dz, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight);
}

void fillElMuHistCollection(HistCollection &histCol, double el1pt, double mu1pt, double el1eta,  double mu1eta, double el1phi, double mu1phi, double el1pTRatio, double mu1pTRatio, double el1ip3D, double mu1ip3D, double el1d0, double mu1d0, double el1dz, double mu1dz, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_el_pTRatio_leading->Fill(el1pTRatio, weight);
  histCol.h_mu_pTRatio_leading->Fill(mu1pTRatio, weight);
  histCol.h_el_ip3d_leading->Fill(el1ip3D, weight);
  histCol.h_mu_ip3d_leading->Fill(mu1ip3D, weight);
  histCol.h_el_d0_leading->Fill(el1d0, weight);
  histCol.h_mu_d0_leading->Fill(mu1d0, weight);
  histCol.h_el_dz_leading->Fill(el1dz, weight);
  histCol.h_mu_dz_leading->Fill(mu1dz, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight);
}

void fillMuHistCollectionWithFatJet(HistCollection &histCol, double mu1pt, double mu2pt, double mu1eta,  double mu2eta, double mu1phi, double mu2phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi, double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_mu_pt_second->Fill(mu2pt, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_mu_eta_second->Fill(mu2eta, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_mu_phi_second->Fill(mu2phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight); 
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight); 
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight); 
  if(nfatjets > 0) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);   
  if(nfatjets > 0) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets > 0) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight); 
}

void fillElHistCollectionWithFatJet(HistCollection &histCol, double el1pt, double el2pt, double el1eta,  double el2eta, double el1phi, double el2phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi, double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_el_pt_second->Fill(el2pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_el_eta_second->Fill(el2eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_el_phi_second->Fill(el2phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight);
  if(nfatjets >= 1) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);
  if(nfatjets >= 1) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets >= 1) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight);
}

void fillElMuHistCollectionWithFatJet(HistCollection &histCol, double el1pt, double mu1pt, double el1eta,  double mu1eta, double el1phi, double mu1phi, double invmass, double MET, int njets, int nbjets, int nfatjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double fatJetPt, double fatJetEta, double fatJetPhi,  double prunMass, double trimMass, double sd0, double tau21, double weight)
{
  histCol.h_el_pt_leading->Fill(el1pt, weight);
  histCol.h_mu_pt_leading->Fill(mu1pt, weight);
  histCol.h_el_eta_leading->Fill(el1eta, weight);
  histCol.h_mu_eta_leading->Fill(mu1eta, weight);
  histCol.h_el_phi_leading->Fill(el1phi, weight);
  histCol.h_mu_phi_leading->Fill(mu1phi, weight);
  histCol.h_InvariantMass->Fill(invmass, weight);
  histCol.h_MET->Fill(MET, weight);
  histCol.h_nJets->Fill(njets, weight);
  histCol.h_nbJets->Fill(nbjets, weight);
  histCol.h_nfatJets->Fill(nfatjets, weight);
  if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
  if(njets >= 2) histCol.h_jet_pt_second->Fill(jet2pt, weight);
  if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
  if(njets >= 2) histCol.h_jet_eta_second->Fill(jet2eta, weight);
  if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
  if(njets >= 2) histCol.h_jet_phi_second->Fill(jet2phi, weight);
  if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
  if(njets >= 2) histCol.h_jet_csv_second->Fill(jet2csv, weight);
  if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
  if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
  if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
  if(nbjets >= 2) histCol.h_bjet_pt_second->Fill(bjet2pt, weight);
  if(nbjets >= 2) histCol.h_bjet_eta_second->Fill(bjet2eta, weight);
  if(nbjets >= 2) histCol.h_bjet_phi_second->Fill(bjet2phi, weight);
  if(nfatjets >= 1) histCol.h_fatjet_pt_leading->Fill(fatJetPt, weight);
  if(nfatjets >= 1) histCol.h_fatjet_eta_leading->Fill(fatJetEta, weight);
  if(nfatjets >= 1) histCol.h_fatjet_phi_leading->Fill(fatJetPhi, weight);
  if(nfatjets > 0) histCol.h_fatjet_prunMass_leading->Fill(prunMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_trimMass_leading->Fill(trimMass, weight);
  if(nfatjets > 0) histCol.h_fatjet_sd0Mass_leading->Fill(sd0, weight);
  if(nfatjets > 0) histCol.h_fatjet_tau21_leading->Fill(tau21, weight);
}

bool passVetoMuonID(std::vector<float> *lep_eta, std::vector<float> *lep_dz, std::vector<float> *lep_dxy, std::vector<float> *lep_relIso03EA, std::vector<float> *lep_sip3d, int idx)
{
    // - POG Loose
    // Already applied at WVZBabyMaker

    // - |eta| < 2.4
    if (not (fabs(lep_eta->at(idx)) < 2.4)) return false;

    // - dz < 0.1
    if (not (fabs(lep_dz->at(idx)) < 0.1)) return false;

    // - dxy < 0.05
    if (not (fabs(lep_dxy->at(idx)) < 0.05)) return false;

    // - RelIso03EA < 0.2 if tagged as Z
    if (not (fabs(lep_relIso03EA->at(idx)) < 0.4)) return false;

    // - |sip3d| < 5
    if (not (fabs(lep_sip3d->at(idx)) < 5)) return false;

    return true;
}


bool passVetoElectronID(std::vector<float> *lep_eta, std::vector<float> *lep_dz, std::vector<float> *lep_dxy, std::vector<float> *lep_relIso03EA, std::vector<float> *lep_sip3d, int idx)
{

    // - POG MVA Veto (i.e. HZZ w.p.)
    // Already applied at WVZBabyMaker

    // - |eta| < 2.5
    if (not (fabs(lep_eta->at(idx)) < 2.5)) return false;

    // - dz < 0.1
    if (not (fabs(lep_dz->at(idx)) < 0.1)) return false;

    // - dxy < 0.05
    if (not (fabs(lep_dxy->at(idx)) < 0.05)) return false;

    // - RelIso03EA < 0.2 if tagged as Z
    if (not (fabs(lep_relIso03EA->at(idx)) < 0.4)) return false;

    // - |sip3d| < 5
    if (not (fabs(lep_sip3d->at(idx)) < 5)) return false;

    return true;

}

bool passVetoLeptonID(std::vector<int> *lep_id, int idx, std::vector<float> *lep_eta, std::vector<float> *lep_dz, std::vector<float> *lep_dxy, std::vector<float> *lep_relIso03EA, std::vector<float> *lep_sip3d)
{
   if (abs(lep_id->at(idx)) == 11)
     return  passVetoElectronID(lep_eta, lep_dz, lep_dxy, lep_relIso03EA, lep_sip3d, idx);
     else return passVetoMuonID(lep_eta, lep_dz, lep_dxy, lep_relIso03EA, lep_sip3d, idx);
}

int selectVetoLeptons(std::vector<float> *lep_pt, std::vector<int> *lep_id, std::vector<float> *lep_eta, std::vector<float> *lep_dz, std::vector<float> *lep_dxy, std::vector<float> *lep_relIso03EA, std::vector<float> *lep_sip3d)
{
  int nVetoLeptons = 0;
  for (unsigned int kk = 0 ; kk < lep_pt->size() ; kk++)
  {
    if (not passVetoLeptonID(lep_id, kk, lep_eta, lep_dz, lep_dxy, lep_relIso03EA, lep_sip3d)) continue;
      nVetoLeptons++;
  }
  return nVetoLeptons;
}

/*
bool Analysis::Is4LeptonEvent()
{
  return nVetoLeptons == 4;
}

bool Is5LeptonEvent()
{
  return nVetoLeptons == 5;
}*/

double mva_electron_2018(leptonInfo &electron)
{
  float notraw = electron.mva;
  if (notraw >  1.0-1.e-6) notraw =  1.0-1.e-6; // protect against inf, -inf due to FP rounding issues
  if (notraw < -1.0+1.e-6) notraw = -1.0+1.e-6;
  return -0.5*log((2.0/(notraw+1))-1.0);
  /*double c, tau, A;
  c = tau = A = 0.0 
  if(fabs(electron.eta) < 0.8) 
  {
    c   = 7.1336238874;
    tau = 16.5605268797;
    A   = 8.22531222391;
    return 
  }*/
}

