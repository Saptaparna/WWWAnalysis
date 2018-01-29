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

//double MUON_MASS = 105.6583745*10e-03;

double MUON_MASS = 0;

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

typedef struct
{
  float pt; 
  float phi;
  float eta;
  float iso;
  int charge;
  bool id;
  float recoEW;
  float triggerEW;
  float trigger;
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
  TH1F *h_mu_pt_trailing;
  TH1F *h_mu_eta_leading;
  TH1F *h_mu_eta_trailing;
  TH1F *h_mu_phi_leading;
  TH1F *h_mu_phi_trailing;
  TH1F *h_InvariantMass;
  TH1F *h_MET;
  TH1F *h_nJets;
  TH1F *h_nbJets;
  TH1F *h_InvariantMassJJ;
  TH1F *h_jet_pt_leading;
  TH1F *h_jet_pt_trailing;
  TH1F *h_jet_eta_leading;
  TH1F *h_jet_eta_trailing;
  TH1F *h_jet_phi_leading;
  TH1F *h_jet_phi_trailing;
} HistCollection;

void initializeHistCollection(HistCollection &histCol, std::string suffix)
{
   histCol.h_mu_pt_leading = new TH1F(("h_mu_pt_leading_"+suffix).c_str(), "Leading muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_leading->Sumw2();
   histCol.h_mu_pt_trailing = new TH1F(("h_mu_pt_trailing_"+suffix).c_str(), "Trailing muon pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_mu_pt_trailing->Sumw2();
   histCol.h_mu_eta_leading = new TH1F(("h_mu_eta_leading_"+suffix).c_str(), "Leading muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_leading->Sumw2();
   histCol.h_mu_eta_trailing = new TH1F(("h_mu_eta_trailing_"+suffix).c_str(), "Trailing muon #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_mu_eta_trailing->Sumw2();
   histCol.h_mu_phi_leading = new TH1F(("h_mu_phi_leading_"+suffix).c_str(), "Leading muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_leading->Sumw2();
   histCol.h_mu_phi_trailing = new TH1F(("h_mu_phi_trailing_"+suffix).c_str(), "Trailing muon #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_mu_phi_trailing->Sumw2();
   histCol.h_InvariantMass=new TH1F(("h_InvariantMass_"+suffix).c_str(), "Di-lepton invariant mass; M_{ll} [GeV]; Events/GeV", 10000, 0.0, 1000.0); histCol.h_InvariantMass->Sumw2();
   histCol.h_MET=new TH1F(("h_MET_"+suffix).c_str(), "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); histCol.h_MET->Sumw2();
   histCol.h_nJets = new TH1F(("h_nJets_"+suffix).c_str(), "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);histCol.h_nJets->Sumw2();
   histCol.h_nbJets = new TH1F(("h_nbJets_"+suffix).c_str(), "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);histCol.h_nbJets->Sumw2();
   histCol.h_InvariantMassJJ = new TH1F(("h_InvariantMassJJ_"+suffix).c_str(), "Di-jet invariant mass M_{jj} [GeV]; Events/GeV", 10000.0, 0.0, 1000.0);histCol.h_InvariantMassJJ->Sumw2();
   histCol.h_jet_pt_leading=new TH1F(("h_jet_pt_leading_"+suffix).c_str(), "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_leading->Sumw2();
   histCol.h_jet_pt_trailing=new TH1F(("h_jet_pt_trailing_"+suffix).c_str(), "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); histCol.h_jet_pt_trailing->Sumw2();
   histCol.h_jet_eta_leading=new TH1F(("h_jet_eta_leading_"+suffix).c_str(), "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_jet_eta_leading->Sumw2();
   histCol.h_jet_eta_trailing=new TH1F(("h_jet_eta_trailing_"+suffix).c_str(), "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); histCol.h_jet_eta_trailing->Sumw2();
   histCol.h_jet_phi_leading=new TH1F(("h_jet_phi_leading_"+suffix).c_str(), "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_leading->Sumw2();
   histCol.h_jet_phi_trailing=new TH1F(("h_jet_phi_trailing_"+suffix).c_str(), "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_trailing->Sumw2();
}

void writeHistCollection(HistCollection &histCol)
{
   histCol.h_mu_pt_leading->Write();
   histCol.h_mu_pt_trailing->Write();
   histCol.h_mu_eta_leading->Write();
   histCol.h_mu_eta_trailing->Write();
   histCol.h_mu_phi_leading->Write();
   histCol.h_mu_phi_trailing->Write();
   histCol.h_InvariantMass->Write();
   histCol.h_MET->Write();
   histCol.h_nJets->Write();
   histCol.h_nbJets->Write();
   histCol.h_InvariantMassJJ->Write();
   histCol.h_jet_pt_leading->Write();
   histCol.h_jet_pt_trailing->Write();
   histCol.h_jet_eta_leading->Write();
   histCol.h_jet_eta_trailing->Write();
   histCol.h_jet_phi_leading->Write();
   histCol.h_jet_phi_trailing->Write();
}
