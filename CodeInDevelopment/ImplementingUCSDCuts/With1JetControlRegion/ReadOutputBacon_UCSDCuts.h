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
  TH1F *h_jet_csv_leading;
  TH1F *h_jet_csv_trailing;
  TH1F *h_bjet_pt_leading;
  TH1F *h_bjet_eta_leading;
  TH1F *h_bjet_phi_leading;
  TH1F *h_bjet_pt_trailing;
  TH1F *h_bjet_eta_trailing;
  TH1F *h_bjet_phi_trailing;
  TH1F *h_el_pt_leading;
  TH1F *h_el_pt_trailing;
  TH1F *h_el_eta_leading;
  TH1F *h_el_eta_trailing;
  TH1F *h_el_phi_leading;
  TH1F *h_el_phi_trailing;
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
   histCol.h_jet_eta_leading=new TH1F(("h_jet_eta_leading_"+suffix).c_str(), "Leading jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_leading->Sumw2();
   histCol.h_jet_eta_trailing=new TH1F(("h_jet_eta_trailing_"+suffix).c_str(), "Trailing jet #eta; #eta; Events", 1200, -6.0, 6.0); histCol.h_jet_eta_trailing->Sumw2();
   histCol.h_jet_phi_leading=new TH1F(("h_jet_phi_leading_"+suffix).c_str(), "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_leading->Sumw2();
   histCol.h_jet_phi_trailing=new TH1F(("h_jet_phi_trailing_"+suffix).c_str(), "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); histCol.h_jet_phi_trailing->Sumw2();
   histCol.h_jet_csv_leading=new TH1F(("h_jet_csv_leading_"+suffix).c_str(), "Leading jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_leading->Sumw2();
   histCol.h_jet_csv_trailing=new TH1F(("h_jet_csv_trailing_"+suffix).c_str(), "Trailing jet CSV discriminator; CSV discriminator; Events", 100, 0, 1.0);histCol.h_jet_csv_trailing->Sumw2();
   histCol.h_bjet_pt_leading=new TH1F(("h_bjet_pt_leading_"+suffix).c_str(), "Leading b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_leading->Sumw2();
   histCol.h_bjet_eta_leading=new TH1F(("h_bjet_eta_leading_"+suffix).c_str(), "Leading b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_leading->Sumw2();
   histCol.h_bjet_phi_leading=new TH1F(("h_bjet_phi_leading_"+suffix).c_str(), "Leading b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_leading->Sumw2();
   histCol.h_bjet_pt_trailing=new TH1F(("h_bjet_pt_trailing_"+suffix).c_str(), "Trailing b-jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000);histCol.h_bjet_pt_trailing->Sumw2();
   histCol.h_bjet_eta_trailing=new TH1F(("h_bjet_eta_trailing_"+suffix).c_str(), "Trailing b-jet |#eta|; |#eta|; Events", 1200, -6.0, 6.0);histCol.h_bjet_eta_trailing->Sumw2();
   histCol.h_bjet_phi_trailing=new TH1F(("h_bjet_phi_trailing_"+suffix).c_str(), "Trailing b-jet #phi;  #phi; Events", 800, -4.0, 4.0);histCol.h_bjet_phi_trailing->Sumw2();
   histCol.h_el_pt_leading = new TH1F(("h_el_pt_leading_"+suffix).c_str(), "Leading electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_leading->Sumw2();     histCol.h_el_pt_trailing = new TH1F(("h_el_pt_trailing_"+suffix).c_str(), "Trailing electron pT; pT [GeV]; Events/GeV", 1000, 0, 1000); histCol.h_el_pt_trailing->Sumw2();
   histCol.h_el_eta_leading = new TH1F(("h_el_eta_leading_"+suffix).c_str(), "Leading electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_leading->Sumw2();
   histCol.h_el_eta_trailing = new TH1F(("h_el_eta_trailing_"+suffix).c_str(), "Trailing electron #eta ; #eta ; Events", 600, -3.0, 3.0); histCol.h_el_eta_trailing->Sumw2();
   histCol.h_el_phi_leading = new TH1F(("h_el_phi_leading_"+suffix).c_str(), "Leading electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_leading->Sumw2();
   histCol.h_el_phi_trailing = new TH1F(("h_el_phi_trailing_"+suffix).c_str(), "Trailing electron #phi ; #phi ; Events", 800, -4.0, 4.0); histCol.h_el_phi_trailing->Sumw2();
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
   histCol.h_jet_csv_leading->Write();
   histCol.h_jet_csv_trailing->Write();
   histCol.h_bjet_pt_leading->Write();
   histCol.h_bjet_eta_leading->Write();
   histCol.h_bjet_phi_leading->Write();
   histCol.h_bjet_pt_trailing->Write();
   histCol.h_bjet_eta_trailing->Write();
   histCol.h_bjet_phi_trailing->Write();
   histCol.h_el_pt_leading->Write();
   histCol.h_el_pt_trailing->Write();
   histCol.h_el_eta_leading->Write();
   histCol.h_el_eta_trailing->Write();
   histCol.h_el_phi_leading->Write();
   histCol.h_el_phi_trailing->Write();
}

void fillMuHistCollection(HistCollection &histCol, double mu1pt, double mu2pt, double mu1eta,  double mu2eta, double mu1phi, double mu2phi, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
{
   histCol.h_mu_pt_leading->Fill(mu1pt, weight);
   histCol.h_mu_pt_trailing->Fill(mu2pt, weight);
   histCol.h_mu_eta_leading->Fill(mu1eta, weight);
   histCol.h_mu_eta_trailing->Fill(mu2eta, weight);
   histCol.h_mu_phi_leading->Fill(mu1phi, weight);
   histCol.h_mu_phi_trailing->Fill(mu2phi, weight);
   histCol.h_InvariantMass->Fill(invmass, weight);
   histCol.h_MET->Fill(MET, weight);
   histCol.h_nJets->Fill(njets, weight);
   histCol.h_nbJets->Fill(nbjets, weight);
   if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
   if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);     
   if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
   if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
   if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
   if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
   if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
   if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
   if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
   if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
   if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight); 
   if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
   if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
   if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
}

void fillElHistCollection(HistCollection &histCol, double el1pt, double el2pt, double el1eta,  double el2eta, double el1phi, double el2phi, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
{
   histCol.h_el_pt_leading->Fill(el1pt, weight);
   histCol.h_el_pt_trailing->Fill(el2pt, weight);
   histCol.h_el_eta_leading->Fill(el1eta, weight);
   histCol.h_el_eta_trailing->Fill(el2eta, weight);
   histCol.h_el_phi_leading->Fill(el1phi, weight);
   histCol.h_el_phi_trailing->Fill(el2phi, weight);
   histCol.h_InvariantMass->Fill(invmass, weight);
   histCol.h_MET->Fill(MET, weight);
   histCol.h_nJets->Fill(njets, weight);
   histCol.h_nbJets->Fill(nbjets, weight);
   if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
   if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
   if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
   if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
   if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
   if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
   if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
   if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
   if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
   if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
   if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
   if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
   if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
   if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
}

void fillElMuHistCollection(HistCollection &histCol, double el1pt, double mu1pt, double el1eta,  double mu1eta, double el1phi, double mu1phi, double invmass, double MET, int njets, int nbjets, double jet1pt, double jet2pt, double jet1eta, double jet2eta, double jet1phi, double jet2phi, double jet1csv, double jet2csv, double bjet1pt, double bjet1eta, double bjet1phi, double bjet2pt, double bjet2eta, double bjet2phi, double weight)
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
    if(njets >= 2) histCol.h_jet_pt_leading->Fill(jet1pt, weight);
    if(njets >= 2) histCol.h_jet_pt_trailing->Fill(jet2pt, weight);
    if(njets >= 2) histCol.h_jet_eta_leading->Fill(jet1eta, weight);
    if(njets >= 2) histCol.h_jet_eta_trailing->Fill(jet2eta, weight);
    if(njets >= 2) histCol.h_jet_phi_leading->Fill(jet1phi, weight);
    if(njets >= 2) histCol.h_jet_phi_trailing->Fill(jet2phi, weight);
    if(njets >= 2) histCol.h_jet_csv_leading->Fill(jet1csv, weight);
    if(njets >= 2) histCol.h_jet_csv_trailing->Fill(jet2csv, weight);
    if(nbjets >= 1) histCol.h_bjet_pt_leading->Fill(bjet1pt, weight);
    if(nbjets >= 1) histCol.h_bjet_eta_leading->Fill(bjet1eta, weight);
    if(nbjets >= 1) histCol.h_bjet_phi_leading->Fill(bjet1phi, weight);
    if(nbjets >= 2) histCol.h_bjet_pt_trailing->Fill(bjet2pt, weight);
    if(nbjets >= 2) histCol.h_bjet_eta_trailing->Fill(bjet2eta, weight);
    if(nbjets >= 2) histCol.h_bjet_phi_trailing->Fill(bjet2phi, weight);
}

double weightMu(std::vector<leptonInfo> v_mu)
{
  return v_mu.at(0).recoEW*v_mu.at(1).recoEW*(1 - (1 - v_mu.at(0).triggerEW)*(1 - v_mu.at(1).triggerEW));
}

double weightEl(std::vector<leptonInfo> v_el)
{
  return v_el.at(0).recoEW*v_el.at(1).recoEW*(1 - (1 - v_el.at(0).triggerEW)*(1 - v_el.at(1).triggerEW));
}

double weightElMu(std::vector<leptonInfo> v_el, std::vector<leptonInfo> v_mu)
{
  return v_el.at(0).recoEW*v_mu.at(0).recoEW*(1 - (1 - v_el.at(0).triggerEW)*(1 - v_mu.at(0).triggerEW));
}
