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

//double MUON_MASS = 0;

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
  float pt_corr;
  float phi_corr;
  float eta_corr;
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

bool sortLeptonsInDescendingpT(leptonInfo lep1, leptonInfo lep2)
{
  return (lep1.pt > lep2.pt);
}

int ReadOutputBacon(std::string infile, std::string treeStr, std::string outfile)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  UInt_t          runNumber;
  ULong64_t       evtNumber;
  UInt_t          lumiSection;
  Bool_t          triggerStatus;
  Float_t         eventWeight;
  UInt_t          nPU;
  Float_t         met;
  Float_t         metPhi;
  vector<float>   *muon_pt;
  vector<float>   *muon_phi;
  vector<float>   *muon_eta;
  vector<float>   *muon_pt_corr;
  vector<float>   *muon_phi_corr;
  vector<float>   *muon_eta_corr;
  vector<float>   *muon_iso;
  vector<float>     *muon_charge;
  vector<bool>    *muon_id;
  vector<float>   *muon_recoEW;
  vector<float>   *muon_triggerEW;
  vector<float>   *muon_trigger;
  vector<float>   *jet_pt;
  vector<float>   *jet_eta;
  vector<float>   *jet_phi;
  vector<float>   *jet_mass;
  vector<float>   *jet_csv;
  vector<float>   *ak8jet_prunMass;
  vector<float>   *ak8jet_trimMass;
  vector<float>   *ak8jet_sd0;
  vector<float>   *ak8jet_pt;
  vector<float>   *ak8jet_phi;
  vector<float>   *ak8jet_eta;
  vector<float>   *ak8jet_tau1;
  vector<float>   *ak8jet_tau2;

  muon_pt = 0;
  muon_phi = 0;
  muon_eta = 0;
  muon_pt_corr = 0;
  muon_phi_corr = 0;
  muon_eta_corr = 0;
  muon_iso = 0;
  muon_charge = 0;
  muon_id = 0;
  muon_recoEW = 0;
  muon_triggerEW = 0;
  muon_trigger = 0;
  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_mass = 0;
  jet_csv = 0;
  ak8jet_prunMass = 0;
  ak8jet_trimMass = 0;
  ak8jet_sd0 = 0;
  ak8jet_pt = 0;
  ak8jet_phi = 0;
  ak8jet_eta = 0;
  ak8jet_tau1 = 0;
  ak8jet_tau2 = 0;

  tree->SetBranchAddress("runNumber", &(runNumber));
  tree->SetBranchAddress("evtNumber", &(evtNumber));
  tree->SetBranchAddress("lumiSection", &(lumiSection));
  tree->SetBranchAddress("eventWeight", &(eventWeight));
  tree->SetBranchAddress("nPU", &(nPU));
  tree->SetBranchAddress("met", &(met)); 
  tree->SetBranchAddress("metPhi", &(metPhi));
  tree->SetBranchAddress("muon_pt", &(muon_pt));
  tree->SetBranchAddress("muon_phi", &(muon_phi));
  tree->SetBranchAddress("muon_eta", &(muon_eta));
  tree->SetBranchAddress("muon_pt_corr", &(muon_pt_corr));
  tree->SetBranchAddress("muon_phi_corr", &(muon_phi_corr));
  tree->SetBranchAddress("muon_eta_corr", &(muon_eta_corr));
  tree->SetBranchAddress("muon_iso", &(muon_iso));
  tree->SetBranchAddress("muon_charge", &(muon_charge));
  tree->SetBranchAddress("muon_id", &(muon_id));
  tree->SetBranchAddress("muon_recoEW", &(muon_recoEW));
  tree->SetBranchAddress("muon_triggerEW", &(muon_triggerEW));
  tree->SetBranchAddress("muon_trigger", &(muon_trigger));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_mass", &(jet_mass));
  tree->SetBranchAddress("jet_csv", &(jet_csv));
  tree->SetBranchAddress("ak8jet_prunMass", &(ak8jet_prunMass));
  tree->SetBranchAddress("ak8jet_trimMass", &(ak8jet_trimMass));
  tree->SetBranchAddress("ak8jet_sd0", &(ak8jet_sd0));
  tree->SetBranchAddress("ak8jet_pt", &(ak8jet_pt));
  tree->SetBranchAddress("ak8jet_phi", &(ak8jet_phi));
  tree->SetBranchAddress("ak8jet_eta", &(ak8jet_eta));
  tree->SetBranchAddress("ak8jet_tau1", &(ak8jet_tau1));
  tree->SetBranchAddress("ak8jet_tau2", &(ak8jet_tau2));
  
  TH1F *h_met = new TH1F("h_met", "h_met", 500, 0.0, 1000.0); h_met->Sumw2();
  TH1F *h_Mll = new TH1F("h_Mll", "h_Mll", 500, 0, 1000.0); h_Mll->Sumw2();
  TH1F *h_nMuons = new TH1F("h_nMuons", "h_nMuons", 10, -0.5, 9.5);h_nMuons->Sumw2();
  TH1F *h_nJets = new TH1F("h_nJets", "h_nJets", 20, -0.5, 19.5);h_nJets->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents = " << nEvents << std::endl;
  double eventWith2SSMuons = 0.0;
  //nEvents = 100;
  
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    //if(not(evtNumber==56 or evtNumber==80)) continue;

    //if(not(evtNumber==6577)) continue;

    h_met->Fill(met);

    std::vector<leptonInfo> v_muons; 
    std::vector<leptonInfo> v_looseMuons;
    for (unsigned int imuon=0; imuon<muon_pt->size(); imuon++)
    {
      leptonInfo muon;
      muon.pt = muon_pt->at(imuon); 
      //std::cout << "muon.pt = " << muon.pt << std::endl;
      muon.phi = muon_phi->at(imuon);
      muon.eta = muon_eta->at(imuon);
      muon.pt_corr = muon_pt_corr->at(imuon);
      muon.phi_corr = muon_phi_corr->at(imuon);
      muon.eta_corr = muon_eta_corr->at(imuon);  
      muon.iso = muon_iso->at(imuon);
      muon.charge = (int)muon_charge->at(imuon);
      muon.id = muon_id->at(imuon);
      muon.recoEW = muon_recoEW->at(imuon);
      muon.trigger = muon_trigger->at(imuon);
      muon.triggerEW = muon_triggerEW->at(imuon);
      //if(muon.id==1 and fabs(muon.eta) < 2.4 and muon.pt > 5.0)
      //{
        //v_muons.push_back(muon);
      //}
      if(muon.id==1 and muon.iso/muon.pt_corr < 0.1 and muon.pt_corr > 25.0 and fabs(muon.eta) < 2.4) 
      {
        /*std::cout << muon.pt_corr << ", " << muon.pt
        << ", " << muon.eta << ", " << muon.phi
        << ", " << muon.id 
        << ", " << muon.iso
        << ", " << muon.charge
        << std::endl;*/
        v_muons.push_back(muon);
      }
      if(muon.id==0 and muon.pt > 10 and fabs(muon.eta) < 2.4) v_looseMuons.push_back(muon);
    }//muon loop

    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);  

    std::cout << "nMuons = " << v_muons.size() << std::endl;
    std::cout << "runNumber = " << runNumber << std::endl;
    std::cout << "evtNumber = " << evtNumber << std::endl;
    std::cout << "lumiSection = " << lumiSection << std::endl;    

    h_nMuons->Fill(v_muons.size());

    std::vector<jetInfo> v_jets;
    for (unsigned int ijet=0; ijet<jet_pt->size(); ijet++)
    {
      jetInfo jet;
      jet.jetPt = jet_pt->at(ijet);
      jet.jetEta = jet_eta->at(ijet);
      jet.jetPhi = jet_phi->at(ijet);
      jet.jetMass = jet_mass->at(ijet);
      jet.jetCSV = jet_csv->at(ijet);
      v_jets.push_back(jet);
    }//jet loop

    //std::cout << "v_jets.size() = " << v_jets.size() << std::endl;

    std::vector<jetInfo> v_selectedJets;
    for(unsigned int iselJet=0; iselJet<v_jets.size(); ++iselJet)
    {
      if(fabs(v_jets.at(iselJet).jetEta)<5.0 and v_jets.at(iselJet).jetPt>30.0)
      {
        bool isGoodJet=true;
        TLorentzVector Jet = fillTLorentzVector(v_jets.at(iselJet).jetPt, v_jets.at(iselJet).jetEta, v_jets.at(iselJet).jetPhi, v_jets.at(iselJet).jetMass);
        for(unsigned int iselMuon=0; iselMuon<v_muons.size(); ++iselMuon)
        {
          TLorentzVector Muon = fillTLorentzVector(v_muons.at(iselMuon).pt, v_muons.at(iselMuon).eta, v_muons.at(iselMuon).phi, MUON_MASS);
          double DRjet_mu = Jet.DeltaR(Muon);
          if(DRjet_mu<0.5) isGoodJet=false;
        }//muon loop closed
        if(isGoodJet) v_selectedJets.push_back(v_jets.at(iselJet));  
      }//jet 4 vector closed
    }//selected jet loop closed

    h_nJets->Fill(v_selectedJets.size());

    std::vector<fatJetInfo> v_fatJets;
    for (unsigned int ifatjet=0; ifatjet<ak8jet_pt->size(); ifatjet++)
    {
      fatJetInfo fatjet;
      fatjet.ak8JetPrunmass = ak8jet_prunMass->at(ifatjet);
      fatjet.ak8JetTrimmass = ak8jet_trimMass->at(ifatjet); 
      fatjet.ak8Jetsd0 = ak8jet_sd0->at(ifatjet);
      fatjet.ak8JetPt = ak8jet_pt->at(ifatjet);
      fatjet.ak8JetEta = ak8jet_eta->at(ifatjet);       
      fatjet.ak8JetPhi = ak8jet_phi->at(ifatjet);
      fatjet.ak8JetTau1 = ak8jet_tau1->at(ifatjet);
      fatjet.ak8JetTau2 = ak8jet_tau2->at(ifatjet);
      v_fatJets.push_back(fatjet);
    }
    
    if(v_muons.size()==2)
    {
      //if(sameValHighPrecision(v_muons.at(0).charge*v_muons.at(1).charge, 1) and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0 and v_selectedJets.size() >= 2)
      if(v_muons.at(0).charge*v_muons.at(1).charge==1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0)
      {
        //std::cout << "v_muons.at(0).charge = " << v_muons.at(0).charge << std::endl;
        //std::cout << "v_muons.at(1).charge = " << v_muons.at(1).charge << std::endl;
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        h_Mll->Fill((mu1+mu2).M());
        eventWith2SSMuons++;
        //eventWith2SSMuons += v_muons.at(0).recoEW*v_muons.at(1).recoEW*v_muons.at(0).trigger*v_muons.at(1).trigger*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW));
      }
    }
  }//end of event loop

  std::cout << "eventWith2SSMuons = " << eventWith2SSMuons << std::endl;
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_met->Write();
  h_Mll->Write();
  h_nMuons->Write();
  h_nJets->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl; 


  return 0;

}
