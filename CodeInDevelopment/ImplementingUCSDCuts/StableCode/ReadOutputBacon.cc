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

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  TLorentzVector lep1;
  float lep1Isolation; 
  float lep1Charge;
  float lep1Flavor;
  float lep1Trigger;
  TLorentzVector lep2;
  float lep2Isolation;
  float lep2Charge;
  float lep2Flavor;
  float lep2Trigger;
} leptonInfo;

typedef struct
{ 
  float jetpt;
  float jeteta;
  float jetphi;
  float jetmass;
  float jetd0;
  float jettag;
  float jetflavor;
} JetInfo;

typedef struct
{
  float prunMass;
  float trimMass;
  float sd0;
  float pt;
  float tau1;
  float tau2;
  float tau3;
} FatJetInfo;

typedef struct
{  
  float genpt;
  float geneta;
  float genphi;
  float genmass;
  float genpdgId;
  float genstatus; 
  float genparent;
} GenInfo;

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.jetpt > jet2.jetpt);
}

bool sortFatJetsInDescendingpT(FatJetInfo fatjet1, FatJetInfo fatjet2)
{
  return (fatjet1.pt > fatjet2.pt);
}

bool sortGenPartsInDescendingpT(GenInfo gen1, GenInfo gen2)
{
  return (gen1.genpt > gen2.genpt);
}

int ReadOutputBacon(std::string infile, std::string outfile)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("tree_TTJets");
  //TChain *tree=new TChain("tree_DYJets");
  //TChain *tree=new TChain("tree_gb2XB_fatjet");
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
  TLorentzVector  *leptonOneP4;
  Float_t         leptonOneIso;
  Int_t           leptonOneQ;
  Int_t           leptonOneFlavor;
  Bool_t          leptonOneTrigger;
  TLorentzVector  *leptonTwoP4;
  Float_t         leptonTwoIso;
  Int_t           leptonTwoQ;
  Int_t           leptonTwoFlavor;
  Bool_t          leptonTwoTrigger;
  TLorentzVector  *jetOneP4;
  Float_t         jetOneD0;
  Float_t         jetOneTag;
  Int_t           jetOneFlavor;
  TLorentzVector  *jetTwoP4;
  Float_t         jetTwoD0;
  Float_t         jetTwoTag;
  Int_t           jetTwoFlavor;
  TLorentzVector  *bjetP4;
  Float_t         bjetD0;
  Float_t         bjetTag;
  Int_t           bjetFlavor;
  TLorentzVector  *genBJetP4;
  Float_t         genBJetTag;
  TLorentzVector  *genJetOneP4;
  TLorentzVector  *genJetTwoP4;
  TLorentzVector  *genJetThreeP4;
  Float_t         genJetOneTag;
  Float_t         genJetTwoTag;
  Float_t         genJetThreeTag;
  UInt_t          nMuons;
  UInt_t          nElectrons;
  UInt_t          nJets;
  UInt_t          nFwdJets;
  UInt_t          nBJets;
  std::vector<float>   *jet_pt;
  std::vector<float>   *jet_eta;
  std::vector<float>   *jet_phi;
  std::vector<float>   *jet_mass;
  std::vector<float>   *jet_D0;
  std::vector<float>   *jet_csv_data;
  std::vector<float>   *jet_csv_mc;
  std::vector<float>   *jet_csv;
  std::vector<float>   *jet_flavor;
  std::vector<float>   *ak8jet_prunMass;
  std::vector<float>   *ak8jet_trimMass;
  std::vector<float>   *ak8jet_sd0;
  std::vector<float>   *ak8jet_pt;
  std::vector<float>   *ak8jet_tau1;
  std::vector<float>   *ak8jet_tau2;
  std::vector<float>   *ak8jet_tau3;
  std::vector<float>   *ca15jet_prunMass;
  std::vector<float>   *ca15jet_trimMass;
  std::vector<float>   *ca15jet_sd0;
  std::vector<float>   *ca15jet_pt;
  std::vector<float>   *ca15jet_tau1;
  std::vector<float>   *ca15jet_tau2;
  std::vector<float>   *ca15jet_tau3;
  std::vector<float>   *gen_pdgId;
  std::vector<float>   *gen_status;
  std::vector<float>   *gen_parent;
  std::vector<float>   *gen_pt;
  std::vector<float>   *gen_eta;
  std::vector<float>   *gen_phi;
  std::vector<float>   *gen_mass;

  leptonOneP4 = 0;
  leptonTwoP4 = 0;
  jetOneP4 = 0;
  jetTwoP4 = 0;
  bjetP4 = 0;
  genBJetP4 = 0;
  genJetOneP4 = 0;
  genJetTwoP4 = 0;
  genJetThreeP4 = 0;
  jet_pt = 0;
  jet_eta = 0;
  jet_phi = 0;
  jet_mass = 0;
  jet_D0 = 0;
  jet_csv_data = 0;
  jet_csv_mc = 0;
  jet_csv = 0;
  jet_flavor = 0;
  ak8jet_prunMass = 0;
  ak8jet_trimMass = 0;
  ak8jet_sd0 = 0;
  ak8jet_pt = 0;
  ak8jet_tau1 = 0;
  ak8jet_tau2 = 0;
  ak8jet_tau3 = 0;
  ca15jet_prunMass = 0;
  ca15jet_trimMass = 0;
  ca15jet_sd0 = 0;
  ca15jet_pt = 0;
  ca15jet_tau1 = 0;
  ca15jet_tau2 = 0;
  ca15jet_tau3 = 0;
  gen_pdgId = 0;
  gen_status = 0;
  gen_parent = 0;
  gen_pt = 0;
  gen_eta = 0;
  gen_phi = 0;
  gen_mass = 0;

  tree->SetBranchAddress("runNumber", &(runNumber));
  tree->SetBranchAddress("evtNumber", &(evtNumber));
  tree->SetBranchAddress("lumiSection", &(lumiSection));
  tree->SetBranchAddress("eventWeight", &(eventWeight));
  tree->SetBranchAddress("nPU", &(nPU));
  tree->SetBranchAddress("met", &(met)); 
  tree->SetBranchAddress("metPhi", &(metPhi));
  tree->SetBranchAddress("leptonOneP4", &(leptonOneP4));
  tree->SetBranchAddress("leptonOneIso", &(leptonOneIso));
  tree->SetBranchAddress("leptonOneQ", &(leptonOneQ));
  tree->SetBranchAddress("leptonOneFlavor", &(leptonOneFlavor));
  tree->SetBranchAddress("leptonOneTrigger", &(leptonOneTrigger));
  tree->SetBranchAddress("leptonTwoP4", &(leptonTwoP4));
  tree->SetBranchAddress("leptonTwoIso", &(leptonTwoIso));
  tree->SetBranchAddress("leptonTwoQ", &(leptonTwoQ));
  tree->SetBranchAddress("leptonTwoFlavor", &(leptonTwoFlavor));
  tree->SetBranchAddress("leptonTwoTrigger", &(leptonTwoTrigger));
  tree->SetBranchAddress("jetOneP4", &(jetOneP4));
  tree->SetBranchAddress("jetOneD0", &(jetOneD0));
  tree->SetBranchAddress("jetOneTag", &(jetOneTag));
  tree->SetBranchAddress("jetOneFlavor", &(jetOneFlavor));
  tree->SetBranchAddress("jetTwoP4", &(jetTwoP4));
  tree->SetBranchAddress("jetTwoD0", &(jetTwoD0));
  tree->SetBranchAddress("jetTwoTag", &(jetTwoTag));
  tree->SetBranchAddress("jetTwoFlavor", &(jetTwoFlavor));
  tree->SetBranchAddress("bjetP4", &(bjetP4));
  tree->SetBranchAddress("bjetD0", &(bjetD0));
  tree->SetBranchAddress("bjetTag", &(bjetTag));
  tree->SetBranchAddress("bjetFlavor", &(bjetFlavor));
  tree->SetBranchAddress("genBJetP4", &(genBJetP4));
  tree->SetBranchAddress("genBJetTag", &(genBJetTag));
  tree->SetBranchAddress("genJetOneP4", &(genJetOneP4));
  tree->SetBranchAddress("genJetOneTag", &(genJetOneTag));
  tree->SetBranchAddress("genJetTwoP4", &(genJetTwoP4));
  tree->SetBranchAddress("genJetTwoTag", &(genJetTwoTag));
  tree->SetBranchAddress("genJetThreeP4", &(genJetThreeP4));
  tree->SetBranchAddress("genJetThreeTag", &(genJetThreeTag));
  tree->SetBranchAddress("nMuons", &(nMuons));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("nFwdJets", &(nFwdJets));
  tree->SetBranchAddress("nBJets", &(nBJets));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_mass", &(jet_mass));
  tree->SetBranchAddress("jet_csv_data", &(jet_csv_data));
  tree->SetBranchAddress("jet_csv_mc", &(jet_csv_mc));
  tree->SetBranchAddress("jet_csv", &(jet_csv));
  tree->SetBranchAddress("jet_D0", &(jet_D0));
  tree->SetBranchAddress("jet_flavor", &(jet_flavor));
  tree->SetBranchAddress("ak8jet_prunMass", &(ak8jet_prunMass));
  tree->SetBranchAddress("ak8jet_trimMass", &(ak8jet_trimMass));
  tree->SetBranchAddress("ak8jet_sd0", &(ak8jet_sd0));
  tree->SetBranchAddress("ak8jet_pt", &(ak8jet_pt));
  tree->SetBranchAddress("ak8jet_tau1", &(ak8jet_tau1));
  tree->SetBranchAddress("ak8jet_tau2", &(ak8jet_tau2));
  tree->SetBranchAddress("ak8jet_tau3", &(ak8jet_tau3));
  tree->SetBranchAddress("ca15jet_prunMass", &(ca15jet_prunMass));
  tree->SetBranchAddress("ca15jet_trimMass", &(ca15jet_trimMass));
  tree->SetBranchAddress("ca15jet_sd0", &(ca15jet_sd0));
  tree->SetBranchAddress("ca15jet_pt", &(ca15jet_pt));
  tree->SetBranchAddress("ca15jet_tau1", &(ca15jet_tau1));
  tree->SetBranchAddress("ca15jet_tau2", &(ca15jet_tau2));
  tree->SetBranchAddress("ca15jet_tau3", &(ca15jet_tau3));
  tree->SetBranchAddress("gen_pdgId", &(gen_pdgId));
  tree->SetBranchAddress("gen_status", &(gen_status));
  tree->SetBranchAddress("gen_parent", &(gen_parent));
  tree->SetBranchAddress("gen_pt", &(gen_pt));
  tree->SetBranchAddress("gen_eta", &(gen_eta));
  tree->SetBranchAddress("gen_phi", &(gen_phi));
  tree->SetBranchAddress("gen_mass", &(gen_mass));

  TH1F *h_deltaR = new TH1F("h_deltaR", "h_deltaR", 10000, 0.0, 10.0); h_deltaR->Sumw2();
  TH1F *h_deltaR1 = new TH1F("h_deltaR1", "h_deltaR1", 10000, 0.0, 10.0); h_deltaR1->Sumw2();
  TH1F *h_deltaR2 = new TH1F("h_deltaR2", "h_deltaR2", 10000, 0.0, 10.0); h_deltaR2->Sumw2();
  TH1F *h_Mqq = new TH1F("h_Mqq", "h_Mqq", 500, 0.0, 100.0); h_Mqq->Sumw2();
  TH1F *h_Mmumu = new TH1F("h_Mmumu", "h_Mmumu", 500, 0, 100); h_Mmumu->Sumw2();
  TH1F *h_Mjj = new TH1F("h_Mjj", "h_Mjj", 500, 0.0, 100.0); h_Mjj->Sumw2();
  TH1F *h_njets = new TH1F("h_njets", "h_njets", 10, -0.5, 9.5); h_njets->Sumw2();
  TH1F *h_jetMass1 = new TH1F("h_jetMass1", "h_jetMass1", 500, 0.0, 100.0);h_jetMass1->Sumw2();
  TH1F *h_jetMass2 = new TH1F("h_jetMass2", "h_jetMass2", 500, 0.0, 100.0);h_jetMass2->Sumw2();
  TH1F *h_jetMass = new TH1F("h_jetMass", "h_jetMass", 500, 0.0, 100.0);h_jetMass->Sumw2();
  TH1F *h_ak8jetPrunMass = new TH1F("h_ak8jetPrunMass", "h_ak8jetPrunMass", 500, 0.0, 100.0);h_ak8jetPrunMass->Sumw2();
  TH1F *h_ak8jetTrimMass = new TH1F("h_ak8jetTrimMass", "h_ak8jetTrimMass", 500, 0.0, 100.0);h_ak8jetTrimMass->Sumw2();
  TH1F *h_ak8jetsd0 = new TH1F("h_ak8jetsd0", "h_ak8jetsd0", 500, 0.0, 100.0); h_ak8jetsd0->Sumw2();
  TH1F *h_ca15jetPrunMass = new TH1F("h_ca15jetPrunMass", "h_ca15jetPrunMass", 500, 0.0, 100.0);h_ca15jetPrunMass->Sumw2();
  TH1F *h_ca15jetTrimMass = new TH1F("h_ca15jetTrimMass", "h_ca15jetTrimMass", 500, 0.0, 100.0);h_ca15jetTrimMass->Sumw2();
  TH1F *h_ca15jetsd0 = new TH1F("h_ca15jetsd0", "h_ca15jetsd0", 500, 0.0, 100.0); h_ca15jetsd0->Sumw2();
  TH1F *h_nsubJettiness = new TH1F("h_nsubJettiness", "h_nsubJettiness; #tau21; Events", 1000, 0.0, 1.0);h_nsubJettiness->Sumw2();

  int n_bothMatched=0;
  int n_oneMatched=0;

  int nEvents=tree->GetEntries();
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    std::vector<JetInfo> jets;
    for (unsigned int k=0; k<jet_pt->size(); k++)
    {
      JetInfo jet;
      jet.jetpt = jet_pt->at(k);
      jet.jeteta = jet_eta->at(k);
      jet.jetphi = jet_phi->at(k);
      jet.jetmass = jet_mass->at(k);
      jet.jetd0 = jet_D0->at(k);
      jet.jettag = jet_csv_mc->at(k);
      jet.jetflavor = jet_flavor->at(k);
      if(fabs(jet.jetflavor) < 5) jets.push_back(jet);
    }
    std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);

    h_njets->Fill(jets.size());
    TLorentzVector jet1, jet2;
    if(jets.size() > 1)
    {
      jet1.SetPtEtaPhiM(jets.at(0).jetpt, jets.at(0).jeteta, jets.at(0).jetphi, jets.at(0).jetmass);
      jet2.SetPtEtaPhiM(jets.at(1).jetpt, jets.at(1).jeteta, jets.at(1).jetphi, jets.at(1).jetmass);
    } 
    
    std::vector<GenInfo> genParts; 
    for (unsigned int m=0; m<gen_pt->size(); m++)  
    {
      GenInfo genPart; 
      genPart.genpt = gen_pt->at(m);
      genPart.geneta = gen_eta->at(m);
      genPart.genphi = gen_phi->at(m);
      genPart.genmass = gen_mass->at(m);
      genPart.genpdgId = gen_pdgId->at(m);
      genPart.genstatus = gen_status->at(m);
      genPart.genparent = gen_parent->at(m);
      if(abs(genPart.genpdgId) < 5 and (genPart.genstatus>= 23 and genPart.genstatus<=30)) genParts.push_back(genPart);
    }
    std::sort (genParts.begin(), genParts.end(), sortGenPartsInDescendingpT);

    TLorentzVector gen1, gen2;
    if(genParts.size() > 1)
    {
      gen1.SetPtEtaPhiM(genParts.at(0).genpt, genParts.at(0).geneta, genParts.at(0).genphi, genParts.at(0).genmass);
      gen2.SetPtEtaPhiM(genParts.at(1).genpt, genParts.at(1).geneta, genParts.at(1).genphi, genParts.at(1).genmass);
      h_Mqq->Fill((gen1+gen2).M());
    }
    int i_matchedJet1 = -1;
    int i_matchedJet2 = -1;
    double mindR1 = 0.4;
    double mindR2 = 0.4;
    for(unsigned int p=0; p<jets.size(); p++)
    {
      TLorentzVector l_recoJet1;
      l_recoJet1.SetPtEtaPhiM(jets.at(p).jetpt, jets.at(p).jeteta, jets.at(p).jetphi, jets.at(p).jetmass); 
      if(gen1.DeltaR(l_recoJet1) < mindR1)
      {
        mindR1 = gen1.DeltaR(l_recoJet1);
        i_matchedJet1 = p;
      }
     }
     h_deltaR1->Fill(mindR1);

     //if(i_matchedJet1!=-1) h_jetMass1->Fill(jets.at(i_matchedJet1).jetmass);

     for(unsigned int q=0; q<jets.size(); q++)
     {  
       TLorentzVector l_recoJet2;
       l_recoJet2.SetPtEtaPhiM(jets.at(q).jetpt, jets.at(q).jeteta, jets.at(q).jetphi, jets.at(q).jetmass);
       if(gen2.DeltaR(l_recoJet2) < mindR2)
       {
         mindR2 = gen2.DeltaR(l_recoJet2);
         i_matchedJet2 = q;
         //h_deltaR2->Fill(mindR2);
       }
     }
     h_deltaR2->Fill(mindR2);

     //if(i_matchedJet2!=-1) h_jetMass2->Fill(jets.at(i_matchedJet2).jetmass);
     
     if((i_matchedJet1 != -1 or i_matchedJet2 != -1) and i_matchedJet2!=i_matchedJet1) 
     {
       if(i_matchedJet1 != -1) h_jetMass1->Fill(jets.at(i_matchedJet1).jetmass);
       if(i_matchedJet2 != -1) h_jetMass2->Fill(jets.at(i_matchedJet2).jetmass);
       if(i_matchedJet1 != -1) h_jetMass->Fill(jets.at(i_matchedJet1).jetmass);
       else if(i_matchedJet2 != -1) h_jetMass->Fill(jets.at(i_matchedJet2).jetmass);
       n_oneMatched++;    
     }

     if(i_matchedJet1 != -1 and i_matchedJet2 != -1 and i_matchedJet2!=i_matchedJet1)
     { 
       n_bothMatched++; 
       TLorentzVector matchedJet1;
       TLorentzVector matchedJet2;
       matchedJet1.SetPtEtaPhiM(jets.at(i_matchedJet1).jetpt, jets.at(i_matchedJet1).jeteta, jets.at(i_matchedJet1).jetphi, jets.at(i_matchedJet1).jetmass);
       matchedJet2.SetPtEtaPhiM(jets.at(i_matchedJet2).jetpt, jets.at(i_matchedJet2).jeteta, jets.at(i_matchedJet2).jetphi, jets.at(i_matchedJet2).jetmass);
       h_Mjj->Fill((matchedJet1+matchedJet2).M());
       h_deltaR->Fill(matchedJet1.DeltaR(matchedJet2));
     }

    std::vector<FatJetInfo> ak8fatjets;
    for(unsigned int i=0; i<ak8jet_prunMass->size(); i++)
    {
      FatJetInfo ak8fatjet;
      ak8fatjet.prunMass = ak8jet_prunMass->at(i);
      ak8fatjet.trimMass = ak8jet_trimMass->at(i);
      ak8fatjet.sd0 = ak8jet_sd0->at(i);
      ak8fatjet.pt = ak8jet_pt->at(i);
      ak8fatjet.tau1 = ak8jet_tau1->at(i);
      ak8fatjet.tau2 = ak8jet_tau2->at(i);
      ak8fatjet.tau3 = ak8jet_tau3->at(i);
      ak8fatjets.push_back(ak8fatjet);
    }
    std::sort(ak8fatjets.begin(), ak8fatjets.end(), sortFatJetsInDescendingpT); 

    if(ak8fatjets.size() > 0) h_ak8jetPrunMass->Fill(ak8fatjets.at(0).prunMass);
    if(ak8fatjets.size() > 0) h_ak8jetTrimMass->Fill(ak8fatjets.at(0).trimMass);  
    if(ak8fatjets.size() > 0) h_ak8jetsd0->Fill(ak8fatjets.at(0).sd0);
    if(ak8fatjets.size() > 0 and ak8fatjets.at(0).tau1 > 0.0) h_nsubJettiness->Fill((double)ak8fatjets.at(0).tau2/(double)ak8fatjets.at(0).tau1);

    std::vector<FatJetInfo> ca15fatjets;
    for(unsigned int i=0; i<ca15jet_prunMass->size(); i++)
    {
      FatJetInfo ca15fatjet;
      ca15fatjet.prunMass = ca15jet_prunMass->at(i);
      ca15fatjet.trimMass = ca15jet_trimMass->at(i);
      ca15fatjet.sd0 = ca15jet_sd0->at(i);
      //ca15fatjet.pt = ca15jet_pt->at(i);
      ca15fatjets.push_back(ca15fatjet);
    }
    //std::sort(ca15fatjets.begin(), ca15fatjets.end(), sortFatJetsInDescendingpT);    
 
    if(ca15fatjets.size() > 0) h_ca15jetPrunMass->Fill(ca15fatjets.at(0).prunMass);
    if(ca15fatjets.size() > 0) h_ca15jetTrimMass->Fill(ca15fatjets.at(0).trimMass);
    if(ca15fatjets.size() > 0) h_ca15jetsd0->Fill(ca15fatjets.at(0).sd0);

    TLorentzVector mu1;   
    TLorentzVector mu2;
    mu1.SetPtEtaPhiE(leptonOneP4->Pt(), leptonOneP4->Eta(), leptonOneP4->Phi(), leptonOneP4->E());
    mu2.SetPtEtaPhiE(leptonTwoP4->Pt(), leptonTwoP4->Eta(), leptonTwoP4->Phi(), leptonTwoP4->E());
    
    h_Mmumu->Fill((mu1+mu2).M());

  }//end of event loop
  std::cout << "n_bothMatched = " << n_bothMatched << std::endl;
  std::cout << "n_oneMatched = " << n_oneMatched << std::endl;

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_njets->Write();
  h_deltaR->Write();
  h_deltaR1->Write();
  h_deltaR2->Write();
  h_Mqq->Write();
  h_Mmumu->Write();
  h_Mjj->Write();
  h_jetMass1->Write();
  h_jetMass2->Write();
  h_jetMass->Write(); 
  h_ak8jetPrunMass->Write();
  h_ak8jetTrimMass->Write();
  h_ak8jetsd0->Write();
  h_nsubJettiness->Write();
  h_ca15jetPrunMass->Write();
  h_ca15jetTrimMass->Write();
  h_ca15jetsd0->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  
  return 0; 
}
