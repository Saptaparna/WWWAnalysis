#include <TH1D.h>
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
#include <TVector2.h>
#include <TF1.h>

using namespace std;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double M)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, M);
  return object_p4;
}

typedef struct
{
  double pT;
  double eta;
  double phi;
  int charge;
} LeptonInfo;

typedef struct
{
  double pT;
  double eta;
  double phi;
} PhotonInfo;

typedef struct
{
  double pT;
  double eta;
  double phi;
  double mass;
  int btag;
} JetInfo;

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

int ReadPPZFiles_Delphes(std::string infile, std::string outfile){

  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("LowPtSUSY_Tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  std::vector<double>   *ph_pt;
  std::vector<double>   *ph_phi;
  vector<double>   *ph_eta;
  Int_t           nPhotons;
  vector<double>   *el_pt;
  vector<double>   *el_phi;
  vector<double>   *el_eta;
  vector<int>     *el_charge;
  Int_t           nElectrons;
  vector<double>   *mu_pt;
  vector<double>   *mu_phi;
  vector<double>   *mu_eta;
  vector<int>     *mu_charge;
  Int_t           nMuons;
  vector<double>   *jet_pt;
  vector<double>   *jet_phi;
  vector<double>   *jet_eta;
  vector<double>   *jet_mass;
  vector<int>     *jet_btag;
  Int_t           nJets;
  Float_t         MET;
  Float_t         MET_Phi; 

  ph_pt = 0;
  ph_phi = 0;
  ph_eta = 0;
  el_pt = 0;
  el_phi = 0;
  el_eta = 0;
  el_charge = 0;
  mu_pt = 0;
  mu_phi = 0;
  mu_eta = 0;
  mu_charge = 0;
  jet_pt = 0;
  jet_phi = 0;
  jet_eta = 0;
  jet_mass = 0;
  jet_btag = 0;
  tree->SetBranchAddress("ph_pt", &(ph_pt));
  tree->SetBranchAddress("ph_phi", &(ph_phi));
  tree->SetBranchAddress("ph_eta", &(ph_eta));
  tree->SetBranchAddress("nPhotons", &(nPhotons));
  tree->SetBranchAddress("el_pt", &(el_pt));
  tree->SetBranchAddress("el_eta", &(el_eta));
  tree->SetBranchAddress("el_phi", &(el_phi));
  tree->SetBranchAddress("el_charge", &(el_charge));
  tree->SetBranchAddress("nElectrons", &(nElectrons));
  tree->SetBranchAddress("mu_pt", &(mu_pt));
  tree->SetBranchAddress("mu_eta", &(mu_eta));
  tree->SetBranchAddress("mu_phi", &(mu_phi));
  tree->SetBranchAddress("mu_charge", &(mu_charge));
  tree->SetBranchAddress("nMuons", &(nMuons));
  tree->SetBranchAddress("jet_pt", &(jet_pt));
  tree->SetBranchAddress("jet_phi", &(jet_phi));
  tree->SetBranchAddress("jet_eta", &(jet_eta));
  tree->SetBranchAddress("jet_mass", &(jet_mass));
  tree->SetBranchAddress("jet_btag", &(jet_btag));
  tree->SetBranchAddress("nJets", &(nJets));
  tree->SetBranchAddress("MET", &(MET));
  tree->SetBranchAddress("MET_Phi", &(MET_Phi));

  TH1D *h_jet_pt_leading=new TH1D("h_jet_pt_leading", "Leading jet pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_jet_pt_leading->Sumw2();
  TH1D *h_jet_pt_trailing=new TH1D("h_jet_pt_trailing", "Trailing jet pT; pT [GeV]; Events/GeV", 1000, 0, 1000); h_jet_pt_trailing->Sumw2();

  TH1D *h_jet_eta_leading=new TH1D("h_jet_eta_leading", "Leading jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_leading->Sumw2();
  TH1D *h_jet_eta_trailing=new TH1D("h_jet_eta_trailing", "Trailing jet #eta; #eta; Events", 600, -3.0, 3.0); h_jet_eta_trailing->Sumw2();

  TH1D *h_jet_phi_leading=new TH1D("h_jet_phi_leading", "Leading jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_leading->Sumw2();
  TH1D *h_jet_phi_trailing=new TH1D("h_jet_phi_trailing", "Trailing jet #phi; #phi; Events/GeV", 800, -4.0, 4.0); h_jet_phi_trailing->Sumw2();

  TH1D *h_nJets = new TH1D("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();

  TH1D *h_WJ_PT_ElEl = new TH1D("h_WJ_PT_ElEl", "W_{J} pT; W_{J} pT [GeV] (Electron Channel); Events", 1000, 0, 1000); h_WJ_PT_ElEl->Sumw2();
  TH1D *h_WJ_PT_MuMu = new TH1D("h_WJ_PT_MuMu", "W_{J} pT; W_{J} pT [GeV] (Muon Channel); Events", 1000, 0, 1000); h_WJ_PT_MuMu->Sumw2();
  TH1D *h_WJ_PT_ElMu = new TH1D("h_WJ_PT_ElMu", "W_{J} pT; W_{J} pT [GeV] (Electron Muon Channel); Events", 1000, 0, 1000); h_WJ_PT_ElMu->Sumw2();

  TH1D *h_Electron0_Charge = new TH1D("h_Electron0_Charge", "h_Electron0_Charge; Charge of the leading electron; Events", 10, -5.5, 4.5); h_Electron0_Charge->Sumw2();
  TH1D *h_Electron1_Charge = new TH1D("h_Electron1_Charge", "h_Electron1_Charge; Charge of the sub-leading electron; Events", 10, -5.5, 4.5); h_Electron1_Charge->Sumw2();
  TH1D *h_Muon0_Charge = new TH1D("h_Muon0_Charge", "h_Muon0_Charge; Charge of the leading muon; Events", 10, -5.5, 4.5); h_Muon0_Charge->Sumw2();
  TH1D *h_Muon1_Charge = new TH1D("h_Muon1_Charge", "h_Muon1_Charge; Charge of the sub-leading muon; Events", 10, -5.5, 4.5); h_Muon1_Charge->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  int n_elec = 0;
  int n_muon = 0;
  for (int i=0; i<nEvents ; ++i)
  {
     tree->GetEvent(i);

     std::vector<JetInfo> jets;
     for (unsigned int j=0; j<jet_pt->size(); ++j)
     {
       JetInfo jet;
       jet.pT = jet_pt->at(j);
       jet.eta = jet_eta->at(j);
       jet.phi = jet_phi->at(j);
       jet.mass = jet_mass->at(j);
       jet.btag = jet_btag->at(j);
       jets.push_back(jet);
     }

     // Now sorting this vector of structs
     std::sort (jets.begin(), jets.end(), sortJetsInDescendingpT);   

     std::vector<PhotonInfo> photons;
     for (unsigned int j=0; j<ph_pt->size(); ++j)
     {
       PhotonInfo photon;
       photon.pT=ph_pt->at(j);
       photon.eta=ph_eta->at(j);
       photon.phi=ph_phi->at(j);
       photons.push_back(photon);
     }
     // Now sorting this vector of structs
     std::sort (photons.begin(), photons.end(), sortPhotonsInDescendingpT);  

     std::vector<LeptonInfo> electrons;
     for (unsigned int j=0; j<el_pt->size(); ++j)
     {
        LeptonInfo electron;
        electron.pT=el_pt->at(j);
        electron.eta=el_eta->at(j);
        electron.phi=el_phi->at(j);
        electron.charge=el_charge->at(j);
        electrons.push_back(electron);
     }

     // Now sorting this vector of structs
     std::sort (electrons.begin(), electrons.end(), sortLeptonsInDescendingpT); 

     std::vector<LeptonInfo> muons;
     for (unsigned int j=0; j<mu_pt->size(); ++j)
     {
       LeptonInfo muon;
       muon.pT=mu_pt->at(j);
       muon.eta=mu_eta->at(j);
       muon.phi=mu_phi->at(j);
       muon.charge=mu_charge->at(j);
       muons.push_back(muon);
     }
     // Now sorting this vector of structs
     std::sort (muons.begin(), muons.end(), sortLeptonsInDescendingpT);

     double HT = 0.0; //jets are sorted. Don't care as far as HT is concerned.
     vector<TLorentzVector> Jet_vector;
     Jet_vector.clear();
     for(unsigned int k=0; k<jets.size(); ++k)
     {
       TLorentzVector Jet;
       if(fabs(jets.at(k).eta)<5.2 and jets.at(k).pT>50.0)
       {
         Jet.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);

         bool isGoodJet=true;
         for(unsigned int j=0; j<electrons.size(); ++j)
         {
           TLorentzVector Electron;
           Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
           double DRjet_el = Jet.DeltaR(Electron);
           if(DRjet_el<0.5) isGoodJet=false;
         }

         for(unsigned int j=0; j<muons.size(); ++j)
         {
           TLorentzVector Muon;
           Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
           double DRjet_mu = Jet.DeltaR(Muon);
           if(DRjet_mu<0.5) isGoodJet=false;
         }
         for(unsigned int l=0; l<photons.size(); ++l)
         {
           TLorentzVector Photon;
           Photon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Photon.SetPtEtaPhiM(photons.at(l).pT, photons.at(l).eta, photons.at(l).phi, 0.0);
           double DRjet_ph = Jet.DeltaR(Photon);
           if(DRjet_ph<0.5) isGoodJet=false;
         }
         if(isGoodJet) Jet_vector.push_back(Jet);
       }//close four vector if
     }

   //Now sorting this vector of structs
   std::sort (Jet_vector.begin(), Jet_vector.end(), sortJetVectorsInDescendingpT);


   //if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
   if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).Pt());

   if(Jet_vector.size()>0 ) h_jet_eta_leading->Fill(Jet_vector.at(0).Eta());
   if(Jet_vector.size()>1 ) h_jet_eta_trailing->Fill(Jet_vector.at(1).Eta());

   if(Jet_vector.size()>0 ) h_jet_phi_leading->Fill(Jet_vector.at(0).Phi());
   if(Jet_vector.size()>1 ) h_jet_phi_trailing->Fill(Jet_vector.at(1).Phi());

   h_nJets->Fill(Jet_vector.size());

   TLorentzVector Jet1;
   TLorentzVector Jet2;

   if(Jet_vector.size()>0 ) Jet1.SetPtEtaPhiE(Jet_vector.at(0).Pt(), Jet_vector.at(0).Eta(), Jet_vector.at(0).Phi(), Jet_vector.at(0).E());
   else Jet1.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);

   if(Jet_vector.size()>1 ) Jet2.SetPtEtaPhiE(Jet_vector.at(1).Pt(), Jet_vector.at(1).Eta(), Jet_vector.at(1).Phi(), Jet_vector.at(1).E());
   else Jet2.SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
   
   int nleptons = electrons.size() + muons.size();

   if(electrons.size()>=2 and electrons.at(0).charge*electrons.at(1).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_WJ_PT_ElEl->Fill((Jet1+Jet2).Pt());
   if(muons.size()>=2 and muons.at(0).charge*muons.at(1).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_WJ_PT_MuMu->Fill((Jet1+Jet2).Pt()); 
   if(muons.size()>=1 and electrons.size()>=1 and muons.at(0).charge*electrons.at(0).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_WJ_PT_ElMu->Fill((Jet1+Jet2).Pt());   

   if(electrons.size()>=1) h_Electron0_Charge->Fill(electrons.at(0).charge);
   if(electrons.size()>=2) h_Electron1_Charge->Fill(electrons.at(1).charge);

   if(muons.size()>=1) h_Muon0_Charge->Fill(muons.at(0).charge);
   if(muons.size()>=2) h_Muon1_Charge->Fill(muons.at(1).charge);

   if(Jet_vector.size()>0 and electrons.size()>=2 and electrons.at(0).charge*electrons.at(1).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
   if(Jet_vector.size()>0 and muons.size()>=2 and muons.at(0).charge*muons.at(1).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
   if(Jet_vector.size()>0 and muons.size()>=1 and electrons.size()>=1 and muons.at(0).charge*electrons.at(0).charge==1 and (Jet1+Jet2).Pt() > 50.0) h_jet_pt_leading->Fill(Jet_vector.at(0).Pt());
  }//event loop closed

  std::cout << "n_elec = " << n_elec << std::endl;
  std::cout << "n_muon = " << n_muon << std::endl;
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_WJ_PT_ElEl->Write();
  h_WJ_PT_MuMu->Write();
  h_WJ_PT_ElMu->Write();
  h_Electron0_Charge->Write();
  h_Electron1_Charge->Write();
  h_Muon0_Charge->Write();
  h_Muon1_Charge->Write();
  h_nJets->Write();
  h_jet_pt_leading->Write();
  h_jet_pt_trailing->Write();

  h_jet_eta_leading->Write();
  h_jet_eta_trailing->Write();

  h_jet_phi_leading->Write();
  h_jet_phi_trailing->Write();

  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;

}
