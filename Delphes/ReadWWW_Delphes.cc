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

typedef struct
{
  TLorentzVector JetLV;
  float BTag_CSV;
} AnalysisJetInfo; 

bool sortLeptonsInDescendingpT(LeptonInfo lep1, LeptonInfo lep2)
{
  return (lep1.pT > lep2.pT);
}

bool sortJetsInDescendingpT(JetInfo jet1, JetInfo jet2)
{
  return (jet1.pT > jet2.pT);
}

bool sortJetVectorsInDescendingpT(AnalysisJetInfo jet1, AnalysisJetInfo jet2)
{
  return (jet1.JetLV.Pt() > jet2.JetLV.Pt());
}
/*
bool sortJetVectorsInDescendingpT(TLorentzVector jet1, TLorentzVector jet2)
{
  return (jet1.Pt() > jet2.Pt());
}*/

bool sortPhotonsInDescendingpT(PhotonInfo pho1, PhotonInfo pho2)
{
  return (pho1.pT > pho2.pT);
}

double mT(TLorentzVector p4, float met, float met_phi)
{
  float phi1 = p4.Phi();
  float phi2 = met_phi;
  float Et1  = p4.Et();
  float Et2  = met;
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

int ReadHGGFiles_Delphes(std::string infile, std::string outfile){

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


  TH1F *h_MET=new TH1F("h_MET", "Missing ET; MET [GeV]; Events/GeV", 600, 0, 600); h_MET->Sumw2();
  TH1F *h_ST = new TH1F("h_ST", "ST (scalar sum of jet + lepton pT); S_T [GeV]; Events/GeV", 5000, 0, 5000.0);h_ST->Sumw2();
  TH1F *h_HTb = new TH1F("h_HTb", "HTb (scalar sum of b-jet pT); H_Tb [GeV]; Events/GeV", 5000, 0, 5000.0);h_HTb->Sumw2();
  TH1F *h_jet_pt_leading=new TH1F("h_jet_pt_leading", "Leading jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_leading->Sumw2();
  TH1F *h_jet_pt_trailing=new TH1F("h_jet_pt_trailing", "Trailing jet pT; pT [GeV]; Events/GeV", 10000, 0, 1000); h_jet_pt_trailing->Sumw2();
  TH1F *h_nJets = new TH1F("h_nJets", "Number of Jets; Number of Jets; Events", 20, -0.5, 19.5);h_nJets->Sumw2();
  TH1F *h_nbJets = new TH1F("h_nbJets", "Number of b-Jets; Number of b-Jets; Events", 20, -0.5, 19.5);h_nbJets->Sumw2();
  TH1F *h_InvariantMass_PhPh=new TH1F("h_InvariantMass_PhPh", "Di-photon invariant mass; m_{#gamma#gamma} [GeV]; Events/GeV", 9000, 0, 300); h_InvariantMass_PhPh->Sumw2();
  TH1F *h_InvariantMass_jj=new TH1F("h_InvariantMass_jj", "jj invariant mass; m_{jj} [GeV]; Events/GeV", 1500, 0, 150); h_InvariantMass_jj->Sumw2();
  TH1F *h_W_pT=new TH1F("h_W_pT", "W boson p_{T}; W p_{T} [GeV]; Events/GeV", 300, 0, 3000); h_W_pT->Sumw2();
  TH1F *h_ST_AllCuts_MuMu = new TH1F("h_ST_AllCuts_MuMu", "ST (scalar sum of jet + lepton pT); S_T [GeV]; Events/GeV", 5000, 0, 5000.0); h_ST_AllCuts_MuMu->Sumw2();
  TH1F *h_ST_AllCuts_ElMu = new TH1F("h_ST_AllCuts_ElMu", "ST (scalar sum of jet + lepton pT); S_T [GeV]; Events/GeV", 5000, 0, 5000.0); h_ST_AllCuts_ElMu->Sumw2();
  TH1F *h_ST_AllCuts_ElEl = new TH1F("h_ST_AllCuts_ElEl", "ST (scalar sum of jet + lepton pT); S_T [GeV]; Events/GeV", 5000, 0, 5000.0); h_ST_AllCuts_ElEl->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  //nEvents = 1;
  int nevents_MuMu = 0.0;
  int nevents_ElMu = 0.0;
  int nevents_ElEl = 0.0;
  for (int i=0; i<nEvents ; ++i)
  {
     tree->GetEvent(i);

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
   
     // filling the muon's properties into a vector of struct
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
   
     double ST = 0.0; //jets are sorted. Don't care as far as ST is concerned.
     double HTb = 0.0; 
     vector<AnalysisJetInfo> Jet_vector;
     Jet_vector.clear();
     vector<AnalysisJetInfo> bJet_vector;
     bJet_vector.clear();
     for(unsigned int k=0; k<jets.size(); ++k)
     {
       AnalysisJetInfo Jet;
       AnalysisJetInfo bJet;
       if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0)
       {
         Jet.JetLV.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);
         Jet.BTag_CSV = jets.at(k).btag;
         bool isGoodJet=true;
         for(unsigned int j=0; j<electrons.size(); ++j)
         {
           TLorentzVector Electron;
           Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
           double DRjet_el = Jet.JetLV.DeltaR(Electron);
           if(DRjet_el<0.5) isGoodJet=false;
         }
         for(unsigned int j=0; j<muons.size(); ++j)
         {
           TLorentzVector Muon;
           Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
           double DRjet_mu = Jet.JetLV.DeltaR(Muon);
           if(DRjet_mu<0.5) isGoodJet=false;
         }  
         if(isGoodJet) Jet_vector.push_back(Jet);
       } 
       if(fabs(jets.at(k).eta)<2.4 and jets.at(k).pT>25.0 and jets.at(k).btag>0)
       {
         bJet.JetLV.SetPtEtaPhiM(jets.at(k).pT, jets.at(k).eta, jets.at(k).phi, jets.at(k).mass);
         bJet.BTag_CSV = jets.at(k).btag;
         bool isGoodbJet=true;
         for(unsigned int j=0; j<electrons.size(); ++j)
         {
           TLorentzVector Electron;
           Electron.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Electron.SetPtEtaPhiM(electrons.at(j).pT, electrons.at(j).eta, electrons.at(j).phi, 0.0);
           double DRjet_el = bJet.JetLV.DeltaR(Electron);
           if(DRjet_el<0.5) isGoodbJet=false;
         }
         for(unsigned int j=0; j<muons.size(); ++j)
         {
           TLorentzVector Muon;
           Muon.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
           Muon.SetPtEtaPhiM(muons.at(j).pT, muons.at(j).eta, muons.at(j).phi, 0.0);
           double DRjet_mu = bJet.JetLV.DeltaR(Muon);
           if(DRjet_mu<0.5) isGoodbJet=false;
         } 
         if(isGoodbJet) bJet_vector.push_back(bJet);
      }//close four vector if
    }//close jet loop

  // Now sorting this vector of structs
  std::sort (Jet_vector.begin(), Jet_vector.end(), sortJetVectorsInDescendingpT);

  for(unsigned int m=0; m<Jet_vector.size(); m++)
  {
    ST += Jet_vector.at(m).JetLV.Pt();
    //std::cout << "ST = " << ST << std::endl;
    //std::cout << "Jet_vector.at(m).Pt() = " << Jet_vector.at(m).Pt() << std::endl;
  }

  for(unsigned int imu=0; imu<muons.size(); imu++)
  {
    ST += muons.at(imu).pT;
    //std::cout << "ST mu = " << ST << std::endl;
    //std::cout << "muons.at(imu).pT = " << muons.at(imu).pT << std::endl; 
  }

  for(unsigned int iel=0; iel<electrons.size(); iel++)
  {
    ST += electrons.at(iel).pT;
    //std::cout << "ST el = " << ST << std::endl;
    //std::cout << "electrons.at(iel).pT = " << electrons.at(iel).pT << std::endl;  
  }

  for(unsigned int m=0; m<bJet_vector.size(); m++)
  {
    HTb += bJet_vector.at(m).JetLV.Pt();
  }

  TLorentzVector jet1;
  jet1.SetPtEtaPhiE(0, 0, 0, 0);

  TLorentzVector jet2;
  jet2.SetPtEtaPhiE(0, 0, 0, 0);

  if(Jet_vector.size()>1)
  {
    jet1.SetPtEtaPhiE(Jet_vector.at(0).JetLV.Pt(), Jet_vector.at(0).JetLV.Eta(), Jet_vector.at(0).JetLV.Phi(), Jet_vector.at(0).JetLV.E());
    jet2.SetPtEtaPhiE(Jet_vector.at(1).JetLV.Pt(), Jet_vector.at(1).JetLV.Eta(), Jet_vector.at(1).JetLV.Phi(), Jet_vector.at(1).JetLV.E());
    h_W_pT->Fill((jet1+jet2).Pt());
  }

  if((jet1+jet2).M()>50.0 and (jet1+jet2).M()<120.0)  h_InvariantMass_jj->Fill((jet1+jet2).M());

  h_MET->Fill(MET);
 
  if(Jet_vector.size()>0 ) h_jet_pt_leading->Fill(Jet_vector.at(0).JetLV.Pt());
  if(Jet_vector.size()>1 ) h_jet_pt_trailing->Fill(Jet_vector.at(1).JetLV.Pt());

  h_nJets->Fill(Jet_vector.size());
  h_nbJets->Fill(bJet_vector.size());
  h_ST->Fill(ST);
  h_HTb->Fill(HTb);

  int nj = 0;
  int nj30 = 0;
  int nb = 0;
  float Mjj = 0;
  float MjjL = -999;
  float DetajjL = -999;
  float MjjVBF = -999;
  float MjjDR1 = -999;
  float DetajjVBF = -999;
  float tmpDR = 9999;
  float tmpDR_DR1 = 9999;
  float DRjj = -999;
  float DRjjDR1 = -999;
  TLorentzVector j0_p4;
  TLorentzVector j1_p4;
  TLorentzVector j0_p4_DR1;
  TLorentzVector j1_p4_DR1;
  for (unsigned int i = 0; i < Jet_vector.size(); ++i)
  {
    // nb jets
    if (not (Jet_vector.at(i).JetLV.Pt() > 20. ))
    continue;
    if (fabs(Jet_vector.at(i).JetLV.Eta()) < 2.4 && Jet_vector.at(i).BTag_CSV==1)
    nb++;
    // njets across all eta
    if (not(Jet_vector.at(i).JetLV.Pt() > 30. ))
    continue;
    if (fabs(Jet_vector.at(i).JetLV.Eta()) < 5.0)
    nj++;

    // njets for central for SS
    if (fabs(Jet_vector.at(i).JetLV.Eta()) < 2.5)
    nj30++;

    // Compute Mjj using the closest two jets
    for (unsigned int j = i + 1; j < Jet_vector.size(); ++j)
    {
       if (not ( Jet_vector.at(j).JetLV.Pt() > 30. ))
       continue;

       // central eta
       if (fabs(Jet_vector.at(i).JetLV.Eta()) < 2.5 and fabs(Jet_vector.at(j).JetLV.Eta()) < 2.5)
       {
         // Choose the closest two jets
         float this_dR = Jet_vector.at(i).JetLV.DeltaR(Jet_vector.at(j).JetLV);
         if (this_dR < tmpDR)
         {
           tmpDR = this_dR;
           Mjj = (Jet_vector.at(i).JetLV + Jet_vector.at(j).JetLV).M();
           DRjj = tmpDR;
           j0_p4 = Jet_vector.at(i).JetLV.Pt() > Jet_vector.at(j).JetLV.Pt() ? Jet_vector.at(i).JetLV : Jet_vector.at(j).JetLV;
           j1_p4 = Jet_vector.at(i).JetLV.Pt() > Jet_vector.at(j).JetLV.Pt() ? Jet_vector.at(j).JetLV : Jet_vector.at(i).JetLV;
         }

         // Choose the jets with angle closest to dR = 1
         if (abs(this_dR - 1.) < abs(tmpDR_DR1 - 1.))
         {
           tmpDR_DR1 = this_dR;
           MjjDR1 = (Jet_vector.at(i).JetLV + Jet_vector.at(j).JetLV).M();
           DRjjDR1 = tmpDR_DR1;
           j0_p4_DR1 = Jet_vector.at(i).JetLV.Pt() > Jet_vector.at(j).JetLV.Pt() ? Jet_vector.at(i).JetLV : Jet_vector.at(j).JetLV;
           j1_p4_DR1 = Jet_vector.at(i).JetLV.Pt() > Jet_vector.at(j).JetLV.Pt() ? Jet_vector.at(j).JetLV : Jet_vector.at(i).JetLV;
         }

         // If they were not set then set (this makes it choose two leading ones)
         if (MjjL < 0)
         {
           MjjL = (Jet_vector.at(i).JetLV + Jet_vector.at(j).JetLV).M();
           DetajjL = fabs(Jet_vector.at(i).JetLV.Eta() - Jet_vector.at(j).JetLV.Eta());
         }
       }
     }
   }
   //std::cout << "Mjj = " << Mjj << std::endl;
   //std::cout << "nb = " << nb << std::endl;

   if(muons.size()==2 and muons.at(0).pT > 25.0 and muons.at(1).pT > 25.0 and muons.at(0).charge*muons.at(1).charge==1)//exactly 2 muons
   {
     //std::cout << "dimuon events" << std::endl;
     TLorentzVector mu1, mu2;
     mu1.SetPtEtaPhiM(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0); 
     mu2.SetPtEtaPhiM(muons.at(1).pT, muons.at(1).eta, muons.at(1).phi, 0.0); 
     if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and (mu1+mu2).M() > 40) nevents_MuMu++;
     //if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and (mu1+mu2).M() > 40) h_ST_AllCuts_MuMu->Fill(ST);
     if(nj30>=2 and nb==0 and DetajjL < 1.5 and (mu1+mu2).M() > 40) h_ST_AllCuts_MuMu->Fill(ST);
   }
   //else if(electrons.size()==2)
   //else if(electrons.size()==2 and electrons.at(0).pT > 25.0 and electrons.at(1).pT > 25.0)
   else if(electrons.size()==2 and electrons.at(0).pT > 25.0 and electrons.at(1).pT > 25.0 and electrons.at(0).charge*electrons.at(1).charge==1)//exactly 2 electrons
   {
     //std::cout << electrons.at(0).pT << std::endl;
     //std::cout << electrons.at(1).pT << std::endl;
     TLorentzVector el1, el2;
     el1.SetPtEtaPhiM(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, 0.0);
     el2.SetPtEtaPhiM(electrons.at(1).pT, electrons.at(1).eta, electrons.at(1).phi, 0.0);
     if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MET > 60.0 and abs((el1+el2).M()-91.1876)>10.0 and (el1+el2).M() > 40) nevents_ElEl++;
     //if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MET > 60.0 and abs((el1+el2).M()-91.1876)>10.0 and (el1+el2).M() > 40) h_ST_AllCuts_ElEl->Fill(ST);
     if(nj30>=2 and nb==0 and DetajjL < 1.5 and MET > 60.0 and abs((el1+el2).M()-91.1876)>10.0 and (el1+el2).M() > 40) h_ST_AllCuts_ElEl->Fill(ST);
   }  
   else if(muons.size()==1 and electrons.size()==1 and electrons.at(0).pT > 25.0 and  muons.at(0).pT > 25.0)//exactly 1 electron and 1 muon
   {
     TLorentzVector el1, mu1;
     mu1.SetPtEtaPhiM(muons.at(0).pT, muons.at(0).eta, muons.at(0).phi, 0.0);
     el1.SetPtEtaPhiM(electrons.at(0).pT, electrons.at(0).eta, electrons.at(0).phi, 0.0);
     double MT0 = mT(el1, MET, MET_Phi);
     double MT1 = mT(mu1, MET, MET_Phi);
     double MTmax = TMath::Max(MT0, MT1);
     if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and (el1+mu1).M() > 30 and MET > 60.0 and MTmax > 90.0) nevents_ElMu++;    
     //if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and (el1+mu1).M() > 30 and MET > 60.0 and MTmax > 90.0) h_ST_AllCuts_ElMu->Fill(ST); 
     if(nj30>=2 and nb==0 and DetajjL < 1.5 and (el1+mu1).M() > 30 and MET > 60.0 and MTmax > 90.0) h_ST_AllCuts_ElMu->Fill(ST);
   }   
  

  }//event loop closed

  std::cout << "nevents_MuMu = " << nevents_MuMu << std::endl;
  std::cout << "nevents_ElEl = " << nevents_ElEl << std::endl;
  std::cout << "nevents_ElMu = " << nevents_ElMu << std::endl;
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_nJets->Write();
  h_HTb->Write();
  h_nbJets->Write();
  h_jet_pt_leading->Write();
  h_jet_pt_trailing->Write();
  h_ST->Write();
  h_MET->Write();
  h_InvariantMass_jj->Write();
  h_W_pT->Write();
  h_ST_AllCuts_MuMu->Write();
  h_ST_AllCuts_ElMu->Write();
  h_ST_AllCuts_ElEl->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;

}
