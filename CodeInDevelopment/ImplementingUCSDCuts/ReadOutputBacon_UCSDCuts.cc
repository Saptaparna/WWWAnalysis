#include "ReadOutputBacon_UCSDCuts.h"


int ReadOutputBacon_UCSDCuts(std::string infile, std::string treeStr, std::string outfile)
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
  Float_t         eventWeightPU;
  Float_t         puWeight;
  UInt_t          nPU;
  Float_t         met;
  Float_t         metPhi;
  vector<float>   *muon_pt;
  vector<float>   *muon_phi;
  vector<float>   *muon_eta;
  vector<float>   *muon_trkIso;
  vector<float>   *muon_pfIso;
  vector<float>   *muon_charge;
  vector<bool>    *muon_id;
  vector<bool>    *muon_id_alternate;
  vector<float>   *muon_recoEW;
  vector<float>   *muon_triggerEW;
  vector<float>   *muon_trigger;
  vector<float>   *electron_pt;
  vector<float>   *electron_phi;
  vector<float>   *electron_eta;
  vector<float>   *electron_pfIso;
  vector<float>   *electron_trkIso;
  vector<float>   *electron_charge;
  vector<bool>    *electron_id;
  vector<float>   *electron_recoEW;
  vector<float>   *electron_triggerEW;
  vector<float>   *electron_trigger;
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
  muon_trkIso = 0;
  muon_pfIso = 0;
  muon_charge = 0;
  muon_id = 0;
  muon_id_alternate = 0;
  muon_recoEW = 0;
  muon_triggerEW = 0;
  muon_trigger = 0;
  electron_pt = 0;
  electron_phi = 0;
  electron_eta = 0;
  electron_pfIso = 0;
  electron_trkIso = 0;
  electron_charge = 0;
  electron_id = 0;
  electron_recoEW = 0;
  electron_triggerEW = 0;
  electron_trigger = 0;
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
  tree->SetBranchAddress("muon_trkIso", &(muon_trkIso));
  tree->SetBranchAddress("muon_pfIso", &(muon_pfIso));
  tree->SetBranchAddress("muon_charge", &(muon_charge));
  tree->SetBranchAddress("muon_id", &(muon_id));
  tree->SetBranchAddress("muon_id_alternate", &(muon_id_alternate));
  tree->SetBranchAddress("muon_recoEW", &(muon_recoEW));
  tree->SetBranchAddress("muon_triggerEW", &(muon_triggerEW));
  tree->SetBranchAddress("muon_trigger", &(muon_trigger));
  tree->SetBranchAddress("electron_pt", &(electron_pt));
  tree->SetBranchAddress("electron_phi", &(electron_phi));
  tree->SetBranchAddress("electron_eta", &(electron_eta));
  tree->SetBranchAddress("electron_pfIso", &(electron_pfIso));
  tree->SetBranchAddress("electron_trkIso", &(electron_trkIso));
  tree->SetBranchAddress("electron_charge", &(electron_charge));
  tree->SetBranchAddress("electron_id", &(electron_id));
  tree->SetBranchAddress("electron_recoEW", &(electron_recoEW));
  tree->SetBranchAddress("electron_triggerEW", &(electron_triggerEW));
  tree->SetBranchAddress("electron_trigger", &(electron_trigger));
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
  
  int nEvents=tree->GetEntries();
  std::cout << "nEvents = " << nEvents << std::endl;
  double cut1_mumu, cut2_mumu, cut3_mumu, cut4_mumu, cut5_mumu, cut6_mumu, cut7_mumu, cut8_mumu;
  cut1_mumu = cut2_mumu = cut3_mumu = cut4_mumu = cut5_mumu = cut6_mumu = cut7_mumu = cut8_mumu = 0.0;
  double cut1_elel, cut2_elel, cut3_elel, cut4_elel, cut5_elel, cut6_elel, cut7_elel, cut8_elel;
  cut1_elel = cut2_elel = cut3_elel = cut4_elel = cut5_elel = cut6_elel = cut7_elel = cut8_elel = 0.0;
  HistCollection mumuHistCut1;
  initializeHistCollection(mumuHistCut1, "2SSTL_MuMu");
  HistCollection mumuHistCut2;
  initializeHistCollection(mumuHistCut2, "2SSTLLV_MuMu");
  HistCollection mumuHistCut3;
  initializeHistCollection(mumuHistCut3, "2SSTLLV2J_MuMu");
  HistCollection mumuHistCut4;
  initializeHistCollection(mumuHistCut4, "2SSTLLV2JBV_MuMu");
  HistCollection mumuHistCut5;
  initializeHistCollection(mumuHistCut5, "2SSTLLV2JBVCMjj_MuMu");
  HistCollection mumuHistCut6;
  initializeHistCollection(mumuHistCut6, "2SSTLLV2JBVCMjjLMjj_MuMu");
  HistCollection mumuHistCut7;
  initializeHistCollection(mumuHistCut7, "2SSTLLV2JBVCMjjLMjjMet_MuMu");
  HistCollection mumuHistCut8;
  initializeHistCollection(mumuHistCut8, "2SSTLLV2JBVCMjjLMjjMetMll_MuMu");
  HistCollection mumuHistCut9;
  initializeHistCollection(mumuHistCut9, "TTbarCR_MuMu");
  HistCollection mumuHistCut10;
  initializeHistCollection(mumuHistCut10, "TTbarCR1b_MuMu");
  HistCollection mumuHistCut11;
  initializeHistCollection(mumuHistCut11, "TTbarCR2b_MuMu");

  HistCollection elelHistCut1;
  initializeHistCollection(elelHistCut1, "2SSTL_ElEl");
  HistCollection elelHistCut2;
  initializeHistCollection(elelHistCut2, "2SSTLLV_ElEl");
  HistCollection elelHistCut3;
  initializeHistCollection(elelHistCut3, "2SSTLLV2J_ElEl");
  HistCollection elelHistCut4;
  initializeHistCollection(elelHistCut4, "2SSTLLV2JBV_ElEl");
  HistCollection elelHistCut5;
  initializeHistCollection(elelHistCut5, "2SSTLLV2JBVCMjj_ElEl");
  HistCollection elelHistCut6;
  initializeHistCollection(elelHistCut6, "2SSTLLV2JBVCMjjLMjj_ElEl");
  HistCollection elelHistCut7;
  initializeHistCollection(elelHistCut7, "2SSTLLV2JBVCMjjLMjjMet_ElEl");
  HistCollection elelHistCut8;
  initializeHistCollection(elelHistCut8, "2SSTLLV2JBVCMjjLMjjMetMll_ElEl");
  HistCollection elelHistCut9;
  initializeHistCollection(elelHistCut9, "TTbarCR_ElEl");
  HistCollection elelHistCut10;
  initializeHistCollection(elelHistCut10, "TTbarCR1b_ElEl");
  HistCollection elelHistCut11;
  initializeHistCollection(elelHistCut11, "TTbarCR2b_ElEl");

  TH1D *h_TotalEvents = new TH1D("h_TotalEvents", "h_TotalEvents", 15, 0.5, 15.5);
  //nEvents = 1000000;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);

    std::vector<leptonInfo> v_muons;
    std::vector<leptonInfo> v_looseMuons;
    for (unsigned int imuon=0; imuon<muon_pt->size(); imuon++)
    {
      leptonInfo muon;
      muon.pt = muon_pt->at(imuon); 
      muon.phi = muon_phi->at(imuon);
      muon.eta = muon_eta->at(imuon);
      muon.pfiso = muon_pfIso->at(imuon);
      muon.trkiso = muon_trkIso->at(imuon);
      muon.charge = (int)muon_charge->at(imuon);
      muon.id = muon_id->at(imuon);
      muon.recoEW = muon_recoEW->at(imuon);
      muon.trigger = muon_trigger->at(imuon);
      muon.triggerEW = muon_triggerEW->at(imuon);
      if(muon.id==1 and muon.trkiso/muon.pt < 0.1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4)
      {
        v_muons.push_back(muon);
      }
      if(muon.pt > 10.0 and muon.id==0) v_looseMuons.push_back(muon);
    }//muon loop

    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);
    std::sort (v_looseMuons.begin(), v_looseMuons.end(), sortLeptonsInDescendingpT);
     
    std::vector<leptonInfo> v_electrons;
    std::vector<leptonInfo> v_looseElectrons;
    for (unsigned int ielectron=0; ielectron<electron_pt->size(); ielectron++)
    {
      leptonInfo electron;
      electron.pt = electron_pt->at(ielectron);
      electron.phi = electron_phi->at(ielectron);
      electron.eta = electron_eta->at(ielectron);
      electron.pfiso = electron_pfIso->at(ielectron);
      electron.trkiso = electron_trkIso->at(ielectron);
      electron.charge = (int)electron_charge->at(ielectron);
      electron.id = electron_id->at(ielectron);
      electron.recoEW = electron_recoEW->at(ielectron);
      electron.trigger = electron_trigger->at(ielectron);
      electron.triggerEW = electron_triggerEW->at(ielectron);
      if(electron.id==1 and electron.trkiso/electron.pt < 0.0588 and electron.pt > 30.0 and fabs(electron.eta) < 2.4)
      {
        v_electrons.push_back(electron);
      }
      if(electron.pt > 10.0 and electron.id==0) v_looseElectrons.push_back(electron);
    }//electron loop

    std::sort (v_electrons.begin(), v_electrons.end(), sortLeptonsInDescendingpT);
    std::sort (v_looseElectrons.begin(), v_looseElectrons.end(), sortLeptonsInDescendingpT);
 
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

    std::vector<AnalysisJetInfo> v_selectedJets;
    for(unsigned int iselJet=0; iselJet<v_jets.size(); ++iselJet)
    {
      AnalysisJetInfo Jet;
  
      if(fabs(v_jets.at(iselJet).jetEta)<3.0 and v_jets.at(iselJet).jetPt>30.0)
      {
        bool isGoodJet=true;
        Jet.JetLV.SetPtEtaPhiM(v_jets.at(iselJet).jetPt, v_jets.at(iselJet).jetEta, v_jets.at(iselJet).jetPhi, v_jets.at(iselJet).jetMass);
        Jet.BTag_CSV = v_jets.at(iselJet).jetCSV;
        //std::cout << "Jet.BTag_CSV = " << Jet.BTag_CSV << std::endl;
        for(unsigned int iselMuon=0; iselMuon<v_muons.size(); ++iselMuon)
        {
          TLorentzVector Muon = fillTLorentzVector(v_muons.at(iselMuon).pt, v_muons.at(iselMuon).eta, v_muons.at(iselMuon).phi, MUON_MASS);
          double DRjet_mu = Jet.JetLV.DeltaR(Muon);
          if(DRjet_mu<0.5) isGoodJet=false;
        }//muon loop closed
        for(unsigned int iselElectron=0; iselElectron<v_electrons.size(); ++iselElectron)
        {
          TLorentzVector Electron = fillTLorentzVector(v_electrons.at(iselElectron).pt, v_electrons.at(iselElectron).eta, v_electrons.at(iselElectron).phi, ELECTRON_MASS);
          double DRjet_el = Jet.JetLV.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodJet=false;
        }//electron loop closed
        if(isGoodJet) v_selectedJets.push_back(Jet);
      //if(isGoodJet and Jet.BTag_CSV > 0.460) nbJets++;
      }//jet 4 vector closed
    }//selected jet loop closed

    std::sort (v_selectedJets.begin(), v_selectedJets.end(), sortJetVectorsInDescendingpT);

    int nbJets_loose = 0;
    int nbJets_tight  = 0;
    std::vector<AnalysisJetInfo> v_selectedBJets;
    for(unsigned int iselbJet=0; iselbJet<v_selectedJets.size(); ++iselbJet)
    {
      AnalysisJetInfo bJet;
      if(v_selectedJets.at(iselbJet).BTag_CSV > 0.460) 
      {
        nbJets_loose++;
        bJet.JetLV.SetPtEtaPhiM(v_selectedJets.at(iselbJet).JetLV.Pt(), v_selectedJets.at(iselbJet).JetLV.Eta(), v_selectedJets.at(iselbJet).JetLV.Phi(), v_selectedJets.at(iselbJet).JetLV.M()); 
        v_selectedBJets.push_back(bJet);
      }
      if(v_selectedJets.at(iselbJet).BTag_CSV > 0.8484) nbJets_tight++;
    }
    if(nbJets_loose != (int)v_selectedBJets.size())
    {
      std::cout << "nbJets_loose = " << nbJets_loose << std::endl;
      std::cout << "v_selectedBJets.size() = " << v_selectedBJets.size() << std::endl;
    }
    double mindR = 0.8;
    //double mindR = 6.77; //sqrt(3.14^2 + 6^2)
    int i_Jet1 = -1;
    int i_Jet2 = -1;
    for(unsigned int iJet=0; iJet<v_selectedJets.size(); iJet++)
    {
      AnalysisJetInfo Jet1;
      Jet1.JetLV.SetPtEtaPhiM(v_selectedJets.at(iJet).JetLV.Pt(), v_selectedJets.at(iJet).JetLV.Eta(), v_selectedJets.at(iJet).JetLV.Phi(), v_selectedJets.at(iJet).JetLV.M());
      for(unsigned int iJet2=iJet+1; iJet2<v_selectedJets.size(); iJet2++)
      {
        AnalysisJetInfo Jet2;
        Jet2.JetLV.SetPtEtaPhiM(v_selectedJets.at(iJet2).JetLV.Pt(), v_selectedJets.at(iJet2).JetLV.Eta(), v_selectedJets.at(iJet2).JetLV.Phi(), v_selectedJets.at(iJet2).JetLV.M());
        double deltaR = Jet1.JetLV.DeltaR(Jet2.JetLV);
        if(deltaR < mindR)
        {
          mindR = deltaR;
          i_Jet1 = iJet;
          i_Jet2 = iJet2;
        }
      }
    }
    double invMassJJ = 0;
    if(i_Jet1 >= 0 and i_Jet2 >= 0)
    {
      TLorentzVector closeJet1 = fillTLorentzVector(v_selectedJets.at(i_Jet1).JetLV.Pt(), v_selectedJets.at(i_Jet1).JetLV.Eta(), v_selectedJets.at(i_Jet1).JetLV.Phi(), v_selectedJets.at(i_Jet1).JetLV.M());
      TLorentzVector closeJet2 = fillTLorentzVector(v_selectedJets.at(i_Jet2).JetLV.Pt(), v_selectedJets.at(i_Jet2).JetLV.Eta(), v_selectedJets.at(i_Jet2).JetLV.Phi(), v_selectedJets.at(i_Jet2).JetLV.M());
      invMassJJ = (closeJet1+closeJet2).M();
    }
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
    
    //if(v_electrons.size() > 0) std::cout << "v_electrons.size() = " << v_electrons.size() << std::endl;
    //if(v_muons.size() > 0) std::cout << "v_muons.size() = " << v_muons.size() << std::endl;

    if(v_muons.size() == 2 and v_electrons.size()==0)
    {
      double eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, jet1csv, jet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = jet1csv = jet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muons.at(0).charge*v_muons.at(1).charge==1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        fillMuHistCollection(mumuHistCut1, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);  
        cut1_mumu+=eventWeightMu;
        h_TotalEvents->Fill(1);
        if(v_looseMuons.size()==0) 
        {
          fillMuHistCollection(mumuHistCut2, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);  
          cut2_mumu+=eventWeightMu;
          h_TotalEvents->Fill(2);
          if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
          {
            cut3_mumu+=eventWeightMu;
            if(v_selectedBJets.size()!=0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
            if(v_selectedBJets.size()!=0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
            if(v_selectedBJets.size()!=0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi(); 
            fillMuHistCollection(mumuHistCut3, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
            h_TotalEvents->Fill(3);
            if(nbJets_loose!=0) continue;
            
            if(fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5 and v_selectedJets.at(0).BTag_CSV > 0.0 and v_selectedJets.at(0).BTag_CSV < 0.460 and v_selectedJets.at(1).BTag_CSV > 0.0 and v_selectedJets.at(1).BTag_CSV < 0.460)
            {
              cut4_mumu+=eventWeightMu;
              fillMuHistCollection(mumuHistCut4, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
              double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
              double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
              h_TotalEvents->Fill(4);
              if(invMassJJ > 60 and invMassJJ < 100) 
              {
                cut5_mumu+=eventWeightMu;
                fillMuHistCollection(mumuHistCut5, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                h_TotalEvents->Fill(5);
                if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                {
                  cut6_mumu+=eventWeightMu;
                  fillMuHistCollection(mumuHistCut6, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);   
                  h_TotalEvents->Fill(6);
                  if(met > 40.0) 
                  {
                    cut7_mumu+=eventWeightMu;
                    fillMuHistCollection(mumuHistCut7, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                    h_TotalEvents->Fill(7);
                    if((mu1+mu2).M() > 40.0) 
                    {
                      cut8_mumu+=eventWeightMu;
                      fillMuHistCollection(mumuHistCut8, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                      h_TotalEvents->Fill(8);
                    }//Mll > 40.0
                  }//met > 40.0 
                }//Leading invmass and delta eta
              }//closest invmass
            }//b-veto
          }//2 jets with pT > 30
        }//loose lepton veto
      }//exactly 2 muons with pT > 30.0
    }//exactly 2 muons
    else if(v_electrons.size() == 2 and v_muons.size()==0)
    {
      double eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, jet1csv, jet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = jet1csv = jet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_electrons.at(0).charge*v_electrons.at(1).charge==1 and v_electrons.at(0).pt > 30.0 and v_electrons.at(1).pt > 30.0)
      {
        TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
        TLorentzVector el2 = fillTLorentzVector(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, ELECTRON_MASS);
        fillElHistCollection(elelHistCut1, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
        cut1_elel+=eventWeightEl;
        h_TotalEvents->Fill(1);
        if(v_looseElectrons.size()==0)
        {
          fillElHistCollection(elelHistCut2, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
          cut2_elel+=eventWeightEl;
          h_TotalEvents->Fill(2);
          if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
          {
            cut3_elel+=eventWeightEl;
            if(v_selectedBJets.size()!=0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
            if(v_selectedBJets.size()!=0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
            if(v_selectedBJets.size()!=0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
            fillElHistCollection(elelHistCut3, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
            h_TotalEvents->Fill(3);
            if(nbJets_loose!=0) continue;
            if(v_looseElectrons.size()==0 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5 and v_selectedJets.at(0).BTag_CSV > 0.0 and v_selectedJets.at(0).BTag_CSV < 0.460 and v_selectedJets.at(1).BTag_CSV > 0.0 and v_selectedJets.at(1).BTag_CSV < 0.460)
            {
              cut4_elel+=eventWeightEl;
              fillElHistCollection(elelHistCut4, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
              double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
              double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
              h_TotalEvents->Fill(4);
              if(invMassJJ > 60 and invMassJJ < 100)
              {
                cut5_elel+=eventWeightEl;
                fillElHistCollection(elelHistCut5, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                h_TotalEvents->Fill(5);
                if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                {
                  cut6_elel+=eventWeightEl;
                  fillElHistCollection(elelHistCut6, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                  h_TotalEvents->Fill(6);
                  if(met > 40.0)
                  {
                    cut7_elel+=eventWeightEl;
                    fillElHistCollection(elelHistCut7, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                    h_TotalEvents->Fill(7);
                    if((el1+el2).M() > 40.0)
                    {
                      cut8_elel+=eventWeightEl;
                      fillElHistCollection(elelHistCut8, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                       h_TotalEvents->Fill(8);
                    }//Mll > 40.0
                  }//met > 40.0 
                }//Leading invmass and delta eta
              }//closest invmass
            }//b-veto
          }//2 jets with pT > 30
        }//loose lepton veto
      }//exactly 2 electrons with pT > 30.0
    }//exactly 2 electrons

  //using bools instead of nested if
  
    bool mut = (v_muons.size()==2 and v_muons.at(0).charge*v_muons.at(1).charge==1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0);
    bool mutlv = mut and v_looseMuons.size()==0;
    bool mutlv2j = mutlv and (v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5); 
    bool mutlv2jbV = mutlv2j and (nbJets_loose==0);
/*
  if(2mut)
  if(2mut) h_MET_mumu->Fill(met); 
  if(2mutlv) h_MET_mumu->Fill(met); 
*/

    //ttbar control region
    if(v_muons.size() >= 2 and v_electrons.size()==0)
    {
      double eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, jet1csv, jet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = jet1csv = jet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muons.at(0).charge*v_muons.at(1).charge==-1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0 and met > 30)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        fillMuHistCollection(mumuHistCut9, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
        if(nbJets_loose > 0 and met > 50 and v_selectedJets.size()>=2)
        {
          fillMuHistCollection(mumuHistCut10, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
          if(nbJets_loose > 1)
          {
            fillMuHistCollection(mumuHistCut11, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), v_selectedBJets.at(1).JetLV.Pt(), v_selectedBJets.at(1).JetLV.Eta(), v_selectedBJets.at(1).JetLV.Phi(), eventWeightMu);
          }//2b-jets  
        }//1b-jet
      }//OS muons
    }//2 muons
    else if(v_electrons.size() >= 2 and v_muons.size()==0)
    {
      double eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, jet1csv, jet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = jet1csv = jet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_electrons.at(0).charge*v_electrons.at(1).charge==-1 and v_electrons.at(0).pt > 30.0 and v_electrons.at(1).pt > 30.0 and met > 30)
      {
        TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
        TLorentzVector el2 = fillTLorentzVector(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, ELECTRON_MASS);
        fillElHistCollection(elelHistCut9, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, jet1csv, jet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
        if(nbJets_loose > 0 and met > 50 and v_selectedJets.size()>=2)
        {
          fillElHistCollection(elelHistCut10, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
          if(nbJets_loose > 1)
          {
            fillElHistCollection(elelHistCut11, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), nbJets_loose, v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), v_selectedJets.at(0).BTag_CSV, v_selectedJets.at(1).BTag_CSV, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), v_selectedBJets.at(1).JetLV.Pt(), v_selectedBJets.at(1).JetLV.Eta(), v_selectedBJets.at(1).JetLV.Phi(), eventWeightEl);
          }//2b-jets
        }//1b-jet
      }//OS electrons
    }//2 electrons
  }//end of event loop

  std::cout << "2SSMuons = " << cut1_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto = " << cut2_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets = " << cut3_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto = " << cut4_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow = " << cut5_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow = " << cut6_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut = " << cut7_mumu << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut+MllCut = " << cut8_mumu << std::endl;
  
  std::cout << "2SSElectrons = " << cut1_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto = " << cut2_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets = " << cut3_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets+BVeto = " << cut4_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets+BVeto+ClosestInvMassWindow = " << cut5_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow = " << cut6_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut = " << cut7_elel << std::endl;
  std::cout << "2SSElectrons+LooseElectronVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut+MllCut = " << cut8_elel << std::endl;
  
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  tFile->mkdir("2SSTL");
  tFile->mkdir("2SSTLLV");
  tFile->mkdir("2SSTLLV2J");
  tFile->mkdir("2SSTLLV2JBV");
  tFile->mkdir("2SSTLLV2JBVCMjj");
  tFile->mkdir("2SSTLLV2JBVCMjjLMjj");
  tFile->mkdir("2SSTLLV2JBVCMjjLMjjMet");
  tFile->mkdir("2SSTLLV2JBVCMjjLMjjMetMll");
  tFile->mkdir("TTbarCR");
  tFile->mkdir("TTbarCR1b");
  tFile->mkdir("TTbarCR2b");
  tFile->cd("2SSTL");
  writeHistCollection(mumuHistCut1);
  writeHistCollection(elelHistCut1);
  tFile->cd("2SSTLLV");
  writeHistCollection(mumuHistCut2);
  writeHistCollection(elelHistCut2);
  tFile->cd("2SSTLLV2J");
  writeHistCollection(mumuHistCut3);
  writeHistCollection(elelHistCut3);
  tFile->cd("2SSTLLV2JBV");
  writeHistCollection(mumuHistCut4);
  writeHistCollection(elelHistCut4);
  tFile->cd("2SSTLLV2JBVCMjj");
  writeHistCollection(mumuHistCut5);
  writeHistCollection(elelHistCut5);
  tFile->cd("2SSTLLV2JBVCMjjLMjj");
  writeHistCollection(mumuHistCut6);
  writeHistCollection(elelHistCut6);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMet");
  writeHistCollection(mumuHistCut7);
  writeHistCollection(elelHistCut7);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMetMll");
  writeHistCollection(mumuHistCut8);
  writeHistCollection(elelHistCut8);
  tFile->cd("TTbarCR");
  writeHistCollection(mumuHistCut9);
  writeHistCollection(elelHistCut9);
  tFile->cd("TTbarCR1b");
  writeHistCollection(mumuHistCut10);
  writeHistCollection(elelHistCut10);
  tFile->cd("TTbarCR2b");
  writeHistCollection(mumuHistCut11);
  writeHistCollection(elelHistCut11);
  tFile->cd();
  h_TotalEvents->Write();
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl; 


  return 0;

}
