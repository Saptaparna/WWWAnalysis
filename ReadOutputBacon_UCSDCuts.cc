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
  UInt_t          nPU;
  Float_t         met;
  Float_t         metPhi;
  vector<float>   *muon_pt;
  vector<float>   *muon_phi;
  vector<float>   *muon_eta;
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
  
  int nEvents=tree->GetEntries();
  std::cout << "nEvents = " << nEvents << std::endl;
  double cut1, cut2, cut3, cut4, cut5, cut6, cut7, cut8;
  cut1 = cut2 = cut3 = cut4 = cut5 = cut6 = cut7 = cut8 = 0.0;
  HistCollection mumuHistCut1;
  initializeHistCollection(mumuHistCut1, "2SSTL");
  HistCollection mumuHistCut2;
  initializeHistCollection(mumuHistCut2, "2SSTLLV");
  HistCollection mumuHistCut3;
  initializeHistCollection(mumuHistCut3, "2SSTLLV2J");
  HistCollection mumuHistCut4;
  initializeHistCollection(mumuHistCut4, "2SSTLLV2JBV");
  HistCollection mumuHistCut5;
  initializeHistCollection(mumuHistCut5, "2SSTLLV2JBVCMjj");
  HistCollection mumuHistCut6;
  initializeHistCollection(mumuHistCut6, "2SSTLLV2JBVCMjjLMjj");
  HistCollection mumuHistCut7;
  initializeHistCollection(mumuHistCut7, "2SSTLLV2JBVCMjjLMjjMet");
  HistCollection mumuHistCut8;
  initializeHistCollection(mumuHistCut8, "2SSTLLV2JBVCMjjLMjjMetMll");
  //nEvents = 10000;
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
      muon.iso = muon_iso->at(imuon);
      muon.charge = (int)muon_charge->at(imuon);
      muon.id = muon_id->at(imuon);
      muon.recoEW = muon_recoEW->at(imuon);
      muon.trigger = muon_trigger->at(imuon);
      muon.triggerEW = muon_triggerEW->at(imuon);
      if(muon.id==1 and muon.iso/muon.pt < 0.1 and muon.pt > 25.0 and fabs(muon.eta) < 2.4)
      {
        v_muons.push_back(muon);
      }
        if(muon.pt > 10.0 and muon.id==0) v_looseMuons.push_back(muon);
    }//muon loop

    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);
    std::sort (v_looseMuons.begin(), v_looseMuons.end(), sortLeptonsInDescendingpT);
      
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
      if(isGoodJet) v_selectedJets.push_back(Jet);
      //if(isGoodJet and Jet.BTag_CSV > 0.460) nbJets++;
      }//jet 4 vector closed
    }//selected jet loop closed

    std::sort (v_selectedJets.begin(), v_selectedJets.end(), sortJetVectorsInDescendingpT);

    int nbJets_loose = 0;
    int nbJets_tight  = 0;
    for(unsigned int iselbJet=0; iselbJet<v_selectedJets.size(); ++iselbJet)
    {
      if(v_selectedJets.at(iselbJet).BTag_CSV > 0.460) nbJets_loose++;
      if(v_selectedJets.at(iselbJet).BTag_CSV > 0.8484) nbJets_tight++;
    }

    double mindR = 0.8;
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
    
    if(v_muons.size() == 2)
    {
      double eventWeight = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW));
      if(v_muons.at(0).charge*v_muons.at(1).charge==1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        mumuHistCut1.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
        mumuHistCut1.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
        mumuHistCut1.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
        mumuHistCut1.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
        mumuHistCut1.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
        mumuHistCut1.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
        mumuHistCut1.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
        mumuHistCut1.h_MET->Fill(met, eventWeight);
        mumuHistCut1.h_nJets->Fill(v_selectedJets.size(), eventWeight);
        mumuHistCut1.h_nbJets->Fill(nbJets_loose, eventWeight);
        cut1+=eventWeight;
        if(v_looseMuons.size()==0) 
        {
          mumuHistCut2.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
          mumuHistCut2.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
          mumuHistCut2.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
          mumuHistCut2.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
          mumuHistCut2.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
          mumuHistCut2.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
          mumuHistCut2.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
          mumuHistCut2.h_MET->Fill(met, eventWeight);
          mumuHistCut2.h_nJets->Fill(v_selectedJets.size(), eventWeight);
          mumuHistCut2.h_nbJets->Fill(nbJets_loose, eventWeight);
          cut2+=eventWeight;
          if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
          {
            cut3+=eventWeight;
            mumuHistCut3.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
            mumuHistCut3.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
            mumuHistCut3.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
            mumuHistCut3.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
            mumuHistCut3.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
            mumuHistCut3.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
            mumuHistCut3.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
            mumuHistCut3.h_MET->Fill(met, eventWeight);
            mumuHistCut3.h_nJets->Fill(v_selectedJets.size(), eventWeight);
            mumuHistCut3.h_nbJets->Fill(nbJets_loose, eventWeight);     
            mumuHistCut3.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
            mumuHistCut3.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight); 
            mumuHistCut3.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
            mumuHistCut3.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
            mumuHistCut3.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
            mumuHistCut3.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
            if(nbJets_loose!=0) continue;
            
            if(v_looseMuons.size()==0 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5 and v_selectedJets.at(0).BTag_CSV > 0.0 and v_selectedJets.at(0).BTag_CSV < 0.460 and v_selectedJets.at(1).BTag_CSV > 0.0 and v_selectedJets.at(1).BTag_CSV < 0.460)
            {
              cut4+=eventWeight;
              mumuHistCut4.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
              mumuHistCut4.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
              mumuHistCut4.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
              mumuHistCut4.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
              mumuHistCut4.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
              mumuHistCut4.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
              mumuHistCut4.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
              mumuHistCut4.h_MET->Fill(met, eventWeight);
              mumuHistCut4.h_nJets->Fill(v_selectedJets.size(), eventWeight);
              mumuHistCut4.h_nbJets->Fill(nbJets_loose, eventWeight);
              mumuHistCut4.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
              mumuHistCut4.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight);
              mumuHistCut4.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
              mumuHistCut4.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
              mumuHistCut4.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
              mumuHistCut4.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
              double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
              double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
              if(invMassJJ > 60 and invMassJJ < 100) 
              {
                cut5+=eventWeight;
                mumuHistCut5.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
                mumuHistCut5.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
                mumuHistCut5.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
                mumuHistCut5.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
                mumuHistCut5.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
                mumuHistCut5.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
                mumuHistCut5.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
                mumuHistCut5.h_MET->Fill(met, eventWeight);
                mumuHistCut5.h_nJets->Fill(v_selectedJets.size(), eventWeight);
                mumuHistCut5.h_nbJets->Fill(nbJets_loose, eventWeight);
                mumuHistCut5.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
                mumuHistCut5.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight);
                mumuHistCut5.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
                mumuHistCut5.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
                mumuHistCut5.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
                mumuHistCut5.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                {
                  cut6+=eventWeight;
                  mumuHistCut6.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
                  mumuHistCut6.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
                  mumuHistCut6.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
                  mumuHistCut6.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
                  mumuHistCut6.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
                  mumuHistCut6.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
                  mumuHistCut6.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
                  mumuHistCut6.h_MET->Fill(met, eventWeight);
                  mumuHistCut6.h_nJets->Fill(v_selectedJets.size(), eventWeight);
                  mumuHistCut6.h_nbJets->Fill(nbJets_loose, eventWeight);
                  mumuHistCut6.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
                  mumuHistCut6.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight);
                  mumuHistCut6.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
                  mumuHistCut6.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
                  mumuHistCut6.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
                  mumuHistCut6.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
                  if(met > 40.0) 
                  {
                    cut7+=eventWeight;
                    mumuHistCut7.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
                    mumuHistCut7.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
                    mumuHistCut7.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
                    mumuHistCut7.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
                    mumuHistCut7.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
                    mumuHistCut7.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
                    mumuHistCut7.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
                    mumuHistCut7.h_MET->Fill(met, eventWeight);
                    mumuHistCut7.h_nJets->Fill(v_selectedJets.size(), eventWeight);
                    mumuHistCut7.h_nbJets->Fill(nbJets_loose, eventWeight);
                    mumuHistCut7.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
                    mumuHistCut7.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight);
                    mumuHistCut7.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
                    mumuHistCut7.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
                    mumuHistCut7.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
                    mumuHistCut7.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
                    if((mu1+mu2).M() > 40.0) 
                    {
                      cut8+=eventWeight;
                      mumuHistCut8.h_mu_pt_leading->Fill(mu1.Pt(), eventWeight);
                      mumuHistCut8.h_mu_pt_trailing->Fill(mu2.Pt(), eventWeight);
                      mumuHistCut8.h_mu_phi_leading->Fill(mu1.Phi(), eventWeight);
                      mumuHistCut8.h_mu_phi_trailing->Fill(mu2.Phi(), eventWeight);
                      mumuHistCut8.h_mu_eta_leading->Fill(mu1.Eta(), eventWeight);
                      mumuHistCut8.h_mu_eta_trailing->Fill(mu2.Eta(), eventWeight);
                      mumuHistCut8.h_InvariantMass->Fill((mu1+mu2).M(), eventWeight);
                      mumuHistCut8.h_MET->Fill(met, eventWeight);
                      mumuHistCut8.h_nJets->Fill(v_selectedJets.size(), eventWeight);
                      mumuHistCut8.h_nbJets->Fill(nbJets_loose, eventWeight);
                      mumuHistCut8.h_jet_pt_leading->Fill(v_selectedJets.at(0).JetLV.Pt(), eventWeight);
                      mumuHistCut8.h_jet_pt_trailing->Fill(v_selectedJets.at(1).JetLV.Pt(), eventWeight);
                      mumuHistCut8.h_jet_eta_leading->Fill(v_selectedJets.at(0).JetLV.Eta(), eventWeight);
                      mumuHistCut8.h_jet_eta_trailing->Fill(v_selectedJets.at(1).JetLV.Eta(), eventWeight);
                      mumuHistCut8.h_jet_phi_leading->Fill(v_selectedJets.at(0).JetLV.Phi(), eventWeight);
                      mumuHistCut8.h_jet_phi_trailing->Fill(v_selectedJets.at(1).JetLV.Phi(), eventWeight);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }//end of event loop

  std::cout << "2SSMuons = " << cut1 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto = " << cut2 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets = " << cut3 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto = " << cut4 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow = " << cut5 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow = " << cut6 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut = " << cut7 << std::endl;
  std::cout << "2SSMuons+LooseMuonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut+MllCut = " << cut8 << std::endl;
    
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
  tFile->cd("2SSTL");
  writeHistCollection(mumuHistCut1);
  tFile->cd("2SSTLLV");
  writeHistCollection(mumuHistCut2);
  tFile->cd("2SSTLLV2J");
  writeHistCollection(mumuHistCut3);
  tFile->cd("2SSTLLV2JBV");
  writeHistCollection(mumuHistCut4);
  tFile->cd("2SSTLLV2JBVCMjj");
  writeHistCollection(mumuHistCut5);
  tFile->cd("2SSTLLV2JBVCMjjLMjj");
  writeHistCollection(mumuHistCut6);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMet");
  writeHistCollection(mumuHistCut7);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMetMll");
  writeHistCollection(mumuHistCut8);
  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl; 


  return 0;

}
