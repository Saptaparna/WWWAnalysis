#include "FatJetAnalyzer.h"

using namespace baconhep;
using namespace std;
double ELECTRON_MASS  = 0.000511;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

FatJetAnalyzer::FatJetAnalyzer() : BLTSelector()
{
}

FatJetAnalyzer::~FatJetAnalyzer()
{
}

void FatJetAnalyzer::Begin(TTree *tree)
{
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
    std::sregex_token_iterator(), std::back_inserter(options));

    params.reset(new Parameters());
    params->setup(options);
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    triggerNames.push_back("HLT_IsoMu22_v*");
    triggerNames.push_back("HLT_IsoTkMu22_v*");
    triggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    triggerNames.push_back("HLT_Ele23_WPLoose_Gsf_v*");

    weights.reset(new WeightUtils(params->period, params->selection, false));
    lumiMask = RunLumiRangeMap();
    if (true) 
    { 
      string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
      //string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt";
      lumiMask.AddJSONFile(jsonFileName);
    }
    
    muonCorr = new rochcor2016();
    
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("eventWeightPU", &eventWeightPU);
    outTree->Branch("eventWeightGen", &eventWeightGen);
    outTree->Branch("puWeight", &puWeight);
    outTree->Branch("nPartons", &nPartons);
    outTree->Branch("topPtWeight", &topPtWeight);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPV", &nPV);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    //leptons
    //muons
    outTree->Branch("muon_pt", &muon_pt);
    outTree->Branch("muon_phi", &muon_phi); 
    outTree->Branch("muon_eta", &muon_eta);
    outTree->Branch("muon_pt_corr", &muon_pt_corr);
    outTree->Branch("muon_phi_corr", &muon_phi_corr);
    outTree->Branch("muon_eta_corr", &muon_eta_corr);
    outTree->Branch("muon_trkIso", &muon_trkIso);
    outTree->Branch("muon_pfIso", &muon_pfIso);
    outTree->Branch("muon_charge", &muon_charge);
    outTree->Branch("muon_d0", &muon_d0);
    outTree->Branch("muon_dz", &muon_dz);
    outTree->Branch("muon_sip3d", &muon_sip3d);
    outTree->Branch("muon_id", &muon_id);
    outTree->Branch("muon_id_alternate", &muon_id_alternate);
    outTree->Branch("muon_id_tightUCSD", &muon_id_tightUCSD);
    outTree->Branch("muon_id_tightMIT", &muon_id_tightMIT);
    outTree->Branch("muon_id_veryLooseUCSD", &muon_id_veryLooseUCSD);
    outTree->Branch("muon_recoEW", &muon_recoEW);
    outTree->Branch("muon_triggerEW", &muon_triggerEW);
    outTree->Branch("muon_trigger", &muon_trigger);
    //electron
    outTree->Branch("electron_pt", &electron_pt);
    outTree->Branch("electron_phi", &electron_phi);
    outTree->Branch("electron_eta", &electron_eta);
    outTree->Branch("electron_pfIso", &electron_pfIso);
    outTree->Branch("electron_trkIso", &electron_trkIso);
    outTree->Branch("electron_charge", &electron_charge);
    outTree->Branch("electron_d0", &electron_d0);
    outTree->Branch("electron_dz", &electron_dz);
    outTree->Branch("electron_sip3d", &electron_sip3d);
    outTree->Branch("electron_id", &electron_id);
    outTree->Branch("electron_id_tightUCSD", &electron_id_tightUCSD);
    outTree->Branch("electron_id_tightMIT", &electron_id_tightMIT);
    outTree->Branch("electron_id_veryLooseUCSD", &electron_id_veryLooseUCSD);
    outTree->Branch("electron_id_HLTsafeMIT", &electron_id_HLTsafeMIT);
    outTree->Branch("electron_recoEW", &electron_recoEW);
    outTree->Branch("electron_triggerEW", &electron_triggerEW);
    outTree->Branch("electron_trigger", &electron_trigger);
    //jets
    outTree->Branch("jet_pt", &jet_pt);
    outTree->Branch("jet_eta", &jet_eta);
    outTree->Branch("jet_phi", &jet_phi);
    outTree->Branch("jet_mass", &jet_mass);
    outTree->Branch("jet_csv", &jet_csv);
    //fat jets
    outTree->Branch("ak8jet_prunMass", &ak8jet_prunMass);
    outTree->Branch("ak8jet_trimMass", &ak8jet_trimMass);
    outTree->Branch("ak8jet_sd0", &ak8jet_sd0);
    outTree->Branch("ak8jet_pt", &ak8jet_pt);
    outTree->Branch("ak8jet_phi", &ak8jet_phi);
    outTree->Branch("ak8jet_eta", &ak8jet_eta);
    outTree->Branch("ak8jet_tau1", &ak8jet_tau1);
    outTree->Branch("ak8jet_tau2", &ak8jet_tau2);

    // event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    ReportPostBegin();
}

Bool_t FatJetAnalyzer::Process(Long64_t entry)
{

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;
    //if(not(fInfo->evtNum==56 or fInfo->evtNum==80)) return kTRUE;

    //muons
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    muon_pt_corr.clear();
    muon_eta_corr.clear();
    muon_phi_corr.clear();
    muon_trkIso.clear();
    muon_pfIso.clear();
    muon_charge.clear();
    muon_d0.clear();
    muon_dz.clear();
    muon_sip3d.clear();
    muon_id.clear();
    muon_id_alternate.clear();
    muon_id_tightUCSD.clear();
    muon_id_tightMIT.clear();
    muon_id_veryLooseUCSD.clear();
    muon_recoEW.clear();
    muon_triggerEW.clear();
    muon_trigger.clear();
    //electrons
    electron_pt.clear();
    electron_eta.clear();
    electron_phi.clear();
    electron_trkIso.clear();
    electron_pfIso.clear();
    electron_charge.clear();
    electron_d0.clear();
    electron_dz.clear();
    electron_sip3d.clear();
    electron_id.clear();
    electron_id_tightUCSD.clear();
    electron_id_tightMIT.clear();
    electron_id_veryLooseUCSD.clear();
    electron_id_HLTsafeMIT.clear();
    electron_recoEW.clear();
    electron_triggerEW.clear();
    electron_trigger.clear();
    //jets
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_mass.clear();
    jet_csv.clear();
    ak8jet_prunMass.clear();
    ak8jet_trimMass.clear();
    ak8jet_sd0.clear();
    ak8jet_pt.clear();
    ak8jet_phi.clear();
    ak8jet_eta.clear();
    ak8jet_tau1.clear();
    ak8jet_tau2.clear();

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    if (isData) 
    {
      RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
    }
    hTotalEvents->Fill(2);

    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) 
    {
      passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }
    if (!passTrigger && isData)
      return kTRUE;
    hTotalEvents->Fill(3);

    eventWeight   = 1;
    eventWeightPU = 1.0;
    eventWeightGen = 1.0;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPV           = fPVArr->GetEntries();
    if (isData) 
    {
      eventWeight *= 1.;
      eventWeightPU = 1.0;
    }
    else
    {
      nPU = fInfo->nPUmean;
      puWeight = weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
      eventWeight *= puWeight;
      eventWeightGen = fGenEvtInfo->weight;
      //std::cout << "puWeight = " << puWeight << std::endl; 
    }

    float topSF = 1.;
    if (!isData) 
    {
      unsigned count = 0;
      for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) 
      {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
        // parton counting for jet-binned Drell-Yan samples
        if (particle->status == 23 
         && (abs(particle->pdgId) < 6 || particle->pdgId == 21) 
         && particle->parent != -2)
        {
          ++count;
        }
        // top pt reweighting: get the scale factor based on the top quark pt
        if (abs(particle->pdgId) == 6 && particle->status == 62) 
        {
           topSF *= exp(0.0615 - 0.0005*particle->pt);
        }
      }
      nPartons = count;
      topPtWeight = sqrt(topSF);
    }
    else
    {
      nPartons = 0.0;
      topPtWeight = 1.0;
    } 
 
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    vector<TMuon*> tmp_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
/*
        if (    muon->pt > 5 
                &&
                fabs(muon->eta) < 2.4
                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
           ) {*/
            tmp_muons.push_back(muon);
        //}
    }
    sort(tmp_muons.begin(), tmp_muons.end(), sort_by_higher_pt<TMuon>); 

    vector<TMuon*> muons;
    for (unsigned i = 0; i < tmp_muons.size(); i++) 
    //for (int i=0; i < fMuonArr->GetEntries(); i++) 
    {
      TMuon* muon = tmp_muons[i];
      //TMuon* muon = (TMuon*) fMuonArr->At(i);
      //assert(muon);
      TLorentzVector muonP4;
      copy_p4(tmp_muons[i], MUON_MASS, muonP4);
      //copy_p4(muon, MUON_MASS, muonP4);
      float qter = 1.;
      if (isData) 
      {
        muonCorr->momcor_data(muonP4, muon->q, 0, qter);
      } else 
      {
        muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
      }
      if (muonP4.Pt() > 25 and muon->trkIso/muonP4.Pt() < 0.1) 
      {
        muon_pt.push_back(muon->pt);
        muon_phi.push_back(muon->phi);
        muon_eta.push_back(muon->eta);
        muon_pt_corr.push_back(muonP4.Pt());
        muon_eta_corr.push_back(muonP4.Eta());
        muon_phi_corr.push_back(muonP4.Phi());
        muon_trkIso.push_back(muon->trkIso);
        muon_pfIso.push_back(muonIso(muon));
        muon_charge.push_back(muon->q);
        muon_d0.push_back(muon->d0);
        muon_dz.push_back(muon->dz);
        muon_sip3d.push_back(muon->sip3d);
        bool isTight = false;
        if(fabs(muon->eta) < 2.4 and (muon->typeBits & baconhep::kGlobal) and (muon->typeBits & baconhep::kGlobal) and muon->muNchi2 < 10. and muon->nMatchStn  > 1 and muon->nPixHits > 0 and fabs(muon->d0) < 0.2 and fabs(muon->dz) < 0.5 and muon->nTkLayers > 5 and muon->nValidHits > 0) isTight = true;
        else isTight = false;
        muon_id.push_back(isTight); 
        muon_id_alternate.push_back(isTightMuon(muon, muon->eta, muonP4.Pt())); 
        muon_id_tightUCSD.push_back(isTightMuonUCSD(muon, muon->eta, muonP4.Pt()));
        muon_id_tightMIT.push_back(isTightMuonMIT(muon, muon->eta, muonP4.Pt()));
        muon_id_veryLooseUCSD.push_back(isVeryLooseMuonUCSD(muon, muon->eta, muonP4.Pt()));
        TLorentzVector muon_p4;
        muon_p4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
        muon_recoEW.push_back(weights->GetMuonRecoEff(muon_p4));    
        pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu22_v*", muon_p4);
        muon_triggerEW.push_back(trigEff.first);
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) 
        {
          triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muon_trigger.push_back(triggered);
      }
    }

    std::vector<TLorentzVector> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) 
    {
      TElectron* electron = (TElectron*) fElectronArr->At(i);
      assert(electron);
      electron_pt.push_back(electron->pt);
      electron_phi.push_back(electron->phi);
      electron_eta.push_back(electron->eta);
      electron_trkIso.push_back(electron->trkIso);
      electron_pfIso.push_back(electronIso(electron, fInfo->rhoJet));//change it with the correct defn of iso
      electron_charge.push_back(electron->q);
      electron_d0.push_back(electron->d0);
      electron_dz.push_back(electron->dz);
      electron_sip3d.push_back(electron->sip3d);
      electron_id.push_back(isTightElectron(electron));
      electron_id_tightUCSD.push_back(isTightElectronUCSD(electron)); 
      electron_id_tightMIT.push_back(isTightElectronMIT(electron));
      electron_id_veryLooseUCSD.push_back(isVeryLooseElectronUCSD(electron));
      electron_id_HLTsafeMIT.push_back(isHLTsafeElectronMIT(electron));
      TLorentzVector electron_p4;
      electron_p4.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, ELECTRON_MASS);
      electron_recoEW.push_back(weights->GetMuonRecoEff(electron_p4));
      pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_Ele27_WPTight_Gsf_v*", electron_p4); 
      electron_triggerEW.push_back(trigEff.first); 
      bool triggered = false;
      for (unsigned i = 0; i < triggerNames.size(); ++i)
      {
        triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
      }
      electron_trigger.push_back(triggered);
    }
 
    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    nJets    = 0;

    //Fat jet info//
    TClonesArray* jetAK8Collection;
    jetAK8Collection = fAddAK8CHSArr;

    for (int i=0; i < jetAK8Collection->GetEntries(); i++)
    {   
      TAddJet* ak8Jet = (TAddJet*) jetAK8Collection->At(i);
      assert(ak8Jet);
      TJet *jet    = (TJet*)fAK4CHSArr->At(ak8Jet->index);
      if(jet->pt > 20 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
      {
        ak8jet_prunMass.push_back(ak8Jet->mass_prun);
        ak8jet_trimMass.push_back(ak8Jet->mass_trim);
        ak8jet_sd0.push_back(ak8Jet->mass_sd0);
        ak8jet_tau1.push_back(ak8Jet->tau1);
        ak8jet_tau2.push_back(ak8Jet->tau2);
        ak8jet_pt.push_back(jet->pt);
        ak8jet_phi.push_back(jet->phi);
        ak8jet_eta.push_back(jet->eta);
      }
    }
    //ak4 jet info//    
    for (int i=0; i < jetCollection->GetEntries(); i++) 
    {
      TJet* jet = (TJet*) jetCollection->At(i);
      assert(jet);
      if(jet->pt > 20 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
      {
        jet_pt.push_back(jet->pt);
        jet_eta.push_back(jet->eta);
        jet_phi.push_back(jet->phi);
        jet_mass.push_back(jet->mass);
        jet_csv.push_back(jet->csv);
      }  
    }
  
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    hTotalEvents->Fill(5);
    
    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}


void FatJetAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void FatJetAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void FatJetAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

Bool_t FatJetAnalyzer::isTightMuon(baconhep::TMuon* muon, float muEta, float muPt)
{
  if (muPt < 5)                                   return false;
  if (fabs(muEta) > 2.4 )                         return false;       
  if (fabs(muEta) < 2.4)
  {
    if(not (muon->typeBits & baconhep::kGlobal) ) return false;  
    if(not (muon->typeBits & baconhep::kPFMuon) ) return false;
    if(muon->muNchi2          >= 10)              return false;
    if(muon->nMatchStn <= 1)                      return false;
    if(muon->nPixHits  <= 0)                      return false; 
    if(muon->nTkLayers <= 5)                      return false;
    if(muon->nValidHits <= 0)                     return false;
    //if(fabs(muon->d0) > 0.2)                      return false;
    //if(fabs(muon->dz) > 0.5)                      return false;  
  }
  return true; 
}

Bool_t FatJetAnalyzer::isTightMuonMIT(baconhep::TMuon* muon, float muEta, float muPt)
{
  if (muPt < 5)                                   return false;
  if (fabs(muEta) > 2.4 )                         return false;
  if (fabs(muEta) < 2.4)
  {
    if(not (muon->typeBits & baconhep::kGlobal) ) return false;
    if(not (muon->typeBits & baconhep::kPFMuon) ) return false;
    if(muon->muNchi2          >= 10)              return false;
    if(muon->nMatchStn <= 1)                      return false;
    if(muon->nPixHits  <= 0)                      return false;
    if(muon->nTkLayers <= 5)                      return false;
    if(muon->nValidHits <= 0)                     return false;
    //if(fabs(muon->d0) > 0.02)                     return false;
    //if(fabs(muon->dz) > 0.1)                      return false;
  }
  return true;
}


Float_t FatJetAnalyzer::muonIso(baconhep::TMuon* muon)
{
  return (muon->chHadIso + std::max( 0., muon->neuHadIso + muon->gammaIso - 0.5*muon->puIso))/muon->pt; //what is muon->pfPt?
}

Float_t FatJetAnalyzer::electronIso(baconhep::TElectron* electron, float rhoFactor)//iso value set at 0.0588
{
  float EA = 0;  
  EAEle[0] = 0.1703;//0.13;
  EAEle[1] = 0.1715;//0.14;
  EAEle[2] = 0.1213;//0.07;
  EAEle[3] = 0.1230;//0.09;
  EAEle[4] = 0.1635;//0.11;
  EAEle[5] = 0.1937;//0.11;
  EAEle[6] = 0.2393;//0.14;
  if (fabs(electron->eta)  <  1.0) EA = EAEle[0];
  else if(fabs(electron->eta) < 1.479) EA = EAEle[1]; 
  else if(fabs(electron->eta) < 2.0) EA = EAEle[2];   
  else if(fabs(electron->eta) < 2.2) EA = EAEle[3];
  else if(fabs(electron->eta) < 2.3) EA = EAEle[4]; 
  else if(fabs(electron->eta) < 2.4) EA = EAEle[5];
  else if(fabs(electron->eta) > 2.4) EA = EAEle[6];
  
  return (electron->chHadIso + std::max(0., (double)electron->neuHadIso + electron->gammaIso - rhoFactor*EA))/electron->pt; 
}

Bool_t FatJetAnalyzer::isTightElectron(baconhep::TElectron* electron)
{
  float invE_invP = fabs(1./ (electron->eoverp * electron->pt) - electron->eoverp * 1./ (electron->eoverp * electron->pt));
  if(fabs(electron->scEta) > 2.5) return false;
  if(fabs(electron->scEta) < 1.479)//barrel
  {
    if(electron->sieie                       > 0.00998) return false;
    if(fabs(electron->dEtaInSeed)            > 0.00308) return false;
    if(fabs(electron->dPhiIn)                > 0.0816 ) return false;
    if(electron->hovere                      > 0.0414 ) return false;
    if(invE_invP                             > 0.0129 ) return false;
    if(electron->nMissingHits                > 1.0    ) return false;
    if(electron->isConv                               ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false; 
  }
  else if(fabs(electron->scEta) > 1.479)//endcap
  {
    if(electron->sieie                 > 0.0292 ) return false; 
    if(fabs(electron->dEtaInSeed)      > 0.00605) return false;
    if(fabs(electron->dPhiIn)          > 0.0394 ) return false;
    if(electron->hovere                > 0.0641 ) return false;
    if(invE_invP                       > 0.0129 ) return false;  
    if(electron->nMissingHits          > 1.0    ) return false;
    if(electron->isConv                         ) return false;
    //if(fabs(electron->d0)              > 0.05   ) return false;//0.10
    //if(fabs(electron->dz)              > 0.1    ) return false;//0.20 
  }
  return true;
}

Bool_t FatJetAnalyzer::isTightElectronMIT(baconhep::TElectron* electron)
{ 
  float invE_invP = fabs(1./ (electron->eoverp * electron->pt) - electron->eoverp * 1./ (electron->eoverp * electron->pt));
  if(fabs(electron->scEta) > 2.5) return false;
  if(fabs(electron->scEta) < 1.479)//barrel
  { 
    if(electron->sieie                       > 0.00998) return false;
    if(fabs(electron->dEtaInSeed)            > 0.00308) return false;
    if(fabs(electron->dPhiIn)                > 0.0816 ) return false;
    if(electron->hovere                      > 0.0414 ) return false;
    if(invE_invP                             > 0.0129 ) return false;
    if(electron->nMissingHits                >= 1     ) return false;
    if(electron->isConv                               ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false;
  }
  else if(fabs(electron->scEta) > 1.479)//endcap
  { 
    if(electron->sieie                       > 0.0292 ) return false;
    if(fabs(electron->dEtaInSeed)            > 0.00605) return false;
    if(fabs(electron->dPhiIn)                > 0.0394 ) return false;
    if(electron->hovere                      > 0.0641 ) return false;
    if(invE_invP                             > 0.0129 ) return false;
    if(electron->nMissingHits                >= 1     ) return false;
    if(electron->isConv                               ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;//0.10
    //if(fabs(electron->dz)                    > 0.1    ) return false;//0.20 
  }
  return true;
}

Bool_t FatJetAnalyzer::isHLTsafeElectronMIT(baconhep::TElectron* electron)
{
  float invE_invP = fabs(1./ (electron->eoverp * electron->pt) - electron->eoverp * 1./ (electron->eoverp * electron->pt));
  if(fabs(electron->scEta) > 2.5) return false;
  if(fabs(electron->scEta) < 1.479)//barrel
  {
    if(electron->sieie                       > 0.011  ) return false;
    if(fabs(electron->dEtaInSeed)            > 0.004  ) return false;
    if(fabs(electron->dPhiIn)                > 0.020  ) return false;
    if(electron->hovere                      > 0.060  ) return false;
    if(invE_invP                             > 0.013  ) return false;
    if(electron->nMissingHits                > 1      ) return false;
    if(electron->isConv                               ) return false;
    if(electron->ecalPFClusIso               > 0.160  ) return false;
    if(electron->hcalPFClusIso               > 0.120  ) return false; 
  }
  else if(fabs(electron->scEta) > 1.479)//endcap
  {
    if(electron->sieie                       > 0.0292 ) return false;
    if(fabs(electron->dEtaInSeed)            > 0.00605) return false;
    if(fabs(electron->dPhiIn)                > 0.0394 ) return false;
    if(electron->hovere                      > 0.0641 ) return false;
    if(invE_invP                             > 0.0129 ) return false;
    if(electron->nMissingHits                >= 1     ) return false;
    if(electron->isConv                               ) return false;
    if(electron->ecalPFClusIso               > 0.120  ) return false;
    if(electron->hcalPFClusIso               > 0.120  ) return false;
  }
  return true;
}


Bool_t FatJetAnalyzer::isTightElectronUCSD(baconhep::TElectron* electron)
{
  //double ip3d = sqrt((electron->d0)*(electron->d0) + (electron->dz)*(electron->dz));
  if(fabs(electron->scEta) > 2.5) return false;  
  if(fabs(electron->scEta) < 1.479)//barrel
  {
    if(electron->nMissingHits                >= 1     ) return false;
    //if(ip3d                                  > 0.015  ) return false;    
    //if(electron->sip3d                       > 4      ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false;
  }
  else if(fabs(electron->scEta) > 1.479)//endcap
  {
    if(electron->nMissingHits                >= 1     ) return false;
    //if(ip3d                                  > 0.015  ) return false;
    //if(electron->sip3d                       > 4      ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false;
  }
  return true;
}

//UCSD electron id loose defintion is tight + loosened isolation 

Bool_t FatJetAnalyzer::isVeryLooseElectronUCSD(baconhep::TElectron* electron)
{
  if(fabs(electron->scEta) > 2.5)                       return false;
  if(fabs(electron->scEta) < 1.479)//barrel
  {
    if(electron->nMissingHits                > 2      ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false;
  }
  else if(fabs(electron->scEta) > 1.479)//endcap
  {
    if(electron->nMissingHits                > 3      ) return false;
    //if(fabs(electron->d0)                    > 0.05   ) return false;
    //if(fabs(electron->dz)                    > 0.1    ) return false;
  }
  return true;
}  

Bool_t FatJetAnalyzer::isTightMuonUCSD(baconhep::TMuon* muon, float muEta, float muPt)
{
  double ip3d = sqrt(muon->d0*muon->d0 + muon->dz*muon->dz);
  if(muPt < 5)                                       return false;
  if(fabs(muEta) > 2.4  )                            return false;  
  if(fabs(muEta) < 2.4)
  {
    if(ip3d        > 0.015  )                        return false; 
    if(muon->sip3d > 4      )                        return false;
    if(fabs(muon->d0) > 0.05)                        return false;
    if(fabs(muon->dz) > 0.1 )                        return false; 
  }
  return true;
}

//UCSD muon id loose defintion is tight + loosened isolation 

Bool_t FatJetAnalyzer::isVeryLooseMuonUCSD(baconhep::TMuon* muon, float muEta, float muPt)
{
  if(muPt < 5)                                       return false;
  if (fabs(muEta) > 2.4  )                           return false; 
  if (fabs(muEta) < 2.4)
  { 
    if(fabs(muon->d0) > 0.05)                        return false;
    if(fabs(muon->dz) > 0.1 )                        return false;
  }
  return true;
}

// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<FatJetAnalyzer> selector(new FatJetAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
