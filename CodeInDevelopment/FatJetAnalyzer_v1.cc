#include "FatJetAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

FatJetAnalyzer::FatJetAnalyzer() : BLTSelector()
{

}

FatJetAnalyzer::~FatJetAnalyzer()
{

}

void FatJetAnalyzer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
    std::sregex_token_iterator(), std::back_inserter(options));

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);
    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "mumu" || params->selection == "emu")
    {
      triggerNames.push_back("HLT_IsoMu22_v*");
      triggerNames.push_back("HLT_IsoTkMu22_v*");
      //triggerNames.push_back("HLT_IsoMu24_eta2p1_v*");
      //triggerNames.push_back("HLT_IsoTkMu24_eta2p1_v*");
    }
    else if (params->selection == "ee")
    {
      triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }

    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false));
    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
      string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt";
      lumiMask.AddJSONFile(jsonFileName);
    }
    // muon momentum corrections
    muonCorr = new rochcor2016();
    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("eventWeightPU", &eventWeightPU);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    // leptons
    outTree->Branch("muon_pt", &muon_pt);
    outTree->Branch("muon_phi", &muon_phi); 
    outTree->Branch("muon_eta", &muon_eta);
    outTree->Branch("muon_iso", &muon_iso);
    outTree->Branch("muon_charge", &muon_charge);
    outTree->Branch("muon_trigger", &muon_trigger);
    outTree->Branch("muon_isoByPt", &muon_isoByPt);
    outTree->Branch("muon_id", &muon_id);
    outTree->Branch("muon_recoEW", &muon_recoEW);
    outTree->Branch("muon_triggerEW", &muon_triggerEW);
    outTree->Branch("electron_pt", &electron_pt);
    outTree->Branch("electron_phi", &electron_phi);
    outTree->Branch("electron_eta", &electron_eta);
    outTree->Branch("electron_charge", &electron_charge);
    outTree->Branch("electron_trigger", &electron_trigger);
    outTree->Branch("electron_id", &electron_id);
    // jets
    outTree->Branch("jet_pt", &jet_pt);
    outTree->Branch("jet_eta", &jet_eta);
    outTree->Branch("jet_phi", &jet_phi);
    outTree->Branch("jet_mass", &jet_mass);
    outTree->Branch("jet_D0", &jet_D0);
    outTree->Branch("jet_csv_data", &jet_csv_data);
    outTree->Branch("jet_csv_mc", &jet_csv_mc);
    outTree->Branch("jet_csv", &jet_csv);
    outTree->Branch("jet_flavor", &jet_flavor); 
    outTree->Branch("ak8jet_prunMass", &ak8jet_prunMass);
    outTree->Branch("ak8jet_trimMass", &ak8jet_trimMass);
    outTree->Branch("ak8jet_sd0", &ak8jet_sd0);
    outTree->Branch("ak8jet_pt", &ak8jet_pt);
    outTree->Branch("ak8jet_phi", &ak8jet_phi);
    outTree->Branch("ak8jet_eta", &ak8jet_eta);
    outTree->Branch("ak8jet_tau1", &ak8jet_tau1);
    outTree->Branch("ak8jet_tau2", &ak8jet_tau2);
    outTree->Branch("ak8jet_tau3", &ak8jet_tau3);
    outTree->Branch("gen_pdgId", &gen_pdgId);
    outTree->Branch("gen_status", &gen_status);
    outTree->Branch("gen_parent", &gen_parent);
    outTree->Branch("gen_pt", &gen_pt);
    outTree->Branch("gen_eta", &gen_eta);
    outTree->Branch("gen_phi", &gen_phi);
    outTree->Branch("gen_mass", &gen_mass);
    outTree->Branch("nPartons", &nPartons);
    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);

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
    electron_pt.clear();
    electron_eta.clear();
    electron_phi.clear();
    electron_charge.clear();
    electron_trigger.clear(); 
    electron_id.clear();
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    muon_iso.clear();
    muon_charge.clear();
    muon_trigger.clear();
    muon_isoByPt.clear();
    muon_id.clear();
    muon_recoEW.clear();
    muon_triggerEW.clear();
    jet_pt.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_mass.clear();
    jet_D0.clear();
    jet_csv_data.clear();
    jet_csv_mc.clear();
    jet_csv.clear();
    jet_flavor.clear();
    gen_pdgId.clear();
    gen_status.clear();
    gen_parent.clear();
    gen_pt.clear();
    gen_eta.clear();
    gen_phi.clear();
    gen_mass.clear();
    ak8jet_prunMass.clear();
    ak8jet_trimMass.clear();
    ak8jet_sd0.clear();
    ak8jet_pt.clear();
    ak8jet_phi.clear();
    ak8jet_eta.clear();
    ak8jet_tau1.clear();
    ak8jet_tau2.clear();
    ak8jet_tau3.clear();

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    // Apply lumi mask
    if (isData) 
    {
      RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
      //if(!lumiMask.HasRunLumi(rl)) 
      //return kTRUE;
    }
    hTotalEvents->Fill(2);
    /* Trigger selection */
    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) 
    {
      passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }
    if (!passTrigger && isData)
      return kTRUE;
    hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    eventWeightPU = 1.0;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPU           = fPVArr->GetEntries();
    if (isData) {
      eventWeight *= 1.;//weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
      eventWeightPU = 1.0;
    }
    else
    {
      eventWeightPU = weights->GetPUWeight(fInfo->nPUmean);
    }

    ///////////////////
    // Select objects//
    ///////////////////

    /* Vertices */
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

    /* MUONS */
    /* Apply a preselection so we can make a collection of muons to clean against */
    vector<TMuon*> muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) 
    {
      TMuon* muon = (TMuon*) fMuonArr->At(i);
      assert(muon);
      TLorentzVector muonP4;
      copy_p4(muon, MUON_MASS, muonP4);
      // Remove muon track pt from muon track isolation variable
      for (int k = i+1; k < fMuonArr->GetEntries(); k++) 
      {
        TLorentzVector muon_k;
        copy_p4(muon, MUON_MASS, muon_k);
        if (muonP4.DeltaR(muon_k) < 0.3) 
        {
          muon->trkIso = max(0., muon->trkIso - muon_k.Pt());
        }
      }
      // Apply rochester muon momentum corrections
      float qter = 1.;
      if (isData) 
      {
        muonCorr->momcor_data(muonP4, muon->q, 0, qter);
      } else 
      {
        muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
      }
      if (muonP4.Pt() > 25) 
      {
        muon_pt.push_back(muon->pt);
        muon_phi.push_back(muon->phi);
        muon_eta.push_back(muon->eta);
        muon_iso.push_back(muon->trkIso);
        muon_charge.push_back(muon->q);
        muon_isoByPt.push_back(muon->trkIso/muonP4.Pt());
        muon_id.push_back(isTightMuon(muon));
        TLorentzVector muon_p4;
        muon_p4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);
        muon_recoEW.push_back(weights->GetMuonRecoEff(muon_p4));    
        pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu22_v*", muon_p4);
        muon_triggerEW.push_back(trigEff.first);
        // trigger matching
        bool triggered = false;
        for (unsigned i = 0; i < triggerNames.size(); ++i) 
        {
          triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
        }
        muon_trigger.push_back(triggered);
      }
    }
    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
      TElectron* electron = (TElectron*) fElectronArr->At(i);
      assert(electron);
      TLorentzVector electronP4;
      copy_p4(electron, ELE_MASS, electronP4);
      electron_pt.push_back(electron->pt);
      electron_phi.push_back(electron->phi);
      electron_eta.push_back(electron->eta);
      electron_charge.push_back(electron->q);
      electron_id.push_back(isTightElectron(electron));
      // trigger matching
      bool triggered = false;
      for (unsigned i = 0; i < triggerNames.size(); ++i) 
      {
        triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
      }
      electron_trigger.push_back(triggered);
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);
    if (not isData)
    {
      /*GenParticles*/
      TClonesArray* genPartCollection;
      genPartCollection = fGenParticleArr; 
      //saving gen info
      unsigned int count = 0;
      for (int i=0; i < genPartCollection->GetEntries(); i++)
      {
        TGenParticle* genPart = (TGenParticle*) genPartCollection->At(i);   
        assert(genPart);
        gen_pdgId.push_back(genPart->pdgId);  
        gen_status.push_back(genPart->status);
        gen_parent.push_back(genPart->parent);
        gen_pt.push_back(genPart->pt);
        gen_eta.push_back(genPart->eta);
        gen_phi.push_back(genPart->phi);
        gen_mass.push_back(genPart->mass);
        if(genPart->status == 23 && (abs(genPart->pdgId) < 6 || genPart->pdgId == 21) && genPart->parent != -2) 
        {
          ++count;
        } 
      }
      nPartons = count; 
    }
    else 
    {
      nPartons = 0;
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
      if(jet->pt > 30 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
      {
        ak8jet_prunMass.push_back(ak8Jet->mass_prun);
        ak8jet_trimMass.push_back(ak8Jet->mass_trim);
        ak8jet_sd0.push_back(ak8Jet->mass_sd0);
        ak8jet_tau1.push_back(ak8Jet->tau1);
        ak8jet_tau2.push_back(ak8Jet->tau2);
        ak8jet_tau3.push_back(ak8Jet->tau3);
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
      if(jet->pt > 30 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
      {
        jet_pt.push_back(jet->pt);
        jet_eta.push_back(jet->eta);
        jet_phi.push_back(jet->phi);
        jet_mass.push_back(jet->mass);
        jet_D0.push_back(jet->d0);
        if(isData)
        {
          jet_csv_data.push_back(jet->csv);
        }
        else
        {
          jet_csv_mc.push_back(particleSelector->BTagModifier(jet, "CSVT"));
          jet_csv.push_back(jet->csv);
          jet_flavor.push_back(jet->hadronFlavor);
        }
      }  
    }

    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    /* MET */
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    /* Apply analysis selections */
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") 
    {
      hTotalEvents->Fill(5);
    } 
    else if (params->selection == "ee") 
    {
      hTotalEvents->Fill(5);
        
    } 
    else if (params->selection == "emu") 
    {
      hTotalEvents->Fill(5);
      /*if (!isData && false) 
      {
        eventWeight *= weights->GetMuonRecoEff(muons[0]);

        // trigger efficiency
        pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muons[0]);
        eventWeight *= trigEff.first/trigEff.second;
      }*/
    } 

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


Bool_t FatJetAnalyzer::isTightMuon(baconhep::TMuon* muon)
{
  if (muon->pt   < 5)          return false;
  if (fabs(muon->eta) < 2.4 and muon->typeBits & baconhep::kGlobal) 
  {
    if(muon->muNchi2   >= 10)  return false;
    if(muon->nMatchStn <= 1)   return false;
    if(muon->nPixHits  <= 0)   return false;
    if(fabs(muon->d0)  >= 0.2) return false;
    if(fabs(muon->dz)  >= 0.5) return false;
    if(muon->nTkLayers <= 5)   return false;
    if(muon->nValidHits <= 0)  return false;
  }
  return true;
}

Bool_t FatJetAnalyzer::isTightElectron(baconhep::TElectron* electron)
{
  if (electron->pt < 20)   return false;
  if (fabs(electron->eta) < 2.5)
  {
    if (not (particleSelector->PassElectronID(electron, cuts->tightElID)))               return false;
    if (not (particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl))) return false;
  }
  return true;
}
// _____________________________________________________________________________
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
