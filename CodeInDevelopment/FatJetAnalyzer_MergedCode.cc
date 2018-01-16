#include "FatJetAnalyzer.h"

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

    if (params->selection == "mumu" || params->selection == "emu")
    {
      triggerNames.push_back("HLT_IsoMu22_v*");
      triggerNames.push_back("HLT_IsoTkMu22_v*");
    }
    else if (params->selection == "ee")
    {
      triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }

    weights.reset(new WeightUtils(params->period, params->selection, false));
    lumiMask = RunLumiRangeMap();
    if (true) 
    { 
      string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt";
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
    outTree->Branch("nPU", &nPU);
    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);
    //leptons
    outTree->Branch("muon_pt", &muon_pt);
    outTree->Branch("muon_phi", &muon_phi); 
    outTree->Branch("muon_eta", &muon_eta);
    outTree->Branch("muon_iso", &muon_iso);
    outTree->Branch("muon_charge", &muon_charge);
    outTree->Branch("muon_id", &muon_id);
    outTree->Branch("muon_recoEW", &muon_recoEW);
    outTree->Branch("muon_triggerEW", &muon_triggerEW);
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
    muon_pt.clear();
    muon_eta.clear();
    muon_phi.clear();
    muon_iso.clear();
    muon_charge.clear();
    muon_id.clear();
    muon_recoEW.clear();
    muon_triggerEW.clear();
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
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPU           = fPVArr->GetEntries();
    if (isData) 
    {
      eventWeight *= 1.;//weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
      eventWeightPU = 1.0;
    }
    else
    {
      eventWeightPU = weights->GetPUWeight(fInfo->nPUmean);
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

    vector<TMuon*> muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) 
    {
      TMuon* muon = (TMuon*) fMuonArr->At(i);
      assert(muon);
      TLorentzVector muonP4;
      copy_p4(muon, MUON_MASS, muonP4);
      for (int k = i+1; k < fMuonArr->GetEntries(); k++) 
      {
        TLorentzVector muon_k;
        copy_p4(muon, MUON_MASS, muon_k);
        if (muonP4.DeltaR(muon_k) < 0.3) 
        {
          muon->trkIso = max(0., muon->trkIso - muon_k.Pt());
        }
      }
      float qter = 1.;
      if (isData) 
      {
        muonCorr->momcor_data(muonP4, muon->q, 0, qter);
      } else 
      {
        muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
      }
      if (muonP4.Pt() > 10) 
      {
        muon_pt.push_back(muon->pt);
        muon_phi.push_back(muon->phi);
        muon_eta.push_back(muon->eta);
        muon_iso.push_back(muon->trkIso);
        muon_charge.push_back(muon->q);
        muon_id.push_back(isTightMuon(muon));
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
        jet_csv.push_back(jet->csv);
      }  
    }
  
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

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
