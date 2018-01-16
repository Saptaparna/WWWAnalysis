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

    if (params->selection == "mumu" || params->selection == "emu") {
        triggerNames.push_back("HLT_IsoMu22_v*");
        triggerNames.push_back("HLT_IsoTkMu22_v*");
    } else if (params->selection == "ee") {
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
    outTree->Branch("nPU", &nPU);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);
    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);

    // jets
    outTree->Branch("jetOneP4", &jetOneP4);
    outTree->Branch("jetOneD0", &jetOneD0);
    outTree->Branch("jetOneTag", &jetOneTag);
    outTree->Branch("jetOneFlavor", &jetOneFlavor);

    outTree->Branch("jetTwoP4", &jetTwoP4);
    outTree->Branch("jetTwoD0", &jetTwoD0);
    outTree->Branch("jetTwoTag", &jetTwoTag);
    outTree->Branch("jetTwoFlavor", &jetTwoFlavor);

    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("bjetTag", &bjetTag);
    outTree->Branch("bjetFlavor", &bjetFlavor);

    outTree->Branch("genBJetP4", &genBJetP4);
    outTree->Branch("genBJetTag", &genBJetTag);
    outTree->Branch("genJetOneP4", &genJetOneP4);
    outTree->Branch("genJetOneTag", &genJetOneTag);
    outTree->Branch("genJetTwoP4", &genJetTwoP4);
    outTree->Branch("genJetTwoTag", &genJetTwoTag);
    outTree->Branch("genJetThreeP4", &genJetThreeP4);
    outTree->Branch("genJetThreeTag", &genJetThreeTag);

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
    outTree->Branch("ak8jet_tau1", &ak8jet_tau1);
    outTree->Branch("ak8jet_tau2", &ak8jet_tau2);
    outTree->Branch("ak8jet_tau3", &ak8jet_tau3);
    outTree->Branch("ca15jet_prunMass", &ca15jet_prunMass);
    outTree->Branch("ca15jet_trimMass", &ca15jet_trimMass);
    outTree->Branch("ca15jet_sd0", &ca15jet_sd0);
    outTree->Branch("ca15jet_pt", &ca15jet_pt);
    outTree->Branch("ca15jet_tau1", &ca15jet_tau1);
    outTree->Branch("ca15jet_tau2", &ca15jet_tau2);
    outTree->Branch("ca15jet_tau3", &ca15jet_tau3);
    outTree->Branch("gen_pdgId", &gen_pdgId);
    outTree->Branch("gen_status", &gen_status);
    outTree->Branch("gen_parent", &gen_parent);
    outTree->Branch("gen_pt", &gen_pt);
    outTree->Branch("gen_eta", &gen_eta);
    outTree->Branch("gen_phi", &gen_phi);
    outTree->Branch("gen_mass", &gen_mass);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

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
    ak8jet_tau1.clear();
    ak8jet_tau2.clear();
    ak8jet_tau3.clear();
    ca15jet_prunMass.clear();
    ca15jet_trimMass.clear();
    ca15jet_sd0.clear();
    ca15jet_pt.clear();
    ca15jet_tau1.clear();
    ca15jet_tau2.clear();
    ca15jet_tau3.clear();
    //if (fInfo->runNum != 275963 || fInfo->evtNum != 88834144)
    //    return kTRUE;

    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);

    // Apply lumi mask
    if (isData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }
    if (!passTrigger && isData)
        return kTRUE;

    hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPU           = fPVArr->GetEntries();
    if (!isData) {
        eventWeight *= 1.;//weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
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
    vector<TMuon*> tmp_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 5 
                && fabs(muon->eta) < 2.4
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
           ) {
            tmp_muons.push_back(muon);
        }
    }
    sort(tmp_muons.begin(), tmp_muons.end(), sort_by_higher_pt<TMuon>);

    // Second pass
    //int trigger_muon_index = -1;
    vector<TLorentzVector> muons;
    vector<TLorentzVector> veto_muons;
    vector<float> muons_iso;
    vector<float> muons_q;
    vector<bool> muons_trigger;
    for (unsigned i = 0; i < tmp_muons.size(); i++) {
        TMuon* muon = tmp_muons[i];
        TLorentzVector muonP4;
        copy_p4(tmp_muons[i], MUON_MASS, muonP4);

        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i+1; j < tmp_muons.size(); j++) {
            TLorentzVector muon_j;
            copy_p4(tmp_muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3) {
                muon->trkIso = max(0., muon->trkIso - muon_j.Pt());
                tmp_muons[j]->trkIso = max(0., tmp_muons[j]->trkIso - muonP4.Pt());
            }
        }

        // Apply rochester muon momentum corrections
        float qter = 1.;
        if (isData) {
            muonCorr->momcor_data(muonP4, muon->q, 0, qter);
        } else {
            muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
        }

        // Fill containers
        if (muon->trkIso/muonP4.Pt() < 0.1) {
            // For synchronization
            //cout << muonP4.Pt() << ", " << muon->pt 
            //     << ", " << muon->eta << ", " << muon->phi 
            //     << endl;

            if (muonP4.Pt() > 20) {
                veto_muons.push_back(muonP4);
            } 

            if (muonP4.Pt() > 25) {
                muons.push_back(muonP4);
                muons_iso.push_back(muon->trkIso);
                muons_q.push_back(muon->q);

                // trigger matching
                bool triggered = false;
                for (unsigned i = 0; i < triggerNames.size(); ++i) {
                    triggered |= trigger->passObj(triggerNames[i], 1, muon->hltMatchBits);
                }
                muons_trigger.push_back(triggered);
            }
        }
    }

    /* ELECTRONS */
    std::vector<TLorentzVector> electrons;
    vector<float> electrons_iso;
    vector<float> electrons_q;
    vector<bool> electrons_trigger;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 20 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            TLorentzVector electronP4;
            copy_p4(electron, ELE_MASS, electronP4);
            electrons.push_back(electronP4);
            electrons_iso.push_back(0.);
            electrons_q.push_back(electron->q);

            // trigger matching
            bool triggered = false;
            for (unsigned i = 0; i < triggerNames.size(); ++i) {
                triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            }
            electrons_trigger.push_back(triggered);
        }
    }

    std::sort(electrons.begin(), electrons.end(), P4SortCondition);

    /*GenParticles*/
    TClonesArray* genPartCollection;
    genPartCollection = fGenParticleArr; 
    //saving gen info
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
    }

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK4CHSArr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    std::vector<TJet*> genbjets;
    std::vector<TJet*> genjets;
    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;

    //saving all jet info
    for (int i=0; i < jetCollection->GetEntries(); i++)
    {   
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);
         
         // Prevent overlap of muons and jets
         TLorentzVector vJet;
         vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
         bool muOverlap = false;
         for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        if(jet->pt > 30 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
        { 
          if(fabs(jet->eta) <= 2.4  && !muOverlap && !elOverlap)//&& particleSelector->PassJetPUID(jet, cuts->looseJetID) && !muOverlap && !elOverlap)
          { 
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
            jet_pt.push_back(jet->pt);
            jet_eta.push_back(jet->eta);
            jet_phi.push_back(jet->phi);
            jet_mass.push_back(jet->mass);
            jet_D0.push_back(jet->d0);
          }
        }
    }
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
       }
    }
    
    TClonesArray* jetCA15Collection;
    jetCA15Collection = fAddCA15CHSArr;

    for (int i=0; i < jetCA15Collection->GetEntries(); i++)
    {   
       TAddJet* ca15Jet = (TAddJet*) jetCA15Collection->At(i);
       assert(ca15Jet);
       //TJet *jet    = (TJet*)fAK4CHSArr->At(ca15Jet->index);
       //if(jet->pt > 30 && fabs(jet->eta) < 4.7 && particleSelector->PassJetID(jet, cuts->looseJetID))
       //{
         ca15jet_prunMass.push_back(ca15Jet->mass_prun);
         ca15jet_trimMass.push_back(ca15Jet->mass_trim);
         ca15jet_sd0.push_back(ca15Jet->mass_sd0);
         //ca15jet_pt.push_back(jet->pt);
         ca15jet_tau1.push_back(ca15Jet->tau1);
         ca15jet_tau2.push_back(ca15Jet->tau2);
         ca15jet_tau3.push_back(ca15Jet->tau3);
       //}
    }     

    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

        //cout << jet->genpt << ", " << jet->genm << ", " 
        //     << jet->geneta << ", " << jet->genphi << ", " 
        //     << jet->partonFlavor << ", " << jet->hadronFlavor << ", " 
        //     << endl;
        //cout << jet->pt << ", " << jet->ptRaw << ", " 
        //     << jet->eta << ", " << jet->phi << ", " 
        //     << muOverlap << ", " << jet->csv << ", " 
        //     << particleSelector->PassJetID(jet, cuts->looseJetID) << endl;
        //cout << endl;

        if (!isData) {
            if (abs(jet->hadronFlavor) == 5) {
                genbjets.push_back(jet);
            } else {
                genjets.push_back(jet);
            }
        }

        if (
                jet->pt > 30 
                && fabs(jet->eta < 4.7)
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        jet->pt > 30 
                        //&& particleSelector->PassJetPUID(jet, cuts->looseJetID)
                        && !muOverlap 
                        && !elOverlap
                   ) { 
                    if (isData) {
                        if (jet->csv > 0.935) {
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "CSVT")) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            jets.push_back(jet);
                            ++nJets;
                        }
                    }
                }
            } else {
                if (jet->pt > 30) {
                    fwdjets.push_back(jet);
                    ++nFwdJets;
                }
            }
        }
    }

    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genjets.begin(), genjets.end(), sort_by_higher_pt<TJet>);
    std::sort(genbjets.begin(), genbjets.end(), sort_by_higher_pt<TJet>);

    // Add additional b jets to the central jet collection
    if (bjets.size() > 1) {
        jets.insert(jets.end(), bjets.begin()+1, bjets.end());
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

    if (params->selection == "mumu") {
        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        //if (muons[0].Pt() < 25 || fabs(muons[0].Eta()) > 2.1)
            //return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muons[0] + muons[1];
        //if (dimuon.M() < 8. || dimuon.M() > 70.)
            //return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4      = muons[0];
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 1;

        leptonTwoP4      = muons[1];
        leptonTwoIso     = muons_iso[1];
        leptonTwoQ       = muons_q[1];
        leptonTwoTrigger = muons_trigger[1];
        leptonTwoFlavor  = 1;

        if (!isData) {
            eventWeight *= weights->GetMuonRecoEff(muons[0]);
            eventWeight *= weights->GetMuonRecoEff(muons[1]);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu22_v*", muons[0]);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu22_v*", muons[1]);
            eventWeight *= 1 - (1 - trigEff1.first)*(1 - trigEff2.first);
        }

    } else if (params->selection == "ee") {

        if (electrons.size() != 2)
            //return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector dielectron;
        dielectron = electrons[0] + electrons[1];
        //if (dielectron.M() < 8. || dielectron.M() > 70.)
            //return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electrons[0];
        leptonOneIso     = electrons_iso[0];
        leptonOneQ       = electrons_q[0];
        leptonOneTrigger = electrons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[1];
        leptonTwoIso     = electrons_iso[1];
        leptonTwoQ       = electrons_q[1];
        leptonTwoTrigger = electrons_trigger[1];
        leptonTwoFlavor  = 0;
    } else if (params->selection == "emu") {

        //if (muons.size() != 1 || electrons.size() != 1)
            //return kTRUE;
        hTotalEvents->Fill(5);

        TLorentzVector dilepton;
        dilepton = muons[0] + electrons[0];
        //if (dilepton.M() < 12 || dilepton.M() > 70)
            //return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = muons[0];
        leptonOneIso     = muons_iso[0];
        leptonOneQ       = muons_q[0];
        leptonOneTrigger = muons_trigger[0];
        leptonOneFlavor  = 0;

        leptonTwoP4      = electrons[0];
        leptonTwoIso     = electrons_iso[0];
        leptonTwoQ       = electrons_q[0];
        leptonTwoTrigger = electrons_trigger[0];
        leptonTwoFlavor  = 1;

        if (!isData && false) {
            eventWeight *= weights->GetMuonRecoEff(muons[0]);

            // trigger efficiency
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muons[0]);
            eventWeight *= trigEff.first/trigEff.second;
        }
    } 


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0     = bjets[0]->d0;
        bjetTag    = bjets[0]->csv;
        bjetFlavor = bjets[0]->hadronFlavor;
    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0     = 0.;
        bjetTag    = 0.;
        bjetFlavor = 0.;
    }

    if (fwdjets.size() > 0) {
        jetOneP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetOneD0     = fwdjets[0]->d0;
        jetOneTag    = 0.;
        jetTwoTag    = 0.;
        jetOneFlavor = fwdjets[0]->hadronFlavor;
    } else if (jets.size() > 1) {
        jetOneP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetOneD0     = jets[0]->d0;
        jetOneTag    = jets[0]->csv;
        jetOneFlavor = jets[0]->hadronFlavor;
        jetTwoP4.SetPtEtaPhiM(jets[1]->pt, jets[1]->eta, jets[1]->phi, jets[1]->mass);
        jetTwoD0     = jets[1]->d0;
        jetTwoTag    = jets[1]->csv;
        jetTwoFlavor = jets[1]->hadronFlavor; 
    } else {
        jetOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetOneD0     = 0.;
        jetTwoD0     = 0.;
        jetOneTag    = 0.;
        jetTwoTag = 0.;
        jetOneFlavor = 0.;
        jetTwoFlavor = 0.;
    } 

    if (genbjets.size() > 0 && !isData) {
        genBJetP4.SetPtEtaPhiM(genbjets[0]->genpt, genbjets[0]->geneta, genbjets[0]->genphi, genbjets[0]->genm);
        genBJetTag = genbjets[0]->csv;
    } else {
        genBJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genBJetTag = 0;
    }

    if (genjets.size() > 2 && !isData) {
        genJetOneP4.SetPtEtaPhiM(genjets[0]->genpt, genjets[0]->geneta, genjets[0]->genphi, genjets[0]->genm);
        genJetOneTag = genjets[0]->csv;
        genJetTwoP4.SetPtEtaPhiM(genjets[1]->genpt, genjets[1]->geneta, genjets[1]->genphi, genjets[1]->genm);
        genJetTwoTag = genjets[1]->csv;
        genJetThreeP4.SetPtEtaPhiM(genjets[2]->genpt, genjets[2]->geneta, genjets[2]->genphi, genjets[2]->genm);
        genJetThreeTag = genjets[2]->csv;
    } else {
        genJetOneP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetOneTag = 0;
        genJetTwoP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetTwoTag = 0;
        genJetThreeP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genJetThreeTag = 0;
    }

    // Synchronization printout
    /*if (nBJets == 1 && nJets == 1 && nFwdJets == 1 && met < 40 && muons.size() >= 2) {
      cout << "Run: " << fInfo->runNum  
      << " Lumi: " << fInfo->lumiSec 
      << " Event: " << fInfo->evtNum << endl;

      cout << nBJets << ", " << nJets << ", " << nFwdJets << endl;

      cout << muons[0].Pt() << ", " << ", " << muons[0].Eta() << ", " 
      << muons[0].Phi() << endl; 
      cout << muons[1].Pt() << ", " << ", " << muons[1].Eta() << ", " 
      << muons[1].Phi() << endl; 

      cout << jetP4.Pt() << ", " << ", " << jetP4.Eta() << ", " 
      << jetP4.Phi() << ", " << jetTag << endl; 
      cout << bjetP4.Pt() << ", " << ", " << bjetP4.Eta() << ", " 
      << jetP4.Phi() << ", " << jetTag << endl; 

      cout << met << endl;

      TLorentzVector dijet = bjetP4 + jetP4;
      TLorentzVector dimuon = muons[0] + muons[1];
      cout << dimuon.Phi() << ", " << dijet.Phi() << ", " << fabs(dimuon.DeltaPhi(dijet)) << endl;
      }*/


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
