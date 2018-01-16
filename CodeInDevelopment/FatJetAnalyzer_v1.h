// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef FATJETANALYZER_HH
#define FATJETANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"
#include "BLT/BLTAnalysis/interface/rochcor2016.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class FatJetAnalyzer: public BLTSelector {
public:
    FatJetAnalyzer();
    ~FatJetAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();
    Bool_t isTightMuon(baconhep::TMuon* muon);   
    Bool_t isTightElectron(baconhep::TElectron* electron);
    TFile *outFile;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    rochcor2016 *muonCorr;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
    UInt_t runNumber, lumiSection, nPU;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight;
    Float_t eventWeightPU;
    Float_t nPartons;
    Bool_t leptonOneTrigger, leptonTwoTrigger;

    Float_t met, metPhi;
    std::vector<float> electron_pt;
    std::vector<float> electron_phi;
    std::vector<float> electron_eta;
    std::vector<float> electron_charge;
    std::vector<float> electron_trigger;
    std::vector<bool>  electron_id;
    std::vector<float> muon_pt;
    std::vector<float> muon_phi;
    std::vector<float> muon_eta;
    std::vector<float> muon_iso;
    std::vector<float> muon_charge;
    std::vector<bool>  muon_trigger;
    std::vector<bool>  muon_id;
    std::vector<float>  muon_isoByPt;
    std::vector<float> muon_recoEW;
    std::vector<float> muon_triggerEW;
    std::vector<float> jet_pt;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_mass;
    std::vector<float> jet_D0;
    std::vector<float> jet_csv_data;
    std::vector<float> jet_csv_mc;
    std::vector<float> jet_csv;
    std::vector<float> jet_flavor;
    std::vector<float> ak8jet_prunMass;
    std::vector<float> ak8jet_trimMass;
    std::vector<float> ak8jet_sd0;
    std::vector<float> ak8jet_pt;
    std::vector<float> ak8jet_phi;
    std::vector<float> ak8jet_eta;
    std::vector<float> ak8jet_tau1;
    std::vector<float> ak8jet_tau2;
    std::vector<float> ak8jet_tau3;
    std::vector<float> gen_pdgId;
    std::vector<float> gen_status;
    std::vector<float> gen_parent;
    std::vector<float> gen_pt;
    std::vector<float> gen_eta;
    std::vector<float> gen_phi;
    std::vector<float> gen_mass;
    UInt_t nJets, nMuons, nElectrons;
    //ClassDef(FatJetAnalyzer,0);
};


#endif  // FATJETANALYZER_HH
