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
    Bool_t isTightMuon(baconhep::TMuon* muon, float muEta, float muPt);   
    Bool_t isTightMuonUCSD(baconhep::TMuon* muon, float muEta, float muPt);
    Bool_t isTightMuonMIT(baconhep::TMuon* muon, float muEta, float muPt);
    Bool_t isVeryLooseMuonUCSD(baconhep::TMuon* muon, float muEta, float muPt);
    Bool_t isTightElectron(baconhep::TElectron* electron);
    Bool_t isTightElectronUCSD(baconhep::TElectron* electron);
    Bool_t isVeryLooseElectronUCSD(baconhep::TElectron* electron);
    Bool_t isTightElectronMIT(baconhep::TElectron* electron);
    Bool_t isHLTsafeElectronMIT(baconhep::TElectron* electron);
    Float_t muonIso(baconhep::TMuon* muon);
    Float_t electronIso(baconhep::TElectron* electron, float rhoFactor);     
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
    UInt_t runNumber, lumiSection, nPU, nPV;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight;
    Float_t eventWeightPU;
    Float_t puWeight;
    Float_t eventWeightGen;
    Float_t topPtWeight;
    Float_t nPartons;
    Bool_t leptonOneTrigger, leptonTwoTrigger;

    Float_t met, metPhi;
    std::vector<float> muon_pt;
    std::vector<float> muon_phi;
    std::vector<float> muon_eta;
    std::vector<float> muon_pt_corr;
    std::vector<float> muon_phi_corr;
    std::vector<float> muon_eta_corr;
    std::vector<float> muon_trkIso;
    std::vector<float> muon_pfIso;
    std::vector<float> muon_charge;
    std::vector<bool>  muon_id;
    std::vector<bool>  muon_id_alternate;
    std::vector<bool>  muon_id_tightUCSD;
    std::vector<bool>  muon_id_tightMIT;
    std::vector<bool>  muon_id_veryLooseUCSD;
    std::vector<float> muon_recoEW;
    std::vector<float> muon_triggerEW;
    std::vector<float> muon_trigger;
    std::vector<float> muon_d0;
    std::vector<float> muon_dz;
    std::vector<float> muon_sip3d;
    std::vector<float> electron_pt;
    std::vector<float> electron_phi;
    std::vector<float> electron_eta;
    std::vector<float> electron_trkIso;
    std::vector<float> electron_pfIso;
    std::vector<float> electron_charge;
    std::vector<bool>  electron_id;
    std::vector<bool>  electron_id_tightUCSD;
    std::vector<bool>  electron_id_tightMIT;
    std::vector<bool>  electron_id_veryLooseUCSD;
    std::vector<bool>  electron_id_HLTsafeMIT;
    std::vector<float> electron_d0;
    std::vector<float> electron_dz;  
    std::vector<float> electron_sip3d;
    std::vector<float> electron_recoEW;
    std::vector<float> electron_triggerEW;
    std::vector<float> electron_trigger;
    std::vector<float> jet_pt;
    std::vector<float> jet_eta;
    std::vector<float> jet_phi;
    std::vector<float> jet_mass;
    std::vector<float> jet_csv;
    std::vector<float> ak8jet_prunMass;
    std::vector<float> ak8jet_trimMass;
    std::vector<float> ak8jet_sd0;
    std::vector<float> ak8jet_pt;
    std::vector<float> ak8jet_phi;
    std::vector<float> ak8jet_eta;
    std::vector<float> ak8jet_tau1;
    std::vector<float> ak8jet_tau2;
    float EAEle[7];
    UInt_t nJets, nMuons, nElectrons;
    //ClassDef(FatJetAnalyzer,0);
};


#endif  // FATJETANALYZER_HH

