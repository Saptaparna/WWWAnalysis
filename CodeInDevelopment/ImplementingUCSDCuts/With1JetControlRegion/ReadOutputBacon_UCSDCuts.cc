#include "ReadOutputBacon_UCSDCuts.h"


int ReadOutputBacon_UCSDCuts(std::string infile, std::string treeStr, std::string outfile, std::string MCSample)
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
  UInt_t          nPV;
  Float_t         met;
  Float_t         metPhi;
  vector<float>   *muon_pt;
  vector<float>   *muon_phi;
  vector<float>   *muon_eta;
  vector<float>   *muon_trkIso;
  vector<float>   *muon_pfIso;
  vector<float>   *muon_charge;
  vector<float>   *muon_d0;
  vector<float>   *muon_dz;
  vector<float>   *muon_sip3d;
  vector<bool>    *muon_id;
  vector<bool>    *muon_id_alternate;
  vector<bool>    *muon_id_tightUCSD;
  vector<bool>    *muon_id_tightMIT;
  vector<bool>    *muon_id_veryLooseUCSD;
  vector<float>   *muon_recoEW;
  vector<float>   *muon_triggerEW;
  vector<float>   *muon_trigger;
  vector<float>   *electron_pt;
  vector<float>   *electron_phi;
  vector<float>   *electron_eta;
  vector<float>   *electron_pfIso;
  vector<float>   *electron_trkIso;
  vector<float>   *electron_charge;
  vector<float>   *electron_d0;
  vector<float>   *electron_dz;
  vector<float>   *electron_sip3d;
  vector<bool>    *electron_id;
  vector<bool>    *electron_id_tightUCSD;
  vector<bool>    *electron_id_tightMIT;
  vector<bool>    *electron_id_veryLooseUCSD;
  vector<bool>    *electron_id_HLTsafeMIT;
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
  muon_d0 = 0;
  muon_dz = 0;
  muon_sip3d = 0;
  muon_id = 0;
  muon_id_alternate = 0;
  muon_id_tightUCSD = 0;
  muon_id_tightMIT = 0;
  muon_id_veryLooseUCSD = 0;
  muon_recoEW = 0;
  muon_triggerEW = 0;
  muon_trigger = 0;
  electron_pt = 0;
  electron_phi = 0;
  electron_eta = 0;
  electron_pfIso = 0;
  electron_trkIso = 0;
  electron_charge = 0;
  electron_d0 = 0;
  electron_dz = 0;
  electron_sip3d = 0;
  electron_id = 0;
  electron_id_tightUCSD = 0;
  electron_id_tightMIT = 0;
  electron_id_veryLooseUCSD = 0;
  electron_id_HLTsafeMIT = 0;
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
  tree->SetBranchAddress("nPV", &(nPV));
  tree->SetBranchAddress("met", &(met)); 
  tree->SetBranchAddress("metPhi", &(metPhi));
  tree->SetBranchAddress("muon_pt", &(muon_pt));
  tree->SetBranchAddress("muon_phi", &(muon_phi));
  tree->SetBranchAddress("muon_eta", &(muon_eta));
  tree->SetBranchAddress("muon_trkIso", &(muon_trkIso));
  tree->SetBranchAddress("muon_pfIso", &(muon_pfIso));
  tree->SetBranchAddress("muon_charge", &(muon_charge));
  tree->SetBranchAddress("muon_d0", &(muon_d0));
  tree->SetBranchAddress("muon_dz", &(muon_dz));
  tree->SetBranchAddress("muon_sip3d", &(muon_sip3d));
  tree->SetBranchAddress("muon_id", &(muon_id));
  tree->SetBranchAddress("muon_id_alternate", &(muon_id_alternate));
  tree->SetBranchAddress("muon_id_tightUCSD", &(muon_id_tightUCSD));
  tree->SetBranchAddress("muon_id_tightMIT", &(muon_id_tightMIT));
  tree->SetBranchAddress("muon_id_veryLooseUCSD", &(muon_id_veryLooseUCSD));
  tree->SetBranchAddress("muon_recoEW", &(muon_recoEW));
  tree->SetBranchAddress("muon_triggerEW", &(muon_triggerEW));
  tree->SetBranchAddress("muon_trigger", &(muon_trigger));
  tree->SetBranchAddress("electron_pt", &(electron_pt));
  tree->SetBranchAddress("electron_phi", &(electron_phi));
  tree->SetBranchAddress("electron_eta", &(electron_eta));
  tree->SetBranchAddress("electron_pfIso", &(electron_pfIso));
  tree->SetBranchAddress("electron_trkIso", &(electron_trkIso));
  tree->SetBranchAddress("electron_charge", &(electron_charge));
  tree->SetBranchAddress("electron_d0", &(electron_d0));
  tree->SetBranchAddress("electron_dz", &(electron_dz));
  tree->SetBranchAddress("electron_sip3d", &(electron_sip3d));
  tree->SetBranchAddress("electron_id", &(electron_id));
  tree->SetBranchAddress("electron_id_tightUCSD", &(electron_id_tightUCSD));
  tree->SetBranchAddress("electron_id_tightMIT", &(electron_id_tightMIT));
  tree->SetBranchAddress("electron_id_veryLooseUCSD", &(electron_id_veryLooseUCSD));
  tree->SetBranchAddress("electron_id_HLTsafeMIT", &(electron_id_HLTsafeMIT));
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
  double cut1_elmu, cut2_elmu, cut3_elmu, cut4_elmu, cut5_elmu, cut6_elmu, cut7_elmu, cut8_elmu, cut9_elmu; 
  cut1_elmu = cut2_elmu = cut3_elmu = cut4_elmu = cut5_elmu = cut6_elmu = cut7_elmu = cut8_elmu = cut9_elmu = 0.0;
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
  HistCollection mumuHistCut12;
  initializeHistCollection(mumuHistCut12, "IpIsoTuning_MuMu");

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
  HistCollection elelHistCut12;
  initializeHistCollection(elelHistCut12, "IpIsoTuning_ElEl");

  HistCollection elmuHistCut1;
  initializeHistCollection(elmuHistCut1, "2SSTL_ElMu");
  HistCollection elmuHistCut2;
  initializeHistCollection(elmuHistCut2, "2SSTLLV_ElMu");
  HistCollection elmuHistCut3;
  initializeHistCollection(elmuHistCut3, "2SSTLLV2J_ElMu");
  HistCollection elmuHistCut4;
  initializeHistCollection(elmuHistCut4, "2SSTLLV2JBV_ElMu");
  HistCollection elmuHistCut5;
  initializeHistCollection(elmuHistCut5, "2SSTLLV2JBVCMjj_ElMu");
  HistCollection elmuHistCut6;
  initializeHistCollection(elmuHistCut6, "2SSTLLV2JBVCMjjLMjj_ElMu");
  HistCollection elmuHistCut7;
  initializeHistCollection(elmuHistCut7, "2SSTLLV2JBVCMjjLMjjMet_ElMu");
  HistCollection elmuHistCut8;
  initializeHistCollection(elmuHistCut8, "2SSTLLV2JBVCMjjLMjjMetMll_ElMu");
  HistCollection elmuHistCut12;
  initializeHistCollection(elmuHistCut12, "2SSTLLV2JBVCMjjLMjjMetMllMtmax_ElMu");
  HistCollection elmuHistCut9;
  initializeHistCollection(elmuHistCut9, "TTbarCR_ElMu");
  HistCollection elmuHistCut10;
  initializeHistCollection(elmuHistCut10, "TTbarCR1b_ElMu");
  HistCollection elmuHistCut11;
  initializeHistCollection(elmuHistCut11, "TTbarCR2b_ElMu");
  HistCollection elmuHistCut13;
  initializeHistCollection(elmuHistCut13, "IpIsoTuning_ElMu");


  TH1D *h_TotalEvents_MuMu = new TH1D("h_TotalEvents_MuMu", "h_TotalEvents_MuMu", 15, -0.5, 14.5); h_TotalEvents_MuMu->Sumw2();
  TH1D *h_TotalEvents_ElMu = new TH1D("h_TotalEvents_ElMu", "h_TotalEvents_ElMu", 15, -0.5, 14.5); h_TotalEvents_ElMu->Sumw2();
  TH1D *h_TotalEvents_ElEl = new TH1D("h_TotalEvents_ElEl", "h_TotalEvents_ElEl", 15, -0.5, 14.5); h_TotalEvents_ElEl->Sumw2();
  TH1D *h_nPU_MuMu = new TH1D("h_nPU_MuMu", "h_nPU_MuMu", 100, -0.5, 99.5); h_nPU_MuMu->Sumw2();
  TH1D *h_nPU_ElEl = new TH1D("h_nPU_ElEl", "h_nPU_ElEl", 100, -0.5, 99.5); h_nPU_ElEl->Sumw2();
  TH1D *h_nPU_ElMu = new TH1D("h_nPU_ElMu", "h_nPU_ElMu", 100, -0.5, 99.5); h_nPU_ElMu->Sumw2();
  
  TH1D *h_nPV_MuMu = new TH1D("h_nPV_MuMu", "h_nPV_MuMu", 100, -0.5, 99.5); h_nPV_MuMu->Sumw2();
  TH1D *h_nPV_ElEl = new TH1D("h_nPV_ElEl", "h_nPV_ElEl", 100, -0.5, 99.5); h_nPV_ElEl->Sumw2();
  TH1D *h_nPV_ElMu = new TH1D("h_nPV_ElMu", "h_nPV_ElMu", 100, -0.5, 99.5); h_nPV_ElMu->Sumw2();

  TH1D *h_leading_mu_IdpfIso_MuMu = new TH1D("h_leading_mu_IdpfIso_MuMu", "h_leading_mu_IdpfIso_MuMu", 10000.0, 0.0, 10.0);h_leading_mu_IdpfIso_MuMu->Sumw2();
  TH1D *h_trailing_mu_IdpfIso_MuMu = new TH1D("h_trailing_mu_IdpfIso_MuMu", "h_trailing_mu_IdpfIso_MuMu", 10000.0, 0.0, 10.0);h_trailing_mu_IdpfIso_MuMu->Sumw2();
  TH1D *h_leading_mu_IdtrkIso_MuMu = new TH1D("h_leading_mu_IdtrkIso_MuMu", "h_leading_mu_IdtrkIso_MuMu", 10000.0, 0.0, 10.0);h_leading_mu_IdtrkIso_MuMu->Sumw2();
  TH1D *h_trailing_mu_IdtrkIso_MuMu = new TH1D("h_trailing_mu_IdtrkIso_MuMu", "h_trailing_mu_IdtrkIso_MuMu", 10000.0, 0.0, 10.0);h_trailing_mu_IdtrkIso_MuMu->Sumw2();

  TH1D *h_leading_mu_IdpfIso_ElMu = new TH1D("h_leading_mu_IdpfIso_ElMu", "h_leading_mu_IdpfIso_ElMu", 10000.0, 0.0, 10.0);h_leading_mu_IdpfIso_ElMu->Sumw2();
  TH1D *h_leading_el_IdpfIso_ElMu = new TH1D("h_leading_el_IdpfIso_ElMu", "h_leading_el_IdpfIso_ElMu", 10000.0, 0.0, 10.0);h_leading_el_IdpfIso_ElMu->Sumw2();
  TH1D *h_leading_mu_IdtrkIso_ElMu = new TH1D("h_leading_mu_IdtrkIso_ElMu", "h_leading_mu_IdtrkIso_ElMu", 10000.0, 0.0, 10.0);h_leading_mu_IdtrkIso_ElMu->Sumw2();
  TH1D *h_leading_el_IdtrkIso_ElMu = new TH1D("h_leading_el_IdtrkIso_ElMu", "h_leading_el_IdtrkIso_ElMu", 10000.0, 0.0, 10.0);h_leading_el_IdtrkIso_ElMu->Sumw2();

  TH1D *h_leading_el_IdpfIso_ElEl = new TH1D("h_leading_el_IdpfIso_ElEl", "h_leading_el_IdpfIso_ElEl", 10000.0, 0.0, 10.0);h_leading_el_IdpfIso_ElEl->Sumw2();
  TH1D *h_trailing_el_IdpfIso_ElEl = new TH1D("h_trailing_el_IdpfIso_ElEl", "h_trailing_el_IdpfIso_ElEl", 10000.0, 0.0, 10.0);h_trailing_el_IdpfIso_ElEl->Sumw2();
  TH1D *h_leading_el_IdtrkIso_ElEl = new TH1D("h_leading_el_IdtrkIso_ElEl", "h_leading_el_IdtrkIso_ElEl", 10000.0, 0.0, 10.0);h_leading_el_IdtrkIso_ElEl->Sumw2();
  TH1D *h_trailing_el_IdtrkIso_ElEl = new TH1D("h_trailing_el_IdtrkIso_ElEl", "h_trailing_el_IdtrkIso_ElEl", 10000.0, 0.0, 10.0);h_trailing_el_IdtrkIso_ElEl->Sumw2();

  TH1D *h_leading_mu_pfIso_MuMu = new TH1D("h_leading_mu_pfIso_MuMu", "h_leading_mu_pfIso_MuMu", 10000.0, 0.0, 10.0);h_leading_mu_pfIso_MuMu->Sumw2();
  TH1D *h_trailing_mu_pfIso_MuMu = new TH1D("h_trailing_mu_pfIso_MuMu", "h_trailing_mu_pfIso_MuMu", 10000.0, 0.0, 10.0);h_trailing_mu_pfIso_MuMu->Sumw2();
  TH1D *h_leading_mu_trkIso_MuMu = new TH1D("h_leading_mu_trkIso_MuMu", "h_leading_mu_trkIso_MuMu", 10000.0, 0.0, 10.0);h_leading_mu_trkIso_MuMu->Sumw2();
  TH1D *h_trailing_mu_trkIso_MuMu = new TH1D("h_trailing_mu_trkIso_MuMu", "h_trailing_mu_trkIso_MuMu", 10000.0, 0.0, 10.0);h_trailing_mu_trkIso_MuMu->Sumw2();

  TH1D *h_leading_mu_pfIso_ElMu = new TH1D("h_leading_mu_pfIso_ElMu", "h_leading_mu_pfIso_ElMu", 10000.0, 0.0, 10.0);h_leading_mu_pfIso_ElMu->Sumw2();
  TH1D *h_leading_el_pfIso_ElMu = new TH1D("h_leading_el_pfIso_ElMu", "h_leading_el_pfIso_ElMu", 10000.0, 0.0, 10.0);h_leading_el_pfIso_ElMu->Sumw2();
  TH1D *h_leading_mu_trkIso_ElMu = new TH1D("h_leading_mu_trkIso_ElMu", "h_leading_mu_trkIso_ElMu", 10000.0, 0.0, 10.0);h_leading_mu_trkIso_ElMu->Sumw2();
  TH1D *h_leading_el_trkIso_ElMu = new TH1D("h_leading_el_trkIso_ElMu", "h_leading_el_trkIso_ElMu", 10000.0, 0.0, 10.0);h_leading_el_trkIso_ElMu->Sumw2();

  TH1D *h_leading_el_pfIso_ElEl = new TH1D("h_leading_el_pfIso_ElEl", "h_leading_el_pfIso_ElEl", 10000.0, 0.0, 10.0);h_leading_el_pfIso_ElEl->Sumw2();
  TH1D *h_trailing_el_pfIso_ElEl = new TH1D("h_trailing_el_pfIso_ElEl", "h_trailing_el_pfIso_ElEl", 10000.0, 0.0, 10.0);h_trailing_el_pfIso_ElEl->Sumw2();
  TH1D *h_leading_el_trkIso_ElEl = new TH1D("h_leading_el_trkIso_ElEl", "h_leading_el_trkIso_ElEl", 10000.0, 0.0, 10.0);h_leading_el_trkIso_ElEl->Sumw2();
  TH1D *h_trailing_el_trkIso_ElEl = new TH1D("h_trailing_el_trkIso_ElEl", "h_trailing_el_trkIso_ElEl", 10000.0, 0.0, 10.0);h_trailing_el_trkIso_ElEl->Sumw2();

  TH1D *h_leading_el_d0_IsoIp_ElEl = new TH1D("h_leading_el_d0_IsoIp_ElEl", "h_leading_el_d0_IsoIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_d0_IsoIp_ElEl->Sumw2();
  TH1D *h_trailing_el_d0_IsoIp_ElEl = new TH1D("h_trailing_el_d0_IsoIp_ElEl", "h_trailing_el_d0_IsoIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_d0_IsoIp_ElEl->Sumw2();

  TH1D *h_leading_el_dz_IsoIp_ElEl = new TH1D("h_leading_el_dz_IsoIp_ElEl", "h_leading_el_dz_IsoIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_dz_IsoIp_ElEl->Sumw2();
  TH1D *h_trailing_el_dz_IsoIp_ElEl = new TH1D("h_trailing_el_dz_IsoIp_ElEl", "h_trailing_el_dz_IsoIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_dz_IsoIp_ElEl->Sumw2();

  TH1D *h_leading_mu_d0_IdIsoIp_MuMu = new TH1D("h_leading_mu_d0_IdIsoIp_MuMu", "h_leading_mu_d0_IdIsoIp_MuMu", 20000.0, -10.0, 10.0);h_leading_mu_d0_IdIsoIp_MuMu->Sumw2();
  TH1D *h_trailing_mu_d0_IdIsoIp_MuMu = new TH1D("h_trailing_mu_d0_IdIsoIp_MuMu", "h_trailing_mu_d0_IdIsoIp_MuMu", 20000.0, -10.0, 10.0);h_trailing_mu_d0_IdIsoIp_MuMu->Sumw2();

  TH1D *h_leading_mu_dz_IdIsoIp_MuMu = new TH1D("h_leading_mu_dz_IdIsoIp_MuMu", "h_leading_mu_dz_IdIsoIp_MuMu", 20000.0, -10.0, 10.0);h_leading_mu_dz_IdIsoIp_MuMu->Sumw2();
  TH1D *h_trailing_mu_dz_IdIsoIp_MuMu = new TH1D("h_trailing_mu_dz_IdIsoIp_MuMu", "h_trailing_mu_dz_IdIsoIp_MuMu", 20000.0, -10.0, 10.0);h_trailing_mu_dz_IdIsoIp_MuMu->Sumw2();

  TH1D *h_leading_mu_d0_IdIsoIp_ElMu = new TH1D("h_leading_mu_d0_IdIsoIp_ElMu", "h_leading_mu_d0_IdIsoIp_ElMu", 20000.0, -10.0, 10.0);h_leading_mu_d0_IdIsoIp_ElMu->Sumw2();
  TH1D *h_leading_el_d0_IdIsoIp_ElMu = new TH1D("h_leading_el_d0_IdIsoIp_ElMu", "h_leading_el_d0_IdIsoIp_ElMu", 20000.0, -10.0, 10.0);h_leading_el_d0_IdIsoIp_ElMu->Sumw2();

  TH1D *h_leading_mu_dz_IdIsoIp_ElMu = new TH1D("h_leading_mu_dz_IdIsoIp_ElMu", "h_leading_mu_dz_IdIsoIp_ElMu", 20000.0, -10.0, 10.0);h_leading_mu_dz_IdIsoIp_ElMu->Sumw2();
  TH1D *h_leading_el_dz_IdIsoIp_ElMu = new TH1D("h_leading_el_dz_IdIsoIp_ElMu", "h_leading_el_dz_IdIsoIp_ElMu", 20000.0, -10.0, 10.0);h_leading_el_dz_IdIsoIp_ElMu->Sumw2();

  TH1D *h_leading_el_d0_IdIsoIp_ElEl = new TH1D("h_leading_el_d0_IdIsoIp_ElEl", "h_leading_el_d0_IdIsoIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_d0_IdIsoIp_ElEl->Sumw2();
  TH1D *h_trailing_el_d0_IdIsoIp_ElEl = new TH1D("h_trailing_el_d0_IdIsoIp_ElEl", "h_trailing_el_d0_IdIsoIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_d0_IdIsoIp_ElEl->Sumw2();

  TH1D *h_leading_el_dz_IdIsoIp_ElEl = new TH1D("h_leading_el_dz_IdIsoIp_ElEl", "h_leading_el_dz_IdIsoIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_dz_IdIsoIp_ElEl->Sumw2();
  TH1D *h_trailing_el_dz_IdIsoIp_ElEl = new TH1D("h_trailing_el_dz_IdIsoIp_ElEl", "h_trailing_el_dz_IdIsoIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_dz_IdIsoIp_ElEl->Sumw2();

  TH1D *h_leading_mu_d0_IdIp_MuMu = new TH1D("h_leading_mu_d0_IdIp_MuMu", "h_leading_mu_d0_IdIp_MuMu", 20000.0, -10.0, 10.0);h_leading_mu_d0_IdIp_MuMu->Sumw2();
  TH1D *h_trailing_mu_d0_IdIp_MuMu = new TH1D("h_trailing_mu_d0_IdIp_MuMu", "h_trailing_mu_d0_IdIp_MuMu", 20000.0, -10.0, 10.0);h_trailing_mu_d0_IdIp_MuMu->Sumw2();

  TH1D *h_leading_mu_dz_IdIp_MuMu = new TH1D("h_leading_mu_dz_IdIp_MuMu", "h_leading_mu_dz_IdIp_MuMu", 20000.0, -10.0, 10.0);h_leading_mu_dz_IdIp_MuMu->Sumw2();
  TH1D *h_trailing_mu_dz_IdIp_MuMu = new TH1D("h_trailing_mu_dz_IdIp_MuMu", "h_trailing_mu_dz_IdIp_MuMu", 20000.0, -10.0, 10.0);h_trailing_mu_dz_IdIp_MuMu->Sumw2();

  TH1D *h_leading_mu_d0_IdIp_ElMu = new TH1D("h_leading_mu_d0_IdIp_ElMu", "h_leading_mu_d0_IdIp_ElMu", 20000.0, -10.0, 10.0);h_leading_mu_d0_IdIp_ElMu->Sumw2();
  TH1D *h_leading_el_d0_IdIp_ElMu = new TH1D("h_leading_el_d0_IdIp_ElMu", "h_leading_el_d0_IdIp_ElMu", 20000.0, -10.0, 10.0);h_leading_el_d0_IdIp_ElMu->Sumw2();

  TH1D *h_leading_mu_dz_IdIp_ElMu = new TH1D("h_leading_mu_dz_IdIp_ElMu", "h_leading_mu_dz_IdIp_ElMu", 20000.0, -10.0, 10.0);h_leading_mu_dz_IdIp_ElMu->Sumw2();
  TH1D *h_leading_el_dz_IdIp_ElMu = new TH1D("h_leading_el_dz_IdIp_ElMu", "h_leading_el_dz_IdIp_ElMu", 20000.0, -10.0, 10.0);h_leading_el_dz_IdIp_ElMu->Sumw2();

  TH1D *h_leading_el_d0_IdIp_ElEl = new TH1D("h_leading_el_d0_IdIp_ElEl", "h_leading_el_d0_IdIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_d0_IdIp_ElEl->Sumw2();
  TH1D *h_trailing_el_d0_IdIp_ElEl = new TH1D("h_trailing_el_d0_IdIp_ElEl", "h_trailing_el_d0_IdIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_d0_IdIp_ElEl->Sumw2();

  TH1D *h_leading_el_dz_IdIp_ElEl = new TH1D("h_leading_el_dz_IdIp_ElEl", "h_leading_el_dz_IdIp_ElEl", 20000.0, -10.0, 10.0);h_leading_el_dz_IdIp_ElEl->Sumw2();
  TH1D *h_trailing_el_dz_IdIp_ElEl = new TH1D("h_trailing_el_dz_IdIp_ElEl", "h_trailing_el_dz_IdIp_ElEl", 20000.0, -10.0, 10.0);h_trailing_el_dz_IdIp_ElEl->Sumw2();
  
  TH1D *h_leading_mu_Iso_Tune_MuMu = new TH1D("h_leading_mu_Iso_Tune_MuMu", "h_leading_mu_Iso_Tune_MuMu", 10000.0, 0.0, 10.0); h_leading_mu_Iso_Tune_MuMu->Sumw2();
  TH1D *h_trailing_mu_Iso_Tune_MuMu = new TH1D("h_trailing_mu_Iso_Tune_MuMu", "h_trailing_mu_Iso_Tune_MuMu", 10000.0, 0.0, 10.0); h_trailing_mu_Iso_Tune_MuMu->Sumw2();
  
  TH1D *h_leading_el_Iso_Tune_ElEl = new TH1D("h_leading_el_Iso_Tune_ElEl", "h_leading_el_Iso_Tune_ElEl", 10000.0, 0.0, 10.0); h_leading_el_Iso_Tune_ElEl->Sumw2();
  TH1D *h_trailing_el_Iso_Tune_ElEl = new TH1D("h_trailing_el_Iso_Tune_ElEl", "h_trailing_el_Iso_Tune_ElEl", 10000.0, 0.0, 10.0); h_trailing_el_Iso_Tune_ElEl->Sumw2();

  TH1D *h_leading_mu_Iso_Tune_ElMu = new TH1D("h_leading_mu_Iso_Tune_ElMu", "h_leading_mu_Iso_Tune_ElMu", 10000.0, 0.0, 10.0); h_leading_mu_Iso_Tune_ElMu->Sumw2();
  TH1D *h_leading_el_Iso_Tune_ElMu = new TH1D("h_leading_el_Iso_Tune_ElMu", "h_leading_el_Iso_Tune_ElMu", 10000.0, 0.0, 10.0); h_leading_el_Iso_Tune_ElMu->Sumw2();  
 
  TH1D *h_leading_mu_d0_Tune_MuMu = new TH1D("h_leading_mu_d0_Tune_MuMu", "h_leading_mu_d0_Tune_MuMu", 20000.0, -10.0, 10.0); h_leading_mu_d0_Tune_MuMu->Sumw2();
  TH1D *h_trailing_mu_d0_Tune_MuMu = new TH1D("h_trailing_mu_d0_Tune_MuMu", "h_trailing_mu_d0_Tune_MuMu", 20000.0, -10.0, 10.0); h_trailing_mu_d0_Tune_MuMu->Sumw2();

  TH1D *h_leading_el_d0_Tune_ElEl = new TH1D("h_leading_el_d0_Tune_ElEl", "h_leading_el_d0_Tune_ElEl", 20000.0, -10.0, 10.0); h_leading_el_d0_Tune_ElEl->Sumw2();
  TH1D *h_trailing_el_d0_Tune_ElEl = new TH1D("h_trailing_el_d0_Tune_ElEl", "h_trailing_el_d0_Tune_ElEl", 20000.0, -10.0, 10.0); h_trailing_el_d0_Tune_ElEl->Sumw2();

  TH1D *h_leading_mu_d0_Tune_ElMu = new TH1D("h_leading_mu_d0_Tune_ElMu", "h_leading_mu_d0_Tune_ElMu", 20000.0, -10.0, 10.0); h_leading_mu_d0_Tune_ElMu->Sumw2();
  TH1D *h_leading_el_d0_Tune_ElMu = new TH1D("h_leading_el_d0_Tune_ElMu", "h_leading_el_d0_Tune_ElMu", 20000.0, -10.0, 10.0); h_leading_el_d0_Tune_ElMu->Sumw2();

  TH1D *h_leading_mu_dz_Tune_MuMu = new TH1D("h_leading_mu_dz_Tune_MuMu", "h_leading_mu_dz_Tune_MuMu", 20000.0, -10.0, 10.0); h_leading_mu_dz_Tune_MuMu->Sumw2();
  TH1D *h_trailing_mu_dz_Tune_MuMu = new TH1D("h_trailing_mu_dz_Tune_MuMu", "h_trailing_mu_dz_Tune_MuMu", 20000.0, -10.0, 10.0); h_trailing_mu_dz_Tune_MuMu->Sumw2();

  TH1D *h_leading_el_dz_Tune_ElEl = new TH1D("h_leading_el_dz_Tune_ElEl", "h_leading_el_dz_Tune_ElEl", 20000.0, -10.0, 10.0); h_leading_el_dz_Tune_ElEl->Sumw2();
  TH1D *h_trailing_el_dz_Tune_ElEl = new TH1D("h_trailing_el_dz_Tune_ElEl", "h_trailing_el_dz_Tune_ElEl", 20000.0, -10.0, 10.0); h_trailing_el_dz_Tune_ElEl->Sumw2();

  TH1D *h_leading_mu_dz_Tune_ElMu = new TH1D("h_leading_mu_dz_Tune_ElMu", "h_leading_mu_dz_Tune_ElMu", 20000.0, -10.0, 10.0); h_leading_mu_dz_Tune_ElMu->Sumw2();
  TH1D *h_leading_el_dz_Tune_ElMu = new TH1D("h_leading_el_dz_Tune_ElMu", "h_leading_el_dz_Tune_ElMu", 20000.0, -10.0, 10.0); h_leading_el_dz_Tune_ElMu->Sumw2();
 
  //nEvents = 1000000;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);

    std::vector<leptonInfo> v_muons;
    std::vector<leptonInfo> v_looseMuons;
    std::vector<leptonInfo> v_muplotIso;
    std::vector<leptonInfo> v_muplotIdIso;
    std::vector<leptonInfo> v_muplotIdIp;
    std::vector<leptonInfo> v_muplotIdIsoIp;
    std::vector<leptonInfo> v_muIpIsoTune;
    for (unsigned int imuon=0; imuon<muon_pt->size(); imuon++)
    {
      leptonInfo muon;
      muon.pt = muon_pt->at(imuon); 
      muon.phi = muon_phi->at(imuon);
      muon.eta = muon_eta->at(imuon);
      muon.pfiso = muon_pfIso->at(imuon);
      muon.trkiso = muon_trkIso->at(imuon);
      muon.charge = (int)muon_charge->at(imuon);
      muon.d0  = muon_d0->at(imuon);
      muon.dz  = muon_dz->at(imuon);
      muon.sip3d  = muon_sip3d->at(imuon);
      muon.id = muon_id->at(imuon);
      muon.recoEW = muon_recoEW->at(imuon);
      muon.trigger = muon_trigger->at(imuon);
      muon.triggerEW = muon_triggerEW->at(imuon);
      muon.isTightUCSD = muon_id_tightUCSD->at(imuon);
      muon.isTightMIT = muon_id_tightMIT->at(imuon);
      muon.isveryLooseUCSD = muon_id_veryLooseUCSD->at(imuon);
      if(muon.id==1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4) v_muplotIdIso.push_back(muon);
      if(muon.pt > 30.0 and fabs(muon.eta) < 2.4) v_muplotIso.push_back(muon); 
      if(muon.isTightMIT==1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4) v_muplotIdIp.push_back(muon);
      if(muon.isTightMIT==1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4 and muon.pfiso < 0.06) v_muplotIdIsoIp.push_back(muon);
      if(muon.isTightMIT==1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4) v_muIpIsoTune.push_back(muon);
      //if(muon.isTightMIT==1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4 and muon.pfiso > 0.06) std::cout << "muon.pfiso = " << muon.pfiso << std::endl;
      //if(muon.id==1 and muon.trkiso/muon.pt < 0.1 and muon.pt > 30.0 and fabs(muon.eta) < 2.4) 
      //if(muon.id==1 and muon.pfiso < 0.15 and muon.pt > 30.0 and fabs(muon.eta) < 2.4)
      if(muon.id==1 and muon.pfiso < 0.06 and muon.pt > 30.0 and fabs(muon.eta) < 2.4)
      //if(muon.isTightUCSD==1 and muon.pfiso < 0.06 and muon.pt > 30.0 and fabs(muon.eta) < 2.4)//UCSD muon cuts
      {
        v_muons.push_back(muon);
      }
      //if(muon.pt > 10.0 and muon.id==1 and muon.pfiso < 0.40) v_looseMuons.push_back(muon);
      if(muon.pt > 10.0 and muon.id==0) v_looseMuons.push_back(muon);
    }//muon loop

    //if(v_muIpIsoTune.size() > 0 and v_muIpIsoTune.at(0).pfiso > 0.06) std::cout << "muon.at(0).pfiso = " << v_muIpIsoTune.at(0).pfiso << std::endl;

    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);
    std::sort (v_looseMuons.begin(), v_looseMuons.end(), sortLeptonsInDescendingpT);
    std::sort (v_muplotIdIso.begin(), v_muplotIdIso.end(), sortLeptonsInDescendingpT);
    std::sort (v_muplotIso.begin(), v_muplotIso.end(), sortLeptonsInDescendingpT);     
    std::sort (v_muplotIdIp.begin(), v_muplotIdIp.end(), sortLeptonsInDescendingpT);
    std::sort (v_muplotIdIsoIp.begin(), v_muplotIdIsoIp.end(), sortLeptonsInDescendingpT);
    std::sort (v_muIpIsoTune.begin(), v_muIpIsoTune.end(), sortLeptonsInDescendingpT);

    std::vector<leptonInfo> v_electrons;
    std::vector<leptonInfo> v_looseElectrons;
    std::vector<leptonInfo> v_elplotIdIso;
    std::vector<leptonInfo> v_elplotIso;
    std::vector<leptonInfo> v_elplotIdIp;
    std::vector<leptonInfo> v_elplotIdIsoIp;
    std::vector<leptonInfo> v_elIpIsoTune;
    for (unsigned int ielectron=0; ielectron<electron_pt->size(); ielectron++)
    {
      leptonInfo electron;
      electron.pt = electron_pt->at(ielectron);
      electron.phi = electron_phi->at(ielectron);
      electron.eta = electron_eta->at(ielectron);
      electron.pfiso = electron_pfIso->at(ielectron);
      electron.trkiso = electron_trkIso->at(ielectron);
      electron.charge = (int)electron_charge->at(ielectron);
      electron.d0  = electron_d0->at(ielectron);
      electron.dz  = electron_dz->at(ielectron);
      electron.sip3d  = electron_sip3d->at(ielectron);
      electron.id = electron_id->at(ielectron);
      electron.isTightUCSD = electron_id_tightUCSD->at(ielectron);
      electron.isTightMIT = electron_id_tightMIT->at(ielectron);
      electron.isveryLooseUCSD = electron_id_veryLooseUCSD->at(ielectron);
      electron.isHLTsafeMIT = electron_id_HLTsafeMIT->at(ielectron);
      electron.recoEW = electron_recoEW->at(ielectron);
      electron.trigger = electron_trigger->at(ielectron);
      electron.triggerEW = electron_triggerEW->at(ielectron);
      if(electron.isTightMIT==1 and electron.pt > 30.0 and fabs(electron.eta) < 2.4) v_elIpIsoTune.push_back(electron);
      if(electron.isTightMIT==1 and electron.pt > 30.0 and fabs(electron.eta) < 2.4) v_elplotIdIp.push_back(electron);
      if(fabs(electron.eta) < 1.479)
      {
        if(electron.isTightMIT==1 and electron.pt > 30.0 and electron.pfiso < 0.0588) v_elplotIdIsoIp.push_back(electron);
      }
      else if(fabs(electron.eta) > 1.479 and fabs(electron.eta) < 2.4)
      {
        if(electron.isTightMIT==1 and electron.pt > 30.0 and electron.pfiso < 0.0571) v_elplotIdIsoIp.push_back(electron);
      }
      if(electron.id==1 and electron.pt > 30.0 and fabs(electron.eta) < 2.4) v_elplotIdIso.push_back(electron); 
      if(electron.pt > 30.0 and fabs(electron.eta) < 2.4) v_elplotIso.push_back(electron);
      if(fabs(electron.eta) < 1.479)
      {
        if(electron.id==1 and electron.pfiso < 0.0588 and electron.pt > 30.0 and fabs(electron.d0) < 0.05 and fabs(electron.dz) < 0.1) v_electrons.push_back(electron);
      }
      else if(fabs(electron.eta) > 1.479 and fabs(electron.eta) < 2.4)
      {
        if(electron.id==1 and electron.pfiso < 0.0571 and electron.pt > 30.0 and fabs(electron.d0) < 0.05 and fabs(electron.dz) < 0.1) v_electrons.push_back(electron);
      }
      if(electron.pt > 10.0 and electron.id==0) v_looseElectrons.push_back(electron);
    }//electron loop

    std::sort (v_electrons.begin(), v_electrons.end(), sortLeptonsInDescendingpT);
    std::sort (v_looseElectrons.begin(), v_looseElectrons.end(), sortLeptonsInDescendingpT);
    std::sort (v_elplotIdIso.begin(), v_elplotIdIso.end(), sortLeptonsInDescendingpT);
    std::sort (v_elplotIso.begin(), v_elplotIso.end(), sortLeptonsInDescendingpT); 
    std::sort (v_elplotIdIp.begin(), v_elplotIdIp.end(), sortLeptonsInDescendingpT);
    std::sort (v_elplotIdIsoIp.begin(), v_elplotIdIsoIp.end(), sortLeptonsInDescendingpT);

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
  
      if(fabs(v_jets.at(iselJet).jetEta)<5.0 and v_jets.at(iselJet).jetPt>30.0)
      {
        bool isGoodJet=true;
        Jet.JetLV.SetPtEtaPhiM(v_jets.at(iselJet).jetPt, v_jets.at(iselJet).jetEta, v_jets.at(iselJet).jetPhi, v_jets.at(iselJet).jetMass);
        Jet.BTag_CSV = v_jets.at(iselJet).jetCSV;
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
      //if(isGoodJet and Jet.BTag_CSV > 0.5426) nbJets++;
      }//jet 4 vector closed
    }//selected jet loop closed

    std::sort (v_selectedJets.begin(), v_selectedJets.end(), sortJetVectorsInDescendingpT);
    int nbJets_loose  = 0;
    int nbJets_tight  = 0;
    std::vector<AnalysisJetInfo> v_selectedBJets;
    for(unsigned int iselbJet=0; iselbJet<v_jets.size(); ++iselbJet)
    {
      AnalysisJetInfo bJet;
      //if(v_jets.at(iselbJet).jetCSV==-10) continue;
      if(fabs(v_jets.at(iselbJet).jetEta)<2.4 and v_jets.at(iselbJet).jetPt > 20.0) bJet.BTag_CSV = v_jets.at(iselbJet).jetCSV;
      if(fabs(v_jets.at(iselbJet).jetEta)<2.4 and v_jets.at(iselbJet).jetPt > 20.0 and v_jets.at(iselbJet).jetCSV > 0.5426) 
      {
        bool isGoodbJet=true;
        bJet.JetLV.SetPtEtaPhiM(v_jets.at(iselbJet).jetPt, v_jets.at(iselbJet).jetEta, v_jets.at(iselbJet).jetPhi, v_jets.at(iselbJet).jetMass); 
        for(unsigned int iselMuon=0; iselMuon<v_muons.size(); ++iselMuon)
        {
          TLorentzVector Muon = fillTLorentzVector(v_muons.at(iselMuon).pt, v_muons.at(iselMuon).eta, v_muons.at(iselMuon).phi, MUON_MASS);
          double DRjet_mu = bJet.JetLV.DeltaR(Muon);
          if(DRjet_mu<0.5) isGoodbJet=false;
        }//muon loop closed
        for(unsigned int iselElectron=0; iselElectron<v_electrons.size(); ++iselElectron)
        {
          TLorentzVector Electron = fillTLorentzVector(v_electrons.at(iselElectron).pt, v_electrons.at(iselElectron).eta, v_electrons.at(iselElectron).phi, ELECTRON_MASS);
          double DRjet_el = bJet.JetLV.DeltaR(Electron);
          if(DRjet_el<0.5) isGoodbJet=false;
        }//electron loop closed
        if(isGoodbJet) v_selectedBJets.push_back(bJet);
        if(isGoodbJet) nbJets_loose++;
      }//bjet 4 vector
      if(v_jets.at(iselbJet).jetCSV > 0.8484) nbJets_tight++;
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
    //control region plots 2 muons without iso or ip cuts and 0 jet bin
    //if(v_muIpIsoTune.size()==2 and v_muIpIsoTune.at(0).charge*v_muIpIsoTune.at(1).charge==1 and v_elIpIsoTune.size()==0 and v_muIpIsoTune.at(0).pfiso > 0.06) std::cout << "muon.at(0).pfiso = " << v_muIpIsoTune.at(0).pfiso << std::endl;
    if(v_muIpIsoTune.size()==2 and v_electrons.size()==0 and v_selectedJets.size()==1)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muIpIsoTune.at(0).recoEW*v_muIpIsoTune.at(1).recoEW*(1 - (1 - v_muIpIsoTune.at(0).triggerEW)*(1 - v_muIpIsoTune.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muIpIsoTune.at(0).recoEW*v_muIpIsoTune.at(1).recoEW*(1 - (1 - v_muIpIsoTune.at(0).triggerEW)*(1 - v_muIpIsoTune.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muIpIsoTune.at(0).charge*v_muIpIsoTune.at(1).charge==1 and v_muIpIsoTune.at(0).pt > 30.0 and v_muIpIsoTune.at(1).pt > 30.0)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muIpIsoTune.at(0).pt, v_muIpIsoTune.at(0).eta, v_muIpIsoTune.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muIpIsoTune.at(1).pt, v_muIpIsoTune.at(1).eta, v_muIpIsoTune.at(1).phi, MUON_MASS);
        fillMuHistCollection(mumuHistCut12, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
        h_leading_mu_Iso_Tune_MuMu->Fill(v_muIpIsoTune.at(0).pfiso, eventWeightMu);
        if(v_muIpIsoTune.at(0).pfiso > 0.06) std::cout << "v_muIpIsoTune.at(0).pfiso 2 = " << v_muIpIsoTune.at(0).pfiso << std::endl;
        h_trailing_mu_Iso_Tune_MuMu->Fill(v_muIpIsoTune.at(1).pfiso, eventWeightMu);
        h_leading_mu_d0_Tune_MuMu->Fill(fabs(v_muIpIsoTune.at(0).d0), eventWeightMu);
        h_trailing_mu_d0_Tune_MuMu->Fill(fabs(v_muIpIsoTune.at(1).d0), eventWeightMu);
        h_leading_mu_dz_Tune_MuMu->Fill(fabs(v_muIpIsoTune.at(0).dz), eventWeightMu);
        h_trailing_mu_dz_Tune_MuMu->Fill(fabs(v_muIpIsoTune.at(1).dz), eventWeightMu);
      }
    }
    else if(v_muIpIsoTune.size()==1 and v_elIpIsoTune.size()==1 and v_selectedJets.size()==1)
    {   
      double eventWeightElMu = 1.0;     
      if(MCSample=="Signal")eventWeightElMu = v_muIpIsoTune.at(0).recoEW*v_elIpIsoTune.at(0).recoEW*(1 - (1 - v_muIpIsoTune.at(0).triggerEW)*(1 - v_elIpIsoTune.at(0).triggerEW));
      else if(MCSample=="MC")eventWeightElMu = v_muIpIsoTune.at(0).recoEW*v_elIpIsoTune.at(0).recoEW*(1 - (1 - v_muIpIsoTune.at(0).triggerEW)*(1 - v_elIpIsoTune.at(0).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muIpIsoTune.at(0).charge*v_elIpIsoTune.at(0).charge==1 and v_muIpIsoTune.at(0).pt > 30.0 and v_elIpIsoTune.at(0).pt > 30.0)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muIpIsoTune.at(0).pt, v_muIpIsoTune.at(0).eta, v_muIpIsoTune.at(0).phi, MUON_MASS);
        TLorentzVector el1 = fillTLorentzVector(v_elIpIsoTune.at(0).pt, v_elIpIsoTune.at(0).eta, v_elIpIsoTune.at(0).phi, ELECTRON_MASS);
        fillElMuHistCollection(elmuHistCut13, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
        h_leading_mu_Iso_Tune_ElMu->Fill(v_muIpIsoTune.at(0).pfiso, eventWeightElMu);
        h_leading_el_Iso_Tune_ElMu->Fill(v_elIpIsoTune.at(0).pfiso, eventWeightElMu);
        h_leading_mu_d0_Tune_ElMu->Fill(fabs(v_muIpIsoTune.at(0).d0), eventWeightElMu);
        h_leading_el_d0_Tune_ElMu->Fill(fabs(v_elIpIsoTune.at(0).d0), eventWeightElMu);
        h_leading_mu_dz_Tune_ElMu->Fill(fabs(v_muIpIsoTune.at(0).dz), eventWeightElMu);
        h_leading_el_dz_Tune_ElMu->Fill(fabs(v_elIpIsoTune.at(0).dz), eventWeightElMu);
      }
    }
    else if(v_elIpIsoTune.size()==2 and v_muons.size()==0 and v_selectedJets.size()==1)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal")eventWeightEl = v_elIpIsoTune.at(0).recoEW*v_elIpIsoTune.at(1).recoEW*(1 - (1 - v_elIpIsoTune.at(0).triggerEW)*(1 - v_elIpIsoTune.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightEl = v_elIpIsoTune.at(0).recoEW*v_elIpIsoTune.at(1).recoEW*(1 - (1 - v_elIpIsoTune.at(0).triggerEW)*(1 - v_elIpIsoTune.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_elIpIsoTune.at(0).charge*v_elIpIsoTune.at(1).charge==1 and v_elIpIsoTune.at(0).pt > 30.0 and v_elIpIsoTune.at(1).pt > 30.0)
      { 
        TLorentzVector el1 = fillTLorentzVector(v_elIpIsoTune.at(0).pt, v_elIpIsoTune.at(0).eta, v_elIpIsoTune.at(0).phi, ELECTRON_MASS);
        TLorentzVector el2 = fillTLorentzVector(v_elIpIsoTune.at(1).pt, v_elIpIsoTune.at(1).eta, v_elIpIsoTune.at(1).phi, ELECTRON_MASS);
        fillElHistCollection(elelHistCut12, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
        h_leading_el_Iso_Tune_ElEl->Fill(v_elIpIsoTune.at(0).pfiso, eventWeightEl);
        h_trailing_el_Iso_Tune_ElEl->Fill(v_elIpIsoTune.at(1).pfiso, eventWeightEl);
        h_leading_el_d0_Tune_ElEl->Fill(fabs(v_elIpIsoTune.at(0).d0), eventWeightEl);
        h_trailing_el_d0_Tune_ElEl->Fill(fabs(v_elIpIsoTune.at(1).d0), eventWeightEl);
        h_leading_el_dz_Tune_ElEl->Fill(fabs(v_elIpIsoTune.at(0).dz), eventWeightEl);
        h_trailing_el_dz_Tune_ElEl->Fill(fabs(v_elIpIsoTune.at(1).dz), eventWeightEl);
      }
    }

    //begin of debug plots
    if(v_muplotIdIso.size()==2 and v_electrons.size()==0)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muplotIdIso.at(0).recoEW*v_muplotIdIso.at(1).recoEW*(1 - (1 - v_muplotIdIso.at(0).triggerEW)*(1 - v_muplotIdIso.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muplotIdIso.at(0).recoEW*v_muplotIdIso.at(1).recoEW*(1 - (1 - v_muplotIdIso.at(0).triggerEW)*(1 - v_muplotIdIso.at(1).triggerEW))*eventWeight;
      if(v_muplotIdIso.at(0).charge*v_muplotIdIso.at(1).charge==1 and v_muplotIdIso.at(0).pt > 30.0 and v_muplotIdIso.at(1).pt > 30.0) 
      {
        h_leading_mu_IdpfIso_MuMu->Fill(v_muplotIdIso.at(0).pfiso, eventWeightMu);
        h_leading_mu_IdtrkIso_MuMu->Fill(v_muplotIdIso.at(0).trkiso/v_muplotIdIso.at(0).pt, eventWeightMu);
        h_trailing_mu_IdpfIso_MuMu->Fill(v_muplotIdIso.at(1).pfiso, eventWeightMu);
        h_trailing_mu_IdtrkIso_MuMu->Fill(v_muplotIdIso.at(1).trkiso/v_muplotIdIso.at(1).pt, eventWeightMu);
      }
    }
    else if(v_muplotIdIso.size()==1 and v_elplotIdIso.size()==1)
    { 
      double eventWeightElMu = 1.0;
      if(MCSample=="Signal")eventWeightElMu = v_elplotIdIso.at(0).recoEW*v_muplotIdIso.at(0).recoEW*(1 - (1 - v_elplotIdIso.at(0).triggerEW)*(1 - v_muplotIdIso.at(0).triggerEW));
      else if(MCSample=="MC")eventWeightElMu = v_elplotIdIso.at(0).recoEW*v_muplotIdIso.at(0).recoEW*(1 - (1 - v_elplotIdIso.at(0).triggerEW)*(1 - v_muplotIdIso.at(0).triggerEW))*eventWeight;
      if(v_elplotIdIso.at(0).charge*v_muplotIdIso.at(0).charge==1 and v_elplotIdIso.at(0).pt > 30.0 and v_muplotIdIso.at(0).pt > 30.0)
      {
        h_leading_mu_IdpfIso_ElMu->Fill(v_muplotIdIso.at(0).pfiso, eventWeightElMu);
        h_leading_mu_IdtrkIso_ElMu->Fill(v_muplotIdIso.at(0).trkiso/v_muplotIdIso.at(0).pt, eventWeightElMu);
        h_leading_el_IdpfIso_ElMu->Fill(v_elplotIdIso.at(0).pfiso, eventWeightElMu);
        h_leading_el_IdtrkIso_ElMu->Fill(v_elplotIdIso.at(0).trkiso/v_elplotIdIso.at(0).pt, eventWeightElMu);
      }
    }
    else if(v_elplotIdIso.size()==2 and v_muons.size()==0)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal")eventWeightEl = v_elplotIdIso.at(0).recoEW*v_elplotIdIso.at(1).recoEW*(1 - (1 - v_elplotIdIso.at(0).triggerEW)*(1 - v_elplotIdIso.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightEl = v_elplotIdIso.at(0).recoEW*v_elplotIdIso.at(1).recoEW*(1 - (1 - v_elplotIdIso.at(0).triggerEW)*(1 - v_elplotIdIso.at(1).triggerEW))*eventWeight;
      if(v_elplotIdIso.at(0).charge*v_elplotIdIso.at(1).charge==1 and v_elplotIdIso.at(0).pt > 30.0 and v_elplotIdIso.at(1).pt > 30.0)
      {
        h_leading_el_IdpfIso_ElEl->Fill(v_elplotIdIso.at(0).pfiso, eventWeightEl);
        h_leading_el_IdtrkIso_ElEl->Fill(v_elplotIdIso.at(0).trkiso/v_elplotIdIso.at(0).pt, eventWeightEl);
        h_trailing_el_IdpfIso_ElEl->Fill(v_elplotIdIso.at(1).pfiso, eventWeightEl);
        h_trailing_el_IdtrkIso_ElEl->Fill(v_elplotIdIso.at(1).trkiso/v_elplotIdIso.at(1).pt, eventWeightEl);
      }
    }
    
    if(v_muplotIso.size()==2 and v_electrons.size()==0)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muplotIso.at(0).recoEW*v_muplotIso.at(1).recoEW*(1 - (1 - v_muplotIso.at(0).triggerEW)*(1 - v_muplotIso.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muplotIso.at(0).recoEW*v_muplotIso.at(1).recoEW*(1 - (1 - v_muplotIso.at(0).triggerEW)*(1 - v_muplotIso.at(1).triggerEW))*eventWeight;
      if(v_muplotIso.at(0).charge*v_muplotIso.at(1).charge==1 and v_muplotIso.at(0).pt > 30.0 and v_muplotIso.at(1).pt > 30.0)
      {
        h_leading_mu_pfIso_MuMu->Fill(v_muplotIso.at(0).pfiso, eventWeightMu);
        h_leading_mu_trkIso_MuMu->Fill(v_muplotIso.at(0).trkiso/v_muplotIso.at(0).pt, eventWeightMu);
        h_trailing_mu_pfIso_MuMu->Fill(v_muplotIso.at(1).pfiso, eventWeightMu);
        h_trailing_mu_trkIso_MuMu->Fill(v_muplotIso.at(1).trkiso/v_muplotIso.at(1).pt, eventWeightMu);
      }
    }
    else if(v_muplotIso.size()==1 and v_elplotIso.size()==1)
    {
      double eventWeightElMu = 1.0;
      if(MCSample=="Signal")eventWeightElMu = v_elplotIso.at(0).recoEW*v_muplotIso.at(0).recoEW*(1 - (1 - v_elplotIso.at(0).triggerEW)*(1 - v_muplotIso.at(0).triggerEW));
      else if(MCSample=="MC")eventWeightElMu = v_elplotIso.at(0).recoEW*v_muplotIso.at(0).recoEW*(1 - (1 - v_elplotIso.at(0).triggerEW)*(1 - v_muplotIso.at(0).triggerEW))*eventWeight;
      if(v_elplotIso.at(0).charge*v_muplotIso.at(0).charge==1 and v_elplotIso.at(0).pt > 30.0 and v_muplotIso.at(0).pt > 30.0)
      {
        h_leading_mu_pfIso_ElMu->Fill(v_muplotIso.at(0).pfiso, eventWeightElMu);
        h_leading_mu_trkIso_ElMu->Fill(v_muplotIso.at(0).trkiso/v_muplotIso.at(0).pt, eventWeightElMu);
        h_leading_el_pfIso_ElMu->Fill(v_elplotIso.at(0).pfiso, eventWeightElMu);
        h_leading_el_trkIso_ElMu->Fill(v_elplotIso.at(0).trkiso/v_elplotIso.at(0).pt, eventWeightElMu);
      }
    }
    else if(v_elplotIso.size()==2 and v_muons.size()==0)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal")eventWeightEl = v_elplotIso.at(0).recoEW*v_elplotIso.at(1).recoEW*(1 - (1 - v_elplotIso.at(0).triggerEW)*(1 - v_elplotIso.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightEl = v_elplotIso.at(0).recoEW*v_elplotIso.at(1).recoEW*(1 - (1 - v_elplotIso.at(0).triggerEW)*(1 - v_elplotIso.at(1).triggerEW))*eventWeight;
      if(v_elplotIso.at(0).charge*v_elplotIso.at(1).charge==1 and v_elplotIso.at(0).pt > 30.0 and v_elplotIso.at(1).pt > 30.0)
      {
        h_leading_el_pfIso_ElEl->Fill(v_elplotIso.at(0).pfiso, eventWeightEl);
        h_leading_el_trkIso_ElEl->Fill(v_elplotIso.at(0).trkiso/v_elplotIso.at(0).pt, eventWeightEl);
        h_trailing_el_pfIso_ElEl->Fill(v_elplotIso.at(1).pfiso, eventWeightEl);
        h_trailing_el_trkIso_ElEl->Fill(v_elplotIso.at(1).trkiso/v_elplotIso.at(1).pt, eventWeightEl);
      }
    }
    //Ip plots
    if(v_muplotIdIsoIp.size()==2 and v_electrons.size()==0)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muplotIdIsoIp.at(0).recoEW*v_muplotIdIsoIp.at(1).recoEW*(1 - (1 - v_muplotIdIsoIp.at(0).triggerEW)*(1 - v_muplotIdIsoIp.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muplotIdIsoIp.at(0).recoEW*v_muplotIdIsoIp.at(1).recoEW*(1 - (1 - v_muplotIdIsoIp.at(0).triggerEW)*(1 - v_muplotIdIsoIp.at(1).triggerEW))*eventWeight;
      if(v_muplotIdIsoIp.at(0).charge*v_muplotIdIsoIp.at(1).charge==1 and v_muplotIdIsoIp.at(0).pt > 30.0 and v_muplotIdIsoIp.at(1).pt > 30.0)
      {
        h_leading_mu_d0_IdIsoIp_MuMu->Fill(v_muplotIdIsoIp.at(0).d0,  eventWeightMu);
        h_trailing_mu_d0_IdIsoIp_MuMu->Fill(v_muplotIdIsoIp.at(1).d0, eventWeightMu);
        h_leading_mu_dz_IdIsoIp_MuMu->Fill(v_muplotIdIsoIp.at(0).dz, eventWeightMu);
        h_trailing_mu_dz_IdIsoIp_MuMu->Fill(v_muplotIdIsoIp.at(1).dz,eventWeightMu);
      }
    }
    else if(v_muplotIdIsoIp.size()==1 and v_elplotIdIsoIp.size()==1)
    {
      double eventWeightElMu = 1.0;
      if(MCSample=="Signal")eventWeightElMu = v_elplotIdIsoIp.at(0).recoEW*v_muplotIdIsoIp.at(0).recoEW*(1 - (1 - v_elplotIdIsoIp.at(0).triggerEW)*(1 - v_muplotIdIsoIp.at(0).triggerEW));
      else if(MCSample=="MC")eventWeightElMu = v_elplotIdIsoIp.at(0).recoEW*v_muplotIdIsoIp.at(0).recoEW*(1 - (1 - v_elplotIdIsoIp.at(0).triggerEW)*(1 - v_muplotIdIsoIp.at(0).triggerEW))*eventWeight;
      if(v_elplotIdIsoIp.at(0).charge*v_muplotIdIsoIp.at(0).charge==1 and v_elplotIdIsoIp.at(0).pt > 30.0 and v_muplotIdIsoIp.at(0).pt > 30.0)
      {
        h_leading_mu_d0_IdIsoIp_ElMu->Fill(v_muplotIdIsoIp.at(0).d0, eventWeightElMu);
        h_leading_el_d0_IdIsoIp_ElMu->Fill(v_elplotIdIsoIp.at(0).d0, eventWeightElMu);
        h_leading_mu_dz_IdIsoIp_ElMu->Fill(v_muplotIdIsoIp.at(0).dz, eventWeightElMu);
        h_leading_el_dz_IdIsoIp_ElMu->Fill(v_elplotIdIsoIp.at(0).dz, eventWeightElMu);
      }
    }
    else if(v_elplotIdIsoIp.size()==2 and v_muons.size()==0)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal")eventWeightEl = v_elplotIdIsoIp.at(0).recoEW*v_elplotIdIsoIp.at(1).recoEW*(1 - (1 - v_elplotIdIsoIp.at(0).triggerEW)*(1 - v_elplotIdIsoIp.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightEl = v_elplotIdIsoIp.at(0).recoEW*v_elplotIdIsoIp.at(1).recoEW*(1 - (1 - v_elplotIdIsoIp.at(0).triggerEW)*(1 - v_elplotIdIsoIp.at(1).triggerEW))*eventWeight;
      if(v_elplotIdIsoIp.at(0).charge*v_elplotIdIsoIp.at(1).charge==1 and v_elplotIdIsoIp.at(0).pt > 30.0 and v_elplotIdIsoIp.at(1).pt > 30.0)
      {
        h_leading_el_d0_IdIsoIp_ElEl->Fill(v_elplotIdIsoIp.at(0).d0, eventWeightEl);
        h_trailing_el_d0_IdIsoIp_ElEl->Fill(v_elplotIdIsoIp.at(1).d0, eventWeightEl);
        h_leading_el_dz_IdIsoIp_ElEl->Fill(v_elplotIdIsoIp.at(0).dz, eventWeightEl);
        h_trailing_el_dz_IdIsoIp_ElEl->Fill(v_elplotIdIsoIp.at(1).dz, eventWeightEl);
      }   
    }

    if(v_muplotIdIp.size()==2 and v_electrons.size()==0)
    { 
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muplotIdIp.at(0).recoEW*v_muplotIdIp.at(1).recoEW*(1 - (1 - v_muplotIdIp.at(0).triggerEW)*(1 - v_muplotIdIp.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muplotIdIp.at(0).recoEW*v_muplotIdIp.at(1).recoEW*(1 - (1 - v_muplotIdIp.at(0).triggerEW)*(1 - v_muplotIdIp.at(1).triggerEW))*eventWeight;
      if(v_muplotIdIp.at(0).charge*v_muplotIdIp.at(1).charge==1 and v_muplotIdIp.at(0).pt > 30.0 and v_muplotIdIp.at(1).pt > 30.0)
      { 
        h_leading_mu_d0_IdIp_MuMu->Fill(v_muplotIdIp.at(0).d0, eventWeightMu);
        h_trailing_mu_d0_IdIp_MuMu->Fill(v_muplotIdIp.at(1).d0, eventWeightMu);
        h_leading_mu_dz_IdIp_MuMu->Fill(v_muplotIdIp.at(0).dz, eventWeightMu);
        h_trailing_mu_dz_IdIp_MuMu->Fill(v_muplotIdIp.at(1).dz, eventWeightMu);
      }
    }
    else if(v_muplotIdIp.size()==1 and v_elplotIdIp.size()==1)
    { 
      double eventWeightElMu = 1.0;
      if(MCSample=="Signal")eventWeightElMu = v_elplotIdIp.at(0).recoEW*v_muplotIdIp.at(0).recoEW*(1 - (1 - v_elplotIdIp.at(0).triggerEW)*(1 - v_muplotIdIp.at(0).triggerEW));
      else if(MCSample=="MC")eventWeightElMu = v_elplotIdIp.at(0).recoEW*v_muplotIdIp.at(0).recoEW*(1 - (1 - v_elplotIdIp.at(0).triggerEW)*(1 - v_muplotIdIp.at(0).triggerEW))*eventWeight;
      if(v_elplotIdIp.at(0).charge*v_muplotIdIp.at(0).charge==1 and v_elplotIdIp.at(0).pt > 30.0 and v_muplotIdIp.at(0).pt > 30.0)
      { 
        h_leading_mu_d0_IdIp_ElMu->Fill(v_muplotIdIp.at(0).d0, eventWeightElMu);
        h_leading_el_d0_IdIp_ElMu->Fill(v_elplotIdIp.at(0).d0, eventWeightElMu);
        h_leading_mu_dz_IdIp_ElMu->Fill(v_muplotIdIp.at(0).dz, eventWeightElMu);
        h_leading_el_dz_IdIp_ElMu->Fill(v_elplotIdIp.at(0).dz, eventWeightElMu);
      }
    }
    else if(v_elplotIdIp.size()==2 and v_muons.size()==0)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal")eventWeightEl = v_elplotIdIp.at(0).recoEW*v_elplotIdIp.at(1).recoEW*(1 - (1 - v_elplotIdIp.at(0).triggerEW)*(1 - v_elplotIdIp.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightEl = v_elplotIdIp.at(0).recoEW*v_elplotIdIp.at(1).recoEW*(1 - (1 - v_elplotIdIp.at(0).triggerEW)*(1 - v_elplotIdIp.at(1).triggerEW))*eventWeight;
      if(v_elplotIdIp.at(0).charge*v_elplotIdIp.at(1).charge==1 and v_elplotIdIp.at(0).pt > 30.0 and v_elplotIdIp.at(1).pt > 30.0)
      {
        h_leading_el_d0_IdIp_ElEl->Fill(v_elplotIdIp.at(0).d0, eventWeightEl);
        h_trailing_el_d0_IdIp_ElEl->Fill(v_elplotIdIp.at(1).d0, eventWeightEl);
        h_leading_el_dz_IdIp_ElEl->Fill(v_elplotIdIp.at(0).dz, eventWeightEl);
        h_trailing_el_dz_IdIp_ElEl->Fill(v_elplotIdIp.at(1).dz, eventWeightEl);
      }
    }
    //end of debug plots
 
    if(v_muons.size() == 2 and v_electrons.size()==0)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal")eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW));
      else if(MCSample=="MC")eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muons.at(0).charge*v_muons.at(1).charge==1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        fillMuHistCollection(mumuHistCut1, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);  
        cut1_mumu+=eventWeightMu;
        h_TotalEvents_MuMu->Fill(1, eventWeightMu);
        h_nPU_MuMu->Fill(nPU, eventWeightMu);
        h_nPV_MuMu->Fill(nPV, eventWeightMu);
        if(v_looseElectrons.size()==0 and v_looseMuons.size()==0) 
        {
          fillMuHistCollection(mumuHistCut2, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);  
          cut2_mumu+=eventWeightMu;
          h_TotalEvents_MuMu->Fill(2, eventWeightMu);
          if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
          {
            cut3_mumu+=eventWeightMu;
            if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
            if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
            if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi(); 
            if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
            if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV;
            fillMuHistCollection(mumuHistCut3, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
            h_TotalEvents_MuMu->Fill(3, eventWeightMu);
            if(v_selectedBJets.size()==0)
            {
              cut4_mumu+=eventWeightMu;
              fillMuHistCollection(mumuHistCut4, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
              double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
              double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
              h_TotalEvents_MuMu->Fill(4, eventWeightMu);
              if(invMassJJ > 60 and invMassJJ < 100) 
              {
                cut5_mumu+=eventWeightMu;
                fillMuHistCollection(mumuHistCut5, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                h_TotalEvents_MuMu->Fill(5, eventWeightMu);
                if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                {
                  cut6_mumu+=eventWeightMu;
                  fillMuHistCollection(mumuHistCut6, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);   
                  h_TotalEvents_MuMu->Fill(6, eventWeightMu);
                  if(met > 40.0) 
                  {
                    cut7_mumu+=eventWeightMu;
                    fillMuHistCollection(mumuHistCut7, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                    h_TotalEvents_MuMu->Fill(7, eventWeightMu);
                    if((mu1+mu2).M() > 40.0) 
                    {
                      cut8_mumu+=eventWeightMu;
                      fillMuHistCollection(mumuHistCut8, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
                      h_TotalEvents_MuMu->Fill(8, eventWeightMu);
                    }//Mll > 40.0
                  }//met > 40.0 
                }//Leading invmass and delta eta
              }//closest invmass
            }//b-veto
          }//2 jets with pT > 30
        }//loose lepton veto
      }//exactly 2 muons with pT > 30.0
    }//exactly 2 muons
    else if(v_electrons.size() == 1 and v_muons.size()==1)
    {
      double eventWeightElMu = 1.0; 
      if(MCSample=="Signal") eventWeightElMu = v_electrons.at(0).recoEW*v_muons.at(0).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_muons.at(0).triggerEW));
      else if(MCSample=="MC") eventWeightElMu =  v_electrons.at(0).recoEW*v_muons.at(0).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_muons.at(0).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_electrons.at(0).charge*v_muons.at(0).charge==1 and v_electrons.at(0).pt > 30.0 and v_muons.at(0).pt > 30.0)
      {
        TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector v_met;
        v_met.SetPtEtaPhiM(met, 0.0, metPhi, 0.0);
        fillElMuHistCollection(elmuHistCut1, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
        cut1_elmu+=eventWeightElMu;
        h_TotalEvents_ElMu->Fill(1, eventWeightElMu);
        h_nPU_ElMu->Fill(nPU, eventWeightElMu);
        h_nPV_ElMu->Fill(nPV, eventWeightElMu);
        if(v_looseElectrons.size()==0 and v_looseMuons.size()==0)
        {
          fillElMuHistCollection(elmuHistCut2, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
          cut2_elmu+=eventWeightElMu;
          h_TotalEvents_ElMu->Fill(2, eventWeightElMu);
          if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
          {
            cut3_elmu+=eventWeightElMu;
            if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
            if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
            if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
            if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
            if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV;
            fillElMuHistCollection(elmuHistCut3, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
            h_TotalEvents_ElMu->Fill(3, eventWeightElMu);
            if(v_selectedBJets.size()==0)
            {
              cut4_elmu+=eventWeightElMu;
              fillElMuHistCollection(elmuHistCut4, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
              double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
              double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
              h_TotalEvents_ElMu->Fill(4, eventWeightElMu);
              if(invMassJJ > 60 and invMassJJ < 100)
              {
                cut5_elmu+=eventWeightElMu;
                fillElMuHistCollection(elmuHistCut5, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                h_TotalEvents_ElMu->Fill(5, eventWeightElMu);
                if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                {
                  cut6_elmu+=eventWeightElMu;
                  fillElMuHistCollection(elmuHistCut6, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
                  h_TotalEvents_ElMu->Fill(6, eventWeightElMu);
                  if(met > 40.0)
                  {
                    cut7_elmu+=eventWeightElMu;
                    fillElMuHistCollection(elmuHistCut7, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
                    h_TotalEvents_ElMu->Fill(7, eventWeightElMu);
                    if((el1+mu1).M() > 40.0)
                    {
                      cut8_elmu+=eventWeightElMu;
                      fillElMuHistCollection(elmuHistCut8, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
                       h_TotalEvents_ElMu->Fill(8, eventWeightElMu);
                       double mt1 = (mu1+v_met).Mt();
                       double mt2 = (el1+v_met).Mt();
                       double mtMax = std::max(mt1, mt2);
                       if(mtMax > 90.0)
                       {
                         cut9_elmu+=eventWeightElMu;
                         fillElMuHistCollection(elmuHistCut12, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
                         h_TotalEvents_ElMu->Fill(9, eventWeightElMu);
                       }//mtmax > 90.0 
                    }//Mll > 40.0
                  }//met > 40.0 
                }//Leading invmass and delta eta
              }//closest invmass
            }//b-veto
          }//2 jets with pT > 30
        }//loose lepton veto
      }//exactly 1 electron with pT > 30.0 and 1 muon with pT > 30.0
    }//exactly 2 electrons
    else if(v_electrons.size() == 2 and v_muons.size()==0)
    {
        double eventWeightEl = 1.0;
        if(MCSample=="MC") eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW))*eventWeight;
        else if(MCSample=="Signal") eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW));
        double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
        jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
        if(v_electrons.at(0).charge*v_electrons.at(1).charge==1 and v_electrons.at(0).pt > 30.0 and v_electrons.at(1).pt > 30.0)
        {
          TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
          TLorentzVector el2 = fillTLorentzVector(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, ELECTRON_MASS);
          if(not ((el1+el2).M() < 80.0 or (el1+el2).M() > 100.0)) continue;
          fillElHistCollection(elelHistCut1, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
          cut1_elel+=eventWeightEl;
          h_TotalEvents_ElEl->Fill(1, eventWeightEl);
          h_nPU_ElEl->Fill(nPU, eventWeightEl);
          h_nPV_ElEl->Fill(nPV, eventWeightEl);
          if(v_looseElectrons.size()==0 and v_looseMuons.size()==0)
          {
            fillElHistCollection(elelHistCut2, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
            cut2_elel+=eventWeightEl;
            h_TotalEvents_ElEl->Fill(2, eventWeightEl);
            if(v_selectedJets.size()>=2 and fabs(v_selectedJets.at(0).JetLV.Eta()) < 2.5 and fabs(v_selectedJets.at(1).JetLV.Eta()) < 2.5)
            {
              cut3_elel+=eventWeightEl;
              if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
              if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
              if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
              if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
              if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV; 
              fillElHistCollection(elelHistCut3, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
              h_TotalEvents_ElEl->Fill(3, eventWeightEl);
              if(v_selectedBJets.size()==0)
              {
                cut4_elel+=eventWeightEl;
                fillElHistCollection(elelHistCut4, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                h_TotalEvents_ElEl->Fill(4, eventWeightEl);
                if(invMassJJ > 60 and invMassJJ < 100)
                {
                  cut5_elel+=eventWeightEl;
                  fillElHistCollection(elelHistCut5, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                  double invMassLeadingJJ = (v_selectedJets.at(0).JetLV + v_selectedJets.at(1).JetLV).M();
                  double deltaEta = fabs(v_selectedJets.at(0).JetLV.Eta()-v_selectedJets.at(1).JetLV.Eta());
                  h_TotalEvents_ElEl->Fill(5, eventWeightEl);
                  if(invMassLeadingJJ < 400.0 and deltaEta < 1.5)
                  {
                    cut6_elel+=eventWeightEl;
                    fillElHistCollection(elelHistCut6, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                    h_TotalEvents_ElEl->Fill(6, eventWeightEl);
                    if(met > 40.0)
                    {
                      cut7_elel+=eventWeightEl;
                      fillElHistCollection(elelHistCut7, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                      h_TotalEvents_ElEl->Fill(7, eventWeightEl);
                      if((el1+el2).M() > 40.0)
                      {
                        cut8_elel+=eventWeightEl;
                        fillElHistCollection(elelHistCut8, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
                        h_TotalEvents_ElEl->Fill(8, eventWeightEl);
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
    bool mutlv2jbV = mutlv2j and (v_selectedBJets.size()==0);

    //ttbar control region
    if(v_muons.size() >= 2 and v_electrons.size()==0)
    {
      double eventWeightMu = 1.0;
      if(MCSample=="Signal") eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW));
      else if(MCSample=="MC") eventWeightMu = v_muons.at(0).recoEW*v_muons.at(1).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_muons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muons.at(0).charge*v_muons.at(1).charge==-1 and v_muons.at(0).pt > 30.0 and v_muons.at(1).pt > 30.0 and met > 30)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector mu2 = fillTLorentzVector(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, MUON_MASS);
        if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
        if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
        if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
        if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
        if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV;
        fillMuHistCollection(mumuHistCut9, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
        if(v_selectedBJets.size() > 0 and met > 50 and v_selectedJets.size()>=2)
        {
          fillMuHistCollection(mumuHistCut10, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), bjet2pt, bjet2eta, bjet2phi, eventWeightMu);
          if(v_selectedBJets.size() > 1)
          {
            fillMuHistCollection(mumuHistCut11, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), (mu1+mu2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), v_selectedBJets.at(1).JetLV.Pt(), v_selectedBJets.at(1).JetLV.Eta(), v_selectedBJets.at(1).JetLV.Phi(), eventWeightMu);
          }//2b-jets  
        }//1b-jet
      }//OS muons
    }//2 muons
    if(v_muons.size() >= 1 and v_electrons.size() >= 1)
    {
      double eventWeightElMu = 1.0;
      if(MCSample=="Signal") eventWeightElMu = v_muons.at(0).recoEW*v_electrons.at(0).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_electrons.at(0).triggerEW));
      else if(MCSample=="MC") eventWeightElMu = v_muons.at(0).recoEW*v_electrons.at(0).recoEW*(1 - (1 - v_muons.at(0).triggerEW)*(1 - v_electrons.at(0).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_muons.at(0).charge*v_electrons.at(0).charge==-1 and v_muons.at(0).pt > 30.0 and v_electrons.at(0).pt > 30.0 and met > 30)
      {
        TLorentzVector mu1 = fillTLorentzVector(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, MUON_MASS);
        TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
        if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
        if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
        if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
        if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
        if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV;
        fillElMuHistCollection(elmuHistCut9, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
        if(v_selectedBJets.size() > 0 and met > 50 and v_selectedJets.size()>=2)
        {
          fillElMuHistCollection(mumuHistCut10, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), bjet2pt, bjet2eta, bjet2phi, eventWeightElMu);
          if(v_selectedBJets.size() > 1)
          {
            fillElMuHistCollection(mumuHistCut11, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), (el1+mu1).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), v_selectedBJets.at(1).JetLV.Pt(), v_selectedBJets.at(1).JetLV.Eta(), v_selectedBJets.at(1).JetLV.Phi(), eventWeightElMu);
          }//2b-jets
        }//1b-jet
      }//OS muons
    }//2 muons
    else if(v_electrons.size() >= 2 and v_muons.size()==0)
    {
      double eventWeightEl = 1.0;
      if(MCSample=="Signal") eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW));
      else if (MCSample=="MC") eventWeightEl = v_electrons.at(0).recoEW*v_electrons.at(1).recoEW*(1 - (1 - v_electrons.at(0).triggerEW)*(1 - v_electrons.at(1).triggerEW))*eventWeight;
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      if(v_electrons.at(0).charge*v_electrons.at(1).charge==-1 and v_electrons.at(0).pt > 30.0 and v_electrons.at(1).pt > 30.0 and met > 30)
      {
        TLorentzVector el1 = fillTLorentzVector(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, ELECTRON_MASS);
        TLorentzVector el2 = fillTLorentzVector(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, ELECTRON_MASS);
        if(v_selectedBJets.size()>0) bjet1pt = v_selectedBJets.at(0).JetLV.Pt();
        if(v_selectedBJets.size()>0) bjet1eta = v_selectedBJets.at(0).JetLV.Eta();
        if(v_selectedBJets.size()>0) bjet1phi = v_selectedBJets.at(0).JetLV.Phi();
        if(v_selectedBJets.size()>0) bjet1csv = v_selectedBJets.at(0).BTag_CSV;
        if(v_selectedBJets.size()>1) bjet2csv = v_selectedBJets.at(1).BTag_CSV;
        fillElHistCollection(elelHistCut9, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
        if(v_selectedBJets.size() > 0 and met > 50 and v_selectedJets.size()>=2)
        {
          fillElHistCollection(elelHistCut10, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), bjet2pt, bjet2eta, bjet2phi, eventWeightEl);
          if(v_selectedBJets.size() > 1)
          {
            fillElHistCollection(elelHistCut11, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), (el1+el2).M(), met, v_selectedJets.size(), v_selectedBJets.size(), v_selectedJets.at(0).JetLV.Pt(), v_selectedJets.at(1).JetLV.Pt(), v_selectedJets.at(0).JetLV.Eta(), v_selectedJets.at(1).JetLV.Eta(), v_selectedJets.at(0).JetLV.Phi(), v_selectedJets.at(1).JetLV.Phi(), bjet1csv, bjet2csv, v_selectedBJets.at(0).JetLV.Pt(), v_selectedBJets.at(0).JetLV.Eta(), v_selectedBJets.at(0).JetLV.Phi(), v_selectedBJets.at(1).JetLV.Pt(), v_selectedBJets.at(1).JetLV.Eta(), v_selectedBJets.at(1).JetLV.Phi(), eventWeightEl);
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
 
  std::cout << "1Electron+1Muon = " << cut1_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto = " << cut2_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets = " << cut3_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto = " << cut4_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto+ClosestInvMassWindow = " << cut5_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow = " << cut6_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut = " << cut7_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut+MllCut = " << cut8_elmu << std::endl;
  std::cout << "1Electron+1Muon+LooseLeptonVeto+2CentralJets+BVeto+ClosestInvMassWindow+LeadingInvMassWindow+MetCut+MllCut+MtMaxCut = " << cut9_elmu << std::endl;
 
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
  tFile->mkdir("2SSTLLV2JBVCMjjLMjjMetMllMtmax");
  tFile->mkdir("TTbarCR");
  tFile->mkdir("TTbarCR1b");
  tFile->mkdir("TTbarCR2b");
  tFile->mkdir("IpIsoTuning");
  tFile->cd("2SSTL");
  writeHistCollection(mumuHistCut1);
  writeHistCollection(elelHistCut1);
  writeHistCollection(elmuHistCut1);
  tFile->cd("2SSTLLV");
  writeHistCollection(mumuHistCut2);
  writeHistCollection(elelHistCut2);
  writeHistCollection(elmuHistCut2);
  tFile->cd("2SSTLLV2J");
  writeHistCollection(mumuHistCut3);
  writeHistCollection(elelHistCut3);
  writeHistCollection(elmuHistCut3);
  tFile->cd("2SSTLLV2JBV");
  writeHistCollection(mumuHistCut4);
  writeHistCollection(elelHistCut4);
  writeHistCollection(elmuHistCut4);
  tFile->cd("2SSTLLV2JBVCMjj");
  writeHistCollection(mumuHistCut5);
  writeHistCollection(elelHistCut5);
  writeHistCollection(elmuHistCut5);
  tFile->cd("2SSTLLV2JBVCMjjLMjj");
  writeHistCollection(mumuHistCut6);
  writeHistCollection(elelHistCut6);
  writeHistCollection(elmuHistCut6);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMet");
  writeHistCollection(mumuHistCut7);
  writeHistCollection(elelHistCut7);
  writeHistCollection(elmuHistCut7);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMetMll");
  writeHistCollection(mumuHistCut8);
  writeHistCollection(elelHistCut8);
  writeHistCollection(elmuHistCut8);
  tFile->cd("2SSTLLV2JBVCMjjLMjjMetMllMtmax");
  writeHistCollection(elmuHistCut12);
  tFile->cd("TTbarCR");
  writeHistCollection(mumuHistCut9);
  writeHistCollection(elelHistCut9);
  writeHistCollection(elmuHistCut9);
  tFile->cd("TTbarCR1b");
  writeHistCollection(mumuHistCut10);
  writeHistCollection(elelHistCut10);
  writeHistCollection(elmuHistCut10);
  tFile->cd("TTbarCR2b");
  writeHistCollection(mumuHistCut11);
  writeHistCollection(elelHistCut11);
  writeHistCollection(elmuHistCut11);
  tFile->cd("IpIsoTuning");
  writeHistCollection(mumuHistCut12);
  writeHistCollection(elelHistCut12);
  writeHistCollection(elmuHistCut13);
  tFile->cd();
  h_TotalEvents_ElEl->Write();
  h_TotalEvents_ElMu->Write();
  h_TotalEvents_MuMu->Write();
  h_nPU_MuMu->Write();
  h_nPU_ElMu->Write();
  h_nPU_ElEl->Write();
  h_nPV_MuMu->Write();
  h_nPV_ElMu->Write();
  h_nPV_ElEl->Write();
  h_leading_mu_IdpfIso_MuMu->Write();
  h_trailing_mu_IdpfIso_MuMu->Write();
  h_leading_mu_IdtrkIso_MuMu->Write();
  h_trailing_mu_IdtrkIso_MuMu->Write();
  h_leading_mu_IdpfIso_ElMu->Write();
  h_leading_el_IdpfIso_ElMu->Write();
  h_leading_mu_IdtrkIso_ElMu->Write();
  h_leading_el_IdtrkIso_ElMu->Write();
  h_leading_el_IdpfIso_ElEl->Write();
  h_trailing_el_IdpfIso_ElEl->Write();
  h_leading_el_IdtrkIso_ElEl->Write();
  h_trailing_el_IdtrkIso_ElEl->Write();
  
  h_leading_mu_pfIso_MuMu->Write();
  h_trailing_mu_pfIso_MuMu->Write();
  h_leading_mu_trkIso_MuMu->Write();
  h_trailing_mu_trkIso_MuMu->Write();
  h_leading_mu_pfIso_ElMu->Write();
  h_leading_el_pfIso_ElMu->Write();
  h_leading_mu_trkIso_ElMu->Write();
  h_leading_el_trkIso_ElMu->Write();
  h_leading_el_pfIso_ElEl->Write();
  h_trailing_el_pfIso_ElEl->Write();
  h_leading_el_trkIso_ElEl->Write();
  h_trailing_el_trkIso_ElEl->Write();
  
  h_leading_mu_d0_IdIsoIp_MuMu->Write();
  h_trailing_mu_d0_IdIsoIp_MuMu->Write();
  h_leading_mu_dz_IdIsoIp_MuMu->Write();
  h_trailing_mu_dz_IdIsoIp_MuMu->Write();

  h_leading_mu_d0_IdIsoIp_ElMu->Write();
  h_leading_el_d0_IdIsoIp_ElMu->Write();
  h_leading_mu_dz_IdIsoIp_ElMu->Write();
  h_leading_el_dz_IdIsoIp_ElMu->Write();

  h_leading_el_d0_IdIsoIp_ElEl->Write();
  h_trailing_el_d0_IdIsoIp_ElEl->Write();
  h_leading_el_dz_IdIsoIp_ElEl->Write();
  h_trailing_el_dz_IdIsoIp_ElEl->Write();

  h_leading_mu_d0_IdIp_MuMu->Write();
  h_trailing_mu_d0_IdIp_MuMu->Write();
  h_leading_mu_dz_IdIp_MuMu->Write();
  h_trailing_mu_dz_IdIp_MuMu->Write();

  h_leading_mu_d0_IdIp_ElMu->Write();
  h_leading_el_d0_IdIp_ElMu->Write();
  h_leading_mu_dz_IdIp_ElMu->Write();
  h_leading_el_dz_IdIp_ElMu->Write();

  h_leading_el_d0_IdIp_ElEl->Write();
  h_trailing_el_d0_IdIp_ElEl->Write();
  h_leading_el_dz_IdIp_ElEl->Write();
  h_trailing_el_dz_IdIp_ElEl->Write();
      
  h_leading_mu_Iso_Tune_MuMu->Write();
  h_trailing_mu_Iso_Tune_MuMu->Write();

  h_leading_mu_Iso_Tune_ElMu->Write();
  h_leading_el_Iso_Tune_ElMu->Write();

  h_leading_el_Iso_Tune_ElEl->Write();
  h_trailing_el_Iso_Tune_ElEl->Write();
  
  h_leading_mu_d0_Tune_MuMu->Write();
  h_trailing_mu_d0_Tune_MuMu->Write();

  h_leading_mu_dz_Tune_MuMu->Write();
  h_trailing_mu_dz_Tune_MuMu->Write();

  h_leading_mu_d0_Tune_ElMu->Write();
  h_leading_el_d0_Tune_ElMu->Write();

  h_leading_mu_dz_Tune_ElMu->Write();
  h_leading_el_dz_Tune_ElMu->Write();

  h_leading_el_d0_Tune_ElEl->Write();
  h_trailing_el_d0_Tune_ElEl->Write();

  h_leading_el_dz_Tune_ElEl->Write();
  h_trailing_el_dz_Tune_ElEl->Write();

  tFile->Close(); 
  std::cout<<"Wrote output file "<<histfilename<<std::endl; 


  return 0;

}
