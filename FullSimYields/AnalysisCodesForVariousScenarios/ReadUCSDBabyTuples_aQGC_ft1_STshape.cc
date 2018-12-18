#include "ReadUCSDBabyTuples.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include <fstream>

using std::string;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool sameVal(double a, double b)
{
   return fabs(a - b) < 1.000e-01;
}


int closestJet(const TLorentzVector& lep_p4, std::vector<TLorentzVector>& v_Jets, float dRmin, float maxAbsEta, TLorentzVector& jet_p4)
{
  bool returnValue = -1;
  for (unsigned int pfjidx=0; pfjidx<v_Jets.size(); ++pfjidx)
  {
    TLorentzVector tmp_jet_p4 = v_Jets.at(pfjidx);
    if (fabs(tmp_jet_p4.Eta()) < maxAbsEta)
    {
      float tmp_dRmin = tmp_jet_p4.DeltaR(lep_p4);
      if (tmp_dRmin < dRmin)
      {
        returnValue = pfjidx;
	dRmin = tmp_dRmin;
	jet_p4 = tmp_jet_p4;
      }
    }
  }
  return  returnValue;
}

int ReadUCSDBabyTuples_aQGC(std::string infile, std::string treeStr, std::string type, std::string sample, std::string outTree, std::string Trigger="None")
{
  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  bool debug=false;

  TFile *outputFile;
  TTree *outputTree;
  double ST_rwgt77_MuMu, ST_rwgt78_MuMu, ST_rwgt79_MuMu, ST_rwgt80_MuMu, ST_rwgt81_MuMu, ST_rwgt82_MuMu, ST_rwgt83_MuMu, ST_rwgt84_MuMu, ST_rwgt85_MuMu, ST_rwgt86_MuMu, ST_rwgt87_MuMu, ST_bkg_MuMu;
  double ST_rwgt77_ElMu, ST_rwgt78_ElMu, ST_rwgt79_ElMu, ST_rwgt80_ElMu, ST_rwgt81_ElMu, ST_rwgt82_ElMu, ST_rwgt83_ElMu, ST_rwgt84_ElMu, ST_rwgt85_ElMu, ST_rwgt86_ElMu, ST_rwgt87_ElMu, ST_bkg_ElMu;
  double ST_rwgt77_ElEl, ST_rwgt78_ElEl, ST_rwgt79_ElEl, ST_rwgt80_ElEl, ST_rwgt81_ElEl, ST_rwgt82_ElEl, ST_rwgt83_ElEl, ST_rwgt84_ElEl, ST_rwgt85_ElEl, ST_rwgt86_ElEl, ST_rwgt87_ElEl, ST_bkg_ElEl;
  double weight_rwgt77_MuMu, weight_rwgt78_MuMu, weight_rwgt79_MuMu, weight_rwgt80_MuMu, weight_rwgt81_MuMu, weight_rwgt82_MuMu, weight_rwgt83_MuMu, weight_rwgt84_MuMu, weight_rwgt85_MuMu, weight_rwgt86_MuMu, weight_rwgt87_MuMu, weight_bkg_MuMu;
  double weight_rwgt77_ElMu, weight_rwgt78_ElMu, weight_rwgt79_ElMu, weight_rwgt80_ElMu, weight_rwgt81_ElMu, weight_rwgt82_ElMu, weight_rwgt83_ElMu, weight_rwgt84_ElMu, weight_rwgt85_ElMu, weight_rwgt86_ElMu, weight_rwgt87_ElMu, weight_bkg_ElMu;
  double weight_rwgt77_ElEl, weight_rwgt78_ElEl, weight_rwgt79_ElEl, weight_rwgt80_ElEl, weight_rwgt81_ElEl, weight_rwgt82_ElEl, weight_rwgt83_ElEl, weight_rwgt84_ElEl, weight_rwgt85_ElEl, weight_rwgt86_ElEl, weight_rwgt87_ElEl, weight_bkg_ElEl; 

  std::string outputtreename=(outTree+".root").c_str();
  outputFile = new TFile((outputtreename).c_str(),"RECREATE");
  outputTree=new TTree("ST_Tree", "ST_Tree");
  outputTree->Branch("ST_rwgt77_MuMu", &ST_rwgt77_MuMu, "ST_rwgt77_MuMu/D");
  outputTree->Branch("ST_rwgt78_MuMu", &ST_rwgt78_MuMu, "ST_rwgt78_MuMu/D");
  outputTree->Branch("ST_rwgt79_MuMu", &ST_rwgt79_MuMu, "ST_rwgt79_MuMu/D");
  outputTree->Branch("ST_rwgt80_MuMu", &ST_rwgt80_MuMu, "ST_rwgt80_MuMu/D");
  outputTree->Branch("ST_rwgt81_MuMu", &ST_rwgt81_MuMu, "ST_rwgt81_MuMu/D");
  outputTree->Branch("ST_rwgt82_MuMu", &ST_rwgt82_MuMu, "ST_rwgt82_MuMu/D");
  outputTree->Branch("ST_rwgt83_MuMu", &ST_rwgt83_MuMu, "ST_rwgt83_MuMu/D");
  outputTree->Branch("ST_rwgt84_MuMu", &ST_rwgt84_MuMu, "ST_rwgt84_MuMu/D");
  outputTree->Branch("ST_rwgt85_MuMu", &ST_rwgt85_MuMu, "ST_rwgt85_MuMu/D");
  outputTree->Branch("ST_rwgt86_MuMu", &ST_rwgt86_MuMu, "ST_rwgt86_MuMu/D");
  outputTree->Branch("ST_rwgt87_MuMu", &ST_rwgt87_MuMu, "ST_rwgt87_MuMu/D");
  outputTree->Branch("ST_bkg_MuMu", &ST_bkg_MuMu, "ST_bkg_MuMu/D");
  outputTree->Branch("ST_rwgt77_ElMu", &ST_rwgt77_ElMu, "ST_rwgt77_ElMu/D");
  outputTree->Branch("ST_rwgt78_ElMu", &ST_rwgt78_ElMu, "ST_rwgt78_ElMu/D");
  outputTree->Branch("ST_rwgt79_ElMu", &ST_rwgt79_ElMu, "ST_rwgt79_ElMu/D");
  outputTree->Branch("ST_rwgt80_ElMu", &ST_rwgt80_ElMu, "ST_rwgt80_ElMu/D");
  outputTree->Branch("ST_rwgt81_ElMu", &ST_rwgt81_ElMu, "ST_rwgt81_ElMu/D");
  outputTree->Branch("ST_rwgt82_ElMu", &ST_rwgt82_ElMu, "ST_rwgt82_ElMu/D");
  outputTree->Branch("ST_rwgt83_ElMu", &ST_rwgt83_ElMu, "ST_rwgt83_ElMu/D");
  outputTree->Branch("ST_rwgt84_ElMu", &ST_rwgt84_ElMu, "ST_rwgt84_ElMu/D");
  outputTree->Branch("ST_rwgt85_ElMu", &ST_rwgt85_ElMu, "ST_rwgt85_ElMu/D");
  outputTree->Branch("ST_rwgt86_ElMu", &ST_rwgt86_ElMu, "ST_rwgt86_ElMu/D");
  outputTree->Branch("ST_rwgt87_ElMu", &ST_rwgt87_ElMu, "ST_rwgt87_ElMu/D");
  outputTree->Branch("ST_bkg_ElMu", &ST_bkg_ElMu, "ST_bkg_ElMu/D");
  outputTree->Branch("ST_rwgt77_ElEl", &ST_rwgt77_ElEl, "ST_rwgt77_ElEl/D");
  outputTree->Branch("ST_rwgt78_ElEl", &ST_rwgt78_ElEl, "ST_rwgt78_ElEl/D");
  outputTree->Branch("ST_rwgt79_ElEl", &ST_rwgt79_ElEl, "ST_rwgt79_ElEl/D");
  outputTree->Branch("ST_rwgt80_ElEl", &ST_rwgt80_ElEl, "ST_rwgt80_ElEl/D");
  outputTree->Branch("ST_rwgt81_ElEl", &ST_rwgt81_ElEl, "ST_rwgt81_ElEl/D");
  outputTree->Branch("ST_rwgt82_ElEl", &ST_rwgt82_ElEl, "ST_rwgt82_ElEl/D");
  outputTree->Branch("ST_rwgt83_ElEl", &ST_rwgt83_ElEl, "ST_rwgt83_ElEl/D");
  outputTree->Branch("ST_rwgt84_ElEl", &ST_rwgt84_ElEl, "ST_rwgt84_ElEl/D");
  outputTree->Branch("ST_rwgt85_ElEl", &ST_rwgt85_ElEl, "ST_rwgt85_ElEl/D");
  outputTree->Branch("ST_rwgt86_ElEl", &ST_rwgt86_ElEl, "ST_rwgt86_ElEl/D");
  outputTree->Branch("ST_rwgt87_ElEl", &ST_rwgt87_ElEl, "ST_rwgt87_ElEl/D");
  outputTree->Branch("ST_bkg_ElEl", &ST_bkg_ElEl, "ST_bkg_ElEl/D");
  outputTree->Branch("weight_rwgt77_MuMu", &weight_rwgt77_MuMu, "weight_rwgt77_MuMu/D");
  outputTree->Branch("weight_rwgt78_MuMu", &weight_rwgt78_MuMu, "weight_rwgt78_MuMu/D");
  outputTree->Branch("weight_rwgt79_MuMu", &weight_rwgt79_MuMu, "weight_rwgt79_MuMu/D");
  outputTree->Branch("weight_rwgt80_MuMu", &weight_rwgt80_MuMu, "weight_rwgt80_MuMu/D");
  outputTree->Branch("weight_rwgt81_MuMu", &weight_rwgt81_MuMu, "weight_rwgt81_MuMu/D");
  outputTree->Branch("weight_rwgt82_MuMu", &weight_rwgt82_MuMu, "weight_rwgt82_MuMu/D");
  outputTree->Branch("weight_rwgt83_MuMu", &weight_rwgt83_MuMu, "weight_rwgt83_MuMu/D");
  outputTree->Branch("weight_rwgt84_MuMu", &weight_rwgt84_MuMu, "weight_rwgt84_MuMu/D");
  outputTree->Branch("weight_rwgt85_MuMu", &weight_rwgt85_MuMu, "weight_rwgt85_MuMu/D");
  outputTree->Branch("weight_rwgt86_MuMu", &weight_rwgt86_MuMu, "weight_rwgt86_MuMu/D");
  outputTree->Branch("weight_rwgt87_MuMu", &weight_rwgt87_MuMu, "weight_rwgt87_MuMu/D");
  outputTree->Branch("weight_bkg_MuMu", &weight_bkg_MuMu, "weight_bkg_MuMu/D");
  outputTree->Branch("weight_rwgt77_ElMu", &weight_rwgt77_ElMu, "weight_rwgt77_ElMu/D");
  outputTree->Branch("weight_rwgt78_ElMu", &weight_rwgt78_ElMu, "weight_rwgt78_ElMu/D");
  outputTree->Branch("weight_rwgt79_ElMu", &weight_rwgt79_ElMu, "weight_rwgt79_ElMu/D");
  outputTree->Branch("weight_rwgt80_ElMu", &weight_rwgt80_ElMu, "weight_rwgt80_ElMu/D");
  outputTree->Branch("weight_rwgt81_ElMu", &weight_rwgt81_ElMu, "weight_rwgt81_ElMu/D");
  outputTree->Branch("weight_rwgt82_ElMu", &weight_rwgt82_ElMu, "weight_rwgt82_ElMu/D");
  outputTree->Branch("weight_rwgt83_ElMu", &weight_rwgt83_ElMu, "weight_rwgt83_ElMu/D");
  outputTree->Branch("weight_rwgt84_ElMu", &weight_rwgt84_ElMu, "weight_rwgt84_ElMu/D");
  outputTree->Branch("weight_rwgt85_ElMu", &weight_rwgt85_ElMu, "weight_rwgt85_ElMu/D");
  outputTree->Branch("weight_rwgt86_ElMu", &weight_rwgt86_ElMu, "weight_rwgt86_ElMu/D");
  outputTree->Branch("weight_rwgt87_ElMu", &weight_rwgt87_ElMu, "weight_rwgt87_ElMu/D");
  outputTree->Branch("weight_bkg_ElMu", &weight_bkg_ElMu, "weight_bkg_ElMu/D");
  outputTree->Branch("weight_rwgt77_ElEl", &weight_rwgt77_ElEl, "weight_rwgt77_ElEl/D");
  outputTree->Branch("weight_rwgt78_ElEl", &weight_rwgt78_ElEl, "weight_rwgt78_ElEl/D");
  outputTree->Branch("weight_rwgt79_ElEl", &weight_rwgt79_ElEl, "weight_rwgt79_ElEl/D");
  outputTree->Branch("weight_rwgt80_ElEl", &weight_rwgt80_ElEl, "weight_rwgt80_ElEl/D");
  outputTree->Branch("weight_rwgt81_ElEl", &weight_rwgt81_ElEl, "weight_rwgt81_ElEl/D");
  outputTree->Branch("weight_rwgt82_ElEl", &weight_rwgt82_ElEl, "weight_rwgt82_ElEl/D");
  outputTree->Branch("weight_rwgt83_ElEl", &weight_rwgt83_ElEl, "weight_rwgt83_ElEl/D");
  outputTree->Branch("weight_rwgt84_ElEl", &weight_rwgt84_ElEl, "weight_rwgt84_ElEl/D");
  outputTree->Branch("weight_rwgt85_ElEl", &weight_rwgt85_ElEl, "weight_rwgt85_ElEl/D");
  outputTree->Branch("weight_rwgt86_ElEl", &weight_rwgt86_ElEl, "weight_rwgt86_ElEl/D");
  outputTree->Branch("weight_rwgt87_ElEl", &weight_rwgt87_ElEl, "weight_rwgt87_ElEl/D");
  outputTree->Branch("weight_bkg_ElEl", &weight_bkg_ElEl, "weight_bkg_ElEl/D");

  Int_t           run;
  Int_t           lumi;
  ULong64_t       evt;
  Int_t           isData;
  Float_t         evt_scale1fb;
  Int_t           evt_passgoodrunlist;
  Int_t           HLT_DoubleMu;
  Int_t           HLT_DoubleEl;
  Int_t           HLT_DoubleEl_DZ;
  Int_t           HLT_DoubleEl_DZ_2;
  Int_t           HLT_MuEG;
  Int_t           HLT_SingleIsoEl8;
  Int_t           HLT_SingleIsoEl17;
  Int_t           HLT_SingleIsoMu8;
  Int_t           HLT_SingleIsoMu17;
  Int_t           mc_HLT_DoubleMu;
  Int_t           mc_HLT_DoubleEl;
  Int_t           mc_HLT_DoubleEl_DZ;
  Int_t           mc_HLT_DoubleEl_DZ_2;
  Int_t           mc_HLT_MuEG;
  Int_t           mc_HLT_SingleIsoEl8;
  Int_t           mc_HLT_SingleIsoEl17;
  Int_t           mc_HLT_SingleIsoMu8;
  Int_t           mc_HLT_SingleIsoMu17;
  vector<float>   *lep_pt;
  vector<float>   *lep_eta;
  vector<float>   *lep_phi;
  vector<float>   *lep_coneCorrPt;
  vector<float>   *lep_ip3d;
  vector<float>   *lep_ip3derr;
  vector<int>     *lep_isTriggerSafe_v1;
  vector<int>     *lep_lostHits;
  vector<int>     *lep_convVeto;
  vector<int>     *lep_motherIdSS;
  vector<int>     *lep_pass_VVV_cutbased_3l_fo;
  vector<int>     *lep_pass_VVV_cutbased_3l_tight;
  vector<int>     *lep_pass_VVV_cutbased_fo;
  vector<int>     *lep_pass_VVV_cutbased_tight;
  vector<int>     *lep_pass_VVV_cutbased_veto;
  vector<int>     *lep_pass_VVV_cutbased_fo_noiso;
  vector<int>     *lep_pass_VVV_cutbased_tight_noiso;
  vector<int>     *lep_pass_VVV_cutbased_veto_noiso;
  vector<int>     *lep_pdgId;
  vector<float>   *lep_dxy;
  vector<float>   *lep_dz;
  vector<float>   *lep_ptRatio;
  vector<float>   *lep_ptRel;
  vector<float>   *lep_pterr;
  vector<float>   *lep_relIso03EAv2;
  vector<float>   *lep_relIso04EAv2;
  vector<int>     *lep_tightCharge;
  vector<float>   *lep_trk_pt;
  vector<int>     *lep_charge;
  vector<float>   *lep_etaSC;
  vector<float>   *lep_MVA;
  vector<int>     *lep_isFromW;
  vector<int>     *lep_isFromZ;
  vector<int>     *lep_isFromB;
  vector<int>     *lep_isFromC;
  vector<int>     *lep_isFromL;
  vector<int>     *lep_isFromLF;
  vector<int>     *lep_genPart_index;
  vector<float>   *jets_csv;
  vector<float>   *jets_up_csv;
  vector<float>   *jets_dn_csv;
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_up_pt;
  Float_t         met_up_phi;
  Float_t         met_dn_pt;
  Float_t         met_dn_phi;
  Int_t           firstgoodvertex;
  Int_t           nTrueInt;
  Int_t           nVert;
  Int_t           nisoTrack_mt2_cleaned_VVV_cutbased_veto;
  Float_t         weight_btagsf;
  Float_t         weight_btagsf_heavy_DN;
  Float_t         weight_btagsf_heavy_UP;
  Float_t         weight_btagsf_light_DN;
  Float_t         weight_btagsf_light_UP;
  Float_t         gen_ht;
  vector<int>     *genPart_motherId;
  vector<int>     *genPart_pdgId;
  vector<int>     *genPart_charge;
  vector<int>     *genPart_status;
  Int_t           ngenLep;
  Int_t           ngenLepFromTau;
  Int_t           Flag_AllEventFilters;
  Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
  Int_t           Flag_HBHEIsoNoiseFilter;
  Int_t           Flag_HBHENoiseFilter;
  Int_t           Flag_badChargedCandidateFilter;
  Int_t           Flag_badMuonFilter;
  Int_t           Flag_badMuonFilterv2;
  Int_t           Flag_badChargedCandidateFilterv2;
  Int_t           Flag_eeBadScFilter;
  Int_t           Flag_globalTightHalo2016;
  Int_t           Flag_goodVertices;
  Int_t           Flag_ecalLaserCorrFilter;
  Int_t           Flag_hcalLaserEventFilter;
  Int_t           Flag_trackingFailureFilter;
  Int_t           Flag_CSCTightHaloFilter;
  Int_t           Flag_CSCTightHalo2015Filter;
  Int_t           Flag_badMuons;
  Int_t           Flag_duplicateMuons;
  Int_t           Flag_noBadMuons;
  Int_t           nVlep;
  Int_t           nTlep;
  Int_t           nTlepSS;
  Int_t           nLlep;
  Int_t           nSFOS;
  Int_t           nSFOSinZ;
  Int_t           nj;
  Int_t           nj_up;
  Int_t           nj_dn;
  Int_t           nj30;
  Int_t           nj30_up;
  Int_t           nj30_dn;
  Int_t           nb;
  Int_t           nb_up;
  Int_t           nb_dn;
  Float_t         Mjj;
  Float_t         Mjj_up;
  Float_t         Mjj_dn;
  Float_t         MjjVBF;
  Float_t         MjjVBF_up;
  Float_t         MjjVBF_dn;
  Float_t         DetajjVBF;
  Float_t         DetajjVBF_up;
  Float_t         DetajjVBF_dn;
  Float_t         MjjL;
  Float_t         MjjL_up;
  Float_t         MjjL_dn;
  Float_t         DetajjL;
  Float_t         DetajjL_up;
  Float_t         DetajjL_dn;
  Float_t         MllSS;
  Float_t         MeeSS;
  Float_t         Mll3L;
  Float_t         Mee3L;
  Float_t         Mll3L1;
  Float_t         M3l;
  Float_t         Pt3l;
  Float_t         M01;
  Float_t         M02;
  Float_t         M12;
  Int_t           isSFOS01;
  Int_t           isSFOS02;
  Int_t           isSFOS12;
  Float_t         DPhi3lMET;
  Float_t         DPhi3lMET_up;
  Float_t         DPhi3lMET_dn;
  Float_t         MTmax;
  Float_t         MTmax_up;
  Float_t         MTmax_dn;
  Float_t         MTmin;
  Float_t         MTmin_up;
  Float_t         MTmin_dn;
  Float_t         MT3rd;
  Float_t         MT3rd_up;
  Float_t         MT3rd_dn;
  Float_t         MTmax3L;
  Float_t         MTmax3L_up;
  Float_t         MTmax3L_dn;
  Int_t           passSSee;
  Int_t           passSSem;
  Int_t           passSSmm;
  Int_t           lep_idx0_SS;
  Int_t           lep_idx1_SS;
  TString         *bkgtype;
  Int_t           vetophoton;
  Float_t         purewgt;
  Float_t         purewgt_up;
  Float_t         purewgt_dn;
  Float_t         ffwgt;
  Float_t         ffwgt_up;
  Float_t         ffwgt_dn;
  Float_t         ffwgt_el_up;
  Float_t         ffwgt_el_dn;
  Float_t         ffwgt_mu_up;
  Float_t         ffwgt_mu_dn;
  Float_t         ffwgt_closure_up;
  Float_t         ffwgt_closure_dn;
  Float_t         ffwgt_closure_el_up;
  Float_t         ffwgt_closure_el_dn;
  Float_t         ffwgt_closure_mu_up;
  Float_t         ffwgt_closure_mu_dn;
  Float_t         ffwgt_full_up;
  Float_t         ffwgt_full_dn;
  Float_t         ffwgtqcd;
  Float_t         ffwgtqcd_up;
  Float_t         ffwgtqcd_dn;
  Float_t         lepsf;
  Float_t         lepsf_up;
  Float_t         lepsf_dn;
  Float_t         trigeff;
  Float_t         trigeff_up;
  Float_t         trigeff_dn;
  Float_t         trigsf;
  Float_t         trigsf_up;
  Float_t         trigsf_dn;
  std::vector<LorentzVector>  *jets_p4;
  std::vector<LorentzVector>  *ak8jets_p4;
  vector<float>   *ak8jets_softdropMass;
  vector<float>   *ak8jets_prunedMass;
  vector<float>   *ak8jets_mass;
  vector<float>   *ak8jets_nJettinessTau1;
  vector<float>   *ak8jets_nJettinessTau2;
  vector<float>   *ak8jets_softdropPuppiSubjet1;
  vector<float>   *ak8jets_softdropPuppiSubjet2;
  vector<float>   *ak8jets_puppi_softdropMass;
  vector<float>   *ak8jets_puppi_nJettinessTau1;
  vector<float>   *ak8jets_puppi_nJettinessTau2;
  vector<float>   *ak8jets_puppi_eta;
  vector<float>   *ak8jets_puppi_phi;
  vector<float>   *ak8jets_puppi_pt;
  vector<float>   *ak8jets_puppi_mass;
  //aqgc variables
  vector<float>   *genweights;
  vector<TString> *genweightsID;
  Float_t         genps_origxwgtup;

  lep_pt = 0;
  lep_eta = 0; 
  lep_pt = 0;
  lep_eta = 0;
  lep_phi = 0;
  lep_coneCorrPt = 0;
  lep_ip3d = 0;
  lep_ip3derr = 0;
  lep_isTriggerSafe_v1 = 0;
  lep_lostHits = 0;
  lep_convVeto = 0;
  lep_motherIdSS = 0;
  lep_pass_VVV_cutbased_3l_fo = 0;
  lep_pass_VVV_cutbased_3l_tight = 0;
  lep_pass_VVV_cutbased_fo = 0;
  lep_pass_VVV_cutbased_tight = 0;
  lep_pass_VVV_cutbased_veto = 0;
  lep_pass_VVV_cutbased_fo_noiso = 0;
  lep_pass_VVV_cutbased_tight_noiso = 0;
  lep_pass_VVV_cutbased_veto_noiso = 0;
  lep_pdgId = 0;
  lep_dxy = 0;
  lep_dz = 0;
  lep_ptRatio = 0;
  lep_ptRel = 0;
  lep_pterr = 0;
  lep_relIso03EAv2 = 0;
  lep_relIso04EAv2 = 0;
  lep_tightCharge = 0;
  lep_trk_pt = 0;
  lep_charge = 0;
  lep_etaSC = 0;
  lep_MVA = 0;
  lep_isFromW = 0;
  lep_isFromZ = 0;
  lep_isFromB = 0;
  lep_isFromC = 0;
  lep_isFromL = 0;
  lep_isFromLF = 0;
  lep_genPart_index = 0;
  jets_csv = 0;
  jets_up_csv = 0;
  jets_dn_csv = 0; 
  genPart_motherId = 0;
  genPart_pdgId = 0;
  genPart_charge = 0;
  genPart_status = 0;
  jets_p4 = 0;
  ak8jets_p4 = 0;
  ak8jets_softdropMass = 0;
  ak8jets_prunedMass = 0;
  ak8jets_mass = 0;
  ak8jets_nJettinessTau1 = 0;
  ak8jets_nJettinessTau2 = 0;
  ak8jets_softdropPuppiSubjet1 = 0;
  ak8jets_softdropPuppiSubjet2 = 0;
  ak8jets_puppi_softdropMass = 0;
  ak8jets_puppi_nJettinessTau1 = 0;
  ak8jets_puppi_nJettinessTau2 = 0;
  ak8jets_puppi_eta = 0;
  ak8jets_puppi_phi = 0;
  ak8jets_puppi_pt = 0;
  ak8jets_puppi_mass = 0;
  genweights = 0;
  genweightsID = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("evt", &(evt));
  tree->SetBranchAddress("isData", &(isData));
  tree->SetBranchAddress("evt_scale1fb", &(evt_scale1fb));
  tree->SetBranchAddress("evt_passgoodrunlist", &(evt_passgoodrunlist));
  tree->SetBranchAddress("firstgoodvertex", &(firstgoodvertex));
  tree->SetBranchAddress("HLT_DoubleMu", &(HLT_DoubleMu));
  tree->SetBranchAddress("HLT_DoubleEl", &(HLT_DoubleEl));
  tree->SetBranchAddress("HLT_DoubleEl_DZ", &(HLT_DoubleEl_DZ));
  tree->SetBranchAddress("HLT_DoubleEl_DZ_2", &(HLT_DoubleEl_DZ_2));
  tree->SetBranchAddress("HLT_MuEG", &(HLT_MuEG));
  tree->SetBranchAddress("mc_HLT_DoubleMu", &(mc_HLT_DoubleMu));
  tree->SetBranchAddress("mc_HLT_DoubleEl", &(mc_HLT_DoubleEl));
  tree->SetBranchAddress("mc_HLT_DoubleEl_DZ", &(mc_HLT_DoubleEl_DZ));
  tree->SetBranchAddress("mc_HLT_DoubleEl_DZ_2", &(mc_HLT_DoubleEl_DZ_2));
  tree->SetBranchAddress("mc_HLT_MuEG", &(mc_HLT_MuEG));
  tree->SetBranchAddress("lep_pt", &(lep_pt));
  tree->SetBranchAddress("lep_eta", &(lep_eta));
  tree->SetBranchAddress("lep_phi", &(lep_phi));
  tree->SetBranchAddress("lep_coneCorrPt", &(lep_coneCorrPt));
  tree->SetBranchAddress("lep_ip3d", &(lep_ip3d));
  tree->SetBranchAddress("lep_ip3derr", &(lep_ip3derr));
  tree->SetBranchAddress("lep_isTriggerSafe_v1", &(lep_isTriggerSafe_v1));
  tree->SetBranchAddress("lep_lostHits", &(lep_lostHits));
  tree->SetBranchAddress("lep_convVeto", &(lep_convVeto));
  tree->SetBranchAddress("lep_motherIdSS", &(lep_motherIdSS));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_3l_fo", &(lep_pass_VVV_cutbased_3l_fo));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_3l_tight", &(lep_pass_VVV_cutbased_3l_tight));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_fo", &(lep_pass_VVV_cutbased_fo));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_tight", &(lep_pass_VVV_cutbased_tight));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_veto", &(lep_pass_VVV_cutbased_veto));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_fo_noiso", &(lep_pass_VVV_cutbased_fo_noiso));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_tight_noiso", &(lep_pass_VVV_cutbased_tight_noiso));
  tree->SetBranchAddress("lep_pass_VVV_cutbased_veto_noiso", &(lep_pass_VVV_cutbased_veto_noiso)); 
  tree->SetBranchAddress("lep_pdgId", &(lep_pdgId));
  tree->SetBranchAddress("lep_dxy", &(lep_dxy));
  tree->SetBranchAddress("lep_dz", &(lep_dz));
  tree->SetBranchAddress("lep_ptRatio", &(lep_ptRatio));
  tree->SetBranchAddress("lep_ptRel", &(lep_ptRel));
  tree->SetBranchAddress("lep_pterr", &(lep_pterr));
  tree->SetBranchAddress("lep_relIso03EAv2", &(lep_relIso03EAv2));
  tree->SetBranchAddress("lep_relIso04EAv2", &(lep_relIso04EAv2));
  tree->SetBranchAddress("lep_tightCharge", &(lep_tightCharge));
  tree->SetBranchAddress("lep_trk_pt", &(lep_trk_pt));
  tree->SetBranchAddress("lep_charge", &(lep_charge));
  tree->SetBranchAddress("lep_etaSC", &(lep_etaSC));
  tree->SetBranchAddress("lep_MVA", &(lep_MVA));
  tree->SetBranchAddress("lep_isFromW", &(lep_isFromW));
  tree->SetBranchAddress("lep_isFromZ", &(lep_isFromZ));
  tree->SetBranchAddress("lep_isFromB", &(lep_isFromB));
  tree->SetBranchAddress("lep_isFromC", &(lep_isFromC));
  tree->SetBranchAddress("lep_isFromL", &(lep_isFromL));
  tree->SetBranchAddress("lep_isFromLF", &(lep_isFromLF)); 
  tree->SetBranchAddress("passSSee", &(passSSee));
  tree->SetBranchAddress("passSSem", &(passSSem));
  tree->SetBranchAddress("passSSmm", &(passSSmm));
  tree->SetBranchAddress("met_pt", &(met_pt));
  tree->SetBranchAddress("met_phi", &(met_phi));
  tree->SetBranchAddress("met_up_pt", &(met_up_pt));
  tree->SetBranchAddress("met_up_phi", &(met_up_phi));
  tree->SetBranchAddress("met_dn_pt", &(met_dn_pt));
  tree->SetBranchAddress("met_dn_phi", &(met_dn_phi));
  tree->SetBranchAddress("purewgt", &(purewgt));
  tree->SetBranchAddress("purewgt_up", &(purewgt_up));
  tree->SetBranchAddress("purewgt_dn", &(purewgt_dn));
  tree->SetBranchAddress("weight_btagsf", &(weight_btagsf));
  tree->SetBranchAddress("weight_btagsf_heavy_DN", &(weight_btagsf_heavy_DN));
  tree->SetBranchAddress("weight_btagsf_heavy_UP", &(weight_btagsf_heavy_UP));
  tree->SetBranchAddress("weight_btagsf_light_DN", &(weight_btagsf_light_DN));
  tree->SetBranchAddress("weight_btagsf_light_UP", &(weight_btagsf_light_UP));
  tree->SetBranchAddress("lepsf", &(lepsf));
  tree->SetBranchAddress("lepsf_up", &(lepsf_up));
  tree->SetBranchAddress("lepsf_dn", &(lepsf_dn));
  tree->SetBranchAddress("trigeff", &(trigeff));
  tree->SetBranchAddress("trigeff_up", &(trigeff_up));
  tree->SetBranchAddress("trigeff_dn", &(trigeff_dn));
  tree->SetBranchAddress("trigsf", &(trigsf));
  tree->SetBranchAddress("trigsf_up", &(trigsf_up));
  tree->SetBranchAddress("trigsf_dn", &(trigsf_dn));
  tree->SetBranchAddress("nVlep", &(nVlep));
  tree->SetBranchAddress("nTlep", &(nTlep));
  tree->SetBranchAddress("nTlepSS", &(nTlepSS));
  tree->SetBranchAddress("nLlep", &(nLlep));
  tree->SetBranchAddress("nj30", &(nj30));
  tree->SetBranchAddress("nj30_up", &(nj30_up));
  tree->SetBranchAddress("nj30_dn", &(nj30_dn));
  tree->SetBranchAddress("nb", &(nb));
  tree->SetBranchAddress("nb_up", &(nb_up));
  tree->SetBranchAddress("nb_dn", &(nb_dn));
  tree->SetBranchAddress("Mjj", &(Mjj));
  tree->SetBranchAddress("Mjj_up", &(Mjj_up));
  tree->SetBranchAddress("Mjj_dn", &(Mjj_dn));
  tree->SetBranchAddress("DetajjL", &(DetajjL));
  tree->SetBranchAddress("DetajjL_up", &(DetajjL_up));
  tree->SetBranchAddress("DetajjL_dn", &(DetajjL_dn));
  tree->SetBranchAddress("MjjL", &(MjjL));
  tree->SetBranchAddress("MjjL_up", &(MjjL_up));
  tree->SetBranchAddress("MjjL_dn", &(MjjL_dn));
  tree->SetBranchAddress("MllSS", &(MllSS));
  tree->SetBranchAddress("vetophoton", &(vetophoton));
  tree->SetBranchAddress("MTmax", &(MTmax));
  tree->SetBranchAddress("MTmax_up", &(MTmax_up));
  tree->SetBranchAddress("MTmax_dn", &(MTmax_dn));
  tree->SetBranchAddress("nisoTrack_mt2_cleaned_VVV_cutbased_veto", &(nisoTrack_mt2_cleaned_VVV_cutbased_veto));
  tree->SetBranchAddress("jets_p4", &(jets_p4));
  tree->SetBranchAddress("ak8jets_p4", &(ak8jets_p4));   
  tree->SetBranchAddress("ak8jets_softdropMass", &(ak8jets_softdropMass));
  tree->SetBranchAddress("ak8jets_prunedMass", &(ak8jets_prunedMass));
  tree->SetBranchAddress("ak8jets_mass", &(ak8jets_mass));
  tree->SetBranchAddress("ak8jets_nJettinessTau1", &(ak8jets_nJettinessTau1));
  tree->SetBranchAddress("ak8jets_nJettinessTau2", &(ak8jets_nJettinessTau2));
  tree->SetBranchAddress("ak8jets_softdropPuppiSubjet1", &(ak8jets_softdropPuppiSubjet1));
  tree->SetBranchAddress("ak8jets_softdropPuppiSubjet2", &(ak8jets_softdropPuppiSubjet2));
  tree->SetBranchAddress("ak8jets_puppi_softdropMass", &(ak8jets_puppi_softdropMass));
  tree->SetBranchAddress("ak8jets_puppi_nJettinessTau1", &(ak8jets_puppi_nJettinessTau1));
  tree->SetBranchAddress("ak8jets_puppi_nJettinessTau2", &(ak8jets_puppi_nJettinessTau2));
  tree->SetBranchAddress("ak8jets_puppi_eta", &(ak8jets_puppi_eta));
  tree->SetBranchAddress("ak8jets_puppi_phi", &(ak8jets_puppi_phi));
  tree->SetBranchAddress("ak8jets_puppi_pt", &(ak8jets_puppi_pt));
  tree->SetBranchAddress("ak8jets_puppi_mass", &(ak8jets_puppi_mass));
  if(sample=="aQGC") tree->SetBranchAddress("genweights", &(genweights));
  if(sample=="aQGC") tree->SetBranchAddress("genweightsID", &(genweightsID));
  if(sample=="aQGC") tree->SetBranchAddress("genps_origxwgtup", &(genps_origxwgtup));

  HistCollection mumuHistCut1;
  initializeHistCollection(mumuHistCut1, "2SSTL_MuMu");
  HistCollection elelHistCut1;
  initializeHistCollection(elelHistCut1, "2SSTL_ElEl");
  HistCollection elmuHistCut1;
  initializeHistCollection(elmuHistCut1, "2SSTL_ElMu");

  TH1D *h_TotalEvents_MuMu_ST_SS = new TH1D("h_TotalEvents_MuMu_ST_SS", "h_TotalEvents_MuMu_ST_SS", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS = new TH1D("h_TotalEvents_ElMu_ST_SS", "h_TotalEvents_ElMu_ST_SS", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS = new TH1D("h_TotalEvents_ElEl_ST_SS", "h_TotalEvents_ElEl_ST_SS", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS->Sumw2();

  //ft0 histograms
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt77 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt77", "h_TotalEvents_MuMu_ST_SS_rwgt77", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt77->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt78 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt78", "h_TotalEvents_MuMu_ST_SS_rwgt78", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt78->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt79 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt79", "h_TotalEvents_MuMu_ST_SS_rwgt79", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt79->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt80 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt80", "h_TotalEvents_MuMu_ST_SS_rwgt80", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt80->Sumw2();  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt81 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt81", "h_TotalEvents_MuMu_ST_SS_rwgt81", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt81->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt82 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt82", "h_TotalEvents_MuMu_ST_SS_rwgt82", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt82->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt83 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt83", "h_TotalEvents_MuMu_ST_SS_rwgt83", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt83->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt84 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt84", "h_TotalEvents_MuMu_ST_SS_rwgt84", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt84->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt85 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt85", "h_TotalEvents_MuMu_ST_SS_rwgt85", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt85->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt86 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt86", "h_TotalEvents_MuMu_ST_SS_rwgt86", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt86->Sumw2();
  TH1D *h_TotalEvents_MuMu_ST_SS_rwgt87 = new TH1D("h_TotalEvents_MuMu_ST_SS_rwgt87", "h_TotalEvents_MuMu_ST_SS_rwgt87", 15, -0.5, 14.5); h_TotalEvents_MuMu_ST_SS_rwgt87->Sumw2();

  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt77 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt77", "h_TotalEvents_ElMu_ST_SS_rwgt77", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt77->Sumw2();  
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt78 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt78", "h_TotalEvents_ElMu_ST_SS_rwgt78", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt78->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt79 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt79", "h_TotalEvents_ElMu_ST_SS_rwgt79", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt79->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt80 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt80", "h_TotalEvents_ElMu_ST_SS_rwgt80", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt80->Sumw2();  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt81 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt81", "h_TotalEvents_ElMu_ST_SS_rwgt81", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt81->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt82 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt82", "h_TotalEvents_ElMu_ST_SS_rwgt82", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt82->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt83 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt83", "h_TotalEvents_ElMu_ST_SS_rwgt83", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt83->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt84 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt84", "h_TotalEvents_ElMu_ST_SS_rwgt84", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt84->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt85 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt85", "h_TotalEvents_ElMu_ST_SS_rwgt85", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt85->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt86 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt86", "h_TotalEvents_ElMu_ST_SS_rwgt86", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt86->Sumw2();
  TH1D *h_TotalEvents_ElMu_ST_SS_rwgt87 = new TH1D("h_TotalEvents_ElMu_ST_SS_rwgt87", "h_TotalEvents_ElMu_ST_SS_rwgt87", 15, -0.5, 14.5); h_TotalEvents_ElMu_ST_SS_rwgt87->Sumw2();

  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt77 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt77", "h_TotalEvents_ElEl_ST_SS_rwgt77", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt77->Sumw2();  
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt78 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt78", "h_TotalEvents_ElEl_ST_SS_rwgt78", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt78->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt79 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt79", "h_TotalEvents_ElEl_ST_SS_rwgt79", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt79->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt80 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt80", "h_TotalEvents_ElEl_ST_SS_rwgt80", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt80->Sumw2();  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt81 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt81", "h_TotalEvents_ElEl_ST_SS_rwgt81", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt81->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt82 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt82", "h_TotalEvents_ElEl_ST_SS_rwgt82", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt82->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt83 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt83", "h_TotalEvents_ElEl_ST_SS_rwgt83", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt83->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt84 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt84", "h_TotalEvents_ElEl_ST_SS_rwgt84", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt84->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt85 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt85", "h_TotalEvents_ElEl_ST_SS_rwgt85", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt85->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt86 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt86", "h_TotalEvents_ElEl_ST_SS_rwgt86", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt86->Sumw2();
  TH1D *h_TotalEvents_ElEl_ST_SS_rwgt87 = new TH1D("h_TotalEvents_ElEl_ST_SS_rwgt87", "h_TotalEvents_ElEl_ST_SS_rwgt87", 15, -0.5, 14.5); h_TotalEvents_ElEl_ST_SS_rwgt87->Sumw2();

  TH1D *h_ST_rwgt77 = new TH1D("h_ST_rwgt77", "h_ST_rwgt77", 100, 0.0, 5000.0); h_ST_rwgt77->Sumw2();
  TH1D *h_ST_rwgt78 = new TH1D("h_ST_rwgt78", "h_ST_rwgt78", 100, 0.0, 5000.0); h_ST_rwgt78->Sumw2();
  TH1D *h_ST_rwgt79 = new TH1D("h_ST_rwgt79", "h_ST_rwgt79", 100, 0.0, 5000.0); h_ST_rwgt79->Sumw2();
  TH1D *h_ST_rwgt80 = new TH1D("h_ST_rwgt80", "h_ST_rwgt80", 100, 0.0, 5000.0); h_ST_rwgt80->Sumw2();
  TH1D *h_ST_rwgt81 = new TH1D("h_ST_rwgt81", "h_ST_rwgt81", 100, 0.0, 5000.0); h_ST_rwgt81->Sumw2();
  TH1D *h_ST_rwgt82 = new TH1D("h_ST_rwgt82", "h_ST_rwgt82", 100, 0.0, 5000.0); h_ST_rwgt82->Sumw2();
  TH1D *h_ST_rwgt83 = new TH1D("h_ST_rwgt83", "h_ST_rwgt83", 100, 0.0, 5000.0); h_ST_rwgt83->Sumw2();
  TH1D *h_ST_rwgt84 = new TH1D("h_ST_rwgt84", "h_ST_rwgt84", 100, 0.0, 5000.0); h_ST_rwgt84->Sumw2();
  TH1D *h_ST_rwgt85 = new TH1D("h_ST_rwgt85", "h_ST_rwgt85", 100, 0.0, 5000.0); h_ST_rwgt85->Sumw2();
  TH1D *h_ST_rwgt86 = new TH1D("h_ST_rwgt86", "h_ST_rwgt86", 100, 0.0, 5000.0); h_ST_rwgt86->Sumw2();
  TH1D *h_ST_rwgt87 = new TH1D("h_ST_rwgt87", "h_ST_rwgt87", 100, 0.0, 5000.0); h_ST_rwgt87->Sumw2();
  TH1D *h_ST_ElEl = new TH1D("h_ST_ElEl", "h_ST_ElEl", 100, 0.0, 5000.0); h_ST_ElEl->Sumw2();
  TH1D *h_ST_ElMu = new TH1D("h_ST_ElMu", "h_ST_ElMu", 100, 0.0, 5000.0); h_ST_ElMu->Sumw2();
  TH1D *h_ST_MuMu = new TH1D("h_ST_MuMu", "h_ST_MuMu", 100, 0.0, 5000.0); h_ST_MuMu->Sumw2();

  int nEvents=tree->GetEntries();
  bool passgenfilterList = false;
  double dilepton_noST = 0.0;
  double dilepton_noST_rawevents = 0.0;
  double dilepton_ST1000 = 0.0;
  double dilepton_ST1000_rawevents = 0.0;
  double dilepton_ST1500 = 0.0;
  double dilepton_ST1500_rawevents = 0.0;
  double dilepton_ST2000 = 0.0;
  double dilepton_ST2000_rawevents = 0.0;

  double dilepton_noST_SM = 0.0;
  double dilepton_noST_rawevents_SM = 0.0;
  double dilepton_ST1000_SM = 0.0;
  double dilepton_ST1000_rawevents_SM = 0.0;
  double dilepton_ST1500_SM = 0.0;
  double dilepton_ST1500_rawevents_SM = 0.0;
  double dilepton_ST2000_SM = 0.0;
  double dilepton_ST2000_rawevents_SM = 0.0;  

  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i); 

    ST_rwgt77_MuMu=ST_rwgt78_MuMu=ST_rwgt79_MuMu=ST_rwgt80_MuMu=ST_rwgt81_MuMu=ST_rwgt82_MuMu=ST_rwgt83_MuMu=ST_rwgt84_MuMu=ST_rwgt85_MuMu=ST_rwgt86_MuMu=ST_rwgt87_MuMu=ST_bkg_MuMu=0.0;
    ST_rwgt77_ElMu=ST_rwgt78_ElMu=ST_rwgt79_ElMu=ST_rwgt80_ElMu=ST_rwgt81_ElMu=ST_rwgt82_ElMu=ST_rwgt83_ElMu=ST_rwgt84_ElMu=ST_rwgt85_ElMu=ST_rwgt86_ElMu=ST_rwgt87_ElMu=ST_bkg_ElMu=0.0;
    ST_rwgt77_ElEl=ST_rwgt78_ElEl=ST_rwgt79_ElEl=ST_rwgt80_ElEl=ST_rwgt81_ElEl=ST_rwgt82_ElEl=ST_rwgt83_ElEl=ST_rwgt84_ElEl=ST_rwgt85_ElEl=ST_rwgt86_ElEl=ST_rwgt87_ElEl=ST_bkg_ElEl=0.0;
    weight_rwgt77_MuMu=weight_rwgt78_MuMu=weight_rwgt79_MuMu=weight_rwgt80_MuMu=weight_rwgt81_MuMu=weight_rwgt82_MuMu=weight_rwgt83_MuMu=weight_rwgt84_MuMu=weight_rwgt85_MuMu=weight_rwgt86_MuMu=weight_rwgt87_MuMu=weight_bkg_MuMu=0.0;
    weight_rwgt77_ElMu=weight_rwgt78_ElMu=weight_rwgt79_ElMu=weight_rwgt80_ElMu=weight_rwgt81_ElMu=weight_rwgt82_ElMu=weight_rwgt83_ElMu=weight_rwgt84_ElMu=weight_rwgt85_ElMu=weight_rwgt86_ElMu=weight_rwgt87_ElMu=weight_bkg_ElMu=0.0;
    weight_rwgt77_ElEl=weight_rwgt78_ElEl=weight_rwgt79_ElEl=weight_rwgt80_ElEl=weight_rwgt81_ElEl=weight_rwgt82_ElEl=weight_rwgt83_ElEl=weight_rwgt84_ElEl=weight_rwgt85_ElEl=weight_rwgt86_ElEl=weight_rwgt87_ElEl=weight_bkg_ElEl=0.0;
 
    double weight = evt_scale1fb*purewgt;
    if(evt_passgoodrunlist==0) continue;
    if(firstgoodvertex!=0) continue;
    if(vetophoton!=0) continue; 

    if(type=="Data" and Trigger=="HLT_MuEG" and HLT_MuEG==0) continue;
    if(type=="Data" and Trigger=="HLT_DoubleEl" and HLT_DoubleEl==0) continue;
    if(type=="Data" and Trigger=="HLT_DoubleMu" and HLT_DoubleMu==0) continue;

    std::vector<TLorentzVector> v_selectedJets;
    for(unsigned int iselJet=0; iselJet<jets_p4->size(); ++iselJet)
    {
      TLorentzVector Jet;

      if(fabs(jets_p4->at(iselJet).Eta())<2.5 and jets_p4->at(iselJet).Pt()>30.0)
      {
        Jet.SetPtEtaPhiE(jets_p4->at(iselJet).Pt(), jets_p4->at(iselJet).Eta(), jets_p4->at(iselJet).Phi(), jets_p4->at(iselJet).E());
        v_selectedJets.push_back(Jet);
      }
    }

    if(nj30!=(int)v_selectedJets.size()) std::cout << "size mismatch found" << std::endl;

    std::vector<fatJetInfo> v_fatJets;
    for (unsigned int ifatjet=0; ifatjet<ak8jets_puppi_pt->size(); ifatjet++)
    {
      fatJetInfo fatjet;
      fatjet.ak8JetPrunmass = ak8jets_prunedMass->at(ifatjet);
      fatjet.ak8Jetsd0 = ak8jets_puppi_softdropMass->at(ifatjet);
      fatjet.ak8JetPt = ak8jets_puppi_pt->at(ifatjet);
      fatjet.ak8JetEta = ak8jets_puppi_eta->at(ifatjet);
      fatjet.ak8JetPhi = ak8jets_puppi_phi->at(ifatjet);
      fatjet.ak8JetTau1 = ak8jets_puppi_nJettinessTau1->at(ifatjet);
      fatjet.ak8JetTau2 = ak8jets_puppi_nJettinessTau2->at(ifatjet);
      v_fatJets.push_back(fatjet);
    }
   
    std::sort (v_fatJets.begin(), v_fatJets.end(), sortFatJetVectorsInDescendingpT);

    std::vector<genInfo> v_genweights;
    float originalWeight = 0.0;
    if(sample=="aQGC")  originalWeight = genps_origxwgtup;  
 
    if(sample=="aQGC")
    {
      for (unsigned int igen=0; igen<genweights->size(); igen++)
      {
        genInfo genW;
        genW.genWeight = genweights->at(igen);
        genW.genWeightID = genweightsID->at(igen);
        v_genweights.push_back(genW);
      }
    }
    //SS 
    if((lep_pdgId->at(0)*lep_pdgId->at(1)==169) and (nVlep==2)*(nLlep==2)*(nTlep==2) and passSSmm==1)//mm
    { 
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, ptrel0, ptrel1, deltaR1, deltaR2;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi =  ptrel0 = ptrel1 = deltaR1 = deltaR2 = 0.0;
      TLorentzVector mu1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), MUON_MASS);
      TLorentzVector mu2 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), MUON_MASS); 
      weight *= lepsf*trigsf; 
      if(type=="MC" and mc_HLT_DoubleMu==1)
      {
        if(nisoTrack_mt2_cleaned_VVV_cutbased_veto==0) 
        {
          if(v_selectedJets.size()>=2)
          //if(nj30>=2)
          {
            if(nb==0)
            { 
              double ST = mu1.Pt() + mu2.Pt() + met_pt;
              for(unsigned int ir=0; ir<v_selectedJets.size(); ir++) ST += v_selectedJets.at(ir).Pt();
              fillMuHistCollection(mumuHistCut1, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), ptrel0, ptrel1, lep_relIso04EAv2->at(0), lep_relIso04EAv2->at(1), deltaR1, deltaR2, (mu1+mu2).M(), met_pt, nj30, nb, MjjL, DetajjL, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, ST, weight);
              if(sample=="aQGC")
              {
                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                {
                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                  /*if(v_genweights.at(ig).genWeightID=="rwgt_77") h_ST_rwgt77->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight);
                  */
                }
              }
              //couts to compute background
              if(DetajjL < 1.5 and MllSS > 40 and met_pt > 60)
              {
                if(sample=="aQGC") 
                {
                  for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                  {
                    double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(1, genWeight); 
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(1, genWeight); 
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_noST += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_noST_rawevents++;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") ST_rwgt87_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") ST_rwgt78_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") ST_rwgt79_MuMu=ST; 
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") ST_rwgt80_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") ST_rwgt81_MuMu=ST; 
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") ST_rwgt82_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") ST_rwgt83_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") ST_rwgt84_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") ST_rwgt85_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") ST_rwgt86_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") ST_rwgt87_MuMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") weight_rwgt87_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") weight_rwgt78_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") weight_rwgt79_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") weight_rwgt80_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") weight_rwgt81_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") weight_rwgt82_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") weight_rwgt83_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") weight_rwgt84_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") weight_rwgt85_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") weight_rwgt86_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") weight_rwgt87_MuMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_SM += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_rawevents_SM++;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") h_ST_rwgt87->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight);
                  }
                }
                ST_bkg_MuMu=ST;
                weight_bkg_MuMu=weight; 
                h_ST_MuMu->Fill(ST, weight);
                h_TotalEvents_MuMu_ST_SS->Fill(1, weight);
                if(ST > 250)
                {
                  if(sample=="aQGC") 
                  {
                    for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                    {
                      double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                      if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(2, genWeight);
                    }
                  }
                  h_TotalEvents_MuMu_ST_SS->Fill(2, weight);
                  if(ST > 500)
                  {
                    if(sample=="aQGC")
                    {
                      for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                      {
                        double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                        if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(3, genWeight);
                      }
                    }
                    h_TotalEvents_MuMu_ST_SS->Fill(3, weight); 
                    if(ST > 750)
                    { 
                      if(sample=="aQGC")
                      {
                        for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                        {
                          double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                          if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(4, genWeight);
                        }
                      } 
                      h_TotalEvents_MuMu_ST_SS->Fill(4, weight);
                      if(ST > 1000)
                      {
                        if(sample=="aQGC")
                        {
                          for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                          {
                            double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000 += genWeight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000_rawevents++;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_SM += genWeight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_rawevents_SM++;
                          }
                        } 
                        h_TotalEvents_MuMu_ST_SS->Fill(5, weight);
                        if(ST > 1500)
                        {
                          if(sample=="aQGC")
                          {
                            for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                            {
                              double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500 += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500_rawevents++;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_SM += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_rawevents_SM++;
                            }
                          }
                          h_TotalEvents_MuMu_ST_SS->Fill(6, weight);
                          if(ST > 2000)
                          {
                            if(sample=="aQGC")
                            {
                              for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                              {
                                double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000 += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000_rawevents++;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_SM += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_rawevents_SM++;
                              }
                            }
                            h_TotalEvents_MuMu_ST_SS->Fill(7, weight);
                            if(ST > 2500)
                            {
                              if(sample=="aQGC")
                              {
                                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                                {
                                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                  if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_MuMu_ST_SS_rwgt77->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_MuMu_ST_SS_rwgt78->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_MuMu_ST_SS_rwgt79->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_MuMu_ST_SS_rwgt80->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_MuMu_ST_SS_rwgt81->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_MuMu_ST_SS_rwgt82->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_MuMu_ST_SS_rwgt83->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_MuMu_ST_SS_rwgt84->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_MuMu_ST_SS_rwgt85->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_MuMu_ST_SS_rwgt86->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_MuMu_ST_SS_rwgt87->Fill(8, genWeight);
                                }
                              }
                              h_TotalEvents_MuMu_ST_SS->Fill(8, weight);
                            }//2500
                          }//2000
                        }//1500 
                      }//1000 
                    }//750
                  }//500
                }//250  
              }//no cut                
            }//deltaeta MllSS met cuts
          }//b-veto
        }//nj>=2
      }//isolated track cut
    }//trigger 
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==143) and (nVlep==2)*(nLlep==2)*(nTlep==2) and passSSem==1)//em
    { 
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector mu1, el1;
      double muip3d, mudxy, mudz, eip3d, edxy, edz, mupTRatio, elpTRatio;
      muip3d=mudxy=mudz=eip3d=edxy=edz=mupTRatio=elpTRatio=0.0; 
      if(abs(lep_pdgId->at(0))==13)
      {      
        mu1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), MUON_MASS); 
        mupTRatio = lep_ptRatio->at(0);
        muip3d = fabs(lep_ip3d->at(0));
        mudxy = fabs(lep_dxy->at(0));
        mudz = fabs(lep_dz->at(0));  
      }
      else if(abs(lep_pdgId->at(0))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
        elpTRatio = lep_ptRatio->at(0);
        eip3d = fabs(lep_ip3d->at(0));
        edxy = fabs(lep_dxy->at(0));
        edz = fabs(lep_dz->at(0));  
      }
      if(abs(lep_pdgId->at(1))==13)
      {
        mu1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), MUON_MASS);
        mupTRatio = lep_ptRatio->at(1);
        muip3d = fabs(lep_ip3d->at(1));
        mudxy = fabs(lep_dxy->at(1));
        mudz = fabs(lep_dz->at(1));  
      }
      else if(abs(lep_pdgId->at(1))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
        elpTRatio = lep_ptRatio->at(1);
        eip3d = fabs(lep_ip3d->at(1));
        edxy = fabs(lep_dxy->at(1));
        edz = fabs(lep_dz->at(1));
      }
      weight *= lepsf*trigsf;
      if(type=="MC" and mc_HLT_MuEG==1)
      {
        if(nisoTrack_mt2_cleaned_VVV_cutbased_veto==0)
        {
          if(v_selectedJets.size()>=2)
          //if(nj30>=2) 
          {
            if(nb==0) 
            {
              double ST = el1.Pt() + mu1.Pt() + met_pt;
              for(unsigned int ir=0; ir<v_selectedJets.size(); ir++) ST += v_selectedJets.at(ir).Pt();
              fillElMuHistCollection(elmuHistCut1, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), mupTRatio, elpTRatio, muip3d, eip3d, mudxy, edxy, mudz, edz, (el1+mu1).M(), met_pt, nj30, nb, MjjL, DetajjL, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, ST, weight);
              if(sample=="aQGC")
              {
                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                {
                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                  /*if(v_genweights.at(ig).genWeightID=="rwgt_77") h_ST_rwgt87->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight);*/
                }
              }
              if(DetajjL < 1.5 and MllSS > 30 and met_pt > 60 and MTmax > 90.0)
              {
                if(sample=="aQGC")
                {
                  for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                  {
                    double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_noST += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_noST_rawevents++;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") ST_rwgt87_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") ST_rwgt78_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") ST_rwgt79_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") ST_rwgt80_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") ST_rwgt81_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") ST_rwgt82_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") ST_rwgt83_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") ST_rwgt84_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") ST_rwgt85_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") ST_rwgt86_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") ST_rwgt87_ElMu=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") weight_rwgt87_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") weight_rwgt78_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") weight_rwgt79_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") weight_rwgt80_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") weight_rwgt81_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") weight_rwgt82_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") weight_rwgt83_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") weight_rwgt84_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") weight_rwgt85_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") weight_rwgt86_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") weight_rwgt87_ElMu=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_SM += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_rawevents_SM++;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight);
                  }
                }
                ST_bkg_ElMu=ST;
                weight_bkg_ElMu=weight;
                h_ST_ElMu->Fill(ST, weight);
                h_TotalEvents_ElMu_ST_SS->Fill(1, weight);
                if(ST > 250)
                {
                  if(sample=="aQGC")
                  {
                    for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                    {
                      double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                      if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(2, genWeight);
                    }
                  }
                  h_TotalEvents_ElMu_ST_SS->Fill(2, weight);
                  if(ST > 500)
                  {
                    if(sample=="aQGC")
                    {
                      for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                      {
                        double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                        if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(3, genWeight);
                      }
                    }
                    h_TotalEvents_ElMu_ST_SS->Fill(3, weight);
                    if(ST > 750)
                    {
                      if(sample=="aQGC")
                      {
                        for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                        {
                          double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                          if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(4, genWeight);
                        }
                      }
                      h_TotalEvents_ElMu_ST_SS->Fill(4, weight);
                      if(ST > 1000)
                      {
                        if(sample=="aQGC")
                        {
                          for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                          {
                            double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000 += genWeight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000_rawevents++;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_SM += genWeight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_rawevents_SM++;
                          }
                        }
                        h_TotalEvents_ElMu_ST_SS->Fill(5, weight);
                        if(ST > 1500)
                        {
                          if(sample=="aQGC")
                          {
                            for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                            {
                              double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500 += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500_rawevents++;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_SM += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_rawevents_SM++;
                            }
                          }
                          h_TotalEvents_ElMu_ST_SS->Fill(6, weight);
                          if(ST > 2000)
                          {
                            if(sample=="aQGC")
                            {
                              for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                              {
                                double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000 += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000_rawevents++;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_SM += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_rawevents_SM++;
                              }
                            }
                            h_TotalEvents_ElMu_ST_SS->Fill(7, weight);
                            if(ST > 2500)
                            {
                              if(sample=="aQGC")
                              {
                                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                                {
                                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                  if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElMu_ST_SS_rwgt77->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElMu_ST_SS_rwgt78->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElMu_ST_SS_rwgt79->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElMu_ST_SS_rwgt80->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElMu_ST_SS_rwgt81->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElMu_ST_SS_rwgt82->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElMu_ST_SS_rwgt83->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElMu_ST_SS_rwgt84->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElMu_ST_SS_rwgt85->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElMu_ST_SS_rwgt86->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElMu_ST_SS_rwgt87->Fill(8, genWeight);
                                }
                              }
                              h_TotalEvents_ElMu_ST_SS->Fill(8, weight);
                            }//2500
                          }//2000
                        }//1500 
                      }//1000 
                    }//750
                  }//500
                }//250  
              }//dletaetajj and mllss and met and mtmax cuts
            }//nb==0  
          }//nb==0 
        }//njets >= 2
      }//isolated tracks
    }//trigger
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==121) and (nVlep==2)*(nLlep==2)*(nTlep==2) and passSSee==1)//ee
    {
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
      TLorentzVector el2 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
      weight *= lepsf*trigsf;
      if(type=="MC" and mc_HLT_DoubleEl_DZ_2==1)
      {
        if(nisoTrack_mt2_cleaned_VVV_cutbased_veto==0)
        {
          if(v_selectedJets.size()>=2)
          //if(nj30>=2)
          {
            if(nb==0) 
            {
              double ST = el1.Pt() + el2.Pt() + met_pt;
              for(unsigned int ir=0; ir<v_selectedJets.size(); ir++) ST += v_selectedJets.at(ir).Pt();
              fillElHistCollection(elelHistCut1, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), (el1+el2).M(), met_pt, nj30, nb, MjjL, DetajjL, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, ST, weight);
              if(sample=="aQGC")
              {
                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                {
                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                  /*if(v_genweights.at(ig).genWeightID=="rwgt_77") h_ST_rwgt87->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight);*/
                }
              }
              if(DetajjL < 1.5 and MllSS > 40 and met_pt > 60 and abs(MllSS-91.1876)>10.0)
              {
                if(sample=="aQGC")
                {
                  for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                  {
                    double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(1, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") ST_rwgt87_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") ST_rwgt78_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") ST_rwgt79_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") ST_rwgt80_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") ST_rwgt81_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") ST_rwgt82_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") ST_rwgt83_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") ST_rwgt84_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") ST_rwgt85_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") ST_rwgt86_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") ST_rwgt87_ElEl=ST;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") weight_rwgt87_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") weight_rwgt78_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") weight_rwgt79_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") weight_rwgt80_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") weight_rwgt81_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") weight_rwgt82_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") weight_rwgt83_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") weight_rwgt84_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") weight_rwgt85_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") weight_rwgt86_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") weight_rwgt87_ElEl=genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_noST += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_SM += genWeight;
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_noST_rawevents_SM++;
                    if(v_genweights.at(ig).genWeightID=="rwgt_78") h_ST_rwgt78->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_79") h_ST_rwgt79->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_80") h_ST_rwgt80->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_81") h_ST_rwgt81->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_82") h_ST_rwgt82->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_83") h_ST_rwgt83->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_84") h_ST_rwgt84->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_85") h_ST_rwgt85->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_86") h_ST_rwgt86->Fill(ST, genWeight);
                    if(v_genweights.at(ig).genWeightID=="rwgt_87") h_ST_rwgt87->Fill(ST, genWeight); 
                  }
                }
                ST_bkg_ElEl=ST;
                weight_bkg_ElEl=weight;
                h_ST_ElEl->Fill(ST, weight);
                h_TotalEvents_ElEl_ST_SS->Fill(1, weight);
                if(ST > 250)
                {
                  if(sample=="aQGC")
                  {
                    for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                    {
                      double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                      if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(2, genWeight);
                      if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(2, genWeight);
                    }
                  } 
                  h_TotalEvents_ElEl_ST_SS->Fill(2, weight);
                  if(ST > 500)
                  {
                    if(sample=="aQGC")
                    {
                      for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                      {
                        double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                        if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(3, genWeight);
                        if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(3, genWeight);
                      }
                    }
                    h_TotalEvents_ElEl_ST_SS->Fill(3, weight);
                    if(ST > 750)
                    {
                      if(sample=="aQGC")
                      {
                        for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                        {
                          double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                          if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(4, genWeight);
                          if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(4, genWeight);
                        }
                      } 
                      h_TotalEvents_ElEl_ST_SS->Fill(4, weight);
                      if(ST > 1000)
                      { 
                        if(sample=="aQGC")
                        {
                          for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                          {
                            double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(5, genWeight);
                            if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000 += genWeight;
                 	    if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1000_rawevents++;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_SM += genWeight;
                            if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1000_rawevents_SM++;
                          }
                        } 
                        h_TotalEvents_ElEl_ST_SS->Fill(5, weight);
                        if(ST > 1500)
                        {
                          if(sample=="aQGC")
                          {
                            for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                            {
                              double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(6, genWeight);
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500 += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST1500_rawevents++;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_SM += genWeight;
                              if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST1500_rawevents_SM++;
                            }
                          } 
                          h_TotalEvents_ElEl_ST_SS->Fill(6, weight);
                          if(ST > 2000)
                          {
                            if(sample=="aQGC")
                            {
                              for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                              {
                                double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(7, genWeight);
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000 += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_77") dilepton_ST2000_rawevents++;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_SM += genWeight;
                                if(v_genweights.at(ig).genWeightID=="rwgt_82") dilepton_ST2000_rawevents_SM++;
                              }
                            } 
                            h_TotalEvents_ElEl_ST_SS->Fill(7, weight);
                            if(ST > 2500)
                            {
                              if(sample=="aQGC")
                              {
                                for(unsigned int ig=0; ig<v_genweights.size(); ig++)
                                {
                                  double genWeight = ((float)v_genweights.at(ig).genWeight/(float)originalWeight)*weight;
                                  if(v_genweights.at(ig).genWeightID=="rwgt_77") h_TotalEvents_ElEl_ST_SS_rwgt77->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_78") h_TotalEvents_ElEl_ST_SS_rwgt78->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_79") h_TotalEvents_ElEl_ST_SS_rwgt79->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_80") h_TotalEvents_ElEl_ST_SS_rwgt80->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_81") h_TotalEvents_ElEl_ST_SS_rwgt81->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_82") h_TotalEvents_ElEl_ST_SS_rwgt82->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_83") h_TotalEvents_ElEl_ST_SS_rwgt83->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_84") h_TotalEvents_ElEl_ST_SS_rwgt84->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_85") h_TotalEvents_ElEl_ST_SS_rwgt85->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_86") h_TotalEvents_ElEl_ST_SS_rwgt86->Fill(8, genWeight);
                                  if(v_genweights.at(ig).genWeightID=="rwgt_87") h_TotalEvents_ElEl_ST_SS_rwgt87->Fill(8, genWeight);
                                }
                              } 
                              h_TotalEvents_ElEl_ST_SS->Fill(8, weight);
                            }//2500
                          }//2000
                        }//1500
                      }//1000
                    }//750
                  }//500 
                }//2500   
              }//no cut  
            }//b-veto
          }//njet>=2 
        }//track veto
      }//trigger
    }//lep selection
    outputTree->Fill();
  }//event loop
  outputTree->Write();
  std::cout << "dilepton_noST = " << dilepton_noST << std::endl;
  std::cout << "dilepton_noST_rawevents = " << dilepton_noST_rawevents << std::endl;
  std::cout << "dilepton_ST1000 = " << dilepton_ST1000 << std::endl;
  std::cout << "dilepton_ST1000_rawevents = " << dilepton_ST1000_rawevents << std::endl;
  std::cout << "dilepton_ST1500 = " << dilepton_ST1500 << std::endl;
  std::cout << "dilepton_ST1500_rawevents = " << dilepton_ST1500_rawevents << std::endl;
  std::cout << "dilepton_ST2000 = " << dilepton_ST2000 << std::endl;
  std::cout << "dilepton_ST2000_rawevents = " << dilepton_ST2000_rawevents << std::endl;

  std::cout << "dilepton_noST_SM = " << dilepton_noST_SM << std::endl;
  std::cout << "dilepton_noST_SM_rawevents = " << dilepton_noST_rawevents_SM << std::endl;
  std::cout << "dilepton_ST1000_SM = " << dilepton_ST1000_SM << std::endl;
  std::cout << "dilepton_ST1000_SM_rawevents = " << dilepton_ST1000_rawevents_SM << std::endl;
  std::cout << "dilepton_ST1500_SM = " << dilepton_ST1500_SM << std::endl;
  std::cout << "dilepton_ST1500_SM_rawevents = " << dilepton_ST1500_rawevents_SM << std::endl;
  std::cout << "dilepton_ST2000_SM = " << dilepton_ST2000_SM << std::endl;
  std::cout << "dilepton_ST2000_SM_rawevents = " << dilepton_ST2000_rawevents_SM << std::endl;  

  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  tFile->mkdir("2SSTL"); 
  tFile->cd("2SSTL");
  writeHistCollection(mumuHistCut1);
  writeHistCollection(elelHistCut1);
  writeHistCollection(elmuHistCut1);
  tFile->cd();
  h_TotalEvents_MuMu_ST_SS->Write();
  h_TotalEvents_ElMu_ST_SS->Write();
  h_TotalEvents_ElEl_ST_SS->Write();
  h_ST_rwgt87->Write();
  h_ST_rwgt78->Write();
  h_ST_rwgt79->Write();
  h_ST_rwgt80->Write();
  h_ST_rwgt81->Write();
  h_ST_rwgt82->Write();
  h_ST_rwgt83->Write();
  h_ST_rwgt84->Write();
  h_ST_rwgt85->Write();
  h_ST_rwgt86->Write();
  h_ST_rwgt87->Write();
  h_ST_ElEl->Write();
  h_ST_ElMu->Write();
  h_ST_MuMu->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt87->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt78->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt79->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt80->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt81->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt82->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt83->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt84->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt85->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt86->Write();
  h_TotalEvents_MuMu_ST_SS_rwgt87->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt87->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt78->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt79->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt80->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt81->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt82->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt83->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt84->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt85->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt86->Write();
  h_TotalEvents_ElMu_ST_SS_rwgt87->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt87->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt78->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt79->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt80->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt81->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt82->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt83->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt84->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt85->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt86->Write();
  h_TotalEvents_ElEl_ST_SS_rwgt87->Write();
  tFile->Close();
  inputFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;
}

