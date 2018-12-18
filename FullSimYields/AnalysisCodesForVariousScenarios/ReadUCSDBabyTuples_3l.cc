#include "ReadUCSDBabyTuples_3l.h"
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

std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName);

bool triggerRequirement(std::vector<int> *lep_pdgID, int mc_El, int mc_El_DZ, int mc_ME, int mc_Mu)
{
  int lepprod01 = lep_pdgID->at(0)*lep_pdgID->at(1);
  if (abs(lepprod01) == 121 && (mc_El==1 || mc_El_DZ==1))
    return true;
  else if (abs(lepprod01) == 143 && mc_ME==1)
    return true;
  else if (abs(lepprod01) == 169 && mc_Mu==1)
    return true;

  int lepprod02 = lep_pdgID->at(0)*lep_pdgID->at(2);
  if (abs(lepprod02) == 121 && (mc_El==1 || mc_El_DZ==1))
     return true;
  else if (abs(lepprod02) == 143 && mc_ME==1)
     return true;
  else if (abs(lepprod02) == 169 && mc_Mu==1)
     return true;

  int lepprod12 = lep_pdgID->at(1)*lep_pdgID->at(2);
  if (abs(lepprod12) == 121 && (mc_El==1 || mc_El_DZ==1))
     return true;
  else if (abs(lepprod12) == 143 && mc_ME==1)
     return true;
  else if (abs(lepprod12) == 169 && mc_Mu==1)
     return true;

  return false;

}


int ReadUCSDBabyTuples_3l(std::string infile, std::string treeStr, std::string Sample, std::string trig, std::string Trigger="None")
{
  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  bool debug=false;
  //debug = true;
  std::string trigfilename=(trig+".root").c_str();
  TFile *trigFile = new TFile((trigfilename).c_str());

  TFile *sfTrkMuFile = new TFile("muon_trk_sf.root");
  TFile *sfidMuFile = new TFile("muon_id_sf.root");
  TFile *sfMuFile = new TFile("muon_sf.root");

  TFile *sfElIDFile = new TFile("egammaEffi.txt_EGM2D.root");
  TH2F  *idSFElHist = (TH2F*) sfElIDFile->Get("EGamma_SF2D");
  TFile *sfElMVAFile = new TFile("egammaEffi.txt_EGM2D.MVA80.root");
  TH2F  *idSFElMVAHist = (TH2F*) sfElMVAFile->Get("EGamma_SF2D");
  TFile *sfElFile = new TFile("elec_sf_iso.root");
  TH2F  *sfElHist = (TH2F*) sfElFile->Get("sf_pt_vs_eta");

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
  tree->SetBranchAddress("Flag_AllEventFilters", &(Flag_AllEventFilters));
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
  tree->SetBranchAddress("DPhi3lMET", &(DPhi3lMET));
  tree->SetBranchAddress("DPhi3lMET_up", &(DPhi3lMET_up));
  tree->SetBranchAddress("DPhi3lMET_dn", &(DPhi3lMET_dn));
  tree->SetBranchAddress("Mll3L", &(Mll3L));
  tree->SetBranchAddress("M3l", &(M3l));
  tree->SetBranchAddress("Mee3L", &(Mee3L)); 
  tree->SetBranchAddress("MTmax3L", &(MTmax3L));
  tree->SetBranchAddress("nSFOS", &(nSFOS)); 
  tree->SetBranchAddress("nj", &(nj));
  tree->SetBranchAddress("Pt3l", &(Pt3l));
  tree->SetBranchAddress("MjjL", &(MjjL));
  tree->SetBranchAddress("Mll3L1", &(Mll3L1));
  tree->SetBranchAddress("MjjL_up", &(MjjL_up));
  tree->SetBranchAddress("MjjL_dn", &(MjjL_dn));
  tree->SetBranchAddress("MllSS", &(MllSS));
  tree->SetBranchAddress("MT3rd", &(MT3rd));
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
  tree->SetBranchAddress("nSFOSinZ", &(nSFOSinZ));

  TH1D *h_ST_0SFOS = new TH1D("h_ST_0SFOS", "h_ST_0SFOS", 5000.0, 0.0, 10000.0); h_ST_0SFOS->Sumw2();
  TH1D *h_ST_1SFOS = new TH1D("h_ST_1SFOS", "h_ST_1SFOS", 5000.0, 0.0, 10000.0); h_ST_1SFOS->Sumw2();
  TH1D *h_ST_2SFOS = new TH1D("h_ST_2SFOS", "h_ST_2SFOS", 5000.0, 0.0, 10000.0); h_ST_2SFOS->Sumw2();

  TH1D *h_RawEvents_0SFOS = new TH1D("h_RawEvents_0SFOS", "h_RawEvents_0SFOS", 15, -0.5, 14.5); h_RawEvents_0SFOS->Sumw2();
  TH1D *h_RawEvents_1SFOS = new TH1D("h_RawEvents_1SFOS", "h_RawEvents_1SFOS", 15, -0.5, 14.5); h_RawEvents_1SFOS->Sumw2();
  TH1D *h_RawEvents_2SFOS = new TH1D("h_RawEvents_2SFOS", "h_RawEvents_2SFOS", 15, -0.5, 14.5); h_RawEvents_2SFOS->Sumw2();

  TH1D *h_TotalEvents_0SFOS = new TH1D("h_TotalEvents_0SFOS", "h_TotalEvents_0SFOS", 15, -0.5, 14.5); h_TotalEvents_0SFOS->Sumw2();
  TH1D *h_TotalEvents_1SFOS = new TH1D("h_TotalEvents_1SFOS", "h_TotalEvents_1SFOS", 15, -0.5, 14.5); h_TotalEvents_1SFOS->Sumw2();
  TH1D *h_TotalEvents_2SFOS = new TH1D("h_TotalEvents_2SFOS", "h_TotalEvents_2SFOS", 15, -0.5, 14.5); h_TotalEvents_2SFOS->Sumw2();

  TH1D *h_TotalEvents_0SFOS_ST = new TH1D("h_TotalEvents_0SFOS_ST", "h_TotalEvents_0SFOS_ST", 15, -0.5, 14.5); h_TotalEvents_0SFOS_ST->Sumw2();
  TH1D *h_TotalEvents_1SFOS_ST = new TH1D("h_TotalEvents_1SFOS_ST", "h_TotalEvents_1SFOS_ST", 15, -0.5, 14.5); h_TotalEvents_1SFOS_ST->Sumw2();
  TH1D *h_TotalEvents_2SFOS_ST = new TH1D("h_TotalEvents_2SFOS_ST", "h_TotalEvents_2SFOS_ST", 15, -0.5, 14.5); h_TotalEvents_2SFOS_ST->Sumw2(); 

  int nEvents=tree->GetEntries();
  bool passgenfilterList = false;
  int n_events_nSFOS = 0.0;
  double n_events_nSFOS_weighted = 0.0;
  int n_events_nSFOS1 = 0.0;
  double n_events_nSFOS1_weighted = 0.0;
  int n_events_nSFOS2 = 0.0;
  double n_events_nSFOS2_weighted = 0.0;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);  
    double weight = evt_scale1fb*purewgt;
    if(evt_passgoodrunlist==0) continue;
    if(firstgoodvertex!=0) continue;
    if(vetophoton!=0) continue; 
    if(Flag_AllEventFilters!=1) continue;
    if(Sample=="Data" and Trigger=="HLT_MuEG" and HLT_MuEG==0) continue;
    if(Sample=="Data" and Trigger=="HLT_DoubleEl" and HLT_DoubleEl==0) continue;
    if(Sample=="Data" and Trigger=="HLT_DoubleMu" and HLT_DoubleMu==0) continue;

    std::vector<TLorentzVector> v_selectedJets;
    for(unsigned int iselJet=0; iselJet<jets_p4->size(); ++iselJet)
    {
      TLorentzVector Jet;

      if(fabs(jets_p4->at(iselJet).Eta())<5.0 and jets_p4->at(iselJet).Pt()>30.0)
      {
        Jet.SetPtEtaPhiE(jets_p4->at(iselJet).Pt(), jets_p4->at(iselJet).Eta(), jets_p4->at(iselJet).Phi(), jets_p4->at(iselJet).E());
        v_selectedJets.push_back(Jet);
      }
    }

    if(nj!=(int)v_selectedJets.size()) std::cout << "size mismatch found" << std::endl;

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

    //3l
    if((nVlep==3)*(nLlep==3)*(nTlep==3))
    {
      bool trig = triggerRequirement(lep_pdgId, mc_HLT_DoubleEl, mc_HLT_DoubleEl_DZ, mc_HLT_MuEG, mc_HLT_DoubleMu);
      if(not trig) continue;
      weight *= lepsf*trigsf;
      //0SFOS
      double ST = 0.0;
      if(lep_pt->at(0) > 25. and lep_pt->at(1) > 20. and lep_pt->at(2) > 20. and nSFOS==0 and nj<=1 and nb==0)
      {
        ST = lep_pt->at(0) + lep_pt->at(1) + lep_pt->at(2);
        ST += met_pt;
        for(unsigned int i=0; i<v_selectedJets.size(); i++) ST += v_selectedJets.at(i).Pt();
        h_ST_0SFOS->Fill(ST, weight);
        h_RawEvents_0SFOS->Fill(1);
        h_TotalEvents_0SFOS->Fill(1, weight);
        if(nj<=1 and nb==0 and DPhi3lMET>2.5 and met_pt>30.0 and Mll3L > 20. and abs(M3l-91.1876) > 15. and abs(Mee3L-91.1876) > 15. and MTmax3L > 90.0)
        {
          n_events_nSFOS++;
          n_events_nSFOS_weighted += weight;
          h_RawEvents_0SFOS->Fill(2);
          h_TotalEvents_0SFOS->Fill(2, weight);
          h_TotalEvents_0SFOS_ST->Fill(1, weight);
          if(ST > 1000)
          {
            h_TotalEvents_0SFOS_ST->Fill(2, weight);
            if(ST > 1500)
            {
              h_TotalEvents_0SFOS_ST->Fill(3, weight);
              if(ST > 2000)
              {
                h_TotalEvents_0SFOS_ST->Fill(4, weight);
                if(ST > 2500)
                {
                  h_TotalEvents_0SFOS_ST->Fill(5, weight);
                }
              }
            }
          }
        } 
      } 
      else if(lep_pt->at(0) > 25. and nSFOS==1 and nj<=1 and nb==0)
      {
        ST = lep_pt->at(0) + lep_pt->at(1) + lep_pt->at(2);
        ST += met_pt;
        for(unsigned int i=0; i<v_selectedJets.size(); i++) ST += v_selectedJets.at(i).Pt();
        h_ST_1SFOS->Fill(ST, weight);
        h_RawEvents_1SFOS->Fill(1);
        h_TotalEvents_1SFOS->Fill(1, weight);
        if(Pt3l>60.0 and DPhi3lMET>2.5 and met_pt>40.0 and Mll3L > 20. and abs(M3l-91.1876) > 10. and nSFOSinZ==0 and MT3rd > 90.0)
        {
          n_events_nSFOS1++;
          n_events_nSFOS1_weighted += weight;
          h_RawEvents_1SFOS->Fill(2);
          h_TotalEvents_1SFOS->Fill(2, weight);
          h_TotalEvents_1SFOS_ST->Fill(1, weight);
          if(ST > 1000)
          {
            h_TotalEvents_1SFOS_ST->Fill(2, weight);
            if(ST > 1500)
            {
              h_TotalEvents_1SFOS_ST->Fill(3, weight);
              if(ST > 2000)
              {
                h_TotalEvents_1SFOS_ST->Fill(4, weight);
                if(ST > 2500)
                {
                  h_TotalEvents_1SFOS_ST->Fill(5, weight);
                }
              }
            }
          } 
        }
      } 
      else if(lep_pt->at(0) > 25. and nSFOS==2 and nj<=1 and nb==0)
      {
        ST = lep_pt->at(0) + lep_pt->at(1) + lep_pt->at(2);
        ST += met_pt;
        for(unsigned int i=0; i<v_selectedJets.size(); i++) ST += v_selectedJets.at(i).Pt();
        h_ST_2SFOS->Fill(ST, weight);
        h_RawEvents_2SFOS->Fill(1);
        h_TotalEvents_2SFOS->Fill(1, weight);
        if(Pt3l>60.0 and DPhi3lMET>2.5 and met_pt>55.0 and (Mll3L > 20. && Mll3L1 > 20.)and abs(M3l-91.1876) > 10. and nSFOSinZ == 0)
        {
          n_events_nSFOS2++;
          n_events_nSFOS2_weighted += weight;
          h_RawEvents_2SFOS->Fill(2);
          h_TotalEvents_2SFOS->Fill(2, weight);
          h_TotalEvents_2SFOS_ST->Fill(1, weight);
          if(ST > 1000)
          {
            h_TotalEvents_2SFOS_ST->Fill(2, weight);
            if(ST > 1500)
            {
              h_TotalEvents_2SFOS_ST->Fill(3, weight);
              if(ST > 2000)
              {
                h_TotalEvents_2SFOS_ST->Fill(4, weight);
                if(ST > 2500)
                {
                  h_TotalEvents_2SFOS_ST->Fill(5, weight);
                }
              }
            }
          }
        }
      } 
    }
    if(debug)
    {
      std::cout << "lep_pt->at(0) = " << lep_pt->at(0) << std::endl;
      std::cout << "nSFOS = " << nSFOS << std::endl;
      std::cout << "nj = " << nj << std::endl;
      std::cout << "nb = " << nb << std::endl;
      std::cout << "Pt3l = " << Pt3l << std::endl;
      std::cout << "DPhi3lMET = " << DPhi3lMET << std::endl;
      std::cout << "met_pt = " << met_pt  << std::endl;
      std::cout << "Mll3L = " << Mll3L << std::endl;
      std::cout << "Mll3L1 = " << Mll3L1 << std::endl; 
      std::cout << "M3l = " << M3l << std::endl;
      std::cout << "nSFOSinZ = " << nSFOSinZ << std::endl;
    }
  }//event loop

  std::cout << "n_events_nSFOS = " << n_events_nSFOS << std::endl;
  std::cout << "n_events_nSFOS_weighted = " << n_events_nSFOS_weighted << std::endl;
  std::cout << "n_events_nSFOS1 = " << n_events_nSFOS1 << std::endl;
  std::cout << "n_events_nSFOS1_weighted = " << n_events_nSFOS1_weighted << std::endl;
  std::cout << "n_events_nSFOS2 = " << n_events_nSFOS2 << std::endl;
  std::cout << "n_events_nSFOS2_weighted = " << n_events_nSFOS2_weighted << std::endl;
  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  h_ST_0SFOS->Write();
  h_ST_1SFOS->Write();
  h_ST_2SFOS->Write();
  h_TotalEvents_0SFOS->Write();
  h_TotalEvents_1SFOS->Write();
  h_TotalEvents_2SFOS->Write();
  h_RawEvents_0SFOS->Write();
  h_RawEvents_1SFOS->Write();
  h_RawEvents_2SFOS->Write();
  h_TotalEvents_0SFOS_ST->Write();
  h_TotalEvents_1SFOS_ST->Write();
  h_TotalEvents_2SFOS_ST->Write();
  tFile->Close();
  inputFile->Close();
  trigFile->Close();
  sfTrkMuFile->Close();
  sfidMuFile->Close();
  sfMuFile->Close();
  sfElIDFile->Close();
  sfElMVAFile->Close();
  sfElFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;
}

// function to make an event list object for MET filtering
std::map<unsigned, std::set<unsigned> > readEventList(char const* _fileName) {
  std::map<unsigned, std::set<unsigned> > list;
  ifstream listFile(_fileName);
  std::cout << _fileName << std::endl;
  if (!listFile.is_open())
    throw std::runtime_error(_fileName);

  unsigned iL(0);
  std::string line;
  while (true) {
    std::getline(listFile, line);
    if (!listFile.good())
      break;
  

    if (line.find(":") == std::string::npos || line.find(":") == line.rfind(":"))
      continue;

    unsigned run(std::atoi(line.substr(0, line.find(":")).c_str()));
    unsigned event(std::atoi(line.substr(line.rfind(":") + 1).c_str()));

    //std::cout << "run = " << run << std::endl;
    //std::cout << "event = " << event << std::endl;

    list[run].insert(event);

    ++iL;
  }

  std::cout << "Loaded " << iL << " events" << std::endl;

  return list;
}
