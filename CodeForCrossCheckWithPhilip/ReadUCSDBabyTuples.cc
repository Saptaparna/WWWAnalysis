#include "ReadUCSDBabyTuples.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

using std::string;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

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

int ReadUCSDBabyTuples(std::string infile, std::string treeStr, std::string Sample, std::string trig, std::string Trigger="None")
{

  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  bool debug=false;

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
  //TLorentzVector *jets_p4; 

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

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("evt", &(evt));
  tree->SetBranchAddress("isData", &(isData));
  tree->SetBranchAddress("evt_scale1fb", &(evt_scale1fb));
  tree->SetBranchAddress("evt_passgoodrunlist", &(evt_passgoodrunlist));
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
   
  HistCollection mumuHistCut1;
  initializeHistCollection(mumuHistCut1, "2OSTL_MuMu");
  HistCollection elelHistCut1;
  initializeHistCollection(elelHistCut1, "2OSTL_ElEl");
  HistCollection elmuHistCut1;
  initializeHistCollection(elmuHistCut1, "2OSTL_ElMu");

  HistCollection mumuHistCut2;
  initializeHistCollection(mumuHistCut2, "2SSTL_MuMu");
  HistCollection elelHistCut2;
  initializeHistCollection(elelHistCut2, "2SSTL_ElEl");
  HistCollection elmuHistCut2;
  initializeHistCollection(elmuHistCut2, "2SSTL_ElMu");

  HistCollection mumuHistCut3;
  initializeHistCollection(mumuHistCut3, "2SS1JetCR_MuMu");
  HistCollection elelHistCut3;
  initializeHistCollection(elelHistCut3, "2SS1JetCR_ElEl");
  HistCollection elmuHistCut3;
  initializeHistCollection(elmuHistCut3, "2SS1JetCR_ElMu");

  TH1D *h_TotalEvents_MuMu_OS = new TH1D("h_TotalEvents_MuMu_OS", "h_TotalEvents_MuMu_OS", 15, -0.5, 14.5); h_TotalEvents_MuMu_OS->Sumw2();
  TH1D *h_TotalEvents_MuMu_SS = new TH1D("h_TotalEvents_MuMu_SS", "h_TotalEvents_MuMu_SS", 15, -0.5, 14.5); h_TotalEvents_MuMu_SS->Sumw2();
  TH1D *h_TotalEvents_ElMu_OS = new TH1D("h_TotalEvents_ElMu_OS", "h_TotalEvents_ElMu_OS", 15, -0.5, 14.5); h_TotalEvents_ElMu_OS->Sumw2();
  TH1D *h_TotalEvents_ElMu_SS = new TH1D("h_TotalEvents_ElMu_SS", "h_TotalEvents_ElMu_SS", 15, -0.5, 14.5); h_TotalEvents_ElMu_SS->Sumw2();
  TH1D *h_TotalEvents_ElEl_OS = new TH1D("h_TotalEvents_ElEl_OS", "h_TotalEvents_ElEl_OS", 15, -0.5, 14.5); h_TotalEvents_ElEl_OS->Sumw2();
  TH1D *h_TotalEvents_ElEl_SS = new TH1D("h_TotalEvents_ElEl_SS", "h_TotalEvents_ElEl_SS", 15, -0.5, 14.5); h_TotalEvents_ElEl_SS->Sumw2();

  TH1D *h_RawEvents_MuMu_SS = new TH1D("h_RawEvents_MuMu_SS", "h_RawEvents_MuMu_SS", 15, -0.5, 14.5); h_RawEvents_MuMu_SS->Sumw2();
  TH1D *h_RawEvents_ElMu_SS = new TH1D("h_RawEvents_ElMu_SS", "h_RawEvents_ElMu_SS", 15, -0.5, 14.5); h_RawEvents_ElMu_SS->Sumw2();
  TH1D *h_RawEvents_ElEl_SS = new TH1D("h_RawEvents_ElEl_SS", "h_RawEvents_ElEl_SS", 15, -0.5, 14.5); h_RawEvents_ElEl_SS->Sumw2();

  int nEvents=tree->GetEntries();
  //nEvents=1000;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);  
    
    double weight = evt_scale1fb*purewgt;
 
    if(evt_passgoodrunlist==0) continue;
 
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
    //OS
    if((lep_pdgId->at(0)*lep_pdgId->at(1)==-169)) //mm
    //if((lep_pdgId->at(0)*lep_pdgId->at(1)==-169) and (nVlep==2)*(nLlep==2)*(nTlep==2))//mm
    {
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, ptrel0, ptrel1, deltaR1, deltaR2;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = ptrel0 = ptrel1 = deltaR1 = deltaR2 = 0.0;
      TLorentzVector mu1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), MUON_MASS);
      TLorentzVector mu2 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), MUON_MASS); 
      weight *= trigSFMuLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*trigSFMuTrail(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(0)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(1)))*idSFMu(sfidMuFile, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFMu(sfidMuFile, lep_pt->at(1), fabs(lep_eta->at(1)))*sfMu(sfMuFile, lep_pt->at(0), fabs(lep_eta->at(0)))*sfMu(sfMuFile, lep_pt->at(1), fabs(lep_eta->at(1))); 
      TLorentzVector jet0_p4, jet1_p4;
      if(Sample=="MC" and mc_HLT_DoubleMu==1)
      {
        fillMuHistCollection(mumuHistCut1, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), ptrel0, ptrel1, lep_relIso04EAv2->at(0), lep_relIso04EAv2->at(1), deltaR1, deltaR2, (mu1+mu2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_MuMu_OS->Fill(1, weight); 
      }
      else if(Sample=="Data")
      {
        weight = 1.0;
        fillMuHistCollection(mumuHistCut1, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), ptrel0, ptrel1, lep_relIso04EAv2->at(0), lep_relIso04EAv2->at(1), deltaR1, deltaR2, (mu1+mu2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_MuMu_OS->Fill(1, weight);
      }
    }
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==-143))//em 
    //else if((lep_pdgId->at(0)*lep_pdgId->at(1)==-143) and (nVlep==2)*(nLlep==2)*(nTlep==2))//em
    {
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector mu1, el1;
      double weightMu1, weightEl1, muip3d, mudxy, mudz, eip3d, edxy, edz, mupTRatio, elpTRatio;
      muip3d=mudxy=mudz=eip3d=edxy=edz=mupTRatio=elpTRatio=0.0;
      weightMu1=weightEl1=1.0;
      if(abs(lep_pdgId->at(0))==13)
      {
        mu1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), MUON_MASS);
        weightMu1 = trigSFMuLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(0)))*idSFMu(sfidMuFile, lep_pt->at(0), fabs(lep_eta->at(0)))*sfMu(sfMuFile, lep_pt->at(0), fabs(lep_eta->at(0)));
        mupTRatio = lep_ptRatio->at(0);
        muip3d = fabs(lep_ip3d->at(0));
        mudxy = fabs(lep_dxy->at(0));
        mudz = fabs(lep_dz->at(0));
      }
      else if(abs(lep_pdgId->at(0))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
        weightEl1 = trigSFElLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(0), fabs(lep_eta->at(0)))*sfEl(sfElFile, sfElHist, lep_pt->at(0), fabs(lep_eta->at(0)));
        elpTRatio = lep_ptRatio->at(0);
        eip3d = fabs(lep_ip3d->at(0));
        edxy = fabs(lep_dxy->at(0));
        edz = fabs(lep_dz->at(0));
      }
      if(abs(lep_pdgId->at(1))==13)
      {
        mu1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), MUON_MASS);
        weightMu1 = trigSFMuLead(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(1)))*idSFMu(sfidMuFile, lep_pt->at(1), fabs(lep_eta->at(1)))*sfMu(sfMuFile, lep_pt->at(1), fabs(lep_eta->at(1)));
        mupTRatio = lep_ptRatio->at(1);
        muip3d = fabs(lep_ip3d->at(1));
        mudxy = fabs(lep_dxy->at(1));
        mudz = fabs(lep_dz->at(1));
      }
      else if(abs(lep_pdgId->at(1))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
        weightEl1 = trigSFElLead(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(1), fabs(lep_eta->at(1)))*sfEl(sfElFile, sfElHist, lep_pt->at(1), fabs(lep_eta->at(1)));
        elpTRatio = lep_ptRatio->at(1);
        eip3d = fabs(lep_ip3d->at(1));
        edxy = fabs(lep_dxy->at(1));
        edz = fabs(lep_dz->at(1));
      }
      weight *= weightMu1*weightEl1;
      if(Sample=="MC" and mc_HLT_MuEG==1)
      {
        fillElMuHistCollection(elmuHistCut1, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), mupTRatio, elpTRatio, muip3d, eip3d, mudxy, edxy, mudz, edz, (el1+mu1).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElMu_OS->Fill(1, weight);
      }
      else if(Sample=="Data")
      {
        weight = 1.0;
        fillElMuHistCollection(elmuHistCut1, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), mupTRatio, elpTRatio, muip3d, eip3d, mudxy, edxy, mudz, edz, (el1+mu1).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElMu_OS->Fill(1, weight);
      }
    }
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==-121))//ee
    //else if((lep_pdgId->at(0)*lep_pdgId->at(1)==-121) and (nVlep==2)*(nLlep==2)*(nTlep==2))//ee
    {
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
      TLorentzVector el2 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
      weight *= trigSFElLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*trigSFElLead(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(1), fabs(lep_eta->at(1)))*sfEl(sfElFile, sfElHist, lep_pt->at(0), fabs(lep_eta->at(0)))*sfEl(sfElFile, sfElHist, lep_pt->at(1), fabs(lep_eta->at(1)));
      if(Sample=="MC" and mc_HLT_DoubleEl_DZ_2==1)
      {
        fillElHistCollection(elelHistCut1, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), (el1+el2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElEl_OS->Fill(1, weight);
      }
      else if(Sample=="Data")
      {
        weight = 1.0;
        fillElHistCollection(elelHistCut1, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), (el1+el2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight); 
        h_TotalEvents_ElEl_OS->Fill(1, weight); 
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
      if(Sample=="MC" and mc_HLT_DoubleMu==1)
      {
        fillMuHistCollection(mumuHistCut2, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), ptrel0, ptrel1, lep_relIso04EAv2->at(0), lep_relIso04EAv2->at(1), deltaR1, deltaR2, (mu1+mu2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_MuMu_SS->Fill(1, weight);
        h_RawEvents_MuMu_SS->Fill(1);
        if(nj30>=2) h_TotalEvents_MuMu_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_MuMu_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_MuMu_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_MuMu_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_MuMu_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_MuMu_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_MuMu_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_MuMu_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_MuMu_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_MuMu_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_TotalEvents_MuMu_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_RawEvents_MuMu_SS->Fill(7);
        if(debug)
        {
          std::cout << "lepsf = " << lepsf << std::endl;
          std::cout << "trigsf = " << trigsf << std::endl;
          std::cout << "weight = " << weight << std::endl;
          std::cout << "nj30 = " << nj30 << std::endl;
          std::cout << "nb = " << nb << std::endl;
          std::cout << "Mjj = " << Mjj << std::endl;
          std::cout << "MjjL = " << MjjL << std::endl;
          std::cout << "DetajjL = " << DetajjL << std::endl;
          std::cout << "MllSS = " << MllSS << std::endl;
        }
      }
      else if(Sample=="Data")
      { 
        weight = 1.0;
        fillMuHistCollection(mumuHistCut2, mu1.Pt(), mu2.Pt(), mu1.Eta(), mu2.Eta(), mu1.Phi(), mu2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), ptrel0, ptrel1, lep_relIso04EAv2->at(0), lep_relIso04EAv2->at(1), deltaR1, deltaR2, (mu1+mu2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_MuMu_SS->Fill(1, weight);
        h_RawEvents_MuMu_SS->Fill(1);
        if(nj30>=2) h_TotalEvents_MuMu_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_MuMu_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_MuMu_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_MuMu_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_MuMu_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_MuMu_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_MuMu_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_MuMu_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_MuMu_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_MuMu_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_TotalEvents_MuMu_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_RawEvents_MuMu_SS->Fill(7);
      }
    }
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==143) and (nVlep==2)*(nLlep==2)*(nTlep==2) and passSSem==1)//em
    { 
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector mu1, el1;
      double weightMu1, weightEl1, muip3d, mudxy, mudz, eip3d, edxy, edz, mupTRatio, elpTRatio;
      muip3d=mudxy=mudz=eip3d=edxy=edz=mupTRatio=elpTRatio=0.0; 
      if(abs(lep_pdgId->at(0))==13)
      {      
        mu1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), MUON_MASS); 
        weightMu1 = trigSFMuLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(0)))*idSFMu(sfidMuFile, lep_pt->at(0), fabs(lep_eta->at(0)))*sfMu(sfMuFile, lep_pt->at(0), fabs(lep_eta->at(0)));
        mupTRatio = lep_ptRatio->at(0);
        muip3d = fabs(lep_ip3d->at(0));
        mudxy = fabs(lep_dxy->at(0));
        mudz = fabs(lep_dz->at(0));  
      }
      else if(abs(lep_pdgId->at(0))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
        weightEl1 = trigSFElLead(trigFile, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(0), fabs(lep_eta->at(0)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(0), fabs(lep_eta->at(0)))*sfEl(sfElFile, sfElHist, lep_pt->at(0), fabs(lep_eta->at(0)));
        elpTRatio = lep_ptRatio->at(0);
        eip3d = fabs(lep_ip3d->at(0));
        edxy = fabs(lep_dxy->at(0));
        edz = fabs(lep_dz->at(0));  
      }
      if(abs(lep_pdgId->at(1))==13)
      {
        mu1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), MUON_MASS);
        weightMu1 = trigSFMuLead(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFMuTrk(sfTrkMuFile, fabs(lep_eta->at(1)))*idSFMu(sfidMuFile, lep_pt->at(1), fabs(lep_eta->at(1)))*sfMu(sfMuFile, lep_pt->at(1), fabs(lep_eta->at(1)));
        mupTRatio = lep_ptRatio->at(1);
        muip3d = fabs(lep_ip3d->at(1));
        mudxy = fabs(lep_dxy->at(1));
        mudz = fabs(lep_dz->at(1));  
      }
      else if(abs(lep_pdgId->at(1))==11)
      {
        el1 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
        weightEl1 = trigSFElLead(trigFile, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFEl(sfElIDFile, idSFElHist, lep_pt->at(1), fabs(lep_eta->at(1)))*idSFElMVA(sfElMVAFile, idSFElMVAHist, lep_pt->at(1), fabs(lep_eta->at(1)))*sfEl(sfElFile, sfElHist, lep_pt->at(1), fabs(lep_eta->at(1)));
        elpTRatio = lep_ptRatio->at(1);
        eip3d = fabs(lep_ip3d->at(1));
        edxy = fabs(lep_dxy->at(1));
        edz = fabs(lep_dz->at(1));
      }
      weight *= lepsf*trigsf;
      if(Sample=="MC" and mc_HLT_MuEG==1)
      {
        fillElMuHistCollection(elmuHistCut2, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), mupTRatio, elpTRatio, muip3d, eip3d, mudxy, edxy, mudz, edz, (el1+mu1).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElMu_SS->Fill(1, weight);
        h_RawEvents_ElMu_SS->Fill(1); 
        if(nj30>=2) h_TotalEvents_ElMu_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_ElMu_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_ElMu_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_ElMu_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_ElMu_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_ElMu_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_ElMu_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_ElMu_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_ElMu_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_ElMu_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30) h_TotalEvents_ElMu_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30) h_RawEvents_ElMu_SS->Fill(7);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0) h_TotalEvents_ElMu_SS->Fill(8, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0) h_RawEvents_ElMu_SS->Fill(8);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0 and MTmax > 90.0) h_TotalEvents_ElMu_SS->Fill(9, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0 and MTmax > 90.0) h_RawEvents_ElMu_SS->Fill(9);
      }
      else if(Sample=="Data")
      {
        weight = 1.0;
        fillElMuHistCollection(elmuHistCut2, el1.Pt(), mu1.Pt(), el1.Eta(), mu1.Eta(), el1.Phi(), mu1.Phi(), mupTRatio, elpTRatio, muip3d, eip3d, mudxy, edxy, mudz, edz, (el1+mu1).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElMu_SS->Fill(1, weight);
        h_RawEvents_ElMu_SS->Fill(1);          
        if(nj30>=2) h_TotalEvents_ElMu_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_ElMu_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_ElMu_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_ElMu_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_ElMu_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_ElMu_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_ElMu_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_ElMu_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_ElMu_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_ElMu_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30) h_TotalEvents_ElMu_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30) h_RawEvents_ElMu_SS->Fill(7);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0) h_TotalEvents_ElMu_SS->Fill(8, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0) h_RawEvents_ElMu_SS->Fill(8);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0 and MTmax > 90.0) h_TotalEvents_ElMu_SS->Fill(9, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 30 and met_pt > 60.0 and MTmax > 90.0) h_RawEvents_ElMu_SS->Fill(9);
      }
    }
    else if((lep_pdgId->at(0)*lep_pdgId->at(1)==121) and (nVlep==2)*(nLlep==2)*(nTlep==2) and passSSee==1)//ee
    {
      double jet1pt, jet1eta, jet1phi, jet2pt, jet2eta, jet2phi, bjet1pt, bjet1csv, bjet2csv, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi;
      jet1pt = jet1eta = jet1phi = jet2pt = jet2eta = jet2phi = bjet1csv = bjet2csv = bjet1pt = bjet1eta = bjet1phi = bjet2pt = bjet2eta = bjet2phi = 0.0;
      TLorentzVector el1 = fillTLorentzVector(lep_pt->at(0), lep_eta->at(0), lep_phi->at(0), ELECTRON_MASS);
      TLorentzVector el2 = fillTLorentzVector(lep_pt->at(1), lep_eta->at(1), lep_phi->at(1), ELECTRON_MASS);
      weight *= lepsf*trigsf;
      if(Sample=="MC" and mc_HLT_DoubleEl_DZ_2==1)
      {
        fillElHistCollection(elelHistCut2, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), (el1+el2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElEl_SS->Fill(1, weight);
        h_RawEvents_ElEl_SS->Fill(1);
        if(nj30>=2) h_TotalEvents_ElEl_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_ElEl_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_ElEl_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_ElEl_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_ElEl_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_ElEl_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_ElEl_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_ElEl_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_ElEl_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_ElEl_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_TotalEvents_ElEl_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_RawEvents_ElEl_SS->Fill(7);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40 and met_pt > 60.0) h_TotalEvents_ElEl_SS->Fill(8, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40 and met_pt > 60.0) h_RawEvents_ElEl_SS->Fill(8);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and met_pt > 60.0 and abs(MllSS-91.1876)>10.0) h_TotalEvents_ElEl_SS->Fill(9, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and met_pt > 60.0 and abs(MllSS-91.1876)>10.0) h_RawEvents_ElEl_SS->Fill(9);
      }
      else if(Sample=="Data")
      {
        weight = 1.0;
        fillElHistCollection(elelHistCut2, el1.Pt(), el2.Pt(), el1.Eta(), el2.Eta(), el1.Phi(), el2.Phi(), lep_ptRatio->at(0), lep_ptRatio->at(1), fabs(lep_ip3d->at(0)), fabs(lep_ip3d->at(1)), fabs(lep_dxy->at(0)), fabs(lep_dxy->at(1)), fabs(lep_dz->at(0)), fabs(lep_dz->at(1)), (el1+el2).M(), met_pt, nj30, nb, jet1pt, jet2pt, jet1eta, jet2eta, jet1phi, jet2phi, bjet1csv, bjet2csv, bjet1pt, bjet1eta, bjet1phi, bjet2pt, bjet2eta, bjet2phi, weight);
        h_TotalEvents_ElEl_SS->Fill(1, weight);
        h_RawEvents_ElEl_SS->Fill(1);
        if(nj30>=2) h_TotalEvents_ElEl_SS->Fill(2, weight);
        if(nj30>=2) h_RawEvents_ElEl_SS->Fill(2);
        if(nj30>=2 and nb==0) h_TotalEvents_ElEl_SS->Fill(3, weight);
        if(nj30>=2 and nb==0) h_RawEvents_ElEl_SS->Fill(3);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_TotalEvents_ElEl_SS->Fill(4, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0) h_RawEvents_ElEl_SS->Fill(4);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_TotalEvents_ElEl_SS->Fill(5, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0) h_RawEvents_ElEl_SS->Fill(5);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_TotalEvents_ElEl_SS->Fill(6, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5) h_RawEvents_ElEl_SS->Fill(6);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_TotalEvents_ElEl_SS->Fill(7, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40) h_RawEvents_ElEl_SS->Fill(7);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40 and met_pt > 60.0) h_TotalEvents_ElEl_SS->Fill(8, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and MllSS > 40 and met_pt > 60.0) h_RawEvents_ElEl_SS->Fill(8);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and met_pt > 60.0 and abs(MllSS-91.1876)>10.0) h_TotalEvents_ElEl_SS->Fill(9, weight);
        if(nj30>=2 and nb==0 and abs(Mjj-80.)<15.0 and MjjL<400.0 and DetajjL < 1.5 and met_pt > 60.0 and abs(MllSS-91.1876)>10.0) h_RawEvents_ElEl_SS->Fill(9);
      }
    }
  }//event loop


  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  tFile->mkdir("2OSTL"); 
  tFile->mkdir("2SSTL"); 
  tFile->mkdir("2SS1JetCR");
  tFile->cd("2OSTL");
  writeHistCollection(mumuHistCut1);
  writeHistCollection(elelHistCut1);
  writeHistCollection(elmuHistCut1);
  tFile->cd("2SSTL");
  writeHistCollection(mumuHistCut2);
  writeHistCollection(elelHistCut2);
  writeHistCollection(elmuHistCut2);
  tFile->cd("2SS1JetCR");
  writeHistCollection(mumuHistCut3);
  writeHistCollection(elelHistCut3);
  writeHistCollection(elmuHistCut3);
  tFile->cd();
  h_TotalEvents_ElEl_OS->Write();
  h_TotalEvents_ElEl_SS->Write();
  h_TotalEvents_ElMu_OS->Write();
  h_TotalEvents_ElMu_SS->Write();
  h_TotalEvents_MuMu_OS->Write();
  h_TotalEvents_MuMu_SS->Write();
  h_RawEvents_ElEl_SS->Write();
  h_RawEvents_ElMu_SS->Write();
  h_RawEvents_MuMu_SS->Write();
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