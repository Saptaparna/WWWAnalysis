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

bool CutHLT(std::vector<leptonInfo> v_leptons, int HLT_DoubleEl_temp, int HLT_MuEG_temp, int HLT_DoubleMu_temp)
{
  bool passTrigger = false;
  for (unsigned int i=0; i<v_leptons.size(); i++)
  {
    for (unsigned int j=0; j<v_leptons.size(); j++)
    {
      if (i==j) continue;
      // Check if any of the combination of leptons pass the trigger thresholds
      // Ele 23 12
      // El23 Mu8
      // Mu23 El12
      // Mu 17 8
      // The thresholds are rounded up to 25, 15, or 10
      if (abs(v_leptons.at(i).id) == 11 and abs(v_leptons.at(j).id) == 11)
        passTrigger |= (HLT_DoubleEl_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 15);
      else if (abs(v_leptons.at(i).id) == 13 and abs(v_leptons.at(j).id) == 11)
        passTrigger |= (HLT_MuEG_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 15);
      else if (abs(v_leptons.at(i).id) == 11 and abs(v_leptons.at(j).id) == 13)
        passTrigger |= (HLT_MuEG_temp and v_leptons.at(i).pt > 25 and v_leptons.at(j).pt > 10);
      else if (abs(v_leptons.at(i).id) == 13 and abs(v_leptons.at(j).id) == 13)
        passTrigger |= (HLT_DoubleMu_temp and v_leptons.at(i).pt > 20 and v_leptons.at(j).pt > 10);
    }
  }
  return passTrigger;
}


void is5leptonZandWtag(std::vector<leptonInfo> v_leptons, double met_pt, double met_phi, int &z_lep1, int &z_lep2, int &z_lep3, int &z_lep4, int &w_lep)
{
  double chi1_sq, chi2_sq, mT;
  z_lep1=z_lep2=z_lep3=z_lep4=w_lep=-999;
  chi1_sq=chi2_sq=0.0;
  double Mz = 91.1876;
  double Mw = 80.379;
  double pair1massDiff, pair2massDiff;
  pair1massDiff=pair2massDiff=0.0;
  double compare1 = 20;
  double compare2 = 20;
  for(unsigned int i=0; i<v_leptons.size(); i++)
  {
    for(unsigned int j=0; j<v_leptons.size(); j++)
    {
      if(i!=j)//make sure not checking same lepton
      {
        if(v_leptons.at(i).id*v_leptons.at(j).id==-121 or v_leptons.at(i).id*v_leptons.at(j).id==-169) //check opposite sign pair 
        {
          chi1_sq = pow(((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz), 2);
          pair1massDiff = ((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz);
          if(fabs(pair1massDiff) <= compare1) compare1 = fabs(pair1massDiff);
          for(unsigned int k=0; k<v_leptons.size(); k++)
          {
            for(unsigned int l=0; l<v_leptons.size(); l++)
            {
              if(j!=l and j!=k and i!=k and i!=l)//make sure not checking same lepton
              {
                if(v_leptons.at(k).id*v_leptons.at(l).id==-121 or v_leptons.at(k).id*v_leptons.at(l).id==-169) //check second opposite sign pair 
                {
                  chi2_sq = pow(((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz), 2);
                  pair2massDiff = ((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz);
                  if(fabs(pair2massDiff) <= compare2) compare2 = fabs(pair2massDiff);
                  if(compare1 < 20.0 and compare2 < 20.0 )//zz tag
                  {
                    z_lep1=i;
                    z_lep2=j;
                    z_lep3=k;
                    z_lep4=l;
                    for(unsigned int m=0; m<v_leptons.size(); m++)
                    {
                      if(i!=m and j!=m and k!=m and l!=m)//make sure not checking same lepton
                      {
                        mT = sqrt(2*met_pt*v_leptons.at(m).pt*(1.0 - cos(v_leptons.at(m).phi - met_phi)));
                        //w_lep=m;
                        if(mT > 50.0) w_lep=m;
                        //if(fabs(mT-Mw)<=20) w_lep=m;
                      }//w-tag
                    }//m-loop  
                  }//zz-tag
                }//second opposite sign pair
              }//lepton index check
            }//l-loop 
          }//k-loop
        }//first opposite sign pair
      }//lepton index check
    }//j-loop
  }//i-loop
}

int ReadUCSDBabyTuples_new(std::string infile, std::string treeStr, std::string Sample, std::string Trigger="None")
{

  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  TFile *sfElRecoFile = new TFile("EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
  TH2F  *idSFElRecoHist = (TH2F*) sfElRecoFile->Get("EGamma_SF2D");
  TFile *sfElRecoLowPtFile = new TFile("EGM2D_BtoH_low_RecoSF_Legacy2016.root");
  TH2F  *idSFElRecoHist = (TH2F*) sfElRecoLowPtFile->Get("EGamma_SF2D");
  TFile *sfElVetoFile = new TFile("2016_ElectronWPVeto_Fall17V2.root");
  TH2F  *idSFElVetoHist = (TH2F*) sfElVetoFile->Get("EGamma_SF2D"); 

  TFile *sfMuRecoFile = new TFile("EfficiencyStudies_2016_rootfiles_RunBCDEF_SF_ID.root"); 
  TH2F  *sfMuRecoHist = (TH2F*) sfMuRecoFile->Get("NUM_MediumID_DEN_genTracks_eta_pt"); 
  TFile *sfMuRecoFile2 = new TFile("EfficiencyStudies_2016_rootfiles_RunGH_SF_ID.root");
  TH2F  *sfMuRecoHist2 = (TH2F*) sfMuRecoFile2->Get("NUM_MediumID_DEN_genTracks_eta_pt");      
  //TFile *sfMuRecoLowPtFile = new TFile("EfficiencyStudies_2018_rootfiles_lowpt_RunABCD_SF_ID.root");
  //TH2F  *sfMuRecoLowPtHist = (TH2F*) sfMuRecoLowPtFile->Get("NUM_MediumID_DEN_genTracks_pt_abseta"); 
  //TFile *sfMuIsoFile = new TFile("EfficiencyStudies_2018_rootfiles_RunABCD_SF_ISO.root");
  //TH2F  *sfMuIsoHist = (TH2F*) sfMuIsoFile->Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta");

  bool debug=false;

  Int_t           run;
  Int_t           lumi;
  ULong64_t       evt;
  Int_t           isData;
  Float_t         evt_scale1fb;
  Float_t         genps_weight;
  Float_t         xsec_br;
  Int_t           evt_passgoodrunlist;
  TString         *CMS4path;
  Int_t           CMS4index;
  Float_t         weight_fr_r1_f1;
  Float_t         weight_fr_r1_f2;
  Float_t         weight_fr_r1_f0p5;
  Float_t         weight_fr_r2_f1;
  Float_t         weight_fr_r2_f2;
  Float_t         weight_fr_r2_f0p5;
  Float_t         weight_fr_r0p5_f1;
  Float_t         weight_fr_r0p5_f2;
  Float_t         weight_fr_r0p5_f0p5;
  Float_t         weight_pdf_up;
  Float_t         weight_pdf_down;
  Float_t         weight_alphas_down;
  Float_t         weight_alphas_up;
  Int_t           HLT_DoubleMu;
  Int_t           HLT_DoubleEl;
  Int_t           HLT_MuEG;
  Int_t           pass_duplicate_ee_em_mm;
  Int_t           pass_duplicate_mm_em_ee;
  Float_t         gen_ht;
  std::vector<LorentzVector>  *gen_V_p4;
  std::vector<float>   *gen_V_pt;
  std::vector<float>   *gen_V_eta;
  std::vector<float>   *gen_V_phi;
  std::vector<float>   *gen_V_mass;
  std::vector<int>     *gen_V_id;
  std::vector<LorentzVector>   *gen_lep_p4;
  std::vector<float>   *gen_lep_pt;
  std::vector<float>   *gen_lep_eta;
  std::vector<float>   *gen_lep_phi;
  std::vector<float>   *gen_lep_mass;
  std::vector<int>     *gen_lep_id;
  Int_t           VHchannel;
  Int_t           Higgschannel;
  Int_t           firstgoodvertex;
  Int_t           nvtx;
  Int_t           nTrueInt;
  std::vector<LorentzVector>  *lep_p4;
  std::vector<float>   *lep_pt;
  std::vector<float>   *lep_eta;
  std::vector<float>   *lep_phi;
  std::vector<float>   *lep_energy;
  std::vector<float>   *lep_mva;
  std::vector<float>   *lep_relIso04DB;
  std::vector<float>   *lep_relIso03EA;
  std::vector<float>   *lep_relIso03EAwLep;
  std::vector<float>   *lep_ip3d;
  std::vector<float>   *lep_sip3d;
  std::vector<float>   *lep_dxy;
  std::vector<float>   *lep_dz;
  std::vector<int>     *lep_mc_id;
  std::vector<int>     *lep_motherIdv2;
  std::vector<int>     *lep_idx;
  std::vector<int>     *lep_id;
  std::vector<int>     *lep_isTightPOG;
  std::vector<int>     *lep_isMediumPOG;
  std::vector<int>     *lep_isCutBasedNoIsoVetoPOG;
  std::vector<int>     *lep_isCutBasedNoIsoLoosePOG;
  std::vector<int>     *lep_isCutBasedNoIsoMediumPOG;
  std::vector<int>     *lep_isCutBasedNoIsoTightPOG;
  std::vector<int>     *lep_isCutBasedIsoVetoPOG;
  std::vector<int>     *lep_isCutBasedIsoLoosePOG;
  std::vector<int>     *lep_isCutBasedIsoMediumPOG;
  std::vector<int>     *lep_isCutBasedIsoTightPOG;  
  Float_t         met_pt;
  Float_t         met_phi;
  Float_t         met_up_pt;
  Float_t         met_up_phi;
  Float_t         met_dn_pt;
  Float_t         met_dn_phi;
  Float_t         met_gen_pt;
  Float_t         met_gen_phi;
  std::vector<LorentzVector>  *jets_p4;
  std::vector<float>   *jets_pt;
  std::vector<float>   *jets_eta;
  std::vector<float>   *jets_phi;
  std::vector<float>   *jets_mass;
  std::vector<LorentzVector>  *jets_cen_p4;
  std::vector<float>   *jets_cen_pt;
  std::vector<float>   *jets_cen_eta;
  std::vector<float>   *jets_cen_phi;
  std::vector<float>   *jets_cen_mass;
  Int_t           nj;
  Int_t           nb;
  Int_t           nbmed;
  Float_t         ht;
  Int_t           nj_cen;
  Float_t         weight_btagsf;
  Float_t         weight_btagsf_heavy_DN;
  Float_t         weight_btagsf_heavy_UP;
  Float_t         weight_btagsf_light_DN;
  Float_t         weight_btagsf_light_UP;

  CMS4path     = 0; 
  gen_V_p4     = 0;
  gen_V_pt     = 0;
  gen_V_eta    = 0;
  gen_V_phi    = 0;
  gen_V_mass   = 0;
  gen_V_id     = 0;
  gen_lep_p4  = 0;
  gen_lep_pt   = 0;
  gen_lep_eta  = 0;
  gen_lep_phi  = 0;
  gen_lep_mass = 0;
  gen_lep_id   = 0;
  lep_p4       = 0;
  lep_pt       = 0;
  lep_eta      = 0;
  lep_phi      = 0;
  lep_energy   = 0;
  lep_mva      = 0;
  lep_relIso04DB = 0;
  lep_relIso03EA = 0;
  lep_relIso03EAwLep = 0;
  lep_ip3d = 0;
  lep_sip3d = 0;
  lep_dxy = 0;
  lep_dz  = 0;
  lep_mc_id = 0;
  lep_motherIdv2 = 0;
  lep_idx = 0;
  lep_id  = 0;
  lep_isTightPOG = 0;
  lep_isMediumPOG = 0;
  lep_isCutBasedNoIsoVetoPOG = 0;
  lep_isCutBasedNoIsoLoosePOG = 0;
  lep_isCutBasedNoIsoMediumPOG = 0;
  lep_isCutBasedNoIsoTightPOG = 0;
  lep_isCutBasedIsoVetoPOG = 0;
  lep_isCutBasedIsoLoosePOG = 0;
  lep_isCutBasedIsoMediumPOG = 0;
  lep_isCutBasedIsoTightPOG = 0;
  jets_p4 = 0;
  jets_pt = 0;
  jets_eta = 0;
  jets_phi = 0;
  jets_mass = 0;
  jets_cen_p4 = 0;
  jets_cen_pt = 0;
  jets_cen_eta = 0;
  jets_cen_phi = 0;
  jets_cen_mass = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("evt", &(evt));
  tree->SetBranchAddress("isData", &(isData));
  tree->SetBranchAddress("evt_scale1fb", &(evt_scale1fb));
  tree->SetBranchAddress("genps_weight", &(genps_weight));
  tree->SetBranchAddress("xsec_br", &(xsec_br));
  tree->SetBranchAddress("evt_passgoodrunlist", &(evt_passgoodrunlist));
  tree->SetBranchAddress("CMS4path", &(CMS4path));
  tree->SetBranchAddress("CMS4index", &(CMS4index));
  tree->SetBranchAddress("weight_fr_r1_f1", &(weight_fr_r1_f1));
  tree->SetBranchAddress("weight_fr_r1_f2", &(weight_fr_r1_f2));
  tree->SetBranchAddress("weight_fr_r1_f0p5", &(weight_fr_r1_f0p5));
  tree->SetBranchAddress("weight_fr_r2_f1", &(weight_fr_r2_f1));
  tree->SetBranchAddress("weight_fr_r2_f2", &(weight_fr_r2_f2));
  tree->SetBranchAddress("weight_fr_r2_f0p5", &(weight_fr_r2_f0p5));
  tree->SetBranchAddress("weight_fr_r0p5_f1", &(weight_fr_r0p5_f1));
  tree->SetBranchAddress("weight_fr_r0p5_f2", &(weight_fr_r0p5_f2));
  tree->SetBranchAddress("weight_fr_r0p5_f0p5", &(weight_fr_r0p5_f0p5));
  tree->SetBranchAddress("weight_pdf_up", &(weight_pdf_up));
  tree->SetBranchAddress("weight_pdf_down", &(weight_pdf_down));
  tree->SetBranchAddress("weight_alphas_down", &(weight_alphas_down));
  tree->SetBranchAddress("weight_alphas_up", &(weight_alphas_up));
  tree->SetBranchAddress("HLT_DoubleMu", &(HLT_DoubleMu));
  tree->SetBranchAddress("HLT_DoubleEl", &(HLT_DoubleEl));
  tree->SetBranchAddress("HLT_MuEG", &(HLT_MuEG));
  tree->SetBranchAddress("pass_duplicate_ee_em_mm", &(pass_duplicate_ee_em_mm));
  tree->SetBranchAddress("pass_duplicate_mm_em_ee", &(pass_duplicate_mm_em_ee));
  tree->SetBranchAddress("gen_ht", &(gen_ht));
  tree->SetBranchAddress("gen_V_p4", &(gen_V_p4));
  tree->SetBranchAddress("gen_V_pt", &(gen_V_pt));
  tree->SetBranchAddress("gen_V_eta", &(gen_V_eta));
  tree->SetBranchAddress("gen_V_phi", &(gen_V_phi));
  tree->SetBranchAddress("gen_V_mass", &(gen_V_mass));
  tree->SetBranchAddress("gen_V_id", &(gen_V_id));
  tree->SetBranchAddress("gen_lep_p4", &(gen_lep_p4));
  tree->SetBranchAddress("gen_lep_pt", &(gen_lep_pt));
  tree->SetBranchAddress("gen_lep_eta", &(gen_lep_eta));
  tree->SetBranchAddress("gen_lep_phi", &(gen_lep_phi));
  tree->SetBranchAddress("gen_lep_mass", &(gen_lep_mass));
  tree->SetBranchAddress("gen_lep_id", &(gen_lep_id));
  tree->SetBranchAddress("VHchannel", &(VHchannel));
  tree->SetBranchAddress("Higgschannel", &(Higgschannel));
  tree->SetBranchAddress("firstgoodvertex", &(firstgoodvertex));
  tree->SetBranchAddress("nvtx", &(nvtx));
  tree->SetBranchAddress("nTrueInt", &(nTrueInt));
  tree->SetBranchAddress("lep_p4", &(lep_p4));
  tree->SetBranchAddress("lep_pt", &(lep_pt));
  tree->SetBranchAddress("lep_eta", &(lep_eta));
  tree->SetBranchAddress("lep_phi", &(lep_phi));
  tree->SetBranchAddress("lep_energy", &(lep_energy));
  tree->SetBranchAddress("lep_mva", &(lep_mva));
  tree->SetBranchAddress("lep_relIso04DB", &(lep_relIso04DB));
  tree->SetBranchAddress("lep_relIso03EA", &(lep_relIso03EA));
  tree->SetBranchAddress("lep_relIso03EAwLep", &(lep_relIso03EAwLep));
  tree->SetBranchAddress("lep_ip3d", &(lep_ip3d));
  tree->SetBranchAddress("lep_sip3d", &(lep_sip3d));
  tree->SetBranchAddress("lep_dxy", &(lep_dxy));
  tree->SetBranchAddress("lep_dz", &(lep_dz));
  tree->SetBranchAddress("lep_mc_id", &(lep_mc_id));
  tree->SetBranchAddress("lep_motherIdv2", &(lep_motherIdv2));
  tree->SetBranchAddress("lep_idx", &(lep_idx));
  tree->SetBranchAddress("lep_id", &(lep_id));
  tree->SetBranchAddress("lep_isTightPOG", &(lep_isTightPOG));
  tree->SetBranchAddress("lep_isMediumPOG", &(lep_isMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoVetoPOG", &(lep_isCutBasedNoIsoVetoPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoLoosePOG", &(lep_isCutBasedNoIsoLoosePOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoMediumPOG", &(lep_isCutBasedNoIsoMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedNoIsoTightPOG", &(lep_isCutBasedNoIsoTightPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoVetoPOG", &(lep_isCutBasedIsoVetoPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoLoosePOG", &(lep_isCutBasedIsoLoosePOG));
  tree->SetBranchAddress("lep_isCutBasedIsoMediumPOG", &(lep_isCutBasedIsoMediumPOG));
  tree->SetBranchAddress("lep_isCutBasedIsoTightPOG", &(lep_isCutBasedIsoTightPOG));
  tree->SetBranchAddress("met_pt", &(met_pt));
  tree->SetBranchAddress("met_phi", &(met_phi));
  tree->SetBranchAddress("met_up_pt", &(met_up_pt));
  tree->SetBranchAddress("met_up_phi", &(met_up_phi));
  tree->SetBranchAddress("met_dn_pt", &(met_dn_pt));
  tree->SetBranchAddress("met_dn_phi", &(met_dn_phi));
  tree->SetBranchAddress("met_gen_pt", &(met_gen_pt));
  tree->SetBranchAddress("met_gen_phi", &(met_gen_phi));
  tree->SetBranchAddress("jets_p4", &(jets_p4));
  tree->SetBranchAddress("jets_pt", &(jets_pt));
  tree->SetBranchAddress("jets_eta", &(jets_eta));
  tree->SetBranchAddress("jets_phi", &(jets_phi));
  tree->SetBranchAddress("jets_mass", &(jets_mass));
  tree->SetBranchAddress("jets_cen_p4", &(jets_cen_p4));
  tree->SetBranchAddress("jets_cen_pt", &(jets_cen_pt));
  tree->SetBranchAddress("jets_cen_eta", &(jets_cen_eta));
  tree->SetBranchAddress("jets_cen_phi", &(jets_cen_phi));
  tree->SetBranchAddress("jets_cen_mass", &(jets_cen_mass));
  tree->SetBranchAddress("nj", &(nj));
  tree->SetBranchAddress("nb", &(nb));
  tree->SetBranchAddress("nbmed", &(nbmed));
  tree->SetBranchAddress("ht", &(ht));
  tree->SetBranchAddress("nj_cen", &(nj_cen));
  tree->SetBranchAddress("weight_btagsf", &(weight_btagsf));
  tree->SetBranchAddress("weight_btagsf_heavy_DN", &(weight_btagsf_heavy_DN));
  tree->SetBranchAddress("weight_btagsf_heavy_UP", &(weight_btagsf_heavy_UP));
  tree->SetBranchAddress("weight_btagsf_light_DN", &(weight_btagsf_light_DN));
  tree->SetBranchAddress("weight_btagsf_light_UP", &(weight_btagsf_light_UP));

  TH1D *h_nLeptons = new TH1D("h_nLeptons", "h_nLeptons", 10, -0.5, 9.5);h_nLeptons->Sumw2();

  TH2D *h_ee_mm_ElMu = new TH2D("h_ee_mm_ElMu", "h_ee_mm_ElMu", 300, 0.0, 300.0, 300, 0.0, 300.0); h_ee_mm_ElMu->Sumw2();
  TH2D *h_chi1_chi2 = new TH2D("h_chi1_chi2", "h_chi1_chi2", 1000, 0.0, 1000.0, 1000, 0.0, 1000.0); h_chi1_chi2->Sumw2();
  TH1D *h_min_chi = new TH1D("h_min_chi", "h_min_chi", 10000, 0.0, 1000.0); h_min_chi->Sumw2(); 
  TH1D *h_min_chi_weight = new TH1D("h_min_chi_weight", "h_min_chi_weight", 10000, 0.0, 1000.0); h_min_chi_weight->Sumw2(); 
  TH2D *h_ee_ee_4El = new TH2D("h_ee_ee_4El", "h_ee_ee_4El", 300, 0.0, 300.0, 300, 0.0, 300.0); h_ee_ee_4El->Sumw2();
  TH2D *h_mm_mm_4Mu = new TH2D("h_mm_mm_4Mu", "h_mm_mm_4Mu", 300, 0.0, 300.0, 300, 0.0, 300.0); h_mm_mm_4Mu->Sumw2();

  HistCollection Cut1_4l;
  initializeHistCollection(Cut1_4l, "4l");

  HistCollection Cut1_4lElMu;
  initializeHistCollection(Cut1_4lElMu, "4lElMu");

  HistCollection Cut1_5l;
  initializeHistCollection(Cut1_5l, "5l");
  
  HistCollection Cut1_5emu;
  initializeHistCollection(Cut1_5emu, "emu");

  TH1D *h_TotalEvents_4Mu = new TH1D("h_TotalEvents_4Mu", "h_TotalEvents_4Mu", 15, -0.5, 14.5); h_TotalEvents_4Mu->Sumw2();
  TH1D *h_TotalEvents_2El2Mu = new TH1D("h_TotalEvents_2El2Mu", "h_TotalEvents_2El2Mu", 15, -0.5, 14.5); h_TotalEvents_2El2Mu->Sumw2();
  TH1D *h_TotalEvents_4El = new TH1D("h_TotalEvents_4El", "h_TotalEvents_4El", 15, -0.5, 14.5); h_TotalEvents_4El->Sumw2();

  TH1D *h_TotalEvents_5l = new TH1D("h_TotalEvents_5l", "h_TotalEvents_5l", 15, -0.5, 14.5); h_TotalEvents_5l->Sumw2();
  TH1D *h_TotalEvents_6l = new TH1D("h_TotalEvents_6l", "h_TotalEvents_6l", 15, -0.5, 14.5); h_TotalEvents_6l->Sumw2();

  TH1D *h_deltaPhi_comb1 = new TH1D("h_deltaPhi_comb1", "h_deltaPhi_comb1", 800, -4.0, 4.0); h_deltaPhi_comb1->Sumw2();
  TH1D *h_deltaPhi_comb2 = new TH1D("h_deltaPhi_comb2", "h_deltaPhi_comb2", 800, -4.0, 4.0); h_deltaPhi_comb2->Sumw2();
  TH1D *h_deltaPhi_max = new TH1D("h_deltaPhi_max", "h_deltaPhi_max", 800, -4.0, 4.0); h_deltaPhi_max->Sumw2();

  TH1D *h_min_chi_sq_2el3mu = new TH1D("h_min_chi_sq_2el3mu", "h_min_chi_sq_2el3mu", 10000, 0.0, 1000.0); h_min_chi_sq_2el3mu->Sumw2();
  TH1D *h_min_chi_sq_3el2mu = new TH1D("h_min_chi_sq_3el2mu", "h_min_chi_sq_3el2mu", 10000, 0.0, 1000.0); h_min_chi_sq_3el2mu->Sumw2();

  TH1D *h_sq_5mu = new TH1D("h_sq_5mu", "h_sq_5mu", 10000, 0.0, 1000.0); h_sq_5mu->Sumw2();
  TH1D *h_sq_5el = new TH1D("h_sq_5el", "h_sq_5el", 10000, 0.0, 1000.0); h_sq_5el->Sumw2();

  TH1D *h_4mupT = new TH1D("h_4mupT", "h_4mupT", 1000, 0.0, 1000.0); h_4mupT->Sumw2();
  TH1D *h_4elpT = new TH1D("h_4elpT", "h_4elpT", 1000, 0.0, 1000.0); h_4elpT->Sumw2();

  TH1D *h_zmass_1 = new TH1D("h_zmass_1", "h_zmass_1", 500, 0.0, 500.0);h_zmass_1->Sumw2();
  TH1D *h_zmass_2 = new TH1D("h_zmass_2", "h_zmass_2", 500, 0.0, 500.0);h_zmass_2->Sumw2();
  TH1D *h_mT = new TH1D("h_mT", "h_mT", 500, 0.0, 500.0);h_mT->Sumw2();
  TH1D *h_min_chisq_5l = new TH1D("h_min_chisq_5l", "h_min_chisq_5l", 10000, 0.0, 1000.0);h_min_chisq_5l->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << nEvents << std::endl;
  //nEvents=10;
  int lep_ZCand_idx1 = 0;
  unsigned int lep_ZCand_idx2 = 0;
  int lep_WCand_idx1 = 0;
  int lep_2ndZCand_idx1 = 0;
  int lep_2ndZCand_idx2 = 0;
  vector <long long> checkDuplicates;
  int nDup = 0;
  long long RUNPREF = 1000 * 1000;
  RUNPREF *= 1000 * 1000;
  int n_5l;
  n_5l=0;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);  
    double weight = evt_scale1fb*getTruePUw2016(nTrueInt);
    if(evt_passgoodrunlist==0) continue;
    if(firstgoodvertex!=0) continue;

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

    std::vector<leptonInfo> v_leptons;
    v_leptons.clear();

    std::vector<leptonInfo> v_muons;
    v_muons.clear();

    for (unsigned int imuon=0; imuon<lep_p4->size(); imuon++)
    {
       leptonInfo muon;
       muon.pt = lep_pt->at(imuon); 
       muon.phi = lep_phi->at(imuon);
       muon.eta = lep_eta->at(imuon);
       muon.sip3d = lep_sip3d->at(imuon);
       muon.dxy  = lep_dxy->at(imuon);
       muon.dz = lep_dz->at(imuon);
       muon.id = lep_id->at(imuon);
       muon.lep_lv.SetPtEtaPhiM(lep_pt->at(imuon), lep_eta->at(imuon), lep_phi->at(imuon), MUON_MASS);
       muon.iso = lep_relIso04DB->at(imuon);
       double recoSF = 0.550*recoSFMu(sfMuRecoFile, sfMuRecoHist, lep_pt->at(imuon), lep_eta->at(imuon)) + 0.450*recoSFMu(sfMuRecoFile2, sfMuRecoHist2, lep_pt->at(imuon), lep_eta->at(imuon)); 
       double recoSFlowPt = 0.550*recoLowPtSFMu(sfMuRecoLowPtFile, sfMuRecoLowPtHist, lep_pt->at(imuon), lep_eta->at(imuon)) + 0.450*recoLowPtSFMu(sfMuRecoLowPtFile2, sfMuRecoLowPtHist2, lep_pt->at(imuon), lep_eta->at(imuon));
       double isoSF = 0.550*isoSFMu(sfMuIsoFile, sfMuIsoHist, lep_pt->at(imuon), lep_eta->at(imuon)) + 0.450*isoSFMu(sfMuIsoFile2, sfMuIsoHist2, lep_pt->at(imuon), lep_eta->at(imuon)); 
       muon.sf = recoSF*recoSFlowPt*isoSF; 
       if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_relIso04DB->at(imuon)) < 0.25 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1) 
       {
         v_muons.push_back(muon);
         v_leptons.push_back(muon);
       } 
    }
    std::sort (v_muons.begin(), v_muons.end(), sortLeptonsInDescendingpT);

    std::vector<leptonInfo> v_electrons;
    v_electrons.clear();    

    for (unsigned int ielectron=0; ielectron<lep_p4->size(); ielectron++)
    {
       leptonInfo electron;
       electron.pt = lep_pt->at(ielectron);
       electron.phi = lep_phi->at(ielectron);
       electron.eta = lep_eta->at(ielectron);
       electron.sip3d = lep_sip3d->at(ielectron);
       electron.dxy  = lep_dxy->at(ielectron);
       electron.dz = lep_dz->at(ielectron);
       electron.mva = lep_mva->at(ielectron);
       electron.id = lep_id->at(ielectron);
       electron.lep_lv.SetPtEtaPhiM(lep_pt->at(ielectron), lep_eta->at(ielectron), lep_phi->at(ielectron), ELECTRON_MASS);
       electron.iso = lep_relIso03EA->at(ielectron); 
       electron.sf = recoSFEl(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*recoLowPtSFEl(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*vetoSFEl(sfElVetoFile, idSFElVetoHist, lep_pt->at(ielectron), lep_eta->at(ielectron)); 
       if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5 and lep_isCutBasedIsoVetoPOG->at(ielectron)==1 and fabs(lep_sip3d->at(ielectron)) < 4) 
       {
         v_electrons.push_back(electron);
         v_leptons.push_back(electron);
       }
    }
    std::sort (v_electrons.begin(), v_electrons.end(), sortLeptonsInDescendingpT);
    
    std::sort (v_leptons.begin(), v_leptons.end(), sortLeptonsInDescendingpT);

    unsigned int nv_leptons = v_muons.size() + v_electrons.size();
    if(nv_leptons != v_leptons.size()) std::cout << "size mismatch" << std::endl;   

    if(not CutHLT(v_leptons, HLT_DoubleEl, HLT_MuEG, HLT_DoubleMu)) continue;
 
    if(Sample=="MC") h_nLeptons->Fill(nv_leptons, weight);
    else if(Sample=="Data") 
    {
      long long dupCheck = run*RUNPREF + evt;
      bool bDuplicate = false;
      for (unsigned int uid = 0; uid < checkDuplicates.size(); uid++)
      {
        if (checkDuplicates[uid] == dupCheck)
        {
          cout<<dupCheck<<endl;
          bDuplicate = true;
          nDup++;
          break;
        }
      }
      if (bDuplicate) continue;
      else checkDuplicates.push_back(dupCheck);
      h_nLeptons->Fill(nv_leptons, 1.0); 
    }
 
    if(v_muons.size()==4)
    {
      TLorentzVector mu1, mu2, mu3, mu4;
      mu1.SetPtEtaPhiM(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, 0.105);
      mu2.SetPtEtaPhiM(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, 0.105);
      mu3.SetPtEtaPhiM(v_muons.at(2).pt, v_muons.at(2).eta, v_muons.at(2).phi, 0.105);
      mu4.SetPtEtaPhiM(v_muons.at(3).pt, v_muons.at(3).eta, v_muons.at(3).phi, 0.105);

      h_4mupT->Fill((mu1+mu2+mu3+mu4).Pt());     
 
      bool ifpass;
      double compare = 15; //Min value of |Mll - MZ|
      int lep_zcand_idx1, lep_zcand_idx2, lep_nonzcand_idx1, lep_nonzcand_idx2;
      lep_zcand_idx1 = lep_zcand_idx2 = lep_nonzcand_idx1 = lep_nonzcand_idx2 = -999; 
      for (unsigned int jj = 0 ; jj < (v_muons.size() - 1) ; jj ++)
      {
        for (unsigned int kk = 0 ; kk < v_muons.size() ; kk ++)
        {
          if (v_muons.at(jj).id == v_muons.at(kk).id) continue;
          TLorentzVector zcand = v_muons.at(jj).lep_lv + v_muons.at(kk).lep_lv; 
          ifpass = true; //based on Philip's classification scheme
          if (fabs(zcand.M() - 91.1875) > 15)  ifpass = false; // within Z mass window
          if (abs(v_muons.at(jj).id) != abs(v_muons.at(kk).id)) ifpass = false; // same-flavor
          if (v_muons.at(jj).id == v_muons.at(kk).id) ifpass = false; // opposite-sign
          if (ifpass && fabs(zcand.M() - 91.1875) > compare) ifpass = false;  
          if (ifpass)
          {
            compare = fabs(zcand.M() - 91.1876);
            lep_zcand_idx1 = jj;
            lep_zcand_idx2 = kk;
          }
        }  
      }
      double Mz = 91.1876;
      double chi1_sq = pow(((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M() - Mz), 2) + pow(((v_muons.at(2).lep_lv+v_muons.at(3).lep_lv).M() - Mz), 2);  
      double chi2_sq = pow(((v_muons.at(0).lep_lv+v_muons.at(2).lep_lv).M() - Mz), 2) + pow(((v_muons.at(1).lep_lv+v_muons.at(3).lep_lv).M() - Mz), 2);
      
      h_chi1_chi2->Fill(chi1_sq, chi2_sq);       
      h_min_chi->Fill(std::min(chi1_sq, chi2_sq));
      TLorentzVector lv_comb01 = (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv);
      TLorentzVector lv_comb23 = (v_muons.at(2).lep_lv+v_muons.at(3).lep_lv);
      TLorentzVector lv_comb02 = (v_muons.at(0).lep_lv+v_muons.at(2).lep_lv);
      TLorentzVector lv_comb13 = (v_muons.at(1).lep_lv+v_muons.at(3).lep_lv);

      h_deltaPhi_comb1->Fill(lv_comb01.DeltaPhi(lv_comb23), weight);
      h_deltaPhi_comb2->Fill(lv_comb02.DeltaPhi(lv_comb13), weight);
      h_deltaPhi_max->Fill(std::max(lv_comb01.DeltaPhi(lv_comb23), lv_comb02.DeltaPhi(lv_comb13)), weight);

      h_mm_mm_4Mu->Fill((v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), (v_muons.at(2).lep_lv+v_muons.at(3).lep_lv).M());

      if(Sample=="MC") h_TotalEvents_4Mu->Fill(1, weight);
      else if(Sample=="Data") h_TotalEvents_4Mu->Fill(1, 1.0);
      double minChi2 = std::min(chi1_sq, chi2_sq); 
      
      if(nb==0) h_TotalEvents_4Mu->Fill(2, weight);

      if(minChi2 < 100.0 and nb==0 and nj==0 and fabs(lv_comb01.DeltaPhi(lv_comb23)) < 2.5 and fabs(lv_comb02.DeltaPhi(lv_comb13)) < 2.5)
      {
        if(Sample=="MC") 
        {
          h_min_chi_weight->Fill(minChi2, weight);
          fillMuHistCollection4l(Cut1_4l, v_muons.at(0).pt, v_muons.at(1).pt, v_muons.at(2).pt, v_muons.at(3).pt, v_muons.at(0).eta, v_muons.at(1).eta, v_muons.at(2).eta, v_muons.at(3).eta, v_muons.at(0).phi, v_muons.at(1).phi, v_muons.at(2).phi, v_muons.at(3).phi, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_muons.at(2).sip3d, v_muons.at(3).sip3d, v_muons.at(0).dxy, v_muons.at(1).dxy, v_muons.at(2).dxy, v_muons.at(3).dxy, v_muons.at(0).dz, v_muons.at(1).dz, v_muons.at(2).dz, v_muons.at(3).dz, met_pt, nj, nb, (mu1+mu2+mu3+mu4).M(), (mu1+mu2).M(), (mu2+mu3).M(), (mu3+mu4).M(),  (mu4+mu1).M(), weight);
          h_TotalEvents_4Mu->Fill(3, weight); 
          if(met_pt > 50.0)
          { 
            h_TotalEvents_4Mu->Fill(4, weight);
          }
        }
        else if(Sample=="Data") 
        {
          h_min_chi_weight->Fill(minChi2, 1.0);
          fillMuHistCollection4l(Cut1_4l, v_muons.at(0).pt, v_muons.at(1).pt, v_muons.at(2).pt, v_muons.at(3).pt, v_muons.at(0).eta, v_muons.at(1).eta, v_muons.at(2).eta, v_muons.at(3).eta, v_muons.at(0).phi, v_muons.at(1).phi, v_muons.at(2).phi, v_muons.at(3).phi, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_muons.at(2).sip3d, v_muons.at(3).sip3d, v_muons.at(0).dxy, v_muons.at(1).dxy, v_muons.at(2).dxy, v_muons.at(3).dxy, v_muons.at(0).dz, v_muons.at(1).dz, v_muons.at(2).dz, v_muons.at(3).dz, met_pt, nj, nb, (mu1+mu2+mu3+mu4).M(), (mu1+mu2).M(), (mu2+mu3).M(), (mu3+mu4).M(),  (mu4+mu1).M(), 1.0);
          h_TotalEvents_4Mu->Fill(3, 1.0);
          if(met_pt > 50.0)
          { 
            h_TotalEvents_4Mu->Fill(4, 1.0);
          }
        }
      }
    }
    else if(v_electrons.size()==4)
    {
      TLorentzVector el1, el2, el3, el4;
      el1.SetPtEtaPhiM(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, 0.000511);
      el2.SetPtEtaPhiM(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, 0.000511);
      el3.SetPtEtaPhiM(v_electrons.at(2).pt, v_electrons.at(2).eta, v_electrons.at(2).phi, 0.000511);
      el4.SetPtEtaPhiM(v_electrons.at(3).pt, v_electrons.at(3).eta, v_electrons.at(3).phi, 0.000511);     

      bool ifpass;
      double compare = 15; //Min value of |Mll - MZ|
      int lep_zcand_idx1, lep_zcand_idx2, lep_nonzcand_idx1, lep_nonzcand_idx2;
      lep_zcand_idx1 = lep_zcand_idx2 = lep_nonzcand_idx1 = lep_nonzcand_idx2 = -999;
      for (unsigned int jj = 0 ; jj < (v_electrons.size() - 1) ; jj ++)
      {
        for (unsigned int kk = 0 ; kk < v_electrons.size() ; kk ++)
        {
          TLorentzVector zcand = v_electrons.at(jj).lep_lv + v_electrons.at(kk).lep_lv;
          ifpass = true; //based on Philip's classification scheme
          if (fabs(zcand.M() - 91.1875) > 15)  ifpass = false; // within Z mass window
          if (abs(v_electrons.at(jj).id) != abs(v_electrons.at(kk).id)) ifpass = false; // same-flavor
          if (v_electrons.at(jj).id == v_electrons.at(kk).id) ifpass = false; // opposite-sign
          if (ifpass && fabs(zcand.M() - 91.1875) > compare) ifpass = false;
          if (ifpass)
          {
            compare = fabs(zcand.M() - 91.1876);
            lep_zcand_idx1 = jj;
            lep_zcand_idx2 = kk;
          }
        }
      }

      double Mz = 91.1876;
      double chi1_sq = pow(((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M() - Mz), 2) + pow(((v_electrons.at(2).lep_lv+v_electrons.at(3).lep_lv).M() - Mz), 2);
      double chi2_sq = pow(((v_electrons.at(0).lep_lv+v_electrons.at(2).lep_lv).M() - Mz), 2) + pow(((v_electrons.at(1).lep_lv+v_electrons.at(3).lep_lv).M() - Mz), 2);
      h_chi1_chi2->Fill(chi1_sq, chi2_sq);
      h_min_chi->Fill(std::min(chi1_sq, chi2_sq));

      TLorentzVector lv_comb01 = (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv);
      TLorentzVector lv_comb23 = (v_electrons.at(2).lep_lv+v_electrons.at(3).lep_lv);
      TLorentzVector lv_comb02 = (v_electrons.at(0).lep_lv+v_electrons.at(2).lep_lv);
      TLorentzVector lv_comb13 = (v_electrons.at(1).lep_lv+v_electrons.at(3).lep_lv);

      h_deltaPhi_comb1->Fill(lv_comb01.DeltaPhi(lv_comb23), weight);
      h_deltaPhi_comb2->Fill(lv_comb02.DeltaPhi(lv_comb13), weight);
      h_deltaPhi_max->Fill(std::max(lv_comb01.DeltaPhi(lv_comb23), lv_comb02.DeltaPhi(lv_comb13)), weight);

      h_ee_ee_4El->Fill((v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M(), (v_electrons.at(2).lep_lv+v_electrons.at(3).lep_lv).M());

      if(Sample=="MC") h_TotalEvents_4El->Fill(1, weight);
      else if(Sample=="Data") h_TotalEvents_4El->Fill(1, 1.0);
      
      if(nb==0) h_TotalEvents_4El->Fill(2, weight);
      
      double minChi2 = std::min(chi1_sq, chi2_sq); 
      if(minChi2 < 100.0 and nb==0 and nj==0 and fabs(lv_comb01.DeltaPhi(lv_comb23)) < 2.5 and fabs(lv_comb02.DeltaPhi(lv_comb13)) < 2.5)
      //if(minChi2 < 100.0 and nb==0)
      {
        if(Sample=="MC") 
        {
          h_min_chi_weight->Fill(minChi2, weight);
          fillElHistCollection4l(Cut1_4l, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_electrons.at(3).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_electrons.at(3).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_electrons.at(3).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_electrons.at(3).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_electrons.at(3).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_electrons.at(3).dz, met_pt, nj, nb, (el1+el2+el3+el4).M(), weight);
          h_TotalEvents_4El->Fill(3, weight);
          if(met_pt > 50.0) 
          {
            h_TotalEvents_4El->Fill(4, weight);
          } 
        }
        else if(Sample=="Data") 
        {
          h_min_chi_weight->Fill(minChi2, 1.0);
          fillElHistCollection4l(Cut1_4l, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_electrons.at(3).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_electrons.at(3).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_electrons.at(3).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_electrons.at(3).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_electrons.at(3).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_electrons.at(3).dz, met_pt, nj, nb, (el1+el2+el3+el4).M(), 1.0);
          h_TotalEvents_4El->Fill(3, 1.0);
          if(met_pt > 50.0)
          { 
            h_TotalEvents_4El->Fill(4, 1.0);
          }
        }
      }
    }
    if(v_muons.size()==2 and v_electrons.size()==2)
    {
      TLorentzVector el1, el2, mu1, mu2;
      mu1.SetPtEtaPhiM(v_muons.at(0).pt, v_muons.at(0).eta, v_muons.at(0).phi, 0.105);
      mu2.SetPtEtaPhiM(v_muons.at(1).pt, v_muons.at(1).eta, v_muons.at(1).phi, 0.105);
      el1.SetPtEtaPhiM(v_electrons.at(0).pt, v_electrons.at(0).eta, v_electrons.at(0).phi, 0.000511);
      el2.SetPtEtaPhiM(v_electrons.at(1).pt, v_electrons.at(1).eta, v_electrons.at(1).phi, 0.000511);
     
      if(Sample=="MC") h_TotalEvents_2El2Mu->Fill(1, weight);
      else if(Sample=="Data") h_TotalEvents_2El2Mu->Fill(1, 1.0); 
 
      if(nb==0) h_TotalEvents_2El2Mu->Fill(2, weight);

      if(fabs((el1+el2).M() - 91.1875) > 10) continue;
      if(fabs((mu1+mu2).M() - 91.1875) > 10) continue;

      if(nb!=0) continue; 

      if(Sample=="MC") 
      {
        fillElMuHistCollection4l(Cut1_4lElMu, v_electrons.at(0).pt, v_electrons.at(1).pt, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_muons.at(0).dz, v_muons.at(1).dz, met_pt, nj, nb, (el1+el2+mu1+mu2).M(), (el1+el2).M(), (mu1+mu2).M(), weight);
        h_TotalEvents_2El2Mu->Fill(3, weight);
        if(met_pt > 40.0)
        { 
          h_TotalEvents_2El2Mu->Fill(4, weight);
        }
      } 
      else if(Sample=="Data") 
      {
        fillElMuHistCollection4l(Cut1_4lElMu, v_electrons.at(0).pt, v_electrons.at(1).pt, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_muons.at(0).dz, v_muons.at(1).dz, met_pt, nj, nb, (el1+el2+mu1+mu2).M(), (el1+el2).M(), (mu1+mu2).M(), 1.0);
        h_TotalEvents_2El2Mu->Fill(3, 1.0);
        if(met_pt > 40.0)
        { 
          h_TotalEvents_2El2Mu->Fill(4, 1.0);
        }
      }
      h_ee_mm_ElMu->Fill((el1+el2).M(), (mu1+mu2).M());
    }
    if(v_leptons.size()==5)
    {
      //n_5l++;
      //weight *= v_leptons.at(0).sf*v_leptons.at(1).sf*v_leptons.at(2).sf*v_leptons.at(3).sf*v_leptons.at(4).sf;
      h_TotalEvents_5l->Fill(1, weight);
      int z_lep1_temp, z_lep2_temp, z_lep3_temp, z_lep4_temp, w_lep_temp;
      z_lep1_temp=z_lep2_temp=z_lep3_temp=z_lep4_temp=w_lep_temp=-999; 
      is5leptonZandWtag(v_leptons, met_pt, met_phi, z_lep1_temp, z_lep2_temp, z_lep3_temp, z_lep4_temp, w_lep_temp);
      if(debug)
      {
        std::cout << "z_lep1 = " << z_lep1_temp << std::endl;
        std::cout << "z_lep2 = " << z_lep2_temp << std::endl;
        std::cout << "z_lep3 = " << z_lep3_temp << std::endl;
        std::cout << "z_lep4 = " << z_lep4_temp << std::endl;
        std::cout << "w_lep = " << w_lep_temp << std::endl;
      }
      double Mw = 80.379;
      int wcand_idx = -999;
      double compare = 20.0;
      for (unsigned int ilep=0; ilep < v_leptons.size(); ilep++)
      {
        double mT = sqrt(2*met_pt*v_leptons.at(ilep).pt*(1.0 - cos(v_leptons.at(ilep).phi - met_phi)));
        if(fabs(mT-Mw) <= compare) 
        {
          compare = fabs(mT-Mw);
          wcand_idx = ilep;  
        }
      }

      if(wcand_idx>=0) h_TotalEvents_5l->Fill(2, weight);

      if(wcand_idx>=0 and fabs(v_leptons.at(wcand_idx).iso) < 0.10) h_TotalEvents_5l->Fill(3, weight);

      if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0)
      {
        weight *= v_leptons.at(z_lep1_temp).sf*v_leptons.at(z_lep2_temp).sf*v_leptons.at(z_lep3_temp).sf*v_leptons.at(z_lep4_temp).sf*v_leptons.at(w_lep_temp).sf;
        h_zmass_1->Fill((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).M(), weight);
        h_zmass_2->Fill((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).M(), weight);
        h_mT->Fill(sqrt(2*met_pt*v_leptons.at(w_lep_temp).pt*(1.0 - cos(v_leptons.at(w_lep_temp).phi - met_phi))), weight);
        double Mz = 91.1876;
        double chi1_sq = pow(((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).M() - Mz), 2);
        double chi2_sq = pow(((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).M() - Mz), 2);
        h_min_chisq_5l->Fill(std::min(chi1_sq, chi2_sq), weight);
        h_TotalEvents_5l->Fill(4, weight);
      }
      
      //if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0) h_TotalEvents_5l->Fill(4, weight); 
      
      if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0) n_5l++; 
      
      if(Sample=="MC") fillLepHistCollection5l(Cut1_5l, v_leptons.at(0).pt, v_leptons.at(1).pt, v_leptons.at(2).pt, v_leptons.at(3).pt, v_leptons.at(4).pt, v_leptons.at(0).eta, v_leptons.at(1).eta, v_leptons.at(2).eta, v_leptons.at(3).eta, v_leptons.at(4).eta, v_leptons.at(0).phi, v_leptons.at(1).phi, v_leptons.at(2).phi, v_leptons.at(3).phi, v_leptons.at(4).phi, v_leptons.at(0).sip3d, v_leptons.at(1).sip3d, v_leptons.at(2).sip3d, v_leptons.at(3).sip3d, v_leptons.at(4).sip3d, v_leptons.at(0).dxy, v_leptons.at(1).dxy, v_leptons.at(2).dxy, v_leptons.at(3).dxy, v_leptons.at(4).dxy, v_leptons.at(0).dz, v_leptons.at(1).dz, v_leptons.at(2).dz, v_leptons.at(3).dz, v_leptons.at(4).dz, met_pt, nj, nb, weight);
      else if(Sample=="Data") fillLepHistCollection5l(Cut1_5l, v_leptons.at(0).pt, v_leptons.at(1).pt, v_leptons.at(2).pt, v_leptons.at(3).pt, v_leptons.at(4).pt, v_leptons.at(0).eta, v_leptons.at(1).eta, v_leptons.at(2).eta, v_leptons.at(3).eta, v_leptons.at(4).eta, v_leptons.at(0).phi, v_leptons.at(1).phi, v_leptons.at(2).phi, v_leptons.at(3).phi, v_leptons.at(4).phi, v_leptons.at(0).sip3d, v_leptons.at(1).sip3d, v_leptons.at(2).sip3d, v_leptons.at(3).sip3d, v_leptons.at(4).sip3d, v_leptons.at(0).dxy, v_leptons.at(1).dxy, v_leptons.at(2).dxy, v_leptons.at(3).dxy, v_leptons.at(4).dxy, v_leptons.at(0).dz, v_leptons.at(1).dz, v_leptons.at(2).dz, v_leptons.at(3).dz, v_leptons.at(4).dz, met_pt, nj, nb, 1.0);
    }
    if(v_leptons.size()>=6) h_TotalEvents_6l->Fill(1, weight); 
    
  }//event loop
  std::cout << "n_5l = " << n_5l << std::endl;
  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  tFile->cd();
  tFile->mkdir("Cut1_4l");
  tFile->cd("Cut1_4l");
  writeHistCollection(Cut1_4l);
  tFile->cd();
  tFile->mkdir("Cut1_4lElMu");
  tFile->cd("Cut1_4lElMu");
  writeHistCollection(Cut1_4lElMu);
  tFile->cd();
  tFile->mkdir("Cut1_5l");
  tFile->cd("Cut1_5l");
  writeHistCollection(Cut1_5l);
  tFile->cd();
  h_zmass_1->Write();
  h_zmass_2->Write();
  h_mT->Write();
  h_min_chisq_5l->Write();
  h_deltaPhi_comb1->Write();
  h_deltaPhi_comb2->Write();
  h_deltaPhi_max->Write();
  h_nLeptons->Write();
  h_ee_mm_ElMu->Write(); 
  h_chi1_chi2->Write();
  h_min_chi->Write();
  h_min_chi_weight->Write();
  h_TotalEvents_4Mu->Write();
  h_TotalEvents_2El2Mu->Write(); 
  h_TotalEvents_4El->Write(); 
  h_TotalEvents_5l->Write();
  h_TotalEvents_6l->Write();
  h_ee_ee_4El->Write();
  h_mm_mm_4Mu->Write();
  h_min_chi_sq_2el3mu->Write();
  h_min_chi_sq_3el2mu->Write();
  h_sq_5mu->Write();
  h_sq_5el->Write();
  h_4mupT->Write();
  h_4elpT->Write();
  tFile->Close();
  sfElRecoFile->Close();
  sfElVetoFile->Close();
  sfMuRecoFile->Close();
  sfMuRecoLowPtFile->Close();
  sfMuIsoFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;
}
