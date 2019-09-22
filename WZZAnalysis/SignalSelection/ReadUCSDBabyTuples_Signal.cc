#include "ReadUCSDBabyTuples_2018.h"
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
  double compare1 = 15;
  double compare2 = 15;
  for(unsigned int i=0; i<v_leptons.size(); i++)
  {
    for(unsigned int j=0; j<v_leptons.size(); j++)
    {
      if(i!=j)//make sure not checking same lepton
      {
        if(v_leptons.at(i).id*v_leptons.at(j).id==-121 or v_leptons.at(i).id*v_leptons.at(j).id==-169) //check opposite sign pair 
        {
          //chi1_sq = pow(((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz), 2);
          pair1massDiff = ((v_leptons.at(i).lep_lv+v_leptons.at(j).lep_lv).M() - Mz);
          for(unsigned int k=0; k<v_leptons.size(); k++)
          {
            for(unsigned int l=0; l<v_leptons.size(); l++)
            {
              if(j!=l and j!=k and i!=k and i!=l)//make sure not checking same lepton
              {
                if(v_leptons.at(k).id*v_leptons.at(l).id==-121 or v_leptons.at(k).id*v_leptons.at(l).id==-169) //check second opposite sign pair 
                {
                  //chi2_sq = pow(((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz), 2);
                  pair2massDiff = ((v_leptons.at(k).lep_lv+v_leptons.at(l).lep_lv).M() - Mz);
                  if(fabs(pair1massDiff) < compare1 and fabs(pair2massDiff) < compare2) 
                  {
                    compare1 = fabs(pair1massDiff);
                    compare2 = fabs(pair2massDiff);
                    z_lep1=i;
                    z_lep2=j;
                    z_lep3=k;
                    z_lep4=l;

                    /*int idx1 = std::max(z_lep1, z_lep2);
                    int idx2 = std::min(z_lep1, z_lep2);
                    int idx3 = std::max(z_lep3, z_lep4);
                    int idx4 = std::min(z_lep3, z_lep4);

                    if(v_leptons.at(idx1).lep_lv.Pt() < 25.0 and v_leptons.at(idx2).lep_lv.Pt() < 10.0) continue;

                    if(v_leptons.at(idx3).lep_lv.Pt() < 25.0 and v_leptons.at(idx4).lep_lv.Pt() < 10.0) continue;
                    */

                    for(unsigned int m=0; m<v_leptons.size(); m++)
                    {
                      if(i!=m and j!=m and k!=m and l!=m)//make sure not checking same lepton
                      {
                        mT = sqrt(2*met_pt*v_leptons.at(m).pt*(1.0 - cos(v_leptons.at(m).phi - met_phi)));
                        w_lep=m;
                        //if(mT > 50.0) w_lep=m;
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

int ReadUCSDBabyTuples_Signal(std::string infile, std::string treeStr, std::string Sample, std::string Trigger="None")
{

  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  TFile *sfElRecoFile = new TFile("egammaEffi.txt_EGM2D_updatedAll.root");
  TH2F  *idSFElRecoHist = (TH2F*) sfElRecoFile->Get("EGamma_SF2D");
  TFile *sfElVetoFile = new TFile("2018_ElectronWPVeto_Fall17V2.root");
  TH2F  *idSFElVetoHist = (TH2F*) sfElVetoFile->Get("EGamma_SF2D"); 

  TFile *sfMuRecoFile = new TFile("EfficiencyStudies_2018_rootfiles_RunABCD_SF_ID.root"); 
  TH2F  *sfMuRecoHist = (TH2F*) sfMuRecoFile->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta"); 
  TFile *sfMuRecoLowPtFile = new TFile("EfficiencyStudies_2018_rootfiles_lowpt_RunABCD_SF_ID.root");
  TH2F  *sfMuRecoLowPtHist = (TH2F*) sfMuRecoLowPtFile->Get("NUM_MediumID_DEN_genTracks_pt_abseta"); 
  TFile *sfMuIsoFile = new TFile("EfficiencyStudies_2018_rootfiles_RunABCD_SF_ISO.root");
  TH2F  *sfMuIsoHist = (TH2F*) sfMuIsoFile->Get("NUM_LooseRelIso_DEN_MediumID_pt_abseta");

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
  std::vector<int>     *lep_mc_motherid;
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
  lep_mc_motherid = 0;
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
  tree->SetBranchAddress("lep_mc_motherid", &(lep_mc_motherid));
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

  TH1D *h_TotalEvents_5l = new TH1D("h_TotalEvents_5l", "h_TotalEvents_5l", 20, -0.5, 19.5); h_TotalEvents_5l->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut1 = new TH1D("h_TotalEvents_5l_Optimize_Cut1", "h_TotalEvents_5l_Optimize_Cut1", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut1->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut2 = new TH1D("h_TotalEvents_5l_Optimize_Cut2", "h_TotalEvents_5l_Optimize_Cut2", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut2->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut3 = new TH1D("h_TotalEvents_5l_Optimize_Cut3", "h_TotalEvents_5l_Optimize_Cut3", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut3->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut4 = new TH1D("h_TotalEvents_5l_Optimize_Cut4", "h_TotalEvents_5l_Optimize_Cut4", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut4->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut5 = new TH1D("h_TotalEvents_5l_Optimize_Cut5", "h_TotalEvents_5l_Optimize_Cut5", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut5->Sumw2();
  TH1D *h_TotalEvents_5l_Optimize_Cut6 = new TH1D("h_TotalEvents_5l_Optimize_Cut6", "h_TotalEvents_5l_Optimize_Cut6", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_Cut6->Sumw2();

  TH1D *h_TotalEvents_5l_Optimize_2D = new TH1D("h_TotalEvents_5l_Optimize_2D", "h_TotalEvents_5l_Optimize_2D", 20, -0.5, 19.5); h_TotalEvents_5l_Optimize_2D->Sumw2();

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
  TH1D *h_nJets = new TH1D("h_nJets", "h_nJets", 15, -0.5, 14.5); h_nJets->Sumw2();

  TH1D *h_pT_zCand_1 = new TH1D("h_pT_zCand_1", "h_pT_zCand_1", 1000, 0.0, 1000.0);h_pT_zCand_1->Sumw2();
  TH1D *h_pT_zCand_2 = new TH1D("h_pT_zCand_2", "h_pT_zCand_2", 1000, 0.0, 1000.0);h_pT_zCand_2->Sumw2();
  TH1D *h_pT_FifthLepton = new TH1D("h_pT_FifthLepton", "h_pT_FifthLepton", 1000, 0.0, 1000.0);h_pT_FifthLepton->Sumw2();
  TH1D *h_pT_FifthLepton_Iso = new TH1D("h_pT_FifthLepton_Iso", "h_pT_FifthLepton_Iso", 5000, 0.0, 5.0);h_pT_FifthLepton_Iso->Sumw2();
  TH1D *h_deltaR_FifthLepton_Jet = new TH1D("h_deltaR_FifthLepton_Jet", "h_deltaR_FifthLepton_Jet", 5000, 0.0, 5.0);h_deltaR_FifthLepton_Jet->Sumw2();

  TH1D *h_Iso_mT50_FifthLepton_pT30 = new TH1D("h_Iso_mT50_FifthLepton_pT30", "h_Iso_mT50_FifthLepton_pT30", 5000, 0.0, 5.0);h_Iso_mT50_FifthLepton_pT30->Sumw2();
  TH1D *h_pT_FifthLepton_mT50_Iso01 = new TH1D("h_pT_FifthLepton_mT50_Iso01", "h_pT_FifthLepton_mT50_Iso01", 1000, 0.0, 1000.0);h_pT_FifthLepton_mT50_Iso01->Sumw2();
  TH1D *h_mT_FifthLepton_pT30_Iso01 = new TH1D("h_mT_FifthLepton_pT30_Iso01", "h_mT_FifthLepton_pT30_Iso01", 500, 0.0, 500.0);h_mT_FifthLepton_pT30_Iso01->Sumw2();

  TH1D *h_mother_z_lep1 = new TH1D("h_mother_z_lep1", "h_mother_z_lep1", 1000, -0.5, 999.5); h_mother_z_lep1->Sumw2();
  TH1D *h_mother_z_lep2 = new TH1D("h_mother_z_lep2", "h_mother_z_lep2", 1000, -0.5, 999.5); h_mother_z_lep2->Sumw2();
  TH1D *h_mother_z_lep3 = new TH1D("h_mother_z_lep3", "h_mother_z_lep3", 1000, -0.5, 999.5); h_mother_z_lep3->Sumw2();
  TH1D *h_mother_z_lep4 = new TH1D("h_mother_z_lep4", "h_mother_z_lep4", 1000, -0.5, 999.5); h_mother_z_lep4->Sumw2();
  TH1D *h_mother_w_lep = new TH1D("h_mother_w_lep", "h_mother_w_lep", 1000, -0.5, 999.5); h_mother_w_lep->Sumw2();

  TH1D *h_leading_JetPt_5l = new TH1D("h_leading_JetPt_5l", "h_leading_JetPt_5l", 1000, 0.0, 1000.0); h_leading_JetPt_5l->Sumw2();
  TH1D *h_sumPt_6l = new TH1D("h_sumPt_6l", "h_sumPt_6l", 5000, 0.0, 5000.0);h_sumPt_6l->Sumw2();
  TH1D *h_triboson_Pt_6l = new TH1D("h_triboson_Pt_6l", "h_triboson_Pt_6l", 5000, 0.0, 5000.0);h_triboson_Pt_6l->Sumw2(); 

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
  int n_5l, n_pos_5l, n_neg_5l;
  n_5l=n_pos_5l=n_neg_5l=0;
  int mother_heavy_flavor, mother_pi0, mother_pi, mother_lep, mother_other;
  mother_heavy_flavor=mother_pi0=mother_pi=mother_lep=mother_other=0;
  double nEvnt_weighted = 0.0; 
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);  
   
    double weight = evt_scale1fb*getTruePUw2018(nTrueInt);
    if(evt_passgoodrunlist==0) continue;
    if(firstgoodvertex!=0) continue;
    nEvnt_weighted += weight;

    //std::cout << "nEvents_weights = " << nEvents*weight << std::endl;
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
       //muon.sf = recoSFMu_Up(sfMuRecoFile, sfMuRecoHist, lep_pt->at(imuon), lep_eta->at(imuon))*recoLowPtSFMu_Up(sfMuRecoLowPtFile, sfMuRecoLowPtHist, lep_pt->at(imuon), lep_eta->at(imuon))*isoSFMu_Up(sfMuIsoFile, sfMuIsoHist, lep_pt->at(imuon), lep_eta->at(imuon));
       muon.sf = recoSFMu_Down(sfMuRecoFile, sfMuRecoHist, lep_pt->at(imuon), lep_eta->at(imuon))*recoLowPtSFMu_Down(sfMuRecoLowPtFile, sfMuRecoLowPtHist, lep_pt->at(imuon), lep_eta->at(imuon))*isoSFMu_Down(sfMuIsoFile, sfMuIsoHist, lep_pt->at(imuon), lep_eta->at(imuon));
       //muon.sf = recoSFMu(sfMuRecoFile, sfMuRecoHist, lep_pt->at(imuon), lep_eta->at(imuon))*recoLowPtSFMu(sfMuRecoLowPtFile, sfMuRecoLowPtHist, lep_pt->at(imuon), lep_eta->at(imuon))*isoSFMu(sfMuIsoFile, sfMuIsoHist, lep_pt->at(imuon), lep_eta->at(imuon));
       muon.charge = lep_id->at(imuon)/13;
       //if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1)
       //if(abs(lep_id->at(imuon))==13 and muon.pt>10.0 and fabs(lep_eta->at(imuon)) < 2.4 and fabs(lep_relIso04DB->at(imuon)) > 0.25 and fabs(lep_sip3d->at(imuon)) < 4 and lep_isMediumPOG->at(imuon)==1) 
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
       electron.charge = lep_id->at(ielectron)/11; 
       //electron.sf = recoSFEl(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*vetoSFEl(sfElVetoFile, idSFElVetoHist, lep_pt->at(ielectron), lep_eta->at(ielectron)); 
       //electron.sf = recoSFEl_Up(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*vetoSFEl_Up(sfElVetoFile, idSFElVetoHist, lep_pt->at(ielectron), lep_eta->at(ielectron));
       electron.sf = recoSFEl_Down(sfElRecoFile, idSFElRecoHist, lep_pt->at(ielectron), lep_eta->at(ielectron))*vetoSFEl_Down(sfElVetoFile, idSFElVetoHist, lep_pt->at(ielectron), lep_eta->at(ielectron)); 
       //if(abs(lep_id->at(ielectron))==11 and electron.pt>10.0 and fabs(lep_eta->at(ielectron)) < 2.5 and lep_isCutBasedNoIsoVetoPOG->at(ielectron)==1 and fabs(lep_sip3d->at(ielectron)) < 4)  
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
      if(Sample=="MC") weight *= v_muons.at(0).sf*v_muons.at(1).sf*v_muons.at(2).sf*v_muons.at(3).sf; 
      else weight = 1.0;      

      h_TotalEvents_4Mu->Fill(1, weight);
 
      fillMuHistCollection4l(Cut1_4l, v_muons.at(0).pt, v_muons.at(1).pt, v_muons.at(2).pt, v_muons.at(3).pt, v_muons.at(0).eta, v_muons.at(1).eta, v_muons.at(2).eta, v_muons.at(3).eta, v_muons.at(0).phi, v_muons.at(1).phi, v_muons.at(2).phi, v_muons.at(3).phi, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_muons.at(2).sip3d, v_muons.at(3).sip3d, v_muons.at(0).dxy, v_muons.at(1).dxy, v_muons.at(2).dxy, v_muons.at(3).dxy, v_muons.at(0).dz, v_muons.at(1).dz, v_muons.at(2).dz, v_muons.at(3).dz, v_muons.at(0).iso, v_muons.at(1).iso, v_muons.at(2).iso, v_muons.at(3).iso, met_pt, nj, nb, (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv+v_muons.at(2).lep_lv+v_muons.at(3).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), (v_muons.at(1).lep_lv+v_muons.at(2).lep_lv).M(), (v_muons.at(2).lep_lv+v_muons.at(3).lep_lv).M(),  (v_muons.at(0).lep_lv+v_muons.at(3).lep_lv).M(), weight);
      if(nb==0) h_TotalEvents_4Mu->Fill(2, weight); 
      if(nb==0 and met_pt > 50.0) h_TotalEvents_4Mu->Fill(3, weight);
    }
    else if(v_electrons.size()==4)
    {

      if(Sample=="MC") weight *= v_electrons.at(0).sf*v_electrons.at(1).sf*v_electrons.at(2).sf*v_electrons.at(3).sf;
      else weight = 1.0;

      h_TotalEvents_4El->Fill(1, weight);
      
      fillElHistCollection4l(Cut1_4l, v_electrons.at(0).pt, v_electrons.at(1).pt, v_electrons.at(2).pt, v_electrons.at(3).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_electrons.at(2).eta, v_electrons.at(3).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_electrons.at(2).phi, v_electrons.at(3).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_electrons.at(2).sip3d, v_electrons.at(3).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_electrons.at(2).dxy, v_electrons.at(3).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_electrons.at(2).dz, v_electrons.at(3).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_electrons.at(2).iso, v_electrons.at(3).iso, met_pt, nj, nb, (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv+v_electrons.at(2).lep_lv+v_electrons.at(3).lep_lv).M(), weight);
      if(nb==0) h_TotalEvents_4El->Fill(2, weight);
      if(nb==0 and met_pt > 50.0) h_TotalEvents_4El->Fill(3, weight);
    }
    else if(v_muons.size()==2 and v_electrons.size()==2)
    {

      if(Sample=="MC") weight *= v_electrons.at(0).sf*v_electrons.at(1).sf*v_muons.at(0).sf*v_muons.at(1).sf;
      else weight = 1.0;

      int nlep = v_muons.size()+v_electrons.size();

      h_TotalEvents_2El2Mu->Fill(1.0, weight);

      fillElMuHistCollection4l(Cut1_4lElMu, v_electrons.at(0).pt, v_electrons.at(1).pt, v_muons.at(0).pt, v_muons.at(1).pt, v_electrons.at(0).eta, v_electrons.at(1).eta, v_muons.at(0).eta, v_muons.at(1).eta, v_electrons.at(0).phi, v_electrons.at(1).phi, v_muons.at(0).phi, v_muons.at(1).phi, v_electrons.at(0).sip3d, v_electrons.at(1).sip3d, v_muons.at(0).sip3d, v_muons.at(1).sip3d, v_electrons.at(0).dxy, v_electrons.at(1).dxy, v_muons.at(0).dxy, v_muons.at(1).dxy, v_electrons.at(0).dz, v_electrons.at(1).dz, v_muons.at(0).dz, v_muons.at(1).dz, v_electrons.at(0).iso, v_electrons.at(1).iso, v_muons.at(0).iso, v_muons.at(1).iso, met_pt, nj, nb, nlep, nlep, 0.0, 0.0, (v_electrons.at(0).lep_lv+v_electrons.at(1).lep_lv).M(), (v_muons.at(0).lep_lv+v_muons.at(1).lep_lv).M(), weight);
      if(nb==0) h_TotalEvents_2El2Mu->Fill(2.0, weight);
      if(nb==0 and met_pt > 50.0) h_TotalEvents_2El2Mu->Fill(3.0, weight);
    }
    if(v_leptons.size()==5)
    {
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

      //if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.10)
      //if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.06)
      if(z_lep1_temp >= 0 and z_lep2_temp >= 0 and z_lep3_temp >= 0 and z_lep4_temp >= 0 and w_lep_temp >= 0) //and nb==0)
      {
        n_5l++;
        weight *= v_leptons.at(z_lep1_temp).sf*v_leptons.at(z_lep2_temp).sf*v_leptons.at(z_lep3_temp).sf*v_leptons.at(z_lep4_temp).sf*v_leptons.at(w_lep_temp).sf;
        //rewrite code to sort 2 z candidates in pT
        /*double Zpt1 = (v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).Pt();
        double Zpt2 = (v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).Pt();
        if(Zpt1 > Zpt2) h_zmass_1->Fill((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).M(), weight);
        else h_zmass_2->Fill((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).M(), weight);
        */
        h_zmass_1->Fill((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).M(), weight);
        h_zmass_2->Fill((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).M(), weight);
        h_mT->Fill(sqrt(2*met_pt*v_leptons.at(w_lep_temp).lep_lv.Pt()*(1.0 - cos(v_leptons.at(w_lep_temp).phi - met_phi))), weight);
        double mT = sqrt(2*met_pt*v_leptons.at(w_lep_temp).lep_lv.Pt()*(1.0 - cos(v_leptons.at(w_lep_temp).phi - met_phi)));
        double iso = fabs(v_leptons.at(w_lep_temp).iso);
        h_pT_zCand_1->Fill((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).Pt(), weight);
        h_pT_zCand_2->Fill((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).Pt(), weight);
        h_pT_FifthLepton->Fill(v_leptons.at(w_lep_temp).lep_lv.Pt(), weight);
        h_pT_FifthLepton_Iso->Fill(v_leptons.at(w_lep_temp).iso, weight);
        if(v_selectedJets.size() > 0.0) h_deltaR_FifthLepton_Jet->Fill(v_selectedJets.at(0).DeltaR(v_leptons.at(z_lep1_temp).lep_lv), weight);
        if(mT > 50.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1) h_pT_FifthLepton_mT50_Iso01->Fill(v_leptons.at(w_lep_temp).lep_lv.Pt(), weight);
        //if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 30.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1) h_mT_FifthLepton_pT30_Iso01->Fill(mT, weight);
        //if(mT > 50.0 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 30.0) h_Iso_mT50_FifthLepton_pT30->Fill(v_leptons.at(w_lep_temp).iso, weight);
        
        if(fabs(v_leptons.at(w_lep_temp).iso) < 0.1) h_mT_FifthLepton_pT30_Iso01->Fill(mT, weight);
        if(mT > 50.0) h_Iso_mT50_FifthLepton_pT30->Fill(v_leptons.at(w_lep_temp).iso, weight);
        //optimization 2D
        if(mT > 10.0 and iso < 0.20) h_TotalEvents_5l_Optimize_2D->Fill(1, weight);
        if(mT > 30.0 and iso < 0.20) h_TotalEvents_5l_Optimize_2D->Fill(2, weight);
        if(mT > 50.0 and iso < 0.20) h_TotalEvents_5l_Optimize_2D->Fill(3, weight);
        //
        if(mT > 10.0 and iso < 0.15) h_TotalEvents_5l_Optimize_2D->Fill(4, weight);
        if(mT > 30.0 and iso < 0.15) h_TotalEvents_5l_Optimize_2D->Fill(5, weight);
        if(mT > 50.0 and iso < 0.15) h_TotalEvents_5l_Optimize_2D->Fill(6, weight); 
        //
        if(mT > 10.0 and iso < 0.10) h_TotalEvents_5l_Optimize_2D->Fill(7, weight);
        if(mT > 30.0 and iso < 0.10) h_TotalEvents_5l_Optimize_2D->Fill(8, weight);
        if(mT > 50.0 and iso < 0.10) h_TotalEvents_5l_Optimize_2D->Fill(9, weight);
        //
        if(mT > 10.0 and iso < 0.05) h_TotalEvents_5l_Optimize_2D->Fill(10, weight);
        if(mT > 30.0 and iso < 0.05) h_TotalEvents_5l_Optimize_2D->Fill(11, weight);
        if(mT > 50.0 and iso < 0.05) h_TotalEvents_5l_Optimize_2D->Fill(12, weight);
        //close
 
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 10) h_TotalEvents_5l_Optimize_Cut1->Fill(1, weight); 
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 20) h_TotalEvents_5l_Optimize_Cut1->Fill(2, weight);
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 30) h_TotalEvents_5l_Optimize_Cut1->Fill(3, weight);
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 40) h_TotalEvents_5l_Optimize_Cut1->Fill(4, weight); 
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 50) h_TotalEvents_5l_Optimize_Cut1->Fill(5, weight);
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 60) h_TotalEvents_5l_Optimize_Cut1->Fill(6, weight);
        if(mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 70) h_TotalEvents_5l_Optimize_Cut1->Fill(7, weight);
        //second var
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 10) h_TotalEvents_5l_Optimize_Cut2->Fill(1, weight);  
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 20) h_TotalEvents_5l_Optimize_Cut2->Fill(2, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 30) h_TotalEvents_5l_Optimize_Cut2->Fill(3, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 40) h_TotalEvents_5l_Optimize_Cut2->Fill(4, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 50) h_TotalEvents_5l_Optimize_Cut2->Fill(5, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 60) h_TotalEvents_5l_Optimize_Cut2->Fill(6, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 70) h_TotalEvents_5l_Optimize_Cut2->Fill(7, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2 and mT > 80) h_TotalEvents_5l_Optimize_Cut2->Fill(8, weight);
        //third var
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.3) h_TotalEvents_5l_Optimize_Cut3->Fill(1, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.25) h_TotalEvents_5l_Optimize_Cut3->Fill(2, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2) h_TotalEvents_5l_Optimize_Cut3->Fill(3, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.15) h_TotalEvents_5l_Optimize_Cut3->Fill(4, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1) h_TotalEvents_5l_Optimize_Cut3->Fill(5, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 10.0 and mT > 10.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.05) h_TotalEvents_5l_Optimize_Cut3->Fill(6, weight);
        //second iteration
        /*if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 5) h_TotalEvents_5l_Optimize_Cut4->Fill(1, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 10) h_TotalEvents_5l_Optimize_Cut4->Fill(2, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 15) h_TotalEvents_5l_Optimize_Cut4->Fill(3, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 20) h_TotalEvents_5l_Optimize_Cut4->Fill(4, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 25) h_TotalEvents_5l_Optimize_Cut4->Fill(5, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 30) h_TotalEvents_5l_Optimize_Cut4->Fill(6, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 35) h_TotalEvents_5l_Optimize_Cut4->Fill(7, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 40) h_TotalEvents_5l_Optimize_Cut4->Fill(8, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 45) h_TotalEvents_5l_Optimize_Cut4->Fill(9, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 50) h_TotalEvents_5l_Optimize_Cut4->Fill(10, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 55) h_TotalEvents_5l_Optimize_Cut4->Fill(11, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 60) h_TotalEvents_5l_Optimize_Cut4->Fill(12, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 65) h_TotalEvents_5l_Optimize_Cut4->Fill(13, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 70) h_TotalEvents_5l_Optimize_Cut4->Fill(14, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 80) h_TotalEvents_5l_Optimize_Cut4->Fill(15, weight);
        if(mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and v_leptons.at(w_lep_temp).lep_lv.Pt() > 100) h_TotalEvents_5l_Optimize_Cut4->Fill(16, weight);
        //second var
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 10) h_TotalEvents_5l_Optimize_Cut5->Fill(1, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 20) h_TotalEvents_5l_Optimize_Cut5->Fill(2, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 30) h_TotalEvents_5l_Optimize_Cut5->Fill(3, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 40) h_TotalEvents_5l_Optimize_Cut5->Fill(4, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 50) h_TotalEvents_5l_Optimize_Cut5->Fill(5, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 60) h_TotalEvents_5l_Optimize_Cut5->Fill(6, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 70) h_TotalEvents_5l_Optimize_Cut5->Fill(7, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1 and mT > 80) h_TotalEvents_5l_Optimize_Cut5->Fill(8, weight);
        //third var
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.3) h_TotalEvents_5l_Optimize_Cut6->Fill(1, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.25) h_TotalEvents_5l_Optimize_Cut6->Fill(2, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.2) h_TotalEvents_5l_Optimize_Cut6->Fill(3, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.15) h_TotalEvents_5l_Optimize_Cut6->Fill(4, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.1) h_TotalEvents_5l_Optimize_Cut6->Fill(5, weight);
        if(v_leptons.at(w_lep_temp).lep_lv.Pt() > 20.0 and mT > 40.0 and fabs(v_leptons.at(w_lep_temp).iso) < 0.05) h_TotalEvents_5l_Optimize_Cut6->Fill(6, weight); 
        */
        double Mz = 91.1876;
        double chi1_sq = pow(((v_leptons.at(z_lep1_temp).lep_lv+v_leptons.at(z_lep2_temp).lep_lv).M() - Mz), 2);
        double chi2_sq = pow(((v_leptons.at(z_lep3_temp).lep_lv+v_leptons.at(z_lep4_temp).lep_lv).M() - Mz), 2);
        h_min_chisq_5l->Fill(std::min(chi1_sq, chi2_sq), weight);
        h_TotalEvents_5l->Fill(4, weight);
        if(mT > 50.0 and iso < 0.1) h_TotalEvents_5l->Fill(5, weight);
        if(debug)
        {
          std::cout << "lep_mc_motherid->at(0) = " << lep_mc_motherid->at(z_lep1_temp) << std::endl;
          std::cout << "lep_mc_motherid->at(1) = " << lep_mc_motherid->at(z_lep2_temp) << std::endl;
          std::cout << "lep_mc_motherid->at(2) = " << lep_mc_motherid->at(z_lep3_temp) << std::endl;
          std::cout << "lep_mc_motherid->at(3) = " << lep_mc_motherid->at(z_lep4_temp) << std::endl;
          std::cout << "lep_mc_motherid->at(4) = " << lep_mc_motherid->at(w_lep_temp) << std::endl;
          std::cout << "event:run:lumi = " << evt << ":" << run << ":" << lumi << std::endl;
          std::cout << "*CMS4path = " << *CMS4path << std::endl;
        }
        h_mother_z_lep1->Fill(abs(lep_mc_motherid->at(z_lep1_temp)));
        h_mother_z_lep2->Fill(abs(lep_mc_motherid->at(z_lep2_temp)));
        h_mother_z_lep3->Fill(abs(lep_mc_motherid->at(z_lep3_temp)));
        h_mother_z_lep4->Fill(abs(lep_mc_motherid->at(z_lep4_temp)));
        h_mother_w_lep->Fill(abs(lep_mc_motherid->at(w_lep_temp)));
        if(abs(lep_mc_motherid->at(w_lep_temp)) > 300 and abs(lep_mc_motherid->at(w_lep_temp)) < 600) mother_heavy_flavor++;
        else if(abs(lep_mc_motherid->at(w_lep_temp)) > 100 and abs(lep_mc_motherid->at(w_lep_temp)) < 200) mother_pi0++;
        else if(abs(lep_mc_motherid->at(w_lep_temp)) > 200 and abs(lep_mc_motherid->at(w_lep_temp)) < 400) mother_pi++;
        else if(abs(lep_mc_motherid->at(w_lep_temp))==21) mother_heavy_flavor++;
        else if(abs(lep_mc_motherid->at(w_lep_temp))==11 or abs(lep_mc_motherid->at(w_lep_temp))==13) mother_lep++;
        else if(abs(lep_mc_motherid->at(w_lep_temp))==23)
        //else if(abs(lep_mc_motherid->at(w_lep_temp)) >= 23  and abs(lep_mc_motherid->at(w_lep_temp)) < 100)
        {
          mother_lep++;
          //std::cout << "abs(lep_mc_motherid->at(w_lep_temp)) = " << abs(lep_mc_motherid->at(w_lep_temp)) << std::endl;
        }
        else mother_other++;

        if(debug)
        {
          std::cout << "lep_motherIdv2->at(0) = " << lep_motherIdv2->at(z_lep1_temp) << std::endl;
          std::cout << "lep_motherIdv2->at(1) = " << lep_motherIdv2->at(z_lep2_temp) << std::endl;
          std::cout << "lep_motherIdv2->at(2) = " << lep_motherIdv2->at(z_lep3_temp) << std::endl;
          std::cout << "lep_motherIdv2->at(3) = " << lep_motherIdv2->at(z_lep4_temp) << std::endl;
          std::cout << "lep_motherIdv2->at(4) = " << lep_motherIdv2->at(w_lep_temp) << std::endl;
 
          std::cout << "v_leptons.at(z_lep1_temp).lep_lv.Pt() = " << v_leptons.at(z_lep1_temp).lep_lv.Pt() << ", v_leptons.at(z_lep1_temp).lep_lv.Eta() = " << v_leptons.at(z_lep1_temp).lep_lv.Eta() << ", v_leptons.at(z_lep1_temp).lep_lv.Phi() = " << v_leptons.at(z_lep1_temp).lep_lv.Phi() << ", v_leptons.at(z_lep1_temp).id = " << v_leptons.at(z_lep1_temp).id << std::endl;
          std::cout << "v_leptons.at(z_lep2_temp).lep_lv.Pt() = " << v_leptons.at(z_lep2_temp).lep_lv.Pt() << ", v_leptons.at(z_lep2_temp).lep_lv.Eta() = " << v_leptons.at(z_lep2_temp).lep_lv.Eta() << ", v_leptons.at(z_lep2_temp).lep_lv.Phi() = " << v_leptons.at(z_lep2_temp).lep_lv.Phi() << ", v_leptons.at(z_lep2_temp).id = " << v_leptons.at(z_lep2_temp).id << std::endl;
          std::cout << "v_leptons.at(z_lep3_temp).lep_lv.Pt() = " << v_leptons.at(z_lep3_temp).lep_lv.Pt() << ", v_leptons.at(z_lep3_temp).lep_lv.Eta() = " << v_leptons.at(z_lep3_temp).lep_lv.Eta() << ", v_leptons.at(z_lep3_temp).lep_lv.Phi() = " << v_leptons.at(z_lep3_temp).lep_lv.Phi() << ", v_leptons.at(z_lep3_temp).id = " << v_leptons.at(z_lep3_temp).id << std::endl;
          std::cout << "v_leptons.at(z_lep4_temp).lep_lv.Pt() = " << v_leptons.at(z_lep4_temp).lep_lv.Pt() << ", v_leptons.at(z_lep4_temp).lep_lv.Eta() = " << v_leptons.at(z_lep4_temp).lep_lv.Eta() << ", v_leptons.at(z_lep4_temp).lep_lv.Phi() = " << v_leptons.at(z_lep4_temp).lep_lv.Phi() << ", v_leptons.at(z_lep4_temp).id = " << v_leptons.at(z_lep4_temp).id << std::endl;
          std::cout << "v_leptons.at(w_lep_temp).lep_lv.Pt() = " << v_leptons.at(w_lep_temp).lep_lv.Pt() << ", v_leptons.at(w_lep_temp).lep_lv.Eta() = " << v_leptons.at(w_lep_temp).lep_lv.Eta() << ", v_leptons.at(w_lep_temp).lep_lv.Phi() = " << v_leptons.at(w_lep_temp).lep_lv.Phi()   << ", v_leptons.at(w_lep_temp).id = " << v_leptons.at(w_lep_temp).id << std::endl;
          std::cout << "event:run:lumi = " << evt << ":" << run << ":" << lumi << std::endl;
          std::cout << "*CMS4path = " << *CMS4path << std::endl;
          std::cout << "v_leptons.at(z_lep1_temp).charge = " << v_leptons.at(z_lep1_temp).charge << std::endl;
          std::cout << "v_leptons.at(z_lep2_temp).charge = " << v_leptons.at(z_lep2_temp).charge << std::endl;
          std::cout << "v_leptons.at(z_lep3_temp).charge = " << v_leptons.at(z_lep3_temp).charge << std::endl;
          std::cout << "v_leptons.at(z_lep4_temp).charge = " << v_leptons.at(z_lep4_temp).charge << std::endl;
          std::cout << "v_leptons.at(w_lep_temp).charge = " << v_leptons.at(w_lep_temp).charge << std::endl;
        }
        if(v_leptons.at(z_lep1_temp).charge*v_leptons.at(z_lep2_temp).charge*v_leptons.at(z_lep3_temp).charge*v_leptons.at(z_lep4_temp).charge*v_leptons.at(w_lep_temp).charge==1) n_pos_5l++;
        if(v_leptons.at(z_lep1_temp).charge*v_leptons.at(z_lep2_temp).charge*v_leptons.at(z_lep3_temp).charge*v_leptons.at(z_lep4_temp).charge*v_leptons.at(w_lep_temp).charge==-1) n_neg_5l++;
      
        /*for(unsigned int ijet=0; ijet<jets_p4->size(); ijet++)
        {
          std::cout << "Jet Pt = " << jets_p4->at(ijet).Pt() << " , Eta = " << jets_p4->at(ijet).Eta() << " , Phi = "  << jets_p4->at(ijet).Phi() << " , M = " << jets_p4->at(ijet).M() << std::endl;
        }*/
        if(jets_p4->size() > 0) h_leading_JetPt_5l->Fill(jets_p4->at(0).Pt());     
 
        if(Sample=="MC") fillLepHistCollection5l(Cut1_5l, v_leptons.at(0).pt, v_leptons.at(1).pt, v_leptons.at(2).pt, v_leptons.at(3).pt, v_leptons.at(4).pt, v_leptons.at(0).eta, v_leptons.at(1).eta, v_leptons.at(2).eta, v_leptons.at(3).eta, v_leptons.at(4).eta, v_leptons.at(0).phi, v_leptons.at(1).phi, v_leptons.at(2).phi, v_leptons.at(3).phi, v_leptons.at(4).phi, v_leptons.at(0).sip3d, v_leptons.at(1).sip3d, v_leptons.at(2).sip3d, v_leptons.at(3).sip3d, v_leptons.at(4).sip3d, v_leptons.at(0).dxy, v_leptons.at(1).dxy, v_leptons.at(2).dxy, v_leptons.at(3).dxy, v_leptons.at(4).dxy, v_leptons.at(0).dz, v_leptons.at(1).dz, v_leptons.at(2).dz, v_leptons.at(3).dz, v_leptons.at(4).dz, v_leptons.at(0).iso, v_leptons.at(1).iso, v_leptons.at(2).iso, v_leptons.at(3).iso, v_leptons.at(4).iso, met_pt, nj, nb, weight);
        else if(Sample=="Data") fillLepHistCollection5l(Cut1_5l, v_leptons.at(0).pt, v_leptons.at(1).pt, v_leptons.at(2).pt, v_leptons.at(3).pt, v_leptons.at(4).pt, v_leptons.at(0).eta, v_leptons.at(1).eta, v_leptons.at(2).eta, v_leptons.at(3).eta, v_leptons.at(4).eta, v_leptons.at(0).phi, v_leptons.at(1).phi, v_leptons.at(2).phi, v_leptons.at(3).phi, v_leptons.at(4).phi, v_leptons.at(0).sip3d, v_leptons.at(1).sip3d, v_leptons.at(2).sip3d, v_leptons.at(3).sip3d, v_leptons.at(4).sip3d, v_leptons.at(0).dxy, v_leptons.at(1).dxy, v_leptons.at(2).dxy, v_leptons.at(3).dxy, v_leptons.at(4).dxy, v_leptons.at(0).dz, v_leptons.at(1).dz, v_leptons.at(2).dz, v_leptons.at(3).dz, v_leptons.at(4).dz, v_leptons.at(0).iso, v_leptons.at(1).iso, v_leptons.at(2).iso, v_leptons.at(3).iso, v_leptons.at(4).iso, met_pt, nj, nb, 1.0);
      }
    }
    if(v_leptons.size()>=6) 
    {  
       weight *= v_leptons.at(0).sf*v_leptons.at(1).sf*v_leptons.at(2).sf*v_leptons.at(3).sf*v_leptons.at(4).sf*v_leptons.at(5).sf;
       double sumPt = v_leptons.at(0).lep_lv.Pt() + v_leptons.at(1).lep_lv.Pt() + v_leptons.at(2).lep_lv.Pt() + v_leptons.at(3).lep_lv.Pt() + v_leptons.at(4).lep_lv.Pt() + v_leptons.at(5).lep_lv.Pt();
       h_sumPt_6l->Fill(sumPt, weight);
       double tribosonPt = (v_leptons.at(0).lep_lv+v_leptons.at(1).lep_lv+v_leptons.at(2).lep_lv+v_leptons.at(3).lep_lv+v_leptons.at(4).lep_lv+v_leptons.at(5).lep_lv).Pt();
       h_triboson_Pt_6l->Fill(tribosonPt, weight);
       h_TotalEvents_6l->Fill(1, weight);
       if(sumPt>250) h_TotalEvents_6l->Fill(2, weight); 
    }   
  }//event loop
  std::cout << "nEvents_weights = " << nEvnt_weighted << std::endl;
  std::cout << "n_5l = " << n_5l << std::endl;
  std::cout << "n_pos_5l = " << n_pos_5l << std::endl;
  std::cout << "n_neg_5l = " << n_neg_5l << std::endl;
  std::cout << "mother_heavy_flavor = " << mother_heavy_flavor << std::endl;
  std::cout << "mother_pi0 = " << mother_pi0 << std::endl;
  std::cout << "mother_pi = " << mother_pi << std::endl;
  std::cout << "mother_lep = " << mother_lep << std::endl;
  std::cout << "mother_other = " << mother_other << std::endl;
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
  h_pT_zCand_1->Write();
  h_pT_zCand_2->Write();
  h_pT_FifthLepton->Write();
  h_pT_FifthLepton_Iso->Write();
  h_deltaR_FifthLepton_Jet->Write();
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
  h_TotalEvents_5l_Optimize_Cut1->Write();
  h_TotalEvents_5l_Optimize_Cut2->Write();
  h_TotalEvents_5l_Optimize_Cut3->Write();
  h_TotalEvents_5l_Optimize_Cut4->Write();
  h_TotalEvents_5l_Optimize_Cut5->Write();
  h_TotalEvents_5l_Optimize_Cut6->Write();
  h_TotalEvents_5l_Optimize_2D->Write();
  h_TotalEvents_6l->Write();
  h_sumPt_6l->Write();
  h_triboson_Pt_6l->Write();
  h_ee_ee_4El->Write();
  h_mm_mm_4Mu->Write();
  h_min_chi_sq_2el3mu->Write();
  h_min_chi_sq_3el2mu->Write();
  h_sq_5mu->Write();
  h_sq_5el->Write();
  h_4mupT->Write();
  h_4elpT->Write();
  h_nJets->Write();
  h_mother_z_lep1->Write();
  h_mother_z_lep2->Write();
  h_mother_z_lep3->Write();
  h_mother_z_lep4->Write();
  h_mother_w_lep->Write();
  h_leading_JetPt_5l->Write();
  h_Iso_mT50_FifthLepton_pT30->Write();
  h_pT_FifthLepton_mT50_Iso01->Write();
  h_mT_FifthLepton_pT30_Iso01->Write();
  tFile->Close();
  sfElRecoFile->Close();
  sfElVetoFile->Close();
  sfMuRecoFile->Close();
  sfMuRecoLowPtFile->Close();
  sfMuIsoFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0;
}
