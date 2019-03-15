#include "jec.h"

using namespace std;

//####################################################################################
CoreUtil::jec::jec()
{
    jet_corrector_pfL1FastJetL2L3 = NULL;
    jecUnc                        = NULL;
    // corrections for later runs in 2016F
    jet_corrector_pfL1FastJetL2L3_postrun278802 = NULL;
    jecUnc_postrun278802                        = NULL;
    jetcorr_filenames_pfL1FastJetL2L3 = std::vector<std::string>();
    jetcorr_filenames_pfL1FastJetL2L3_postrun278802 = std::vector<std::string>();
    // jet uncertainty filename
    jetunc_filename_pfL1FastJetL2L3 = "";

    TString path = gSystem->Getenv("COREDIR");

    // If path is null
    if (path.IsNull())
    {
        printf("CoreUtil::jec::jec() Must set COREDIR=/path/to/CORE\n");
        abort();
    }

    // Append the fullpath via hardcode
    jecdatapath = path;
    jecdatapath += "/Tools/jetcorr/data/";

    printf("CoreUtil::jec::jec() Loading %s\n", path.Data());

}

//####################################################################################
CoreUtil::jec::~jec()
{
    if (jet_corrector_pfL1FastJetL2L3)               delete jet_corrector_pfL1FastJetL2L3;
    if (jet_corrector_pfL1FastJetL2L3_postrun278802) delete jet_corrector_pfL1FastJetL2L3_postrun278802;
    if (jecUnc)                                      delete jecUnc;
    if (jecUnc_postrun278802)                        delete jecUnc_postrun278802;
}

//####################################################################################
void CoreUtil::jec::resetObjects()
{
    if (jetcorr_filenames_pfL1FastJetL2L3_postrun278802.size() > 0)
        jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
    if (jetcorr_filenames_pfL1FastJetL2L3.size() > 0)
        jetcorr_filenames_pfL1FastJetL2L3.clear();
}

//####################################################################################
FactorizedJetCorrector* CoreUtil::jec::getJECL1(bool isData, int run)
{
    if (isData && run >= 278802 && run <= 278808)
    {
        return jet_corrector_pfL1_postrun278802;
    }
    else
    {
        return jet_corrector_pfL1;
    }
}

//####################################################################################
FactorizedJetCorrector* CoreUtil::jec::getJEC(bool isData, int run)
{
    if (isData && run >= 278802 && run <= 278808)
    {
        return jet_corrector_pfL1FastJetL2L3_postrun278802;
    }
    else
    {
        return jet_corrector_pfL1FastJetL2L3;
    }
}

//####################################################################################
JetCorrectionUncertainty* CoreUtil::jec::getJECUnc(bool isData, int run)
{
    if (isData && run >= 278802 && run <= 278808)
    {
        return jecUnc_postrun278802;;
    }
    else
    {
        return jecUnc;
    }
}

//####################################################################################
void CoreUtil::jec::setJECFor(TString filename, bool isfastsim)
{
    // If current filename same as the one already set to then skip
    if (current_filename == filename)
        return;

    // Reset
    resetObjects();

//    cout << endl;
    cout << "CoreUtil:: Create new JEC for " << filename << endl;

    // Set the current_filename to the provided filename
    current_filename = filename;

    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2015/2016
    //-----------------------------------------------------------------------------------------------------------------------------------
    //JECs for 76X
    if (filename.Contains("16Dec2015") || filename.Contains("76X_mcRun2"))
    {
        if (filename.Contains("Run2015C") || filename.Contains("Run2015D"))
        {
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else
        {
            // files for 76X MC
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt");
        }
    }
    else if (filename.Contains("80MiniAODv") || filename.Contains("RelVal"))
    {
        // files for 80X MC, ICHEP production
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_Uncertainty_AK4PFchs.txt");
    }
    else if (filename.Contains("Summer16") || filename.Contains("TEST") || filename.Contains("PrivateWWW") || filename.Contains("TChiWH_HToVVTauTau_HToBB_mChargino200_mLSP50") || filename.Contains("WWW_BW") || filename.Contains("GluGlu"))
    {
        // files for 80X MC, Summer16 (Moriond17) production
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt");
    }
    else if (filename.Contains("Run2016") || filename.Contains("CMSSW_8_0_11_V08-00-06"))
    {
        // // This scheme will eventually be used for Moriond. For now just use JECs from ICHEP
        if (filename.Contains("Run2016B") ||
            filename.Contains("Run2016C") ||
            filename.Contains("Run2016D"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK4PFchs.txt");
        }
        if (filename.Contains("Run2016E") ||
            filename.Contains("Run2016F"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
        if (filename.Contains("Run2016G"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
        if (filename.Contains("Run2016H"))
        {
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2017
    //-----------------------------------------------------------------------------------------------------------------------------------
    else if (filename.Contains("Run2017"))
    {
        if (filename.Contains("Run2017B"))
        {
            setJECFileNames("Fall17_17Nov2017B_V6", "DATA");
        }
        if (filename.Contains("Run2017C"))
        {
            setJECFileNames("Fall17_17Nov2017C_V6", "DATA");
        }
        if (filename.Contains("Run2017D"))
        {
            setJECFileNames("Fall17_17Nov2017D_V6", "DATA");
        }
        if (filename.Contains("Run2017E"))
        {
            setJECFileNames("Fall17_17Nov2017E_V6", "DATA");
        }
        if (filename.Contains("Run2017F"))
        {
            setJECFileNames("Fall17_17Nov2017F_V6", "DATA");
        }
    }
    else if (filename.Contains("Fall17MiniAOD"))
    {
        setJECFileNames("Fall17_17Nov2017_V6", "MC");
    }
    else if (filename.Contains("PhotophobicALP"))
    {
        setJECFileNames("Fall17_17Nov2017_V6", "MC");
    }

    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2016 fast sim
    //-----------------------------------------------------------------------------------------------------------------------------------
    if (isfastsim)
    {
        // files for 25ns fastsim samples
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_Uncertainty_AK4PFchs.txt");
    }


    //===================================================================================================================================
    // Now create the jet corrector based on what it is set to
    //===================================================================================================================================

    if (jetcorr_filenames_pfL1FastJetL2L3.size() == 0)
    {
        cout << "CoreUtil:: Error, sample not found. Check the JECs." << endl;
        exit(100);
    }

    // Print what is used
    cout << "CoreUtil:: JECs used:" << endl;
    for (size_t jecind = 0; jecind < jetcorr_filenames_pfL1FastJetL2L3.size(); jecind++)
    {
        cout << "CoreUtil:: " << jetcorr_filenames_pfL1FastJetL2L3.at(jecind) << endl;
    }

    // Now create
    jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
    jet_corrector_pfL1  = makeJetCorrector(jetcorr_filenames_pfL1);
    jet_corrector_pfL1_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1_postrun278802);
    // For jet uncertainty for 2016 the code was written in bad ways such that the jetunc was created on the fly above in the function
    // Starting 2017 I tried to clean up a little bit so that it is created at the end
    if (!jecUnc)
        jecUnc = new JetCorrectionUncertainty(jetunc_filename_pfL1FastJetL2L3);
}

//####################################################################################
void CoreUtil::jec::setJECFileNames(std::string jecEra, std::string type)
{
    std::string basepath = jecdatapath.Data(); // A bit messy? Oh what the hell, this whole thing is messy... Besides I like TString generally
    jetcorr_filenames_pfL1FastJetL2L3.clear();
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L1FastJet_AK4PFchs.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L3Absolute_AK4PFchs.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L2L3Residual_AK4PFchs.txt");
    jetunc_filename_pfL1FastJetL2L3 = basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_Uncertainty_AK4PFchs.txt";
}

void CoreUtil::jec::setFatJECFor(TString filename, bool isfastsim)
{
    // If current filename same as the one already set to then skip
    if (current_filename == filename)
        return;
    
    // Reset
    resetObjects();
    
    //    cout << endl;
    cout << "CoreUtil:: Create new JEC for " << filename << endl;
    
    // Set the current_filename to the provided filename
    current_filename = filename;
    
    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2015/2016
    //-----------------------------------------------------------------------------------------------------------------------------------
    //JECs for 76X
    if (filename.Contains("16Dec2015") || filename.Contains("76X_mcRun2"))
    {
        if (filename.Contains("Run2015C") || filename.Contains("Run2015D"))
        {
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_76X/DATA/Fall15_25nsV2_DATA_Uncertainty_AK4PFchs.txt");
        }
        else
        {
            // files for 76X MC
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_76X/MC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt");
        }
    }
    else if (filename.Contains("80MiniAODv") || filename.Contains("RelVal"))
    {
        // files for 80X MC, ICHEP production
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/MC/Spring16_25nsV1_MC_Uncertainty_AK4PFchs.txt");
    }
    else if (filename.Contains("Summer16") || filename.Contains("TEST") || filename.Contains("PrivateWWW") || filename.Contains("TChiWH_HToVVTauTau_HToBB_mChargino200_mLSP50") || filename.Contains("WWW_BW") || filename.Contains("GluGlu"))
    {
        // files for 80X MC, Summer16 (Moriond17) production
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/MC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt");
    }
    else if (filename.Contains("Run2016") || filename.Contains("CMSSW_8_0_11_V08-00-06"))
    {
        // // This scheme will eventually be used for Moriond. For now just use JECs from ICHEP
        if (filename.Contains("Run2016B") ||
            filename.Contains("Run2016C") ||
            filename.Contains("Run2016D"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016BCDV3_DATA_Uncertainty_AK4PFchs.txt");
        }
        if (filename.Contains("Run2016E") ||
            filename.Contains("Run2016F"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016EFV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
        if (filename.Contains("Run2016G"))
        {
            // files for 80X Data
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
        if (filename.Contains("Run2016H"))
        {
            jetcorr_filenames_pfL1.clear();
            jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.clear();
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt");
            jetcorr_filenames_pfL1_postrun278802.clear();
            jetcorr_filenames_pfL1_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.clear();
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt");
            jetcorr_filenames_pfL1FastJetL2L3_postrun278802.push_back("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt");
            jecUnc_postrun278802 = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/DATA/Summer16_23Sep2016HV3_DATA_Uncertainty_AK4PFchs.txt");
            jet_corrector_pfL1FastJetL2L3_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3_postrun278802);
        }
    }
    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2017
    //-----------------------------------------------------------------------------------------------------------------------------------
    else if (filename.Contains("Run2017"))
    {
        if (filename.Contains("Run2017B"))
        {
            setFatJECFileNames("Fall17_17Nov2017B_V6", "DATA");
        }
        if (filename.Contains("Run2017C"))
        {
            setFatJECFileNames("Fall17_17Nov2017C_V6", "DATA");
        }
        if (filename.Contains("Run2017D"))
        {
            setFatJECFileNames("Fall17_17Nov2017D_V6", "DATA");
        }
        if (filename.Contains("Run2017E"))
        {
            setFatJECFileNames("Fall17_17Nov2017E_V6", "DATA");
        }
        if (filename.Contains("Run2017F"))
        {
            setFatJECFileNames("Fall17_17Nov2017F_V6", "DATA");
        }
    }
    else if (filename.Contains("Fall17MiniAOD"))
    {
        setFatJECFileNames("Fall17_17Nov2017_V6", "MC");
    }
    else if (filename.Contains("PhotophobicALP"))
    {
        setFatJECFileNames("Fall17_17Nov2017_V6", "MC");
    }
    
    //-----------------------------------------------------------------------------------------------------------------------------------
    // 2016 fast sim
    //-----------------------------------------------------------------------------------------------------------------------------------
    if (isfastsim)
    {
        // files for 25ns fastsim samples
        jetcorr_filenames_pfL1.clear();
        jetcorr_filenames_pfL1.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.clear();
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L2Relative_AK4PFchs.txt");
        jetcorr_filenames_pfL1FastJetL2L3.push_back("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_L3Absolute_AK4PFchs.txt");
        jecUnc = new JetCorrectionUncertainty("coreutil/data/jetCorrections/source_80X/FASTSIM/Spring16_FastSimV1_Uncertainty_AK4PFchs.txt");
    }
    
    
    //===================================================================================================================================
    // Now create the jet corrector based on what it is set to
    //===================================================================================================================================
    
    if (jetcorr_filenames_pfL1FastJetL2L3.size() == 0)
    {
        cout << "CoreUtil:: Error, sample not found. Check the JECs." << endl;
        exit(100);
    }
    
    // Print what is used
    cout << "CoreUtil:: JECs used:" << endl;
    for (size_t jecind = 0; jecind < jetcorr_filenames_pfL1FastJetL2L3.size(); jecind++)
    {
        cout << "CoreUtil:: " << jetcorr_filenames_pfL1FastJetL2L3.at(jecind) << endl;
    }
    
    // Now create
    jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);
    jet_corrector_pfL1  = makeJetCorrector(jetcorr_filenames_pfL1);
    jet_corrector_pfL1_postrun278802  = makeJetCorrector(jetcorr_filenames_pfL1_postrun278802);
    // For jet uncertainty for 2016 the code was written in bad ways such that the jetunc was created on the fly above in the function
    // Starting 2017 I tried to clean up a little bit so that it is created at the end
    if (!jecUnc)
        jecUnc = new JetCorrectionUncertainty(jetunc_filename_pfL1FastJetL2L3);
}

//####################################################################################
void CoreUtil::jec::setFatJECFileNames(std::string jecEra, std::string type)
{
    std::string basepath = jecdatapath.Data(); // A bit messy? Oh what the hell, this whole thing is messy... Besides I like TString generally
    jetcorr_filenames_pfL1FastJetL2L3.clear();
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L1FastJet_AK8PFPuppi.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L2Relative_AK8PFPuppi.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L3Absolute_AK8PFPuppi.txt");
    jetcorr_filenames_pfL1FastJetL2L3.push_back(basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_L2L3Residual_AK8PFPuppi.txt");
    jetunc_filename_pfL1FastJetL2L3 = basepath+"run2_25ns/"+jecEra+"_"+type+"/"+jecEra+"_"+type+"_Uncertainty_AK8PFPuppi.txt";
}
