#ifndef jec_h
#define jec_h

#include <vector>
#include <string>
#include <iostream>

#include "TString.h"
#include "TSystem.h"

#include "Tools/jetcorr/FactorizedJetCorrector.h"
#include "Tools/JetCorrector.h"
#include "Tools/jetcorr/JetCorrectionUncertainty.h"

namespace CoreUtil
{
    class jec
    {
        public:
            // PU correction
            std::vector<std::string> jetcorr_filenames_pfL1;
            FactorizedJetCorrector   * jet_corrector_pfL1;
            JetCorrectionUncertainty * jecUnc_pfL1;

            // PU corrections for later runs in 2016F
            std::vector<std::string> jetcorr_filenames_pfL1_postrun278802;
            FactorizedJetCorrector   * jet_corrector_pfL1_postrun278802;
            JetCorrectionUncertainty * jecUnc_pfL1_postrun278802;

            // Nominal correction
            std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
            std::string jetunc_filename_pfL1FastJetL2L3;
            FactorizedJetCorrector   * jet_corrector_pfL1FastJetL2L3;
            JetCorrectionUncertainty * jecUnc;

            // corrections for later runs in 2016F
            std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3_postrun278802;
            FactorizedJetCorrector   * jet_corrector_pfL1FastJetL2L3_postrun278802;
            JetCorrectionUncertainty * jecUnc_postrun278802;

            TString current_filename;
            TString jecdatapath;

            jec();
            ~jec();
            void resetObjects();
            FactorizedJetCorrector* getJECL1(bool = false, int = -999);
            FactorizedJetCorrector* getJEC(bool = false, int = -999);
            JetCorrectionUncertainty* getJECUnc(bool = false, int = -999);
            void setJECFor(TString, bool = false);
            void setJECFileNames(std::string, std::string);
            void setFatJECFor(TString, bool = false);
            void setFatJECFileNames(std::string, std::string);
    };
}

#endif
