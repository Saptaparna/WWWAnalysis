#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"
#include "TFile.h"
#include "THStack.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include <fstream>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

void makeDatacard_noST()
{
  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"}; 

  int size = sizeof(coeff)/sizeof(TString);

  std::cout << "size = " << size << std::endl;

  double yield1 = 10.5783;
  double yield2 = 10.0037;
  double yield3 = 9.45216;
  double yield4 = 8.92382;
  double yield5 = 8.41863;
  double yield6 = 7.9366;
  double yield7 = 7.47773;
  double yield8 = 7.04201;
  double yield9 = 6.62945;
  double yield10 = 6.24005;
  double yield11 = 5.87381;
  double yield12 = 5.53072;
  double yield13 = 5.21079;
  double yield14 = 4.91402;
  double yield15 = 4.6404;
  double yield16 = 4.38994;
  double yield17 = 4.16264;
  double yield18 = 3.9585;
  double yield19 = 3.77751;
  double yield20 = 3.61968;
  double yield21 = 3.48501;
  double yield22 = 3.37349;
  double yield23 = 3.28513;
  double yield24 = 3.21993;
  double yield25 = 3.17789;
  double yield26 = 3.159;
  double yield27 = 3.16327;
  double yield28 = 3.19069;
  double yield29 = 3.24128;
  double yield30 = 3.31502;
  double yield31 = 3.41192;
  double yield32 = 3.53197;
  double yield33 = 3.67518;
  double yield34 = 3.84155;
  double yield35 = 4.03108;
  double yield36 = 4.24376;
  double yield37 = 4.4796;
  double yield38 = 4.7386;
  double yield39 = 5.02076;
  double yield40 = 5.32607;
  double yield41 = 5.65454;
  double yield42 = 6.00616;
  double yield43 = 6.38095;
  double yield44 = 6.77889;
  double yield45 = 7.19998;
  double yield46 = 7.64424;
  double yield47 = 8.11165;
  double yield48 = 8.60222;
  double yield49 = 9.11594;
  double yield50 = 9.65283;
  double yield51 = 10.2129;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;
  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6, 32.6};

  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_noST_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}

/*
void makeDatacard_ST250()
{

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {8.52, 5.97, 4.0, 2.6, 1.77, 1.52, 1.84, 2.73, 4.19, 6.23, 6.12};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9, 10.9};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST250_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}

void makeDatacard_ST500()
{

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {7.35, 4.82, 2.87, 1.48, 0.652, 0.395, 0.703, 1.58, 3.02, 5.03, 8.21};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {4.36, 4.36, 4.36, 4.36, 4.36, 4.36, 4.36, 4.36, 4.36, 4.36, 4.36};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST500_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}
*/
void makeDatacard_ST750()
{

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {6.66, 4.29, 2.45, 1.14, 0.363, 0.118, 0.405, 1.22, 2.58, 6.7, 4.45};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST750_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}

void makeDatacard_ST1000()
{

  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};  

  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 6.81848;
  double yield2 = 6.29911;
  double yield3 = 5.80078;
  double yield4 = 5.32351;
  double yield5 = 4.86727;
  double yield6 = 4.43208;
  double yield7 = 4.01794;
  double yield8 = 3.62483;
  double yield9 = 3.25278;
  double yield10 = 2.90177;
  double yield11 = 2.5718;
  double yield12 = 2.26288;
  double yield13 = 1.975;
  double yield14 = 1.70816;
  double yield15 = 1.46237;
  double yield16 = 1.23763;
  double yield17 = 1.03393;
  double yield18 = 0.851271;
  double yield19 = 0.689659;
  double yield20 = 0.549092;
  double yield21 = 0.429569;
  double yield22 = 0.331091;
  double yield23 = 0.253657;
  double yield24 = 0.197268;
  double yield25 = 0.161923;
  double yield26 = 0.147622;
  double yield27 = 0.154366;
  double yield28 = 0.182155;
  double yield29 = 0.230988;
  double yield30 = 0.300865;
  double yield31 = 0.391787;
  double yield32 = 0.503754;
  double yield33 = 0.636765;
  double yield34 = 0.79082;
  double yield35 = 0.96592;
  double yield36 = 1.16206;
  double yield37 = 1.37925;
  double yield38 = 1.61749;
  double yield39 = 1.87676;
  double yield40 = 2.15709;
  double yield41 = 2.45845;
  double yield42 = 2.78086;
  double yield43 = 3.12432;
  double yield44 = 3.48882;
  double yield45 = 3.87436;
  double yield46 = 4.28095;
  double yield47 = 4.70859;
  double yield48 = 5.15727;
  double yield49 = 5.62699;
  double yield50 = 6.11775;
  double yield51 = 6.62957 ;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66, 3.66};

  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST1000_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}


void makeDatacard_ST1500()
{
  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"}; 
  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 5.67921;
  double yield2 = 5.23835;
  double yield3 = 4.81539;
  double yield4 = 4.41034;
  double yield5 = 4.0232;
  double yield6 = 3.65396;
  double yield7 = 3.30263;
  double yield8 = 2.9692;
  double yield9 = 2.65368;
  double yield10 = 2.35607;
  double yield11 = 2.07636;
  double yield12 = 1.81455;
  double yield13 = 1.57065;
  double yield14 = 1.34466;
  double yield15 = 1.13658;
  double yield16 = 0.946397;
  double yield17 = 0.774123;
  double yield18 = 0.619754;
  double yield19 = 0.483292;
  double yield20 = 0.364735;
  double yield21 = 0.264084;
  double yield22 = 0.181339;
  double yield23 = 0.1165;
  double yield24 = 0.0695669;
  double yield25 = 0.0405394;
  double yield26 = 0.0294178;
  double yield27 = 0.036202;
  double yield28 = 0.060892;
  double yield29 = 0.103488;
  double yield30 = 0.16399;
  double yield31 = 0.242397;
  double yield32 = 0.338711;
  double yield33 = 0.45293;
  double yield34 = 0.585055;
  double yield35 = 0.735086;
  double yield36 = 0.903023;
  double yield37 = 1.08887;
  double yield38 = 1.29261;
  double yield39 = 1.51427;
  double yield40 = 1.75383;
  double yield41 = 2.01129;
  double yield42 = 2.28667;
  double yield43 = 2.57994;
  double yield44 = 2.89113;
  double yield45 = 3.22022;
  double yield46 = 3.56721;
  double yield47 = 3.93211;
  double yield48 = 4.31492;
  double yield49 = 4.71563;
  double yield50 = 5.13425;
  double yield51 = 5.57078;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04};

  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST1500_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl; 
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;  
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;  
  }
}

void makeDatacard_ST2000()
{

  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"}; 

  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 4.48827*0.75;
  double yield2 = 4.13778*0.75;
  double yield3 = 3.80155*0.75;
  double yield4 = 3.47958*0.75;
  double yield5 = 3.17188*0.75;
  double yield6 = 2.87843*0.75;
  double yield7 = 2.59925*0.75;
  double yield8 = 2.33433*0.75;
  double yield9 = 2.08366*0.75;
  double yield10 = 1.84726*0.75;
  double yield11 = 1.62512*0.75;
  double yield12 = 1.41725*0.75;
  double yield13 = 1.22363*0.75;
  double yield14 = 1.04427*0.75;
  double yield15 = 0.879178*0.75;
  double yield16 = 0.728344*0.75;
  double yield17 = 0.591771*0.75;
  double yield18 = 0.469459*0.75;
  double yield19 = 0.361408*0.75;
  double yield20 = 0.267619*0.75;
  double yield21 = 0.18809*0.75;
  double yield22 = 0.122823*0.75;
  double yield23 = 0.0718162*0.75;
  double yield24 = 0.0350709*0.75;
  double yield25 = 0.0125866*0.75;
  double yield26 = 0.00436339*0.75;
  double yield27 = 0.0104013*0.75;
  double yield28 = 0.0307002*0.75;
  double yield29 = 0.0652603*0.75;
  double yield30 = 0.114081*0.75;
  double yield31 = 0.177164*0.75;
  double yield32 = 0.254507*0.75;
  double yield33 = 0.346111*0.75;
  double yield34 = 0.451977*0.75;
  double yield35 = 0.572103*0.75;
  double yield36 = 0.706491*0.75;
  double yield37 = 0.85514*0.75;
  double yield38 = 1.01805*0.75;
  double yield39 = 1.19522*0.75;
  double yield40 = 1.38665*0.75;
  double yield41 = 1.59235*0.75;
  double yield42 = 1.8123*0.75;
  double yield43 = 2.04651*0.75;
  double yield44 = 2.29499*0.75;
  double yield45 = 2.55773*0.75;
  double yield46 = 2.83473*0.75;
  double yield47 = 3.12599*0.75;
  double yield48 = 3.43151*0.75;
  double yield49 = 3.75129*0.75;
  double yield50 = 4.08533*0.75;
  double yield51 = 4.43363*0.75;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242, 0.242};

  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST2000_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi      lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_bkg    lnN    1.00      1.30    background uncertainty 30%" <<std::endl;
  }
}

void makeDatacard_ST2500()
{
  
  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
  
  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 1.99993;
  double yield2 = 1.84337;
  double yield3 = 1.69319;
  double yield4 = 1.5494;
  double yield5 = 1.41198;
  double yield6 = 1.28096;
  double yield7 = 1.15632;
  double yield8 = 1.03806;
  double yield9 = 0.926183;
  double yield10 = 0.820694;
  double yield11 = 0.721589;
  double yield12 = 0.628868;
  double yield13 = 0.542532;
  double yield14 = 0.46258;
  double yield15 = 0.389012;
  double yield16 = 0.321829;
  double yield17 = 0.26103;
  double yield18 = 0.206615;
  double yield19 = 0.158585;
  double yield20 = 0.116939;
  double yield21 = 0.0816769;
  double yield22 = 0.0527995;
  double yield23 = 0.0303065;
  double yield24 = 0.0141979;
  double yield25 = 0.00447351;
  double yield26 = 0.0011335;
  double yield27 = 0.00417783;
  double yield28 = 0.0136065;
  double yield29 = 0.0294195;
  double yield30 = 0.0516168;
  double yield31 = 0.0801985;
  double yield32 = 0.115165;
  double yield33 = 0.156515;
  double yield34 = 0.20425;
  double yield35 = 0.258369;
  double yield36 = 0.318872;
  double yield37 = 0.38576;
  double yield38 = 0.459032;
  double yield39 = 0.538688;
  double yield40 = 0.624729;
  double yield41 = 0.717154;
  double yield42 = 0.815963;
  double yield43 = 0.921157;
  double yield44 = 1.03273;
  double yield45 = 1.1507;
  double yield46 = 1.27504;
  double yield47 = 1.40577;
  double yield48 = 1.54289;
  double yield49 = 1.68639;
  double yield50 = 1.83627;
  double yield51 = 1.99254;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504, 0.0504};

  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];

  for(int i=0; i<size; i++)
  {
    std::ofstream outfile("realistic-counting-experiment_ST2500_ft_"+coeff[i]+".txt");
    outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
    outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
    outfile<<"imax 1  number of channels" <<std::endl;
    outfile<<"jmax *  number of backgrounds" <<std::endl;
    outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
    outfile<<"bin 1" <<std::endl;
    outfile<<"observation 0" <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
    outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
    outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
    outfile<<"# on each process and bin" <<std::endl;
    outfile<<"bin         1                1" <<std::endl;
    outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
    outfile<<"process     0                1" <<std::endl;
    outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
    outfile<<"------------" <<std::endl;
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}

void makeDatacard_ST2000_MuMu()
{
    
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    
    int size = sizeof(coeff)/sizeof(TString);
    
    double yield1 = 1.90079;
    double yield2 = 1.75241;
    double yield3 = 1.61006;
    double yield4 = 1.47376;
    double yield5 = 1.34349;
    double yield6 = 1.21926;
    double yield7 = 1.10108;
    double yield8 = 0.988924;
    double yield9 = 0.882812;
    double yield10 = 0.782738;
    double yield11 = 0.688703;
    double yield12 = 0.600706;
    double yield13 = 0.518747;
    double yield14 = 0.442827;
    double yield15 = 0.372945;
    double yield16 = 0.309102;
    double yield17 = 0.251297;
    double yield18 = 0.199531;
    double yield19 = 0.153803;
    double yield20 = 0.114114;
    double yield21 = 0.0804631;
    double yield22 = 0.0528507;
    double yield23 = 0.0312767;
    double yield24 = 0.0157412;
    double yield25 = 0.00624417;
    double yield26 = 0.00278558;
    double yield27 = 0.00536544;
    double yield28 = 0.0139838;
    double yield29 = 0.0286405;
    double yield30 = 0.0493358;
    double yield31 = 0.0760695;
    double yield32 = 0.108842;
    double yield33 = 0.147652;
    double yield34 = 0.192501;
    double yield35 = 0.243389;
    double yield36 = 0.300315;
    double yield37 = 0.363279;
    double yield38 = 0.432282;
    double yield39 = 0.507324;
    double yield40 = 0.588403;
    double yield41 = 0.675522;
    double yield42 = 0.768678;
    double yield43 = 0.867874;
    double yield44 = 0.973107;
    double yield45 = 1.08438;
    double yield46 = 1.20169;
    double yield47 = 1.32504;
    double yield48 = 1.45443;
    double yield49 = 1.58985;
    double yield50 = 1.73132;
    double yield51 = 1.87882;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0299};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST2000_MuMu_ft_"+coeff[i]+".txt");
        outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
        outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
        outfile<<"imax 1  number of channels" <<std::endl;
        outfile<<"jmax *  number of backgrounds" <<std::endl;
        outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
        outfile<<"bin 1" <<std::endl;
        outfile<<"observation 0" <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
        outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
        outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
        outfile<<"# on each process and bin" <<std::endl;
        outfile<<"bin         1                1" <<std::endl;
        outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
        outfile<<"process     0                1" <<std::endl;
        outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
        outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
    }
}

void makeDatacard_ST2000_ElMu()
{
    
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    
    int size = sizeof(coeff)/sizeof(TString);
    
    double yield1 = 2.0596;
    double yield2 = 1.89866;
    double yield3 = 1.74428;
    double yield4 = 1.59645;
    double yield5 = 1.45516;
    double yield6 = 1.32043;
    double yield7 = 1.19226;
    double yield8 = 1.07063;
    double yield9 = 0.955552;
    double yield10 = 0.847028;
    double yield11 = 0.745055;
    double yield12 = 0.649633;
    double yield13 = 0.560763;
    double yield14 = 0.478444;
    double yield15 = 0.402677;
    double yield16 = 0.333461;
    double yield17 = 0.270797;
    double yield18 = 0.214684;
    double yield19 = 0.165122;
    double yield20 = 0.122112;
    double yield21 = 0.0856533;
    double yield22 = 0.055746;
    double yield23 = 0.03239;
    double yield24 = 0.0155855;
    double yield25 = 0.00533247;
    double yield26 = 0.00163082;
    double yield27 = 0.00448058;
    double yield28 = 0.0138818;
    double yield29 = 0.0298344;
    double yield30 = 0.0523384;
    double yield31 = 0.0813938;
    double yield32 = 0.117001;
    double yield33 = 0.159159;
    double yield34 = 0.207869;
    double yield35 = 0.26313;
    double yield36 = 0.324942;
    double yield37 = 0.393306;
    double yield38 = 0.468222;
    double yield39 = 0.549688;
    double yield40 = 0.637707;
    double yield41 = 0.732276;
    double yield42 = 0.833397;
    double yield43 = 0.94107;
    double yield44 = 1.05529;
    double yield45 = 1.17607;
    double yield46 = 1.3034;
    double yield47 = 1.43727;
    double yield48 = 1.5777;
    double yield49 = 1.72468;
    double yield50 = 1.87822;
    double yield51 = 2.0383;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908, 0.0908};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST2000_ElMu_ft_"+coeff[i]+".txt");
        outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
        outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
        outfile<<"imax 1  number of channels" <<std::endl;
        outfile<<"jmax *  number of backgrounds" <<std::endl;
        outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
        outfile<<"bin 1" <<std::endl;
        outfile<<"observation 0" <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
        outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
        outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
        outfile<<"# on each process and bin" <<std::endl;
        outfile<<"bin         1                1" <<std::endl;
        outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
        outfile<<"process     0                1" <<std::endl;
        outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
        outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
    }
}

void makeDatacard_ST2000_ElEl()
{
    
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    
    int size = sizeof(coeff)/sizeof(TString);
    
    double yield1 = 0.523985;
    double yield2 = 0.483107;
    double yield3 = 0.443889;
    double yield4 = 0.406333;
    double yield5 = 0.370437;
    double yield6 = 0.336203;
    double yield7 = 0.303629;
    double yield8 = 0.272717;
    double yield9 = 0.243466;
    double yield10 = 0.215875;
    double yield11 = 0.189946;
    double yield12 = 0.165678;
    double yield13 = 0.143071;
    double yield14 = 0.122125;
    double yield15 = 0.10284;
    double yield16 = 0.0852158;
    double yield17 = 0.0692529;
    double yield18 = 0.0549511;
    double yield19 = 0.0423102;
    double yield20 = 0.0313305;
    double yield21 = 0.0220117;
    double yield22 = 0.0143541;
    double yield23 = 0.00835743;
    double yield24 = 0.00402185;
    double yield25 = 0.00134731;
    double yield26 = 0.000333823;
    double yield27 = 0.000981383;
    double yield28 = 0.00328999;
    double yield29 = 0.00725965;
    double yield30 = 0.0128904;
    double yield31 = 0.0201821;
    double yield32 = 0.0291349;
    double yield33 = 0.0397488;
    double yield34 = 0.0520236;
    double yield35 = 0.0659596;
    double yield36 = 0.0815566;
    double yield37 = 0.0988146;
    double yield38 = 0.117734;
    double yield39 = 0.138314;
    double yield40 = 0.160555;
    double yield41 = 0.184457;
    double yield42 = 0.210021;
    double yield43 = 0.237245;
    double yield44 = 0.26613;
    double yield45 = 0.296677;
    double yield46 = 0.328884;
    double yield47 = 0.362753;
    double yield48 = 0.398282;
    double yield49 = 0.435473;
    double yield50 = 0.474324;
    double yield51 = 0.514837;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388, 0.0388};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST2000_ElEl_ft_"+coeff[i]+".txt");
        outfile<<"# Simple counting experiment, with one signal and a few background processes" <<std::endl;
        outfile<<"# Simplified version of the 35/pb H->WW analysis for mH = 160 GeV" <<std::endl;
        outfile<<"imax 1  number of channels" <<std::endl;
        outfile<<"jmax *  number of backgrounds" <<std::endl;
        outfile<<"kmax *  number of nuisance parameters (sources of systematical uncertainties)"  <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# we have just one channel, in which we observe 0 events" <<std::endl;
        outfile<<"bin 1" <<std::endl;
        outfile<<"observation 0" <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"# now we list the expected events for signal and all backgrounds in that bin" <<std::endl;
        outfile<<"# the second 'process' line must have a positive number for backgrounds, and 0 for signal" <<std::endl;
        outfile<<"# then we list the independent sources of uncertainties, and give their effect (syst. error)" <<std::endl;
        outfile<<"# on each process and bin" <<std::endl;
        outfile<<"bin         1                1" <<std::endl;
        outfile<<"process     ft_"+coeff[i]+"        bkg" <<std::endl;
        outfile<<"process     0                1" <<std::endl;
        outfile<<"rate        "+std::to_string(signal[i])+"         "+std::to_string(bkg[i]) <<std::endl;
        outfile<<"------------" <<std::endl;
        outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
        outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
    }
}

