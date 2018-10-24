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

  double yield1 = 9.30881;
  double yield2 = 8.74075;
  double yield3 = 8.1956;
  double yield4 = 7.67338;
  double yield5 = 7.17408;
  double yield6 = 6.69771;
  double yield7 = 6.24425;
  double yield8 = 5.81372;
  double yield9 = 5.40611;
  double yield10 = 5.02142;
  double yield11 = 4.65966;
  double yield12 = 4.32081;
  double yield13 = 4.00489;
  double yield14 = 3.71189;
  double yield15 = 3.44181;
  double yield16 = 3.19466;
  double yield17 = 2.97042;
  double yield18 = 2.76911;
  double yield19 = 2.59072;
  double yield20 = 2.43526;
  double yield21 = 2.30271;
  double yield22 = 2.19309;
  double yield23 = 2.10639;
  double yield24 = 2.04261;
  double yield25 = 2.00176;
  double yield26 = 1.98382;
  double yield27 = 1.98881;
  double yield28 = 2.01672;
  double yield29 = 2.06755;
  double yield30 = 2.14131;
  double yield31 = 2.23799;
  double yield32 = 2.35759;
  double yield33 = 2.50011;
  double yield34 = 2.66555;
  double yield35 = 2.85392;
  double yield36 = 3.0652;
  double yield37 = 3.29941;
  double yield38 = 3.55654;
  double yield39 = 3.8366;
  double yield40 = 4.13957;
  double yield41 = 4.46547;
  double yield42 = 4.81429;
  double yield43 = 5.18604;
  double yield44 = 5.5807;
  double yield45 = 5.99829;
  double yield46 = 6.4388;
  double yield47 = 6.90223;
  double yield48 = 7.38858;
  double yield49 = 7.89786;
  double yield50 = 8.43006;
  double yield51 = 8.98517;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;
  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2};

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

  double yield1 = 6.12406;
  double yield2 = 5.6499;
  double yield3 = 5.195;
  double yield4 = 4.75934;
  double yield5 = 4.34292;
  double yield6 = 3.94576;
  double yield7 = 3.56785;
  double yield8 = 3.20918;
  double yield9 = 2.86977;
  double yield10 = 2.5496;
  double yield11 = 2.24868;
  double yield12 = 1.96701;
  double yield13 = 1.70459;
  double yield14 = 1.46142;
  double yield15 = 1.2375;
  double yield16 = 1.03282;
  double yield17 = 0.847396;
  double yield18 = 0.681219;
  double yield19 = 0.534292;
  double yield20 = 0.406613;
  double yield21 = 0.298183;
  double yield22 = 0.209001;
  double yield23 = 0.139069;
  double yield24 = 0.0883852;
  double yield25 = 0.0569503;
  double yield26 = 0.0447644;
  double yield27 = 0.0518272;
  double yield28 = 0.0781389;
  double yield29 = 0.123699;
  double yield30 = 0.188509;
  double yield31 = 0.272567;
  double yield32 = 0.375874;
  double yield33 = 0.49843;
  double yield34 = 0.640234;
  double yield35 = 0.801288;
  double yield36 = 0.98159;
  double yield37 = 1.18114;
  double yield38 = 1.39994;
  double yield39 = 1.63799;
  double yield40 = 1.89529;
  double yield41 = 2.17183;
  double yield42 = 2.46763;
  double yield43 = 2.78267;
  double yield44 = 3.11697;
  double yield45 = 3.47051;
  double yield46 = 3.8433;
  double yield47 = 4.23534;
  double yield48 = 4.64663;
  double yield49 = 5.07716;
  double yield50 = 5.52695;
  double yield51 = 5.99598;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258};

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

  double yield1 = 4.4906*0.75;
  double yield2 = 4.14008*0.75;
  double yield3 = 3.80383*0.75;
  double yield4 = 3.48185*0.75;
  double yield5 = 3.17416*0.75;
  double yield6 = 2.88073*0.75;
  double yield7 = 2.60158*0.75;
  double yield8 = 2.3367*0.75;
  double yield9 = 2.08609*0.75;
  double yield10 = 1.84976*0.75;
  double yield11 = 1.6277*0.75;
  double yield12 = 1.41992*0.75;
  double yield13 = 1.22641*0.75;
  double yield14 = 1.04717*0.75;
  double yield15 = 0.88221*0.75;
  double yield16 = 0.731521*0.75;
  double yield17 = 0.595105*0.75;
  double yield18 = 0.472963*0.75;
  double yield19 = 0.365095*0.75;
  double yield20 = 0.271501*0.75;
  double yield21 = 0.19218*0.75;
  double yield22 = 0.127133*0.75;
  double yield23 = 0.0763602*0.75;
  double yield24 = 0.0398606*0.75;
  double yield25 = 0.0176348*0.75;
  double yield26 = 0.00968266*0.75;
  double yield27 = 0.0160042*0.75;
  double yield28 = 0.0365995*0.75;
  double yield29 = 0.0714685*0.75;
  double yield30 = 0.120611*0.75;
  double yield31 = 0.184028*0.75;
  double yield32 = 0.261718*0.75;
  double yield33 = 0.353681*0.75;
  double yield34 = 0.459919*0.75;
  double yield35 = 0.58043*0.75;
  double yield36 = 0.715215*0.75;
  double yield37 = 0.864274*0.75;
  double yield38 = 1.02761*0.75;
  double yield39 = 1.20521*0.75;
  double yield40 = 1.39709*0.75;
  double yield41 = 1.60325*0.75;
  double yield42 = 1.82367*0.75;
  double yield43 = 2.05837*0.75;
  double yield44 = 2.30735*0.75;
  double yield45 = 2.5706*0.75;
  double yield46 = 2.84812*0.75;
  double yield47 = 3.13991*0.75;
  double yield48 = 3.44598*0.75;
  double yield49 = 3.76633*0.75;
  double yield50 = 4.10094*0.75;
  double yield51 = 4.44983*0.75;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321};

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
    outfile<<"lumi      lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_bkg    lnN    1.00      1.30    background uncertainty 30%" <<std::endl; 
     
  }
}

void makeDatacard_ST2000()
{

  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"}; 

  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 2.88123;
  double yield2 = 2.656;
  double yield3 = 2.43994;
  double yield4 = 2.23306;
  double yield5 = 2.03536;
  double yield6 = 1.84684;
  double yield7 = 1.66749;
  double yield8 = 1.49731;
  double yield9 = 1.33631;
  double yield10 = 1.18449;
  double yield11 = 1.04185;
  double yield12 = 0.908376;
  double yield13 = 0.784081;
  double yield14 = 0.668963;
  double yield15 = 0.563021;
  double yield16 = 0.466255;
  double yield17 = 0.378665;
  double yield18 = 0.30025;
  double yield19 = 0.231012;
  double yield20 = 0.17095;
  double yield21 = 0.120063;
  double yield22 = 0.0783531;
  double yield23 = 0.0458187;
  double yield24 = 0.0224604;
  double yield25 = 0.00827796;
  double yield26 = 0.00327155;
  double yield27 = 0.00744112;
  double yield28 = 0.0207867;
  double yield29 = 0.0433082;
  double yield30 = 0.0750057;
  double yield31 = 0.115879;
  double yield32 = 0.165929;
  double yield33 = 0.225154;
  double yield34 = 0.293556;
  double yield35 = 0.371133;
  double yield36 = 0.457886;
  double yield37 = 0.553816;
  double yield38 = 0.658921;
  double yield39 = 0.773202;
  double yield40 = 0.89666;
  double yield41 = 1.02929;
  double yield42 = 1.1711;
  double yield43 = 1.32209;
  double yield44 = 1.48225;
  double yield45 = 1.65159;
  double yield46 = 1.8301;
  double yield47 = 2.01779;
  double yield48 = 2.21465;
  double yield49 = 2.42069;
  double yield50 = 2.63591;
  double yield51 = 2.8603;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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
    outfile<<"lumi    lnN    1.026     1.026   lumi affects both signal and all backgrounds. lnN = lognormal" <<std::endl;
    outfile<<"xs_SS   lnN    -         1.30    background uncertainty 30%" <<std::endl;
  }
}

void makeDatacard_ST2500()
{
  
  TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
  
  int size = sizeof(coeff)/sizeof(TString);

  double yield1 = 1.27015;
  double yield2 = 1.17081;
  double yield3 = 1.07551;
  double yield4 = 0.984255;
  double yield5 = 0.89705;
  double yield6 = 0.813892;
  double yield7 = 0.73478;
  double yield8 = 0.659715;
  double yield9 = 0.588697;
  double yield10 = 0.521725;
  double yield11 = 0.458801;
  double yield12 = 0.399923;
  double yield13 = 0.345092;
  double yield14 = 0.294307;
  double yield15 = 0.24757;
  double yield16 = 0.204879;
  double yield17 = 0.166235;
  double yield18 = 0.131637;
  double yield19 = 0.101087;
  double yield20 = 0.0745827;
  double yield21 = 0.0521256;
  double yield22 = 0.0337152;
  double yield23 = 0.0193516;
  double yield24 = 0.00903479;
  double yield25 = 0.0027647;
  double yield26 = 0.000541373;
  double yield27 = 0.0023648;
  double yield28 = 0.00823498;
  double yield29 = 0.0181519;
  double yield30 = 0.0321156;
  double yield31 = 0.0501261;
  double yield32 = 0.0721833;
  double yield33 = 0.0982872;
  double yield34 = 0.128438;
  double yield35 = 0.162635;
  double yield36 = 0.20088;
  double yield37 = 0.243171;
  double yield38 = 0.289508;
  double yield39 = 0.339893;
  double yield40 = 0.394324;
  double yield41 = 0.452802;
  double yield42 = 0.515327;
  double yield43 = 0.581898;
  double yield44 = 0.652517;
  double yield45 = 0.727182;
  double yield46 = 0.805893;
  double yield47 = 0.888652;
  double yield48 = 0.975457;
  double yield49 = 1.06631;
  double yield50 = 1.16121;
  double yield51 = 1.26015;

  double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};

  double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];

  signal[25] = 0.000001;

  double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

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

void makeDatacard_ST1500_0OSSF()
{
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    int size = sizeof(coeff)/sizeof(TString);
    
    double yield1 = 1.15956;
    double yield2 = 1.06904;
    double yield3 = 0.982218;
    double yield4 = 0.899082;
    double yield5 = 0.819633;
    double yield6 = 0.743873;
    double yield7 = 0.6718;
    double yield8 = 0.603415;
    double yield9 = 0.538719;
    double yield10 = 0.47771;
    double yield11 = 0.420389;
    double yield12 = 0.366757;
    double yield13 = 0.316812;
    double yield14 = 0.270555;
    double yield15 = 0.227987;
    double yield16 = 0.189106;
    double yield17 = 0.153913;
    double yield18 = 0.122409;
    double yield19 = 0.0945919;
    double yield20 = 0.0704631;
    double yield21 = 0.0500224;
    double yield22 = 0.0332696;
    double yield23 = 0.0202049;
    double yield24 = 0.0108281;
    double yield25 = 0.00513928;
    double yield26 = 0.00313848;
    double yield27 = 0.00482567;
    double yield28 = 0.0102008;
    double yield29 = 0.019264;
    double yield30 = 0.0320152;
    double yield31 = 0.0484543;
    double yield32 = 0.0685815;
    double yield33 = 0.0923966;
    double yield34 = 0.1199;
    double yield35 = 0.151091;
    double yield36 = 0.18597;
    double yield37 = 0.224537;
    double yield38 = 0.266792;
    double yield39 = 0.312735;
    double yield40 = 0.362366;
    double yield41 = 0.415685;
    double yield42 = 0.472692;
    double yield43 = 0.533387;
    double yield44 = 0.59777;
    double yield45 = 0.665841;
    double yield46 = 0.7376;
    double yield47 = 0.813047;
    double yield48 = 0.892182;
    double yield49 = 0.975005;
    double yield50 = 1.06152;
    double yield51 = 1.15172;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST1500_0OSSF_ft_"+coeff[i]+".txt");
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

void makeDatacard_ST1500_1OSSF()
{
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    int size = sizeof(coeff)/sizeof(TString);
    
    double yield1 = 1.83963;
    double yield2 = 1.69603;
    double yield3 = 1.55829;
    double yield4 = 1.42638;
    double yield5 = 1.30032;
    double yield6 = 1.1801;
    double yield7 = 1.06573;
    double yield8 = 0.957197;
    double yield9 = 0.85451;
    double yield10 = 0.757667;
    double yield11 = 0.666667;
    double yield12 = 0.581511;
    double yield13 = 0.502198;
    double yield14 = 0.428729;
    double yield15 = 0.361104;
    double yield16 = 0.299322;
    double yield17 = 0.243384;
    double yield18 = 0.19329;
    double yield19 = 0.149039;
    double yield20 = 0.110631;
    double yield21 = 0.0780676;
    double yield22 = 0.0513475;
    double yield23 = 0.030471;
    double yield24 = 0.0154381;
    double yield25 = 0.00624879;
    double yield26 = 0.00290311;
    double yield27 = 0.00540103;
    double yield28 = 0.0137426;
    double yield29 = 0.0279277;
    double yield30 = 0.0479564;
    double yield31 = 0.0738288;
    double yield32 = 0.105545;
    double yield33 = 0.143104;
    double yield34 = 0.186507;
    double yield35 = 0.235754;
    double yield36 = 0.290845;
    double yield37 = 0.351779;
    double yield38 = 0.418556;
    double yield39 = 0.491177;
    double yield40 = 0.569642;
    double yield41 = 0.65395;
    double yield42 = 0.744102;
    double yield43 = 0.840098;
    double yield44 = 0.941937;
    double yield45 = 1.04962;
    double yield46 = 1.16315;
    double yield47 = 1.28252;
    double yield48 = 1.40773;
    double yield49 = 1.53879;
    double yield50 = 1.67569;
    double yield51 = 1.81843;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST1500_1OSSF_ft_"+coeff[i]+".txt");
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

void makeDatacard_ST1500_2OSSF()
{
    TString coeff[51] = {"neg2_5", "neg2_4", "neg2_3", "neg2_2", "neg2_1", "neg2_0", "neg1_9", "neg1_8", "neg1_7", "neg1_6", "neg1_5", "neg1_4", "neg1_3", "neg1_2", "neg1_1", "neg1_0", "neg0_9", "neg0_8", "neg0_7", "neg0_6", "neg0_5", "neg0_4", "neg0_3", "neg0_2", "neg0_1", "0_0", "pos0_1", "pos0_2", "pos0_3", "pos0_4", "pos0_5", "pos0_6", "pos0_7", "pos0_8", "pos0_9", "pos1_0", "pos1_1", "pos1_2", "pos1_3", "pos1_4", "pos1_5", "pos1_6", "pos1_7", "pos1_8", "pos1_9", "pos2_0", "pos2_1", "pos2_2", "pos2_3", "pos2_4", "pos2_5"};
    int size = sizeof(coeff)/sizeof(TString);
    double yield1 = 1.49024;
    double yield2 = 1.3739;
    double yield3 = 1.2623;
    double yield4 = 1.15544;
    double yield5 = 1.05332;
    double yield6 = 0.955941;
    double yield7 = 0.863299;
    double yield8 = 0.775396;
    double yield9 = 0.692233;
    double yield10 = 0.613808;
    double yield11 = 0.540124;
    double yield12 = 0.471178;
    double yield13 = 0.406972;
    double yield14 = 0.347506;
    double yield15 = 0.292779;
    double yield16 = 0.242791;
    double yield17 = 0.197542;
    double yield18 = 0.157033;
    double yield19 = 0.121263;
    double yield20 = 0.0902325;
    double yield21 = 0.0639413;
    double yield22 = 0.0423895;
    double yield23 = 0.0255771;
    double yield24 = 0.013504;
    double yield25 = 0.00617031;
    double yield26 = 0.00357595;
    double yield27 = 0.00572095;
    double yield28 = 0.0126053;
    double yield29 = 0.024229;
    double yield30 = 0.0405921;
    double yield31 = 0.0616945;
    double yield32 = 0.0875363;
    double yield33 = 0.118117;
    double yield34 = 0.153438;
    double yield35 = 0.193498;
    double yield36 = 0.238297;
    double yield37 = 0.287835;
    double yield38 = 0.342113;
    double yield39 = 0.401131;
    double yield40 = 0.464887;
    double yield41 = 0.533383;
    double yield42 = 0.606619;
    double yield43 = 0.684593;
    double yield44 = 0.767307;
    double yield45 = 0.854761;
    double yield46 = 0.946953;
    double yield47 = 1.04389;
    double yield48 = 1.14556;
    double yield49 = 1.25197;
    double yield50 = 1.36312;
    double yield51 = 1.47901;
    
    double signal_beforeSub[51] = {yield1, yield2, yield3, yield4, yield5, yield6, yield7, yield8, yield9, yield10, yield11, yield12, yield13, yield14, yield15, yield16, yield17, yield18, yield19, yield20, yield21, yield22, yield23, yield24, yield25, yield26, yield27, yield28, yield29, yield30, yield31, yield32, yield33, yield34, yield35, yield36, yield37, yield38, yield39, yield40, yield41, yield42, yield43, yield44, yield45, yield46, yield47, yield48, yield49, yield50, yield51};
    
    double signal[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    std::cout << "signal_beforeSub[25] = " << signal_beforeSub[25] << std::endl;
    
    for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[25];
    
    signal[25] = 0.000001;
    
    double bkg[51] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double bkg_beforeAdd[51] = {0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001, 0.00000001};
    
    for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[25];
    
    for(int i=0; i<size; i++)
    {
        std::ofstream outfile("realistic-counting-experiment_ST1500_2OSSF_ft_"+coeff[i]+".txt");
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

