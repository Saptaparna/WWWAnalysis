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

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {8.98, 6.44, 4.47, 3.07, 2.24, 1.98, 2.3, 3.19, 4.66, 6.7, 9.31};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2, 14.2};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

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

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {6.00, 3.84, 2.17, 0.982, 0.273, 0.0455, 0.299, 1.03, 2.25, 3.95, 6.12};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258, 0.258};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

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

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"}; 
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {4.45, 2.85, 1.6, 0.716, 0.184, 0.00985, 0.192, 0.731, 1.63, 2.88, 4.49};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321, 0.00321}; 
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

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

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {2.86, 1.83, 1.03, 0.458, 0.116, 0.00296, 0.12, 0.466, 1.04, 1.85, 2.88};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

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

  TString coeff[11] = {"neg2_5", "neg2_0", "neg1_5", "neg1_0", "neg0_5", "0_0", "pos0_5", "pos1_0", "pos1_5", "pos2_0", "pos2_5"};
  int size = sizeof(coeff)/sizeof(TString);

  double signal_beforeSub[11] = {1.26, 0.806, 0.453, 0.201, 0.05, 0.000344, 0.0519, 0.205, 0.459, 0.814, 1.27};
  double signal[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  std::cout << "signal_beforeSub[5] = " << signal_beforeSub[5] << std::endl;

  for(int i=0; i<size; i++) signal[i] = signal_beforeSub[i] - signal_beforeSub[5];

  signal[5] = 0.000001;
  double bkg[11] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double bkg_beforeAdd[11] = {0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001};
  for(int i=0; i<size; i++) bkg[i] = bkg_beforeAdd[i] + signal_beforeSub[5];

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


