#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TCanvas.h>
using std::string;

void makeGraphs_pValue_All()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
   
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double noST_p_value1 = 0.134554;
    double noST_p_value2 = 0.15191;
    double noST_p_value3 = 0.170217;
    double noST_p_value4 = 0.18934;
    double noST_p_value5 = 0.209124;
    double noST_p_value6 = 0.229401;
    double noST_p_value7 = 0.249992;
    double noST_p_value8 = 0.270709;
    double noST_p_value9 = 0.291366;
    double noST_p_value10 = 0.311782;
    double noST_p_value11 = 0.331771;
    double noST_p_value12 = 0.35117;
    double noST_p_value13 = 0.369823;
    double noST_p_value14 = 0.387588;
    double noST_p_value15 = 0.404341;
    double noST_p_value16 = 0.419969;
    double noST_p_value17 = 0.43438;
    double noST_p_value18 = 0.447493;
    double noST_p_value19 = 0.459242;
    double noST_p_value20 = 0.469572;
    double noST_p_value21 = 0.478443;
    double noST_p_value22 = 0.485821;
    double noST_p_value23 = 0.491655;
    double noST_p_value24 = 0.496442;
    double noST_p_value25 = 0.498914;
    double noST_p_value26 = 0.5;
    double noST_p_value27 = 0.499698;
    double noST_p_value28 = 0.498008;
    double noST_p_value29 = 0.494363;
    double noST_p_value30 = 0.489284;
    double noST_p_value31 = 0.482795;
    double noST_p_value32 = 0.474764;
    double noST_p_value33 = 0.465253;
    double noST_p_value34 = 0.454299;
    double noST_p_value35 = 0.441949;
    double noST_p_value36 = 0.428263;
    double noST_p_value37 = 0.413312;
    double noST_p_value38 = 0.397181;
    double noST_p_value39 = 0.379973;
    double noST_p_value40 = 0.361805;
    double noST_p_value41 = 0.34281;
    double noST_p_value42 = 0.323134;
    double noST_p_value43 = 0.302939;
    double noST_p_value44 = 0.282395;
    double noST_p_value45 = 0.261688;
    double noST_p_value46 = 0.241003;
    double noST_p_value47 = 0.220526;
    double noST_p_value48 = 0.200443;
    double noST_p_value49 = 0.180926;
    double noST_p_value50 = 0.162141;
    double noST_p_value51 = 0.144232;

    double p_value_noST[51] = {noST_p_value1, noST_p_value2, noST_p_value3, noST_p_value4, noST_p_value5, noST_p_value6, noST_p_value7, noST_p_value8, noST_p_value9, noST_p_value10, noST_p_value11, noST_p_value12, noST_p_value13, noST_p_value14, noST_p_value15, noST_p_value16, noST_p_value17, noST_p_value18, noST_p_value19, noST_p_value20, noST_p_value21, noST_p_value22, noST_p_value23, noST_p_value24, noST_p_value25, noST_p_value26, noST_p_value27, noST_p_value28, noST_p_value29, noST_p_value30, noST_p_value31, noST_p_value32, noST_p_value33, noST_p_value34, noST_p_value35, noST_p_value36, noST_p_value37, noST_p_value38, noST_p_value39, noST_p_value40, noST_p_value41, noST_p_value42, noST_p_value43, noST_p_value44, noST_p_value45, noST_p_value46, noST_p_value47, noST_p_value48, noST_p_value49, noST_p_value50, noST_p_value51}; 

    double ST1000_p_value1 = 4.24651e-07;
    double ST1000_p_value2 = 1.54217e-06;
    double ST1000_p_value3 = 5.22373e-06;
    double ST1000_p_value4 = 1.6506e-05;
    double ST1000_p_value5 = 4.86661e-05;
    double ST1000_p_value6 = 0.000133939;
    double ST1000_p_value7 = 0.000344316;
    double ST1000_p_value8 = 0.000827447;
    double ST1000_p_value9 = 0.00186068;
    double ST1000_p_value10 = 0.00392012;
    double ST1000_p_value11 = 0.00774891;
    double ST1000_p_value12 = 0.0143952;
    double ST1000_p_value13 = 0.0251803;
    double ST1000_p_value14 = 0.0415623;
    double ST1000_p_value15 = 0.0648885;
    double ST1000_p_value16 = 0.0960749;
    double ST1000_p_value17 = 0.135281;
    double ST1000_p_value18 = 0.181702;
    double ST1000_p_value19 = 0.23352;
    double ST1000_p_value20 = 0.288077;
    double ST1000_p_value21 = 0.342178;
    double ST1000_p_value22 = 0.392472;
    double ST1000_p_value23 = 0.435804;
    double ST1000_p_value24 = 0.469456;
    double ST1000_p_value25 = 0.491316;
    double ST1000_p_value26 = 0.5;
    double ST1000_p_value27 = 0.495003;
    double ST1000_p_value28 = 0.476497;
    double ST1000_p_value29 = 0.445805;
    double ST1000_p_value30 = 0.404806;
    double ST1000_p_value31 = 0.356071;
    double ST1000_p_value32 = 0.302659;
    double ST1000_p_value33 = 0.247898;
    double ST1000_p_value34 = 0.195063;
    double ST1000_p_value35 = 0.146989;
    double ST1000_p_value36 = 0.105741;
    double ST1000_p_value37 = 0.0723995;
    double ST1000_p_value38 = 0.0470446;
    double ST1000_p_value39 = 0.0289339;
    double ST1000_p_value40 = 0.0168021;
    double ST1000_p_value41 = 0.0091925;
    double ST1000_p_value42 = 0.00472869;
    double ST1000_p_value43 = 0.00228325;
    double ST1000_p_value44 = 0.00103323;
    double ST1000_p_value45 = 0.000437669;
    double ST1000_p_value46 = 0.000173349;
    double ST1000_p_value47 = 6.41436e-05;
    double ST1000_p_value48 = 2.21592e-05;
    double ST1000_p_value49 = 7.14367e-06;
    double ST1000_p_value50 = 2.14834e-06;
    double ST1000_p_value51 = 6.02642e-07;

    double p_value_ST1000[51] = {ST1000_p_value1, ST1000_p_value2, ST1000_p_value3, ST1000_p_value4, ST1000_p_value5, ST1000_p_value6, ST1000_p_value7, ST1000_p_value8, ST1000_p_value9, ST1000_p_value10, ST1000_p_value11, ST1000_p_value12, ST1000_p_value13, ST1000_p_value14, ST1000_p_value15, ST1000_p_value16, ST1000_p_value17, ST1000_p_value18, ST1000_p_value19, ST1000_p_value20, ST1000_p_value21, ST1000_p_value22, ST1000_p_value23, ST1000_p_value24, ST1000_p_value25, ST1000_p_value26, ST1000_p_value27, ST1000_p_value28, ST1000_p_value29, ST1000_p_value30, ST1000_p_value31, ST1000_p_value32, ST1000_p_value33, ST1000_p_value34, ST1000_p_value35, ST1000_p_value36, ST1000_p_value37, ST1000_p_value38, ST1000_p_value39, ST1000_p_value40, ST1000_p_value41, ST1000_p_value42, ST1000_p_value43, ST1000_p_value44, ST1000_p_value45, ST1000_p_value46, ST1000_p_value47, ST1000_p_value48, ST1000_p_value49, ST1000_p_value50, ST1000_p_value51};
   
    double ST1500_p_value1 = 4.00434e-11;
    double ST1500_p_value2 = 2.89704e-10;
    double ST1500_p_value3 = 1.90293e-09;
    double ST1500_p_value4 = 1.13526e-08;
    double ST1500_p_value5 = 6.15367e-08;
    double ST1500_p_value6 = 3.03298e-07;
    double ST1500_p_value7 = 1.3601e-06;
    double ST1500_p_value8 = 5.55436e-06;
    double ST1500_p_value9 = 2.06777e-05;
    double ST1500_p_value10 = 7.02526e-05;
    double ST1500_p_value11 = 0.000218133;
    double ST1500_p_value12 = 0.000619906;
    double ST1500_p_value13 = 0.0016153;
    double ST1500_p_value14 = 0.0038669;
    double ST1500_p_value15 = 0.00852296;
    double ST1500_p_value16 = 0.0173394;
    double ST1500_p_value17 = 0.0326507;
    double ST1500_p_value18 = 0.0570807;
    double ST1500_p_value19 = 0.0929579;
    double ST1500_p_value20 = 0.141534;
    double ST1500_p_value21 = 0.202242;
    double ST1500_p_value22 = 0.272258;
    double ST1500_p_value23 = 0.346448;
    double ST1500_p_value24 = 0.41741;
    double ST1500_p_value25 = 0.474421;
    double ST1500_p_value26 = 0.5;
    double ST1500_p_value27 = 0.47934;
    double ST1500_p_value28 = 0.42488;
    double ST1500_p_value29 = 0.354882;
    double ST1500_p_value30 = 0.28064;
    double ST1500_p_value31 = 0.209845;
    double ST1500_p_value32 = 0.147887;
    double ST1500_p_value33 = 0.0978567;
    double ST1500_p_value34 = 0.0605638;
    double ST1500_p_value35 = 0.034931;
    double ST1500_p_value36 = 0.0187115;
    double ST1500_p_value37 = 0.00928034;
    double ST1500_p_value38 = 0.00424965;
    double ST1500_p_value39 = 0.00179231;
    double ST1500_p_value40 = 0.000694601;
    double ST1500_p_value41 = 0.000246861;
    double ST1500_p_value42 = 8.03233e-05;
    double ST1500_p_value43 = 2.3888e-05;
    double ST1500_p_value44 = 6.48441e-06;
    double ST1500_p_value45 = 1.6048e-06;
    double ST1500_p_value46 = 3.61729e-07;
    double ST1500_p_value47 = 7.41948e-08;
    double ST1500_p_value48 = 1.38371e-08;
    double ST1500_p_value49 = 2.34496e-09;
    double ST1500_p_value50 = 3.60972e-10;
    double ST1500_p_value51 = 5.04486e-11;

    double p_value_ST1500[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51}; 

    double ST2000_p_value1 = 5.26354e-09;
    double ST2000_p_value2 = 2.40529e-08;
    double ST2000_p_value3 = 1.02188e-07;
    double ST2000_p_value4 = 4.03794e-07;
    double ST2000_p_value5 = 1.4849e-06;
    double ST2000_p_value6 = 5.08491e-06;
    double ST2000_p_value7 = 1.62277e-05;
    double ST2000_p_value8 = 4.83017e-05;
    double ST2000_p_value9 = 0.0001342;
    double ST2000_p_value10 = 0.000348381;
    double ST2000_p_value11 = 0.000845935;
    double ST2000_p_value12 = 0.00192378;
    double ST2000_p_value13 = 0.00410254;
    double ST2000_p_value14 = 0.00821606;
    double ST2000_p_value15 = 0.0154771;
    double ST2000_p_value16 = 0.0274723;
    double ST2000_p_value17 = 0.0460383;
    double ST2000_p_value18 = 0.0729928;
    double ST2000_p_value19 = 0.109739;
    double ST2000_p_value20 = 0.156826;
    double ST2000_p_value21 = 0.213577;
    double ST2000_p_value22 = 0.277873;
    double ST2000_p_value23 = 0.346135;
    double ST2000_p_value24 = 0.413172;
    double ST2000_p_value25 = 0.470837;
    double ST2000_p_value26 = 0.5;
    double ST2000_p_value27 = 0.475137;
    double ST2000_p_value28 = 0.418986;
    double ST2000_p_value29 = 0.352391;
    double ST2000_p_value30 = 0.284001;
    double ST2000_p_value31 = 0.219171;
    double ST2000_p_value32 = 0.161622;
    double ST2000_p_value33 = 0.113605;
    double ST2000_p_value34 = 0.0759222;
    double ST2000_p_value35 = 0.0481236;
    double ST2000_p_value36 = 0.0288647;
    double ST2000_p_value37 = 0.0163483;
    double ST2000_p_value38 = 0.00872642;
    double ST2000_p_value39 = 0.00438208;
    double ST2000_p_value40 = 0.00206681;
    double ST2000_p_value41 = 0.00091426;
    double ST2000_p_value42 = 0.000378782;
    double ST2000_p_value43 = 0.000146803;
    double ST2000_p_value44 = 5.31692e-05;
    double ST2000_p_value45 = 1.79766e-05;
    double ST2000_p_value46 = 5.66915e-06;
    double ST2000_p_value47 = 1.66617e-06;
    double ST2000_p_value48 = 4.56063e-07;
    double ST2000_p_value49 = 1.1618e-07;
    double ST2000_p_value50 = 2.75288e-08;
    double ST2000_p_value51 = 6.06462e-09;

    double p_value_ST2000[51] = {ST2000_p_value1, ST2000_p_value2, ST2000_p_value3, ST2000_p_value4, ST2000_p_value5, ST2000_p_value6, ST2000_p_value7, ST2000_p_value8, ST2000_p_value9, ST2000_p_value10, ST2000_p_value11, ST2000_p_value12, ST2000_p_value13, ST2000_p_value14, ST2000_p_value15, ST2000_p_value16, ST2000_p_value17, ST2000_p_value18, ST2000_p_value19, ST2000_p_value20, ST2000_p_value21, ST2000_p_value22, ST2000_p_value23, ST2000_p_value24, ST2000_p_value25, ST2000_p_value26, ST2000_p_value27, ST2000_p_value28, ST2000_p_value29, ST2000_p_value30, ST2000_p_value31, ST2000_p_value32, ST2000_p_value33, ST2000_p_value34, ST2000_p_value35, ST2000_p_value36, ST2000_p_value37, ST2000_p_value38, ST2000_p_value39, ST2000_p_value40, ST2000_p_value41, ST2000_p_value42, ST2000_p_value43, ST2000_p_value44, ST2000_p_value45, ST2000_p_value46, ST2000_p_value47, ST2000_p_value48, ST2000_p_value49, ST2000_p_value50, ST2000_p_value51}; 

    double ST2500_p_value1 = 1.80647e-05;
    double ST2500_p_value2 = 4.02304e-05;
    double ST2500_p_value3 = 8.63723e-05;
    double ST2500_p_value4 = 0.000178825;
    double ST2500_p_value5 = 0.000357163;
    double ST2500_p_value6 = 0.000688458;
    double ST2500_p_value7 = 0.00128132;
    double ST2500_p_value8 = 0.0023036;
    double ST2500_p_value9 = 0.00400267;
    double ST2500_p_value10 = 0.00672553;
    double ST2500_p_value11 = 0.0109341;
    double ST2500_p_value12 = 0.017211;
    double ST2500_p_value13 = 0.0262469;
    double ST2500_p_value14 = 0.038808;
    double ST2500_p_value15 = 0.0556754;
    double ST2500_p_value16 = 0.0775656;
    double ST2500_p_value17 = 0.105031;
    double ST2500_p_value18 = 0.138359;
    double ST2500_p_value19 = 0.177481;
    double ST2500_p_value20 = 0.221917;
    double ST2500_p_value21 = 0.270725;
    double ST2500_p_value22 = 0.322519;
    double ST2500_p_value23 = 0.375439;
    double ST2500_p_value24 = 0.426985;
    double ST2500_p_value25 = 0.473052;
    double ST2500_p_value26 = 0.499999;
    double ST2500_p_value27 = 0.477005;
    double ST2500_p_value28 = 0.431871;
    double ST2500_p_value29 = 0.380639;
    double ST2500_p_value30 = 0.327734;
    double ST2500_p_value31 = 0.275736;
    double ST2500_p_value32 = 0.226564;
    double ST2500_p_value33 = 0.181648;
    double ST2500_p_value34 = 0.141971;
    double ST2500_p_value35 = 0.108062;
    double ST2500_p_value36 = 0.0800244;
    double ST2500_p_value37 = 0.0576044;
    double ST2500_p_value38 = 0.0402707;
    double ST2500_p_value39 = 0.0273185;
    double ST2500_p_value40 = 0.0179691;
    double ST2500_p_value41 = 0.0114519;
    double ST2500_p_value42 = 0.00706667;
    double ST2500_p_value43 = 0.00421957;
    double ST2500_p_value44 = 0.00243654;
    double ST2500_p_value45 = 0.00135987;
    double ST2500_p_value46 = 0.000733188;
    double ST2500_p_value47 = 0.000381694;
    double ST2500_p_value48 = 0.000191783;
    double ST2500_p_value49 = 9.29637e-05;
    double ST2500_p_value50 = 4.34569e-05;
    double ST2500_p_value51 = 1.95846e-05;

    double p_value_ST2500[51] = {ST2500_p_value1, ST2500_p_value2, ST2500_p_value3, ST2500_p_value4, ST2500_p_value5, ST2500_p_value6, ST2500_p_value7, ST2500_p_value8, ST2500_p_value9, ST2500_p_value10, ST2500_p_value11, ST2500_p_value12, ST2500_p_value13, ST2500_p_value14, ST2500_p_value15, ST2500_p_value16, ST2500_p_value17, ST2500_p_value18, ST2500_p_value19, ST2500_p_value20, ST2500_p_value21, ST2500_p_value22, ST2500_p_value23, ST2500_p_value24, ST2500_p_value25, ST2500_p_value26, ST2500_p_value27, ST2500_p_value28, ST2500_p_value29, ST2500_p_value30, ST2500_p_value31, ST2500_p_value32, ST2500_p_value33, ST2500_p_value34, ST2500_p_value35, ST2500_p_value36, ST2500_p_value37, ST2500_p_value38, ST2500_p_value39, ST2500_p_value40, ST2500_p_value41, ST2500_p_value42, ST2500_p_value43, ST2500_p_value44, ST2500_p_value45, ST2500_p_value46, ST2500_p_value47, ST2500_p_value48, ST2500_p_value49, ST2500_p_value50, ST2500_p_value51};

    double ST1500_p_value_Signal10Unc1 = 4.00433e-11;
    double ST1500_p_value_Signal10Unc2 = 2.89704e-10;
    double ST1500_p_value_Signal10Unc3 = 1.90293e-09;
    double ST1500_p_value_Signal10Unc4 = 1.13526e-08;
    double ST1500_p_value_Signal10Unc5 = 6.15367e-08;
    double ST1500_p_value_Signal10Unc6 = 3.03298e-07;
    double ST1500_p_value_Signal10Unc7 = 1.3601e-06;
    double ST1500_p_value_Signal10Unc8 = 5.55436e-06;
    double ST1500_p_value_Signal10Unc9 = 2.06778e-05;
    double ST1500_p_value_Signal10Unc10 = 7.02528e-05;
    double ST1500_p_value_Signal10Unc11 = 0.000218134;
    double ST1500_p_value_Signal10Unc12 = 0.000619907;
    double ST1500_p_value_Signal10Unc13 = 0.00161531;
    double ST1500_p_value_Signal10Unc14 = 0.0038669;
    double ST1500_p_value_Signal10Unc15 = 0.00852297;
    double ST1500_p_value_Signal10Unc16 = 0.0173394;
    double ST1500_p_value_Signal10Unc17 = 0.0326507;
    double ST1500_p_value_Signal10Unc18 = 0.0570807;
    double ST1500_p_value_Signal10Unc19 = 0.0929579;
    double ST1500_p_value_Signal10Unc20 = 0.141534;
    double ST1500_p_value_Signal10Unc21 = 0.202242;
    double ST1500_p_value_Signal10Unc22 = 0.272259;
    double ST1500_p_value_Signal10Unc23 = 0.346449;
    double ST1500_p_value_Signal10Unc24 = 0.41741;
    double ST1500_p_value_Signal10Unc25 = 0.474423;
    double ST1500_p_value_Signal10Unc26 = 0.5;
    double ST1500_p_value_Signal10Unc27 = 0.479342;
    double ST1500_p_value_Signal10Unc28 = 0.42488;
    double ST1500_p_value_Signal10Unc29 = 0.354883;
    double ST1500_p_value_Signal10Unc30 = 0.280641;
    double ST1500_p_value_Signal10Unc31 = 0.209845;
    double ST1500_p_value_Signal10Unc32 = 0.147887;
    double ST1500_p_value_Signal10Unc33 = 0.0978567;
    double ST1500_p_value_Signal10Unc34 = 0.0605638;
    double ST1500_p_value_Signal10Unc35 = 0.034931;
    double ST1500_p_value_Signal10Unc36 = 0.0187115;
    double ST1500_p_value_Signal10Unc37 = 0.00928035;
    double ST1500_p_value_Signal10Unc38 = 0.00424965;
    double ST1500_p_value_Signal10Unc39 = 0.00179231;
    double ST1500_p_value_Signal10Unc40 = 0.000694602;
    double ST1500_p_value_Signal10Unc41 = 0.000246862;
    double ST1500_p_value_Signal10Unc42 = 8.03235e-05;
    double ST1500_p_value_Signal10Unc43 = 2.3888e-05;
    double ST1500_p_value_Signal10Unc44 = 6.48441e-06;
    double ST1500_p_value_Signal10Unc45 = 1.6048e-06;
    double ST1500_p_value_Signal10Unc46 = 3.61729e-07;
    double ST1500_p_value_Signal10Unc47 = 7.41948e-08; 
    double ST1500_p_value_Signal10Unc48 = 1.38371e-08;
    double ST1500_p_value_Signal10Unc49 = 2.34496e-09;
    double ST1500_p_value_Signal10Unc50 = 3.60972e-10;
    double ST1500_p_value_Signal10Unc51 = 5.04486e-11; 

    double ST1500_p_value_Signal10Unc[51] = {ST1500_p_value_Signal10Unc1, ST1500_p_value_Signal10Unc2, ST1500_p_value_Signal10Unc3, ST1500_p_value_Signal10Unc4, ST1500_p_value_Signal10Unc5, ST1500_p_value_Signal10Unc6, ST1500_p_value_Signal10Unc7, ST1500_p_value_Signal10Unc8, ST1500_p_value_Signal10Unc9, ST1500_p_value_Signal10Unc10, ST1500_p_value_Signal10Unc11, ST1500_p_value_Signal10Unc12, ST1500_p_value_Signal10Unc13, ST1500_p_value_Signal10Unc14, ST1500_p_value_Signal10Unc15, ST1500_p_value_Signal10Unc16, ST1500_p_value_Signal10Unc17, ST1500_p_value_Signal10Unc18, ST1500_p_value_Signal10Unc19, ST1500_p_value_Signal10Unc20, ST1500_p_value_Signal10Unc21, ST1500_p_value_Signal10Unc22, ST1500_p_value_Signal10Unc23, ST1500_p_value_Signal10Unc24, ST1500_p_value_Signal10Unc25, ST1500_p_value_Signal10Unc26, ST1500_p_value_Signal10Unc27, ST1500_p_value_Signal10Unc28, ST1500_p_value_Signal10Unc29, ST1500_p_value_Signal10Unc30, ST1500_p_value_Signal10Unc31, ST1500_p_value_Signal10Unc32, ST1500_p_value_Signal10Unc33, ST1500_p_value_Signal10Unc34, ST1500_p_value_Signal10Unc35, ST1500_p_value_Signal10Unc36, ST1500_p_value_Signal10Unc37, ST1500_p_value_Signal10Unc38, ST1500_p_value_Signal10Unc39, ST1500_p_value_Signal10Unc40, ST1500_p_value_Signal10Unc41, ST1500_p_value_Signal10Unc42, ST1500_p_value_Signal10Unc43, ST1500_p_value_Signal10Unc44, ST1500_p_value_Signal10Unc45, ST1500_p_value_Signal10Unc46, ST1500_p_value_Signal10Unc47, ST1500_p_value_Signal10Unc48, ST1500_p_value_Signal10Unc49, ST1500_p_value_Signal10Unc50, ST1500_p_value_Signal10Unc51};

    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(51, param, p_value_noST);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p_value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(51, param, p_value_ST1000);
    gr_1->SetMarkerStyle(20);
    gr_1->SetMarkerColor(kCyan);
    gr_1->SetLineColor(kCyan);
    gr_1->SetMarkerSize(0.9);
    gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(51, param, p_value_ST1500);
    gr_2->SetMarkerStyle(20);
    gr_2->SetMarkerColor(kBlack);
    gr_2->SetLineColor(kBlack);
    gr_2->SetMarkerSize(0.9);
    gr_2->Draw("LP*");

    TGraph *gr_3 = new TGraph(51, param, p_value_ST2000);
    gr_3->SetMarkerStyle(20);
    gr_3->SetMarkerColor(kGreen+3);
    gr_3->SetLineColor(kGreen+3);
    gr_3->SetMarkerSize(0.9);
    gr_3->Draw("LP*");

    TGraph *gr_4 = new TGraph(51, param, p_value_ST2500);
    gr_4->SetMarkerStyle(20);
    gr_4->SetMarkerColor(kMagenta);
    gr_4->SetLineColor(kMagenta);
    gr_4->SetMarkerSize(0.9);
    gr_4->Draw("LP*");

    TGraph *gr_6 = new TGraph(51, param, ST1500_p_value_Signal10Unc);
    gr_6->SetMarkerStyle(20);
    gr_6->SetMarkerColor(kAzure+4);
    gr_6->SetLineColor(kAzure+4);
    gr_6->SetMarkerSize(0.9);
    gr_6->Draw("LP*");

    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.65,0.70,0.89,0.89,NULL,"brNDC");
    leg->AddEntry(gr, "no ST cut", "pl");
    leg->AddEntry(gr_1, "ST 1000", "pl");
    leg->AddEntry(gr_2, "ST 1500", "pl");
    leg->AddEntry(gr_6, "ST 1500: Signal Syst 10%", "pl");
    leg->AddEntry(gr_3, "ST 2000", "pl");
    leg->AddEntry(gr_4, "ST 2500", "pl");
    leg->SetBorderSize(0);
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr_2->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit ST 1500 = " << m << endl;
    
    while (m<=2.0 and gr_2->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit ST 1500 = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_3->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit ST 2000 = " << m1 << endl;

    while (m1<=2.0 and gr_3->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit ST 2000 = " << m1 << endl;

    c1.SaveAs("Parabola_ft0_p_value_Finegrid_All.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Finegrid_All.png");
    
}

void makeGraphs_pValue_combined()
{
   
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
   
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5}; 
 
    double ST_p_value1 = 2.00092e-14;
    double ST_p_value2 = 2.96911e-13;
    double ST_p_value3 = 3.87444e-12;
    double ST_p_value4 = 4.43719e-11;
    double ST_p_value5 = 4.45959e-10;
    double ST_p_value6 = 3.93108e-09;
    double ST_p_value7 = 3.03592e-08;
    double ST_p_value8 = 2.05597e-07;
    double ST_p_value9 = 1.2211e-06;
    double ST_p_value10 = 6.36445e-06;
    double ST_p_value11 = 2.91611e-05;
    double ST_p_value12 = 0.00011762;
    double ST_p_value13 = 0.00041846;
    double ST_p_value14 = 0.00131697;
    double ST_p_value15 = 0.00367947;
    double ST_p_value16 = 0.00915728;
    double ST_p_value17 = 0.0203963;
    double ST_p_value18 = 0.0408639;
    double ST_p_value19 = 0.0741211;
    double ST_p_value20 = 0.122441;
    double ST_p_value21 = 0.185583;
    double ST_p_value22 = 0.260065;
    double ST_p_value23 = 0.339134;
    double ST_p_value24 = 0.41418;
    double ST_p_value25 = 0.473493;
    double ST_p_value26 = 0.5;
    double ST_p_value27 = 0.478901;
    double ST_p_value28 = 0.422322;
    double ST_p_value29 = 0.348458;
    double ST_p_value30 = 0.269423;
    double ST_p_value31 = 0.194096;
    double ST_p_value32 = 0.129374;
    double ST_p_value33 = 0.0791853;
    double ST_p_value34 = 0.0442033;
    double ST_p_value35 = 0.0223485;
    double ST_p_value36 = 0.0101711;
    double ST_p_value37 = 0.00414478;
    double ST_p_value38 = 0.0015063;
    double ST_p_value39 = 0.000485981;
    double ST_p_value40 = 0.000138753;
    double ST_p_value41 = 3.49547e-05;
    double ST_p_value42 = 7.75757e-06;
    double ST_p_value43 = 1.51304e-06;
    double ST_p_value44 = 2.59008e-07;
    double ST_p_value45 = 3.88997e-08;
    double ST_p_value46 = 5.12168e-09;
    double ST_p_value47 = 5.911e-10;  
    double ST_p_value48 = 5.97987e-11;
    double ST_p_value49 = 5.30843e-12;
    double ST_p_value50 = 4.13423e-13;
    double ST_p_value51 = 2.83312e-14;
 
    double ST_p_value[51] = {ST_p_value1, ST_p_value2, ST_p_value3, ST_p_value4, ST_p_value5, ST_p_value6, ST_p_value7, ST_p_value8, ST_p_value9, ST_p_value10, ST_p_value11, ST_p_value12, ST_p_value13, ST_p_value14, ST_p_value15, ST_p_value16, ST_p_value17, ST_p_value18, ST_p_value19, ST_p_value20, ST_p_value21, ST_p_value22, ST_p_value23, ST_p_value24, ST_p_value25, ST_p_value26, ST_p_value27, ST_p_value28, ST_p_value29, ST_p_value30, ST_p_value31, ST_p_value32, ST_p_value33, ST_p_value34, ST_p_value35, ST_p_value36, ST_p_value37, ST_p_value38, ST_p_value39, ST_p_value40, ST_p_value41, ST_p_value42, ST_p_value43, ST_p_value44, ST_p_value45, ST_p_value46, ST_p_value47, ST_p_value48, ST_p_value49, ST_p_value50, ST_p_value51}; 

    double ST_p_value_Syst_401 = 1.08291e-13;
    double ST_p_value_Syst_402 = 1.25412e-12;
    double ST_p_value_Syst_403 = 1.31348e-11;
    double ST_p_value_Syst_404 = 1.23907e-10;
    double ST_p_value_Syst_405 = 1.05059e-09;
    double ST_p_value_Syst_406 = 7.98594e-09;
    double ST_p_value_Syst_407 = 5.42609e-08;
    double ST_p_value_Syst_408 = 3.29223e-07;
    double ST_p_value_Syst_409 = 1.78083e-06;
    double ST_p_value_Syst_4010 = 8.57882e-06;
    double ST_p_value_Syst_4011 = 3.68086e-05;
    double ST_p_value_Syst_4012 = 0.000140645;
    double ST_p_value_Syst_4013 = 0.000478843;
    double ST_p_value_Syst_4014 = 0.00145483;
    double ST_p_value_Syst_4015 = 0.00395347;
    double ST_p_value_Syst_4016 = 0.00963096;
    double ST_p_value_Syst_4017 = 0.021108;
    double ST_p_value_Syst_4018 = 0.0417924;
    double ST_p_value_Syst_4019 = 0.0751704;
    double ST_p_value_Syst_4020 = 0.123464;
    double ST_p_value_Syst_4021 = 0.186434;
    double ST_p_value_Syst_4022 = 0.26066;
    double ST_p_value_Syst_4023 = 0.339472;
    double ST_p_value_Syst_4024 = 0.414321;
    double ST_p_value_Syst_4025 = 0.473526;
    double ST_p_value_Syst_4026 = 0.5;
    double ST_p_value_Syst_4027 = 0.478903;
    double ST_p_value_Syst_4028 = 0.422441;
    double ST_p_value_Syst_4029 = 0.348761;
    double ST_p_value_Syst_4030 = 0.269978;
    double ST_p_value_Syst_4031 = 0.194911;
    double ST_p_value_Syst_4032 = 0.130375;
    double ST_p_value_Syst_4033 = 0.0802326;
    double ST_p_value_Syst_4034 = 0.045148;
    double ST_p_value_Syst_4035 = 0.0230862;
    double ST_p_value_Syst_4036 = 0.010671;
    double ST_p_value_Syst_4037 = 0.00443915;
    double ST_p_value_Syst_4038 = 0.00165713;
    double ST_p_value_Syst_4039 = 0.000553219;
    double ST_p_value_Syst_4040 = 0.000164849;
    double ST_p_value_Syst_4041 = 4.37751e-05;
    double ST_p_value_Syst_4042 = 1.03575e-05;
    double ST_p_value_Syst_4043 = 2.18165e-06;
    double ST_p_value_Syst_4044 = 4.09215e-07;
    double ST_p_value_Syst_4045 = 6.84404e-08;
    double ST_p_value_Syst_4046 = 1.02164e-08;
    double ST_p_value_Syst_4047 = 1.36348e-09;
    double ST_p_value_Syst_4048 = 1.63007e-10;
    double ST_p_value_Syst_4049 = 1.75094e-11;
    double ST_p_value_Syst_4050 = 1.693e-12;
    double ST_p_value_Syst_4051 = 1.48076e-13;

    double ST_p_value_Syst_40[51] = {ST_p_value_Syst_401, ST_p_value_Syst_402, ST_p_value_Syst_403, ST_p_value_Syst_404, ST_p_value_Syst_405, ST_p_value_Syst_406, ST_p_value_Syst_407, ST_p_value_Syst_408, ST_p_value_Syst_409, ST_p_value_Syst_4010, ST_p_value_Syst_4011, ST_p_value_Syst_4012, ST_p_value_Syst_4013, ST_p_value_Syst_4014, ST_p_value_Syst_4015, ST_p_value_Syst_4016, ST_p_value_Syst_4017, ST_p_value_Syst_4018, ST_p_value_Syst_4019, ST_p_value_Syst_4020, ST_p_value_Syst_4021, ST_p_value_Syst_4022, ST_p_value_Syst_4023, ST_p_value_Syst_4024, ST_p_value_Syst_4025, ST_p_value_Syst_4026, ST_p_value_Syst_4027, ST_p_value_Syst_4028, ST_p_value_Syst_4029, ST_p_value_Syst_4030, ST_p_value_Syst_4031, ST_p_value_Syst_4032, ST_p_value_Syst_4033, ST_p_value_Syst_4034, ST_p_value_Syst_4035, ST_p_value_Syst_4036, ST_p_value_Syst_4037, ST_p_value_Syst_4038, ST_p_value_Syst_4039, ST_p_value_Syst_4040, ST_p_value_Syst_4041, ST_p_value_Syst_4042, ST_p_value_Syst_4043, ST_p_value_Syst_4044, ST_p_value_Syst_4045, ST_p_value_Syst_4046, ST_p_value_Syst_4047, ST_p_value_Syst_4048, ST_p_value_Syst_4049, ST_p_value_Syst_4050, ST_p_value_Syst_4051};

    double ST_p_value_Syst_501 = 6.77183e-13;
    double ST_p_value_Syst_502 = 6.03191e-12;
    double ST_p_value_Syst_503 = 4.99621e-11;
    double ST_p_value_Syst_504 = 3.82674e-10;
    double ST_p_value_Syst_505 = 2.69993e-09;
    double ST_p_value_Syst_506 = 1.74724e-08;
    double ST_p_value_Syst_507 = 1.03226e-07;
    double ST_p_value_Syst_508 = 5.55216e-07;
    double ST_p_value_Syst_509 = 2.70962e-06;
    double ST_p_value_Syst_5010 = 1.19657e-05;
    double ST_p_value_Syst_5011 = 4.77383e-05;
    double ST_p_value_Syst_5012 = 0.000171767;
    double ST_p_value_Syst_5013 = 0.000556852;
    double ST_p_value_Syst_5014 = 0.00162663;
    double ST_p_value_Syst_5015 = 0.00428539;
    double ST_p_value_Syst_5016 = 0.0101924;
    double ST_p_value_Syst_5017 = 0.021938;
    double ST_p_value_Syst_5018 = 0.0428625;
    double ST_p_value_Syst_5019 = 0.0763698;
    double ST_p_value_Syst_5020 = 0.124626;
    double ST_p_value_Syst_5021 = 0.187398;
    double ST_p_value_Syst_5022 = 0.261332;
    double ST_p_value_Syst_5023 = 0.339852;
    double ST_p_value_Syst_5024 = 0.41448;
    double ST_p_value_Syst_5025 = 0.473563;
    double ST_p_value_Syst_5026 = 0.5;
    double ST_p_value_Syst_5027 = 0.47893;
    double ST_p_value_Syst_5028 = 0.422576;
    double ST_p_value_Syst_5029 = 0.349103;
    double ST_p_value_Syst_5030 = 0.270604;
    double ST_p_value_Syst_5031 = 0.195833;
    double ST_p_value_Syst_5032 = 0.131511;
    double ST_p_value_Syst_5033 = 0.0814286;
    double ST_p_value_Syst_5034 = 0.0462354;
    double ST_p_value_Syst_5035 = 0.0239448;
    double ST_p_value_Syst_5036 = 0.0112619;
    double ST_p_value_Syst_5037 = 0.00479457;
    double ST_p_value_Syst_5038 = 0.00184428;
    double ST_p_value_Syst_5039 = 0.000639623;
    double ST_p_value_Syst_5040 = 0.000199887;
    double ST_p_value_Syst_5041 = 5.62787e-05;
    double ST_p_value_Syst_5042 = 1.4295e-05;
    double ST_p_value_Syst_5043 = 3.27811e-06;
    double ST_p_value_Syst_5044 = 6.80003e-07;
    double ST_p_value_Syst_5045 = 1.27978e-07;
    double ST_p_value_Syst_5046 = 2.19119e-08;
    double ST_p_value_Syst_5047 = 3.42489e-09;
    double ST_p_value_Syst_5048 = 4.90509e-10;
    double ST_p_value_Syst_5049 = 6.46735e-11;
    double ST_p_value_Syst_5050 = 7.87882e-12;
    double ST_p_value_Syst_5051 = 8.92522e-13;

    double ST_p_value_Syst_50[51] = {ST_p_value_Syst_501, ST_p_value_Syst_502, ST_p_value_Syst_503, ST_p_value_Syst_504, ST_p_value_Syst_505, ST_p_value_Syst_506, ST_p_value_Syst_507, ST_p_value_Syst_508, ST_p_value_Syst_509, ST_p_value_Syst_5010, ST_p_value_Syst_5011, ST_p_value_Syst_5012, ST_p_value_Syst_5013, ST_p_value_Syst_5014, ST_p_value_Syst_5015, ST_p_value_Syst_5016, ST_p_value_Syst_5017, ST_p_value_Syst_5018, ST_p_value_Syst_5019, ST_p_value_Syst_5020, ST_p_value_Syst_5021, ST_p_value_Syst_5022, ST_p_value_Syst_5023, ST_p_value_Syst_5024, ST_p_value_Syst_5025, ST_p_value_Syst_5026, ST_p_value_Syst_5027, ST_p_value_Syst_5028, ST_p_value_Syst_5029, ST_p_value_Syst_5030, ST_p_value_Syst_5031, ST_p_value_Syst_5032, ST_p_value_Syst_5033, ST_p_value_Syst_5034, ST_p_value_Syst_5035, ST_p_value_Syst_5036, ST_p_value_Syst_5037, ST_p_value_Syst_5038, ST_p_value_Syst_5039, ST_p_value_Syst_5040, ST_p_value_Syst_5041, ST_p_value_Syst_5042, ST_p_value_Syst_5043, ST_p_value_Syst_5044, ST_p_value_Syst_5045, ST_p_value_Syst_5046, ST_p_value_Syst_5047, ST_p_value_Syst_5048, ST_p_value_Syst_5049, ST_p_value_Syst_5050, ST_p_value_Syst_5051};

    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};  
 
    TGraph *gr = new TGraph(51, param, ST_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    //gr->SetMinimum(0.02);
    gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");

    TGraph *gr_1 = new TGraph(51, param, ST_p_value_Syst_40);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");

    TGraph *gr_2 = new TGraph(51, param, ST_p_value_Syst_50);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");

    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");

    TLegend *leg = new TLegend(0.38,0.30,0.62,0.50,NULL,"brNDC");
    leg->AddEntry(gr,   "Combined Limit: Syst 30%", "pl");
    leg->AddEntry(gr_1, "Combined Limit: Syst 40%", "pl");
    leg->AddEntry(gr_2, "Combined Limit: Syst 50%", "pl");
    leg->SetBorderSize(0);
    leg->Draw();

    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 30% syst = " << m << endl;

    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 30% syst = " << m << endl;

    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 40% syst = " << m1 << endl;

    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 40% syst = " << m1 << endl;

    double m2=-2.0;
    while (m2<=2.0 and gr_2->Eval(m2)<gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected lower limit combined 50% syst = " << m2 << endl;

    while (m2<=2.0 and gr_2->Eval(m2)>gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected upper limit combined 50% syst = " << m2 << endl;


    c1.SaveAs("Parabola_ft0_p_value_All.pdf");
    c1.SaveAs("Parabola_ft0_p_value_All.png");
}

void makeGraphs_pValue_channelwiseSeparation()
{
   
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
   
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5}; 

    
    double ST_p_value1 = 2.00092e-14;
    double ST_p_value2 = 2.96911e-13;
    double ST_p_value3 = 3.87444e-12;
    double ST_p_value4 = 4.43719e-11;
    double ST_p_value5 = 4.45959e-10;
    double ST_p_value6 = 3.93108e-09;
    double ST_p_value7 = 3.03592e-08;
    double ST_p_value8 = 2.05597e-07;
    double ST_p_value9 = 1.2211e-06;
    double ST_p_value10 = 6.36445e-06;
    double ST_p_value11 = 2.91611e-05;
    double ST_p_value12 = 0.00011762;
    double ST_p_value13 = 0.00041846;
    double ST_p_value14 = 0.00131697;
    double ST_p_value15 = 0.00367947;
    double ST_p_value16 = 0.00915728;
    double ST_p_value17 = 0.0203963;
    double ST_p_value18 = 0.0408639;
    double ST_p_value19 = 0.0741211;
    double ST_p_value20 = 0.122441;
    double ST_p_value21 = 0.185583;
    double ST_p_value22 = 0.260065;
    double ST_p_value23 = 0.339134;
    double ST_p_value24 = 0.41418;
    double ST_p_value25 = 0.473493;
    double ST_p_value26 = 0.5;
    double ST_p_value27 = 0.478901;
    double ST_p_value28 = 0.422322;
    double ST_p_value29 = 0.348458;
    double ST_p_value30 = 0.269423;
    double ST_p_value31 = 0.194096;
    double ST_p_value32 = 0.129374;
    double ST_p_value33 = 0.0791853;
    double ST_p_value34 = 0.0442033;
    double ST_p_value35 = 0.0223485;
    double ST_p_value36 = 0.0101711;
    double ST_p_value37 = 0.00414478;
    double ST_p_value38 = 0.0015063;
    double ST_p_value39 = 0.000485981;
    double ST_p_value40 = 0.000138753;
    double ST_p_value41 = 3.49547e-05;
    double ST_p_value42 = 7.75757e-06;
    double ST_p_value43 = 1.51304e-06;
    double ST_p_value44 = 2.59008e-07;
    double ST_p_value45 = 3.88997e-08;
    double ST_p_value46 = 5.12168e-09;
    double ST_p_value47 = 5.911e-10;  
    double ST_p_value48 = 5.97987e-11;
    double ST_p_value49 = 5.30843e-12;
    double ST_p_value50 = 4.13423e-13;
    double ST_p_value51 = 2.83312e-14;
 
    double ST_p_value[51] = {ST_p_value1, ST_p_value2, ST_p_value3, ST_p_value4, ST_p_value5, ST_p_value6, ST_p_value7, ST_p_value8, ST_p_value9, ST_p_value10, ST_p_value11, ST_p_value12, ST_p_value13, ST_p_value14, ST_p_value15, ST_p_value16, ST_p_value17, ST_p_value18, ST_p_value19, ST_p_value20, ST_p_value21, ST_p_value22, ST_p_value23, ST_p_value24, ST_p_value25, ST_p_value26, ST_p_value27, ST_p_value28, ST_p_value29, ST_p_value30, ST_p_value31, ST_p_value32, ST_p_value33, ST_p_value34, ST_p_value35, ST_p_value36, ST_p_value37, ST_p_value38, ST_p_value39, ST_p_value40, ST_p_value41, ST_p_value42, ST_p_value43, ST_p_value44, ST_p_value45, ST_p_value46, ST_p_value47, ST_p_value48, ST_p_value49, ST_p_value50, ST_p_value51}; 

    double ST_ch_p_value1 = 2.08461e-15;
    double ST_ch_p_value2 = 3.7339e-14;
    double ST_ch_p_value3 = 5.82712e-13;
    double ST_ch_p_value4 = 7.9261e-12;
    double ST_ch_p_value5 = 9.38835e-11;
    double ST_ch_p_value6 = 9.66164e-10;
    double ST_ch_p_value7 = 8.65168e-09;
    double ST_ch_p_value8 = 6.73493e-08;
    double ST_ch_p_value9 = 4.55982e-07;
    double ST_ch_p_value10 = 2.68752e-06;
    double ST_ch_p_value11 = 1.38098e-05;
    double ST_ch_p_value12 = 6.19581e-05;
    double ST_ch_p_value13 = 0.00024301;
    double ST_ch_p_value14 = 0.000836966;
    double ST_ch_p_value15 = 0.00253467;
    double ST_ch_p_value16 = 0.00678216;
    double ST_ch_p_value17 = 0.0161013;
    double ST_ch_p_value18 = 0.0341198;
    double ST_ch_p_value19 = 0.0649089;
    double ST_ch_p_value20 = 0.11158;
    double ST_ch_p_value21 = 0.174508;
    double ST_ch_p_value22 = 0.250724;
    double ST_ch_p_value23 = 0.332603;
    double ST_ch_p_value24 = 0.410993;
    double ST_ch_p_value25 = 0.472776;
    double ST_ch_p_value26 = 0.5;
    double ST_ch_p_value27 = 0.478054;
    double ST_ch_p_value28 = 0.419224;
    double ST_ch_p_value29 = 0.342264;
    double ST_ch_p_value30 = 0.260124;
    double ST_ch_p_value31 = 0.183114;
    double ST_ch_p_value32 = 0.118309;
    double ST_ch_p_value33 = 0.0696629;
    double ST_ch_p_value34 = 0.0370952;
    double ST_ch_p_value35 = 0.0177478;
    double ST_ch_p_value36 = 0.0075878;
    double ST_ch_p_value37 = 0.00287588;
    double ST_ch_p_value38 = 0.000964501;
    double ST_ch_p_value39 = 0.000284746;
    double ST_ch_p_value40 = 7.38322e-05;
    double ST_ch_p_value41 = 1.67463e-05;
    double ST_ch_p_value42 = 3.31281e-06;
    double ST_ch_p_value43 = 5.72072e-07;
    double ST_ch_p_value44 = 8.59723e-08;
    double ST_ch_p_value45 = 1.12348e-08;
    double ST_ch_p_value46 = 1.27646e-09;
    double ST_ch_p_value47 = 1.26203e-10;
    double ST_ch_p_value48 = 1.08505e-11;
    double ST_ch_p_value49 = 8.11093e-13;
    double ST_ch_p_value50 = 5.28601e-14;
    double ST_ch_p_value51 = 3.00391e-15;

    double ST_ch_p_value[51] = {ST_ch_p_value1, ST_ch_p_value2, ST_ch_p_value3, ST_ch_p_value4, ST_ch_p_value5, ST_ch_p_value6, ST_ch_p_value7, ST_ch_p_value8, ST_ch_p_value9, ST_ch_p_value10, ST_ch_p_value11, ST_ch_p_value12, ST_ch_p_value13, ST_ch_p_value14, ST_ch_p_value15, ST_ch_p_value16, ST_ch_p_value17, ST_ch_p_value18, ST_ch_p_value19, ST_ch_p_value20, ST_ch_p_value21, ST_ch_p_value22, ST_ch_p_value23, ST_ch_p_value24, ST_ch_p_value25, ST_ch_p_value26, ST_ch_p_value27, ST_ch_p_value28, ST_ch_p_value29, ST_ch_p_value30, ST_ch_p_value31, ST_ch_p_value32, ST_ch_p_value33, ST_ch_p_value34, ST_ch_p_value35, ST_ch_p_value36, ST_ch_p_value37, ST_ch_p_value38, ST_ch_p_value39, ST_ch_p_value40, ST_ch_p_value41, ST_ch_p_value42, ST_ch_p_value43, ST_ch_p_value44, ST_ch_p_value45, ST_ch_p_value46, ST_ch_p_value47, ST_ch_p_value48, ST_ch_p_value49, ST_ch_p_value50, ST_ch_p_value51};

    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};  
 
    TGraph *gr = new TGraph(51, param, ST_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    //gr->SetMinimum(0.02);
    gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");

    TGraph *gr_1 = new TGraph(51, param, ST_ch_p_value);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");

    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");

    TLegend *leg = new TLegend(0.38,0.30,0.62,0.60,NULL,"brNDC");
    //TLegend *leg = new TLegend(0.38,0.30,0.62,0.55,NULL,"brNDC");
    leg->SetHeader("Combined Limit");
    leg->AddEntry(gr,   "Syst 30%", "pl");
    leg->AddEntry(gr_1,   "Separated by channel", "pl");
    leg->SetBorderSize(0);
    leg->Draw();

    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 30% syst = " << m << endl;

    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 30% syst = " << m << endl;

    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined channel-wise = " << m1 << endl;

    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined channel-wise = " << m1 << endl;

    c1.SaveAs("Parabola_ft0_p_value_SeparateChannelwise.pdf");
    c1.SaveAs("Parabola_ft0_p_value_SeparateChannelwise.png");
}

void makeGraphs_pValue_channelwiseplots_SS()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
   
    double ElEl_p_value1 = 0.0821412;
    double ElEl_p_value2 = 0.0946987;
    double ElEl_p_value3 = 0.108459;
    double ElEl_p_value4 = 0.123416;
    double ElEl_p_value5 = 0.139546;
    double ElEl_p_value6 = 0.156805;
    double ElEl_p_value7 = 0.175127;
    double ElEl_p_value8 = 0.194423;
    double ElEl_p_value9 = 0.214584;
    double ElEl_p_value10 = 0.235483;
    double ElEl_p_value11 = 0.256969;
    double ElEl_p_value12 = 0.278876;
    double ElEl_p_value13 = 0.301019;
    double ElEl_p_value14 = 0.323201;
    double ElEl_p_value15 = 0.345205;
    double ElEl_p_value16 = 0.366806;
    double ElEl_p_value17 = 0.387761;
    double ElEl_p_value18 = 0.407813;
    double ElEl_p_value19 = 0.42669;
    double ElEl_p_value20 = 0.444099;
    double ElEl_p_value21 = 0.459733;
    double ElEl_p_value22 = 0.473264;
    double ElEl_p_value23 = 0.484367;
    double ElEl_p_value24 = 0.492693;
    double ElEl_p_value25 = 0.497973;
    double ElEl_p_value26 = 0.5;
    double ElEl_p_value27 = 0.498702;
    double ElEl_p_value28 = 0.494126;
    double ElEl_p_value29 = 0.486449;
    double ElEl_p_value30 = 0.475945;
    double ElEl_p_value31 = 0.462906;
    double ElEl_p_value32 = 0.447707;
    double ElEl_p_value33 = 0.430661;
    double ElEl_p_value34 = 0.412082;
    double ElEl_p_value35 = 0.392264;
    double ElEl_p_value36 = 0.371486;
    double ElEl_p_value37 = 0.350007;
    double ElEl_p_value38 = 0.328071;
    double ElEl_p_value39 = 0.30591;
    double ElEl_p_value40 = 0.283741;
    double ElEl_p_value41 = 0.261765;
    double ElEl_p_value42 = 0.240171;
    double ElEl_p_value43 = 0.21913;
    double ElEl_p_value44 = 0.198794;
    double ElEl_p_value45 = 0.179297;
    double ElEl_p_value46 = 0.160753;
    double ElEl_p_value47 = 0.143253;
    double ElEl_p_value48 = 0.12687;
    double ElEl_p_value49 = 0.111651;
    double ElEl_p_value50 = 0.0976268;
    double ElEl_p_value51 = 0.0848048;
 
    
    double ElEl_p_value[51] = {ElEl_p_value1, ElEl_p_value2, ElEl_p_value3, ElEl_p_value4, ElEl_p_value5, ElEl_p_value6, ElEl_p_value7, ElEl_p_value8, ElEl_p_value9, ElEl_p_value10, ElEl_p_value11, ElEl_p_value12, ElEl_p_value13, ElEl_p_value14, ElEl_p_value15, ElEl_p_value16, ElEl_p_value17, ElEl_p_value18, ElEl_p_value19, ElEl_p_value20, ElEl_p_value21, ElEl_p_value22, ElEl_p_value23, ElEl_p_value24, ElEl_p_value25, ElEl_p_value26, ElEl_p_value27, ElEl_p_value28, ElEl_p_value29, ElEl_p_value30, ElEl_p_value31, ElEl_p_value32, ElEl_p_value33, ElEl_p_value34, ElEl_p_value35, ElEl_p_value36, ElEl_p_value37, ElEl_p_value38, ElEl_p_value39, ElEl_p_value40, ElEl_p_value41, ElEl_p_value42, ElEl_p_value43, ElEl_p_value44, ElEl_p_value45, ElEl_p_value46, ElEl_p_value47, ElEl_p_value48, ElEl_p_value49, ElEl_p_value50, ElEl_p_value51};
    
    double ElMu_p_value1 = 0.00125994;
    double ElMu_p_value2 = 0.00213009;
    double ElMu_p_value3 = 0.00350014;
    double ElMu_p_value4 = 0.0055926;
    double ElMu_p_value5 = 0.00869374;
    double ElMu_p_value6 = 0.013154;
    double ElMu_p_value7 = 0.0193823;
    double ElMu_p_value8 = 0.0278311;
    double ElMu_p_value9 = 0.0389664;
    double ElMu_p_value10 = 0.0532316;
    double ElMu_p_value11 = 0.071003;
    double ElMu_p_value12 = 0.0925402;
    double ElMu_p_value13 = 0.117941;
    double ElMu_p_value14 = 0.147103;
    double ElMu_p_value15 = 0.179703;
    double ElMu_p_value16 = 0.215192;
    double ElMu_p_value17 = 0.252804;
    double ElMu_p_value18 = 0.291588;
    double ElMu_p_value19 = 0.330447;
    double ElMu_p_value20 = 0.368181;
    double ElMu_p_value21 = 0.403522;
    double ElMu_p_value22 = 0.435165;
    double ElMu_p_value23 = 0.461808;
    double ElMu_p_value24 = 0.482188;
    double ElMu_p_value25 = 0.495205;
    double ElMu_p_value26 = 0.5;
    double ElMu_p_value27 = 0.496304;
    double ElMu_p_value28 = 0.484318;
    double ElMu_p_value29 = 0.464837;
    double ElMu_p_value30 = 0.43894;
    double ElMu_p_value31 = 0.407873;
    double ElMu_p_value32 = 0.372939;
    double ElMu_p_value33 = 0.335444;
    double ElMu_p_value34 = 0.296661;
    double ElMu_p_value35 = 0.257803;
    double ElMu_p_value36 = 0.219981;
    double ElMu_p_value37 = 0.184169;
    double ElMu_p_value38 = 0.151157;
    double ElMu_p_value39 = 0.121525;
    double ElMu_p_value40 = 0.0956245;
    double ElMu_p_value41 = 0.0735867;
    double ElMu_p_value42 = 0.0553373;
    double ElMu_p_value43 = 0.0406354;
    double ElMu_p_value44 = 0.0291177;
    double ElMu_p_value45 = 0.0203456;
    double ElMu_p_value46 = 0.0138543;
    double ElMu_p_value47 = 0.00918865;
    double ElMu_p_value48 = 0.0059321;
    double ElMu_p_value49 = 0.003726;
    double ElMu_p_value50 = 0.00227579;
    double ElMu_p_value51 = 0.00135118;
    
    double ElMu_p_value[51] = {ElMu_p_value1, ElMu_p_value2, ElMu_p_value3, ElMu_p_value4, ElMu_p_value5, ElMu_p_value6, ElMu_p_value7, ElMu_p_value8, ElMu_p_value9, ElMu_p_value10, ElMu_p_value11, ElMu_p_value12, ElMu_p_value13, ElMu_p_value14, ElMu_p_value15, ElMu_p_value16, ElMu_p_value17, ElMu_p_value18, ElMu_p_value19, ElMu_p_value20, ElMu_p_value21, ElMu_p_value22, ElMu_p_value23, ElMu_p_value24, ElMu_p_value25, ElMu_p_value26, ElMu_p_value27, ElMu_p_value28, ElMu_p_value29, ElMu_p_value30, ElMu_p_value31, ElMu_p_value32, ElMu_p_value33, ElMu_p_value34, ElMu_p_value35, ElMu_p_value36, ElMu_p_value37, ElMu_p_value38, ElMu_p_value39, ElMu_p_value40, ElMu_p_value41, ElMu_p_value42, ElMu_p_value43, ElMu_p_value44, ElMu_p_value45, ElMu_p_value46, ElMu_p_value47, ElMu_p_value48, ElMu_p_value49, ElMu_p_value50, ElMu_p_value51};
   
    double MuMu_p_value1 = 0.000311897;
    double MuMu_p_value2 = 0.000583324;
    double MuMu_p_value3 = 0.00105643;
    double MuMu_p_value4 = 0.00185334;
    double MuMu_p_value5 = 0.0031513;
    double MuMu_p_value6 = 0.00519571;
    double MuMu_p_value7 = 0.00831066;
    double MuMu_p_value8 = 0.012905;
    double MuMu_p_value9 = 0.0194646;
    double MuMu_p_value10 = 0.0285362;
    double MuMu_p_value11 = 0.0406923;
    double MuMu_p_value12 = 0.0564832;
    double MuMu_p_value13 = 0.0763762;
    double MuMu_p_value14 = 0.10069;
    double MuMu_p_value15 = 0.129532;
    double MuMu_p_value16 = 0.162746;
    double MuMu_p_value17 = 0.199883;
    double MuMu_p_value18 = 0.240193;
    double MuMu_p_value19 = 0.282634;
    double MuMu_p_value20 = 0.325908;
    double MuMu_p_value21 = 0.368481;
    double MuMu_p_value22 = 0.4086;
    double MuMu_p_value23 = 0.444251;
    double MuMu_p_value24 = 0.473097;
    double MuMu_p_value25 = 0.492512;
    double MuMu_p_value26 = 0.5;
    double MuMu_p_value27 = 0.494391;
    double MuMu_p_value28 = 0.47659;
    double MuMu_p_value29 = 0.448932;
    double MuMu_p_value30 = 0.414118;
    double MuMu_p_value31 = 0.374517;
    double MuMu_p_value32 = 0.332186;
    double MuMu_p_value33 = 0.288914;
    double MuMu_p_value34 = 0.246265;
    double MuMu_p_value35 = 0.205573;
    double MuMu_p_value36 = 0.167921;
    double MuMu_p_value37 = 0.134101;
    double MuMu_p_value38 = 0.104606;
    double MuMu_p_value39 = 0.0796343;
    double MuMu_p_value40 = 0.0591137;
    double MuMu_p_value41 = 0.0427519;
    double MuMu_p_value42 = 0.0300999;
    double MuMu_p_value43 = 0.020615;
    double MuMu_p_value44 = 0.0137249;
    double MuMu_p_value45 = 0.00887665;
    double MuMu_p_value46 = 0.0055736;
    double MuMu_p_value47 = 0.00339556;
    double MuMu_p_value48 = 0.00200604;
    double MuMu_p_value49 = 0.00114873;
    double MuMu_p_value50 = 0.000637248;
    double MuMu_p_value51 = 0.000342341;

    double MuMu_p_value[51] = {MuMu_p_value1, MuMu_p_value2, MuMu_p_value3, MuMu_p_value4, MuMu_p_value5, MuMu_p_value6, MuMu_p_value7, MuMu_p_value8, MuMu_p_value9, MuMu_p_value10, MuMu_p_value11, MuMu_p_value12, MuMu_p_value13, MuMu_p_value14, MuMu_p_value15, MuMu_p_value16, MuMu_p_value17, MuMu_p_value18, MuMu_p_value19, MuMu_p_value20, MuMu_p_value21, MuMu_p_value22, MuMu_p_value23, MuMu_p_value24, MuMu_p_value25, MuMu_p_value26, MuMu_p_value27, MuMu_p_value28, MuMu_p_value29, MuMu_p_value30, MuMu_p_value31, MuMu_p_value32, MuMu_p_value33, MuMu_p_value34, MuMu_p_value35, MuMu_p_value36, MuMu_p_value37, MuMu_p_value38, MuMu_p_value39, MuMu_p_value40, MuMu_p_value41, MuMu_p_value42, MuMu_p_value43, MuMu_p_value44, MuMu_p_value45, MuMu_p_value46, MuMu_p_value47, MuMu_p_value48, MuMu_p_value49, MuMu_p_value50, MuMu_p_value51};

    double ST2000_p_value1 = 1.35132e-05;
    double ST2000_p_value2 = 3.56246e-05;
    double ST2000_p_value3 = 8.89924e-05;
    double ST2000_p_value4 = 0.00021074;
    double ST2000_p_value5 = 0.000473316;
    double ST2000_p_value6 = 0.00100895;
    double ST2000_p_value7 = 0.0020427;
    double ST2000_p_value8 = 0.0039315;
    double ST2000_p_value9 = 0.00720075;
    double ST2000_p_value10 = 0.0125643;
    double ST2000_p_value11 = 0.0209127;
    double ST2000_p_value12 = 0.0332506;
    double ST2000_p_value13 = 0.0505829;
    double ST2000_p_value14 = 0.0737481;
    double ST2000_p_value15 = 0.103235;
    double ST2000_p_value16 = 0.139023;
    double ST2000_p_value17 = 0.180472;
    double ST2000_p_value18 = 0.226316;
    double ST2000_p_value19 = 0.274739;
    double ST2000_p_value20 = 0.323545;
    double ST2000_p_value21 = 0.370347;
    double ST2000_p_value22 = 0.412773;
    double ST2000_p_value23 = 0.448634;
    double ST2000_p_value24 = 0.476034;
    double ST2000_p_value25 = 0.493535;
    double ST2000_p_value26 = 0.5;
    double ST2000_p_value27 = 0.495247;
    double ST2000_p_value28 = 0.479384;
    double ST2000_p_value29 = 0.453424;
    double ST2000_p_value30 = 0.418746;
    double ST2000_p_value31 = 0.377189;
    double ST2000_p_value32 = 0.330907;
    double ST2000_p_value33 = 0.282252;
    double ST2000_p_value34 = 0.23362;
    double ST2000_p_value35 = 0.187252;
    double ST2000_p_value36 = 0.145032;
    double ST2000_p_value37 = 0.108319;
    double ST2000_p_value38 = 0.0778498;
    double ST2000_p_value39 = 0.0537372;
    double ST2000_p_value40 = 0.0355593;
    double ST2000_p_value41 = 0.0225183;
    double ST2000_p_value42 = 0.0136256;
    double ST2000_p_value43 = 0.00786641;
    double ST2000_p_value44 = 0.00432734;
    double ST2000_p_value45 = 0.00226566;
    double ST2000_p_value46 = 0.00112784;
    double ST2000_p_value47 = 0.000533316;
    double ST2000_p_value48 = 0.000239371;
    double ST2000_p_value49 = 0.000101911;
    double ST2000_p_value50 = 4.11333e-05;
    double ST2000_p_value51 = 1.57325e-05;

    double p_value_ST2000[51] = {ST2000_p_value1, ST2000_p_value2, ST2000_p_value3, ST2000_p_value4, ST2000_p_value5, ST2000_p_value6, ST2000_p_value7, ST2000_p_value8, ST2000_p_value9, ST2000_p_value10, ST2000_p_value11, ST2000_p_value12, ST2000_p_value13, ST2000_p_value14, ST2000_p_value15, ST2000_p_value16, ST2000_p_value17, ST2000_p_value18, ST2000_p_value19, ST2000_p_value20, ST2000_p_value21, ST2000_p_value22, ST2000_p_value23, ST2000_p_value24, ST2000_p_value25, ST2000_p_value26, ST2000_p_value27, ST2000_p_value28, ST2000_p_value29, ST2000_p_value30, ST2000_p_value31, ST2000_p_value32, ST2000_p_value33, ST2000_p_value34, ST2000_p_value35, ST2000_p_value36, ST2000_p_value37, ST2000_p_value38, ST2000_p_value39, ST2000_p_value40, ST2000_p_value41, ST2000_p_value42, ST2000_p_value43, ST2000_p_value44, ST2000_p_value45, ST2000_p_value46, ST2000_p_value47, ST2000_p_value48, ST2000_p_value49, ST2000_p_value50, ST2000_p_value51}; 

    double SS_Combined_p_value1 = 1.40007e-06;
    double SS_Combined_p_value2 = 4.43432e-06;
    double SS_Combined_p_value3 = 1.31961e-05;
    double SS_Combined_p_value4 = 3.6936e-05;
    double SS_Combined_p_value5 = 9.72936e-05;
    double SS_Combined_p_value6 = 0.00024118;
    double SS_Combined_p_value7 = 0.000563542;
    double SS_Combined_p_value8 = 0.00124167;
    double SS_Combined_p_value9 = 0.00258269;
    double SS_Combined_p_value10 = 0.00507492;
    double SS_Combined_p_value11 = 0.00944187;
    double SS_Combined_p_value12 = 0.0166445;
    double SS_Combined_p_value13 = 0.027848;
    double SS_Combined_p_value14 = 0.0443196;
    double SS_Combined_p_value15 = 0.0671838;
    double SS_Combined_p_value16 = 0.0972673;
    double SS_Combined_p_value17 = 0.134732;
    double SS_Combined_p_value18 = 0.179042;
    double SS_Combined_p_value19 = 0.228739;
    double SS_Combined_p_value20 = 0.281699;
    double SS_Combined_p_value21 = 0.335105;
    double SS_Combined_p_value22 = 0.385951;
    double SS_Combined_p_value23 = 0.43099;
    double SS_Combined_p_value24 = 0.466923;
    double SS_Combined_p_value25 = 0.490858;
    double SS_Combined_p_value26 = 0.5;
    double SS_Combined_p_value27 = 0.49312;
    double SS_Combined_p_value28 = 0.47124;
    double SS_Combined_p_value29 = 0.436799;
    double SS_Combined_p_value30 = 0.392843;
    double SS_Combined_p_value31 = 0.342689;
    double SS_Combined_p_value32 = 0.289399;
    double SS_Combined_p_value33 = 0.23625;
    double SS_Combined_p_value34 = 0.185955;
    double SS_Combined_p_value35 = 0.140759;
    double SS_Combined_p_value36 = 0.102253;
    double SS_Combined_p_value37 = 0.0710818;
    double SS_Combined_p_value38 = 0.0472097;
    double SS_Combined_p_value39 = 0.0298824;
    double SS_Combined_p_value40 = 0.0179922;
    double SS_Combined_p_value41 = 0.0102878;
    double SS_Combined_p_value42 = 0.00557264;
    double SS_Combined_p_value43 = 0.00285888;
    double SS_Combined_p_value44 = 0.00138595;
    double SS_Combined_p_value45 = 0.000634553;
    double SS_Combined_p_value46 = 0.000273866;
    double SS_Combined_p_value47 = 0.000111439;
    double SS_Combined_p_value48 = 4.26984e-05;
    double SS_Combined_p_value49 = 1.53904e-05;
    double SS_Combined_p_value50 = 5.21724e-06;
    double SS_Combined_p_value51 = 1.66275e-06;

    double SS_Combined_p_value[51] = {SS_Combined_p_value1, SS_Combined_p_value2, SS_Combined_p_value3, SS_Combined_p_value4, SS_Combined_p_value5, SS_Combined_p_value6, SS_Combined_p_value7, SS_Combined_p_value8, SS_Combined_p_value9, SS_Combined_p_value10, SS_Combined_p_value11, SS_Combined_p_value12, SS_Combined_p_value13, SS_Combined_p_value14, SS_Combined_p_value15, SS_Combined_p_value16, SS_Combined_p_value17, SS_Combined_p_value18, SS_Combined_p_value19, SS_Combined_p_value20, SS_Combined_p_value21, SS_Combined_p_value22, SS_Combined_p_value23, SS_Combined_p_value24, SS_Combined_p_value25, SS_Combined_p_value26, SS_Combined_p_value27, SS_Combined_p_value28, SS_Combined_p_value29, SS_Combined_p_value30, SS_Combined_p_value31, SS_Combined_p_value32, SS_Combined_p_value33, SS_Combined_p_value34, SS_Combined_p_value35, SS_Combined_p_value36, SS_Combined_p_value37, SS_Combined_p_value38, SS_Combined_p_value39, SS_Combined_p_value40, SS_Combined_p_value41, SS_Combined_p_value42, SS_Combined_p_value43, SS_Combined_p_value44, SS_Combined_p_value45, SS_Combined_p_value46, SS_Combined_p_value47, SS_Combined_p_value48, SS_Combined_p_value49, SS_Combined_p_value50, SS_Combined_p_value51};
 
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(51, param, ElEl_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(51, param, ElMu_p_value);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
   
    TGraph *gr_2 = new TGraph(51, param, MuMu_p_value);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");

    TGraph *gr_3 = new TGraph(51, param, p_value_ST2000);
    gr_3->SetLineColor(kCyan);
    //gr_2->SetLineWidth(3);
    gr_3->Draw("LP*");

    TGraph *gr_4 = new TGraph(51, param, SS_Combined_p_value);
    gr_4->SetLineColor(kOrange);
    //gr_2->SetLineWidth(3);
    gr_4->Draw("LP*");
 
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.40,0.30,0.60,0.50,NULL,"brNDC");
    leg->AddEntry(gr,   "Dilepton: ElEl", "pl");
    leg->AddEntry(gr_1,   "Dilepton: ElMu", "pl");
    leg->AddEntry(gr_2,   "Dilepton: MuMu", "pl");
    leg->AddEntry(gr_3,   "Dilepton: Summed", "pl");
    leg->AddEntry(gr_4,   "Dilepton: 3 bins", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined ElEl = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined ElEl = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined ElMu = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined ElMu = " << m1 << endl;
    
    double m2=-2.0;
    while (m2<=2.0 and gr_2->Eval(m2)<gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected lower limit combined MuMu = " << m2 << endl;
    
    while (m2<=2.0 and gr_2->Eval(m2)>gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected upper limit combined MuMu = " << m2 << endl;

    double m3=-2.0;
    while (m3<=2.0 and gr_3->Eval(m3)<gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected lower limit combined = " << m3 << endl;
    
    while (m3<=2.0 and gr_3->Eval(m3)>gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected upper limit combined = " << m3 << endl;

    double m4=-2.0;
    while (m4<=2.0 and gr_4->Eval(m4)<gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected lower limit combined = " << m4 << endl;

    while (m4<=2.0 and gr_4->Eval(m4)>gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected upper limit combined = " << m4 << endl;

    c1.SaveAs("Parabola_ft0_p_value_Channelwise.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Channelwise.png");
}

void makeGraphs_pValue_channelwiseplots_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double ZeroOSSF_p_value1 = 0.000386417;
    double ZeroOSSF_p_value2 = 0.00068111;
    double ZeroOSSF_p_value3 = 0.00116893;
    double ZeroOSSF_p_value4 = 0.00195411;
    double ZeroOSSF_p_value5 = 0.00318303;
    double ZeroOSSF_p_value6 = 0.00505382;
    double ZeroOSSF_p_value7 = 0.00782458;
    double ZeroOSSF_p_value8 = 0.011818;
    double ZeroOSSF_p_value9 = 0.0174202;
    double ZeroOSSF_p_value10 = 0.0250727;
    double ZeroOSSF_p_value11 = 0.0352532;
    double ZeroOSSF_p_value12 = 0.0484476;
    double ZeroOSSF_p_value13 = 0.0651127;
    double ZeroOSSF_p_value14 = 0.0856308;
    double ZeroOSSF_p_value15 = 0.110262;
    double ZeroOSSF_p_value16 = 0.139102;
    double ZeroOSSF_p_value17 = 0.172039;
    double ZeroOSSF_p_value18 = 0.208735;
    double ZeroOSSF_p_value19 = 0.248611;
    double ZeroOSSF_p_value20 = 0.29084;
    double ZeroOSSF_p_value21 = 0.334366;
    double ZeroOSSF_p_value22 = 0.377868;
    double ZeroOSSF_p_value23 = 0.419673;
    double ZeroOSSF_p_value24 = 0.457415;
    double ZeroOSSF_p_value25 = 0.486973;
    double ZeroOSSF_p_value26 = 0.5;
    double ZeroOSSF_p_value27 = 0.488883;
    double ZeroOSSF_p_value28 = 0.460331;
    double ZeroOSSF_p_value29 = 0.423079;
    double ZeroOSSF_p_value30 = 0.381515;
    double ZeroOSSF_p_value31 = 0.338087;
    double ZeroOSSF_p_value32 = 0.294508;
    double ZeroOSSF_p_value33 = 0.252123;
    double ZeroOSSF_p_value34 = 0.212012;
    double ZeroOSSF_p_value35 = 0.175019;
    double ZeroOSSF_p_value36 = 0.141746;
    double ZeroOSSF_p_value37 = 0.112551;
    double ZeroOSSF_p_value38 = 0.0875626;
    double ZeroOSSF_p_value39 = 0.066703;
    double ZeroOSSF_p_value40 = 0.0497238;
    double ZeroOSSF_p_value41 = 0.0362512;
    double ZeroOSSF_p_value42 = 0.0258332;
    double ZeroOSSF_p_value43 = 0.0179847;
    double ZeroOSSF_p_value44 = 0.0122259;
    double ZeroOSSF_p_value45 = 0.00811159;
    double ZeroOSSF_p_value46 = 0.00525035;
    double ZeroOSSF_p_value47 = 0.00331394;
    double ZeroOSSF_p_value48 = 0.00203894;
    double ZeroOSSF_p_value49 = 0.00122239;
    double ZeroOSSF_p_value50 = 0.00071383;
    double ZeroOSSF_p_value51 = 0.000405919;    
    
    double ZeroOSSF_p_value[51] = {ZeroOSSF_p_value1, ZeroOSSF_p_value2, ZeroOSSF_p_value3, ZeroOSSF_p_value4, ZeroOSSF_p_value5, ZeroOSSF_p_value6, ZeroOSSF_p_value7, ZeroOSSF_p_value8, ZeroOSSF_p_value9, ZeroOSSF_p_value10, ZeroOSSF_p_value11, ZeroOSSF_p_value12, ZeroOSSF_p_value13, ZeroOSSF_p_value14, ZeroOSSF_p_value15, ZeroOSSF_p_value16, ZeroOSSF_p_value17, ZeroOSSF_p_value18, ZeroOSSF_p_value19, ZeroOSSF_p_value20, ZeroOSSF_p_value21, ZeroOSSF_p_value22, ZeroOSSF_p_value23, ZeroOSSF_p_value24, ZeroOSSF_p_value25, ZeroOSSF_p_value26, ZeroOSSF_p_value27, ZeroOSSF_p_value28, ZeroOSSF_p_value29, ZeroOSSF_p_value30, ZeroOSSF_p_value31, ZeroOSSF_p_value32, ZeroOSSF_p_value33, ZeroOSSF_p_value34, ZeroOSSF_p_value35, ZeroOSSF_p_value36, ZeroOSSF_p_value37, ZeroOSSF_p_value38, ZeroOSSF_p_value39, ZeroOSSF_p_value40, ZeroOSSF_p_value41, ZeroOSSF_p_value42, ZeroOSSF_p_value43, ZeroOSSF_p_value44, ZeroOSSF_p_value45, ZeroOSSF_p_value46, ZeroOSSF_p_value47, ZeroOSSF_p_value48, ZeroOSSF_p_value49, ZeroOSSF_p_value50, ZeroOSSF_p_value51};
    
    double OneOSSF_p_value1 = 1.74308e-05;
    double OneOSSF_p_value2 = 4.03969e-05;
    double OneOSSF_p_value3 = 8.98668e-05;
    double OneOSSF_p_value4 = 0.000192015;
    double OneOSSF_p_value5 = 0.000394194;
    double OneOSSF_p_value6 = 0.000777979;
    double OneOSSF_p_value7 = 0.00147681;
    double OneOSSF_p_value8 = 0.00269816;
    double OneOSSF_p_value9 = 0.0047474;
    double OneOSSF_p_value10 = 0.00804998;
    double OneOSSF_p_value11 = 0.0131647;
    double OneOSSF_p_value12 = 0.0207801;
    double OneOSSF_p_value13 = 0.0316873;
    double OneOSSF_p_value14 = 0.0467217;
    double OneOSSF_p_value15 = 0.0666768;
    double OneOSSF_p_value16 = 0.0921952;
    double OneOSSF_p_value17 = 0.12365;
    double OneOSSF_p_value18 = 0.161036;
    double OneOSSF_p_value19 = 0.203895;
    double OneOSSF_p_value20 = 0.25127;
    double OneOSSF_p_value21 = 0.301707;
    double OneOSSF_p_value22 = 0.353259;
    double OneOSSF_p_value23 = 0.403433;
    double OneOSSF_p_value24 = 0.448875;
    double OneOSSF_p_value25 = 0.484225;
    double OneOSSF_p_value26 = 0.5;
    double OneOSSF_p_value27 = 0.488011;
    double OneOSSF_p_value28 = 0.454796;
    double OneOSSF_p_value29 = 0.410415;
    double OneOSSF_p_value30 = 0.36069;
    double OneOSSF_p_value31 = 0.309167;
    double OneOSSF_p_value32 = 0.258437;
    double OneOSSF_p_value33 = 0.210518;
    double OneOSSF_p_value34 = 0.166932;
    double OneOSSF_p_value35 = 0.128711;
    double OneOSSF_p_value36 = 0.0963855;
    double OneOSSF_p_value37 = 0.0700214;
    double OneOSSF_p_value38 = 0.0492941;
    double OneOSSF_p_value39 = 0.0335928;
    double OneOSSF_p_value40 = 0.0221388;
    double OneOSSF_p_value41 = 0.0140967;
    double OneOSSF_p_value42 = 0.00866471;
    double OneOSSF_p_value43 = 0.00513707;
    double OneOSSF_p_value44 = 0.00293545;
    double OneOSSF_p_value45 = 0.00161557;
    double OneOSSF_p_value46 = 0.000855816;
    double OneOSSF_p_value47 = 0.000436107;   
    double OneOSSF_p_value48 = 0.00021366;
    double OneOSSF_p_value49 = 0.000100582;
    double OneOSSF_p_value50 = 4.54787e-05;
    double OneOSSF_p_value51 = 1.97421e-05;
 
    double OneOSSF_p_value[51] = {OneOSSF_p_value1, OneOSSF_p_value2, OneOSSF_p_value3, OneOSSF_p_value4, OneOSSF_p_value5, OneOSSF_p_value6, OneOSSF_p_value7, OneOSSF_p_value8, OneOSSF_p_value9, OneOSSF_p_value10, OneOSSF_p_value11, OneOSSF_p_value12, OneOSSF_p_value13, OneOSSF_p_value14, OneOSSF_p_value15, OneOSSF_p_value16, OneOSSF_p_value17, OneOSSF_p_value18, OneOSSF_p_value19, OneOSSF_p_value20, OneOSSF_p_value21, OneOSSF_p_value22, OneOSSF_p_value23, OneOSSF_p_value24, OneOSSF_p_value25, OneOSSF_p_value26, OneOSSF_p_value27, OneOSSF_p_value28, OneOSSF_p_value29, OneOSSF_p_value30, OneOSSF_p_value31, OneOSSF_p_value32, OneOSSF_p_value33, OneOSSF_p_value34, OneOSSF_p_value35, OneOSSF_p_value36, OneOSSF_p_value37, OneOSSF_p_value38, OneOSSF_p_value39, OneOSSF_p_value40, OneOSSF_p_value41, OneOSSF_p_value42, OneOSSF_p_value43, OneOSSF_p_value44, OneOSSF_p_value45, OneOSSF_p_value46, OneOSSF_p_value47, OneOSSF_p_value48, OneOSSF_p_value49, OneOSSF_p_value50, OneOSSF_p_value51};
    
    double TwoOSSF_p_value1 = 5.81205e-05;
    double TwoOSSF_p_value2 = 0.000120464;
    double TwoOSSF_p_value3 = 0.000241147;
    double TwoOSSF_p_value4 = 0.00046642;
    double TwoOSSF_p_value5 = 0.000872007;
    double TwoOSSF_p_value6 = 0.00157653;
    double TwoOSSF_p_value7 = 0.00275766;
    double TwoOSSF_p_value8 = 0.00466934;
    double TwoOSSF_p_value9 = 0.00765739;
    double TwoOSSF_p_value10 = 0.0121697;
    double TwoOSSF_p_value11 = 0.018755;
    double TwoOSSF_p_value12 = 0.0280475;
    double TwoOSSF_p_value13 = 0.0407302;
    double TwoOSSF_p_value14 = 0.0574802;
    double TwoOSSF_p_value15 = 0.0788956;
    double TwoOSSF_p_value16 = 0.105412;
    double TwoOSSF_p_value17 = 0.13722;
    double TwoOSSF_p_value18 = 0.174195;
    double TwoOSSF_p_value19 = 0.215848;
    double TwoOSSF_p_value20 = 0.261306;
    double TwoOSSF_p_value21 = 0.309317;
    double TwoOSSF_p_value22 = 0.358235;
    double TwoOSSF_p_value23 = 0.405971;
    double TwoOSSF_p_value24 = 0.449634;
    double TwoOSSF_p_value25 = 0.484342;
    double TwoOSSF_p_value26 = 0.5;
    double TwoOSSF_p_value27 = 0.48686;
    double TwoOSSF_p_value28 = 0.453416;
    double TwoOSSF_p_value29 = 0.41034;
    double TwoOSSF_p_value30 = 0.362851;
    double TwoOSSF_p_value31 = 0.313944;
    double TwoOSSF_p_value32 = 0.265773;
    double TwoOSSF_p_value33 = 0.220013;
    double TwoOSSF_p_value34 = 0.177954;
    double TwoOSSF_p_value35 = 0.140509;
    double TwoOSSF_p_value36 = 0.108201;
    double TwoOSSF_p_value37 = 0.0811863;
    double TwoOSSF_p_value38 = 0.0593025;
    double TwoOSSF_p_value39 = 0.0421335;
    double TwoOSSF_p_value40 = 0.0290936;
    double TwoOSSF_p_value41 = 0.0195095;
    double TwoOSSF_p_value42 = 0.0126957;
    double TwoOSSF_p_value43 = 0.00801199;
    double TwoOSSF_p_value44 = 0.00490024;
    double TwoOSSF_p_value45 = 0.00290288;
    double TwoOSSF_p_value46 = 0.00166472;
    double TwoOSSF_p_value47 = 0.000923654;
    double TwoOSSF_p_value48 = 0.000495629;
    double TwoOSSF_p_value49 = 0.000257081;
    double TwoOSSF_p_value50 = 0.000128845;
    double TwoOSSF_p_value51 = 6.23709e-05;   

    double TwoOSSF_p_value[51] = {TwoOSSF_p_value1, TwoOSSF_p_value2, TwoOSSF_p_value3, TwoOSSF_p_value4, TwoOSSF_p_value5, TwoOSSF_p_value6, TwoOSSF_p_value7, TwoOSSF_p_value8, TwoOSSF_p_value9, TwoOSSF_p_value10, TwoOSSF_p_value11, TwoOSSF_p_value12, TwoOSSF_p_value13, TwoOSSF_p_value14, TwoOSSF_p_value15, TwoOSSF_p_value16, TwoOSSF_p_value17, TwoOSSF_p_value18, TwoOSSF_p_value19, TwoOSSF_p_value20, TwoOSSF_p_value21, TwoOSSF_p_value22, TwoOSSF_p_value23, TwoOSSF_p_value24, TwoOSSF_p_value25, TwoOSSF_p_value26, TwoOSSF_p_value27, TwoOSSF_p_value28, TwoOSSF_p_value29, TwoOSSF_p_value30, TwoOSSF_p_value31, TwoOSSF_p_value32, TwoOSSF_p_value33, TwoOSSF_p_value34, TwoOSSF_p_value35, TwoOSSF_p_value36, TwoOSSF_p_value37, TwoOSSF_p_value38, TwoOSSF_p_value39, TwoOSSF_p_value40, TwoOSSF_p_value41, TwoOSSF_p_value42, TwoOSSF_p_value43, TwoOSSF_p_value44, TwoOSSF_p_value45, TwoOSSF_p_value46, TwoOSSF_p_value47, TwoOSSF_p_value48, TwoOSSF_p_value49, TwoOSSF_p_value50, TwoOSSF_p_value51};
   
    double ST1500_p_value1 = 4.00434e-11;
    double ST1500_p_value2 = 2.89704e-10;
    double ST1500_p_value3 = 1.90293e-09;
    double ST1500_p_value4 = 1.13526e-08;
    double ST1500_p_value5 = 6.15367e-08;
    double ST1500_p_value6 = 3.03298e-07;
    double ST1500_p_value7 = 1.3601e-06;
    double ST1500_p_value8 = 5.55436e-06;
    double ST1500_p_value9 = 2.06777e-05;
    double ST1500_p_value10 = 7.02526e-05;
    double ST1500_p_value11 = 0.000218133;
    double ST1500_p_value12 = 0.000619906;
    double ST1500_p_value13 = 0.0016153;
    double ST1500_p_value14 = 0.0038669;
    double ST1500_p_value15 = 0.00852296;
    double ST1500_p_value16 = 0.0173394;
    double ST1500_p_value17 = 0.0326507;
    double ST1500_p_value18 = 0.0570807;
    double ST1500_p_value19 = 0.0929579;
    double ST1500_p_value20 = 0.141534;
    double ST1500_p_value21 = 0.202242;
    double ST1500_p_value22 = 0.272258;
    double ST1500_p_value23 = 0.346448;
    double ST1500_p_value24 = 0.41741;
    double ST1500_p_value25 = 0.474421;
    double ST1500_p_value26 = 0.5;
    double ST1500_p_value27 = 0.47934;
    double ST1500_p_value28 = 0.42488;
    double ST1500_p_value29 = 0.354882;
    double ST1500_p_value30 = 0.28064;
    double ST1500_p_value31 = 0.209845;
    double ST1500_p_value32 = 0.147887;
    double ST1500_p_value33 = 0.0978567;
    double ST1500_p_value34 = 0.0605638;
    double ST1500_p_value35 = 0.034931;
    double ST1500_p_value36 = 0.0187115;
    double ST1500_p_value37 = 0.00928034;
    double ST1500_p_value38 = 0.00424965;
    double ST1500_p_value39 = 0.00179231;
    double ST1500_p_value40 = 0.000694601;
    double ST1500_p_value41 = 0.000246861;
    double ST1500_p_value42 = 8.03233e-05;
    double ST1500_p_value43 = 2.3888e-05;
    double ST1500_p_value44 = 6.48441e-06;
    double ST1500_p_value45 = 1.6048e-06;
    double ST1500_p_value46 = 3.61729e-07;
    double ST1500_p_value47 = 7.41948e-08;
    double ST1500_p_value48 = 1.38371e-08;
    double ST1500_p_value49 = 2.34496e-09;
    double ST1500_p_value50 = 3.60972e-10;
    double ST1500_p_value51 = 5.04486e-11;

    double p_value_ST1500[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};

    double Tri_Combined_p_value1 = 3.73145e-11;
    double Tri_Combined_p_value2 = 2.71526e-10;
    double Tri_Combined_p_value3 = 1.79212e-09;
    double Tri_Combined_p_value4 = 1.07446e-08;
    double Tri_Combined_p_value5 = 5.85451e-08;
    double Tri_Combined_p_value6 = 2.89765e-07;
    double Tri_Combined_p_value7 = 1.30514e-06;
    double Tri_Combined_p_value8 = 5.35129e-06;
    double Tri_Combined_p_value9 = 1.99934e-05;
    double Tri_Combined_p_value10 = 6.81868e-05;
    double Tri_Combined_p_value11 = 0.000212439;
    double Tri_Combined_p_value12 = 0.000605868;
    double Tri_Combined_p_value13 = 0.00158294;
    double Tri_Combined_p_value14 = 0.00380184;
    double Tri_Combined_p_value15 = 0.0084001;
    double Tri_Combined_p_value16 = 0.017125;
    double Tri_Combined_p_value17 = 0.0323027;
    double Tri_Combined_p_value18 = 0.0565771;
    double Tri_Combined_p_value19 = 0.0923185;
    double Tri_Combined_p_value20 = 0.140789;
    double Tri_Combined_p_value21 = 0.201352;
    double Tri_Combined_p_value22 = 0.271577;
    double Tri_Combined_p_value23 = 0.345665;
    double Tri_Combined_p_value24 = 0.417059;
    double Tri_Combined_p_value25 = 0.474349;
    double Tri_Combined_p_value26 = 0.5;
    double Tri_Combined_p_value27 = 0.479134;
    double Tri_Combined_p_value28 = 0.424286;
    double Tri_Combined_p_value29 = 0.354001;
    double Tri_Combined_p_value30 = 0.27959;
    double Tri_Combined_p_value31 = 0.208917;
    double Tri_Combined_p_value32 = 0.147018;
    double Tri_Combined_p_value33 = 0.0971104;
    double Tri_Combined_p_value34 = 0.0599737;
    double Tri_Combined_p_value35 = 0.0345217;
    double Tri_Combined_p_value36 = 0.0184655;
    double Tri_Combined_p_value37 = 0.00912906;
    double Tri_Combined_p_value38 = 0.00416956;
    double Tri_Combined_p_value39 = 0.00175321;
    double Tri_Combined_p_value40 = 0.000677845;
    double Tri_Combined_p_value41 = 0.00024014;
    double Tri_Combined_p_value42 = 7.78138e-05;
    double Tri_Combined_p_value43 = 2.30617e-05;
    double Tri_Combined_p_value44 = 6.23638e-06;
    double Tri_Combined_p_value45 = 1.53612e-06;
    double Tri_Combined_p_value46 = 3.44656e-07;
    double Tri_Combined_p_value47 = 7.03791e-08;
    double Tri_Combined_p_value48 = 1.30625e-08;
    double Tri_Combined_p_value49 = 2.20106e-09;
    double Tri_Combined_p_value50 = 3.37111e-10;
    double Tri_Combined_p_value51 = 4.68588e-11;

    double Tri_Combined_p_value[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};

    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(51, param, ZeroOSSF_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(51, param, OneOSSF_p_value);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(51, param, TwoOSSF_p_value);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");
    
    TGraph *gr_3 = new TGraph(51, param, p_value_ST1500);
    gr_3->SetLineColor(kCyan);
    //gr_2->SetLineWidth(3);
    gr_3->Draw("LP*");
    
    TGraph *gr_4 = new TGraph(51, param, Tri_Combined_p_value);
    gr_4->SetLineColor(kOrange);
    gr_4->Draw("LP*");

    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
   
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: 0OSSF", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 1OSSF", "pl");
    leg->AddEntry(gr_2,   "Trilepton: 2OSSF", "pl");
    leg->AddEntry(gr_3,   "Trilepton: All", "pl");
    leg->AddEntry(gr_4,   "Trilepton: 3 bins", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 0OSSF = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 0OSSF = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 1OSSF = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 1OSSF = " << m1 << endl;

    double m2=-2.0;
    while (m2<=2.0 and gr_2->Eval(m2)<gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected lower limit combined 2SFOS = " << m2 << endl;

    while (m2<=2.0 and gr_2->Eval(m2)>gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected upper limit combined 2SFOS = " << m2 << endl;

    double m3=-2.0;
    while (m3<=2.0 and gr_3->Eval(m3)<gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected lower limit combined channel = " << m3 << endl;

    while (m3<=2.0 and gr_3->Eval(m3)>gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected upper limit combined channel = " << m3 << endl;
   
    double m4=-2.0;
    while (m4<=2.0 and gr_4->Eval(m4)<gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected lower limit combined channel 3 bins = " << m4 << endl;

    while (m4<=2.0 and gr_4->Eval(m4)>gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected upper limit combined channel 3 bins = " << m4 << endl;
 
    c1.SaveAs("Parabola_ft0_p_value_Channelwise_3l.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Channelwise_3l.png");
}

void makeGraphs_pValue_All_withSignalUnc()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double noST_p_value1 = 0.134554;
    double noST_p_value2 = 0.15191;
    double noST_p_value3 = 0.170217;
    double noST_p_value4 = 0.18934;
    double noST_p_value5 = 0.209124;
    double noST_p_value6 = 0.229401;
    double noST_p_value7 = 0.249992;
    double noST_p_value8 = 0.270709;
    double noST_p_value9 = 0.291366;
    double noST_p_value10 = 0.311782;
    double noST_p_value11 = 0.331771;
    double noST_p_value12 = 0.35117;
    double noST_p_value13 = 0.369823;
    double noST_p_value14 = 0.387588;
    double noST_p_value15 = 0.404341;
    double noST_p_value16 = 0.419969;
    double noST_p_value17 = 0.43438;
    double noST_p_value18 = 0.447493;
    double noST_p_value19 = 0.459242;
    double noST_p_value20 = 0.469572;
    double noST_p_value21 = 0.478443;
    double noST_p_value22 = 0.485821;
    double noST_p_value23 = 0.491655;
    double noST_p_value24 = 0.496442;
    double noST_p_value25 = 0.498914;
    double noST_p_value26 = 0.5;
    double noST_p_value27 = 0.499698;
    double noST_p_value28 = 0.498008;
    double noST_p_value29 = 0.494363;
    double noST_p_value30 = 0.489284;
    double noST_p_value31 = 0.482795;
    double noST_p_value32 = 0.474764;
    double noST_p_value33 = 0.465253;
    double noST_p_value34 = 0.454299;
    double noST_p_value35 = 0.441949;
    double noST_p_value36 = 0.428263;
    double noST_p_value37 = 0.413312;
    double noST_p_value38 = 0.397181;
    double noST_p_value39 = 0.379973;
    double noST_p_value40 = 0.361805;
    double noST_p_value41 = 0.34281;
    double noST_p_value42 = 0.323134;
    double noST_p_value43 = 0.302939;
    double noST_p_value44 = 0.282395;
    double noST_p_value45 = 0.261688;
    double noST_p_value46 = 0.241003;
    double noST_p_value47 = 0.220526;
    double noST_p_value48 = 0.200443;
    double noST_p_value49 = 0.180926;
    double noST_p_value50 = 0.162141;
    double noST_p_value51 = 0.144232;
    
    double p_value_noST[51] = {noST_p_value1, noST_p_value2, noST_p_value3, noST_p_value4, noST_p_value5, noST_p_value6, noST_p_value7, noST_p_value8, noST_p_value9, noST_p_value10, noST_p_value11, noST_p_value12, noST_p_value13, noST_p_value14, noST_p_value15, noST_p_value16, noST_p_value17, noST_p_value18, noST_p_value19, noST_p_value20, noST_p_value21, noST_p_value22, noST_p_value23, noST_p_value24, noST_p_value25, noST_p_value26, noST_p_value27, noST_p_value28, noST_p_value29, noST_p_value30, noST_p_value31, noST_p_value32, noST_p_value33, noST_p_value34, noST_p_value35, noST_p_value36, noST_p_value37, noST_p_value38, noST_p_value39, noST_p_value40, noST_p_value41, noST_p_value42, noST_p_value43, noST_p_value44, noST_p_value45, noST_p_value46, noST_p_value47, noST_p_value48, noST_p_value49, noST_p_value50, noST_p_value51};
    
    double ST1000_p_value1 = 4.24651e-07;
    double ST1000_p_value2 = 1.54217e-06;
    double ST1000_p_value3 = 5.22373e-06;
    double ST1000_p_value4 = 1.6506e-05;
    double ST1000_p_value5 = 4.86661e-05;
    double ST1000_p_value6 = 0.000133939;
    double ST1000_p_value7 = 0.000344316;
    double ST1000_p_value8 = 0.000827447;
    double ST1000_p_value9 = 0.00186068;
    double ST1000_p_value10 = 0.00392012;
    double ST1000_p_value11 = 0.00774891;
    double ST1000_p_value12 = 0.0143952;
    double ST1000_p_value13 = 0.0251803;
    double ST1000_p_value14 = 0.0415623;
    double ST1000_p_value15 = 0.0648885;
    double ST1000_p_value16 = 0.0960749;
    double ST1000_p_value17 = 0.135281;
    double ST1000_p_value18 = 0.181702;
    double ST1000_p_value19 = 0.23352;
    double ST1000_p_value20 = 0.288077;
    double ST1000_p_value21 = 0.342178;
    double ST1000_p_value22 = 0.392472;
    double ST1000_p_value23 = 0.435804;
    double ST1000_p_value24 = 0.469456;
    double ST1000_p_value25 = 0.491316;
    double ST1000_p_value26 = 0.5;
    double ST1000_p_value27 = 0.495003;
    double ST1000_p_value28 = 0.476497;
    double ST1000_p_value29 = 0.445805;
    double ST1000_p_value30 = 0.404806;
    double ST1000_p_value31 = 0.356071;
    double ST1000_p_value32 = 0.302659;
    double ST1000_p_value33 = 0.247898;
    double ST1000_p_value34 = 0.195063;
    double ST1000_p_value35 = 0.146989;
    double ST1000_p_value36 = 0.105741;
    double ST1000_p_value37 = 0.0723995;
    double ST1000_p_value38 = 0.0470446;
    double ST1000_p_value39 = 0.0289339;
    double ST1000_p_value40 = 0.0168021;
    double ST1000_p_value41 = 0.0091925;
    double ST1000_p_value42 = 0.00472869;
    double ST1000_p_value43 = 0.00228325;
    double ST1000_p_value44 = 0.00103323;
    double ST1000_p_value45 = 0.000437669;
    double ST1000_p_value46 = 0.000173349;
    double ST1000_p_value47 = 6.41436e-05;
    double ST1000_p_value48 = 2.21592e-05;
    double ST1000_p_value49 = 7.14367e-06;
    double ST1000_p_value50 = 2.14834e-06;
    double ST1000_p_value51 = 6.02642e-07;
    
    double p_value_ST1000[51] = {ST1000_p_value1, ST1000_p_value2, ST1000_p_value3, ST1000_p_value4, ST1000_p_value5, ST1000_p_value6, ST1000_p_value7, ST1000_p_value8, ST1000_p_value9, ST1000_p_value10, ST1000_p_value11, ST1000_p_value12, ST1000_p_value13, ST1000_p_value14, ST1000_p_value15, ST1000_p_value16, ST1000_p_value17, ST1000_p_value18, ST1000_p_value19, ST1000_p_value20, ST1000_p_value21, ST1000_p_value22, ST1000_p_value23, ST1000_p_value24, ST1000_p_value25, ST1000_p_value26, ST1000_p_value27, ST1000_p_value28, ST1000_p_value29, ST1000_p_value30, ST1000_p_value31, ST1000_p_value32, ST1000_p_value33, ST1000_p_value34, ST1000_p_value35, ST1000_p_value36, ST1000_p_value37, ST1000_p_value38, ST1000_p_value39, ST1000_p_value40, ST1000_p_value41, ST1000_p_value42, ST1000_p_value43, ST1000_p_value44, ST1000_p_value45, ST1000_p_value46, ST1000_p_value47, ST1000_p_value48, ST1000_p_value49, ST1000_p_value50, ST1000_p_value51};
    
    double ST1500_p_value1 = 4.00434e-11;
    double ST1500_p_value2 = 2.89704e-10;
    double ST1500_p_value3 = 1.90293e-09;
    double ST1500_p_value4 = 1.13526e-08;
    double ST1500_p_value5 = 6.15367e-08;
    double ST1500_p_value6 = 3.03298e-07;
    double ST1500_p_value7 = 1.3601e-06;
    double ST1500_p_value8 = 5.55436e-06;
    double ST1500_p_value9 = 2.06777e-05;
    double ST1500_p_value10 = 7.02526e-05;
    double ST1500_p_value11 = 0.000218133;
    double ST1500_p_value12 = 0.000619906;
    double ST1500_p_value13 = 0.0016153;
    double ST1500_p_value14 = 0.0038669;
    double ST1500_p_value15 = 0.00852296;
    double ST1500_p_value16 = 0.0173394;
    double ST1500_p_value17 = 0.0326507;
    double ST1500_p_value18 = 0.0570807;
    double ST1500_p_value19 = 0.0929579;
    double ST1500_p_value20 = 0.141534;
    double ST1500_p_value21 = 0.202242;
    double ST1500_p_value22 = 0.272258;
    double ST1500_p_value23 = 0.346448;
    double ST1500_p_value24 = 0.41741;
    double ST1500_p_value25 = 0.474421;
    double ST1500_p_value26 = 0.5;
    double ST1500_p_value27 = 0.47934;
    double ST1500_p_value28 = 0.42488;
    double ST1500_p_value29 = 0.354882;
    double ST1500_p_value30 = 0.28064;
    double ST1500_p_value31 = 0.209845;
    double ST1500_p_value32 = 0.147887;
    double ST1500_p_value33 = 0.0978567;
    double ST1500_p_value34 = 0.0605638;
    double ST1500_p_value35 = 0.034931;
    double ST1500_p_value36 = 0.0187115;
    double ST1500_p_value37 = 0.00928034;
    double ST1500_p_value38 = 0.00424965;
    double ST1500_p_value39 = 0.00179231;
    double ST1500_p_value40 = 0.000694601;
    double ST1500_p_value41 = 0.000246861;
    double ST1500_p_value42 = 8.03233e-05;
    double ST1500_p_value43 = 2.3888e-05;
    double ST1500_p_value44 = 6.48441e-06;
    double ST1500_p_value45 = 1.6048e-06;
    double ST1500_p_value46 = 3.61729e-07;
    double ST1500_p_value47 = 7.41948e-08;
    double ST1500_p_value48 = 1.38371e-08;
    double ST1500_p_value49 = 2.34496e-09;
    double ST1500_p_value50 = 3.60972e-10;
    double ST1500_p_value51 = 5.04486e-11;
    
    double p_value_ST1500[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};
    
    double ST2000_p_value1 = 5.26354e-09;
    double ST2000_p_value2 = 2.40529e-08;
    double ST2000_p_value3 = 1.02188e-07;
    double ST2000_p_value4 = 4.03794e-07;
    double ST2000_p_value5 = 1.4849e-06;
    double ST2000_p_value6 = 5.08491e-06;
    double ST2000_p_value7 = 1.62277e-05;
    double ST2000_p_value8 = 4.83017e-05;
    double ST2000_p_value9 = 0.0001342;
    double ST2000_p_value10 = 0.000348381;
    double ST2000_p_value11 = 0.000845935;
    double ST2000_p_value12 = 0.00192378;
    double ST2000_p_value13 = 0.00410254;
    double ST2000_p_value14 = 0.00821606;
    double ST2000_p_value15 = 0.0154771;
    double ST2000_p_value16 = 0.0274723;
    double ST2000_p_value17 = 0.0460383;
    double ST2000_p_value18 = 0.0729928;
    double ST2000_p_value19 = 0.109739;
    double ST2000_p_value20 = 0.156826;
    double ST2000_p_value21 = 0.213577;
    double ST2000_p_value22 = 0.277873;
    double ST2000_p_value23 = 0.346135;
    double ST2000_p_value24 = 0.413172;
    double ST2000_p_value25 = 0.470837;
    double ST2000_p_value26 = 0.5;
    double ST2000_p_value27 = 0.475137;
    double ST2000_p_value28 = 0.418986;
    double ST2000_p_value29 = 0.352391;
    double ST2000_p_value30 = 0.284001;
    double ST2000_p_value31 = 0.219171;
    double ST2000_p_value32 = 0.161622;
    double ST2000_p_value33 = 0.113605;
    double ST2000_p_value34 = 0.0759222;
    double ST2000_p_value35 = 0.0481236;
    double ST2000_p_value36 = 0.0288647;
    double ST2000_p_value37 = 0.0163483;
    double ST2000_p_value38 = 0.00872642;
    double ST2000_p_value39 = 0.00438208;
    double ST2000_p_value40 = 0.00206681;
    double ST2000_p_value41 = 0.00091426;
    double ST2000_p_value42 = 0.000378782;
    double ST2000_p_value43 = 0.000146803;
    double ST2000_p_value44 = 5.31692e-05;
    double ST2000_p_value45 = 1.79766e-05;
    double ST2000_p_value46 = 5.66915e-06;
    double ST2000_p_value47 = 1.66617e-06;
    double ST2000_p_value48 = 4.56063e-07;
    double ST2000_p_value49 = 1.1618e-07;
    double ST2000_p_value50 = 2.75288e-08;
    double ST2000_p_value51 = 6.06462e-09;
    
    double p_value_ST2000[51] = {ST2000_p_value1, ST2000_p_value2, ST2000_p_value3, ST2000_p_value4, ST2000_p_value5, ST2000_p_value6, ST2000_p_value7, ST2000_p_value8, ST2000_p_value9, ST2000_p_value10, ST2000_p_value11, ST2000_p_value12, ST2000_p_value13, ST2000_p_value14, ST2000_p_value15, ST2000_p_value16, ST2000_p_value17, ST2000_p_value18, ST2000_p_value19, ST2000_p_value20, ST2000_p_value21, ST2000_p_value22, ST2000_p_value23, ST2000_p_value24, ST2000_p_value25, ST2000_p_value26, ST2000_p_value27, ST2000_p_value28, ST2000_p_value29, ST2000_p_value30, ST2000_p_value31, ST2000_p_value32, ST2000_p_value33, ST2000_p_value34, ST2000_p_value35, ST2000_p_value36, ST2000_p_value37, ST2000_p_value38, ST2000_p_value39, ST2000_p_value40, ST2000_p_value41, ST2000_p_value42, ST2000_p_value43, ST2000_p_value44, ST2000_p_value45, ST2000_p_value46, ST2000_p_value47, ST2000_p_value48, ST2000_p_value49, ST2000_p_value50, ST2000_p_value51};
    
    double ST2500_p_value1 = 1.80647e-05;
    double ST2500_p_value2 = 4.02304e-05;
    double ST2500_p_value3 = 8.63723e-05;
    double ST2500_p_value4 = 0.000178825;
    double ST2500_p_value5 = 0.000357163;
    double ST2500_p_value6 = 0.000688458;
    double ST2500_p_value7 = 0.00128132;
    double ST2500_p_value8 = 0.0023036;
    double ST2500_p_value9 = 0.00400267;
    double ST2500_p_value10 = 0.00672553;
    double ST2500_p_value11 = 0.0109341;
    double ST2500_p_value12 = 0.017211;
    double ST2500_p_value13 = 0.0262469;
    double ST2500_p_value14 = 0.038808;
    double ST2500_p_value15 = 0.0556754;
    double ST2500_p_value16 = 0.0775656;
    double ST2500_p_value17 = 0.105031;
    double ST2500_p_value18 = 0.138359;
    double ST2500_p_value19 = 0.177481;
    double ST2500_p_value20 = 0.221917;
    double ST2500_p_value21 = 0.270725;
    double ST2500_p_value22 = 0.322519;
    double ST2500_p_value23 = 0.375439;
    double ST2500_p_value24 = 0.426985;
    double ST2500_p_value25 = 0.473052;
    double ST2500_p_value26 = 0.499999;
    double ST2500_p_value27 = 0.477005;
    double ST2500_p_value28 = 0.431871;
    double ST2500_p_value29 = 0.380639;
    double ST2500_p_value30 = 0.327734;
    double ST2500_p_value31 = 0.275736;
    double ST2500_p_value32 = 0.226564;
    double ST2500_p_value33 = 0.181648;
    double ST2500_p_value34 = 0.141971;
    double ST2500_p_value35 = 0.108062;
    double ST2500_p_value36 = 0.0800244;
    double ST2500_p_value37 = 0.0576044;
    double ST2500_p_value38 = 0.0402707;
    double ST2500_p_value39 = 0.0273185;
    double ST2500_p_value40 = 0.0179691;
    double ST2500_p_value41 = 0.0114519;
    double ST2500_p_value42 = 0.00706667;
    double ST2500_p_value43 = 0.00421957;
    double ST2500_p_value44 = 0.00243654;
    double ST2500_p_value45 = 0.00135987;
    double ST2500_p_value46 = 0.000733188;
    double ST2500_p_value47 = 0.000381694;
    double ST2500_p_value48 = 0.000191783;
    double ST2500_p_value49 = 9.29637e-05;
    double ST2500_p_value50 = 4.34569e-05;
    double ST2500_p_value51 = 1.95846e-05;
    
    double p_value_ST2500[51] = {ST2500_p_value1, ST2500_p_value2, ST2500_p_value3, ST2500_p_value4, ST2500_p_value5, ST2500_p_value6, ST2500_p_value7, ST2500_p_value8, ST2500_p_value9, ST2500_p_value10, ST2500_p_value11, ST2500_p_value12, ST2500_p_value13, ST2500_p_value14, ST2500_p_value15, ST2500_p_value16, ST2500_p_value17, ST2500_p_value18, ST2500_p_value19, ST2500_p_value20, ST2500_p_value21, ST2500_p_value22, ST2500_p_value23, ST2500_p_value24, ST2500_p_value25, ST2500_p_value26, ST2500_p_value27, ST2500_p_value28, ST2500_p_value29, ST2500_p_value30, ST2500_p_value31, ST2500_p_value32, ST2500_p_value33, ST2500_p_value34, ST2500_p_value35, ST2500_p_value36, ST2500_p_value37, ST2500_p_value38, ST2500_p_value39, ST2500_p_value40, ST2500_p_value41, ST2500_p_value42, ST2500_p_value43, ST2500_p_value44, ST2500_p_value45, ST2500_p_value46, ST2500_p_value47, ST2500_p_value48, ST2500_p_value49, ST2500_p_value50, ST2500_p_value51};
    
    double ST1500_p_value_Signal10Unc1 = 4.00433e-11;
    double ST1500_p_value_Signal10Unc2 = 2.89704e-10;
    double ST1500_p_value_Signal10Unc3 = 1.90293e-09;
    double ST1500_p_value_Signal10Unc4 = 1.13526e-08;
    double ST1500_p_value_Signal10Unc5 = 6.15367e-08;
    double ST1500_p_value_Signal10Unc6 = 3.03298e-07;
    double ST1500_p_value_Signal10Unc7 = 1.3601e-06;
    double ST1500_p_value_Signal10Unc8 = 5.55436e-06;
    double ST1500_p_value_Signal10Unc9 = 2.06778e-05;
    double ST1500_p_value_Signal10Unc10 = 7.02528e-05;
    double ST1500_p_value_Signal10Unc11 = 0.000218134;
    double ST1500_p_value_Signal10Unc12 = 0.000619907;
    double ST1500_p_value_Signal10Unc13 = 0.00161531;
    double ST1500_p_value_Signal10Unc14 = 0.0038669;
    double ST1500_p_value_Signal10Unc15 = 0.00852297;
    double ST1500_p_value_Signal10Unc16 = 0.0173394;
    double ST1500_p_value_Signal10Unc17 = 0.0326507;
    double ST1500_p_value_Signal10Unc18 = 0.0570807;
    double ST1500_p_value_Signal10Unc19 = 0.0929579;
    double ST1500_p_value_Signal10Unc20 = 0.141534;
    double ST1500_p_value_Signal10Unc21 = 0.202242;
    double ST1500_p_value_Signal10Unc22 = 0.272259;
    double ST1500_p_value_Signal10Unc23 = 0.346449;
    double ST1500_p_value_Signal10Unc24 = 0.41741;
    double ST1500_p_value_Signal10Unc25 = 0.474423;
    double ST1500_p_value_Signal10Unc26 = 0.5;
    double ST1500_p_value_Signal10Unc27 = 0.479342;
    double ST1500_p_value_Signal10Unc28 = 0.42488;
    double ST1500_p_value_Signal10Unc29 = 0.354883;
    double ST1500_p_value_Signal10Unc30 = 0.280641;
    double ST1500_p_value_Signal10Unc31 = 0.209845;
    double ST1500_p_value_Signal10Unc32 = 0.147887;
    double ST1500_p_value_Signal10Unc33 = 0.0978567;
    double ST1500_p_value_Signal10Unc34 = 0.0605638;
    double ST1500_p_value_Signal10Unc35 = 0.034931;
    double ST1500_p_value_Signal10Unc36 = 0.0187115;
    double ST1500_p_value_Signal10Unc37 = 0.00928035;
    double ST1500_p_value_Signal10Unc38 = 0.00424965;
    double ST1500_p_value_Signal10Unc39 = 0.00179231;
    double ST1500_p_value_Signal10Unc40 = 0.000694602;
    double ST1500_p_value_Signal10Unc41 = 0.000246862;
    double ST1500_p_value_Signal10Unc42 = 8.03235e-05;
    double ST1500_p_value_Signal10Unc43 = 2.3888e-05;
    double ST1500_p_value_Signal10Unc44 = 6.48441e-06;
    double ST1500_p_value_Signal10Unc45 = 1.6048e-06;
    double ST1500_p_value_Signal10Unc46 = 3.61729e-07;
    double ST1500_p_value_Signal10Unc47 = 7.41948e-08;
    double ST1500_p_value_Signal10Unc48 = 1.38371e-08;
    double ST1500_p_value_Signal10Unc49 = 2.34496e-09;
    double ST1500_p_value_Signal10Unc50 = 3.60972e-10;
    double ST1500_p_value_Signal10Unc51 = 5.04486e-11;
    
    double ST1500_p_value_Signal10Unc[51] = {ST1500_p_value_Signal10Unc1, ST1500_p_value_Signal10Unc2, ST1500_p_value_Signal10Unc3, ST1500_p_value_Signal10Unc4, ST1500_p_value_Signal10Unc5, ST1500_p_value_Signal10Unc6, ST1500_p_value_Signal10Unc7, ST1500_p_value_Signal10Unc8, ST1500_p_value_Signal10Unc9, ST1500_p_value_Signal10Unc10, ST1500_p_value_Signal10Unc11, ST1500_p_value_Signal10Unc12, ST1500_p_value_Signal10Unc13, ST1500_p_value_Signal10Unc14, ST1500_p_value_Signal10Unc15, ST1500_p_value_Signal10Unc16, ST1500_p_value_Signal10Unc17, ST1500_p_value_Signal10Unc18, ST1500_p_value_Signal10Unc19, ST1500_p_value_Signal10Unc20, ST1500_p_value_Signal10Unc21, ST1500_p_value_Signal10Unc22, ST1500_p_value_Signal10Unc23, ST1500_p_value_Signal10Unc24, ST1500_p_value_Signal10Unc25, ST1500_p_value_Signal10Unc26, ST1500_p_value_Signal10Unc27, ST1500_p_value_Signal10Unc28, ST1500_p_value_Signal10Unc29, ST1500_p_value_Signal10Unc30, ST1500_p_value_Signal10Unc31, ST1500_p_value_Signal10Unc32, ST1500_p_value_Signal10Unc33, ST1500_p_value_Signal10Unc34, ST1500_p_value_Signal10Unc35, ST1500_p_value_Signal10Unc36, ST1500_p_value_Signal10Unc37, ST1500_p_value_Signal10Unc38, ST1500_p_value_Signal10Unc39, ST1500_p_value_Signal10Unc40, ST1500_p_value_Signal10Unc41, ST1500_p_value_Signal10Unc42, ST1500_p_value_Signal10Unc43, ST1500_p_value_Signal10Unc44, ST1500_p_value_Signal10Unc45, ST1500_p_value_Signal10Unc46, ST1500_p_value_Signal10Unc47, ST1500_p_value_Signal10Unc48, ST1500_p_value_Signal10Unc49, ST1500_p_value_Signal10Unc50, ST1500_p_value_Signal10Unc51};
    
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(51, param, p_value_noST);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p_value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(51, param, p_value_ST1000);
    gr_1->SetMarkerStyle(20);
    gr_1->SetMarkerColor(kCyan);
    gr_1->SetLineColor(kCyan);
    gr_1->SetMarkerSize(0.9);
    gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(51, param, p_value_ST1500);
    gr_2->SetMarkerStyle(20);
    gr_2->SetMarkerColor(kBlack);
    gr_2->SetLineColor(kBlack);
    gr_2->SetMarkerSize(0.9);
    gr_2->Draw("LP*");
    
    TGraph *gr_3 = new TGraph(51, param, p_value_ST2000);
    gr_3->SetMarkerStyle(20);
    gr_3->SetMarkerColor(kGreen+3);
    gr_3->SetLineColor(kGreen+3);
    gr_3->SetMarkerSize(0.9);
    gr_3->Draw("LP*");
    
    TGraph *gr_4 = new TGraph(51, param, p_value_ST2500);
    gr_4->SetMarkerStyle(20);
    gr_4->SetMarkerColor(kMagenta);
    gr_4->SetLineColor(kMagenta);
    gr_4->SetMarkerSize(0.9);
    gr_4->Draw("LP*");
    
    TGraph *gr_6 = new TGraph(51, param, ST1500_p_value_Signal10Unc);
    gr_6->SetMarkerStyle(20);
    gr_6->SetMarkerColor(kAzure+4);
    gr_6->SetLineColor(kAzure+4);
    gr_6->SetMarkerSize(0.9);
    gr_6->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.65,0.70,0.89,0.89,NULL,"brNDC");
    leg->AddEntry(gr, "no ST cut", "pl");
    leg->AddEntry(gr_1, "ST 1000", "pl");
    leg->AddEntry(gr_2, "ST 1500", "pl");
    leg->AddEntry(gr_6, "ST 1500: Signal Syst 10%", "pl");
    leg->AddEntry(gr_3, "ST 2000", "pl");
    leg->AddEntry(gr_4, "ST 2500", "pl");
    leg->SetBorderSize(0);
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr_2->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit ST 1500 = " << m << endl;
    
    while (m<=2.0 and gr_2->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit ST 1500 = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_3->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit ST 2000 = " << m1 << endl;
    
    while (m1<=2.0 and gr_3->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit ST 2000 = " << m1 << endl;
    
    c1.SaveAs("Parabola_ft0_p_value_Finegrid_All.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Finegrid_All.png");
    
}

void makeGraphs_pValue_combined_withsignalUnc()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double ST_p_value1 = 2.00092e-14;
    double ST_p_value2 = 2.96911e-13;
    double ST_p_value3 = 3.87444e-12;
    double ST_p_value4 = 4.43719e-11;
    double ST_p_value5 = 4.45959e-10;
    double ST_p_value6 = 3.93108e-09;
    double ST_p_value7 = 3.03592e-08;
    double ST_p_value8 = 2.05597e-07;
    double ST_p_value9 = 1.2211e-06;
    double ST_p_value10 = 6.36445e-06;
    double ST_p_value11 = 2.91611e-05;
    double ST_p_value12 = 0.00011762;
    double ST_p_value13 = 0.00041846;
    double ST_p_value14 = 0.00131697;
    double ST_p_value15 = 0.00367947;
    double ST_p_value16 = 0.00915728;
    double ST_p_value17 = 0.0203963;
    double ST_p_value18 = 0.0408639;
    double ST_p_value19 = 0.0741211;
    double ST_p_value20 = 0.122441;
    double ST_p_value21 = 0.185583;
    double ST_p_value22 = 0.260065;
    double ST_p_value23 = 0.339134;
    double ST_p_value24 = 0.41418;
    double ST_p_value25 = 0.473493;
    double ST_p_value26 = 0.5;
    double ST_p_value27 = 0.478901;
    double ST_p_value28 = 0.422322;
    double ST_p_value29 = 0.348458;
    double ST_p_value30 = 0.269423;
    double ST_p_value31 = 0.194096;
    double ST_p_value32 = 0.129374;
    double ST_p_value33 = 0.0791853;
    double ST_p_value34 = 0.0442033;
    double ST_p_value35 = 0.0223485;
    double ST_p_value36 = 0.0101711;
    double ST_p_value37 = 0.00414478;
    double ST_p_value38 = 0.0015063;
    double ST_p_value39 = 0.000485981;
    double ST_p_value40 = 0.000138753;
    double ST_p_value41 = 3.49547e-05;
    double ST_p_value42 = 7.75757e-06;
    double ST_p_value43 = 1.51304e-06;
    double ST_p_value44 = 2.59008e-07;
    double ST_p_value45 = 3.88997e-08;
    double ST_p_value46 = 5.12168e-09;
    double ST_p_value47 = 5.911e-10;
    double ST_p_value48 = 5.97987e-11;
    double ST_p_value49 = 5.30843e-12;
    double ST_p_value50 = 4.13423e-13;
    double ST_p_value51 = 2.83312e-14;
    
    double ST_p_value[51] = {ST_p_value1, ST_p_value2, ST_p_value3, ST_p_value4, ST_p_value5, ST_p_value6, ST_p_value7, ST_p_value8, ST_p_value9, ST_p_value10, ST_p_value11, ST_p_value12, ST_p_value13, ST_p_value14, ST_p_value15, ST_p_value16, ST_p_value17, ST_p_value18, ST_p_value19, ST_p_value20, ST_p_value21, ST_p_value22, ST_p_value23, ST_p_value24, ST_p_value25, ST_p_value26, ST_p_value27, ST_p_value28, ST_p_value29, ST_p_value30, ST_p_value31, ST_p_value32, ST_p_value33, ST_p_value34, ST_p_value35, ST_p_value36, ST_p_value37, ST_p_value38, ST_p_value39, ST_p_value40, ST_p_value41, ST_p_value42, ST_p_value43, ST_p_value44, ST_p_value45, ST_p_value46, ST_p_value47, ST_p_value48, ST_p_value49, ST_p_value50, ST_p_value51};

    double ST_p_value_Signal10Unc1 = 2.00091e-14;
    double ST_p_value_Signal10Unc2 = 2.9691e-13;
    double ST_p_value_Signal10Unc3 = 3.87444e-12;
    double ST_p_value_Signal10Unc4 = 4.43719e-11;
    double ST_p_value_Signal10Unc5 = 4.45958e-10;
    double ST_p_value_Signal10Unc6 = 3.93107e-09;
    double ST_p_value_Signal10Unc7 = 3.03592e-08;
    double ST_p_value_Signal10Unc8 = 2.05597e-07;
    double ST_p_value_Signal10Unc9 = 1.2211e-06;
    double ST_p_value_Signal10Unc10 = 6.36446e-06;
    double ST_p_value_Signal10Unc11 = 2.91611e-05;
    double ST_p_value_Signal10Unc12 = 0.00011762;
    double ST_p_value_Signal10Unc13 = 0.00041846;
    double ST_p_value_Signal10Unc14 = 0.00131697;
    double ST_p_value_Signal10Unc15 = 0.00367947;
    double ST_p_value_Signal10Unc16 = 0.00915728;
    double ST_p_value_Signal10Unc17 = 0.0203963;
    double ST_p_value_Signal10Unc18 = 0.0408639;
    double ST_p_value_Signal10Unc19 = 0.0741211;
    double ST_p_value_Signal10Unc20 = 0.122441;
    double ST_p_value_Signal10Unc21 = 0.185583;
    double ST_p_value_Signal10Unc22 = 0.260065;
    double ST_p_value_Signal10Unc23 = 0.339136;
    double ST_p_value_Signal10Unc24 = 0.41418;
    double ST_p_value_Signal10Unc25 = 0.473493;
    double ST_p_value_Signal10Unc26 = 0.5;
    double ST_p_value_Signal10Unc27 = 0.478879;
    double ST_p_value_Signal10Unc28 = 0.422322;
    double ST_p_value_Signal10Unc29 = 0.348458;
    double ST_p_value_Signal10Unc30 = 0.269423;
    double ST_p_value_Signal10Unc31 = 0.194096;
    double ST_p_value_Signal10Unc32 = 0.129374;
    double ST_p_value_Signal10Unc33 = 0.0791853;
    double ST_p_value_Signal10Unc34 = 0.0442033;
    double ST_p_value_Signal10Unc35 = 0.0223485;
    double ST_p_value_Signal10Unc36 = 0.0101711;
    double ST_p_value_Signal10Unc37 = 0.00414478;
    double ST_p_value_Signal10Unc38 = 0.0015063;
    double ST_p_value_Signal10Unc39 = 0.000485981;
    double ST_p_value_Signal10Unc40 = 0.000138753;
    double ST_p_value_Signal10Unc41 = 3.49547e-05;
    double ST_p_value_Signal10Unc42 = 7.75758e-06;
    double ST_p_value_Signal10Unc43 = 1.51305e-06;
    double ST_p_value_Signal10Unc44 = 2.59008e-07;
    double ST_p_value_Signal10Unc45 = 3.88997e-08;
    double ST_p_value_Signal10Unc46 = 5.12168e-09;
    double ST_p_value_Signal10Unc47 = 5.911e-10;
    double ST_p_value_Signal10Unc48 = 5.97986e-11;
    double ST_p_value_Signal10Unc49 = 5.30842e-12;
    double ST_p_value_Signal10Unc50 = 4.13422e-13;
    double ST_p_value_Signal10Unc51 = 2.83312e-14;
   
    double ST_p_value_Signal10Unc[51] = {ST_p_value_Signal10Unc1, ST_p_value_Signal10Unc2, ST_p_value_Signal10Unc3, ST_p_value_Signal10Unc4, ST_p_value_Signal10Unc5, ST_p_value_Signal10Unc6, ST_p_value_Signal10Unc7, ST_p_value_Signal10Unc8, ST_p_value_Signal10Unc9, ST_p_value_Signal10Unc10, ST_p_value_Signal10Unc11, ST_p_value_Signal10Unc12, ST_p_value_Signal10Unc13, ST_p_value_Signal10Unc14, ST_p_value_Signal10Unc15, ST_p_value_Signal10Unc16, ST_p_value_Signal10Unc17, ST_p_value_Signal10Unc18, ST_p_value_Signal10Unc19, ST_p_value_Signal10Unc20, ST_p_value_Signal10Unc21, ST_p_value_Signal10Unc22, ST_p_value_Signal10Unc23, ST_p_value_Signal10Unc24, ST_p_value_Signal10Unc25, ST_p_value_Signal10Unc26, ST_p_value_Signal10Unc27, ST_p_value_Signal10Unc28, ST_p_value_Signal10Unc29, ST_p_value_Signal10Unc30, ST_p_value_Signal10Unc31, ST_p_value_Signal10Unc32, ST_p_value_Signal10Unc33, ST_p_value_Signal10Unc34, ST_p_value_Signal10Unc35, ST_p_value_Signal10Unc36, ST_p_value_Signal10Unc37, ST_p_value_Signal10Unc38, ST_p_value_Signal10Unc39, ST_p_value_Signal10Unc40, ST_p_value_Signal10Unc41, ST_p_value_Signal10Unc42, ST_p_value_Signal10Unc43, ST_p_value_Signal10Unc44, ST_p_value_Signal10Unc45, ST_p_value_Signal10Unc46, ST_p_value_Signal10Unc47, ST_p_value_Signal10Unc48, ST_p_value_Signal10Unc49, ST_p_value_Signal10Unc50, ST_p_value_Signal10Unc51};

    double ST_p_value_Signal25Unc1 = 5.0779e-11;
    double ST_p_value_Signal25Unc2 = 3.82803e-10;
    double ST_p_value_Signal25Unc3 = 2.60936e-09;
    double ST_p_value_Signal25Unc4 = 1.60768e-08;
    double ST_p_value_Signal25Unc5 = 8.9488e-08;
    double ST_p_value_Signal25Unc6 = 4.50608e-07;
    double ST_p_value_Signal25Unc7 = 2.05207e-06;
    double ST_p_value_Signal25Unc8 = 8.45663e-06;
    double ST_p_value_Signal25Unc9 = 3.15926e-05;
    double ST_p_value_Signal25Unc10 = 0.000107013;
    double ST_p_value_Signal25Unc11 = 0.000329479;
    double ST_p_value_Signal25Unc12 = 0.000922798;
    double ST_p_value_Signal25Unc13 = 0.00235874;
    double ST_p_value_Signal25Unc14 = 0.00550869;
    double ST_p_value_Signal25Unc15 = 0.0117956;
    double ST_p_value_Signal25Unc16 = 0.0232318;
    double ST_p_value_Signal25Unc17 = 0.042238;
    double ST_p_value_Signal25Unc18 = 0.0711818;
    double ST_p_value_Signal25Unc19 = 0.111659;
    double ST_p_value_Signal25Unc20 = 0.16384;
    double ST_p_value_Signal25Unc21 = 0.225961;
    double ST_p_value_Signal25Unc22 = 0.294416;
    double ST_p_value_Signal25Unc23 = 0.363978;
    double ST_p_value_Signal25Unc24 = 0.428187;
    double ST_p_value_Signal25Unc25 = 0.47794;
    double ST_p_value_Signal25Unc26 = 0.5;
    double ST_p_value_Signal25Unc27 = 0.482533;
    double ST_p_value_Signal25Unc28 = 0.434834;
    double ST_p_value_Signal25Unc29 = 0.372009;
    double ST_p_value_Signal25Unc30 = 0.302727;
    double ST_p_value_Signal25Unc31 = 0.233885;
    double ST_p_value_Signal25Unc32 = 0.17089;
    double ST_p_value_Signal25Unc33 = 0.117393;
    double ST_p_value_Signal25Unc34 = 0.0754773;
    double ST_p_value_Signal25Unc35 = 0.0452019;
    double ST_p_value_Signal25Unc36 = 0.025115;
    double ST_p_value_Signal25Unc37 = 0.012882;
    double ST_p_value_Signal25Unc38 = 0.00608018;
    double ST_p_value_Signal25Unc39 = 0.00263262;
    double ST_p_value_Signal25Unc40 = 0.00104187;
    double ST_p_value_Signal25Unc41 = 0.000376346;
    double ST_p_value_Signal25Unc42 = 0.000123758;
    double ST_p_value_Signal25Unc43 = 3.69881e-05;
    double ST_p_value_Signal25Unc44 = 1.00232e-05;
    double ST_p_value_Signal25Unc45 = 2.46309e-06;
    double ST_p_value_Signal25Unc46 = 5.47852e-07;
    double ST_p_value_Signal25Unc47 = 1.10202e-07;
    double ST_p_value_Signal25Unc48 = 2.00501e-08;
    double ST_p_value_Signal25Unc49 = 3.29643e-09;
    double ST_p_value_Signal25Unc50 = 4.9013e-10;
    double ST_p_value_Signal25Unc51 = 6.58268e-11;

    double ST_p_value_Signal25Unc[51] = {ST_p_value_Signal25Unc1, ST_p_value_Signal25Unc2, ST_p_value_Signal25Unc3, ST_p_value_Signal25Unc4, ST_p_value_Signal25Unc5, ST_p_value_Signal25Unc6, ST_p_value_Signal25Unc7, ST_p_value_Signal25Unc8, ST_p_value_Signal25Unc9, ST_p_value_Signal25Unc10, ST_p_value_Signal25Unc11, ST_p_value_Signal25Unc12, ST_p_value_Signal25Unc13, ST_p_value_Signal25Unc14, ST_p_value_Signal25Unc15, ST_p_value_Signal25Unc16, ST_p_value_Signal25Unc17, ST_p_value_Signal25Unc18, ST_p_value_Signal25Unc19, ST_p_value_Signal25Unc20, ST_p_value_Signal25Unc21, ST_p_value_Signal25Unc22, ST_p_value_Signal25Unc23, ST_p_value_Signal25Unc24, ST_p_value_Signal25Unc25, ST_p_value_Signal25Unc26, ST_p_value_Signal25Unc27, ST_p_value_Signal25Unc28, ST_p_value_Signal25Unc29, ST_p_value_Signal25Unc30, ST_p_value_Signal25Unc31, ST_p_value_Signal25Unc32, ST_p_value_Signal25Unc33, ST_p_value_Signal25Unc34, ST_p_value_Signal25Unc35, ST_p_value_Signal25Unc36, ST_p_value_Signal25Unc37, ST_p_value_Signal25Unc38, ST_p_value_Signal25Unc39, ST_p_value_Signal25Unc40, ST_p_value_Signal25Unc41, ST_p_value_Signal25Unc42, ST_p_value_Signal25Unc43, ST_p_value_Signal25Unc44, ST_p_value_Signal25Unc45, ST_p_value_Signal25Unc46, ST_p_value_Signal25Unc47, ST_p_value_Signal25Unc48, ST_p_value_Signal25Unc49, ST_p_value_Signal25Unc50, ST_p_value_Signal25Unc51};
 
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    TGraph *gr = new TGraph(51, param, ST_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    //gr->SetMinimum(0.02);
    gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(51, param, ST_p_value_Signal10Unc);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    //gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(51, param, ST_p_value_Signal25Unc);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.38,0.30,0.62,0.50,NULL,"brNDC");
    leg->SetHeader("Combined Limit");
    leg->AddEntry(gr,   "Signal: no Syst, Bkg: 30%", "pl");
    //leg->AddEntry(gr_1, "Signal: 10% Syst, Bkg: 30%", "pl");
    leg->AddEntry(gr_2, "Signal: 25% reduced", "pl");
    leg->SetBorderSize(0);
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 30% syst = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 30% syst = " << m << endl;
    /*
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 30% syst = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 30% syst = " << m1 << endl;
    */
    double m2=-2.0;
    while (m2<=2.0 and gr_2->Eval(m2)<gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected lower limit combined 25% signal reduction = " << m2 << endl;
    
    while (m2<=2.0 and gr_2->Eval(m2)>gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected upper limit combined 25% signal reduction = " << m2 << endl;
    
    
    c1.SaveAs("Parabola_ft0_p_value_All_SignalSyst_10_Bkg_30.pdf");
    c1.SaveAs("Parabola_ft0_p_value_All_SignalSyst_10_Bkg_30.png");
}

void makeGraphs_pValue_STshape_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double ZeroOSSF_p_value1 = 0.000386417;
    double ZeroOSSF_p_value2 = 0.00068111;
    double ZeroOSSF_p_value3 = 0.00116893;
    double ZeroOSSF_p_value4 = 0.00195411;
    double ZeroOSSF_p_value5 = 0.00318303;
    double ZeroOSSF_p_value6 = 0.00505382;
    double ZeroOSSF_p_value7 = 0.00782458;
    double ZeroOSSF_p_value8 = 0.011818;
    double ZeroOSSF_p_value9 = 0.0174202;
    double ZeroOSSF_p_value10 = 0.0250727;
    double ZeroOSSF_p_value11 = 0.0352532;
    double ZeroOSSF_p_value12 = 0.0484476;
    double ZeroOSSF_p_value13 = 0.0651127;
    double ZeroOSSF_p_value14 = 0.0856308;
    double ZeroOSSF_p_value15 = 0.110262;
    double ZeroOSSF_p_value16 = 0.139102;
    double ZeroOSSF_p_value17 = 0.172039;
    double ZeroOSSF_p_value18 = 0.208735;
    double ZeroOSSF_p_value19 = 0.248611;
    double ZeroOSSF_p_value20 = 0.29084;
    double ZeroOSSF_p_value21 = 0.334366;
    double ZeroOSSF_p_value22 = 0.377868;
    double ZeroOSSF_p_value23 = 0.419673;
    double ZeroOSSF_p_value24 = 0.457415;
    double ZeroOSSF_p_value25 = 0.486973;
    double ZeroOSSF_p_value26 = 0.5;
    double ZeroOSSF_p_value27 = 0.488883;
    double ZeroOSSF_p_value28 = 0.460331;
    double ZeroOSSF_p_value29 = 0.423079;
    double ZeroOSSF_p_value30 = 0.381515;
    double ZeroOSSF_p_value31 = 0.338087;
    double ZeroOSSF_p_value32 = 0.294508;
    double ZeroOSSF_p_value33 = 0.252123;
    double ZeroOSSF_p_value34 = 0.212012;
    double ZeroOSSF_p_value35 = 0.175019;
    double ZeroOSSF_p_value36 = 0.141746;
    double ZeroOSSF_p_value37 = 0.112551;
    double ZeroOSSF_p_value38 = 0.0875626;
    double ZeroOSSF_p_value39 = 0.066703;
    double ZeroOSSF_p_value40 = 0.0497238;
    double ZeroOSSF_p_value41 = 0.0362512;
    double ZeroOSSF_p_value42 = 0.0258332;
    double ZeroOSSF_p_value43 = 0.0179847;
    double ZeroOSSF_p_value44 = 0.0122259;
    double ZeroOSSF_p_value45 = 0.00811159;
    double ZeroOSSF_p_value46 = 0.00525035;
    double ZeroOSSF_p_value47 = 0.00331394;
    double ZeroOSSF_p_value48 = 0.00203894;
    double ZeroOSSF_p_value49 = 0.00122239;
    double ZeroOSSF_p_value50 = 0.00071383;
    double ZeroOSSF_p_value51 = 0.000405919;
    
    double ZeroOSSF_p_value[51] = {ZeroOSSF_p_value1, ZeroOSSF_p_value2, ZeroOSSF_p_value3, ZeroOSSF_p_value4, ZeroOSSF_p_value5, ZeroOSSF_p_value6, ZeroOSSF_p_value7, ZeroOSSF_p_value8, ZeroOSSF_p_value9, ZeroOSSF_p_value10, ZeroOSSF_p_value11, ZeroOSSF_p_value12, ZeroOSSF_p_value13, ZeroOSSF_p_value14, ZeroOSSF_p_value15, ZeroOSSF_p_value16, ZeroOSSF_p_value17, ZeroOSSF_p_value18, ZeroOSSF_p_value19, ZeroOSSF_p_value20, ZeroOSSF_p_value21, ZeroOSSF_p_value22, ZeroOSSF_p_value23, ZeroOSSF_p_value24, ZeroOSSF_p_value25, ZeroOSSF_p_value26, ZeroOSSF_p_value27, ZeroOSSF_p_value28, ZeroOSSF_p_value29, ZeroOSSF_p_value30, ZeroOSSF_p_value31, ZeroOSSF_p_value32, ZeroOSSF_p_value33, ZeroOSSF_p_value34, ZeroOSSF_p_value35, ZeroOSSF_p_value36, ZeroOSSF_p_value37, ZeroOSSF_p_value38, ZeroOSSF_p_value39, ZeroOSSF_p_value40, ZeroOSSF_p_value41, ZeroOSSF_p_value42, ZeroOSSF_p_value43, ZeroOSSF_p_value44, ZeroOSSF_p_value45, ZeroOSSF_p_value46, ZeroOSSF_p_value47, ZeroOSSF_p_value48, ZeroOSSF_p_value49, ZeroOSSF_p_value50, ZeroOSSF_p_value51};
    /*
    double OneOSSF_p_value1 = 1.74308e-05;
    double OneOSSF_p_value2 = 4.03969e-05;
    double OneOSSF_p_value3 = 8.98668e-05;
    double OneOSSF_p_value4 = 0.000192015;
    double OneOSSF_p_value5 = 0.000394194;
    double OneOSSF_p_value6 = 0.000777979;
    double OneOSSF_p_value7 = 0.00147681;
    double OneOSSF_p_value8 = 0.00269816;
    double OneOSSF_p_value9 = 0.0047474;
    double OneOSSF_p_value10 = 0.00804998;
    double OneOSSF_p_value11 = 0.0131647;
    double OneOSSF_p_value12 = 0.0207801;
    double OneOSSF_p_value13 = 0.0316873;
    double OneOSSF_p_value14 = 0.0467217;
    double OneOSSF_p_value15 = 0.0666768;
    double OneOSSF_p_value16 = 0.0921952;
    double OneOSSF_p_value17 = 0.12365;
    double OneOSSF_p_value18 = 0.161036;
    double OneOSSF_p_value19 = 0.203895;
    double OneOSSF_p_value20 = 0.25127;
    double OneOSSF_p_value21 = 0.301707;
    double OneOSSF_p_value22 = 0.353259;
    double OneOSSF_p_value23 = 0.403433;
    double OneOSSF_p_value24 = 0.448875;
    double OneOSSF_p_value25 = 0.484225;
    double OneOSSF_p_value26 = 0.5;
    double OneOSSF_p_value27 = 0.488011;
    double OneOSSF_p_value28 = 0.454796;
    double OneOSSF_p_value29 = 0.410415;
    double OneOSSF_p_value30 = 0.36069;
    double OneOSSF_p_value31 = 0.309167;
    double OneOSSF_p_value32 = 0.258437;
    double OneOSSF_p_value33 = 0.210518;
    double OneOSSF_p_value34 = 0.166932;
    double OneOSSF_p_value35 = 0.128711;
    double OneOSSF_p_value36 = 0.0963855;
    double OneOSSF_p_value37 = 0.0700214;
    double OneOSSF_p_value38 = 0.0492941;
    double OneOSSF_p_value39 = 0.0335928;
    double OneOSSF_p_value40 = 0.0221388;
    double OneOSSF_p_value41 = 0.0140967;
    double OneOSSF_p_value42 = 0.00866471;
    double OneOSSF_p_value43 = 0.00513707;
    double OneOSSF_p_value44 = 0.00293545;
    double OneOSSF_p_value45 = 0.00161557;
    double OneOSSF_p_value46 = 0.000855816;
    double OneOSSF_p_value47 = 0.000436107;
    double OneOSSF_p_value48 = 0.00021366;
    double OneOSSF_p_value49 = 0.000100582;
    double OneOSSF_p_value50 = 4.54787e-05;
    double OneOSSF_p_value51 = 1.97421e-05;
    
    double OneOSSF_p_value[51] = {OneOSSF_p_value1, OneOSSF_p_value2, OneOSSF_p_value3, OneOSSF_p_value4, OneOSSF_p_value5, OneOSSF_p_value6, OneOSSF_p_value7, OneOSSF_p_value8, OneOSSF_p_value9, OneOSSF_p_value10, OneOSSF_p_value11, OneOSSF_p_value12, OneOSSF_p_value13, OneOSSF_p_value14, OneOSSF_p_value15, OneOSSF_p_value16, OneOSSF_p_value17, OneOSSF_p_value18, OneOSSF_p_value19, OneOSSF_p_value20, OneOSSF_p_value21, OneOSSF_p_value22, OneOSSF_p_value23, OneOSSF_p_value24, OneOSSF_p_value25, OneOSSF_p_value26, OneOSSF_p_value27, OneOSSF_p_value28, OneOSSF_p_value29, OneOSSF_p_value30, OneOSSF_p_value31, OneOSSF_p_value32, OneOSSF_p_value33, OneOSSF_p_value34, OneOSSF_p_value35, OneOSSF_p_value36, OneOSSF_p_value37, OneOSSF_p_value38, OneOSSF_p_value39, OneOSSF_p_value40, OneOSSF_p_value41, OneOSSF_p_value42, OneOSSF_p_value43, OneOSSF_p_value44, OneOSSF_p_value45, OneOSSF_p_value46, OneOSSF_p_value47, OneOSSF_p_value48, OneOSSF_p_value49, OneOSSF_p_value50, OneOSSF_p_value51};
    
    double TwoOSSF_p_value1 = 5.81205e-05;
    double TwoOSSF_p_value2 = 0.000120464;
    double TwoOSSF_p_value3 = 0.000241147;
    double TwoOSSF_p_value4 = 0.00046642;
    double TwoOSSF_p_value5 = 0.000872007;
    double TwoOSSF_p_value6 = 0.00157653;
    double TwoOSSF_p_value7 = 0.00275766;
    double TwoOSSF_p_value8 = 0.00466934;
    double TwoOSSF_p_value9 = 0.00765739;
    double TwoOSSF_p_value10 = 0.0121697;
    double TwoOSSF_p_value11 = 0.018755;
    double TwoOSSF_p_value12 = 0.0280475;
    double TwoOSSF_p_value13 = 0.0407302;
    double TwoOSSF_p_value14 = 0.0574802;
    double TwoOSSF_p_value15 = 0.0788956;
    double TwoOSSF_p_value16 = 0.105412;
    double TwoOSSF_p_value17 = 0.13722;
    double TwoOSSF_p_value18 = 0.174195;
    double TwoOSSF_p_value19 = 0.215848;
    double TwoOSSF_p_value20 = 0.261306;
    double TwoOSSF_p_value21 = 0.309317;
    double TwoOSSF_p_value22 = 0.358235;
    double TwoOSSF_p_value23 = 0.405971;
    double TwoOSSF_p_value24 = 0.449634;
    double TwoOSSF_p_value25 = 0.484342;
    double TwoOSSF_p_value26 = 0.5;
    double TwoOSSF_p_value27 = 0.48686;
    double TwoOSSF_p_value28 = 0.453416;
    double TwoOSSF_p_value29 = 0.41034;
    double TwoOSSF_p_value30 = 0.362851;
    double TwoOSSF_p_value31 = 0.313944;
    double TwoOSSF_p_value32 = 0.265773;
    double TwoOSSF_p_value33 = 0.220013;
    double TwoOSSF_p_value34 = 0.177954;
    double TwoOSSF_p_value35 = 0.140509;
    double TwoOSSF_p_value36 = 0.108201;
    double TwoOSSF_p_value37 = 0.0811863;
    double TwoOSSF_p_value38 = 0.0593025;
    double TwoOSSF_p_value39 = 0.0421335;
    double TwoOSSF_p_value40 = 0.0290936;
    double TwoOSSF_p_value41 = 0.0195095;
    double TwoOSSF_p_value42 = 0.0126957;
    double TwoOSSF_p_value43 = 0.00801199;
    double TwoOSSF_p_value44 = 0.00490024;
    double TwoOSSF_p_value45 = 0.00290288;
    double TwoOSSF_p_value46 = 0.00166472;
    double TwoOSSF_p_value47 = 0.000923654;
    double TwoOSSF_p_value48 = 0.000495629;
    double TwoOSSF_p_value49 = 0.000257081;
    double TwoOSSF_p_value50 = 0.000128845;
    double TwoOSSF_p_value51 = 6.23709e-05;
    
    double TwoOSSF_p_value[51] = {TwoOSSF_p_value1, TwoOSSF_p_value2, TwoOSSF_p_value3, TwoOSSF_p_value4, TwoOSSF_p_value5, TwoOSSF_p_value6, TwoOSSF_p_value7, TwoOSSF_p_value8, TwoOSSF_p_value9, TwoOSSF_p_value10, TwoOSSF_p_value11, TwoOSSF_p_value12, TwoOSSF_p_value13, TwoOSSF_p_value14, TwoOSSF_p_value15, TwoOSSF_p_value16, TwoOSSF_p_value17, TwoOSSF_p_value18, TwoOSSF_p_value19, TwoOSSF_p_value20, TwoOSSF_p_value21, TwoOSSF_p_value22, TwoOSSF_p_value23, TwoOSSF_p_value24, TwoOSSF_p_value25, TwoOSSF_p_value26, TwoOSSF_p_value27, TwoOSSF_p_value28, TwoOSSF_p_value29, TwoOSSF_p_value30, TwoOSSF_p_value31, TwoOSSF_p_value32, TwoOSSF_p_value33, TwoOSSF_p_value34, TwoOSSF_p_value35, TwoOSSF_p_value36, TwoOSSF_p_value37, TwoOSSF_p_value38, TwoOSSF_p_value39, TwoOSSF_p_value40, TwoOSSF_p_value41, TwoOSSF_p_value42, TwoOSSF_p_value43, TwoOSSF_p_value44, TwoOSSF_p_value45, TwoOSSF_p_value46, TwoOSSF_p_value47, TwoOSSF_p_value48, TwoOSSF_p_value49, TwoOSSF_p_value50, TwoOSSF_p_value51};
    
    double ST1500_p_value1 = 4.00434e-11;
    double ST1500_p_value2 = 2.89704e-10;
    double ST1500_p_value3 = 1.90293e-09;
    double ST1500_p_value4 = 1.13526e-08;
    double ST1500_p_value5 = 6.15367e-08;
    double ST1500_p_value6 = 3.03298e-07;
    double ST1500_p_value7 = 1.3601e-06;
    double ST1500_p_value8 = 5.55436e-06;
    double ST1500_p_value9 = 2.06777e-05;
    double ST1500_p_value10 = 7.02526e-05;
    double ST1500_p_value11 = 0.000218133;
    double ST1500_p_value12 = 0.000619906;
    double ST1500_p_value13 = 0.0016153;
    double ST1500_p_value14 = 0.0038669;
    double ST1500_p_value15 = 0.00852296;
    double ST1500_p_value16 = 0.0173394;
    double ST1500_p_value17 = 0.0326507;
    double ST1500_p_value18 = 0.0570807;
    double ST1500_p_value19 = 0.0929579;
    double ST1500_p_value20 = 0.141534;
    double ST1500_p_value21 = 0.202242;
    double ST1500_p_value22 = 0.272258;
    double ST1500_p_value23 = 0.346448;
    double ST1500_p_value24 = 0.41741;
    double ST1500_p_value25 = 0.474421;
    double ST1500_p_value26 = 0.5;
    double ST1500_p_value27 = 0.47934;
    double ST1500_p_value28 = 0.42488;
    double ST1500_p_value29 = 0.354882;
    double ST1500_p_value30 = 0.28064;
    double ST1500_p_value31 = 0.209845;
    double ST1500_p_value32 = 0.147887;
    double ST1500_p_value33 = 0.0978567;
    double ST1500_p_value34 = 0.0605638;
    double ST1500_p_value35 = 0.034931;
    double ST1500_p_value36 = 0.0187115;
    double ST1500_p_value37 = 0.00928034;
    double ST1500_p_value38 = 0.00424965;
    double ST1500_p_value39 = 0.00179231;
    double ST1500_p_value40 = 0.000694601;
    double ST1500_p_value41 = 0.000246861;
    double ST1500_p_value42 = 8.03233e-05;
    double ST1500_p_value43 = 2.3888e-05;
    double ST1500_p_value44 = 6.48441e-06;
    double ST1500_p_value45 = 1.6048e-06;
    double ST1500_p_value46 = 3.61729e-07;
    double ST1500_p_value47 = 7.41948e-08;
    double ST1500_p_value48 = 1.38371e-08;
    double ST1500_p_value49 = 2.34496e-09;
    double ST1500_p_value50 = 3.60972e-10;
    double ST1500_p_value51 = 5.04486e-11;
    
    double p_value_ST1500[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};
    
    double Tri_Combined_p_value1 = 3.73145e-11;
    double Tri_Combined_p_value2 = 2.71526e-10;
    double Tri_Combined_p_value3 = 1.79212e-09;
    double Tri_Combined_p_value4 = 1.07446e-08;
    double Tri_Combined_p_value5 = 5.85451e-08;
    double Tri_Combined_p_value6 = 2.89765e-07;
    double Tri_Combined_p_value7 = 1.30514e-06;
    double Tri_Combined_p_value8 = 5.35129e-06;
    double Tri_Combined_p_value9 = 1.99934e-05;
    double Tri_Combined_p_value10 = 6.81868e-05;
    double Tri_Combined_p_value11 = 0.000212439;
    double Tri_Combined_p_value12 = 0.000605868;
    double Tri_Combined_p_value13 = 0.00158294;
    double Tri_Combined_p_value14 = 0.00380184;
    double Tri_Combined_p_value15 = 0.0084001;
    double Tri_Combined_p_value16 = 0.017125;
    double Tri_Combined_p_value17 = 0.0323027;
    double Tri_Combined_p_value18 = 0.0565771;
    double Tri_Combined_p_value19 = 0.0923185;
    double Tri_Combined_p_value20 = 0.140789;
    double Tri_Combined_p_value21 = 0.201352;
    double Tri_Combined_p_value22 = 0.271577;
    double Tri_Combined_p_value23 = 0.345665;
    double Tri_Combined_p_value24 = 0.417059;
    double Tri_Combined_p_value25 = 0.474349;
    double Tri_Combined_p_value26 = 0.5;
    double Tri_Combined_p_value27 = 0.479134;
    double Tri_Combined_p_value28 = 0.424286;
    double Tri_Combined_p_value29 = 0.354001;
    double Tri_Combined_p_value30 = 0.27959;
    double Tri_Combined_p_value31 = 0.208917;
    double Tri_Combined_p_value32 = 0.147018;
    double Tri_Combined_p_value33 = 0.0971104;
    double Tri_Combined_p_value34 = 0.0599737;
    double Tri_Combined_p_value35 = 0.0345217;
    double Tri_Combined_p_value36 = 0.0184655;
    double Tri_Combined_p_value37 = 0.00912906;
    double Tri_Combined_p_value38 = 0.00416956;
    double Tri_Combined_p_value39 = 0.00175321;
    double Tri_Combined_p_value40 = 0.000677845;
    double Tri_Combined_p_value41 = 0.00024014;
    double Tri_Combined_p_value42 = 7.78138e-05;
    double Tri_Combined_p_value43 = 2.30617e-05;
    double Tri_Combined_p_value44 = 6.23638e-06;
    double Tri_Combined_p_value45 = 1.53612e-06;
    double Tri_Combined_p_value46 = 3.44656e-07;
    double Tri_Combined_p_value47 = 7.03791e-08;
    double Tri_Combined_p_value48 = 1.30625e-08;
    double Tri_Combined_p_value49 = 2.20106e-09;
    double Tri_Combined_p_value50 = 3.37111e-10;
    double Tri_Combined_p_value51 = 4.68588e-11;
    
    double Tri_Combined_p_value[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};
    */
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
   
    double ST_shape_ZeroSFOS[11] = {1.03951e-05, 0.000556546, 0.0107108, 0.0818906, 0.279586, 0.499397, 0.271035, 0.0761098, 0.00973483, 0.000494683, 9.03808e-06};
    
    double ST_shape_OneSFOS[11] = {9.79566e-08, 3.14837e-05, 0.00221108, 0.0400688, 0.227616, 0.499486, 0.219845, 0.0377864,  0.00196048, 2.68776e-05, 8.08226e-08};
    
    double ST_shape_TwoSFOS[11] = {3.28425e-06, 0.000262321, 0.00673895, 0.0638775, 0.256836, 0.49955, 0.252657, 0.0617097, 0.00635205, 0.000238666, 2.93993e-06};
    
    double param11[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    TGraph *gr = new TGraph(51, param, ZeroOSSF_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param11, ST_shape_ZeroSFOS);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_2 = new TGraph(11, param11, ST_shape_OneSFOS);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");
    
    TGraph *gr_3 = new TGraph(11, param11, ST_shape_TwoSFOS);
    gr_3->SetLineColor(kCyan);
    //gr_1->SetLineWidth(3);
    gr_3->Draw("LP*");
    /*
    TGraph *gr_2 = new TGraph(51, param, TwoOSSF_p_value);
    gr_2->SetLineColor(kGreen+3);
    //gr_2->SetLineWidth(3);
    gr_2->Draw("LP*");
    
    TGraph *gr_3 = new TGraph(51, param, p_value_ST1500);
    gr_3->SetLineColor(kCyan);
    //gr_2->SetLineWidth(3);
    gr_3->Draw("LP*");
    
    TGraph *gr_4 = new TGraph(51, param, Tri_Combined_p_value);
    gr_4->SetLineColor(kOrange);
    gr_4->Draw("LP*");
    */
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: 0OSSF", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 0OSSF Shape", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 1OSSF Shape", "pl");
    leg->AddEntry(gr_2,   "Trilepton: 2OSSF Shape", "pl");
    //leg->AddEntry(gr_3,   "Trilepton: All", "pl");
    //leg->AddEntry(gr_4,   "Trilepton: 3 bins", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 0OSSF = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 0OSSF = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 0OSSF shape based = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 0OSSF shape based = " << m1 << endl;
    
    double m2=-2.0;
    while (m2<=2.0 and gr_2->Eval(m2)<gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected lower limit combined 2SFOS = " << m2 << endl;
    
    while (m2<=2.0 and gr_2->Eval(m2)>gr_5->Eval(m2)) m2+=0.001;
    std::cout << "expected upper limit combined 2SFOS = " << m2 << endl;
    
    double m3=-2.0;
    while (m3<=2.0 and gr_3->Eval(m3)<gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected lower limit combined channel = " << m3 << endl;
    
    while (m3<=2.0 and gr_3->Eval(m3)>gr_5->Eval(m3)) m3+=0.001;
    std::cout << "expected upper limit combined channel = " << m3 << endl;
    /*
    double m4=-2.0;
    while (m4<=2.0 and gr_4->Eval(m4)<gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected lower limit combined channel 3 bins = " << m4 << endl;
    
    while (m4<=2.0 and gr_4->Eval(m4)>gr_5->Eval(m4)) m4+=0.001;
    std::cout << "expected upper limit combined channel 3 bins = " << m4 << endl;
    */
    c1.SaveAs("Parabola_ft0_p_value_Channelwise_3l.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Channelwise_3l.png");
}

void makeGraphs_pValue_STshape_ZeroOSSF_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double ZeroOSSF_p_value1 = 0.000386417;
    double ZeroOSSF_p_value2 = 0.00068111;
    double ZeroOSSF_p_value3 = 0.00116893;
    double ZeroOSSF_p_value4 = 0.00195411;
    double ZeroOSSF_p_value5 = 0.00318303;
    double ZeroOSSF_p_value6 = 0.00505382;
    double ZeroOSSF_p_value7 = 0.00782458;
    double ZeroOSSF_p_value8 = 0.011818;
    double ZeroOSSF_p_value9 = 0.0174202;
    double ZeroOSSF_p_value10 = 0.0250727;
    double ZeroOSSF_p_value11 = 0.0352532;
    double ZeroOSSF_p_value12 = 0.0484476;
    double ZeroOSSF_p_value13 = 0.0651127;
    double ZeroOSSF_p_value14 = 0.0856308;
    double ZeroOSSF_p_value15 = 0.110262;
    double ZeroOSSF_p_value16 = 0.139102;
    double ZeroOSSF_p_value17 = 0.172039;
    double ZeroOSSF_p_value18 = 0.208735;
    double ZeroOSSF_p_value19 = 0.248611;
    double ZeroOSSF_p_value20 = 0.29084;
    double ZeroOSSF_p_value21 = 0.334366;
    double ZeroOSSF_p_value22 = 0.377868;
    double ZeroOSSF_p_value23 = 0.419673;
    double ZeroOSSF_p_value24 = 0.457415;
    double ZeroOSSF_p_value25 = 0.486973;
    double ZeroOSSF_p_value26 = 0.5;
    double ZeroOSSF_p_value27 = 0.488883;
    double ZeroOSSF_p_value28 = 0.460331;
    double ZeroOSSF_p_value29 = 0.423079;
    double ZeroOSSF_p_value30 = 0.381515;
    double ZeroOSSF_p_value31 = 0.338087;
    double ZeroOSSF_p_value32 = 0.294508;
    double ZeroOSSF_p_value33 = 0.252123;
    double ZeroOSSF_p_value34 = 0.212012;
    double ZeroOSSF_p_value35 = 0.175019;
    double ZeroOSSF_p_value36 = 0.141746;
    double ZeroOSSF_p_value37 = 0.112551;
    double ZeroOSSF_p_value38 = 0.0875626;
    double ZeroOSSF_p_value39 = 0.066703;
    double ZeroOSSF_p_value40 = 0.0497238;
    double ZeroOSSF_p_value41 = 0.0362512;
    double ZeroOSSF_p_value42 = 0.0258332;
    double ZeroOSSF_p_value43 = 0.0179847;
    double ZeroOSSF_p_value44 = 0.0122259;
    double ZeroOSSF_p_value45 = 0.00811159;
    double ZeroOSSF_p_value46 = 0.00525035;
    double ZeroOSSF_p_value47 = 0.00331394;
    double ZeroOSSF_p_value48 = 0.00203894;
    double ZeroOSSF_p_value49 = 0.00122239;
    double ZeroOSSF_p_value50 = 0.00071383;
    double ZeroOSSF_p_value51 = 0.000405919;
    
    double ZeroOSSF_p_value[51] = {ZeroOSSF_p_value1, ZeroOSSF_p_value2, ZeroOSSF_p_value3, ZeroOSSF_p_value4, ZeroOSSF_p_value5, ZeroOSSF_p_value6, ZeroOSSF_p_value7, ZeroOSSF_p_value8, ZeroOSSF_p_value9, ZeroOSSF_p_value10, ZeroOSSF_p_value11, ZeroOSSF_p_value12, ZeroOSSF_p_value13, ZeroOSSF_p_value14, ZeroOSSF_p_value15, ZeroOSSF_p_value16, ZeroOSSF_p_value17, ZeroOSSF_p_value18, ZeroOSSF_p_value19, ZeroOSSF_p_value20, ZeroOSSF_p_value21, ZeroOSSF_p_value22, ZeroOSSF_p_value23, ZeroOSSF_p_value24, ZeroOSSF_p_value25, ZeroOSSF_p_value26, ZeroOSSF_p_value27, ZeroOSSF_p_value28, ZeroOSSF_p_value29, ZeroOSSF_p_value30, ZeroOSSF_p_value31, ZeroOSSF_p_value32, ZeroOSSF_p_value33, ZeroOSSF_p_value34, ZeroOSSF_p_value35, ZeroOSSF_p_value36, ZeroOSSF_p_value37, ZeroOSSF_p_value38, ZeroOSSF_p_value39, ZeroOSSF_p_value40, ZeroOSSF_p_value41, ZeroOSSF_p_value42, ZeroOSSF_p_value43, ZeroOSSF_p_value44, ZeroOSSF_p_value45, ZeroOSSF_p_value46, ZeroOSSF_p_value47, ZeroOSSF_p_value48, ZeroOSSF_p_value49, ZeroOSSF_p_value50, ZeroOSSF_p_value51};
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    double ST_shape_ZeroSFOS[11] = {1.03951e-05, 0.000556546, 0.0107108, 0.0818906, 0.279586, 0.499397, 0.271035, 0.0761098, 0.00973483, 0.000494683, 9.03808e-06};
    
    //double ST_shape_OneSFOS[11] = {9.79566e-08, 3.14837e-05, 0.00221108, 0.0400688, 0.227616, 0.499486, 0.219845, 0.0377864,  0.00196048, 2.68776e-05, 8.08226e-08};
    
    //double ST_shape_TwoSFOS[11] = {3.28425e-06, 0.000262321, 0.00673895, 0.0638775, 0.256836, 0.49955, 0.252657, 0.0617097, 0.00635205, 0.000238666, 2.93993e-06};
    
    double param11[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    TGraph *gr = new TGraph(51, param, ZeroOSSF_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param11, ST_shape_ZeroSFOS);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: 0OSSF", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 0OSSF Shape", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 0OSSF = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 0OSSF = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 0OSSF shape based = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 0OSSF shape based = " << m1 << endl;
    
    c1.SaveAs("Parabola_ft0_p_value_Shape_0OSSF.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Shape_0OSSF.png");
}

void makeGraphs_pValue_STshape_OneOSSF_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double OneOSSF_p_value1 = 1.74308e-05;
    double OneOSSF_p_value2 = 4.03969e-05;
    double OneOSSF_p_value3 = 8.98668e-05;
    double OneOSSF_p_value4 = 0.000192015;
    double OneOSSF_p_value5 = 0.000394194;
    double OneOSSF_p_value6 = 0.000777979;
    double OneOSSF_p_value7 = 0.00147681;
    double OneOSSF_p_value8 = 0.00269816;
    double OneOSSF_p_value9 = 0.0047474;
    double OneOSSF_p_value10 = 0.00804998;
    double OneOSSF_p_value11 = 0.0131647;
    double OneOSSF_p_value12 = 0.0207801;
    double OneOSSF_p_value13 = 0.0316873;
    double OneOSSF_p_value14 = 0.0467217;
    double OneOSSF_p_value15 = 0.0666768;
    double OneOSSF_p_value16 = 0.0921952;
    double OneOSSF_p_value17 = 0.12365;
    double OneOSSF_p_value18 = 0.161036;
    double OneOSSF_p_value19 = 0.203895;
    double OneOSSF_p_value20 = 0.25127;
    double OneOSSF_p_value21 = 0.301707;
    double OneOSSF_p_value22 = 0.353259;
    double OneOSSF_p_value23 = 0.403433;
    double OneOSSF_p_value24 = 0.448875;
    double OneOSSF_p_value25 = 0.484225;
    double OneOSSF_p_value26 = 0.5;
    double OneOSSF_p_value27 = 0.488011;
    double OneOSSF_p_value28 = 0.454796;
    double OneOSSF_p_value29 = 0.410415;
    double OneOSSF_p_value30 = 0.36069;
    double OneOSSF_p_value31 = 0.309167;
    double OneOSSF_p_value32 = 0.258437;
    double OneOSSF_p_value33 = 0.210518;
    double OneOSSF_p_value34 = 0.166932;
    double OneOSSF_p_value35 = 0.128711;
    double OneOSSF_p_value36 = 0.0963855;
    double OneOSSF_p_value37 = 0.0700214;
    double OneOSSF_p_value38 = 0.0492941;
    double OneOSSF_p_value39 = 0.0335928;
    double OneOSSF_p_value40 = 0.0221388;
    double OneOSSF_p_value41 = 0.0140967;
    double OneOSSF_p_value42 = 0.00866471;
    double OneOSSF_p_value43 = 0.00513707;
    double OneOSSF_p_value44 = 0.00293545;
    double OneOSSF_p_value45 = 0.00161557;
    double OneOSSF_p_value46 = 0.000855816;
    double OneOSSF_p_value47 = 0.000436107;
    double OneOSSF_p_value48 = 0.00021366;
    double OneOSSF_p_value49 = 0.000100582;
    double OneOSSF_p_value50 = 4.54787e-05;
    double OneOSSF_p_value51 = 1.97421e-05;
    
    double OneOSSF_p_value[51] = {OneOSSF_p_value1, OneOSSF_p_value2, OneOSSF_p_value3, OneOSSF_p_value4, OneOSSF_p_value5, OneOSSF_p_value6, OneOSSF_p_value7, OneOSSF_p_value8, OneOSSF_p_value9, OneOSSF_p_value10, OneOSSF_p_value11, OneOSSF_p_value12, OneOSSF_p_value13, OneOSSF_p_value14, OneOSSF_p_value15, OneOSSF_p_value16, OneOSSF_p_value17, OneOSSF_p_value18, OneOSSF_p_value19, OneOSSF_p_value20, OneOSSF_p_value21, OneOSSF_p_value22, OneOSSF_p_value23, OneOSSF_p_value24, OneOSSF_p_value25, OneOSSF_p_value26, OneOSSF_p_value27, OneOSSF_p_value28, OneOSSF_p_value29, OneOSSF_p_value30, OneOSSF_p_value31, OneOSSF_p_value32, OneOSSF_p_value33, OneOSSF_p_value34, OneOSSF_p_value35, OneOSSF_p_value36, OneOSSF_p_value37, OneOSSF_p_value38, OneOSSF_p_value39, OneOSSF_p_value40, OneOSSF_p_value41, OneOSSF_p_value42, OneOSSF_p_value43, OneOSSF_p_value44, OneOSSF_p_value45, OneOSSF_p_value46, OneOSSF_p_value47, OneOSSF_p_value48, OneOSSF_p_value49, OneOSSF_p_value50, OneOSSF_p_value51};
    
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    double ST_shape_OneSFOS[11] = {9.79566e-08, 3.14837e-05, 0.00221108, 0.0400688, 0.227616, 0.499486, 0.219845, 0.0377864,  0.00196048, 2.68776e-05, 8.08226e-08};
    
    //double ST_shape_TwoSFOS[11] = {3.28425e-06, 0.000262321, 0.00673895, 0.0638775, 0.256836, 0.49955, 0.252657, 0.0617097, 0.00635205, 0.000238666, 2.93993e-06};
    
    double param11[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    TGraph *gr = new TGraph(51, param, OneOSSF_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param11, ST_shape_OneSFOS);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: 1OSSF", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 1OSSF Shape", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 1OSSF = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 1OSSF = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 1OSSF shape based = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 1OSSF shape based = " << m1 << endl;
    
    c1.SaveAs("Parabola_ft0_p_value_Shape_1OSSF.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Shape_1OSSF.png");
}

void makeGraphs_pValue_STshape_TwoOSSF_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double TwoOSSF_p_value1 = 5.81205e-05;
    double TwoOSSF_p_value2 = 0.000120464;
    double TwoOSSF_p_value3 = 0.000241147;
    double TwoOSSF_p_value4 = 0.00046642;
    double TwoOSSF_p_value5 = 0.000872007;
    double TwoOSSF_p_value6 = 0.00157653;
    double TwoOSSF_p_value7 = 0.00275766;
    double TwoOSSF_p_value8 = 0.00466934;
    double TwoOSSF_p_value9 = 0.00765739;
    double TwoOSSF_p_value10 = 0.0121697;
    double TwoOSSF_p_value11 = 0.018755;
    double TwoOSSF_p_value12 = 0.0280475;
    double TwoOSSF_p_value13 = 0.0407302;
    double TwoOSSF_p_value14 = 0.0574802;
    double TwoOSSF_p_value15 = 0.0788956;
    double TwoOSSF_p_value16 = 0.105412;
    double TwoOSSF_p_value17 = 0.13722;
    double TwoOSSF_p_value18 = 0.174195;
    double TwoOSSF_p_value19 = 0.215848;
    double TwoOSSF_p_value20 = 0.261306;
    double TwoOSSF_p_value21 = 0.309317;
    double TwoOSSF_p_value22 = 0.358235;
    double TwoOSSF_p_value23 = 0.405971;
    double TwoOSSF_p_value24 = 0.449634;
    double TwoOSSF_p_value25 = 0.484342;
    double TwoOSSF_p_value26 = 0.5;
    double TwoOSSF_p_value27 = 0.48686;
    double TwoOSSF_p_value28 = 0.453416;
    double TwoOSSF_p_value29 = 0.41034;
    double TwoOSSF_p_value30 = 0.362851;
    double TwoOSSF_p_value31 = 0.313944;
    double TwoOSSF_p_value32 = 0.265773;
    double TwoOSSF_p_value33 = 0.220013;
    double TwoOSSF_p_value34 = 0.177954;
    double TwoOSSF_p_value35 = 0.140509;
    double TwoOSSF_p_value36 = 0.108201;
    double TwoOSSF_p_value37 = 0.0811863;
    double TwoOSSF_p_value38 = 0.0593025;
    double TwoOSSF_p_value39 = 0.0421335;
    double TwoOSSF_p_value40 = 0.0290936;
    double TwoOSSF_p_value41 = 0.0195095;
    double TwoOSSF_p_value42 = 0.0126957;
    double TwoOSSF_p_value43 = 0.00801199;
    double TwoOSSF_p_value44 = 0.00490024;
    double TwoOSSF_p_value45 = 0.00290288;
    double TwoOSSF_p_value46 = 0.00166472;
    double TwoOSSF_p_value47 = 0.000923654;
    double TwoOSSF_p_value48 = 0.000495629;
    double TwoOSSF_p_value49 = 0.000257081;
    double TwoOSSF_p_value50 = 0.000128845;
    double TwoOSSF_p_value51 = 6.23709e-05;
    
    double TwoOSSF_p_value[51] = {TwoOSSF_p_value1, TwoOSSF_p_value2, TwoOSSF_p_value3, TwoOSSF_p_value4, TwoOSSF_p_value5, TwoOSSF_p_value6, TwoOSSF_p_value7, TwoOSSF_p_value8, TwoOSSF_p_value9, TwoOSSF_p_value10, TwoOSSF_p_value11, TwoOSSF_p_value12, TwoOSSF_p_value13, TwoOSSF_p_value14, TwoOSSF_p_value15, TwoOSSF_p_value16, TwoOSSF_p_value17, TwoOSSF_p_value18, TwoOSSF_p_value19, TwoOSSF_p_value20, TwoOSSF_p_value21, TwoOSSF_p_value22, TwoOSSF_p_value23, TwoOSSF_p_value24, TwoOSSF_p_value25, TwoOSSF_p_value26, TwoOSSF_p_value27, TwoOSSF_p_value28, TwoOSSF_p_value29, TwoOSSF_p_value30, TwoOSSF_p_value31, TwoOSSF_p_value32, TwoOSSF_p_value33, TwoOSSF_p_value34, TwoOSSF_p_value35, TwoOSSF_p_value36, TwoOSSF_p_value37, TwoOSSF_p_value38, TwoOSSF_p_value39, TwoOSSF_p_value40, TwoOSSF_p_value41, TwoOSSF_p_value42, TwoOSSF_p_value43, TwoOSSF_p_value44, TwoOSSF_p_value45, TwoOSSF_p_value46, TwoOSSF_p_value47, TwoOSSF_p_value48, TwoOSSF_p_value49, TwoOSSF_p_value50, TwoOSSF_p_value51};

    
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
    
    double ST_shape_TwoSFOS[11] = {3.28425e-06, 0.000262321, 0.00673895, 0.0638775, 0.256836, 0.49955, 0.252657, 0.0617097, 0.00635205, 0.000238666, 2.93993e-06};
    
    double param11[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    TGraph *gr = new TGraph(51, param, TwoOSSF_p_value);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param11, ST_shape_TwoSFOS);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: 2OSSF", "pl");
    leg->AddEntry(gr_1,   "Trilepton: 2OSSF Shape", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined 2OSSF = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined 2OSSF = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined 2OSSF shape based = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined 2OSSF shape based = " << m1 << endl;
    
    c1.SaveAs("Parabola_ft0_p_value_Shape_2OSSF.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Shape_2OSSF.png");
}

void makeGraphs_pValue_STshape_Comb_3l()
{
    
    gROOT->SetStyle("Plain");
    TCanvas c1("c1","Cleaning Plot", 10, 10, 1200, 800);
    c1.SetLogy();
    
    double param[51] = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
    
    double ST1500_p_value1 = 4.00434e-11;
    double ST1500_p_value2 = 2.89704e-10;
    double ST1500_p_value3 = 1.90293e-09;
    double ST1500_p_value4 = 1.13526e-08;
    double ST1500_p_value5 = 6.15367e-08;
    double ST1500_p_value6 = 3.03298e-07;
    double ST1500_p_value7 = 1.3601e-06;
    double ST1500_p_value8 = 5.55436e-06;
    double ST1500_p_value9 = 2.06777e-05;
    double ST1500_p_value10 = 7.02526e-05;
    double ST1500_p_value11 = 0.000218133;
    double ST1500_p_value12 = 0.000619906;
    double ST1500_p_value13 = 0.0016153;
    double ST1500_p_value14 = 0.0038669;
    double ST1500_p_value15 = 0.00852296;
    double ST1500_p_value16 = 0.0173394;
    double ST1500_p_value17 = 0.0326507;
    double ST1500_p_value18 = 0.0570807;
    double ST1500_p_value19 = 0.0929579;
    double ST1500_p_value20 = 0.141534;
    double ST1500_p_value21 = 0.202242;
    double ST1500_p_value22 = 0.272258;
    double ST1500_p_value23 = 0.346448;
    double ST1500_p_value24 = 0.41741;
    double ST1500_p_value25 = 0.474421;
    double ST1500_p_value26 = 0.5;
    double ST1500_p_value27 = 0.47934;
    double ST1500_p_value28 = 0.42488;
    double ST1500_p_value29 = 0.354882;
    double ST1500_p_value30 = 0.28064;
    double ST1500_p_value31 = 0.209845;
    double ST1500_p_value32 = 0.147887;
    double ST1500_p_value33 = 0.0978567;
    double ST1500_p_value34 = 0.0605638;
    double ST1500_p_value35 = 0.034931;
    double ST1500_p_value36 = 0.0187115;
    double ST1500_p_value37 = 0.00928034;
    double ST1500_p_value38 = 0.00424965;
    double ST1500_p_value39 = 0.00179231;
    double ST1500_p_value40 = 0.000694601;
    double ST1500_p_value41 = 0.000246861;
    double ST1500_p_value42 = 8.03233e-05;
    double ST1500_p_value43 = 2.3888e-05;
    double ST1500_p_value44 = 6.48441e-06;
    double ST1500_p_value45 = 1.6048e-06;
    double ST1500_p_value46 = 3.61729e-07;
    double ST1500_p_value47 = 7.41948e-08;
    double ST1500_p_value48 = 1.38371e-08;
    double ST1500_p_value49 = 2.34496e-09;
    double ST1500_p_value50 = 3.60972e-10;
    double ST1500_p_value51 = 5.04486e-11;
    
    double p_value_ST1500[51] = {ST1500_p_value1, ST1500_p_value2, ST1500_p_value3, ST1500_p_value4, ST1500_p_value5, ST1500_p_value6, ST1500_p_value7, ST1500_p_value8, ST1500_p_value9, ST1500_p_value10, ST1500_p_value11, ST1500_p_value12, ST1500_p_value13, ST1500_p_value14, ST1500_p_value15, ST1500_p_value16, ST1500_p_value17, ST1500_p_value18, ST1500_p_value19, ST1500_p_value20, ST1500_p_value21, ST1500_p_value22, ST1500_p_value23, ST1500_p_value24, ST1500_p_value25, ST1500_p_value26, ST1500_p_value27, ST1500_p_value28, ST1500_p_value29, ST1500_p_value30, ST1500_p_value31, ST1500_p_value32, ST1500_p_value33, ST1500_p_value34, ST1500_p_value35, ST1500_p_value36, ST1500_p_value37, ST1500_p_value38, ST1500_p_value39, ST1500_p_value40, ST1500_p_value41, ST1500_p_value42, ST1500_p_value43, ST1500_p_value44, ST1500_p_value45, ST1500_p_value46, ST1500_p_value47, ST1500_p_value48, ST1500_p_value49, ST1500_p_value50, ST1500_p_value51};
    
    
    double p_value_disc[51] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

    
    double ST_shape_Combined[11] = {2.93156e-16, 2.5247e-10, 5.03892e-06, 0.00340927, 0.124823, 0.499089, 0.124823, 0.00340927, 5.03892e-06, 2.5247e-10, 2.93156e-16};
    
    double param11[11] = {-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
    
    TGraph *gr = new TGraph(51, param, p_value_ST1500);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->SetLineColor(kBlue);
    gr->SetMarkerSize(0.9);
    gr->SetMaximum(1.0);
    gr->SetMinimum(0.02);
    //gr->SetMinimum(0.000002);
    gr->SetTitle("");
    gr->SetTitle("p value vs Param (ft0)");
    gr->GetXaxis()->SetTitle("Param (ft0)");
    gr->GetYaxis()->SetTitle("p-value");
    gr->Draw("ALP*");
    
    TGraph *gr_1 = new TGraph(11, param11, ST_shape_Combined);
    gr_1->SetLineColor(kBlack);
    //gr_1->SetLineWidth(3);
    gr_1->Draw("LP*");
    
    TGraph *gr_5 = new TGraph(51, param, p_value_disc);
    gr_5->SetLineColor(kRed);
    gr_5->SetLineWidth(3);
    gr_5->Draw("L");
    
    TLegend *leg = new TLegend(0.60,0.65,0.80,0.89,NULL,"brNDC");
    leg->AddEntry(gr,   "Trilepton: Combined", "pl");
    leg->AddEntry(gr_1,   "Trilepton: Combined Shape", "pl");
    leg->Draw();
    
    double m=-2.0;
    while (m<=2.0 and gr->Eval(m)<gr_5->Eval(m)) m+=0.001;
    std::cout << "expected lower limit combined Comb = " << m << endl;
    
    while (m<=2.0 and gr->Eval(m)>gr_5->Eval(m)) m+=0.001;
    std::cout << "expected upper limit combined Comb = " << m << endl;
    
    double m1=-2.0;
    while (m1<=2.0 and gr_1->Eval(m1)<gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected lower limit combined Comb shape based = " << m1 << endl;
    
    while (m1<=2.0 and gr_1->Eval(m1)>gr_5->Eval(m1)) m1+=0.001;
    std::cout << "expected upper limit combined Comb shape based = " << m1 << endl;
    
    c1.SaveAs("Parabola_ft0_p_value_Shape_Comb.pdf");
    c1.SaveAs("Parabola_ft0_p_value_Shape_Comb.png");
}

