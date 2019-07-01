import sys
import ROOT as rt
import math
from LHEevent import *
from LHEfile import *
import plotTools

if __name__ == '__main__':

    #Bprime histograms
    MW_jj = rt.TH1D("MW_jj", "MW_jj", 500, 0., 500)
    MW_jj.Sumw2()
    MInvariantMass_mumu = rt.TH1F("MInvariantMass_mumu", "MInvariantMass_mumu", 500, 0., 500);
    MInvariantMass_mumu.Sumw2()    
    MInvariantMass_qq = rt.TH1F("MInvariantMass_qq", "MInvariantMass_qq", 500, 0., 1000.0);
    MInvariantMass_qq.Sumw2()
    EW_jj = rt.TH1D("EW_jj", "EW_jj", 500, 0., 1000)
    EW_jj.Sumw2()
    EW_qq = rt.TH1D("EW_qq", "EW_qq", 500, 0., 1000)
    EW_qq.Sumw2()
    pTW_jj = rt.TH1D("pTW_jj", "pTW_jj", 500, 0., 500)
    pTW_jj.Sumw2()
    gamma_qq = rt.TH1D("gamma_qq", "gamma_qq", 100, 0.0, 100.)
    gamma_qq.Sumw2()
    DeltaR_qq = rt.TH1D("DeltaR_qq", "DeltaR_qq", 1000, 0.0, 10.0)
    DeltaR_qq.Sumw2()
    DeltaR_qq_pT = rt.TH2D("DeltaR_qq_pT", "DeltaR_qq_pT", 1000, 0.0, 10.0, 5000, 0.0, 500.0)
    DeltaR_qq_pT.Sumw2()
    DeltaR_qq_pT_Wmass = rt.TH2D("DeltaR_qq_pT_Wmass", "DeltaR_qq_pT_Wmass", 1000, 0.0, 10.0, 5000, 0.0, 500.0)
    DeltaR_qq_pT_Wmass.Sumw2()
    h_boostedpT_deltaR04 = rt.TH1D("h_boostedpT_deltaR04", "h_boostedpT_deltaR04", 5000, 0.0, 5000.0)
    h_boostedpT_deltaR04.Sumw2()
    h_boostedDeltaR_pT150 = rt.TH1D("h_boostedDeltaR_pT150", "h_boostedDeltaR_pT150", 500, 0.0, 5.0)
    h_boostedDeltaR_pT150.Sumw2()
    mu1_lv = rt.TLorentzVector()
    mu2_lv = rt.TLorentzVector()
    q1_lv = rt.TLorentzVector()
    q2_lv = rt.TLorentzVector()
    

    # find events in file
    myLHEfile = LHEfile(sys.argv[1])
    myLHEfile.setMax(100000)
    #myLHEfile.setMax(1)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_mu = 0
        n_q = 0
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            if abs(p['ID'])  == 24: MW_jj.Fill(p['M'])
            if abs(p['ID'])  == 24: EW_jj.Fill(p['E'])
            if (abs(p['ID'])  == 24 and rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']) > 50.0): pTW_jj.Fill(rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']))
            if abs(p['ID']) == 13: 
              n_mu += 1
              if n_mu==1:  mu1_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
              if n_mu==2:  mu2_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
              if n_mu==2:  MInvariantMass_mumu.Fill((mu1_lv+mu2_lv).M())
            if ((abs(p['ID']) == 1 or abs(p['ID']) == 2 or abs(p['ID']) == 3 or abs(p['ID']) == 4)):
              n_q += 1
              print n_q
              #print p['Px'], p['Py'], p['Pz'], p['E']
              if n_q==3: q1_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
              if n_q==4: q2_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
              if n_q==4: MInvariantMass_qq.Fill((q1_lv+q2_lv).M())
              if n_q==4: EW_qq.Fill((q1_lv+q2_lv).E());
              if n_q==4: DeltaR_qq.Fill(q1_lv.DeltaR(q2_lv)); #q1_lv.Pt() > 180.0 and q2_lv.Pt() > 180.0): DeltaR_qq.Fill(q1_lv.DeltaR(q2_lv));
              if n_q==4: DeltaR_qq_pT.Fill(q1_lv.DeltaR(q2_lv), (q1_lv.Pt()+q2_lv.Pt())/2.0)
              if n_q==4 and ((q1_lv+q2_lv).M() > 70.0 and (q1_lv+q2_lv).M() < 90.0): DeltaR_qq_pT_Wmass.Fill(q1_lv.DeltaR(q2_lv), (q1_lv.Pt()+q2_lv.Pt())/2.0) 
              if n_q==4: h_boostedpT_deltaR04.Fill((q1_lv.Pt()+q2_lv.Pt()))
              #if n_q==4 and q1_lv.DeltaR(q2_lv) < 0.4: h_boostedpT_deltaR04.Fill((q1_lv.Pt()+q2_lv.Pt())/2.0)
              if n_q==4 and (q1_lv.Pt()+q2_lv.Pt())/2.0 > 150.0: h_boostedDeltaR_pT150.Fill(q1_lv.DeltaR(q2_lv)) 
              if (q1_lv+q2_lv).M() > 0.0: gamma_qq.Fill((q1_lv+q2_lv).E()/(q1_lv+q2_lv).M());
        del oneEvent, myLHEevent
        
    # write the histograms
    histoFILE = rt.TFile(sys.argv[2],"RECREATE")
    MW_jj.Write()
    EW_jj.Write()
    MInvariantMass_mumu.Write();
    MInvariantMass_qq.Write();
    pTW_jj.Write()
    EW_qq.Write()
    gamma_qq.Write()
    DeltaR_qq.Write()
    DeltaR_qq_pT.Write()
    DeltaR_qq_pT_Wmass.Write()
    h_boostedpT_deltaR04.Write()
    h_boostedDeltaR_pT150.Write()
    histoFILE.Close()
