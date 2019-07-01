import sys
import ROOT as rt
import math
from LHEevent import *
from LHEfile import *
import plotTools

if __name__ == '__main__':

    #Bprime histograms
    h_Lepton_1_Pt = rt.TH1F("h_Lepton_1_Pt", "h_Lepton_1_Pt", 1000, 0., 1000)
    h_Lepton_1_Pt.Sumw2()
    h_Lepton_2_Pt = rt.TH1F("h_Lepton_2_Pt", "h_Lepton_2_Pt", 1000, 0., 1000)
    h_Lepton_2_Pt.Sumw2()
    h_Lepton_3_Pt = rt.TH1F("h_Lepton_3_Pt", "h_Lepton_3_Pt", 1000, 0., 1000)
    h_Lepton_3_Pt.Sumw2()
    h_Lepton_4_Pt = rt.TH1F("h_Lepton_4_Pt", "h_Lepton_4_Pt", 1000, 0., 1000)
    h_Lepton_4_Pt.Sumw2()
    h_Lepton_5_Pt = rt.TH1F("h_Lepton_5_Pt", "h_Lepton_5_Pt", 1000, 0., 1000)
    h_Lepton_5_Pt.Sumw2()
    h_MET = rt.TH1F("h_MET", "h_MET", 1000, 0., 1000)
    h_MET.Sumw2()

    # find events in file
    myLHEfile = LHEfile(sys.argv[1])
    myLHEfile.setMax(100000)
    #myLHEfile.setMax(1)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_l = 0
        transverseMomentum = []
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            if (abs(p['ID']) == 11 or abs(p['ID']) == 13 or abs(p['ID']) == 15): 
              transverseMomentum.append(rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']))
              transverseMomentum.sort(reverse=True)
              #print transverseMomentum
              #print p['Px']
              #print p['Py']
              n_l += 1
              if n_l==5:  print n_l
              if n_l==5:  h_Lepton_1_Pt.Fill(transverseMomentum[0])
              if n_l==5:  h_Lepton_2_Pt.Fill(transverseMomentum[1])
              if n_l==5:  h_Lepton_3_Pt.Fill(transverseMomentum[2])
              if n_l==5:  h_Lepton_4_Pt.Fill(transverseMomentum[3])
              if n_l==5:  h_Lepton_5_Pt.Fill(transverseMomentum[4])
            if (abs(p['ID']) == 12 or abs(p['ID']) == 14 or abs(p['ID']) == 16):
              h_MET.Fill(rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']));
        del oneEvent, myLHEevent
        
    # write the histograms
    histoFILE = rt.TFile(sys.argv[2],"RECREATE")
    h_Lepton_1_Pt.Write()
    h_Lepton_2_Pt.Write()
    h_Lepton_3_Pt.Write()
    h_Lepton_4_Pt.Write()
    h_Lepton_5_Pt.Write()
    h_MET.Write()
    histoFILE.Close()
