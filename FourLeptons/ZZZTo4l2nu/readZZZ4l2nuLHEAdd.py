import sys
import ROOT as rt
import math
from LHEevent import *
from LHEfile import *
import plotTools
#from dataclasses import dataclass
from ROOT import TFile, TTree
from array import array
#from rootpy import stl

"""
class LV:
    px: double
    py: double
    pz: double
    e: double
"""

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
    h_MET = rt.TH1F("h_MET", "h_MET", 1000, 0., 1000)
    h_MET.Sumw2()
    nu1_lv = rt.TLorentzVector()
    nu2_lv = rt.TLorentzVector()
    f = TFile('test_tree.root', 'recreate')
    t = TTree('t_lv', 'tree')
    
    maxn = 4
    n = array( 'i', [ 0 ] )
    t.Branch('nlep', n, 'nlep/I')
    p_x = vector('float')
    p_y = vector('float')
    p_z = vector('float')
    energy = vector('float')
    t.Branch('px', p_x, 'px')
    t.Branch('py', p_y, 'py')
    t.Branch('pz', p_z, 'pz')
    t.Branch('e',  energy, 'e')

    # find events in file
    myLHEfile = LHEfile(sys.argv[1])
    #myLHEfile.setMax(100000)
    myLHEfile.setMax(1)
    eventsReadIn = myLHEfile.readEvents()
    for oneEvent in eventsReadIn:
        myLHEevent = LHEevent()
        myLHEevent.fillEvent(oneEvent)
        n_l = 0
        nu = 0
        n = 0
        p_x = p_y = p_z = energy = 0
        transverseMomentum = []
        fourvector = list()
        for i in range(0,len(myLHEevent.Particles)):
            p = myLHEevent.Particles[i]
            if (abs(p['ID']) == 11 or abs(p['ID']) == 13 or abs(p['ID']) == 15): 
              transverseMomentum.append(rt.TMath.Sqrt(p['Px']*p['Px'] + p['Py']*p['Py']))
              transverseMomentum.sort(reverse=True)
              #print transverseMomentum
              n_l += 1
              n = n_l
              print p['Px'] 
              p_x.push_back(p['Px'])
              p_y.push_back(p['Py'])
              p_z.push_back(p['Pz'])
              energy.push_back(p['E'])
              if n_l==4:  print n_l
              if n_l==4:  h_Lepton_1_Pt.Fill(transverseMomentum[0])
              if n_l==4:  h_Lepton_2_Pt.Fill(transverseMomentum[1])
              if n_l==4:  h_Lepton_3_Pt.Fill(transverseMomentum[2])
              if n_l==4:  h_Lepton_4_Pt.Fill(transverseMomentum[3])
            if (p['ID'] == 12 or p['ID'] == 14 or p['ID'] == 16): nu1_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
            if (p['ID'] == -12 or p['ID'] == -14 or p['ID'] == -16): nu2_lv = rt.TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E'])
            if (abs(p['ID']) == 12 or abs(p['ID']) == 14 or abs(p['ID']) == 16): nu += 1
            if (nu==2): h_MET.Fill((nu1_lv+nu2_lv).Pt())
            t.Fill()    
        #del oneEvent, myLHEevent
        
    # write the histograms
    f.Write()
    f.Close()
    histoFILE = rt.TFile(sys.argv[2],"RECREATE")
    h_Lepton_1_Pt.Write()
    h_Lepton_2_Pt.Write()
    h_Lepton_3_Pt.Write()
    h_Lepton_4_Pt.Write()
    h_MET.Write()
    histoFILE.Close()
