#!/usr/bin/env python
# encoding: utf-8

# File        : JetStudy.py
# Author      : Ben Wu
# Contact     : benwu@fnal.gov
# Date        : 2018 Oct 26
#
# Description :

from rootpy.plotting import Hist,Canvas
from rootpy.interactive import wait
from rootpy.plotting.style import get_style
import numpy as np

class JetStudy():
    def __init__ (self, jettype):
        self.hismap = {}
        self.BookHistogram()
        self.Jets = []
        self.jettype = jettype

    def BookHistogram(self):
        self.hismap["JetPt"]=     Hist(100, 0, 1000, title="JetPt;Jet Pt;Event")
        self.hismap["JetEta"]=     Hist(20, -5, 5, title="JetEta;Jet Eta;Event")
        self.hismap["JetPhi"]=     Hist(20, -5, 5, title="JetPhi;Jet Phi;Event")
        self.hismap["JetMass"]=   Hist(100, 0, 200,  title="JetMass;Jet Mass;Event")
        self.hismap["nJet"]=     Hist(100, 0, 100, title="NJets;NJets;Event")
        self.hismap["nJet30"]=     Hist(100, 0, 100, title="NJets30;NJets30;Event")

    def Run(self, jetmap):
        jets = np.asarray(jetmap[self.jettype])
        self.hismap["nJet"].Fill(jets.shape[0])
        self.hismap["nJet30"].Fill(np.count_nonzero(jets[:, 0] > 30))
        for j in jets:
            self.hismap["JetPt"].Fill(j[0])
            self.hismap["JetEta"].Fill(j[1])
            self.hismap["JetPhi"].Fill(j[2])
            self.hismap["JetMass"].Fill(j[3])


    def GetHist(self):
        return self.hismap


    def Draw(self):
        pass
        get_style('ATLAS')
        canvas = Canvas(width=700, height=500)
        canvas.SetLeftMargin(0.15)
        canvas.SetBottomMargin(0.15)
        canvas.SetTopMargin(0.10)
        canvas.SetRightMargin(0.05)
        for k,v in self.hismap.items():
            canvas.Clear()
            v.Draw()
            wait()

