#!/usr/bin/env python
# encoding: utf-8

# File        : L1JetAlgo.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Oct 25
#
# Description :


import uproot
import numpy as np
import itertools
import ROOT
from pyjet import cluster
import math
from rootpy.plotting import Hist,Canvas
from rootpy.plotting.style import get_style
import ROOT
from JetStudy import JetStudy
from rootpy.io import root_open

IsTest=False


def deltaPhi(x, y):
    dphi = x-y
    if dphi >= math.pi:
        dphi -= 2*math.pi
    if dphi < -1*math.pi:
        dphi += 2*math.pi
    return dphi

def CluMetric(x, y, k, R=0.4):
    minv = min( pow(x[0], k * 2), pow(y[0], k * 2))
    deta = pow(x[1]-y[1], 2)
    dphi = pow(deltaPhi(x[2], y[2]), 2)
    return  minv * (deta+dphi) / R**2

def CludR(x, y):
    deta = pow(x[1]-y[1], 2)
    dphi = pow(deltaPhi(x[2], y[2]), 2)
    return (deta+dphi)

def GetSeparationList(idx, pars, k, R):
    combs = list(itertools.combinations(idx, 2))
    combm = {}
    for i in combs:
        first = list(i)[0]
        second = list(i)[1]
        combm[i]  = CluMetric(list(pars[first]), list(pars[second]), k, R)
    for i in idx:
        combm[(i, i)] = list(pars[i])[0]**(2*k)
    return combm


def GetDistanceList(idx, pars, R):
    ptlist ={list(pars[x])[0] : x for x in idx}
    if len(ptlist.keys()) == 0:
        return {}
    maxidx = ptlist[max(ptlist.keys())]
    combm = {}

    for i in idx:
        if i == maxidx:
            combm[(i, i)] = R**2
        else:
            combm[(maxidx, i)] = CludR(list(pars[maxidx]), list(pars[i]))
    return combm


def EMerging(x, y):
    nX = ROOT.TLorentzVector(0, 0, 0, 0)
    nY = ROOT.TLorentzVector(0, 0, 0, 0)
    nX.SetPtEtaPhiM(*x)
    nY.SetPtEtaPhiM(*y)
    newZ = nX + nY
    return np.array([newZ.Pt(), newZ.Eta(), newZ.Phi(), newZ.M()])

def WTAMerging(x, y):
    newX = []
    if x[0] >= y[0]:
        newX= [x[0]+y[0], x[1], x[2], 0]
    else:
        newX= [x[0]+y[0], y[1], y[2], 0]
    return np.array(newX)

class ClusAK4:
    def __init__(self, pars, cluMed="metric", merMed="EMerge"):
        self.pars = pars
        self.k =  -1
        self.R = 0.4
        self.ObjList = range(0, self.pars.shape[0])
        self.JetList = []
        self.cluMed = cluMed
        self.merMed = merMed
        if self.cluMed  is "metric":
            self.SepList = GetSeparationList(self.ObjList, self.pars, self.k, self.R)
        else:
            self.SepList = GetDistanceList(self.ObjList, self.pars, self.R)
        # print("Dine ini")

    def GetMin(self):
        m = min(self.SepList, key=self.SepList.get)
        # print(m)
        if m[0] == m[1]:
            self.MoveSingleToJet(m[0])
        else:
            self.MergingObjs(m)

    def MoveSingleToJet(self, objIdx):
        self.JetList.append(self.pars[objIdx])
        self.ObjList = [x for x in self.ObjList if x is not objIdx]
        self.UpdateLists()

    def MergingObjs(self, m):
        if self.merMed is "EMerge":
            newobj = EMerging(self.pars[m[0]], self.pars[m[1]])
            # print("newobj", newobj)
        else:
            newobj = WTAMerging(self.pars[m[0]], self.pars[m[1]])
        self.pars=np.concatenate((self.pars, np.array([newobj])), axis=0)
        self.ObjList = [x for x in self.ObjList if x not in m]
        newIdx =len(self.pars)-1
        self.ObjList.append(newIdx)
        self.UpdateLists()
        # print(self.JetList)

    def UpdateLists(self):
        if self.cluMed  is "metric":
            # print("i nupdatelist")
            self.SepList = GetSeparationList(self.ObjList, self.pars, self.k, self.R)
        else:
            self.SepList = GetDistanceList(self.ObjList, self.pars, self.R)


    def Run(self):
        while len(self.ObjList) > 0:
            self.GetMin()
        # print(self.JetList)
        return self.JetList


def StudyInputs(filelist, treetype, outname):

    JStudyMap = {
        "ExAK4"   : JetStudy("ExAK4"),
        "WTAAK4"  : JetStudy("WTAAK4"),
        "AproAK4" : JetStudy("AproAK4"),
        "FjAK4"   : JetStudy("FjAK4")
    }

    for par in uproot.iterate(filelist, 'Par_%s/tree' % treetype, ['pt', 'eta', 'phi', 'mass']):
        for i, x in enumerate(zip( par['pt'], par['eta'], par['phi'], par['mass'])):
            print(i)
            objs = np.array(list(zip(x[0], x[1], x[2],x[3])))
            jetmap = {}

            vv = np.zeros(len(objs), dtype=[('pT', '<f8'), ('eta', '<f8'), ('phi', '<f8'), ('mass', '<f8')])
            transObjs = objs.transpose()
            vv['pT'] = transObjs[0]
            vv['eta'] = transObjs[1]
            vv['phi'] = transObjs[2]
            vv['mass'] = transObjs[3]

            exak4s = ClusAK4(objs, cluMed="metric", merMed="EMerge")
            wtaak4s = ClusAK4(objs, cluMed = "metric", merMed="WTA")
            apcak4s = ClusAK4(objs, cluMed="distance", merMed="WTA")
            jetmap["ExAK4"] = exak4s.Run()
            jetmap["WTAAK4"] = wtaak4s.Run()
            jetmap["AproAK4"] = apcak4s.Run()

            s = cluster(vv, R= 0.4, p=-1)
            fjs = s.inclusive_jets()
            fjj = []
            for j in fjs:
                fjj.append(np.array([j.pt, j.eta, j.phi, j.mass]))
            jetmap["FjAK4"] = fjj

            for k, v in JStudyMap.items():
                v.Run(jetmap)

            if IsTest:
                break


    f = root_open("%s.root" % outname, "recreate")
    f.mkdir("d")
    for k, v in JStudyMap.items():
        f.cd()
        f.mkdir(k)
        f.cd(k)
        for name, h  in v.GetHist().items():
            h.write(name)
    f.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--inputFiles', help='an integer for the accumulator')
    parser.add_argument('--outputFile')
    args = parser.parse_args()

    ifiles = []
    with open(args.inputFiles) as files:
        for j in files.readlines():
            ifiles.append(j.strip())
    StudyInputs(ifiles, "PF", args.outputFile)


