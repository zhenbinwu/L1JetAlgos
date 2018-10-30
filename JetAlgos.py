#!/usr/bin/env python
# encoding: utf-8

# File        : JetAlgos.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Oct 30
#
# Description : 

import uproot
import itertools
import numpy as np
import ROOT
import math

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

# CluMetric vectorization
def CluMetricV(x, y, k, R=0.4):
    Xpt = np.power(x[:, 0], k * 2)
    Ypt = np.power(y[:, 0], k * 2)
    deta = np.power(x[:, 1]-y[:, 1], 2)
    dPhi = x[:, 2] - y[:, 2]
    np.subtract(dPhi, 2*math.pi, out = dPhi, where= (dPhi >=math.pi))
    np.add(dPhi, 2*math.pi,  out =dPhi , where=(dPhi < -1*math.pi))
    np.power(dPhi, 2, out=dPhi)
    return np.multiply(np.divide(np.add(deta, dPhi), R**2), np.minimum(Xpt, Ypt))

def CludR(x, y):
    deta = pow(x[1]-y[1], 2)
    dphi = pow(deltaPhi(x[2], y[2]), 2)
    return (deta+dphi)

# CluR vectorization
def CludRV(x, y):
    deta = np.power(x[:, 1]-y[:, 1], 2)
    dPhi = x[:, 2] - y[:, 2]
    np.subtract(dPhi, 2*math.pi, out = dPhi, where= (dPhi >=math.pi))
    np.add(dPhi, 2*math.pi,  out =dPhi , where=(dPhi < -1*math.pi))
    np.power(dPhi, 2, out=dPhi)
    return  np.add(deta, dPhi)

def GetSeparationList(idx, pars, k, R):
    combm = {}
    if len(idx) >= 2:
        combs  = np.asarray(list(itertools.combinations(idx, 2)))
        vecVal = CluMetricV(pars[combs[:, 0]], pars[combs[:, 1]], k, R)
        combm  = dict(zip(tuple(map(tuple, combs)) , tuple(vecVal)))
    for i in idx:
        combm[(i, i)] = list(pars[i])[0]**(2*k)
    return combm


def GetDistanceList(idx, pars, R):
    ptlist ={list(pars[x])[0] : x for x in idx}
    combm = {}
    if len(ptlist.keys()) == 0:
        return combm

    maxidx = ptlist[max(ptlist.keys())]
    npars  = pars[idx]
    maxv   = np.tile(pars[maxidx], (len(npars),1))
    dis    = CludRV(maxv, npars)
    maxk   = tuple([(maxidx, x) for x in idx])
    combm  = dict(zip(maxk, dis))
    combm[(maxidx, maxidx)] = R**2
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

    def GetMin(self):
        m = min(self.SepList, key=self.SepList.get)
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
        else:
            newobj = WTAMerging(self.pars[m[0]], self.pars[m[1]])
        self.pars=np.concatenate((self.pars, np.array([newobj])), axis=0)
        self.ObjList = [x for x in self.ObjList if x not in m]
        newIdx =len(self.pars)-1
        self.ObjList.append(newIdx)
        self.UpdateLists()

    def UpdateLists(self):
        if self.cluMed  is "metric":
            self.SepList = GetSeparationList(self.ObjList, self.pars, self.k, self.R)
        else:
            self.SepList = GetDistanceList(self.ObjList, self.pars, self.R)


    def Run(self):
        while len(self.ObjList) > 0:
            self.GetMin()
        return self.JetList
