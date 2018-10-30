#!/usr/bin/env python
# encoding: utf-8

# File        : L1JetAlgo.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2018 Oct 25
#
# Description :


import uproot
import itertools
import numpy as np
import ROOT
import math
import rootpy

from rootpy.io import root_open
from rootpy.tree import TreeChain
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

## User
from JetStudy import JetStudy
from JetAlgos import ClusAK4

IsTest=False


def StudyInputs(filelist, treetype, outname):

    JStudyMap = {
        "ExAK4"   : JetStudy("ExAK4"),
        "WTAAK4"  : JetStudy("WTAAK4"),
        "AproAK4" : JetStudy("AproAK4"),
    }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Write the output ~~~~~
    f = root_open("%s.root" % outname, "recreate")
    treemap = {}

    class Jets(TreeModel):
        num_vals = IntCol()
        pts = FloatArrayCol(500, length_name='num_vals')
        etas = FloatArrayCol(500, length_name='num_vals')
        phis = FloatArrayCol(500, length_name='num_vals')
        mass = FloatArrayCol(500, length_name='num_vals')

    for k, v in JStudyMap.items():
        treemap [k] = Tree(k, model=Jets)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running code ~~~~~
    for par_ in uproot.iterate(filelist, 'Par_%s/tree' % treetype, ['pt', 'eta', 'phi', 'mass']):
        par = {key.decode("utf-8"): value for (key, value) in par_.items()}
        for i, x in enumerate(zip( par['pt'], par['eta'], par['phi'], par['mass'])):
            print(i)
            objs = np.array(list(zip(x[0], x[1], x[2],x[3])))
            jetmap = {}


            exak4s = ClusAK4(objs, cluMed="metric", merMed="EMerge")
            wtaak4s = ClusAK4(objs, cluMed = "metric", merMed="WTA")
            apcak4s = ClusAK4(objs, cluMed="distance", merMed="WTA")
            jetmap["ExAK4"] = exak4s.Run()
            jetmap["WTAAK4"] = wtaak4s.Run()
            jetmap["AproAK4"] = apcak4s.Run()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Save jets ~~~~~
            for k, v in JStudyMap.items():
                jets = np.asarray(jetmap[k])
                treemap[k].num_vals = len(jets)
                treemap[k].pts = jets[:, 0]
                treemap[k].etas = jets[:, 1]
                treemap[k].phis = jets[:, 2]
                treemap[k].mass = jets[:, 3]
                treemap[k].fill()

            for k, v in JStudyMap.items():
                v.Run(jetmap)

            if i == 1:
                break

    for k, v in JStudyMap.items():
        treemap[k].Write()

    newf = ROOT.TChain("myak4%s/tree" % treetype )
    [newf.Add(x) for x in ifiles]
    newt=newf.CloneTree(-1, "fast")
    newt.Write()

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
    StudyInputs(ifiles, "Puppi", args.outputFile)
