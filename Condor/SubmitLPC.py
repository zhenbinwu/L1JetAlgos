#!/usr/bin/env python
# coding: utf-8

import os
import re
import time
import subprocess
import glob
import tarfile
import shutil
import getpass

DelExe    = '../L1JetAlgo.py'
# DelExe    = '../FatJetProducer_cfg.py'
OutDir = '/store/user/%s/Phase2L1/JetAlgo' %  getpass.getuser()
# tempdir = '/uscmst1b_scratch/lpc1/lpctrig/benwu/CondorTemp/' % getpass.getuser()
# tempdir = '/uscmst1b_scratch/lpc1/lpctrig/benwu/CondorTemp/' % getpass.getuser()
tempdir = '/uscmst1b_scratch/lpc1/3DayLifetime/benwu/TestCondor/'
# ProjectName = 'Ak4forMurli_v0'
ProjectName = '10XJets_v4'
argument = "--inputFiles=%s.$(Process).list --outputFile=%s_$(Process)"
# argument = "inputFiles=%s.$(Process).list outputFile=%s_$(Process) maxEvents=10 \n Queue %d \n"

Process = {
    'QCD_PU0'     : ['QCD_PU0.list',     3],
    # 'QCD_PU140'   : ['QCD_PU140.list',   3],
    # 'TTbar_PU0'   : ['TTbar_PU0.list',   3],
    # 'TTbar_PU140' : ['TTbar_PU140.list', 3],
    # 'TTbar_PU140' : ['../../TTbar_PU140.list', 10],
    # 'TTbar_PU200' : ['../../TTbar_PU140.list', 10],
}

def Condor_Sub(condor_file):
    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(condor_file))
    print "To submit condor with " + condor_file
    os.system("condor_submit " + condor_file)
    os.chdir(curdir)


def SplitPro(key, file, fraction):
    splitedfiles = []
    filelistdir = tempdir + '/' + "FileList"
    try:
        os.makedirs(filelistdir)
    except OSError:
        pass


    filename = os.path.abspath(file)
    if fraction == 1:
        splitedfiles.append(os.path.abspath(filename))
        shutil.copy2(os.path.abspath(filename), "%s/%s" % (filelistdir, os.path.basename(filename)))
        return splitedfiles

    f = open(filename, 'r')
    lines = f.readlines()
    if len(lines) <= fraction:
        lineperfile = 1
        fraction = len(lines)
    else:
        lineperfile = len(lines) / fraction
        if len(lines) % fraction > 0:
            lineperfile += 1


    for i in range(0, fraction):
        wlines = []
        if i == fraction - 1 :
            wlines = lines[lineperfile*i :]
        else:
            wlines = lines[lineperfile*i : lineperfile*(i+1)]
        if len(wlines) > 0:
            outf = open("%s/%s.%d.list" % (filelistdir, key, i), 'w')
            outf.writelines(wlines)
            splitedfiles.append(os.path.abspath("%s/%s.%d.list" % (filelistdir, key, i)))
        outf.close()

    return splitedfiles

def my_process():
    ## temp dir for submit
    global tempdir
    global Mergeblock
    global ProjectName
    ProjectName = time.strftime('%b%d') + ProjectName
    tempdir = tempdir + os.getlogin() + "/" + ProjectName +  "/"
    print tempdir
    try:
        os.makedirs(tempdir)
    except OSError:
        pass

    ## Create the output directory
    outdir = OutDir +  "/" + ProjectName + "/"
    try:
        os.makedirs("/eos/uscms/%s" % outdir)
    except OSError:
        pass

    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = tempdir + "/" + "RunExe.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("RunExe.csh", "r"):
            #line = line.replace("DELDIR", os.environ['PWD'])
            line = line.replace("DELSCR", os.environ['SCRAM_ARCH'])
            line = line.replace("DELDIR", os.environ['CMSSW_VERSION'])
            line = line.replace("DELEXE", DelExe.split('/')[-1])
            line = line.replace("OUTDIR", outdir)
            outfile.write(line)

    ### Create Tarball
    NewNpro = {}
    Tarfiles = []
    for key, value in Process.items():
        if value[0] == "":
            value[0] = "../FileList/"+key+".list"
        if not os.path.isfile(value[0]):
            continue
        npro = GetProcess(key, value)
        Tarfiles+=npro
        NewNpro[key] = len(npro)
    Tarfiles+=glob.glob("../*py")

    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.environ['CMSSW_BASE'])
    tarballname ="%s/CMSSW.tgz" % tempdir
    with tarfile.open(tarballname, "w:gz", dereference=False) as tar:
        [tar.add(f) for f in TarCMSSW()]
        tar.close()
    os.chdir(curdir)
    process = subprocess.Popen( "xrdcp %s root://cmseos.fnal.gov/%s " % (tarballname, outdir ),
                               shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    tarballnames = []

    Tarfiles.append(os.path.abspath(DelExe))
    # Tarfiles += GetNeededFileList(key)
    tarballname ="%s/%s.tar.gz" % (tempdir, ProjectName)
    with tarfile.open(tarballname, "w:gz", dereference=True) as tar:
        [tar.add(f, arcname=f.split('/')[-1]) for f in Tarfiles]
        tar.close()
    tarballnames.append(tarballname)
    tarballnames.append("/uscms/home/benwu//pyjet.tgz")

    ### Update condor files
    for key, value in Process.items():
        arg = "\nArguments = " + argument % (key, key) + "\n Queue %d \n" % NewNpro[key]
        # if NewNpro[key] > 1:
            # arg = "\nArguments = %s.$(Process).list %s_$(Process).root \nQueue %d \n" % (key, key, NewNpro[key])
        # else:
            # arg = "\nArguments = inputFiles=%s.list outputFile=%s maxEvents=10 \n Queue\n" % (key, key)

        ## Prepare the condor file
        condorfile = tempdir + "/" + "condor_" + ProjectName +"_" + key
        with open(condorfile, "wt") as outfile:
            for line in open("condor_template", "r"):
                line = line.replace("EXECUTABLE", os.path.abspath(RunHTFile))
                line = line.replace("TARFILES", ", ".join(tarballnames))
                line = line.replace("TEMPDIR", tempdir)
                line = line.replace("PROJECTNAME", ProjectName)
                line = line.replace("ARGUMENTS", arg)
                line = line.replace("X509_USER_PROXY", os.environ['X509_USER_PROXY'])
                outfile.write(line)

        Condor_Sub(condorfile)

def GetProcess(key, value):
    if len(value) == 1:
        return SplitPro(key, value[0], 1)
    else :
        return SplitPro(key, value[0], value[1])

def GetNeededFileList(key):
    relist = []
    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.environ['CMSSW_BASE'])
    g = glob.glob("")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*root")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*csv")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*cfg")
    relist += [os.path.abspath(h) for h in g]
    g = glob.glob("../*model")
    relist += [os.path.abspath(h) for h in g]
    process = subprocess.Popen( "ldd %s " % os.path.abspath(DelExe) , shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for l in process.stdout:
        if os.getenv('USER') in l:
            relist.append(l.strip().split(' ')[2])
    os.chdir(curdir)
    return relist

def TarCMSSW():
    from itertools import chain
    relist = []
    # result = [chain.from_iterable(glob(os.path.join(x[0], '*')) for x in os.walk('.'))]
    libs = [y for x in os.walk('lib') for y in glob.glob(os.path.join(x[0], '*'))]
    srcs = [ y for x in os.walk('src') for y in glob.glob(os.path.join(x[0], '*'))]
    # print srcs
    wantsrc = []
    for j in srcs:
        if os.path.splitext(j)[1] == ".py":
            wantsrc.append(j)
        if "/data/" in j:
            wantsrc.append(j)

    return libs + wantsrc


if __name__ == "__main__":
    my_process()

