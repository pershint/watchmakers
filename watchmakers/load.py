from ROOT import TRandom3
from ROOT import TChain,TGraph,TGraphErrors,gSystem,gROOT,TH1D,TH2D,TFile,TCanvas
from ROOT import THStack,Double
from ROOT import kRed,kBlue,kGreen,kCyan,kOrange

import os.path
from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import array as npa
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,sqrt,absolute

from math import pow,exp,log10

try:
    from root_numpy import root2rec
    from rootpy.plotting import Canvas,Hist,Hist2D,Graph
    from rootpy.plotting.style import set_style
    from rootpy.io import root_open
    #from rootpy.interactive import wait

    set_style('ATLAS')

    warnings.simplefilter("ignore")
except:
    print "Could not load in root_numpy or rootpy, they are required to run this module."

defaultValues  = [3,5000,0,'root_files/merged_ntuple_watchman','null', \
                  'processed_watchman.root',10.,2.0,100.0,6.0,0.1,0.1,5.42,6.4,8.0]

docstring = """
    Usage: watchmakers.py [options]
    
    Arguments:
    
    Options:
    -m                  generate macro files
    -N=<N>              Number of MC script that were run [Default: %d]
    -e=<runBeamEntry>   Number of entries per macro [Default: %d]
    -j=<jobType>        Create submision scripts and macros [Default: %d]
    -f=<ifile>          Input file [Default: %s]
    --ntupleout=<outN>  Name of ntuple out [Default: %s]
    -n                  generate ntuple from rat-pac root files
    -o=<outputfile>     Efficiency output file [Default: %s]
    -M                  Merge result files
    -a                  Do the analysis on the merged file
    -r=<rate>           rate of accidentals in hz [Default: %f]
    -d=<distance>       Maximal distance between two events (m) [Default: %f]
    -t=<time>           Maximal time between two events (micro) [Default: %f]
    -T=<tubes>          Minimal number of tubes hit [Default: %f]
    -g=<goodness>       Bonsai position goodness parameter [Default: %f]
    -G=<Goodness>       Bonsai direction goodness parameter [Default: %f]
    -P                  Pick a single physics process to analyis (used for ntup)
    -C                  Pick a single coverage to analyse
    -R                  Read analyisis result
    -D                  Delete all current photocoverage directory.
    --fv=<fidV>         Fiducial Volome [Default: %f]
    --psup=<psupV>      Distance to PMT support, assuming right cylinder [Default: %f]
    --tankDis=<tankV>   Distance to tank wall, assuming right cylinder [Default: %f]
    """ % (defaultValues[0],defaultValues[1],defaultValues[2],defaultValues[3],\
           defaultValues[4],defaultValues[5],defaultValues[6],defaultValues[7],\
           defaultValues[8],defaultValues[9],defaultValues[10],defaultValues[11],\
           defaultValues[12],defaultValues[13],defaultValues[14])

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print 'using docopt as the user control interface'
except ImportError:
    print 'docopt is not a recognized module, it is required to run this module'


gSystem.Load("libRATEvent")
gSystem.AddIncludePath(" -I$RATROOT/include")
gROOT.LoadMacro("watchmakers/goldenFileExtractor.C")
from ROOT import goldenFileExtractor

def loadSimulationParameters():
    #Chain and subsequent isotopes
    d = {}
    # Water and PMT contamination
    d['CHAIN_238U_NA'] =['238U','234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['232Th','228Ac','212Pb','212Bi','208Tl']
    d['CHAIN_222Rn_NA'] = ['222Rn','214Pb','214Bi','210Bi','210Tl']
    # Radioisotope that should have beta-Neutron modes, (beta only generated)
    A = ['16','17','18','17','18','8','9','11']
    Z = ['6','6','6','7','7','2','3','3']
    ZA = A
    for i in range(len(A)):
        ZA[i] = str(int(A[i])*1000 +int(Z[i]))
    d['A_Z'] =  ZA
    #Oscillated spectrum at Boulby and IMB site
    d['ibd'] = ['boulby','imb']
    d['N'] = ['neutron']
    d['IBD'] = ['IBD']
    # Fast neutron contamination
    d['FN'] = ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC','QBBC',\
    'QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    #
    #List of all physics process clumped together
    iso = ['CHAIN_238U_NA','CHAIN_232Th_NA','CHAIN_222Rn_NA','A_Z','ibd',\
    'FN','N','IBD']
    loc = ['PMT','PMT','FV','RN','S','FN','N','I']
    #Photocoverage selected
    coverage = ['10pct','15pct','20pct','25pct','30pct','35pct','40pct']
    coveragePCT = {'10pct':9.86037,'15pct':14.887,'20pct':19.4453,\
    '25pct':24.994,'30pct':28.8925,'35pct':34.3254,'40pct':39.1385}

    return d, iso,loc,coverage,coveragePCT

