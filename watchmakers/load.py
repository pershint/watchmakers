from ROOT import TRandom3
from ROOT import TChain,TGraph,TGraphErrors,gSystem,gROOT,TH1D,TH2D,TFile,TCanvas
from ROOT import THStack,Double
from ROOT import kRed,kBlue,kGreen,kCyan,kOrange

from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
from ROOT import TCanvas,TLine, TLatex
from ROOT import sqrt
from ROOT import gStyle,gPad,TPaletteAxis
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

#    set_style('ATLAS')

    warnings.simplefilter("ignore")
except:
    print "Could not load in root_numpy or rootpy, they are required to run this module."

defaultValues  = [3,2500,'merged_ntuple_watchman','null', \
                  'processed_watchman.root',10.,2.0,100.0,6.0,0.1,0.1,5.42,6.4,8.0,2805.,'day',\
                  'boulby',1.0]

docstring = """
    Usage: watchmakers.py [options]
    
    Arguments:
    
    Options:
    -m                  generate macro files
    -N=<N>              Number of MC script that were run [Default: %d]
    -e=<runBeamEntry>   Number of entries per macro (U/Th event x5) [Default: %d]
    -j=<jobType>        Create submision scripts and macros
    -f=<ifile>          Input file [Default: %s]
    --ntupleout=<outN>  Name of ntuple out [Default: %s]
    -n                  generate ntuple from single rat-pac root files
    --extractNtup       generate ntuple from all rat-pac root files
    -o=<outputfile>     Efficiency output file [Default: %s]
    -M                  Merge result files
    -a                  Do the analysis on the merged file
    -r=<rate>           rate of accidentals in hz [Default: %f]
    -d=<distance>       Maximal distance between two events (m) [Default: %f]
    -t=<time>           Maximal time between two events (micro) [Default: %f]
    -T=<tubes>          Minimal number of tubes hit [Default: %f]
    -g=<goodness>       Bonsai position goodness parameter [Default: %f]
    -G=<Goodness>       Bonsai direction goodness parameter [Default: %f]
    -P=<proc>           Pick a single physics process to analyis/merge (used for ntup)
    -L=<loc>            Pick a single physics location to analyis/merge (used for ntup)
    -C                  Pick a single coverage to analyse
    -R                  Read analyisis result
    -D                  Delete all current photocoverage directory.
    --fv=<fidV>         Fiducial Volome [Default: %f]
    --psup=<psupV>      Distance to PMT support, assuming right cylinder [Default: %f]
    --tankDis=<tankV>   Distance to tank wall, assuming right cylinder [Default: %f]
    --detectorDepth=<DD> Depth of detector [Default: %f]
    --customJob         Custom job for photocoverage 02-2017
    --timeScale=<_ts>   Integration period (sec,day,month,year) [Default: %s]
    --site=<_site>      Site of the experiment (boulby,fairport) [Default: %s]
    --OnOff=<_OOratio>  Ratio of reactor on to reactor off [Default: %d]

    """ % (defaultValues[0],defaultValues[1],defaultValues[2],defaultValues[3],\
           defaultValues[4],defaultValues[5],defaultValues[6],defaultValues[7],\
           defaultValues[8],defaultValues[9],defaultValues[10],defaultValues[11],\
           defaultValues[12],defaultValues[13],defaultValues[14],defaultValues[15],\
           defaultValues[16],defaultValues[17])

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print 'using docopt as the user control interface'
except ImportError:
    print 'docopt is not a recognized module, it is required to run this module'


gSystem.Load("libRATEvent")
gSystem.AddIncludePath(" -I$RATROOT/include")
gROOT.LoadMacro("$WATCHENV/watchmakers/goldenFileExtractor.C")
from ROOT import goldenFileExtractor

def loadSimulationParameters():
    #Chain and subsequent isotopes
    d = {}
    # Water and PMT contamination
    d['CHAIN_238U_NA'] =['238U','234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['232Th','228Ac','212Pb','212Bi','208Tl']
    d['CHAIN_222Rn_NA'] = ['222Rn','214Pb','214Bi','210Bi','210Tl']
    # Radioisotope that should have beta-Neutron modes, (beta only generated)
#    A = ['16','17','18','17','18','8','9','11']
#    Z = ['6','6','6','7','7','2','3','3']
    #Reduced selection
    A = ['9','11']
    Z = ['3','3']
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

def loadAnalysisParameters(timeScale='day'):


    # Default units are in sec. Conversion factor are below
    timeSec     = 1.0/365./24./3600.
    
    # Number of free proton
    if timeScale == 'sec':
        timeS   = 1.0
    if timeScale == 'day':
        timeS   = 24.0*3600.
    if timeScale == 'month':
        timeS   = 365.0/12.*24.0*3600.
    if timeScale == 'year':
        timeS   = 365.0*24.0*3600.
        
    #Mass in kilograms
    mass = 2.0

    #Evaluate FV to total detector volume ratio
    nKiloTons   = 3.22
    FreeProtons = 0.6065
    TNU         = FreeProtons* nKiloTons *timeSec
    FVkTonRatio = pow(float(arguments['--fv']),3)/pow(float(arguments['--tankDis']),3)
    
    
    #Fast neutrons conversion
    #Rock mass
    volumeR         = (2.*22.5*23.8*1.0+2.*17*23.8*1.0+2.*22.5*17.*1.0)
    density         = 2.39 #from McGrath
    rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    
    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    
    
    boulbyIBDRate   = 1120.8*.4/.6 *TNU #//924.48*TNU Taken from website, average corrected
    fairportIBDRate = 7583.*TNU

    inta        = ['si','so','eo','ei']

    #Add the U-238 chain
    proc        = ['234Pa','214Pb','214Bi','210Bi','210Tl']
    loca        = ['PMT',  'PMT',  'PMT',  'PMT',  'PMT']
    acc         = ['acc',  'acc',  'acc',  'acc',  'acc']
    br          = [1.0,     1.0,    1.0,   1.0 ,   0.002]
    site        = ['',      '',     '',     '',     '']
    arr         = empty(5)
    arr[:]      = 0.993
    Activity    = arr
    #Add the Th-232 chain
    proc        +=['232Th','228Ac','212Pb','212Bi','208Tl']
    loca        +=['PMT'  ,'PMT',   'PMT', 'PMT',  'PMT'  ]
    acc         +=['acc'  ,'acc',   'acc', 'acc',  'acc'  ]
    br          += [1.0,     1.0,    1.0,   1.0 ,   1.0]
    site        += ['',      '',     '',     '',     '']
    arr         = empty(5)
    arr[:]      = 0.124
    Activity    = append(   Activity,arr)
    #Add the Rn-222 chain
    proc        +=['214Pb','214Bi','210Bi','210Tl']
    loca        +=['FV',   'FV',   'FV',   'FV']
    acc         +=['acc',  'acc',  'acc',   'acc']
    br          += [1.0,   1.0,   1.0,     0.002]
    site        += ['',     '',     '',     '']
    arr = empty(4)
    arr[:]      = 6.4
    Activity    = append(   Activity,arr)


    #Add the neutrino signal
    proc        +=['imb','imb','boulby','boulby']
    loca        +=['S','S',     'S',  'S']
    acc         +=['di', 'corr', 'di', 'corr']
    br          += [1.0,  1.0, 1.0 , 1.0]
    site        += [''   ,'' , '', '']
    arr         = npa([fairportIBDRate,fairportIBDRate,boulbyIBDRate,boulbyIBDRate])
    Activity    = append(    Activity,arr)

    # Add the neutron
    proc        +=['neutron','neutron']
    loca        +=['N',     'N']
    acc         +=['corr',  'corr']
    br          += [1.0,   1.0]
    arr         = npa([fairportIBDRate,boulbyIBDRate])
#    print "Neutrino activity ",arr*timeS/nKiloTons
    Activity    = append(    Activity,arr)
    site        += [ '','boulby']

    # add a fast neutron at Fairport
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    loca        +=  ['FN','FN','FN','FN','FN','FN','FN','FN']
    acc         +=  ['corr','corr','corr','corr','corr','corr','corr','corr']
    br          +=  [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    arr = empty(8)
    arr[:]      = fastNeutrons[0]
    Activity    = append(Activity,arr)
    site        += ['',      '','',      '','',      '','',      '']
    # add a fast neutron at Boulby
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    loca        +=  ['FN','FN','FN','FN','FN','FN','FN','FN']
    acc         +=  ['corr','corr','corr','corr','corr','corr','corr','corr']
    br          +=  [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]
    arr = empty(8)
    arr[:]      = fastNeutrons[1]
    Activity    = append(Activity,arr)
    site        += ['boulby','boulby','boulby','boulby','boulby','boulby',\
    'boulby','boulby']

#    # Read in the different radionuclide
#    proc        +=  ['16006','17006','18006','17007','18007','8002','9003',\
#    '11003']
#    loca        +=  ['RN','RN','RN','RN','RN','RN','RN','RN']
#    acc         +=  ['di','di','di','di','di','di','di','di']
#    #normalised to 9Li from SK
#    arr         = npa([0.02,0.001,0.001,0.59*0.002,4e-6,0.23,1.9,0.01])/1.9
#    arr         *= radionuclideRate[0]
#    Activity    = append(Activity,arr)
#    br         +=  [.988,1.0,1.0,0.951,0.143,0.16,0.495,0.927]
#    site        += ['','','','','','','','']
#

    # Read in the different radionuclide
    proc        +=  ['9003', '11003']
    loca        +=  ['RN','RN']
    acc         +=  ['di','di']
    #normalised to 9Li from SK
    arr         = npa([1.9,0.01])/1.9
    arr         *= radionuclideRate[0]
    Activity    = append(Activity,arr)
    br         +=  [0.495,0.927]
    site        += ['','']



#    # Read in the different radionuclide
#    proc        +=  ['16006','17006','18006','17007','18007','8002','9003',\
#    '11003']
#    loca        +=  ['RN','RN','RN','RN','RN','RN','RN','RN']
#    acc         +=  ['di','di','di','di','di','di','di','di']
#    #normalised to 9Li from SK
#    arr         = npa([ 0.02,0.001,0.001,0.59*0.002,4e-6,0.23,1.9,0.01])/1.9
#    arr         *= radionuclideRate[1]
#    Activity    = append(Activity,arr)
#    br         +=  [.988,1.0,1.0,0.951,0.143,0.16,0.495,0.927]
#    site        += ['boulby','boulby','boulby','boulby','boulby','boulby',\
#    'boulby','boulby']
#
#

    proc        +=  ['9003', '11003']
    loca        +=  ['RN','RN']
    acc         +=  ['di','di']
    #normalised to 9Li from SK
    arr         = npa([1.9,0.01])/1.9
    arr         *= radionuclideRate[1]
    Activity    = append(Activity,arr)
    br         +=  [0.495,0.927]
    site        += ['boulby','boulby']

    return inta,proc,loca,acc,arr,Activity,br,site,timeS,boulbyIBDRate*FVkTonRatio,mass


