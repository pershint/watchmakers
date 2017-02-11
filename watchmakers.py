#!/usr/bin/python

from watchmakers.load import *
from watchmakers.io_operations import *
from watchmakers.analysis import *
from watchmakers.sensitivity import *
from watchmakers.read import *

######################## Start of main function ###########################

if __name__ == "__main__":
    
    if arguments['-D']:
        deleteAllWorkDirectories()

    if arguments['-m']:
        generateMacros(int(arguments['-N']),int(arguments['-e']))

    if arguments['-j']:
        generateJobs(int(arguments['-N']),arguments)

    if arguments['-n']:
        extractNtuple(arguments)

    if arguments['--extractNtup']:
        extractNtupleALL(arguments)

    if arguments['-M']:
        mergeNtupleFiles(arguments)

    if arguments['-a']:
        g,h = {},{}
        g,h = runAnalysisProcess(arguments["-f"],g,h)
        writeResultsToFile(arguments["-o"],g,h)

    if arguments['-R']:
#        site = 'boulby'
#        hBoulby = TH2D('hBoulby','hBoulby',40,0.5,40.5,40,0.5,40.5)
        sensitivityMap()
#        hBoulby.SaveAs('test.C')
#        runSensitivity()


    if arguments['--customJob']:
        customJob(arguments)


######################## Waba Luba Dub Dub!! ###########################
