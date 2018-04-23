#!/usr/bin/python

from watchmakers.load import *
from watchmakers.io_operations import *
from watchmakers.analysis import *
from watchmakers.sensitivity import *
from watchmakers.data import *

######################## Start of main function ###########################

if __name__ == "__main__":

    print docstring
    print arguments
#     print defaultValues
    print ""

    if arguments['-D']:
        deleteAllWorkDirectories()

    if arguments['-m']:
        if arguments['--newVers']:
            generateMacrosNew(int(arguments['-N']),int(arguments['-e']))
        else:
            generateMacros(int(arguments['-N']),int(arguments['-e']))

    if arguments['-j']:
        if arguments['--newVers']:
            generateJobsNew(int(arguments['-N']),arguments)
        else:
            generateJobs(int(arguments['-N']),arguments)

    if arguments['-n']:
        extractNtuple(arguments)

    if arguments['--extractNtup']:
        extractNtupleALL(arguments)

    if arguments['-M']:
        if arguments['--newVers']:
            mergeNtupleFilesNew(arguments)
        else:
            mergeNtupleFiles(arguments)

    if arguments['-A']:
        g,h = {},{}
        g,h = runAnalysisProcess(arguments["-f"],g,h)
        writeResultsToFile(arguments["-o"],g,h)

    if  arguments["--efficiency"]:
        extractHistogramWitCorrectRate()

    if arguments['--sensitivity']:
        sensitivityMapNew()

    if arguments['--customJob']:
        customJob(arguments)

    if arguments['--pass1Trigger']:
        performPass1(arguments)

    if arguments['--pass2Trigger']:
        performPass2(arguments)

    if arguments['--fileDict']:
        createFileDictionary(arguments)

    if arguments['--PDFs']:
        extractPDFandCorrectRate()

    if arguments['--sensIBD']:
        if arguments['--newVers']:
            sensitivityMapPass2New()
        else:
            sensitivityMapPass2()

    if arguments['--evalRate']:
        if arguments['--newVers']:
            readEfficiencyHistogram()
        else:
            readEfficiencyHistogram()

    if arguments['--findRate']:
        if arguments['--newVers']:
            findRate()
        else:
            findRate()


######################## Waba Luba Dub Dub!! ###########################
