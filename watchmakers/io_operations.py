from load import *



def testCreateDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

def deleteDirectory(directory):
    if os.path.exists(directory):
        rmtree(directory)
#


def macroGenerator(percentage,isotope,location,runs,events):
    
    covPCT = {'10pct':0.1,'15pct':0.15,'20pct':0.2,\
    '25pct':0.25,'30pct':0.30,'35pct':0.35,'40pct':0.40}

    #Part of the macro that is the same for all jobs
    dir = os.getcwd()
    
    header = '''
/glg4debug/glg4param omit_muon_processes  0.0
/glg4debug/glg4param omit_hadronic_processes  0.0

/rat/db/set DETECTOR experiment "Watchman"
/rat/db/set DETECTOR detector_factory "Watchman"
/rat/db/set WATCHMAN_PARAMS photocathode_coverage %4.2f


/run/initialize

# BEGIN EVENT LOOP
/rat/proc lesssimpledaq
/rat/proc fitbonsai
/rat/proc fitcentroid
/rat/proc fitpath
/rat/proc count
/rat/procset update 1000

# Use IO.default_output_filename
/rat/proclast outroot
/rat/procset file "%s/root_files/%s/%s/watchman_%s_%s_%s_%d.root"
#END EVENT LOOP

''' %(covPCT[percentage],dir,isotope,percentage,isotope,percentage,location,runs)
    
    
    #Part of macro that varies with the various conditions
    if location == 'PMT':
        line1 = '''
/generator/add decaychain %s:regexfill
/generator/pos/set inner_pmts[0-9]+

/run/beamOn %d''' %(isotope,events*5)
    elif location == 'FV':
        line1 = '''
/generator/add decaychain %s:fill:poisson
/generator/pos/set  0 0 0
/generator/rate/set 6.43

/run/beamOn %d''' %(isotope,events*5)
    elif location == 'FN':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 1434.0
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(events)
    elif location == 'FNimb':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 1434.0
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(events)
    elif location == 'FNboulby':
        line1 = '''
/generator/add combo fastneutron:regexfill
/generator/pos/set rock_[0-9]+
/generator/vtx/set 0 0 0
/generator/fastneutron/depth 2805.
/generator/fastneutron/enthresh 10.0
/generator/fastneutron/sidewalls 1.0

/run/beamOn %d'''%(events)
    elif    location=='I':
        line1 = '''
/generator/add combo ibd:fill
/generator/vtx/set  1 0 0
/generator/pos/set 0 0 0

/run/beamOn %d'''%(events)
    elif location == 'S':
        line1 ='''
/generator/add combo spectrum:fill
/generator/vtx/set e+ %s
/generator/pos/set 0 0 0

/run/beamOn %d'''%(isotope,events)
    elif location =='N':
        line1 = '''
/generator/add combo gun2:fill
/generator/vtx/set %s  0 0 0 0 0.001 0.20
/generator/pos/set 0 0 0

/run/beamOn %d'''%(isotope,events)
    elif location == 'RN':
        AZ = isotope
        A =  int(int(AZ)/1000)
        Z = int(AZ) - A*1000
        line1 = '''
/generator/add combo isotope:fill
/generator/pos/set 0 0 0
/generator/vtx/set GenericIon 0 0 0
/generator/isotope/A %s.0
/generator/isotope/Z %s.0
/generator/isotope/E 0.0

/run/beamOn %d''' %(A,Z,events)
    else:
        line1 = 'A'
        print location
    return header+line1



def jobString(percentage,j,runs,models,arguments):
#    directory = "/p/lscratche/adg/Watchboy/simplifiedData/rp_sim/wm"
    directory   = os.getcwd()
    softDir     = "/usr/gapps/adg/geant4/rat_pac_and_dependency"
    ratDir      = os.environ['RATROOT']
    rootDir     = os.environ['ROOTSYS']
    g4Dir       =  os.environ['G4INSTALL']
    watchmakersDir = os.environ['WATCHENV']
    
    software    = "%s/bin/rat" %(ratDir)
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    
    ele =  d["%s"%(iso[j])]
    location = loc[j]
    
    goodness     = float(arguments['-g'])
    case = int(arguments['-j'])
#    print defaultValues

    additionalString      = ""
    additionalCommands    = ""

    if float(arguments['-r'])          != defaultValues[6]:
        additionalString += "_rate_%f" %(float(arguments['-r']))
        additionalCommands += " -r %f " %(float(arguments['-r']))
    
    if float(arguments['-d'])          != defaultValues[7]:
        additionalString += "_deltaR_%f" %(float(arguments['-d']))
        additionalCommands += " -d %f" %(float(arguments['-d']))
    
    if float(arguments['-t'])          != defaultValues[8]:
        additionalString += "_deltaT_%f" %(float(arguments['-t']))
        additionalCommands +=  " -t %f" %(float(arguments['-t']))

    if float(arguments['-T'])            != (defaultValues[9]):
        additionalString += "_nhitMin_%d" %(int(arguments['-T']))
        additionalCommands += " -T %d" %(int(arguments['-T']))

    if float(arguments['-g'])          != defaultValues[10]:
        additionalString += "_posGood_%f" %(float(arguments['-g']))
        additionalCommands += " -g %f" %(float(arguments['-g']))

    if float(arguments['-G'])          != defaultValues[11]:
        additionalString += "_dirGood_%f" %(float(arguments['-G']))
        additionalCommands += " -G %f" %(float(arguments['-G']))

    if float(arguments['--fv'])        !=  defaultValues[12]:
        additionalString += "_FVboundary_%f" %(float(arguments['--fv']))
        additionalCommands +=  "--fv %f" %(float(arguments['--fv']))

    if float(arguments['--psup'])      != defaultValues[13]:
        additionalString += "_PMTboundary_%f" %(float(arguments['--psup']))
        additionalCommands += "--psup %f" %(float(arguments['--psup']))

    if float(arguments['--tankDis'])   != defaultValues[14]:
        additionalString += "_Tankboundary_%f" %(float(arguments['--tankDist']))
        additionalCommands +=" --tankDist %f" %(float(arguments['--tankDist']))


    if additionalString != "":
        print additionalString
        print additionalCommands
        print type(additionalCommands)

    line1 = """#!/bin/sh
#MSUB -N WM_%s_%s_%d    #name of job
#MSUB -A adg         # sets bank account
#MSUB -l nodes=1:ppn=1,walltime=23:59:59,partition=borax  # uses 1 node
#MSUB -q pbatch         #pool
#MSUB -o %s/log/wmpc_%s_%s_%d.log
#MSUB -e %s/log/wmpc_%s_%s_%d.err
#MSUB -d %s  # directory to run from
#MSUB -V
#MSUB                     # no more psub commands

source %s/bin/thisroot.sh
source %s/../../../bin/geant4.sh
source %s/geant4make.sh
source %s/env.sh
source %s/env_wm.sh
export G4NEUTRONHP_USE_ONLY_PHOTONEVAPORATION=1\n
""" %(percentage,location,runs,\
directory,percentage,location,runs,\
directory,percentage,location,runs,\
directory,\
rootDir,g4Dir,g4Dir,ratDir,watchmakersDir)

    for mods in models:
        if location == "FN":
            line1 += "export PHYSLIST=%s\n" %(mods)
        if case == 1 or case == 2:
            line1 += "%s -l log/%s/%s/rat.%s_%s_%s_%d.log %s/macro_%s/run%s_%s_%d.mac\n" %(software,\
                                                      mods,percentage,percentage,mods,location,runs,\
                                                                                 directory,percentage,mods,location,runs)
        if case == 1 or case == 3:
            fileN = "root_files/%s/%s/watchman_%s_%s_%s_%d.root" %(mods,percentage,mods,percentage,location,runs)
            if additionalString != "":
                fileNO = "ntuple_root_files/%s/%s/watchman_%s_%s_%s%s_%d.root" %(mods,percentage,mods,percentage,location,additionalString,runs)
                line1 += "watch -n %s -f %s --ntupleout %s\n" %(additionalCommands,fileN,fileNO)

            else:
                line1 += "watch -n -f %s\n" %(fileN)

    return line1



def generateMacros(N,e):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
#    N = int(arguments['-N'])
    print N,e
    ##Clean or create macro directories
    for idx,cover in enumerate(coverage):
        dir = "macro_%s" %(cover)
        testCreateDirectory(dir)
#    for idx,cover in enumerate(coverage):
#        dir = "%s" %(cover)
#        testCreateDirectory(dir)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                for val in range(N):
                    line = macroGenerator(cover,ii,loc[j],val,e )
                    dir = "macro_%s" %(cover)

                    outfile = open("%s/run%s_%s_%d.mac" %(dir,ii,\
                    loc[j],val),"wb")
                    outfile.writelines(line)
                    outfile.close
    return 0

def removeMacrosAndDirectories():
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    for idx,cover in enumerate(coverage):
        dir = "macro_%s" %(cover)
        deleteDirectory(dir)


def generateJobs(N,arguments):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    
    '''Find wheter the jobs folder exist: if no create, if yes clean and recreate'''
    directory = 'jobs'
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

    for ii in loc:
        for idx,cover in enumerate(coverage):
            directory = "jobs/%s/%s" %(ii,cover)
            if not os.path.exists(directory):
                os.makedirs(directory)
            else:
                rmtree(directory)
                os.makedirs(directory)

    '''Find wheter the jobs folder exist: if no create, if yes clean and recreate'''
    directory = 'log'
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        rmtree(directory)
        os.makedirs(directory)

#    directory = 'root_files'
#    if not os.path.exists(directory):
#        os.makedirs(directory)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "root_files/%s/%s" %(ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)


#
#    directory = 'ntuple_root_files'
#    if not os.path.exists(directory):
#        os.makedirs(directory)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "ntuple_root_files/%s/%s" %(ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                directory = "log/%s/%s" %(ii,cover)
                if not os.path.exists(directory):
                    os.makedirs(directory)


    '''Make sure that the softlink are correct for Bonsai input'''

    ratDir      = os.environ['RATROOT']

    src = ratDir+'/fit_param.dat'
    dst = os.getcwd()+'/fit_param.dat'
    if not os.path.exists(dst):
        os.symlink(src,dst)

    src = ratDir+'/like.bin'
    dst = os.getcwd()+'/like.bin'
    if not os.path.exists(dst):
        os.symlink(src,dst)

    job_list = '''#!/bin/sh
'''

    for j in range(len(iso)):
        for idx,cover in enumerate(coverage):
            models  = d["%s" %(iso[j])]
            for index in range(N):
                line = jobString(cover,j,index,models,arguments)
                stringFile = "jobs/%s/%s/jobs%s_%s_%s_%d.sh" %(loc[j],cover,cover,\
                                                            "%s"%(iso[int(j)]),loc[j],index)
                if index == 0:
                    job_list+= '(msub ' + stringFile +') || ./'+ stringFile + '\n'
                
                outfile = open(stringFile,"wb")
                outfile.writelines(line)
                if index < N-1:
                    stringFile1 = "(msub jobs/%s/%s/jobs%s_%s_%s_%d.sh || ./jobs/%s/%s/jobs%s_%s_%s_%d.sh)" %(models,cover,cover,\
                                                                                                 "%s"%(iso[int(j)]),loc[j],index+1,models,cover,cover,\
                                                                                                 "%s"%(iso[int(j)]),loc[j],index+1)
                    outfile.writelines(stringFile1)
                outfile.close
                os.chmod(stringFile,S_IRWXU)


    outfile = open('sub_jobs',"wb")
    outfile.writelines(job_list)
    outfile.close
    os.chmod('sub_jobs',S_IRWXG)
    os.chmod('sub_jobs',S_IRWXU)
    return 0

def deleteAllWorkDirectories():
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()

    dir = "log"
    if os.path.exists(dir):
        rmtree(dir)

    dir = "jobs"
    if os.path.exists(dir):
        rmtree(dir)

    for idx,cover in enumerate(coverage):
        dir = "macro_%s" %(cover)
        if os.path.exists(dir):
            rmtree(dir)

    if os.path.exists('fit_param.dat'):
        os.remove('fit_param.dat')

    if os.path.exists('like.bin'):
        os.remove('like.bin')

    if os.path.exists('sub_jobs'):
        os.remove('sub_jobs')

def writeResultsToFile(s,g,h):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    f_root = TFile(s,"recreate")
    for gE in g:
        g["%s"%(gE)].Write()
    for hE in h:
        h["%s"%(hE)].Write()
    f_root.Close()


def mergeFiles():
    # Read external requirements
    #arguments = docopt.docopt(docstring)
    # Load internal requirements
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    trees = {}
    pathFinal = "root_files/merged_ntuple_watchman"
    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                t_name  = "data_%s_%s_%s"%(ii,cover,loc[j])
                trees[t_name] = TChain("data")
                
                s = "ntuple_root_files/%s/%s/watchman_%s_%s_%s_*.root" %(ii,cover,ii,cover,loc[j])
                sw = "%s_%s_%s_%s.root"%(pathFinal,ii,cover,loc[j])
            
                print "Writing ", sw,"from",s
                trees[t_name].Add(s)
                print "Number of entries ",trees[t_name].GetEntries()
                trees[t_name].Merge(sw)
                del trees[t_name]
    del trees
    return 0


def extractNtuple(arguments):
    N            = int(arguments["-N"])
    rate         = float(arguments["-r"])
    timemask     = float(arguments['-t'])*1000.0
    distancemask = float(arguments['-d'])
    goodness     = float(arguments['-g'])
    dirGoodness  = float(arguments['-G'])
    minNHIT      = float(arguments['-T'])
    file         = arguments["-f"]
    fidV         = float(arguments["--fv"])
    pmtV         = float(arguments["--psup"])
    tankV        = float(arguments["--tankDis"])
    outF         = arguments["--ntupleout"]
    
    print file
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    goldenFileExtractor(file,minNHIT,goodness,dirGoodness,timemask,\
                        rate,distancemask,fidV,pmtV,tankV,outF)
