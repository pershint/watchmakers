#import watchmakers as PR

from watchmakers.load import *
from watchmakers.analysis import *
from io_operations import testEnabledCondition

from numpy import max
t = arguments['--timeScale']

def drange(start, stop, step):
	rii= start
	while rii<stop:
		yield rii
		rii+=step


def accidentalCanvas(Graphs):

    proc    = []
    loca    = []
    type    = []
    color   = []
    lineS   = []
    acc     = []

    proc    += ['234Pa','214Pb','214Bi','210Bi','210Tl']
    loca    += ['PMT','PMT','PMT','PMT','PMT']
    type    += ['si','si','si','si','si']
    color   += [kO+0,  kO+1, kO+2, kO+3, kO-2]

    proc    +=['232Th','228Ac','212Pb','212Bi','208Tl']
    loca    +=['PMT','PMT','PMT','PMT','PMT']
    type    += ['si','si','si','si','si']
    color   += [kM+0,  kM-6, kM+2, kM-3, kM-7]

    proc    +=['214Pb','214Bi','210Bi','210Tl']
    loca    +=['FV','FV','FV','FV']
    type    += ['ei','ei','ei','ei']
    color   += [kB+0,  kB+1, kB+2, kB+3]

    proc    +=['214Pb','214Bi','210Bi','210Tl']
    loca    +=['FV','FV','FV','FV']
    type    += ['si','si','si','si']
    color   += [kG+0,  kG+1, kG+2, kG+3]

    proc    +=['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    loca    +=['FN','FN','FN','FN','FN','FN','FN','FN']
    type    +=['si','si','si','si','si','si','si','si']
    color   += [kO+2,kO+4,kO+2,kO+4,kO+2,kO+4,kO+2,kO+4]


def siteCanvas(Graphs,site,cut,hist):
    proc        = []
    loca        = []
    type        = []
    color       = []
    lineS       = []
    acc         = []
    scale       = []

    #fast neutrons
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]


    lineS       += [1,1,1,1,1,1,1,1]
    scale       += [1.,1.,1.,1.,1.,1.,1.,1.]


    for PC in range(10,41):
        S = 0.
        B = Graphs['all_cut_accidental'].Eval(PC)

        for i,arr in enumerate(proc):
            _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
            if scale[i] > 0.:
                B+=Graphs[_str].Eval(PC)*scale[i]
            elif scale[i] < 0:
                S+=Graphs[_str].Eval(PC)
				#        print cut,S,B,S/B,S/sqrt(B+S)
        hist.Fill(PC,cut,S/sqrt(B+S))

    C1 = TCanvas('C1','%s'%(t),1200,800)

    Graphs['all_cut_accidental'].Draw('AC')
    Graphs['all_cut_accidental'].SetLineColor(2)
    Graphs['all_cut_accidental'].SetLineWidth(3)
    Graphs['all_cut_accidental'].SetTitle('doubles rate (all cuts)')
    Graphs['all_cut_accidental'].GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    Graphs['all_cut_accidental'].GetXaxis().SetTitle('photo-coverage [%]')
    if site == 'boulby':
        Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,70/31.)
    else:
        Graphs['all_cut_accidental'].GetYaxis().SetRangeUser(0,400/31.)

    for i,arr in enumerate(proc):
        _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
        Graphs[_str].Draw('same')
        Graphs[_str].SetLineColor(color[i])
        Graphs[_str].SetLineStyle(lineS[i])
        if loca[i] == 'PMT':
            _descrpt  = '%s %s glass'%(proc[i],loca[i])
            Graphs[_str].SetTitle(_descrpt)
        elif loca[i] == 'FV':
            _descrpt  = '%s %s in water (%s)'%(proc[i],loca[i],type[i])
            Graphs[_str].SetTitle(_descrpt)
        else:
            if proc[i]=='imb':
                Graphs[_str].SetTitle('Prompt-positron at Perry')
            if proc[i]=='boulby':
                Graphs[_str].SetTitle('Prompt-positron at Boulby')
        Graphs[_str].SetMarkerStyle(1)

    t1 = TLatex()
    if site == 'boulby':
        line = TLine(7, 4.,42, 4)
        t1.DrawLatex( 10,65., 'pe > %d.0'%(cut))
    else:
        line = TLine(7, 40,42, 40)
        t1.DrawLatex( 10,363., 'pe > %d.0'%(cut))

    line.SetLineColor(4)
    line.Draw('same')
    C1.SetGridy()
    C1.Update()

	#    C1.SaveAs('gif/pc_%s_%d.gif'%(site,cut))
    C1.SaveAs('gif/pc_%s.gif+'%(site))
    return hist


def sensitivityMap():

    site = arguments["--site"]
    if site != 'boulby':
        site = ''
    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print site,OnOffRatio

    cores = int(arguments["--cores"])
    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print "Value of reduction of radionuclide greater than 1, setting to 0"
            rnRedux = 0.0
    else:
         rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.

    if t == 'day':
        timeAdjustment = 1.0

    if t == 'month':
        timeAdjustment = 1./31.

    if t == 'year':
        timeAdjustment = 1./365.

    maxTime = 400.*timeAdjustment

    proc        = []
    loca        = []
    type        = []
    color       = []
    lineS       = []
    acc         = []
    scale       = []
    #fast neutrons
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]
    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    if additionalString == "":
        additionalString = "_default"

    if site == 'boulby':
        location = 'Boulby '
    else:
        location = 'Fairport '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();

    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)

	   # palette=TPaletteAxis(hist.FindObject("palette"))
	   # #palette=histo.FindObject("palette")
	   # palette.SetLabelFont(42)
	   # palette.SetLabelSize(0.030)
	   # palette.SetX1NDC(0.905)
	   # palette.SetX2NDC(0.918)
	   # axis = palette.GetAxis()
	   # axis.SetDecimals(kTRUE)
	   # axis.SetMaxDigits(3)
	   # axis.SetNoExponent(kFALSE)
	   # gPad.Update()


    fileIN = 'processed_watchman_March%s.root' %(additionalString)
    fileIN = 'ntuple_root_files%s/processed_watchman.root' %(additionalString)

    print 'Will read in', fileIN

    _graph = Graph()
    _graph.SetName('_graph')
    _graph.SetTitle('off-time [%s] required till a 3#sigma discovery'%(t))

    _graphS = Graph()
    _graphS.SetName('_graphS')
    _graphS.SetTitle('IBDs from %d-core at full power'%(cores))

    _graphB = Graph()
    _graphB.SetName('_graphB')
    _graphB.SetTitle('total background rate after optimization')

    _graphBFN = Graph()
    _graphBFN.SetLineColor(kO+6)
    _graphBFN.SetName('_graphBFN')
    _graphBFN.SetTitle('total fast neutron component')

    _graphBRN = Graph()
    _graphBRN.SetLineColor(kG+3)
    _graphBRN.SetName('_graphBFN')
    _graphBRN.SetTitle('total radionuclide component')


    _graphBAC = Graph()
    _graphBAC.SetLineColor(2)
    _graphBAC.SetName('_graphBFN')
    _graphBAC.SetTitle('total accidental component')


    _graph.SetLineWidth(2)
    _graphS.SetLineWidth(2)
    _graphS.SetLineColor(4)
    _graphB.SetLineWidth(2)
    _graphBAC.SetLineWidth(2)
    _graphBFN.SetLineWidth(2)
    _graphBRN.SetLineWidth(2)




    minY    =   0

    gCnt = 0
    for PC in range(10,41):
		#        print 'processing photocoverage ',PC,'% :',
        S_tmp           = 0.
        B_tmp           = 0.0
        BFN_tmp           = 0.0
        BRN_tmp           = 0.0
        BAC_tmp           = 0.0
        SoverS_B        = 0.
        SoverS_B_sys20  = 0.
        T3SIGMA_MIN     = 16000
        cut_tmp     = 0.
        c1.cd(1)
        hist.Draw('colz')
        gPad.SetTicks()

        c1.cd(2)
        _graph.Draw('AL')
        _graph.GetXaxis().SetTitle('photocoverage [%]')
        _graph.GetYaxis().SetTitle('time [%s]'%(t))
        gPad.SetTicks()
        gPad.SetGridy()


        c1.cd(3)
        _graphS.Draw('AL')
        y = _graphS.GetYaxis().GetXmax()
        _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
        _graphS.GetXaxis().SetTitle('photocoverage [%]')
        _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.cd(4)
        _graphB.Draw('AL')
        _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
        _graphB.GetXaxis().SetTitle('photocoverage [%]')
        _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        _graphBAC.Draw('same')
        _graphBFN.Draw('same')
        _graphBRN.Draw('same')
        c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.Update()
        for cut in range(10,30):
            Graphs =obtainAbsoluteEfficiency(fileIN,timeScale=t,cut=cut)

            S       = 0.
            BFN     = 0.
            BRN     = 0.
            BAC     = 0.
            _str    = 'all_cut_accidental'
            BAC = Graphs[_str].Eval(PC)
            for i,arr in enumerate(proc):
                _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
                if scale[i] > 0.:
                    if loca[i] == 'FN%s'%(site):
                        BFN += Graphs[_str].Eval(PC)*scale[i]
                    elif loca[i] == 'RN%s'%(site):
                        BRN += Graphs[_str].Eval(PC)*scale[i]*(1.0-rnRedux)
                elif scale[i] < 0:
                    S =Graphs[_str].Eval(PC)
                    B = BAC+BRN+BFN
                    SSBB = S/sqrt(S+S+B)
                    #Include 29% other reactor
                    if site == 'boulby' and cores ==1 :
                        T3SIGMA = 9*(S+2*(S*1.29+B)/OnOffRatio)/S/S #1.29% other reactor
                    elif site == 'boulby' and cores ==2 :
                        S*=2                            #2-cores
                        T3SIGMA = 9*(S*1.145+B)/S/S #1.29% other reactor
                    else:
                        T3SIGMA = 9*(S+2*(S*0.05+B)/OnOffRatio)/S/S #5% other reactor

                    if  T3SIGMA < maxTime:
                        print PC,cut,S,S/1.3,BAC,BRN,BFN,B,SSBB,T3SIGMA
                        hist.Fill(PC,cut,T3SIGMA)
                    if  T3SIGMA >= maxTime:
                        hist.Fill(PC,cut,maxTime)

                    if T3SIGMA <   T3SIGMA_MIN:
                            S_tmp       = S
                            B_tmp       = B
                            BFN_tmp       = BFN
                            BRN_tmp       = BRN
                            BAC_tmp       = BAC
                            T3SIGMA_MIN    =T3SIGMA
        _graph.SetPoint(gCnt,PC,T3SIGMA_MIN)
        _graphS.SetPoint(gCnt,PC,S_tmp)
        _graphB.SetPoint(gCnt,PC,B_tmp)
        _graphBFN.SetPoint(gCnt,PC,BFN_tmp)
        _graphBRN.SetPoint(gCnt,PC,BRN_tmp)
        _graphBAC.SetPoint(gCnt,PC,BAC_tmp)

        gCnt+=1

    c1.cd(1)
    hist.Draw('colz')
    gPad.SetTicks()

    c1.cd(2)
    _graph.Draw('AL')
    _graph.GetXaxis().SetTitle('photocoverage [%]')
    _graph.GetYaxis().SetTitle('time [%s]'%(t))
    gPad.SetTicks()
    gPad.SetGridy()


    c1.cd(3)
    _graphS.Draw('AL')
    y = _graphS.GetYaxis().GetXmax()
    _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
    _graphS.GetXaxis().SetTitle('photocoverage [%]')
    _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    c1.cd(4)
    _graphB.Draw('AL')
    _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
    _graphB.GetXaxis().SetTitle('photocoverage [%]')
    _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    _graphBAC.Draw('same')
    _graphBFN.Draw('same')
    _graphBRN.Draw('same')
    c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    hist.SaveAs('results/%sSensitivityMap%s.C'%(site,additionalString))
    _graph.SaveAs('results/%sSensitivityMapGraph%s.C'%(site,additionalString))
    _graphS.SaveAs('results/%sSensitivityMapGraphS%s.C'%(site,additionalString))
    _graphB.SaveAs('results/%sSensitivityMapGraphB%s.C'%(site,additionalString))

    c1.SaveAs('results/%sSensitivityMapCanvas%s.C'%(site,additionalString))
    c1.SaveAs('results/%sSensitivityMapCanvas%s.gif'%(site,additionalString))


def sensitivityMapNew():

    site = arguments["--site"]
    if site != 'boulby':
        site = ''
    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print site,OnOffRatio

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print "Value of reduction of radionuclide greater than 1, setting to 0"
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.

    if t == 'day':
        timeAdjustment = 1.0

    if t == 'month':
        timeAdjustment = 1./31.

    if t == 'year':
        timeAdjustment = 1./365.

    maxTime = 400.*timeAdjustment

    proc        = []
    loca        = []
    type        = []
    color       = []
    lineS       = []
    acc         = []
    scale       = []
    #fast neutrons
    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]
    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    if additionalString == "":
        additionalString = "_default"

    if site == 'boulby':
        location = 'Boulby '
    else:
        location = 'Fairport '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();

    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)

    fileIN = 'processed_watchman_March%s.root' %(additionalString)
    fileIN = 'ntuple_root_files%s/processed_watchman.root' %(additionalString)

    print 'Will read in', fileIN
    print "%s %s %10s %15s %15s %15s %13s %13s %13s %13s" %('PC','cut','S','eff','BAC','BRN','BFN','B_tot','SSBB','T3SIGMA')
    _graph = Graph()
    _graph.SetName('_graph')
    _graph.SetTitle('off-time [%s] required till a 3#sigma discovery'%(t))

    _graphS = Graph()
    _graphS.SetName('_graphS')
    _graphS.SetTitle('IBDs from %d-core at full power'%(cores))

    _graphB = Graph()
    _graphB.SetName('_graphB')
    _graphB.SetTitle('total background rate after optimization')

    _graphBFN = Graph()
    _graphBFN.SetLineColor(kO+6)
    _graphBFN.SetName('_graphBFN')
    _graphBFN.SetTitle('total fast neutron component')

    _graphBRN = Graph()
    _graphBRN.SetLineColor(kG+3)
    _graphBRN.SetName('_graphBFN')
    _graphBRN.SetTitle('total radionuclide component')


    _graphBAC = Graph()
    _graphBAC.SetLineColor(2)
    _graphBAC.SetName('_graphBFN')
    _graphBAC.SetTitle('total accidental component')

    _graph.SetLineWidth(2)
    _graphS.SetLineWidth(2)
    _graphS.SetLineColor(4)
    _graphB.SetLineWidth(2)
    _graphBAC.SetLineWidth(2)
    _graphBFN.SetLineWidth(2)
    _graphBRN.SetLineWidth(2)

    minY    =   0

    gCnt = 0
    for PC in range(10,41):
		#        print 'processing photocoverage ',PC,'% :',
        S_tmp           = 0.
        B_tmp           = 0.0
        BFN_tmp           = 0.0
        BRN_tmp           = 0.0
        BAC_tmp           = 0.0
        SoverS_B        = 0.
        SoverS_B_sys20  = 0.
        T3SIGMA_MIN     = 16000
        cut_tmp     = 0.

        c1.cd(1)
        hist.Draw('colz')
        gPad.SetTicks()

        c1.cd(2)
        _graph.Draw('AL')
        _graph.GetXaxis().SetTitle('photocoverage [%]')
        _graph.GetYaxis().SetTitle('time [%s]'%(t))
        gPad.SetTicks()
        gPad.SetGridy()

        c1.cd(3)
        _graphS.Draw('AL')
        y = _graphS.GetYaxis().GetXmax()
        _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
        _graphS.GetXaxis().SetTitle('photocoverage [%]')
        _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.cd(4)
        _graphB.Draw('AL')
        _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
        _graphB.GetXaxis().SetTitle('photocoverage [%]')
        _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
        _graphBAC.Draw('same')
        _graphBFN.Draw('same')
        _graphBRN.Draw('same')
        c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
        gPad.SetGridy()
        gPad.SetTicks()

        c1.Update()

        for cut in range(4,40):
            Graphs =obtainAbsoluteEfficiency(fileIN,timeScale=t,cut=cut)

            S       = 0.
            BFN     = 0.
            BRN     = 0.
            BAC     = 0.
            _str    = 'all_cut_accidental'
            BAC = Graphs[_str].Eval(PC)
            for i,arr in enumerate(proc):
                _str = 'scaled_%s_%s_%s_1_abs_%s' %(type[i],proc[i],loca[i],acc[i])
                if scale[i] > 0.:
                    if loca[i] == 'FN%s'%(site):
                        BFN += Graphs[_str].Eval(PC)*scale[i]
                        # print 'BFN',BFN,Graphs[_str].Eval(PC),scale[i]
                    elif loca[i] == 'RN%s'%(site):
                        BRN += Graphs[_str].Eval(PC)*scale[i]*(1.0-rnRedux)
                        # print 'BRN',BRN,Graphs[_str].Eval(PC),scale[i]

                elif scale[i] < 0:
                    S =Graphs[_str].Eval(PC)
                    B = BAC+BRN+BFN
                    SSBB = S/sqrt(S+S+B)
                    #Include 29% other reactor
                    if site == 'boulby' and cores ==1 :
                        T3SIGMA = 9*(S+2*(S*1.29+B)/OnOffRatio)/S/S #1.29% other reactor
                    elif site == 'boulby' and cores ==2 :
                        S*=2                            #2-cores
                        T3SIGMA = 9*(S*1.145+B)/S/S #1.29% other reactor
                    else:
                        T3SIGMA = 9*(S+2*(S*0.05+B)/OnOffRatio)/S/S #5% other reactor

                    if  T3SIGMA < maxTime:
                        print PC,cut,S,S/1.3,BAC,BRN,BFN,B,SSBB,T3SIGMA
                        hist.Fill(PC,cut,T3SIGMA)
                    if  T3SIGMA >= maxTime:
                        hist.Fill(PC,cut,maxTime)

                    if T3SIGMA <   T3SIGMA_MIN:
                            S_tmp       = S
                            B_tmp       = B
                            BFN_tmp       = BFN
                            BRN_tmp       = BRN
                            BAC_tmp       = BAC
                            T3SIGMA_MIN    =T3SIGMA
        _graph.SetPoint(gCnt,PC,T3SIGMA_MIN)
        _graphS.SetPoint(gCnt,PC,S_tmp)
        _graphB.SetPoint(gCnt,PC,B_tmp)
        _graphBFN.SetPoint(gCnt,PC,BFN_tmp)
        _graphBRN.SetPoint(gCnt,PC,BRN_tmp)
        _graphBAC.SetPoint(gCnt,PC,BAC_tmp)

        gCnt+=1

    c1.cd(1)
    hist.Draw('colz')
    gPad.SetTicks()

    c1.cd(2)
    _graph.Draw('AL')
    _graph.GetXaxis().SetTitle('photocoverage [%]')
    _graph.GetYaxis().SetTitle('time [%s]'%(t))
    gPad.SetTicks()
    gPad.SetGridy()


    c1.cd(3)
    _graphS.Draw('AL')
    y = _graphS.GetYaxis().GetXmax()
    _graphS.GetYaxis().SetRangeUser(0.,_graphS.GetYaxis().GetXmax()*2.)
    _graphS.GetXaxis().SetTitle('photocoverage [%]')
    _graphS.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    c1.cd(3).BuildLegend(0.42, 0.84, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    c1.cd(4)
    _graphB.Draw('AL')
    _graphB.GetYaxis().SetRangeUser(0., _graphB.GetYaxis().GetXmax()*2.)
    _graphB.GetXaxis().SetTitle('photocoverage [%]')
    _graphB.GetYaxis().SetTitle('event rate [%s^{-1}]'%(t))
    _graphBAC.Draw('same')
    _graphBFN.Draw('same')
    _graphBRN.Draw('same')
    c1.cd(4).BuildLegend(0.42, 0.72, 0.85, 0.88)
    gPad.SetGridy()
    gPad.SetTicks()

    hist.SaveAs('results/%sSensitivityMap%s.C'%(site,additionalString))
    _graph.SaveAs('results/%sSensitivityMapGraph%s.C'%(site,additionalString))
    _graphS.SaveAs('results/%sSensitivityMapGraphS%s.C'%(site,additionalString))
    _graphB.SaveAs('results/%sSensitivityMapGraphB%s.C'%(site,additionalString))

    c1.SaveAs('results/%sSensitivityMapCanvas%s.C'%(site,additionalString))
    c1.SaveAs('results/%sSensitivityMapCanvas%s.gif'%(site,additionalString))

    c2 = TCanvas('c2','c2',1618,1000)
    # c2.SetRightMargin(0.53)
    c2.Divide(1,1)
    c2.cd(1)
    _graph.Draw('AL')
    _graph.GetXaxis().SetTitle('photocoverage [%]')
    _graph.GetYaxis().SetTitle('time [%s]'%(t))
    gPad.SetTicks()
    gPad.SetGridy()
    c2.SaveAs('results/%sSensitivityTim%s.C'%(site,additionalString))
    c2.SaveAs('results/%sSensitivityTime%s.gif'%(site,additionalString))



def runSensitivity():
    hBoulby = TH2D('hBoulby','hBoulby',50,0.5,50.5,50,0.5,50.5)
    print 'Boulby:'
    for cut in range(4,40):
        Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',timeScale=t,cut=cut)
        hBoulby =  siteCanvas(Graphs,'boulby',cut,hBoulby)
    hBoulby.SaveAs('hBoulby.C')


    hFairport = TH2D('hFairport','hFairport',50,0.5,50.5,50,0.5,50.5)
    print 'Fairport:'
    for cut in range(4,40):
        Graphs = obtainAbsoluteEfficiency('processed_data_watchman.root',timeScale=t,cut=cut)
        hFairport =  siteCanvas(Graphs,'',cut,hFairport)
    hFairport.SaveAs('hFairport.C')

    return 0


#hBoulby = TH2D('hBoulby','hBoulby',50,0.5,50.5,50,0.5,50.5)
#hBoulby = sensitivityMap('boulby',hBoulby)
