from load import *
from io_operations import testEnabledCondition

def fillHistograms(inFile,a1,t1,h,cover,ii,locj,covPCT):
    # Obtain logarithmic binnings
    nbins, xbins, ybins = logx_logy_array()
    
    additionalString,additionalCommands = testEnabledCondition(arguments)
    
    fiducialVolume = float(arguments["--fv"])
    pmtDist         = float(arguments["--psup"])
#    print "Fiducial volume is ", fiducialVolume
    #Read-in file
    try:
        s =  "ntuple_root_files%s/%s_%s_%s_%s.root"%(additionalString,inFile,ii,cover,locj)
        print "Reading in ",s 
#        data = root2array(s)
        t           = root2rec(s)

        #Apply some analysis
        r           = npa(t.reco_r<fiducialVolume,dtype=bool)
        z           = npa(absolute(t.reco_z)<fiducialVolume,dtype=bool)
        
        #isFV        = logical_and(r,z,dtype=bool)
        isFV        = npa(t.FV==1,dtype=bool)
        notFV       = npa(t.FV!=1,dtype=bool)
        isFV_t      = npa(t.FV_truth==1,dtype=bool)
        notFV_t     = npa(t.FV_truth!=1,dtype=bool)

        iCandidate  = npa(t.candidate==1,dtype=bool)
        totalEvtGen = npa(t.all_ev==t.all_ev_tot,dtype=bool)
        
        tot         = float(sum(totalEvtGen))
        totD        = float(len(t.FV))

        si          = logical_and(isFV,notFV_t,dtype=bool)
        so          = logical_and(notFV,isFV_t,dtype=bool)
        ei          = logical_and(isFV,isFV_t,dtype=bool)
        eo          = logical_and(notFV,notFV_t,dtype=bool)

        ssi         = sum(si)
        sso         = sum(so)
        sei         = sum(ei)
        seo         = sum(eo)
        print 'MC events generated',tot,', number of MC recorded',len(t.FV)
        print ssi,sso,sei,seo,', absolute: ', ssi/tot,sso/tot,sei/tot,seo/tot, \
        'relative',ssi/totD,sso/totD,sei/totD,seo/totD

        #Define set of histograms to fill out
        for subE in range(5):
            if subE ==0:
                subEv0   = npa(t.candidate==1,dtype=bool)
            else:
                subEv0   = npa(t.detected_ev==subE,dtype=bool)
            suma = sum(subEv0)
            if suma:
                print 'Sub event ',subE,' : ',suma
                totRel = float(suma)
                mask = logical_and(subEv0,si,dtype=bool)
                if sum(mask):
                    type = "si"
                    print subE,type,totRel,'(',sum(mask),')'

                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)

                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl       ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDist))
                    s_dl       = "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]
                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDist))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)




                mask = logical_and(subEv0,so)
                if sum(mask):

                    type        = "so"
                    print subE,type,totRel,'(',sum(mask),')'

                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)
                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)


                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl        ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]

                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDist))
                    s_dl        ="%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]

                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDist))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)


                mask = logical_and(subEv0,ei)
                if sum(mask):

                    type        = "ei"
                    print subE,type,totRel,'(',sum(mask),')'


                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)
                    anaVar      = t.pe[mask,...]
                    s_dl        = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)


                    anaVar      = t.nhit[mask,...]
                    s_dl= "%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl= "%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDist))
                    s_dl= "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)

                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]

                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDist))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)

                mask = logical_and(subEv0,eo)
                if sum(mask):

                    type        = "eo"
                    print subE,type,totRel,'(',sum(mask),')'
                    string = "%s_%s_%s_%d_abs"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/tot)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/tot)
                    string = "%s_%s_%s_%d_rel"%(type,ii,locj,subE)
                    if string in h:
                        point = h[string].GetN()
                        h[string].SetPoint(point,covPCT,sum(mask)/totRel)
                    else:
                        h[string] = Graph()
                        h[string].SetName(string)
                        h[string].SetPoint(0,covPCT,sum(mask)/totRel)

                    anaVar      = t.pe[mask,...]
                    s_dl         = "%s_pe_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    anaVar      = t.nhit[mask,...]
                    s_dl        ="%s_nhit_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]     = Hist(5000,0,500,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('photoelectrons')
                    h[s_dl].SetYTitle('counts [a.u.]')
                    h[s_dl].fill_array(anaVar)

                    i_d         = t.inner_dist[mask,...]
                    i_t         = t.inner_time[mask,...]
                    anaVar      = column_stack((i_d,i_t/1e9))
                    s_dl        ="%s_dDdT_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(xbins,ybins,name=s_dl,title=s_dl);
                    h[s_dl].SetXTitle('inter-event distance [m]')
                    h[s_dl].SetYTitle('inter-event time [s]')
                    h[s_dl].fill_array(anaVar)
                    r_r         = t.reco_r[mask,...]
                    r_z         = t.reco_z[mask,...]
                    anaVar      = column_stack((power(r_r/pmtDist,2),r_z/pmtDist))
                    s_dl       = "%s_dRdZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)
                    t_r         = t.true_r[mask,...]
                    t_z         = t.true_z[mask,...]
                    anaVar      = column_stack((power(t_r/pmtDist,2),t_z/pmtDist))
                    s_dl       = "%s_tRtZ_%s_%s_%s_%d"%(type,cover,ii,locj,subE)
                    h[s_dl]= Hist2D(100,0,1.6,200,-1.6,1.6,name=s_dl,title=s_dl)
                    h[s_dl].SetXTitle('(r / R_{pmt})^{2} [unitless]')
                    h[s_dl].SetYTitle('(z / Z_{pmt}) [unitless]')
                    h[s_dl].fill_array(anaVar)


    except:
        print "Could not read file ",s

    print ""
#        f.Close()
    return h

def logx_logy_array(nbins = 500,xmin = 1e-2,xmax = 30.,ymin = 1e-9,ymax = 1e3):
    #x-axis
    logxmin = log10(xmin)
    logxmax = log10(xmax)
    xbinwidth= (logxmax-logxmin)/float(nbins)
    xbins	= zeros(nbins,dtype=float)
    xbins[0] = xmin
    #y-axis
    logymin = log10(ymin)
    logymax = log10(ymax)
    ybinwidth= (logymax-logymin)/float(nbins)
    ybins	= zeros(nbins,dtype=float)
    ybins[0] = ymin
    for i in range(nbins):
        xbins[i] = xmin + pow(10,logxmin+i*xbinwidth)
        ybins[i] = ymin + pow(10,logymin+i*ybinwidth)
    return nbins,xbins,ybins




def obtainAbsoluteEfficiency(f,timeScale='day',cut = 10.0):

    covPCT      = {'9.86037':1432., '14.887':2162.,'19.4453':2824.,'24.994':3558.,'28.8925':4196.,'34.3254':4985.,'39.1385':5684.}
    pct         = npa([9.86037,14.887,19.4453,24.994,28.8925,34.3254,39.1385])
    pctVal      = ['10pct','15pct','20pct','25pct','30pct','35pct','40pct']

    f           = TFile(f,'read')
    EG          = {}
    


    inta,proc,loca,acc,arr,Activity,br,site,timeS,boulbyNeutronFV,mass = loadAnalysisParameters(timeScale)


    x,y         = Double(0.),Double(0.)
    boolFirst   = True
    for _inta in inta:
        for ii in range(len(proc)):
            _str1               = "%s_%s_%s_1_abs" %(_inta,proc[ii],loca[ii])
            _scal_str1          = "scaled_%s_%s_%s%s_1_abs_%s"%(_inta,proc[ii],\
            loca[ii],site[ii],acc[ii])
            _strEff             = "eff_%s_%s_%s%s_1_abs" %(_inta,proc[ii],\
            loca[ii],site[ii])

            EG[_strEff]         = Graph()
            cnter = 0
            EG[_strEff].SetPoint(cnter,0.,0.)
            cnter+=1
            for cc,val in enumerate(pctVal):
                s              = "%s_pe_%s_%s_%s_1"%(_inta,val,proc[ii],loca[ii])
                eff = histIntegral(s,f,cut)
                if eff>0.00:
                    EG[_strEff].SetPoint(cnter,pct[cc],eff)
                    cnter+=1
#                    if loca[ii]=="RN":
#                        print _inta,val,proc[ii],loca[ii],cnter,pct[cc],eff

            try:
                EG[_str1] = f.Get(_str1)
                EG[_scal_str1] = EG[_str1].Clone()

                for i in range(EG[_str1].GetN()):
                    EG[_str1].GetPoint(i,x,y)
                    nY      = y*Activity[ii]*timeS*br[ii]*EG[_strEff].Eval(x)
#                    if _inta == 'ei' and loca[ii] == 'S':
#                        print 'A',_str1,_strEff,y,Activity[ii],timeS,br[ii],EG[_strEff].Eval(x),nY
                    if loca[ii] == 'PMT':
                        nY      *= mass*covPCT["%s"%(x)]
                        EG[_scal_str1].SetPoint(i,x,nY)
                    else:
                        EG[_scal_str1].SetPoint(i,x,nY)
            except:
                a = 0

    _scal_acc,_scal_acc_notFV,_scal_acc1,_scal_acc_notFV1 ="scaled_accidental",\
    "scaled_accidental_notFV","cut_accidental","all_cut_accidental"
    EG[_scal_acc],EG[_scal_acc_notFV],EG[_scal_acc1],EG[_scal_acc_notFV1] = \
    Graph(),Graph(),Graph(),Graph()

    for i,p  in enumerate(pct):
        EG[_scal_acc].SetPoint(i ,p ,0.0)
        EG[_scal_acc_notFV].SetPoint(i,p,0.)
        EG[_scal_acc1].SetPoint(i,p,0.)
        EG[_scal_acc_notFV1].SetPoint(i,p,0.)


    for _inta in inta:
        for ii in range(len(proc)):
            _scal_str1          = "scaled_%s_%s_%s%s_1_abs_%s"%(_inta,\
            proc[ii],loca[ii],site[ii],acc[ii])
            for iii,value in enumerate(pct):
                try:
                    if acc[ii] == 'acc' and (_inta == 'si' or _inta == 'ei'):
                        x   = float(value)
                        oY  = EG[_scal_acc].Eval(x)
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_acc].SetPoint(iii,x,oY+nY)
                    if acc[ii] == 'di' and (_inta == 'si' or _inta == 'ei'):
                        x   = float(value)
                        aY  = EG['scaled_ei_neutron_Nboulby_1_abs_corr'].Eval(x)
                        oY  = aY/(boulbyNeutronFV*timeS)
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_str1].SetPoint(iii,x,oY*nY)
                    if acc[ii] == 'acc' and _inta != 'si' and _inta != 'ei':
                        x   = float(value)
                        oY  = EG[_scal_acc_notFV].Eval(x)
                        nY  = EG[_scal_str1].Eval(float(value))
                        EG[_scal_acc_notFV].SetPoint(iii,x,oY+nY)
                except:
                    a = 1


    distanceEff = 0.004
    for iii,value in enumerate(pct):
        x   = float(value)
        oY  = EG[_scal_acc].Eval(x)
        nY  = 0.0001/timeS*oY*oY
        EG[_scal_acc1].SetPoint(iii,x,nY)
        EG[_scal_acc_notFV1].SetPoint(iii,x,nY*distanceEff)
    return EG



def pickColor(H,_loc,r_c,o_c,b_c,c_c ):
    if _loc=='PMT':
        H.SetLineColor(kOrange+o_c)
        H.SetFillColor(kOrange+o_c)
        H.SetMarkerColor(kOrange+o_c)
        o_c+=1
    if _loc=='RN':
        H.SetLineColor(7)
        H.SetMarkerColor(7)
    if _loc=='FV':
        H.SetLineColor(kRed+r_c)
        H.SetFillColor(kRed+r_c)
        H.SetMarkerColor(kRed+r_c)
        r_c+=1
    if _loc=='S':
        H.SetLineColor(kBlue+b_c)
        H.SetFillColor(kBlue+b_c)
        H.SetMarkerColor(kBlue+b_c)
        b_c+=1
    if _loc=='N':
        H.SetLineColor(kCyan+c_c)
        H.SetMarkerColor(kCyan+c_c)
        c_c+=1
    return r_c,o_c,b_c,c_c, H


def integralCoincidence(R,lowerBound,upperBound):
    low = -exp(-lowerBound) * (1+lowerBound )
    up  = -exp(-upperBound) * (1+upperBound)
    return up - low


def histIntegral(s,f,cut):
    H = {}
    H[s] = f.Get(s)
    if H[s]!= None:
        a = H[s].Integral(0,int(cut*10))
        N = H[s].Integral(0,5000)
    else:
        a = 1.0
        N = 1.0

    return (1.0 - a/N )

def runAnalysisProcess(f,g,h):
    d,iso,loc,coverage,coveragePCT = loadSimulationParameters()
    aDT1 =tDT1 = a1 = t1 = 0

    for j in range(len(iso)):
        for ii in d["%s"%(iso[int(j)])]:
            for idx,cover in enumerate(coverage):
                print  coverage[idx],coveragePCT[cover],ii,loc[j]
                h  = fillHistograms(f,a1,t1,h,cover,ii,loc[j],\
                float(coveragePCT[cover]))
    return g,h
