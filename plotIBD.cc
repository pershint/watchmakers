

{


  TMultiGraph *mg = new TMultiGraph();


  for (int pmt_num=1;pmt_num<4;pmt_num++) {
    double darkrate[4]={1.0,5.0,10.0,15.0};
    double eff[4];
    
    for (int darkcounter=1;darkcounter<16;darkcounter++) {
      double PMTCoverage=10.0*pmt_num;
      if ((darkcounter==1)||(darkcounter==5)||(darkcounter==10)||(darkcounter==15)) {
	double PMTDarkRate=1000*darkcounter;
	
	TString string1="200ns_tts_newLikelihood";
	TString darkRateString="1kHz";
	TString pmtCoverageString="10pct";
	if ((PMTCoverage==10.0)&&(PMTDarkRate==5000))  {darkRateString="5kHz";};
	if ((PMTCoverage==10.0)&&(PMTDarkRate==10000)) {darkRateString="10kHz";};
	if ((PMTCoverage==10.0)&&(PMTDarkRate==15000)) {darkRateString="15kHz";};
	if ((PMTCoverage==20.0)&&(PMTDarkRate==1000)) {pmtCoverageString="20pct";};
	if ((PMTCoverage==20.0)&&(PMTDarkRate==5000)) {pmtCoverageString="20pct";darkRateString="5kHz";}
	if ((PMTCoverage==20.0)&&(PMTDarkRate==10000)) {pmtCoverageString="20pct";darkRateString="10kHz";}
	if ((PMTCoverage==20.0)&&(PMTDarkRate==15000)) {pmtCoverageString="20pct";darkRateString="15kHz";}
	if ((PMTCoverage==30.0)&&(PMTDarkRate==1000)) {pmtCoverageString="30pct";};
	if ((PMTCoverage==30.0)&&(PMTDarkRate==5000)) {pmtCoverageString="30pct";darkRateString="5kHz";}
	if ((PMTCoverage==30.0)&&(PMTDarkRate==10000)) {pmtCoverageString="30pct";darkRateString="10kHz";}
	if ((PMTCoverage==30.0)&&(PMTDarkRate==15000)) {pmtCoverageString="30pct";darkRateString="15kHz";}
	
	
	
	TFile *file_xxpct_darkxxkHz = TFile::Open("../ibd_"+string1+"/dark"+darkRateString+"/ibd_"+pmtCoverageString+"/ntuple_files/ntuple_watchman_ibd_all.root");
	TTree *mytree_xxpct_darkxxkHz = (TTree*) file_xxpct_darkxxkHz->Get("data");
	double pe=0.0;
	double goodness=0.0;
	TVector3 *posTruth;
	TVector3 *posReco;
	int sub_ev_cnt=0;
	int sub_ev=0;
	mytree_xxpct_darkxxkHz->SetBranchAddress("pe",&pe);
	mytree_xxpct_darkxxkHz->SetBranchAddress("pos_goodness",&goodness);
	mytree_xxpct_darkxxkHz->SetBranchAddress("posTruth",&posTruth);
	mytree_xxpct_darkxxkHz->SetBranchAddress("posReco",&posReco);
	mytree_xxpct_darkxxkHz->SetBranchAddress("sub_ev_cnt",&sub_ev_cnt);
	mytree_xxpct_darkxxkHz->SetBranchAddress("sub_ev",&sub_ev);
	
	
	
	TString hpe_name="hpe_"+pmtCoverageString+"_dark"+darkRateString;
	TH1D *hpe_xxpct_darkxxkHz=new TH1D(hpe_name,"Events in Fiducial Vol. ("+pmtCoverageString+", "+darkRateString+")",220,-20,200.0);
	mytree_xxpct_darkxxkHz->Draw("pe >> "+hpe_name,"((sub_ev_cnt==1&&sub_ev==1)||(sub_ev_cnt==2.&&sub_ev==1)||(sub_ev_cnt==0&&sub_ev==0))&&sqrt(posTruth->X()*posTruth->X()+posTruth->Y()*posTruth->Y())<5420.&&posTruth->Z()>-5420.&&posTruth->Z()<5420");
	
	hpe_xxpct_darkxxkHz->SetXTitle("Photo electrons");
	hpe_xxpct_darkxxkHz->SetLineColor(4);
	hpe_xxpct_darkxxkHz->SetLineWidth(2);
	hpe_xxpct_darkxxkHz->GetXaxis()->SetTitleSize(.06);
	gStyle->SetTitleFontSize(.08);
	gStyle->SetLabelSize(.06, "XY");
	TCanvas * c2 = new TCanvas("c", "c", 600, 800);
	c2->Divide(1,3);
	c2->cd(1);
	hpe_xxpct_darkxxkHz->Draw();
	//c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_"+pmtCoverageString+"_dark"+darkRateString+".C");
        //c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_"+pmtCoverageString+"_dark"+darkRateString+".pdf");
        //c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_"+pmtCoverageString+"_dark"+darkRateString+".root");
	
	double numIBD_fiducialTruth=hpe_xxpct_darkxxkHz->Integral();
	cout << "number of events in fiducial volume is " << numIBD_fiducialTruth << "\n";
	
	TH1D *hpe_xxpct_darkxxkHz_detectedPrompt=new TH1D(hpe_name+"_detectedPrompt","Detecte events in Fiducial Vol. ("+pmtCoverageString+", "+darkRateString+")",200,-20,200.0);
	c1->cd();
	mytree_xxpct_darkxxkHz->Draw("pe >> "+hpe_name+"_detectedPrompt","sub_ev_cnt==2.&&sub_ev==1.&&(sqrt(posReco->X()*posReco->X()+posReco->Y()*posReco->Y())<5420.&&posReco->Z()>-5420.&&posReco->Z()<5420)");
	hpe_xxpct_darkxxkHz_detectedPrompt->SetXTitle("Photo electrons");
	hpe_xxpct_darkxxkHz_detectedPrompt->SetTitle("Good reconstructed fiducial prompt events ("+pmtCoverageString+", "+darkRateString+")");
	hpe_xxpct_darkxxkHz_detectedPrompt->SetLineColor(4);
	hpe_xxpct_darkxxkHz_detectedPrompt->SetLineWidth(2);
	hpe_xxpct_darkxxkHz_detectedPrompt->GetXaxis()->SetTitleSize(.06);
	c2->cd(2);
	hpe_xxpct_darkxxkHz_detectedPrompt->Draw();
	//c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detectedPrompt"+pmtCoverageString+"_dark"+darkRateString+".C");
        //c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detectedPrompt"+pmtCoverageString+"_dark"+darkRateString+".pdf");
        //c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detectedPrompt"+pmtCoverageString+"_dark"+darkRateString+".root");
	
	double numIBD_fiducialReco=hpe_xxpct_darkxxkHz_detectedPrompt->Integral();
	cout << "number of prompt events in fiducial volume is " << numIBD_fiducialReco << "\n";
	
	TH1D *hpe_xxpct_darkxxkHz_detectedDelayed=new TH1D(hpe_name+"_detectedDelayed","Detecte events in Fiducial Vol. ("+pmtCoverageString+", "+darkRateString+")",200,-20,200.0);
	c1->cd();
	mytree_xxpct_darkxxkHz->Draw("pe >> "+hpe_name+"_detectedDelayed","sub_ev_cnt==2.&&sub_ev==2.&&(sqrt(posReco->X()*posReco->X()+posReco->Y()*posReco->Y())<5420.&&posReco->Z()>-5420.&&posReco->Z()<5420)");
	hpe_xxpct_darkxxkHz_detectedDelayed->SetXTitle("Photo electrons");
	hpe_xxpct_darkxxkHz_detectedDelayed->SetTitle("Good reconstructed fiducial delayed events ("+pmtCoverageString+", "+darkRateString+")");
	hpe_xxpct_darkxxkHz_detectedDelayed->SetLineColor(4);
	hpe_xxpct_darkxxkHz_detectedDelayed->SetLineWidth(2);
	hpe_xxpct_darkxxkHz_detectedDelayed->GetXaxis()->SetTitleSize(.06);
	c2->cd(3);
	hpe_xxpct_darkxxkHz_detectedDelayed->Draw();
	c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detected"+pmtCoverageString+"_dark"+darkRateString+".C");
        c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detected"+pmtCoverageString+"_dark"+darkRateString+".pdf");
        c2->Print("plots/ibd_200ns_tts_newLikelihood/hpe_detected"+pmtCoverageString+"_dark"+darkRateString+".root");
	
	double numIBD_fiducialReco=hpe_xxpct_darkxxkHz_detectedDelayed->Integral();
	cout << "number of delayed events in fiducial volume is " << numIBD_fiducialReco << "\n";
	
	//now go through each of the events
	int numEvents=mytree_xxpct_darkxxkHz->GetEntries();
	int numFidEvents=0;
	int numTriggeredAndReconEvents=0;
	bool goodPrompt=false;
	TVector3 promptReco(0.0,0.0,0.0);
	TVector3 delayedReco(0.0,0.0,0.0);
	double promptDelayedDist=0.0;
	TH1D *hpromptDelayedDist_xxpct_darkxxkHz = new TH1D("hpromptDelayedDist_"+pmtCoverageString+"_dark"+darkRateString,"",400,0,8000.0);
	hpromptDelayedDist_xxpct_darkxxkHz->SetTitle("Prompt Delayed Distance ("+pmtCoverageString+", "+darkRateString+")");
	hpromptDelayedDist_xxpct_darkxxkHz->SetLineColor(4);
	hpromptDelayedDist_xxpct_darkxxkHz->SetLineWidth(2);
	hpromptDelayedDist_xxpct_darkxxkHz->SetXTitle("mm");
	
	for (int evnum=0;evnum<numEvents;evnum++) {
	  mytree_xxpct_darkxxkHz->GetEntry(evnum);
	  if (sub_ev_cnt==2.&&(sqrt(posReco->X()*posReco->X()+posReco->Y()*posReco->Y())<5420.&&posReco->Z()>-5420.&&posReco->Z()<5420)) {
	    //either a good prompt or delayed event
	    if (sub_ev==1) {
	      promptReco.SetXYZ(posReco->X(),posReco->Y(),posReco->Z());
	      goodPrompt=true;
	    }
	    if (sub_ev==2) {
	      delayedReco.SetXYZ(posReco->X(),posReco->Y(),posReco->Z());
	      if (goodPrompt) {
		promptDelayedDist=sqrt((posReco->X()-promptReco.X())**2+(posReco->Y()-promptReco.Y())**2+(posReco->Z()-promptReco.Z())**2);
		hpromptDelayedDist_xxpct_darkxxkHz->Fill(promptDelayedDist);
		goodPrompt=false;
	      }
	    }
	  } 
	  
	}
	
	c1->cd();
	hpromptDelayedDist_xxpct_darkxxkHz->Draw();
	double efficiency=hpromptDelayedDist_xxpct_darkxxkHz->Integral(1,100)/numIBD_fiducialTruth;
	cout << endl;
	cout << "total number of antineutrino interactions in fiducial " << numIBD_fiducialTruth << endl;
	cout << "total prompt/delayed pairs is " << hpromptDelayedDist_xxpct_darkxxkHz->Integral() << endl;
	cout << "total number of good reconstructed events is " << hpromptDelayedDist_xxpct_darkxxkHz->Integral(1,100) << endl;
	cout << "efficiency is " << efficiency << endl;
	c1->Print("plots/ibd_200ns_tts_newLikelihood/hpromptDelayedDist_"+pmtCoverageString+"_dark"+darkRateString+".C");
	c1->Print("plots/ibd_200ns_tts_newLikelihood/hpromptDelayedDist_"+pmtCoverageString+"_dark"+darkRateString+".pdf");
	c1->Print("plots/ibd_200ns_tts_newLikelihood/hpromptDelayedDist_"+pmtCoverageString+"_dark"+darkRateString+".root");
	
	if (darkcounter==1) eff[0]=efficiency;
	if (darkcounter==5) eff[1]=efficiency;
	if (darkcounter==10) eff[2]=efficiency;
	if (darkcounter==15) eff[3]=efficiency;
	
	
      }
    }
    
    TGraph *g1=new TGraph(4,darkrate,eff);
    g1->SetName(pmtCoverageString);
    g1->SetLineColor(pmt_num+1);
    g1->SetMarkerColor(pmt_num+1);
    g1->SetMarkerStyle(21);
    g1->GetXaxis()->SetTitle("PMT Dark Rate (kHz)");
    if (pmt_num==1) g1->SetTitle("10% PMT coverage");
    if (pmt_num==2) g1->SetTitle("20% PMT coverage");
    if (pmt_num==3) g1->SetTitle("30% PMT coverage");
    g1->Draw("ACP");

    mg->Add(g1,"lp");
    
  }

  mg->Draw("A");
  mg->GetXaxis()->SetTitle("PMT Dark Rate (kHz)");
  mg->GetYaxis()->SetTitle("Efficiency");
  mg->GetYaxis()->SetRangeUser(0.,1.);
  TLegend* leg_frac = new TLegend(0.12, 0.45, 0.35, 0.85);
  leg_frac->AddEntry(g1, "", "f");

  c1->BuildLegend();
  c1->Print("plots/ibd_200ns_tts_newLikelihood/ibdEfficiency.C");
  c1->Print("plots/ibd_200ns_tts_newLikelihood/ibdEfficiency.pdf");
  c1->Print("plots/ibd_200ns_tts_newLikelihood/ibdEfficiency.root");


}

