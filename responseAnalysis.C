///////////////////////////////////
//
// responseAnalysis
//
// - take reco response maps and calculate stuff (JER/JES, closure, whatever...)
//
//
//
//////////////////////////////////

// gaussian fit fxn
double fitfxn(double *x, double *par){
    double pt = x[0];
    double fitval = par[0]*TMath::Exp(-0.5*TMath::Power(((x[0] - par[1])/par[2]),2));
    return fitval;
}

void responseAnalysis(){

	TFile *f1;
	// responses with jet-energy corrections applied
	TH2D *response_corr_allJets, *response_corr_bJets, *response_corr_cJets, *response_corr_udJets, *response_corr_sJets, *response_corr_gJets, *response_corr_xJets;
	TH1D *JER_corr_allJets, *JER_corr_bJets, *JER_corr_cJets, *JER_corr_udJets, *JER_corr_sJets, *JER_corr_gJets, *JER_corr_xJets;
	TH1D *JES_corr_allJets, *JES_corr_bJets, *JES_corr_cJets, *JES_corr_udJets, *JES_corr_sJets, *JES_corr_gJets, *JES_corr_xJets;
	// responses without jet-energy corrections applied
	TH2D *response_raw_allJets, *response_raw_bJets, *response_raw_cJets, *response_raw_udJets, *response_raw_sJets, *response_raw_gJets, *response_raw_xJets;
	TH1D *JER_raw_allJets, *JER_raw_bJets, *JER_raw_cJets, *JER_raw_udJets, *JER_raw_sJets, *JER_raw_gJets, *JER_raw_xJets;
	TH1D *JES_raw_allJets, *JES_raw_bJets, *JES_raw_cJets, *JES_raw_udJets, *JES_raw_sJets, *JES_raw_gJets, *JES_raw_xJets;

	

	TString inputFile;
	inputFile = "/home/clayton/Analysis/code/skimming/PYTHIAHYDJET_scan/rootFiles/response/PYTHIAHYDJET_response_19Aug22.root";
	
	f1 = TFile::Open(inputFile);

	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C3",response_corr_allJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C3",response_corr_bJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C3",response_corr_cJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C3",response_corr_udJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C3",response_corr_sJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C3",response_corr_gJets);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_xJets_C3",response_corr_xJets);

	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_allJets_C3",response_raw_allJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_bJets_C3",response_raw_bJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_cJets_C3",response_raw_cJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_udJets_C3",response_raw_udJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_sJets_C3",response_raw_sJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_gJets_C3",response_raw_gJets);
	f1->GetObject("h_matchedRawJetPtOverGenJetPt_genJetPt_xJets_C3",response_raw_xJets);

    
	
	TH2D *h2_test = (TH2D*) response_corr_allJets->Clone("h2_test");
	
	TH1D *hx_test = h2_test->ProjectionY();
    	TAxis *xaxis = hx_test->GetXaxis();


    	//const int Nedges = 28;
    	//double pt_axis[Nedges] = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100,110,120,140,160,180,210,250,300};
    	const int Nedges = 18;
    	double pt_axis[Nedges] = {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,210,250,300};
    	double x[Nedges-1];
    	double ept_axis[Nedges-1];

	int firstxbin = 0;
    	int lastxbin = 0;


	double JER_corr_allJets_[Nedges-1], JER_corr_bJets_[Nedges-1], JER_corr_cJets_[Nedges-1], JER_corr_udJets_[Nedges-1], JER_corr_sJets_[Nedges-1], JER_corr_gJets_[Nedges-1], JER_corr_xJets_[Nedges-1];	
	double JES_corr_allJets_[Nedges-1], JES_corr_bJets_[Nedges-1], JES_corr_cJets_[Nedges-1], JES_corr_udJets_[Nedges-1], JES_corr_sJets_[Nedges-1], JES_corr_gJets_[Nedges-1], JES_corr_xJets_[Nedges-1];	

	double eJER_corr_allJets_[Nedges-1], eJER_corr_bJets_[Nedges-1], eJER_corr_cJets_[Nedges-1], eJER_corr_udJets_[Nedges-1], eJER_corr_sJets_[Nedges-1], eJER_corr_gJets_[Nedges-1], eJER_corr_xJets_[Nedges-1];	
	double eJES_corr_allJets_[Nedges-1], eJES_corr_bJets_[Nedges-1], eJES_corr_cJets_[Nedges-1], eJES_corr_udJets_[Nedges-1], eJES_corr_sJets_[Nedges-1], eJES_corr_gJets_[Nedges-1], eJES_corr_xJets_[Nedges-1];	
	TF1 *fitFxn_corr_allJets = new TF1("fitFxn_corr_allJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_bJets = new TF1("fitFxn_corr_bJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_cJets = new TF1("fitFxn_corr_cJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_udJets = new TF1("fitFxn_corr_udJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_sJets = new TF1("fitFxn_corr_sJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_gJets = new TF1("fitFxn_corr_gJets",fitfxn,0,300,3);
	TF1 *fitFxn_corr_xJets = new TF1("fitFxn_corr_xJets",fitfxn,0,300,3);


	fitFxn_corr_allJets->SetParameter(1,1.0);
	fitFxn_corr_bJets->SetParameter(1,1.0);
	fitFxn_corr_cJets->SetParameter(1,1.0);
	fitFxn_corr_udJets->SetParameter(1,1.0);
	fitFxn_corr_sJets->SetParameter(1,1.0);
	fitFxn_corr_gJets->SetParameter(1,1.0);
	fitFxn_corr_xJets->SetParameter(1,1.0);
	
	fitFxn_corr_allJets->SetParameter(2,0.3);
	fitFxn_corr_bJets->SetParameter(2,0.3);
	fitFxn_corr_cJets->SetParameter(2,0.3);
	fitFxn_corr_udJets->SetParameter(2,0.3);
	fitFxn_corr_sJets->SetParameter(2,0.3);
	fitFxn_corr_gJets->SetParameter(2,0.3);
	fitFxn_corr_xJets->SetParameter(2,0.3);
	
	fitFxn_corr_allJets->SetParLimits(2,0,1);
	fitFxn_corr_bJets->SetParLimits(2,0,1);
	fitFxn_corr_cJets->SetParLimits(2,0,1);
	fitFxn_corr_udJets->SetParLimits(2,0,1);
	fitFxn_corr_sJets->SetParLimits(2,0,1);
	fitFxn_corr_gJets->SetParLimits(2,0,1);
	fitFxn_corr_xJets->SetParLimits(2,0,1);


	double JER_tempVal_corr_allJets, JER_tempVal_corr_bJets, JER_tempVal_corr_cJets, JER_tempVal_corr_udJets, JER_tempVal_corr_sJets, JER_tempVal_corr_gJets, JER_tempVal_corr_xJets;
	double eJER_tempVal_corr_allJets, eJER_tempVal_corr_bJets, eJER_tempVal_corr_cJets, eJER_tempVal_corr_udJets, eJER_tempVal_corr_sJets, eJER_tempVal_corr_gJets, eJER_tempVal_corr_xJets;

	double JES_tempVal_corr_allJets, JES_tempVal_corr_bJets, JES_tempVal_corr_cJets, JES_tempVal_corr_udJets, JES_tempVal_corr_sJets, JES_tempVal_corr_gJets, JES_tempVal_corr_xJets;
	double eJES_tempVal_corr_allJets, eJES_tempVal_corr_bJets, eJES_tempVal_corr_cJets, eJES_tempVal_corr_udJets, eJES_tempVal_corr_sJets, eJES_tempVal_corr_gJets, eJES_tempVal_corr_xJets;



	double JER_raw_allJets_[Nedges-1], JER_raw_bJets_[Nedges-1], JER_raw_cJets_[Nedges-1], JER_raw_udJets_[Nedges-1], JER_raw_sJets_[Nedges-1], JER_raw_gJets_[Nedges-1], JER_raw_xJets_[Nedges-1];	
	double JES_raw_allJets_[Nedges-1], JES_raw_bJets_[Nedges-1], JES_raw_cJets_[Nedges-1], JES_raw_udJets_[Nedges-1], JES_raw_sJets_[Nedges-1], JES_raw_gJets_[Nedges-1], JES_raw_xJets_[Nedges-1];

	double eJER_raw_allJets_[Nedges-1], eJER_raw_bJets_[Nedges-1], eJER_raw_cJets_[Nedges-1], eJER_raw_udJets_[Nedges-1], eJER_raw_sJets_[Nedges-1], eJER_raw_gJets_[Nedges-1], eJER_raw_xJets_[Nedges-1];	
	double eJES_raw_allJets_[Nedges-1], eJES_raw_bJets_[Nedges-1], eJES_raw_cJets_[Nedges-1], eJES_raw_udJets_[Nedges-1], eJES_raw_sJets_[Nedges-1], eJES_raw_gJets_[Nedges-1], eJES_raw_xJets_[Nedges-1];


	TF1 *fitFxn_raw_allJets = new TF1("fitFxn_raw_allJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_bJets = new TF1("fitFxn_raw_bJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_cJets = new TF1("fitFxn_raw_cJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_udJets = new TF1("fitFxn_raw_udJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_sJets = new TF1("fitFxn_raw_sJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_gJets = new TF1("fitFxn_raw_gJets",fitfxn,0,300,3);
	TF1 *fitFxn_raw_xJets = new TF1("fitFxn_raw_xJets",fitfxn,0,300,3);
	


	fitFxn_raw_allJets->SetParameter(1,1.0);
	fitFxn_raw_bJets->SetParameter(1,1.0);
	fitFxn_raw_cJets->SetParameter(1,1.0);
	fitFxn_raw_udJets->SetParameter(1,1.0);
	fitFxn_raw_sJets->SetParameter(1,1.0);
	fitFxn_raw_gJets->SetParameter(1,1.0);
	fitFxn_raw_xJets->SetParameter(1,1.0);
	
	fitFxn_raw_allJets->SetParameter(2,0.3);
	fitFxn_raw_bJets->SetParameter(2,0.3);
	fitFxn_raw_cJets->SetParameter(2,0.3);
	fitFxn_raw_udJets->SetParameter(2,0.3);
	fitFxn_raw_sJets->SetParameter(2,0.3);
	fitFxn_raw_gJets->SetParameter(2,0.3);
	fitFxn_raw_xJets->SetParameter(2,0.3);

	double JER_tempVal_raw_allJets, JER_tempVal_raw_bJets, JER_tempVal_raw_cJets, JER_tempVal_raw_udJets, JER_tempVal_raw_sJets, JER_tempVal_raw_gJets, JER_tempVal_raw_xJets;
	double eJER_tempVal_raw_allJets, eJER_tempVal_raw_bJets, eJER_tempVal_raw_cJets, eJER_tempVal_raw_udJets, eJER_tempVal_raw_sJets, eJER_tempVal_raw_gJets, eJER_tempVal_raw_xJets;

	double JES_tempVal_raw_allJets, JES_tempVal_raw_bJets, JES_tempVal_raw_cJets, JES_tempVal_raw_udJets, JES_tempVal_raw_sJets, JES_tempVal_raw_gJets, JES_tempVal_raw_xJets;
	double eJES_tempVal_raw_allJets, eJES_tempVal_raw_bJets, eJES_tempVal_raw_cJets, eJES_tempVal_raw_udJets, eJES_tempVal_raw_sJets, eJES_tempVal_raw_gJets, eJES_tempVal_raw_xJets;


	TH1D *projX_corr_allJets, *projX_corr_bJets, *projX_corr_cJets, *projX_corr_udJets, *projX_corr_sJets, *projX_corr_gJets, *projX_corr_xJets;
	TH1D *projX_raw_allJets, *projX_raw_bJets, *projX_raw_cJets, *projX_raw_udJets, *projX_raw_sJets, *projX_raw_gJets, *projX_raw_xJets;


	for(int i = 0; i < Nedges-1; i++){
	
        	x[i] = (pt_axis[i+1]+pt_axis[i])/2;
        	ept_axis[i] = (pt_axis[i+1] - pt_axis[i])/2;

        	cout << "i = " << i << ", x[i] = " << x[i] << endl;

        	firstxbin = xaxis->FindBin(pt_axis[i]+0.001);
        	lastxbin = xaxis->FindBin(pt_axis[i+1]-0.001);

		projX_corr_allJets = response_corr_allJets->ProjectionX("projX_corr_allJets",firstxbin,lastxbin);	
		projX_corr_bJets = response_corr_bJets->ProjectionX("projX_corr_bJets",firstxbin,lastxbin);	
		projX_corr_cJets = response_corr_cJets->ProjectionX("projX_corr_cJets",firstxbin,lastxbin);	
		projX_corr_udJets = response_corr_udJets->ProjectionX("projX_corr_udJets",firstxbin,lastxbin);	
		projX_corr_sJets = response_corr_sJets->ProjectionX("projX_corr_sJets",firstxbin,lastxbin);	
		projX_corr_gJets = response_corr_gJets->ProjectionX("projX_corr_gJets",firstxbin,lastxbin);	
		projX_corr_xJets = response_corr_xJets->ProjectionX("projX_corr_xJets",firstxbin,lastxbin);	
	
		projX_raw_allJets = response_raw_allJets->ProjectionX("projX_raw_allJets",firstxbin,lastxbin);	
		projX_raw_bJets = response_raw_bJets->ProjectionX("projX_raw_bJets",firstxbin,lastxbin);	
		projX_raw_cJets = response_raw_cJets->ProjectionX("projX_raw_cJets",firstxbin,lastxbin);	
		projX_raw_udJets = response_raw_udJets->ProjectionX("projX_raw_udJets",firstxbin,lastxbin);	
		projX_raw_sJets = response_raw_sJets->ProjectionX("projX_raw_sJets",firstxbin,lastxbin);	
		projX_raw_gJets = response_raw_gJets->ProjectionX("projX_raw_gJets",firstxbin,lastxbin);	
		projX_raw_xJets = response_raw_xJets->ProjectionX("projX_raw_xJets",firstxbin,lastxbin);
	
		

		projX_corr_allJets->Fit(fitFxn_corr_allJets,"R","",0,5);
		projX_corr_bJets->Fit(fitFxn_corr_bJets,"R","",0,5);
		projX_corr_cJets->Fit(fitFxn_corr_cJets,"R","",0,5);
		projX_corr_udJets->Fit(fitFxn_corr_udJets,"R","",0,5);
		projX_corr_sJets->Fit(fitFxn_corr_sJets,"R","",0,5);
		projX_corr_gJets->Fit(fitFxn_corr_gJets,"R","",0,5);
		projX_corr_xJets->Fit(fitFxn_corr_xJets,"R","",0,5);

		JES_tempVal_corr_allJets = fitFxn_corr_allJets->GetParameter(1);
		eJES_tempVal_corr_allJets = fitFxn_corr_allJets->GetParError(1);
		JER_tempVal_corr_allJets = fitFxn_corr_allJets->GetParameter(2);
		eJER_tempVal_corr_allJets = fitFxn_corr_allJets->GetParError(2);

		JES_tempVal_corr_bJets = fitFxn_corr_bJets->GetParameter(1);
		eJES_tempVal_corr_bJets = fitFxn_corr_bJets->GetParError(1);
		JER_tempVal_corr_bJets = fitFxn_corr_bJets->GetParameter(2);
		eJER_tempVal_corr_bJets = fitFxn_corr_bJets->GetParError(2);

		JES_tempVal_corr_cJets = fitFxn_corr_cJets->GetParameter(1);
		eJES_tempVal_corr_cJets = fitFxn_corr_cJets->GetParError(1);
		JER_tempVal_corr_cJets = fitFxn_corr_cJets->GetParameter(2);
		eJER_tempVal_corr_cJets = fitFxn_corr_cJets->GetParError(2);

		JES_tempVal_corr_udJets = fitFxn_corr_udJets->GetParameter(1);
		eJES_tempVal_corr_udJets = fitFxn_corr_udJets->GetParError(1);
		JER_tempVal_corr_udJets = fitFxn_corr_udJets->GetParameter(2);
		eJER_tempVal_corr_udJets = fitFxn_corr_udJets->GetParError(2);

		JES_tempVal_corr_sJets = fitFxn_corr_sJets->GetParameter(1);
		eJES_tempVal_corr_sJets = fitFxn_corr_sJets->GetParError(1);
		JER_tempVal_corr_sJets = fitFxn_corr_sJets->GetParameter(2);
		eJER_tempVal_corr_sJets = fitFxn_corr_sJets->GetParError(2);

		JES_tempVal_corr_gJets = fitFxn_corr_gJets->GetParameter(1);
		eJES_tempVal_corr_gJets = fitFxn_corr_gJets->GetParError(1);
		JER_tempVal_corr_gJets = fitFxn_corr_gJets->GetParameter(2);
		eJER_tempVal_corr_gJets = fitFxn_corr_gJets->GetParError(2);

		JES_tempVal_corr_xJets = fitFxn_corr_xJets->GetParameter(1);
		eJES_tempVal_corr_xJets = fitFxn_corr_xJets->GetParError(1);
		JER_tempVal_corr_xJets = fitFxn_corr_xJets->GetParameter(2);
		eJER_tempVal_corr_xJets = fitFxn_corr_xJets->GetParError(2);

	
		projX_raw_allJets->Fit(fitFxn_raw_allJets,"R","",0,5);
                projX_raw_bJets->Fit(fitFxn_raw_bJets,"R","",0,5);
                projX_raw_cJets->Fit(fitFxn_raw_cJets,"R","",0,5);
                projX_raw_udJets->Fit(fitFxn_raw_udJets,"R","",0,5);
                projX_raw_sJets->Fit(fitFxn_raw_sJets,"R","",0,5);
                projX_raw_gJets->Fit(fitFxn_raw_gJets,"R","",0,5);
                projX_raw_xJets->Fit(fitFxn_raw_xJets,"R","",0,5);

                JES_tempVal_raw_allJets = fitFxn_raw_allJets->GetParameter(1);
                eJES_tempVal_raw_allJets = fitFxn_raw_allJets->GetParError(1);
                JER_tempVal_raw_allJets = fitFxn_raw_allJets->GetParameter(2);
                eJER_tempVal_raw_allJets = fitFxn_raw_allJets->GetParError(2);

                JES_tempVal_raw_bJets = fitFxn_raw_bJets->GetParameter(1);
                eJES_tempVal_raw_bJets = fitFxn_raw_bJets->GetParError(1);
                JER_tempVal_raw_bJets = fitFxn_raw_bJets->GetParameter(2);
                eJER_tempVal_raw_bJets = fitFxn_raw_bJets->GetParError(2);

                JES_tempVal_raw_cJets = fitFxn_raw_cJets->GetParameter(1);
                eJES_tempVal_raw_cJets = fitFxn_raw_cJets->GetParError(1);
                JER_tempVal_raw_cJets = fitFxn_raw_cJets->GetParameter(2);
                eJER_tempVal_raw_cJets = fitFxn_raw_cJets->GetParError(2);

                JES_tempVal_raw_udJets = fitFxn_raw_udJets->GetParameter(1);
                eJES_tempVal_raw_udJets = fitFxn_raw_udJets->GetParError(1);
                JER_tempVal_raw_udJets = fitFxn_raw_udJets->GetParameter(2);
                eJER_tempVal_raw_udJets = fitFxn_raw_udJets->GetParError(2);

                JES_tempVal_raw_sJets = fitFxn_raw_sJets->GetParameter(1);
                eJES_tempVal_raw_sJets = fitFxn_raw_sJets->GetParError(1);
                JER_tempVal_raw_sJets = fitFxn_raw_sJets->GetParameter(2);
                eJER_tempVal_raw_sJets = fitFxn_raw_sJets->GetParError(2);

                JES_tempVal_raw_gJets = fitFxn_raw_gJets->GetParameter(1);
                eJES_tempVal_raw_gJets = fitFxn_raw_gJets->GetParError(1);
                JER_tempVal_raw_gJets = fitFxn_raw_gJets->GetParameter(2);
                eJER_tempVal_raw_gJets = fitFxn_raw_gJets->GetParError(2);

                JES_tempVal_raw_xJets = fitFxn_raw_xJets->GetParameter(1);
                eJES_tempVal_raw_xJets = fitFxn_raw_xJets->GetParError(1);
                JER_tempVal_raw_xJets = fitFxn_raw_xJets->GetParameter(2);
                eJER_tempVal_raw_xJets = fitFxn_raw_xJets->GetParError(2);


		// fill input arrays with temp values
		JER_corr_allJets_[i] = JER_tempVal_corr_allJets;
		JER_corr_bJets_[i] = JER_tempVal_corr_bJets;
		JER_corr_cJets_[i] = JER_tempVal_corr_cJets;
		JER_corr_udJets_[i] = JER_tempVal_corr_udJets;
		JER_corr_sJets_[i] = JER_tempVal_corr_sJets;
		JER_corr_gJets_[i] = JER_tempVal_corr_gJets;
		JER_corr_xJets_[i] = JER_tempVal_corr_xJets;

		eJER_corr_allJets_[i] = eJER_tempVal_corr_allJets;
		eJER_corr_bJets_[i] = eJER_tempVal_corr_bJets;
		eJER_corr_cJets_[i] = eJER_tempVal_corr_cJets;
		eJER_corr_udJets_[i] = eJER_tempVal_corr_udJets;
		eJER_corr_sJets_[i] = eJER_tempVal_corr_sJets;
		eJER_corr_gJets_[i] = eJER_tempVal_corr_gJets;
		eJER_corr_xJets_[i] = eJER_tempVal_corr_xJets;

		JES_corr_allJets_[i] = JES_tempVal_corr_allJets;
                JES_corr_bJets_[i] = JES_tempVal_corr_bJets;
                JES_corr_cJets_[i] = JES_tempVal_corr_cJets;
                JES_corr_udJets_[i] = JES_tempVal_corr_udJets;
                JES_corr_sJets_[i] = JES_tempVal_corr_sJets;
                JES_corr_gJets_[i] = JES_tempVal_corr_gJets;
                JES_corr_xJets_[i] = JES_tempVal_corr_xJets;

                eJES_corr_allJets_[i] = eJES_tempVal_corr_allJets;
                eJES_corr_bJets_[i] = eJES_tempVal_corr_bJets;
                eJES_corr_cJets_[i] = eJES_tempVal_corr_cJets;
                eJES_corr_udJets_[i] = eJES_tempVal_corr_udJets;
                eJES_corr_sJets_[i] = eJES_tempVal_corr_sJets;
                eJES_corr_gJets_[i] = eJES_tempVal_corr_gJets;
                eJES_corr_xJets_[i] = eJES_tempVal_corr_xJets;		
		
		
		JER_raw_allJets_[i] = JER_tempVal_raw_allJets;
                JER_raw_bJets_[i] = JER_tempVal_raw_bJets;
                JER_raw_cJets_[i] = JER_tempVal_raw_cJets;
                JER_raw_udJets_[i] = JER_tempVal_raw_udJets;
                JER_raw_sJets_[i] = JER_tempVal_raw_sJets;
                JER_raw_gJets_[i] = JER_tempVal_raw_gJets;
                JER_raw_xJets_[i] = JER_tempVal_raw_xJets;

                eJER_raw_allJets_[i] = eJER_tempVal_raw_allJets;
                eJER_raw_bJets_[i] = eJER_tempVal_raw_bJets;
                eJER_raw_cJets_[i] = eJER_tempVal_raw_cJets;
                eJER_raw_udJets_[i] = eJER_tempVal_raw_udJets;
                eJER_raw_sJets_[i] = eJER_tempVal_raw_sJets;
                eJER_raw_gJets_[i] = eJER_tempVal_raw_gJets;
                eJER_raw_xJets_[i] = eJER_tempVal_raw_xJets;

		JES_raw_allJets_[i] = JES_tempVal_raw_allJets;
                JES_raw_bJets_[i] = JES_tempVal_raw_bJets;
                JES_raw_cJets_[i] = JES_tempVal_raw_cJets;
                JES_raw_udJets_[i] = JES_tempVal_raw_udJets;
                JES_raw_sJets_[i] = JES_tempVal_raw_sJets;
                JES_raw_gJets_[i] = JES_tempVal_raw_gJets;
                JES_raw_xJets_[i] = JES_tempVal_raw_xJets;

                eJES_raw_allJets_[i] = eJES_tempVal_raw_allJets;
                eJES_raw_bJets_[i] = eJES_tempVal_raw_bJets;
                eJES_raw_cJets_[i] = eJES_tempVal_raw_cJets;
                eJES_raw_udJets_[i] = eJES_tempVal_raw_udJets;
                eJES_raw_sJets_[i] = eJES_tempVal_raw_sJets;
                eJES_raw_gJets_[i] = eJES_tempVal_raw_gJets;
                eJES_raw_xJets_[i] = eJES_tempVal_raw_xJets;







	}




	TGraphErrors *gr1 = new TGraphErrors(Nedges-1,x,JER_corr_allJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr2 = new TGraphErrors(Nedges-1,x,JER_corr_bJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr3 = new TGraphErrors(Nedges-1,x,JER_corr_cJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr4 = new TGraphErrors(Nedges-1,x,JER_corr_udJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr5 = new TGraphErrors(Nedges-1,x,JER_corr_gJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr6 = new TGraphErrors(Nedges-1,x,JER_corr_sJets_,ept_axis,eJER_corr_allJets_);
	TGraphErrors *gr7 = new TGraphErrors(Nedges-1,x,JER_corr_xJets_,ept_axis,eJER_corr_allJets_);

	TGraphErrors *Gr1 = new TGraphErrors(Nedges-1,x,JES_corr_allJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr2 = new TGraphErrors(Nedges-1,x,JES_corr_bJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr3 = new TGraphErrors(Nedges-1,x,JES_corr_cJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr4 = new TGraphErrors(Nedges-1,x,JES_corr_udJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr5 = new TGraphErrors(Nedges-1,x,JES_corr_gJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr6 = new TGraphErrors(Nedges-1,x,JES_corr_sJets_,ept_axis,eJES_corr_allJets_);
	TGraphErrors *Gr7 = new TGraphErrors(Nedges-1,x,JES_corr_xJets_,ept_axis,eJES_corr_allJets_);

	gr1->SetLineColor(kBlack);
	gr2->SetLineColor(kBlue);
	gr3->SetLineColor(kRed-2);
	gr4->SetLineColor(kGreen+2);
	gr5->SetLineColor(kMagenta+2);
	gr6->SetLineColor(kCyan+1);
	gr7->SetLineColor(kBlue+6);

	Gr1->SetLineColor(kBlack);
	Gr2->SetLineColor(kBlue);
	Gr3->SetLineColor(kRed-2);
	Gr4->SetLineColor(kGreen+2);
	Gr5->SetLineColor(kMagenta+2);
	Gr6->SetLineColor(kCyan+1);
	Gr7->SetLineColor(kBlue+6);

	
	double legendTextSize = 0.035;
	TLegend *l = new TLegend(0.5,0.6,0.89,0.89);
	l->SetBorderSize(0);
	l->SetTextSize(legendTextSize);
	l->AddEntry(gr2,"b");
	l->AddEntry(gr3,"c");
	l->AddEntry(gr4,"ud");
	l->AddEntry(gr5,"s");
	l->AddEntry(gr6,"g");
	l->AddEntry(gr7,"x");
	
	
	TCanvas *c1 = new TCanvas("c1","c1",700,700);
   	c1->cd();
   	TPad *p1 = new TPad("p1","p1",0,0,1,1);
   	p1->SetLeftMargin(0.2);
   	p1->SetBottomMargin(0.15);
   	p1->Draw();
   	p1->cd();
   	TMultiGraph  *mg1  = new TMultiGraph();

   	mg1->Add(gr2);
   	mg1->Add(gr3);
   	mg1->Add(gr4);
	mg1->Add(gr5);	
	mg1->Add(gr6);	
	mg1->Add(gr7);

	mg1->GetYaxis()->SetRangeUser(0.0,0.4);
   	mg1->GetXaxis()->SetRangeUser(50,300);
   	mg1->GetYaxis()->SetTitle("#sigma (#font[52]{p}_{T}^{recoJet} / #font[52]{p}_{T}^{genJet})");
   	mg1->GetXaxis()->SetTitle("#font[52]{p}_{T}^{genJet} [GeV]");
   	mg1->GetXaxis()->SetLabelSize(0.04);
   	mg1->GetXaxis()->SetTitleSize(0.05);
   	mg1->GetYaxis()->SetLabelSize(0.04);
   	mg1->GetYaxis()->SetTitleSize(0.05);
   	mg1->GetXaxis()->SetTitleFont(42);
   	mg1->GetXaxis()->SetLabelFont(42);
   	mg1->GetYaxis()->SetTitleFont(42);
   	mg1->GetYaxis()->SetLabelFont(42);

	mg1->Draw("AP");	
	l->Draw();





	TCanvas *c2 = new TCanvas("c2","c2",700,700);
   	c2->cd();
   	TPad *p2 = new TPad("p2","p2",0,0,1,1);
   	p2->SetLeftMargin(0.2);
   	p2->SetBottomMargin(0.15);
   	p2->Draw();
   	p2->cd();
   	TMultiGraph  *mg2  = new TMultiGraph();

   	mg2->Add(Gr2);
   	mg2->Add(Gr3);
   	mg2->Add(Gr4);
	mg2->Add(Gr5);	
	mg2->Add(Gr6);	
	mg2->Add(Gr7);

	mg2->GetYaxis()->SetRangeUser(0.8,1.2);
   	mg2->GetXaxis()->SetRangeUser(50,300);
   	mg2->GetYaxis()->SetTitle("#mu (#font[52]{p}_{T}^{recoJet} / #font[52]{p}_{T}^{genJet})");
   	mg2->GetXaxis()->SetTitle("#font[52]{p}_{T}^{genJet} [GeV]");
   	mg2->GetXaxis()->SetLabelSize(0.04);
   	mg2->GetXaxis()->SetTitleSize(0.05);
   	mg2->GetYaxis()->SetLabelSize(0.04);
   	mg2->GetYaxis()->SetTitleSize(0.05);
   	mg2->GetXaxis()->SetTitleFont(42);
   	mg2->GetXaxis()->SetLabelFont(42);
   	mg2->GetYaxis()->SetTitleFont(42);
   	mg2->GetYaxis()->SetLabelFont(42);

	mg2->Draw("AP");	











}
