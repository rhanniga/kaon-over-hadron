void makeUSCorrections(string inputFile){
    TFile* input = new TFile(inputFile.c_str());
    TH2D* hKaon2Dpeak = (TH2D*)input->Get("hKaon2Dpeak");
    TH2D* hKaon2DLside = (TH2D*)input->Get("hKaon2DLside");
    TH2D* hKaon2DRside = (TH2D*)input->Get("hKaon2DRside");
    TH2D* hKaonLS2Dpeak = (TH2D*)input->Get("hKaonLS2Dpeak");
    TH2D* hKaonLS2DLside = (TH2D*)input->Get("hKaonLS2DLside");
    TH2D* hKaonLS2DRside = (TH2D*)input->Get("hKaonLS2DRside");

    TH2D* trigDistSameUS = (TH2D*)input->Get("fTrigSameUSDist");
    TH2D* trigDistSameLS = (TH2D*)input->Get("fTrigSameLSDist");

    hKaon2Dpeak->SetName("uncorrectedhKaon2Dpeak");

    TH2D* hKaonBGPeakRegionL = (TH2D*)hKaon2DLside->Clone("hKaonBGPeakRegionL");
    hKaonBGPeakRegionL->Scale(1.0/(hKaon2DLside->Integral(hKaon2DLside->GetXaxis()->FindBin(-1.2), hKaon2DLside->GetXaxis()->FindBin(1.2), 1, hKaon2DLside->GetYaxis()->GetNbins())));
    hKaonBGPeakRegionL->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hKaonBGPeakRegionL_deta = (TH1D*)hKaonBGPeakRegionL->ProjectionX("hKaonBGPeakRegionL_deta", 1, hKaonBGPeakRegionL->GetYaxis()->GetNbins());
    TH1D* hKaonBGPeakRegionL_dphi = (TH1D*)hKaonBGPeakRegionL->ProjectionY("hKaonBGPeakRegionL_dphi", hKaonBGPeakRegionL->GetXaxis()->FindBin(-1.2), hKaonBGPeakRegionL->GetXaxis()->FindBin(1.2));

    TH2D* hKaonBGPeakRegionR = (TH2D*)hKaon2DRside->Clone("hKaonBGPeakRegionR");
    hKaonBGPeakRegionR->Scale(1.0/(hKaon2DRside->Integral(hKaon2DRside->GetXaxis()->FindBin(-1.2), hKaon2DRside->GetXaxis()->FindBin(1.2), 1, hKaon2DRside->GetYaxis()->GetNbins())));
    hKaonBGPeakRegionR->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hKaonBGPeakRegionR_deta = (TH1D*)hKaonBGPeakRegionR->ProjectionX("hKaonBGPeakRegionR_deta", 1, hKaonBGPeakRegionR->GetYaxis()->GetNbins());
    TH1D* hKaonBGPeakRegionR_dphi = (TH1D*)hKaonBGPeakRegionR->ProjectionY("hKaonBGPeakRegionR_dphi", hKaonBGPeakRegionR->GetXaxis()->FindBin(-1.2), hKaonBGPeakRegionR->GetXaxis()->FindBin(1.2));

    TH2D* hKaonBGPeakRegion = (TH2D*)hKaonBGPeakRegionL->Clone("hKaonBGPeakregion");
    hKaonBGPeakRegion->Add(hKaonBGPeakRegionR);
    hKaonBGPeakRegion->Scale(0.5);
    hKaonBGPeakRegion->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* hKaonBGPeakRegion_deta = (TH1D*)hKaonBGPeakRegion->ProjectionX("hKaonBGPeakRegion_deta", 1, hKaonBGPeakRegion->GetYaxis()->GetNbins());
    TH1D* hKaonBGPeakRegion_dphi = (TH1D*)hKaonBGPeakRegion->ProjectionY("hKaonBGPeakRegion_dphi", hKaonBGPeakRegion->GetXaxis()->FindBin(-1.2), hKaonBGPeakRegion->GetXaxis()->FindBin(1.2));


    //US residual checks between SB average and the Left and Right separately
    TH2D* resLeftVsAvg = (TH2D*)hKaonBGPeakRegionL->Clone("resLeftVsAvg");
    resLeftVsAvg->Add(hKaonBGPeakRegion, -1.0);
    resLeftVsAvg->Divide(hKaonBGPeakRegionL);
    resLeftVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_deta = (TH1D*)hKaonBGPeakRegionL_deta->Clone("resLeftVsAvg_deta");
    resLeftVsAvg_deta->Add(hKaonBGPeakRegion_deta, -1.0);
    resLeftVsAvg_deta->Divide(hKaonBGPeakRegionL_deta);
    resLeftVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resLeftVsAvg_dphi = (TH1D*)hKaonBGPeakRegionL_dphi->Clone("resLeftVsAvg_dphi");
    resLeftVsAvg_dphi->Add(hKaonBGPeakRegion_dphi, -1.0);
    resLeftVsAvg_dphi->Divide(hKaonBGPeakRegionL_dphi);

    TH2D* resRightVsAvg = (TH2D*)hKaonBGPeakRegionR->Clone("resRightVsAbg");
    resRightVsAvg->Add(hKaonBGPeakRegion, -1.0);
    resRightVsAvg->Divide(hKaonBGPeakRegionR);
    resRightVsAvg->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_deta = (TH1D*)hKaonBGPeakRegionR_deta->Clone("resRightVsAvg_deta");
    resRightVsAvg_deta->Add(hKaonBGPeakRegion_deta, -1.0);
    resRightVsAvg_deta->Divide(hKaonBGPeakRegionR_deta);
    resRightVsAvg_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resRightVsAvg_dphi = (TH1D*)hKaonBGPeakRegionR_dphi->Clone("resRightVsAvg_dphi");
    resRightVsAvg_dphi->Add(hKaonBGPeakRegion_dphi, -1.0);
    resRightVsAvg_dphi->Divide(hKaonBGPeakRegionR_dphi);



    Float_t leftscale = hKaon2DLside->Integral(hKaon2DLside->GetXaxis()->FindBin(-1.2), hKaon2DLside->GetXaxis()->FindBin(1.2), 1, hKaon2DLside->GetYaxis()->GetNbins())/hKaonLS2DLside->Integral(hKaon2DLside->GetXaxis()->FindBin(-1.2), hKaon2DLside->GetXaxis()->FindBin(1.2), 1, hKaon2DLside->GetYaxis()->GetNbins());

    TH2D* LLSsubhKaon2DLside = (TH2D*)hKaon2DLside->Clone("LLSsubhKaon2DLside");
    TH2D* LLSsubhKaon2Dpeak = (TH2D*)hKaon2Dpeak->Clone("LLSsubhKaon2Dpeak");
    LLSsubhKaon2DLside->Add(hKaonLS2DLside, -1.0*leftscale);
    //LLSsubhKaon2DLside->Divide(hKaonLS2DLside);
    //LLSsubhKaon2DLside->Scale(1.0/leftscale);
    LLSsubhKaon2Dpeak->Add(hKaonLS2Dpeak, -1.0*leftscale);

    TH1D* LLSsubhKaon2DLside_deta = LLSsubhKaon2DLside->ProjectionX("LLSsubhKaon2DLside_deta", 1, LLSsubhKaon2DLside->GetYaxis()->GetNbins());
    TH1D* LLSsubhKaon2DLside_dphi = LLSsubhKaon2DLside->ProjectionY("LLSsubhKaon2DLside_dphi", LLSsubhKaon2DLside->GetXaxis()->FindBin(-1.2), LLSsubhKaon2DLside->GetXaxis()->FindBin(1.2));

    //Float_t rightscale = hKaon2DRside->Integral(1, hKaon2DRside->GetXaxis()->GetNbins(), 1, hKaon2DRside->GetYaxis()->GetNbins())/hKaonLS2DRside->Integral(1, hKaon2DRside->GetXaxis()->GetNbins(), 1, hKaon2DRside->GetYaxis()->GetNbins());
    gStyle->SetOptStat(0);
    Float_t rightscale = hKaon2DRside->Integral(hKaon2DRside->GetXaxis()->FindBin(-1.2), hKaon2DRside->GetXaxis()->FindBin(1.2), 1, hKaon2DRside->GetYaxis()->GetNbins())/hKaonLS2DRside->Integral(hKaon2DRside->GetXaxis()->FindBin(-1.2), hKaon2DRside->GetXaxis()->FindBin(1.2), 1, hKaon2DRside->GetYaxis()->GetNbins());
    TH2D* RLSsubhKaon2DRside = (TH2D*)hKaon2DRside->Clone("RLSsubhKaon2DRside");
    TH2D* RLSsubhKaon2Dpeak = (TH2D*)hKaon2Dpeak->Clone("RLSsubhKaon2Dpeak");
    TH2D* RLSsubhKaon2DRsideScaled = (TH2D*)hKaon2DRside->Clone("RLSsubhKaon2DRsideScaled");
    RLSsubhKaon2DRsideScaled->GetXaxis()->SetRangeUser(-1.2, 1.2);
    RLSsubhKaon2DRsideScaled->Scale(rightscale);
    TH1D* RLSsubhKaonDPhiRsideScaled = RLSsubhKaon2DRsideScaled->ProjectionY();
    RLSsubhKaonDPhiRsideScaled->SetTitle("h-#Kaon^{0} scaled R-sideband #Delta#varphi distribution");
    RLSsubhKaonDPhiRsideScaled->SetLineColor(6);
    RLSsubhKaonDPhiRsideScaled->SetLineWidth(3);
    RLSsubhKaonDPhiRsideScaled->SetMarkerColor(6);
    TCanvas *presCanvas = new TCanvas("presCanvas", "Presentation Canvas", 0, 10, 1600, 1200);
    presCanvas->cd();
    RLSsubhKaonDPhiRsideScaled->Draw();

    RLSsubhKaon2DRside->Add(hKaonLS2DRside, -1.0*rightscale);
    //RLSsubhKaon2DRside->Divide(hKaonLS2DRside);
    //RLSsubhKaon2DRside->Scale(1.0/rightscale);
    RLSsubhKaon2Dpeak->Add(hKaonLS2Dpeak, -1.0*rightscale);
    RLSsubhKaon2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);

    TH1D* RLSsubhKaon2DRside_deta = RLSsubhKaon2DRside->ProjectionX("RLSsubhKaon2DRside_deta", 1, RLSsubhKaon2DRside->GetYaxis()->GetNbins());
    TH1D* RLSsubhKaon2DRside_dphi = RLSsubhKaon2DRside->ProjectionY("RLSsubhKaon2DRside_dphi", RLSsubhKaon2DRside->GetXaxis()->FindBin(-1.2), RLSsubhKaon2DRside->GetXaxis()->FindBin(1.2));

    TH1D* RLSsubhKaon2Dpeak_deta = RLSsubhKaon2Dpeak->ProjectionX("RLSsubhKaon2Dpeak_deta", 1, RLSsubhKaon2Dpeak->GetYaxis()->GetNbins());
    TH1D* RLSsubhKaon2Dpeak_dphi = RLSsubhKaon2Dpeak->ProjectionY("RLSsubhKaon2Dpeak_dphi", RLSsubhKaon2Dpeak->GetXaxis()->FindBin(-1.2), RLSsubhKaon2Dpeak->GetXaxis()->FindBin(1.2));


    TH1D* scales = new TH1D("scales", "scales", 2, -1, 1);
    scales->SetBinContent(1, leftscale);
    scales->SetBinContent(2, rightscale);

    TH2D* rebinRLSsubhKaon2Dpeak = (TH2D*)RLSsubhKaon2Dpeak->Clone("rebinRLSsubhKaon2Dpeak");
    rebinRLSsubhKaon2Dpeak->Rebin2D(2, 2);

    //Using US estimate for BG to subtract off the from the peak region:

    Float_t scaleUS = (rightscale)*hKaonLS2Dpeak->Integral(hKaonLS2Dpeak->GetXaxis()->FindBin(-1.2), hKaonLS2Dpeak->GetXaxis()->FindBin(1.2), 1, hKaonLS2Dpeak->GetYaxis()->GetNbins());
    Float_t scaletest = (rightscale)*hKaon2Dpeak->Integral(hKaon2Dpeak->GetXaxis()->FindBin(-1.2), hKaon2Dpeak->GetXaxis()->FindBin(1.2), 1, hKaon2Dpeak->GetYaxis()->GetNbins());


    printf("\n\nscaleUS = %e\n\ntestscale = %e \n\n", scaleUS, scaletest);

    //avg of right and left US sideband tests
    TH2D* AvgUSsubhKaon2Dpeak = (TH2D*)hKaon2Dpeak->Clone("AvgUSsubhKaon2Dpeak");
    AvgUSsubhKaon2Dpeak->Add(hKaonBGPeakRegion, -1.0*scaleUS);
    AvgUSsubhKaon2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhKaon2Dpeak_deta = (TH1D*)AvgUSsubhKaon2Dpeak->ProjectionX("AvgUSsubhKaon2Dpeak_deta", 1, AvgUSsubhKaon2Dpeak->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhKaon2Dpeak_dphi = (TH1D*)AvgUSsubhKaon2Dpeak->ProjectionY("AvgUSsubhKaon2Dpeak_dphi", AvgUSsubhKaon2Dpeak->GetXaxis()->FindBin(-1.2), AvgUSsubhKaon2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhKaon2Dpeakleftscale = (TH2D*)hKaon2Dpeak->Clone("AvgUSsubhKaon2Dpeakleftscale");
    AvgUSsubhKaon2Dpeakleftscale->Add(hKaonBGPeakRegion, -1.0*scaleUS*leftscale/rightscale);
    AvgUSsubhKaon2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhKaon2Dpeakleftscale_deta = (TH1D*)AvgUSsubhKaon2Dpeakleftscale->ProjectionX("AvgUSsubhKaon2Dpeakleftscale_deta", 1, AvgUSsubhKaon2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhKaon2Dpeakleftscale_dphi = (TH1D*)AvgUSsubhKaon2Dpeakleftscale->ProjectionY("AvgUSsubhKaon2Dpeakleftscale_dphi", AvgUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(-1.2), AvgUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* AvgUSsubhKaon2Dpeakavgscale = (TH2D*)hKaon2Dpeak->Clone("AvgUSsubhKaon2Dpeakavgscale");
    AvgUSsubhKaon2Dpeakavgscale->Add(hKaonBGPeakRegion, -1.0*scaleUS*(leftscale + rightscale)/(2.0*rightscale));
    AvgUSsubhKaon2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* AvgUSsubhKaon2Dpeakavgscale_deta = (TH1D*)AvgUSsubhKaon2Dpeakavgscale->ProjectionX("AvgUSsubhKaon2Dpeakavgscale_deta", 1, AvgUSsubhKaon2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* AvgUSsubhKaon2Dpeakavgscale_dphi = (TH1D*)AvgUSsubhKaon2Dpeakavgscale->ProjectionY("AvgUSsubhKaon2Dpeakavgscale_dphi", AvgUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(-1.2), AvgUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //right side US sideband tests
    TH2D* RSUSsubhKaon2Dpeak = (TH2D*)hKaon2Dpeak->Clone("RSUSsubhKaon2Dpeak");
    RSUSsubhKaon2Dpeak->Add(hKaonBGPeakRegionR, -1.0*scaleUS);
    RSUSsubhKaon2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhKaon2Dpeak_deta = (TH1D*)RSUSsubhKaon2Dpeak->ProjectionX("RSUSsubhKaon2Dpeak_deta", 1, RSUSsubhKaon2Dpeak->GetYaxis()->GetNbins());
    TH1D* RSUSsubhKaon2Dpeak_dphi = (TH1D*)RSUSsubhKaon2Dpeak->ProjectionY("RSUSsubhKaon2Dpeak_dphi", RSUSsubhKaon2Dpeak->GetXaxis()->FindBin(-1.2), RSUSsubhKaon2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* RSUSsubhKaon2Dpeakleftscale = (TH2D*)hKaon2Dpeak->Clone("RSUSsubhKaon2Dpeakleftscale");
    RSUSsubhKaon2Dpeakleftscale->Add(hKaonBGPeakRegionR, -1.0*scaleUS*leftscale/rightscale);
    RSUSsubhKaon2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhKaon2Dpeakleftscale_deta = (TH1D*)RSUSsubhKaon2Dpeakleftscale->ProjectionX("RSUSsubhKaon2Dpeakleftscale_deta", 1, RSUSsubhKaon2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhKaon2Dpeakleftscale_dphi = (TH1D*)RSUSsubhKaon2Dpeakleftscale->ProjectionY("RSUSsubhKaon2Dpeakleftscale_dphi", RSUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(-1.2), RSUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    //GOTO HERE
    TH2D* RSUSsubhKaon2Dpeakavgscale = (TH2D*)hKaon2Dpeak->Clone("RSUSsubhKaon2Dpeakavgscale");
    RSUSsubhKaon2Dpeakavgscale->Add(hKaonBGPeakRegionR, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    RSUSsubhKaon2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* RSUSsubhKaon2Dpeakavgscale_deta = (TH1D*)RSUSsubhKaon2Dpeakavgscale->ProjectionX("RSUSsubhKaon2Dpeakavgscale_deta", 1, RSUSsubhKaon2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* RSUSsubhKaon2Dpeakavgscale_dphi = (TH1D*)RSUSsubhKaon2Dpeakavgscale->ProjectionY("RSUSsubhKaon2Dpeakavgscale_dphi", RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(1.2));
    //END GOTO

    // double signalOverTotalArray[3] = {0.09745, 0.1191, 0.1697}; // have to change this for each multiplicity bin

    // TH2D* RSUSsubhKaon2Dpeakavgscale = (TH2D*)hKaon2Dpeak->Clone("RSUSsubhKaon2Dpeakavgscale");

    // // Using BG = TOTAL - SIGNAL = TOTAL(1-S/TOTAL)
    // double bgIntegral = RSUSsubhKaon2Dpeakavgscale->Integral(RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(1.2), 1, RSUSsubhKaon2Dpeakavgscale->GetYaxis()->GetNbins())*(1 - signalOverTotalArray[2]);

    // RSUSsubhKaon2Dpeakavgscale->Add(hKaonBGPeakRegionR, -1.0*bgIntegral);
    // RSUSsubhKaon2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    // TH1D* RSUSsubhKaon2Dpeakavgscale_deta = (TH1D*)RSUSsubhKaon2Dpeakavgscale->ProjectionX("RSUSsubhKaon2Dpeakavgscale_deta", 1, RSUSsubhKaon2Dpeakavgscale->GetYaxis()->GetNbins());
    // TH1D* RSUSsubhKaon2Dpeakavgscale_dphi = (TH1D*)RSUSsubhKaon2Dpeakavgscale->ProjectionY("RSUSsubhKaon2Dpeakavgscale_dphi", RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(-1.2), RSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(1.2));

    //left side US sideband tests
    TH2D* LSUSsubhKaon2Dpeak = (TH2D*)hKaon2Dpeak->Clone("LSUSsubhKaon2Dpeak");
    LSUSsubhKaon2Dpeak->Add(hKaonBGPeakRegionL, -1.0*scaleUS);
    LSUSsubhKaon2Dpeak->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhKaon2Dpeak_deta = (TH1D*)LSUSsubhKaon2Dpeak->ProjectionX("LSUSsubhKaon2Dpeak_deta", 1, LSUSsubhKaon2Dpeak->GetYaxis()->GetNbins());
    TH1D* LSUSsubhKaon2Dpeak_dphi = (TH1D*)LSUSsubhKaon2Dpeak->ProjectionY("LSUSsubhKaon2Dpeak_dphi", LSUSsubhKaon2Dpeak->GetXaxis()->FindBin(-1.2), LSUSsubhKaon2Dpeak->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhKaon2Dpeakleftscale = (TH2D*)hKaon2Dpeak->Clone("LSUSsubhKaon2Dpeakleftscale");
    LSUSsubhKaon2Dpeakleftscale->Add(hKaonBGPeakRegionL, -1.0*scaleUS*leftscale/rightscale);
    LSUSsubhKaon2Dpeakleftscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhKaon2Dpeakleftscale_deta = (TH1D*)LSUSsubhKaon2Dpeakleftscale->ProjectionX("LSUSsubhKaon2Dpeakleftscale_deta", 1, LSUSsubhKaon2Dpeakleftscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhKaon2Dpeakleftscale_dphi = (TH1D*)LSUSsubhKaon2Dpeakleftscale->ProjectionY("LSUSsubhKaon2Dpeakleftscale_dphi", LSUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(-1.2), LSUSsubhKaon2Dpeakleftscale->GetXaxis()->FindBin(1.2));

    TH2D* LSUSsubhKaon2Dpeakavgscale = (TH2D*)hKaon2Dpeak->Clone("LSUSsubhKaon2Dpeakavgscale");
    LSUSsubhKaon2Dpeakavgscale->Add(hKaonBGPeakRegionL, -1.0*scaleUS*(leftscale+rightscale)/(2.0*rightscale));
    LSUSsubhKaon2Dpeakavgscale->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* LSUSsubhKaon2Dpeakavgscale_deta = (TH1D*)LSUSsubhKaon2Dpeakavgscale->ProjectionX("LSUSsubhKaon2Dpeakavgscale_deta", 1, LSUSsubhKaon2Dpeakavgscale->GetYaxis()->GetNbins());
    TH1D* LSUSsubhKaon2Dpeakavgscale_dphi = (TH1D*)LSUSsubhKaon2Dpeakavgscale->ProjectionY("LSUSsubhKaon2Dpeakavgscale_dphi", LSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(-1.2), LSUSsubhKaon2Dpeakavgscale->GetXaxis()->FindBin(1.2));


    TH2D* resUSvsLS = (TH2D*)AvgUSsubhKaon2Dpeak->Clone("resUSvsLS");
    resUSvsLS->Add(RLSsubhKaon2Dpeak, -1.0);
    resUSvsLS->Divide(AvgUSsubhKaon2Dpeak);
    resUSvsLS->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_deta = (TH1D*)AvgUSsubhKaon2Dpeak_deta->Clone("resUSvsLS_deta");
    resUSvsLS_deta->Add(RLSsubhKaon2Dpeak_deta, -1.0);
    resUSvsLS_deta->Divide(AvgUSsubhKaon2Dpeak_deta);
    resUSvsLS_deta->GetXaxis()->SetRangeUser(-1.2, 1.2);
    TH1D* resUSvsLS_dphi = (TH1D*)AvgUSsubhKaon2Dpeak_dphi->Clone("resUSvsLS_dphi");
    resUSvsLS_dphi->Add(RLSsubhKaon2Dpeak_dphi, -1.0);
    resUSvsLS_dphi->Divide(AvgUSsubhKaon2Dpeak_dphi);



    TFile* output = new TFile(Form("US_syst_%s", inputFile.c_str()), "RECREATE");
    LLSsubhKaon2DLside->Write();
    LLSsubhKaon2DLside_deta->Write();
    LLSsubhKaon2DLside_dphi->Write();
    LLSsubhKaon2Dpeak->Write();
    RLSsubhKaon2DRside->Write();
    RLSsubhKaon2DRside_deta->Write();
    RLSsubhKaon2DRside_dphi->Write();
    RLSsubhKaon2Dpeak->Write();
    RLSsubhKaon2Dpeak_deta->Write();
    RLSsubhKaon2Dpeak_dphi->Write();
    rebinRLSsubhKaon2Dpeak->Write();
    scales->Write();
    hKaon2Dpeak->Write();
    hKaonBGPeakRegionL->Write();
    hKaonBGPeakRegionL_deta->Write();
    hKaonBGPeakRegionL_dphi->Write();
    hKaonBGPeakRegionR->Write();
    hKaonBGPeakRegionR_deta->Write();
    hKaonBGPeakRegionR_dphi->Write();
    hKaonBGPeakRegion->Write();
    hKaonBGPeakRegion_deta->Write();
    hKaonBGPeakRegion_dphi->Write();
    resLeftVsAvg->Write();
    resLeftVsAvg_deta->Write();
    resLeftVsAvg_dphi->Write();
    resRightVsAvg->Write();
    resRightVsAvg_deta->Write();
    resRightVsAvg_dphi->Write();
    AvgUSsubhKaon2Dpeak->Write();
    AvgUSsubhKaon2Dpeak_deta->Write();
    AvgUSsubhKaon2Dpeak_dphi->Write();
    AvgUSsubhKaon2Dpeakleftscale->Write();
    AvgUSsubhKaon2Dpeakleftscale_deta->Write();
    AvgUSsubhKaon2Dpeakleftscale_dphi->Write();
    AvgUSsubhKaon2Dpeakavgscale->Write();
    AvgUSsubhKaon2Dpeakavgscale_deta->Write();
    AvgUSsubhKaon2Dpeakavgscale_dphi->Write();
    RSUSsubhKaon2Dpeak->Write();
    RSUSsubhKaon2Dpeak_deta->Write();
    RSUSsubhKaon2Dpeak_dphi->Write();
    RSUSsubhKaon2Dpeakleftscale->Write();
    RSUSsubhKaon2Dpeakleftscale_deta->Write();
    RSUSsubhKaon2Dpeakleftscale_dphi->Write();
    RSUSsubhKaon2Dpeakavgscale->Write();
    RSUSsubhKaon2Dpeakavgscale_deta->Write();
    RSUSsubhKaon2Dpeakavgscale_dphi->Write();
    LSUSsubhKaon2Dpeak->Write();
    LSUSsubhKaon2Dpeak_deta->Write();
    LSUSsubhKaon2Dpeak_dphi->Write();
    LSUSsubhKaon2Dpeakleftscale->Write();
    LSUSsubhKaon2Dpeakleftscale_deta->Write();
    LSUSsubhKaon2Dpeakleftscale_dphi->Write();
    LSUSsubhKaon2Dpeakavgscale->Write();
    LSUSsubhKaon2Dpeakavgscale_deta->Write();
    LSUSsubhKaon2Dpeakavgscale_dphi->Write();
    resUSvsLS->Write();
    resUSvsLS_deta->Write();
    resUSvsLS_dphi->Write();
}


