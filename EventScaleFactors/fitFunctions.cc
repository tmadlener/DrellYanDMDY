#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../Include/fitFunctions.hh"
#include "../Include/MyTools.hh"
#include <TEntryList.h>
#endif

const int targetSpecificBin=0;
const int targetEt=0;
const int targetEta=1;

// The following systematic error is derived in a series of studies
// documented in this presentation:
//    https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=214643
// additionally, the numbers are bumped up to 1%, which is order
// of systematics for these fits coming from signal and background
// parameterization, as indicated by several studies presented at EGM. 
// There is no exact number on which EGM group signed off, so
// a generic 1% is used.
//    This set of systematics is applicable to 6 bins of Et.
// For other Et binnings, the array would be 3% below 20 GeV and 1% above. 
const double effRecoSystErr[6] = {0.030, 0.022, 0.01, 0.01, 0.01, 0.01};

// The following systematic error is related to the fact that
// we measure the average trigger efficiency of the leading and trailing
// electrons. The leading electron is always a subset of the
// trailing for present trigers, so the average of leading+trailing
// has the same value as the efficiency for trailing. However, the 
// HLT efficiency for leading, measured separately, is a bit different.
// The assigned systematics is the difference between the dedicated
// measurement of the trailing efficiency and the average, divided by two.
// (The division by two is because really the error is relevant only
// to the efficiency in data. However, the functions at present do not know
// whether it is data or MC. Assining half ot the error to each data and MC
// is an approimation, but will give nearly the desired result on the data/MC
// ratio.
const double effHltSystErrBarrel[6] = {0.0, 0.0, 0.006, 0.002, 0.002, 0.001};
const double effHltSystErrEndcap[6] = {0.0, 0.0, 0.010, 0.002, 0.002, 0.001}; // For endcap, the error is larger

/*
void measurePassAndFail(double &signal, double &signalErr, double &efficiency, double &efficiencyErr,TTree *passTree, TTree *failTree,TCanvas *passCanvas, TCanvas *failCanvas,const char* setBinsType){

  // Define data sets
  RooRealVar mass("mass","mass",60, 120);
  mass.setBins(120);
  RooRealVar pt ("et" ,"et" ,10.0, 1000);
  RooRealVar eta("eta","eta",-10, 10);
  RooDataSet *dataPass = new RooDataSet("dataPass","dataPass",passTree,RooArgSet(mass,pt,eta));
  RooDataSet *dataFail = new RooDataSet("dataFail","dataFail",RooArgList(mass,pt,eta),Import(*failTree));
  RooCategory probeType("probeType","probeType");
  probeType.defineType("pass");
  probeType.defineType("fail");
  RooAbsData *data;
  // If needed do binned fit
  bool unbinnedFit = false;
  if(unbinnedFit){
    data = new RooDataSet("data","data",mass,Index(probeType),
			  Import("pass",*dataPass), Import("fail",*dataFail));
  }else{
    RooDataHist *dataPassBinned = dataPass->binnedClone("dataPassBinned","dataPassBinned");
    RooDataHist *dataFailBinned = dataFail->binnedClone("dataFailBinned","dataFailBinned");
    data = new RooDataHist("data","data",mass,Index(probeType),
			   Import("pass",*dataPassBinned), Import("fail",*dataFailBinned));    
  }
  // Define the PDFs
  //
  // Common pieces
  //
  // True signal model
  RooRealVar zMass ("zMass" ,"zMass" ,91.188);
  RooRealVar zWidth("zWidth","zWidth",2.495);
  RooBreitWigner bwPdf("bwPdf","bwPdf",mass,zMass,zWidth);
  // Total signal events and efficiency
  RooRealVar nsignal("nsignal","nsignal",1000,0.0,1.0e7);
  RooRealVar eff    ("eff"    ,"eff"    ,0.7,0.0,1.0);

  //  PDF for the PASS sample
  // 
  // Background
  RooRealVar lambdaBgPass("lambdaBgPass","lambdaBgPass",-0.1, -0.5, 0.0);
  RooExponential bgPassPdf("bgPassPdf","bgPassPdf",mass,lambdaBgPass);
  // Signal
  //     - resolution function
  RooRealVar cbMeanPass("cbMeanPass","cbMeanPass"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthPass("cbWidthPass","cbWidthPass",1.0,  0.1, 5.0);
  RooRealVar cbAlphaPass("cbAlphaPass","cbAlphaPass",5.0,  0.0,20.0);
  RooRealVar cbNPass("cbNPass","cbNPass"            ,1.0,  0.0,10.0);
  RooCBShape cbPassPdf("cbPassPdf","cbPassPdf",mass,cbMeanPass,cbWidthPass,cbAlphaPass,cbNPass);
  //     - realistic model
  mass.setBins(10000,setBinsType);
  RooFFTConvPdf signalPassPdf("signalPassPdf","signalPassPdf",mass, bwPdf, cbPassPdf);
  // Combine signal and background
  RooFormulaVar nsigPass("nsigPass","nsigPass","@0*@1",RooArgList(nsignal,eff));
  RooRealVar nbgPass ("nbgPass" ,"nbgPass" ,1,0.0,1.0e5);
  RooAddPdf passPdf("passPdf","passPdf",RooArgList(signalPassPdf,bgPassPdf), RooArgList(nsigPass,nbgPass));
  
  //
  //  PDF for the FAIL sample
  // 
  // Background
  RooRealVar lambdaBgFail("lambdaBgFail","lambdaBgFail",-0.1, -0.5, 0.0);
  RooExponential bgFailPdf("bgFailPdf","bgFailPdf",mass,lambdaBgFail);
  // Signal
  //     - resolution function
  RooRealVar cbMeanFail("cbMeanFail","cbMeanFail"   ,0.0, -10.0,10.0);
  RooRealVar cbWidthFail("cbWidthFail","cbWidthFail",1.0,   0.1, 5.0);
  RooRealVar cbAlphaFail("cbAlphaFail","cbAlphaFail",5.0,   0.0,20.0);
  RooRealVar cbNFail("cbNFail","cbNFail"            ,1.0,   0.0,10.0);
  RooCBShape cbFailPdf("cbFailPdf","cbFailPdf",mass,cbMeanFail,cbWidthFail,cbAlphaFail,cbNFail);
  //     - realistic model
  RooFFTConvPdf signalFailPdf("signalFailPdf","signalFailPdf",mass, bwPdf, cbFailPdf);
  // Combine signal and background
  RooFormulaVar nsigFail("nsigFail","nsigFail","@0*(1.0-@1)",RooArgList(nsignal,eff));
  RooRealVar nbgFail ("nbgFail" ,"nbgFail" ,1,0.0,1.0e5);
  RooAddPdf failPdf("failPdf","failPdf",RooArgList(signalFailPdf,bgFailPdf), RooArgList(nsigFail,nbgFail));
  
  // Combine pass and fail
  RooSimultaneous fullPdf("fullPdf","fullPdf",probeType);
  fullPdf.addPdf(passPdf,"pass");
  fullPdf.addPdf(failPdf,"fail");

  // Do the fit

  // Start with a reasonable point and do rough approximation first
  double total = dataPass->numEntries()+dataFail->numEntries();
  nsignal.setVal( 0.99*total);
  eff.setVal(0.90);
  nbgPass.setVal(0.01*total);
  nbgFail.setVal(0.01*total);
  cbAlphaPass.setVal(1.0);
  cbAlphaFail.setVal(0.5);
  cbNPass    .setVal(5.0);
  cbNFail    .setVal(5.0);
  cbAlphaPass.setConstant(kTRUE);
  cbAlphaFail.setConstant(kTRUE);
  cbNPass    .setConstant(kTRUE);
  cbNFail    .setConstant(kTRUE);
  RooFitResult *result = fullPdf.fitTo(*data,Extended(kTRUE),Save(),RooFit::NumCPU(2,true));

  // Release shape parameters and refine the fit
  cbAlphaPass.setConstant(kFALSE);
  cbAlphaFail.setConstant(kFALSE);
  cbNPass    .setConstant(kFALSE);
  cbNFail    .setConstant(kFALSE);
  result = fullPdf.fitTo(*data,Extended(kTRUE),Save(),RooFit::NumCPU(2,true));

  cout << "Fit status 1st iteration " << result->status() << endl;
//   if(!result->status()){
//     result = fullPdf.fitTo(*data,Extended(kTRUE),Save());
//     cout << "Fit status 2d iteration " << result->status() << endl;
//   }

  // Plot
  passCanvas->cd();
  passCanvas->SetWindowPosition(0,0);
  passCanvas->Draw();
  RooPlot *framePass = mass.frame();
  dataPass->plotOn(framePass);
  passPdf.plotOn(framePass);
  passPdf.plotOn(framePass,Components("bgPassPdf"),LineStyle(kDashed));
  framePass->Draw();

  failCanvas->cd();
  failCanvas->SetWindowPosition(0+ failCanvas->GetWindowWidth(),0);
  failCanvas->Draw();
  RooPlot *frameFail = mass.frame();
  dataFail->plotOn(frameFail);
  failPdf.plotOn(frameFail);
  failPdf.plotOn(frameFail,Components("bgFailPdf"),LineStyle(kDashed));
  frameFail->Draw();

  signal        = nsignal.getVal();
  signalErr     = nsignal.getError();
  efficiency    = eff.getVal();
  efficiencyErr = eff.getError();

  return;
}
*/

void measureEfficiency(TTree *passTree, TTree *failTree, 
		       int method, int etBinning, int etaBinning, 
		       TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
		       bool useTemplates, TFile *templatesFile, 
		       TFile *resultsRootFile, TFile *plotsRootFile,
		       int NsetBins, DYTools::TEfficiencyKind_t effType, 
		       const char* setBinsType, 
		       TString dirTag, const TString &picFileExtraTag, 
		       int puBin) {
  // puBin is important for the fit

  // For COUNTnCOUNT method we should write to root file results
  // from measureEfficiencyCountAndCount routine, otherwise
  // from measureEfficiencyWithFit routine.
  bool saveCountingToRootFile = true;
  if( method == DYTools::COUNTnFIT || method == DYTools::FITnFIT )
    saveCountingToRootFile = false;
  
  // Always report counting method results
  measureEfficiencyCountAndCount(passTree, failTree, etBinning, etaBinning, 
				 canvas, effOutput, 
				 saveCountingToRootFile, resultsRootFile, 
				 plotsRootFile, effType);
  
  if( method == DYTools::COUNTnFIT || method == DYTools::FITnFIT ) {
    measureEfficiencyWithFit(passTree, failTree, 
			     method, etBinning, etaBinning, 
			     canvas, effOutput, fitLog,
			     useTemplates, templatesFile, 
			     resultsRootFile, plotsRootFile,
			     NsetBins, effType, setBinsType, dirTag, 
			     picFileExtraTag, puBin);
  }
  

  return;
}

// -------------------------------------------------------------------------

void measureEfficiencyPU(TTree *passTreeFull, TTree *failTreeFull, 
			 int method, int etBinning, int etaBinning, 
			 TCanvas *canvas,ofstream &effOutput, ofstream &fitLog,
			 bool useTemplates, TFile *templatesFile, 
			 const TString &resultRootFileBase,
			 int NsetBins, DYTools::TEfficiencyKind_t effType, 
			 const char* setBinsType, 
			 TString dirTag, const TString &picFileExtraTag,
			 int puDependence
			 ) {

  if (!puDependence) {
    fitLog << "\nmeasureEfficiencyPU: no PU-dependence requested\n";
    TString resRootFName=resultRootFileBase + TString(".root");
    TFile *resultsRootFile=new TFile(resRootFName,"recreate");
    TString resPlotsFName=resultRootFileBase + TString("-plots.root");
    TFile *resultPlotsFile=new TFile(resPlotsFName,"recreate");
    measureEfficiency(passTreeFull,failTreeFull,method,etBinning,etaBinning,
		      canvas,effOutput,fitLog,useTemplates,
		      templatesFile,resultsRootFile,resultPlotsFile,
		      NsetBins,effType,setBinsType,
		      dirTag,picFileExtraTag);
  }
  else {
    // prepare pu-dependent trees
    vector<TTree*> passTreeV, failTreeV;
    passTreeV.reserve(DYTools::nPVBinCount);
    failTreeV.reserve(DYTools::nPVBinCount);
    Double_t mass,et,eta;
    UInt_t nPV;

    // 1. create trees
    for (int pass=0; pass<2; ++pass) {
      for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
	char buf[50];
	UInt_t pvMin=UInt_t(DYTools::nPVLimits[pu_i  ]+0.6);
	UInt_t pvMax=UInt_t(DYTools::nPVLimits[pu_i+1]-0.4);
	sprintf(buf,"%sTreePU_%u_%u",(pass)? "pass":"fail",pvMin,pvMax);
	//std::cout << "buf=" << buf << "\n";
	TTree *tree=new TTree(buf,buf);
	tree->SetDirectory(0);
	tree->Branch("mass",&mass,"mass/D");
	tree->Branch("et",&et,"et/D");
	tree->Branch("eta",&eta,"eta/D");
	if (pass) passTreeV.push_back(tree); else failTreeV.push_back(tree);
      }
    }

    // fill trees
    UInt_t excludedPass=0, excludedFail=0;
    for (int pass=0; pass<2; ++pass) {
      vector<TTree*> *treeV = (pass) ? &passTreeV : &failTreeV;
      TTree *treeFull=(pass) ? passTreeFull : failTreeFull;
      treeFull->SetBranchAddress("mass",&mass);
      treeFull->SetBranchAddress("et",&et);
      treeFull->SetBranchAddress("eta",&eta);
      treeFull->SetBranchAddress("nGoodPV",&nPV);
      for (UInt_t i=0; i<treeFull->GetEntries(); ++i) {
	treeFull->GetEntry(i);
	unsigned int pu_idx=(unsigned int)(DYTools::findPUBin(nPV));
	if (pu_idx<treeV->size()) {
	  (*treeV)[pu_idx]->Fill();
	}
	else {
	  if (pass) excludedPass++; else excludedFail++;
	}
      }
    }
    if (1) {
      effOutput << "\nEvent distribution by PU\n";
      effOutput << "  PU range    passCount    failCount\n";
      char buf[100];
      UInt_t nPass=excludedPass;
      UInt_t nFail=excludedFail;
      for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
	nPass+=passTreeV[pu_i]->GetEntries();
	nFail+=failTreeV[pu_i]->GetEntries();
	sprintf(buf,"  %4.1lf..%4.1lf    %8lld   %8lld", 
		DYTools::nPVLimits[pu_i],DYTools::nPVLimits[pu_i+1],
		passTreeV[pu_i]->GetEntries(), failTreeV[pu_i]->GetEntries());
	effOutput << buf << "\n";
      }
      sprintf(buf,"  excluded      %8u   %8u", excludedPass,excludedFail);
      effOutput << buf << "\n"
		<<"-------------------------------------\n";
      sprintf(buf,"  total         %8u   %8u", nPass,nFail);
      effOutput << buf << "\n";
    }
  
    // measure efficiency
    for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
      std::cout << " pu_i=" << pu_i << "\n";
      char buf[50];
      UInt_t pvMin=UInt_t(DYTools::nPVLimits[pu_i  ]+0.6);
      UInt_t pvMax=UInt_t(DYTools::nPVLimits[pu_i+1]-0.4);
      sprintf(buf,"_%u_%u.root",pvMin,pvMax);
      TString resRootFName=resultRootFileBase + TString(buf);
      TFile *resultsRootFile=new TFile(resRootFName,"recreate");
      sprintf(buf,"-plots_%u_%u.root",pvMin,pvMax);
      TString resPlotFName=resultRootFileBase + TString(buf);
      TFile *resultPlotFile=new TFile(resPlotFName,"recreate");
      fitLog << "\nmeasureEfficiencyPU: PU range " 
	     << pvMin << " - " << pvMax << "\n";
      effOutput << "\nmeasureEfficiencyPU: PU range " 
		<< pvMin << " - " << pvMax << "\n";
      std::cout << "call measure efficiency" << std::endl;
      measureEfficiency(passTreeV[pu_i],failTreeV[pu_i],
			method,etBinning,etaBinning,
			canvas,effOutput,fitLog,useTemplates,templatesFile,
			resultsRootFile,resultPlotFile,
			NsetBins,effType,setBinsType,
			dirTag,picFileExtraTag, pu_i+1);
      std::cout << "done measure efficiency" << std::endl;
    }

    ClearVec(passTreeV); ClearVec(failTreeV);
  }

  return;
}

void measureEfficiencyCountAndCount(TTree *passTree, TTree *failTree, 
			    int etBinning, int etaBinning, 
			    TCanvas *canvas, ofstream &effOutput,
			    bool saveResultsToRootFile, TFile *resultsRootFile,
				    TFile *resultPlotsFile, 
				    DYTools::TEfficiencyKind_t effType){
  
  int nEt                = DYTools::getNEtBins(etBinning);
  const double *limitsEt = DYTools::getEtBinLimits(etBinning);
  
  int nEta                = DYTools::getNEtaBins(etaBinning);
  const double *limitsEta = DYTools::getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);

  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);
  TMatrixD effArray2DWeighted(nEt, nEta);
  TMatrixD effArrayErrLow2DWeighted(nEt, nEta);
  TMatrixD effArrayErrHigh2DWeighted(nEt, nEta);
 
  std::vector<std::string> lines;
  lines.reserve(50);

  effOutput << endl;
  effOutput << "Efficiency, counting method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";

  lines.push_back("\nEfficiency, counting method (weighted):\n");
  lines.push_back("     SC ET         SC eta           efficiency             pass         fail\n");
  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
      double effCount, effErrLowCount, effErrHighCount;
      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCutFormat= DYTools::signedEtaBinning(etaBinning) ?
	" ( eta >= %5.3f && eta < %5.3f ) " :
	" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ";

      double limitsEtaMin = limitsEta[j];
      double limitsEtaMax = limitsEta[j+1];	  
      // For RECO efficiency, to increase fit stability, MERGE barrel
      // eta bins, and separately endcap eta bins, for the binning ETABINS5
      // and Et < 20 GeV
      //
      // While this is count and count method with no fit involved, the change
      // is made here as well to match the corresponding change in the
      // measureEfficiencyWithFit(...) function, so that the averaging over
      // the eta bins is the same.
      //
      // In the way it is implemented now, there is duplication of work.
      // For example, when the eta bins 0.0-0.8 and 0.8-1.4442 are merged,
      // we do the full fitting of the data in 0.0-1.4442 twice,
      // once for j=0, and then again for j=1. The result of this will be identical
      // and some CPU will be wasted. However, the bookkeeping this way is much
      // easier, and the code is easier to read. This can be changed in the 
      // future if needed.
      //   Since the binning is not changed, the entries in the efficiency array
      // will be exactly the same for the merged bins. 
      //   The calculation of the error on the event-level scale factors as 
      // a function of mass takes into account this 100% correlation in
      // the efficiencies of the merged bins. This is done in calcEventEff.C
      bool isRECO=(effType == DYTools::RECO) ? true : false;
      if( isRECO && etaBinning == DYTools::ETABINS5 && limitsEt[i+1]  <= 20.0){
	if ( j == 0 || j == 1 ){
	  // Barrel eta bins
	  limitsEtaMin = limitsEta[0];
	  limitsEtaMax = limitsEta[2];
	  printf("MERGE two barrel eta bins j=0 j=1 into one. The efficiency measurement\n");
	  printf("      is done twice with identical data, this is the instance j=%d\n", j);
	} else if ( j ==3 || j ==4 ) {
	  // Endcap eta bins
	  limitsEtaMin = limitsEta[3];
	  limitsEtaMax = limitsEta[5];
	  printf("MERGE two endcap eta bins j=3 j=4 into one. The efficiency measurement\n");
	  printf("      is done twice with identical data, this is the instance j=%d\n", j);
	} else {
	  // Anything else (really, the rapidity gap, j==2)
	  // No need to set anything, set above in declaration/initialization
	}
      }

      TString etaCut = TString::Format( etaCutFormat,
					limitsEtaMin, limitsEtaMax);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);

      // padIdx = 1 + 2*(i + j*nEt)
      int padIdx= 1 + ( (2*i + 2*j*nEt) % (2*DYTools::maxTnPCanvasDivisions) );
      canvas->cd(padIdx);
      passTree->Draw("weight >> hwPass",cut);
      canvas->cd(padIdx+1);
      failTree->Draw("weight >> hwFail",cut);
      TH1 *hwPass=(TH1*)gDirectory->Get("hwPass");
      TH1 *hwFail=(TH1*)gDirectory->Get("hwFail");
      double probesPassWeighted = hwPass->GetMean() * hwPass->GetEntries();
      double probesFailWeighted = hwFail->GetMean() * hwFail->GetEntries();
      double effCountWeighted, effErrLowCountWeighted, effErrHighCountWeighted;
      //std::cout << " probesPass=" << probesPass << ", probesFail=" << probesFail << "\n";
      //std::cout << " hwPass->GetEntries=" << hwPass->GetEntries() << ", hwFail->GetEntries=" << hwFail->GetEntries() << "\n";
      //std::cout << " probesPassWeighted=" << probesPassWeighted << ", probesFailWeighted=" << probesFailWeighted << "\n";

      DYTools::calcEfficiency( probesPass, probesPass+probesFail, 
			       DYTools::EFF_CLOPPER_PEARSON,
			       effCount, effErrLowCount, effErrHighCount);
      DYTools::calcEfficiency( probesPassWeighted, 
			       probesPassWeighted+probesFailWeighted,
			       DYTools::EFF_CLOPPER_PEARSON,
			       effCountWeighted, 
			       effErrLowCountWeighted,effErrHighCountWeighted);
      const int len=200;
      char strOut[len+1];
      snprintf(strOut,len, "hweightPass_Et_%1.0f-%1.0f__Eta_%5.3f-%5.3f",
	      limitsEt[i],limitsEt[i+1],
	      limitsEta[j],limitsEta[j+1]);
      hwPass->SetName(strOut);
      if (resultPlotsFile) { resultPlotsFile->cd(); hwPass->Write(); }
      snprintf(strOut,len, "hweightFail_Et_%1.0f-%1.0f__Eta_%5.3f-%5.3f",
	      limitsEt[i],limitsEt[i+1],
	      limitsEta[j],limitsEta[j+1]);
      hwFail->SetName(strOut);
      if (resultPlotsFile) { resultPlotsFile->cd(); hwFail->Write(); }

      snprintf(strOut,len, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1],
	     effCount*100, effErrHighCount*100, effErrLowCount*100,
	     probesPass, probesFail);
      effOutput << strOut;

      snprintf(strOut,len, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      effCountWeighted*100, 
	      effErrHighCountWeighted*100, effErrLowCountWeighted*100,
	      probesPassWeighted, probesFailWeighted);
      lines.push_back(strOut);

      // padIdx = 1 + 2*(i + j*nEt)
      padIdx= 1 + ( (2*i + 2*j*nEt) % (2*DYTools::maxTnPCanvasDivisions) );
      canvas->cd(padIdx);
      passTree->Draw("mass",cut);
      canvas->cd(padIdx+1);
      failTree->Draw("mass",cut);
      canvas->Update();
      if (resultPlotsFile) {
	resultPlotsFile->cd();
	TString canvName=TString::Format("canvMassDistr_Et_%2.0f-%2.0f_abs_eta_%5.3f-%5.3f",
					 limitsEt[i],limitsEt[i+1],
					 limitsEta[j],limitsEta[j+1]);
	canvName.ReplaceAll(".","_");
	TCanvas ctemp(canvName,canvName,900,400);
	ctemp.Divide(2,1);
	ctemp.cd(1); passTree->Draw("mass",cut);
	ctemp.cd(2); failTree->Draw("mass",cut);
	ctemp.Write();
      }

      // Add systematics for HLT efficiency (see comments at the top
      // when the arrays are introduced
      bool isHLT = (effType == DYTools::HLT) ? true : false;
      if( isHLT) {
	if(etBinning == DYTools::ETBINS6 && etaBinning == DYTools::ETABINS5){
	  double extraErr = effHltSystErrBarrel[i];
	  if( j > 2 ) 
	    extraErr = effHltSystErrEndcap[i];
	  double extraErrSquared = extraErr*extraErr;
	  // Update the errors
	  effErrLowCount = sqrt(effErrLowCount*effErrLowCount + extraErrSquared);
	  effErrHighCount = sqrt(effErrHighCount*effErrHighCount + extraErrSquared);
	  effErrLowCountWeighted = sqrt(effErrLowCountWeighted*effErrLowCountWeighted + extraErrSquared);
	  effErrHighCountWeighted = sqrt(effErrHighCountWeighted*effErrHighCountWeighted + extraErrSquared);
	}
      }

      effArray2D(i,j) = effCount;
      effArrayErrLow2D(i,j) = effErrLowCount;
      effArrayErrHigh2D(i,j) = effErrHighCount;
      effArray2DWeighted(i,j) = effCountWeighted;
      effArrayErrLow2DWeighted(i,j) = effErrLowCountWeighted;
      effArrayErrHigh2DWeighted(i,j) = effErrHighCountWeighted;
    }
  }
  effOutput << endl;
  for (unsigned int i=0; i<lines.size(); ++i) {
    effOutput << lines[i];
  }
  effOutput << endl;

  if(saveResultsToRootFile){
    if(resultsRootFile && resultsRootFile->IsOpen()){
      resultsRootFile->cd();
      effArray2D.Write("effArray2D");
      effArrayErrLow2D.Write("effArrayErrLow2D");
      effArrayErrHigh2D.Write("effArrayErrHigh2D");    
      effArray2DWeighted.Write("effArray2DWeighted");
      effArrayErrLow2DWeighted.Write("effArrayErrLow2DWeighted");
      effArrayErrHigh2DWeighted.Write("effArrayErrHigh2DWeighted");    
      resultsRootFile->Close();
    }else assert(0);
  }

  if (resultPlotsFile && resultPlotsFile->IsOpen()) {
    resultPlotsFile->cd();
    canvas->Write();
  }

  return;
}

// --------------------------------------------------

void measureEfficiencyWithFit(TTree *passTree, TTree *failTree, 
			      int method, int etBinning, int etaBinning, 
			      TCanvas *canvas, ofstream &effOutput, ofstream &fitLog,
			      bool useTemplates, TFile *templatesFile, 
			      TFile *resultsRootFile, TFile *resultPlotsFile,
			      int NsetBins, DYTools::TEfficiencyKind_t effType,
			      const char* setBinsType, 
			      TString dirTag, const TString &picFileExtraTag,
			      int puBin){
  
  int nEt                = DYTools::getNEtBins(etBinning);
  const double *limitsEt = DYTools::getEtBinLimits(etBinning);

  int nEta                = DYTools::getNEtaBins(etaBinning);
  const double *limitsEta = DYTools::getEtaBinLimits(etaBinning);
  printf("eta bins %d\n", nEta);
  
  bool isRECO=(effType == DYTools::RECO) ? true : false;
      
  TMatrixD effArray2D(nEt, nEta);
  TMatrixD effArrayErrLow2D(nEt, nEta);
  TMatrixD effArrayErrHigh2D(nEt, nEta);

  TString methodStr;
  switch(method) {
  case DYTools::COUNTnCOUNT: methodStr="Count+Count"; break;
  case DYTools::COUNTnFIT: methodStr="Count+Fit"; break;
  case DYTools::FITnFIT: methodStr="Fit+Fit"; break;
  default:
    std::cout << "unknown method of efficiency calculation\n";
    assert(0);
  }
 
  effOutput << endl;
  effOutput << "Efficiency, " << methodStr << " method:\n";  
  effOutput << "     SC ET         SC eta           efficiency             pass         fail\n";
  if (targetSpecificBin) {
    TString msg= Form("\n\ttargetSpecificBin is on. targetEt=%d, targetEta=%d\n",targetEt,targetEta);
    std::cout << msg;
    effOutput << msg;
  }

  for(int j=0; j<nEta; j++){
    for(int i=0; i<nEt; i++){
  //for(int j=0; j<1; j++){
    //for(int i=0; i<1; i++){
      if (targetSpecificBin==1)  {
	if ((targetEt>=0) && (targetEt!=i)) continue;
	if ((targetEta>=0) && (targetEta!=j)) continue;
      }

      TString etCut = TString::Format(" ( et >=%6.1f && et <%6.1f ) ",
				      limitsEt[i], limitsEt[i+1]);
      TString etaCutFormat= DYTools::signedEtaBinning(etaBinning) ?
	" ( eta >= %5.3f && eta < %5.3f ) " :
	" ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ";

      double limitsEtaMin = limitsEta[j];
      double limitsEtaMax = limitsEta[j+1];	  
      // For RECO efficiency, to increase fit stability, MERGE barrel
      // eta bins, and separately endcap eta bins, for the binning ETABINS5,
      // and Et < 20. 
      //
      // In the way it is implemented now, there is duplication of work.
      // For example, when the eta bins 0.0-0.8 and 0.8-1.4442 are merged,
      // we do the full fitting of the data in 0.0-1.4442 twice,
      // once for j=0, and then again for j=1. The result of this will be identical
      // and some CPU will be wasted. However, the bookkeeping this way is much
      // easier, and the code is easier to read. This can be changed in the 
      // future if needed.
      //   Since the binning is not changed, the entries in the efficiency array
      // will be exactly the same for the merged bins. 
      //   The calculation of the error on the event-level scale factors as 
      // a function of mass takes into account this 100% correlation in
      // the efficiencies of the merged bins. This is done in calcEventEff.C
      if( isRECO && etaBinning == DYTools::ETABINS5 && limitsEt[i+1] <= 20.0){
	if ( j == 0 || j == 1 ){
	  // Barrel eta bins
	  limitsEtaMin = limitsEta[0];
	  limitsEtaMax = limitsEta[2];
	  printf("MERGE two barrel eta bins j=0 j=1 into one. The efficiency measurement\n");
	  printf("      is done twice with identical data, this is the instance j=%d for Et bin i=%d\n", j,i);
	} else if ( j ==3 || j ==4 ) {
	  // Endcap eta bins
	  limitsEtaMin = limitsEta[3];
	  limitsEtaMax = limitsEta[5];
	  printf("MERGE two endcap eta bins j=3 j=4 into one. The efficiency measurement\n");
	  printf("      is done twice with identical data, this is the instance j=%d for Et bin i=%d\n", j,i);
	} else {
	  // Anything else (really, the rapidity gap, j==2)
	  // No need to set anything, set above in declaration/initialization
	}
      }

      TString etaCut = TString::Format( etaCutFormat,
					limitsEtaMin, limitsEtaMax);
      TString cut = etCut + TString(" && ") + etaCut;
      //cout << cut.Data() << endl;
      double probesPass = passTree->GetEntries(cut);
      double probesFail = failTree->GetEntries(cut);
      // padIdx = 1 + 2*(i + j*nEt)
      int padIdx= 1 + ( (2*i + 2*j*nEt) % (2*DYTools::maxTnPCanvasDivisions) );
      TPad *passPad = (TPad*)canvas->GetPad(padIdx);
      TPad *failPad = (TPad*)canvas->GetPad(padIdx + 1);
      double efficiency, efficiencyErrHi, efficiencyErrLo;
      printf("\n ==\n");
      char strOut[200];
      sprintf(strOut,
	   " ==   Start fitting Et: %3.0f - %3.0f  and eta:  %5.3f - %5.3f \n",
	     limitsEt[i], limitsEt[i+1],
	     limitsEta[j], limitsEta[j+1]);
      printf("%s",strOut);
      printf(" ==\n\n");
      fitLog << endl << strOut << endl;
      
      if(!useTemplates){
	fitMass(passTree, failTree, cut, method, 
		efficiency,efficiencyErrHi, efficiencyErrLo, passPad, failPad, 
		resultPlotsFile,
		fitLog, NsetBins, isRECO, setBinsType, dirTag);
      }
      else{
	printf("\nMASS TEMPLATES ARE USED IN THE FIT\n\n");
	// In case templates are used, find the right templates
	TH1F *templatePass = 
	  getPassTemplate(i,j,etaBinning, templatesFile, puBin);
	TH1F *templateFail = 
	  getFailTemplate(i,j,etaBinning, templatesFile, puBin);

	// In case if MERGE of the eta bins is needed for RECO efficiency (see above
	// more detailed comments about this, we add the templates of the appropriate
	// bins.
	// Note: the original templates are cloned, so that we do not mess with 
	// the original ones. This is important because each template is used twice
	// in this implementation, once for each of the bins being merged.
	if( isRECO && etaBinning == DYTools::ETABINS5 && limitsEt[i+1] <= 20.0 ){
	  if ( j == 0 || j == 1 ){
	    // Barrel eta bins
	    TH1F *templatePassOne = (TH1F*)getPassTemplate(i,0,etaBinning, templatesFile, puBin)->Clone();
	    TH1F *templatePassTwo = (TH1F*)getPassTemplate(i,1,etaBinning, templatesFile, puBin)->Clone();
	    templatePassOne->Add(templatePassTwo);
	    templatePass = templatePassOne;
	    TH1F *templateFailOne = (TH1F*)getFailTemplate(i,0,etaBinning, templatesFile, puBin)->Clone();
	    TH1F *templateFailTwo = (TH1F*)getFailTemplate(i,1,etaBinning, templatesFile, puBin)->Clone();
	    templateFailOne->Add(templateFailTwo);
	    templateFail = templateFailOne;
	    printf("MERGE templates for two barrel eta bins j=0 j=1 into one. The efficiency measurement\n");
	    printf("      is done twice with identical data and templates, this is the instance j=%d for Et bin i=%d\n", j,i);
	  } else if ( j ==3 || j ==4 ) {
	    // Endcap eta bins
	    TH1F *templatePassOne = (TH1F*)getPassTemplate(i,3,etaBinning, templatesFile, puBin)->Clone();
	    TH1F *templatePassTwo = (TH1F*)getPassTemplate(i,4,etaBinning, templatesFile, puBin)->Clone();
	    templatePassOne->Add(templatePassTwo);
	    templatePass = templatePassOne;
	    TH1F *templateFailOne = (TH1F*)getFailTemplate(i,3,etaBinning, templatesFile, puBin)->Clone();
	    TH1F *templateFailTwo = (TH1F*)getFailTemplate(i,4,etaBinning, templatesFile, puBin)->Clone();
	    templateFailOne->Add(templateFailTwo);
	    templateFail = templateFailOne;
	    limitsEtaMin = limitsEta[3];
	    limitsEtaMax = limitsEta[5];
	    printf("MERGE templates for two endcap eta bins j=3 j=4 into one. The efficiency measurement\n");
	    printf("      is done twice with identical data and templates, this is the instance j=%d for Et bin i=%d\n", j,i);
	  } else {
	    // Anything else (really, the rapidity gap, j==2)
	    // No need to set anything, set above in declaration/initialization
	  }
	}

	if (0) {
	  TCanvas *cx=new TCanvas("cx","cx",600,600);
	  templatePass->SetLineColor(kGreen);
	  templatePass->DrawCopy("hist");
	  templateFail->SetLineColor(kRed);
	  templateFail->DrawCopy("hist same");
	  cx->Update();
	  double xx; std::cin>>xx;
	}
      
	fitMassWithTemplates(passTree, failTree, cut, method, 
			     efficiency, efficiencyErrHi, efficiencyErrLo,
			     passPad, failPad, resultPlotsFile,
			     fitLog, templatePass, templateFail, 
			     isRECO, setBinsType, dirTag, picFileExtraTag);
	//std::cout << "after fitMassWithTemplates: enter a value\n"; double x; std::cin >> x;
	char buf[50];
	sprintf(buf,"tmp_eta%d_et%d.png",j,i);
	canvas->SaveAs(buf);
 
      }
            

      sprintf(strOut, "   %3.0f - %3.0f   %5.3f - %5.3f   %5.1f +%5.1f -%5.1f    %10.0f  %10.0f\n",
	      limitsEt[i], limitsEt[i+1],
	      limitsEta[j], limitsEta[j+1],
	      efficiency*100, efficiencyErrHi*100, efficiencyErrLo*100,
	      probesPass, probesFail);
      effOutput << strOut;
      effArray2D(i,j) = efficiency;
      // Add systematics for RECO efficiency found from a fit
      if( isRECO ){
	if( etBinning == DYTools::ETBINS6 ){
	  efficiencyErrLo = sqrt(efficiencyErrLo*efficiencyErrLo
				 + effRecoSystErr[i]*effRecoSystErr[i]);
	  efficiencyErrHi = sqrt(efficiencyErrHi*efficiencyErrHi
				 + effRecoSystErr[i]*effRecoSystErr[i]);
	  printf("RECO efficiency for (Et,eta) index (%d,%d): %.0f%% of systematic error added\n", i,j, effRecoSystErr[i]*100);
	}else{
	  printf("RECO efficiency for (Et,eta) index (%d,%d): no systematic is added\n",i,j );
	}
      }
      
      effArrayErrLow2D(i,j) = efficiencyErrLo;
      effArrayErrHigh2D(i,j) = efficiencyErrHi;
      
    }
  }

  effOutput << endl;

  if(resultsRootFile && resultsRootFile->IsOpen()){
    resultsRootFile->cd();
    effArray2D.Write("effArray2D");
    effArrayErrLow2D.Write("effArrayErrLow2D");
    effArrayErrHigh2D.Write("effArrayErrHigh2D");    
    resultsRootFile->Close();
  }else assert(0);

  if (resultPlotsFile && resultPlotsFile->IsOpen()) {
    resultPlotsFile->cd();
    canvas->Write();
  }

  //std::cout << "before leaving measureEfficiencyWithFit: enter a value\n"; double x; std::cin >> x;
   
  return;
}

// --------------------------------------------------

int getTemplateBin(int etBin, int etaBin, int etaBinning){

  int templateBin = -1;

  if( etBin != -1 && etaBin != -1)
    templateBin = etBin * DYTools::getNEtaBins(etaBinning) + etaBin;

  return templateBin;

}

// --------------------------------------------------

TString getTemplateName(int etBin, int etaBin, const char *pass_fail_str, 
			int puBin) {
  char buf[50];
  sprintf(buf,"hMassTemplate_Et%d_eta%d",etBin,etaBin);
  int len=0;
  if (puBin>0) {
    len=strlen(buf);
    int puMin=int(DYTools::nPVLimits[puBin-1]+0.6);
    int puMax=int(DYTools::nPVLimits[puBin  ]-0.4);
    sprintf(buf+len,"_%d_%d",puMin,puMax);
  }
  len=strlen(buf);
  sprintf(buf+len,"_%s",pass_fail_str);
  return TString(buf);
}

// --------------------------------------------------

TH1F * getPassTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = getTemplateName(etBin,etaBin,"pass",puBin);

  hist = (TH1F*)file->Get(name);
  return hist;
}

// --------------------------------------------------

TH1F * getFailTemplate(int etBin, int etaBin, int etaBinning, TFile *file, 
		       int puBin){

  TH1F *hist = 0;
  if(file == 0)
    return hist;

  int templateBin = getTemplateBin(etBin, etaBin, etaBinning);
  if( templateBin == -1 )
    return hist;

  TString name = getTemplateName(etBin,etaBin,"fail",puBin);

  hist = (TH1F*)file->Get(name);
  return hist;
}
