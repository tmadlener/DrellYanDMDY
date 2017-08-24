#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

using std::string;
using std::stringstream;

TString tagDirYields = "";
TString tagDirConstants = "";
Double_t luminosity = 0;

//
//  Main code
//

void createDummyUnfoldingSystematics(const TString conf){

  // First, read the configuration file. The configuration file
  // is the same as the one used for the cross section calculation
  // script, we need to know the location of data yields and
  // unfolding matrices

  ifstream ifs;
  ifs.open(conf.Data());
  if(!ifs.is_open())
    assert(0);
  string line;
  int state = 0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state==0){
      stringstream ss1(line); ss1 >> luminosity;
      state++;
    }else if(state==1){
      tagDirYields = TString(line);
      state++;
    }else if(state==2){
      tagDirConstants = TString(line);
      break;
    }
  }
  ifs.close();
  
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();

  TVectorD unfoldedYieldsMean(nUnfoldingBins);
  TVectorD unfoldedYieldsRMS(nUnfoldingBins);
  TVectorD unfoldedYieldsSquaredMean(nUnfoldingBins);
  
  unfoldedYieldsMean = 1e-6; 		//C: Just to see if something changes
  unfoldedYieldsRMS = 1e-6; 		//C: Just to see if something changes
  unfoldedYieldsSquaredMean = 1e-6; 	//C: Just to see if something changes

  TVectorD unfoldingSystPercentSmear(nUnfoldingBins);
  TVectorD unfoldingSystPercentFsr(nUnfoldingBins); 
  TVectorD unfoldingSystPercent(nUnfoldingBins); 

  unfoldingSystPercentSmear = 1e-6; 		//C: Just to see if something changes
  unfoldingSystPercentFsr = 1e-6; 		//C: Just to see if something changes
  unfoldingSystPercent = 1e-6; 			//C: Just to see if something changes

  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/unfolding_systematics") 
				+ DYTools::analysisTag + TString("_dummy.root"));

  TFile fa(unfoldingSystFileName,"recreate");
  unfolding::writeBinningArrays(fa);
  unfoldedYieldsMean.Write("unfoldedYieldsMeanFI");
  unfoldedYieldsRMS.Write("unfoldedYieldsRMSFI");
  unfoldingSystPercentSmear.Write("unfoldingSystPercentSmearFI");
  unfoldingSystPercentFsr.Write("unfoldingSystPercentFsrFI");
  unfoldingSystPercent.Write("unfoldingSystPercentFI");
  fa.Close();

  std::cout << "Write DUMMY unfolding systematics, filled with zeros, into the file " << unfoldingSystFileName << std::endl;

  return;
}

