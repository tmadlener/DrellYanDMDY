#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

// The global variables will be used in several functions:
TString tagDirYields = "";
TString tagDirConstants = "";
Double_t luminosity = 0;

void createDummyAcceptanceSystematics(const TString conf){

  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "calcAcceptanceSystematics: _DebugRun_ detected. Terminating the script\n";
    return;
  }

  // normal calculation

  // First, read the configuration file. The configuration file
  // is the same as the one used for the cross section calculation
  // script, we need to know the location of data yields and
  // unfolding matrices

  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  std::string line;
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

  /////////////////////////////////
  //calculate Fsr systematics 
  /////////////////////////////////

  TMatrixD accFsrMax(DYTools::nMassBins,DYTools::npTBinsMax);
  TMatrixD accFsrMin(DYTools::nMassBins,DYTools::npTBinsMax);
  TMatrixD accSystPercentFsr(DYTools::nMassBins,DYTools::npTBinsMax);
  accFsrMax = 0;
  accFsrMin = 0;
  accSystPercentFsr = 0;

  // Store constants in the file
  TString fileEnd(DYTools::analysisTag + TString("_dummy.root"));
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString acceptanceSystFileName(outputDir+TString("/acceptance_FSR_systematics") + fileEnd);

  TFile fa(acceptanceSystFileName,"recreate");
  accSystPercentFsr.Write("accSystPercent");
  accFsrMin.Write("accFSRSystMin");
  accFsrMax.Write("accFsrSystMax");
  unfolding::writeBinningArrays(fa);
  fa.Close();

  std::cout << "Write DUMMY acceptance systematics file filled with zeros, file name " << acceptanceSystFileName << std::endl;

  return;
}

