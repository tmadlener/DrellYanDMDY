#include "TVectorD.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingTools.hh"

// -----------------------------------------------------------

void createDummyEscaleSystematics(TString lumiTag="DY_j22_19712pb", int saveTexTable=0){

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  
  TVectorD escaleRandomizedSystRelative(nUnfoldingBins);
  TVectorD escaleResidualDiffSystRelative(nUnfoldingBins);
  TVectorD escaleFitShapeSystRelative(nUnfoldingBins);
  TVectorD escaleEtaBinSystRelative(nUnfoldingBins);
  TVectorD escaleSystPercent(nUnfoldingBins);
  escaleRandomizedSystRelative = 1e-6; 		//C: Just to see if something changes
  escaleResidualDiffSystRelative=1e-6; 		//C: Just to see if something changes
  escaleFitShapeSystRelative = 1e-6; 		//C: Just to see if something changes
  escaleEtaBinSystRelative = 1e-6; 		//C: Just to see if something changes
  escaleSystPercent = 1e-6; 			//C: Just to see if something changes

  TString finalFName=TString("../root_files/systematics/") + lumiTag + 
    TString("/escale_systematics") + DYTools::analysisTag + TString("_dummy.root");
  TFile fout(finalFName,"recreate");
  unfolding::writeBinningArrays(fout);
  escaleRandomizedSystRelative.Write("escaleRandomizedSystRelativeFI");
  escaleFitShapeSystRelative.Write("escaleFitShapeSystRelativeFI");
  escaleResidualDiffSystRelative.Write("escaleResidualDiffSystRelativeFI");
  escaleEtaBinSystRelative.Write("escaleEtaBinSystRelativeFI");
  escaleSystPercent.Write("escaleSystPercentFI");

  fout.Close();

  std::cout << "Created DUMMY escale systematics filled with zeros, file name " << finalFName << std::endl;

  return;
}


