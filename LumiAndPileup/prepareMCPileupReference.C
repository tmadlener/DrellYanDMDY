#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

bool useMeanPileup = true;

void prepareMCPileupReference(){

//   TFile f("/home/hep/ikrav/work/ntuples/DrellYan_8TeV_53X_local/s12-zeem20-v7a_tight-loose_skim.root");
  TFile f("/home/hep/ikrav/work/ntuples/DrellYan_8TeV_53X_with_regression/s12-zllm50-v7a_NoReg_ntuple.root");

  TTree *tree = (TTree*)f.Get("Events");

  TH1F *pileup = new TH1F("pileup","",100,0,100);

  if( useMeanPileup ){
    tree->Draw("Info.nPUmean>>pileup");
    printf("Use the MEAN pileup field\n");
  }else{
    tree->Draw("Info.nPU>>pileup");
    printf("Use the GENERATED pileup field\n");
  }

  pileup->SetDirectory(0);
  pileup->Print();

  TFile fout("mcPileupHildreth_full2012.root", "recreate");
  fout.cd();
  pileup->Write();
  fout.Close();
}
