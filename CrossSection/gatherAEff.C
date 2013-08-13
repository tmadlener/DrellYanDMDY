#include <TVectorD.h>
#include <TMatrixD.h>

#include <fstream>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingTools.hh"

typedef TH1F Histo_t; // match the type returned in MyTools.hh

// -----------------------------------------

void AdjustDim(std::string &line) {
  const std::string expect =( DYTools::study2D) ? "2D" : "1D";
  const std::string canHave=(!DYTools::study2D) ? "2D" : "1D";
  if (!PosOk(line,expect)) {
    size_t pos=line.find(canHave);
    if (PosOk(pos)) {
      std::cout << "changing <" << line << ">\n";
      line[pos]=expect[0];
      std::cout << "      to <" << line << ">\n";
    }
  }
  return;
}

// -----------------------------------------

int LoadMatrix(const TString &fname, TMatrixD &M, TMatrixD &Merr, const TString &field, const TString &fieldErr) {
  TFile f(fname,"read");
  TMatrixD *Mptr=(TMatrixD*)f.Get(field);
  TMatrixD *MerrPtr=(TMatrixD*)f.Get(fieldErr);
  f.Close();
  if (!Mptr || !MerrPtr) {
    std::cout << "error in loading from <" << fname << ">\n";
    if (!Mptr) std::cout << " failed to load <" << field << ">\n";
    if (!MerrPtr) std::cout << " failed to load <" << fieldErr << ">\n";
    return 0;
  }
  if ((Mptr->GetNrows() != M.GetNrows()) ||
      (Mptr->GetNcols() != M.GetNcols()) ) {
    std::cout << "dim mismatch in <" << fname << ">\n";
    return 0;
  }
  M = *Mptr;
  Merr= *MerrPtr;
  delete Mptr;
  delete MerrPtr;
  return 1;
}

// -----------------------------------------

int GetMassProfile(const TMatrixD &m, int yIdx, TVectorD &v) {
  v=0;
  for (int i=0; i<m.GetNrows(); ++i) {
    v[i] = m(i,yIdx);
  }
  return 1;
}

// -----------------------------------------

int GetMassProfile1D(const TMatrixD &m, const TMatrixD &mErr,
		     TVectorD &v, TVectorD &vErr) {
  int res=
    GetMassProfile(m,0,v) &&
    GetMassProfile(mErr,0,vErr);
  return res;
}

// -----------------------------------------

// -----------------------------------------
// -----------------------------------------

void gatherAEff(const char *inpFile=NULL) {
  std::string line;
  std::vector<std::string> lines;

  if (!inpFile) {
    std::cout << "\ninpFile=NULL\n";
    std::cout << "Creating sample input file\n";
    const char *sampleFName="aeff_sample.inp";
    std::ofstream fout(sampleFName);
    fout << "../root_files/constants/DY_j22_19789pb/acceptance_constants1D.root\n";
    fout << "../root_files/constants/DY_j22_19789pb/event_efficiency_constants1D.root\n";
    fout << "../root_files/constants/DY_j22_19789pb/scale_factors_1D_Full2012_hltEffOld_PU.root\n";
    fout.close();
    std::cout << "check the file <" << sampleFName << ">\n";
    std::cout << "1D/2D in the name does not matter\n";
    return;
  }
  else { // load input
    std::ifstream finput(inpFile);
    if (!finput) {
      std::cout << "failed to open file <" << inpFile << ">\n";
      return;
    }
    getline(finput,line); lines.push_back(line);
    std::cout << "line1=<" << line << ">\n";
    getline(finput,line); lines.push_back(line);
    std::cout << "line2=<" << line << ">\n";
    getline(finput,line); lines.push_back(line);
    std::cout << "line3=<" << line << ">\n";
    finput.close();
  }
  std::cout << dashline;

  std::string fnameAcc,fnameEff,fnameRho;
  int count=0;
  for (unsigned int i=0; i<lines.size(); ++i) {
    line=lines[i];
    AdjustDim(line);
    if (PosOk(line,"acceptance")) { fnameAcc=line; count++; }
    else if (PosOk(line,"efficiency")) { fnameEff=line; count++; }
    else if (PosOk(line,"scale_factors")) { fnameRho=line; count++; }
    else {
      std::cout << "could not identify the file <" << line << ">\n";
      return;
    }
  }
  if (count!=3) { 
    std::cout << "not all files were identified\n";
    return;
  }
  std::cout << "files were identified ok\n";
  std::cout << dashline;

  TMatrixD effM(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD effErrM(effM);
  TMatrixD accM(effM), accErrM(effM);
  TMatrixD rhoM(effM), rhoErrM(effM);

  if (!LoadMatrix(fnameAcc,accM,accErrM,"acceptanceMatrix","acceptanceErrMatrix") ||
      !LoadMatrix(fnameEff,effM,effErrM,"efficiencyArray","efficiencyErrArray") ||
      !LoadMatrix(fnameRho,rhoM,rhoErrM,"scaleFactor","scaleFactorErr") ) {
    std::cout << "failed to load field\n";
    return;
  }

  TString outFName;
  if (DYTools::study2D) {
    outFName="dyee_aeff_2D.root";
    TFile fout(outFName,"recreate");
    accM.Write("acceptance");
    accErrM.Write("acceptanceErr");
    effM.Write("efficiency");
    effErrM.Write("efficiencyErr");
    rhoM.Write("scaleFactor");
    rhoErrM.Write("scaleFactorErr");
    unfolding::writeBinningArrays(fout);

    for (int i=0; i<DYTools::nMassBins; ++i) {
      TString massRange=Form("_%1.0lf_%2.0lf",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
      TString hAccName=TString("hAcc") + massRange;
      Histo_t *hAcc=extractRapidityDependence(hAccName,"",accM,accErrM,i,0);
      TString hEffName=TString("hEff") + massRange;
      Histo_t *hEff=extractRapidityDependence(hEffName,"",effM,effErrM,i,0);
      TString hRhoName=TString("hRho") + massRange;
      Histo_t *hRho=extractRapidityDependence(hRhoName,"",rhoM,rhoErrM,i,0);
      if (!hAcc || !hEff || !hRho) {
	std::cout << "got unexpected null histo\n";
	break;
      }

      hAcc->Write();
      hEff->Write();
      hRho->Write();
    }
    fout.Close();
  }
  else { // 1D case
    outFName="dyee_aeff_1D.root";
    std::cout << "accM: rows " << accM.GetNrows() << ", cols " << accM.GetNcols() << "\n";
    TVectorD eff(DYTools::nMassBins), effErr(eff);
    TVectorD acc(eff), accErr(eff);
    TVectorD rho(eff), rhoErr(eff);

    GetMassProfile1D(accM,accErrM, acc,accErr);
    GetMassProfile1D(effM,effErrM, eff,effErr);
    GetMassProfile1D(rhoM,rhoErrM, rho,rhoErr);
    
    Histo_t *hAcc=extractMassDependence("hAcc","",accM,accErrM,0,0,0);
    Histo_t *hEff=extractMassDependence("hEff","",effM,effErrM,0,0,0);
    Histo_t *hRho=extractMassDependence("hRho","",rhoM,rhoErrM,0,0,0);

    if (!hAcc || !hEff || !hRho) {
      std::cout << "got unexpected null histo\n";
    }
    else {
      TFile fout(outFName,"recreate");
      acc.Write("acceptance");
      accErr.Write("acceptanceErr");
      eff.Write("efficiency");
      effErr.Write("efficiencyErr");
      rho.Write("scaleFactor");
      rhoErr.Write("scaleFactorErr");
      unfolding::writeBinningArrays(fout);

      hAcc->Write();
      hEff->Write();
      hRho->Write();
      fout.Close();
    }
  }
  std::cout << "file <" << outFName << "> created\n";

  return;
}
