/*!
 *  \file Fitter.h
 *  \brief Header file for different fitting functions
 */

#ifndef FITTER_H
#define FITTER_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <TROOT.h>
#include "ROOT/RDataFrame.hxx"
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveStats.h>
#include <TLorentzVector.h>
#include <TColor.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TFrame.h>
#include <TRatioPlot.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include <RooFit.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooClassFactory.h>
#include <RooExtendPdf.h>
#include <RooProduct.h>
#include <RooProdPdf.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooMinimizer.h>
#include <RooFormulaVar.h>
#include <RooWrapperPdf.h>
#include <RooFitResult.h>
#include <TMinuit.h>
#include <RooConstVar.h>
#include <RooChebychev.h>
#include <RooArgusBG.h>
#include <RooFFTConvPdf.h>
#include <RooCBShape.h>
#include <RooCrystalBall.h>
#include <RooExponential.h>
#include "AnalysisVariable.h"

/*These namespaces can be useful*/
using namespace std;

int FitterMass(ROOT::RDataFrame *input_dataframe,
               AnalysisVariable *analvar_tobefit,
               bool debug = false,
               string closuretest = "");

#endif