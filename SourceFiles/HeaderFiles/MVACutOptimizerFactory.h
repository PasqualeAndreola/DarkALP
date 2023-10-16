/*!
 *  \file FakeNNFactory.h
 *  \brief Header file for \ref FakeNNFactory function
 */

#ifndef MVACUTOPTIMIZERFACTORY_H
#define MVACUTOPTIMIZERFACTORY_H

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
#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include <TMVA/Config.h>
#include <TMVA/Configurable.h>
#include <TMVA/TMVAGui.h>
#include "TMVAMethod.h"
#include "RootFileCreator.h"
#include "AnalysisVariable.h"

/*These namespaces can be useful*/
using namespace std;

int MVACutOptimizerFactory(unordered_map<string, pair<string, string>> mva_cut_factory_files,
                           vector<string> features,
                           vector<AnalysisVariable> *vartobefit,
                           vector<TMVAMethod> tmvamethods,
                           chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                           bool debug_variableranking = false);

int MVACutOptimizerFactoryMultiSigBkg(unordered_map<string, pair<string, string>> mva_cut_sig_factory_files,
                                      unordered_map<string, pair<string, string>> mva_cut_bkg_factory_files,
                                      string outputfilename,
                                      vector<string> features,
                                      vector<TMVAMethod> tmvamethods,
                                      chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                                      bool debug_variableranking = true);

#endif