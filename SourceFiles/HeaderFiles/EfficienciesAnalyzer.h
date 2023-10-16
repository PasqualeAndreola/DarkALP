/*! 
 *  \file EfficienciesAnalyzer.h 
 *  \brief Header file for \ref EfficienciesAnalyzer function
 */

#ifndef EFFICIENCIESANALYZER_H
#define EFFICIENCIESANALYZER_H

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
#include <TVectorD.h>
#include "AnalysisVariable.h"

/*These namespaces can be useful*/
using namespace std;

Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                        unordered_map<string, string> mapofcuts,
                        chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                        bool outputprint = false,
                        bool debug = false);

Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                        string cut,
                        chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                        bool outputprint =false,
                        bool debug = false);

Double_t EfficienciesAnalyzer(pair<string, string> data_holder,
                              string cut,
                              chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                              bool outputprint = false,
                              bool debug = false);

int ROCCurve(pair<string, string> data_holder, AnalysisVariable *var_to_be_analyzed, string truth_cut, chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                        bool outputprint =false,
                        bool debug = false);
#endif