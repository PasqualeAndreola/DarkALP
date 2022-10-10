/*! 
 *  \file VariablePlotter.h 
 *  \brief Header file for \ref VariablePlotter function
 */

#ifndef VARIABLEPLOTTER_H
#define VARIABLEPLOTTER_H

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
#include "AnalysisVariable.h"

/*These namespaces can be useful*/
using namespace std;

int VariablePlotter(unordered_map<string, pair<string, string>> data_holder,
                        vector<AnalysisVariable> *tobediscriminated,
                        string cuts,
                        chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                        bool debug = false);

#endif