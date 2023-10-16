/*!
 *  \file TreeSkimmer.h
 *  \brief Header file for some functions used to create *.root files
 */

#ifndef TREESKIMMER_H
#define TREESKIMMER_H

#include <iostream>
#include <fstream>
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

int TreeSkimmer(string infile, string cuts="");
                          
#endif