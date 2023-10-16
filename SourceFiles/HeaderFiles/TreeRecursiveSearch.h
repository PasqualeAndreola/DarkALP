/*!
 *  \file TreeRecursiveSearch.h
 *  \brief Header file for \ref TreeRecursiveSearch function
 */

#ifndef TREERECURSIVESEARCH_H
#define TREERECURSIVESEARCH_H

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
#include <TF1.h>
#include <TFrame.h>
#include <TRatioPlot.h>
#include <TPaveText.h>
#include <TMVA/DataLoader.h>
#include <TMVA/RReader.hxx>
#include <TMVA/Factory.h>
#include <TMVA/Config.h>
#include <TMVA/Configurable.h>
#include <TMVA/TMVAGui.h>
#include "TMVAMethod.h"
#include <TMVA/Reader.h>
#include <TLine.h>

/*These namespaces can be useful*/
using namespace std;

int TreeRecursiveSearch(TDirectory *directory_to_explore, vector<pair<string, string>> *treeinthefile, TList *treelist);

#endif