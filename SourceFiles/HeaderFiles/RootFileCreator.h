/*!
 *  \file RootFileCreator.h
 *  \brief Header file for some functions used to create *.root files
 */

#ifndef ROOTFILECREATOR_H
#define ROOTFILECREATOR_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <TROOT.h>
#include "ROOT/RDataFrame.hxx"
#include <TChainElement.h>
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

int RootFileCreatorDefine(pair<string, string> input_file_tree_name,
                          pair<string, string> output_file_tree_name,
                          unordered_map<string, vector<string>> anvar_tobedefined_variables,
                          string cut = "");
                          
int RootFileCreatorFilterer(pair<string, string> data_holder,
                            string output_filename,
                            string cut,
                            string output_path = "OutputFiles/",
                            chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                            bool outputprint = false,
                            bool debug = false);

int RootFileCreatorFilterer(TChain *treechain_to_filter,
                            string output_filename,
                            string cut,
                            int nevents=0,
                            chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                            bool outputprint = false,
                            bool debug = false);
                        
int RootFileCreatorTree(TTree *treetobewritten,
                        string output_filename = "",
                        string output_path = "OutputFiles/",
                        vector<string> branchestobewritten = {},
                        chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                        bool outputprint = false,
                        bool debug = false);

int RootFileCreatorTreeBranches(TTree *treetobewritten,
                                vector<string> branchestobewritten,
                                string output_filename = "",
                                string output_path = "OutputFiles/",
                                chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                                bool outputprint = false,
                                bool debug = false);

int RootFileCreatorExtensionPathPurger(string *stringtobepurged);

string RootFileCreatorExtensionPathPurger(string stringtobepurged);

#endif