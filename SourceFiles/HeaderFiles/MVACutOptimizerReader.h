/*!
 *  \file MVACutOptimizerReader.h
 *  \brief Header file for \ref MVACutOptimizerReader function
 */

#ifndef MVACUTOPTIMIZERREADER_H
#define MVACUTOPTIMIZERREADER_H

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

int MVACutOptimizerReader(string mva_cut_factory_filename,
                          pair<string, string> file_to_be_read,
                          vector<TMVAMethod> tmvamethods,
                          string outputfile_name = "",
                          bool testsamplefactoryefficency = false);

int MVACutOptimizerReader(string mva_cut_factory_filename,
                          unordered_map<string, pair<string, string>> files_to_be_read,
                          vector<TMVAMethod> tmvamethods);

#endif