/*!
 *  \file EfficienciesRatioAnalyzer.h
 *  \brief Header file for \ref EfficienciesRatioAnalyzer functions.
 * 
 *  \class These functions are meant to compute the ratio between efficiencies of the cuts in the simulation files 
 *         and efficiencies of the cuts in the data files
 *  
 */

#ifndef EFFICIENCIESRATIOANALYZER_H
#define EFFICIENCIESRATIOANALYZER_H

#include <TROOT.h>
#include <TH1.h>
#include <TLegend.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <ostream>

//It can be useful to use these namespaces
using namespace std;

int EfficienciesRatioAnalyzer(vector<AnalysisVariable> *vec_analysisvariables, ofstream& csv_efficiencies_ratios_file, unordered_map<string, pair<string, string>> data_holder);

#endif