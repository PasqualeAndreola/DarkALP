/*!
 *  \file EfficienciesAnalyzer.cpp
 *  \brief Source file for \ref EfficienciesAnalyzer function implementation
 */

/*!
 *  \fn Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
 *                     unordered_map<string, string> mapofcuts,
 *                      chrono::_V2::system_clock::time_point start,
 *                      bool outputprint,
 *                       bool debug) "";
 *  \brief Function used to import selected variables from the files in input. Efficiencies of the cuts of the variables
 *         are computed and written in an output file defined by the user.
 *
 *  \param data_holder Represents a list of pair of names: the first name of the pair represents the input file,
 *                     whereas the second name of the pair represent the name of the tree that holds the variables.
 *  \param mapofcuts Map of strings that express the cuts being used in the analysis
 *  \param start Starting time of the main function (default is the starting time of the function)
 */
#include "HeaderFiles/EfficienciesAnalyzer.h"
#include "HeaderFiles/PrintFuncInfo.h"
#include "HeaderFiles/EfficienciesAnalyzer.h"

Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                              unordered_map<string, string> mapofcuts,
                              chrono::_V2::system_clock::time_point start,
                              bool outputprint,
                              bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataset = NULL;

    // Loading dataframes which hold the variables needed to compute the efficiencies
    // Looping over chosen variables in order to compute the efficiencies
    for (unordered_map<string, pair<string, string>>::iterator dataset_iterator = data_holder.begin(); dataset_iterator != data_holder.end(); dataset_iterator++)
    {
        const char *datasetname = dataset_iterator->first.data();
        char *filename = dataset_iterator->second.first.data();
        char *treename = dataset_iterator->second.second.data();
        bool fileexist = TFile(filename).IsZombie();
        bool treeexist = TTree(filename, treename).IsZombie();
        if (fileexist == false && treeexist == false)
            dataset = new ROOT::RDataFrame(treename, filename);
        Double_t neventsprecut = dataset->Count().GetValue();
        Double_t neventsminicut = dataset->Filter(mapofcuts["MinCuts"]).Count().GetValue();
        Double_t efficiencymincut = neventsminicut / neventsprecut;
        Double_t neventsaftercut = dataset->Filter(mapofcuts["Cuts"]).Count().GetValue();
        Double_t efficiencycut = neventsaftercut / neventsprecut;
        Double_t efficiency_mincut_vs_cut = neventsaftercut / neventsminicut;
        if (outputprint == true)
        {
            cout << "Number of events before the cuts in " << datasetname << " file are: " << neventsprecut << endl;
            cout << "Number of events after the minimal cuts in " << datasetname << " file are: " << neventsminicut << endl;
            printf("\033[1;33mEfficiency of minimal cuts wrt no cuts in the file %s is: %.4f \033[0m\n", datasetname, efficiencymincut);
            cout << "Number of events after the cuts in " << datasetname << " file are: " << neventsaftercut << endl;
            printf("\033[1;32mEfficiency of cuts wrt no cuts in the file %s is: %.4f \033[0m\n", datasetname, efficiencycut);
            printf("\033[1;36mEfficiency of cuts wrt minimal cuts in the file %s is: %.4f \033[0m\n", datasetname, efficiency_mincut_vs_cut);
            cout << endl;
        }
    }

    cout << endl;
    return 0;
}

/*!
 *  \fn Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
 *                      string cut,
 *                      chrono::_V2::system_clock::time_point start,
 *                      bool outputprint,
 *                       bool debug) "";
 *  \brief Function used to import selected variables from the files in input. Efficiencies of the cuts of the variables
 *         are computed and written in an output file defined by the user.
 *
 *  \param data_holder Represents a list of pair of names: the first name of the pair represents the input file,
 *                     whereas the second name of the pair represent the name of the tree that holds the variables.
 *  \param Cut String that expresses the cuts being used in the analysis
 *  \param start Starting time of the main function (default is the starting time of the function)
 */

Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                              string cut,
                              chrono::_V2::system_clock::time_point start,
                              bool outputprint,
                              bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataset = NULL;

    // Loading dataframes which hold the variables needed to compute the efficiencies
    // Looping over chosen variables in order to compute the efficiencies
    for (unordered_map<string, pair<string, string>>::iterator dataset_iterator = data_holder.begin(); dataset_iterator != data_holder.end(); dataset_iterator++)
    {
        const char *datasetname = dataset_iterator->first.data();
        char *filename = dataset_iterator->second.first.data();
        char *treename = dataset_iterator->second.second.data();
        bool fileexist = TFile(filename).IsZombie();
        bool treeexist = TTree(filename, treename).IsZombie();
        if (fileexist == false && treeexist == false)
            dataset = new ROOT::RDataFrame(treename, filename);
        Double_t neventsprecut = dataset->Count().GetValue();
        Double_t neventsaftercut = dataset->Filter(cut).Count().GetValue();
        Double_t efficiencycut = neventsaftercut / neventsprecut;
        if (outputprint == true)
        {
            cout << "Number of events before the cuts in " << datasetname << " file are: " << neventsprecut << endl;
            cout << "Number of events after the cuts in " << datasetname << " file are: " << neventsaftercut << endl;
            printf("\033[1;32mEfficiency of cuts wrt no cuts in the file %s is: %.4f \033[0m\n", datasetname, efficiencycut);
            cout << endl;
        }
    }

    cout << endl;
    return 0;
}

/*!
 *  \fn Double_t EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
 *                      string cut,
 *                      chrono::_V2::system_clock::time_point start,
 *                      bool outputprint,
 *                       bool debug) "";
 *  \brief Function used to import selected variables from the files in input. Efficiencies of the cuts of the variables
 *         are computed and written in an output file defined by the user.
 *
 *  \param data_holder Represents a pair of strings (the first one represents the file name; the second one represents the tree name)
 *  \param Cut String that expresses the cuts being used in the analysis
 *  \param start Starting time of the main function (default is the starting time of the function)
 */

Double_t EfficienciesAnalyzer(pair<string, string> data_holder,
                              string cut,
                              chrono::_V2::system_clock::time_point start,
                              bool outputprint,
                              bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataset = NULL;

    // Loading dataframe using file and tree names
    char *filename = data_holder.first.data();
    char *treename = data_holder.second.data();
    bool fileexist = TFile(filename).IsZombie();
    bool treeexist = TTree(filename, treename).IsZombie();
    if (fileexist == false && treeexist == false) dataset = new ROOT::RDataFrame(treename, filename);
    Double_t neventsprecut = dataset->Count().GetValue();
    Double_t neventsaftercut = dataset->Filter(cut).Count().GetValue();
    Double_t efficiencycut = neventsaftercut / neventsprecut;
    if (outputprint == true)
    {
        cout << "Number of events before the cuts in " << filename << " file are: " << neventsprecut << endl;
        cout << "Number of events after the cuts in " << filename << " file are: " << neventsaftercut << endl;
        printf("\033[1;32mEfficiency of cuts wrt no cuts in the file %s is: %.4f \033[0m\n", filename, efficiencycut);
        cout << endl;
    }

    return efficiencycut;
}
int ROCCurve(pair<string, string> data_holder, AnalysisVariable *var_to_be_analyzed, string truth_cut, chrono::_V2::system_clock::time_point start, bool outputprint, bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataset = NULL;
    TH1F* histogram_analysisvariable = NULL;
    const char *variable = var_to_be_analyzed->variable_expression;
    const char *var_pretty_name = var_to_be_analyzed->variable_prettyname;
    Double_t varmincut = var_to_be_analyzed->variable_cut_min;
    string kinematic_min_cut = TString::Format("%s>%f", variable, varmincut).Data();
    int bins = var_to_be_analyzed->variable_bins;
    vector<Double_t> *signal_efficiency_points = new vector<Double_t>;
    vector<Double_t> *background_rejection_points = new vector<Double_t>;

    // Loading dataframe using file and tree names
    char *filename = data_holder.first.data();
    char *treename = data_holder.second.data();
    bool fileexist = TFile(filename).IsZombie();
    bool treeexist = TTree(filename, treename).IsZombie();
    if (fileexist == false && treeexist == false) dataset = new ROOT::RDataFrame(treename, filename);
    Double_t n_tot_events = dataset->Count().GetValue(), min = dataset->Min(variable).GetValue(), max = dataset->Max(variable).GetValue();
    histogram_analysisvariable = (TH1F*)dataset->Fill(TH1F(var_pretty_name, var_pretty_name, 
                                               bins, min, max),
                                               {variable}).GetPtr()->Clone();
    Double_t var_hist_stddev = histogram_analysisvariable->GetStdDev();
    Double_t sigeff = (dataset->Filter((kinematic_min_cut+"&&"+truth_cut).data()).Count().GetValue())/(n_tot_events);
    signal_efficiency_points->push_back(sigeff);
    Double_t bkgeff = (dataset->Filter(("!("+kinematic_min_cut+")&& !("+truth_cut+")").data()).Count().GetValue())/(n_tot_events);
    background_rejection_points->push_back(bkgeff);
    while (sigeff > 0)
    {
        if (0.95*sigeff > dataset->Filter((kinematic_min_cut+"&&"+truth_cut).data()).Count().GetValue()/(n_tot_events))
        {
            sigeff = (dataset->Filter((kinematic_min_cut+"&&"+truth_cut).data()).Count().GetValue())/(n_tot_events);
            signal_efficiency_points->push_back(sigeff);
            bkgeff = (dataset->Filter(("!("+kinematic_min_cut+")&& !("+truth_cut+")").data()).Count().GetValue())/(n_tot_events);
            background_rejection_points->push_back(bkgeff);
        }
        cout << "Sigeff " << sigeff << " bkgeff " << bkgeff << endl;
        varmincut = varmincut+var_hist_stddev/10;
        kinematic_min_cut = TString::Format("%s>%f", variable, varmincut);
    }
    TGraph* bkgeff_vs_sigeff = new TGraph(signal_efficiency_points->size(), signal_efficiency_points->data(), background_rejection_points->data());
    bkgeff_vs_sigeff->Draw();
    gPad->SaveAs("bkgeff_vs_sigeff.png");
    return 0;
}
