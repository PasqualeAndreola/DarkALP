/*!
 *  \file EfficienciesAnalyzer.cpp
 *  \brief Source file for \ref EfficienciesAnalyzer function implementation
 */

/*!
 *  \fn int EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                        vector<AnalysisVariable> *var_to_be_analyzed,
                        bool debug) "";
 *  \brief Function used to import selected variables from the files in input. Efficiencies of the cuts of the variables
 *         are computed and written in an output file defined by the user.
 *
 *  \param data_holder Represents a list of pair of names: the first name of the pair represents the input file,
 *                     whereas the second name of the pair represent the name of the tree that holds the variables.
 *  \param *var_to_be_analyzed Vector of instances of the \ref AnalysisVariable class.
 *  \param cuts String that expresses the cut being used in the analysis
 *  \param start Starting time of the main function (default is the starting time of the function)
 */

#include "HeaderFiles/AnalysisVariable.h"
#include "HeaderFiles/PrintFuncInfo.h"
#include "HeaderFiles/RootFileCreator.h"

int AnalysisVariable::VariablePlotter(pair<string, string> *data_holder,
                                      AnalysisVariable *var_to_be_analyzed,
                                      string cut,
                                      chrono::_V2::system_clock::time_point start,
                                      bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    string filelabel = RootFileCreatorExtensionPathPurger(data_holder->first);
    ROOT::RDataFrame *dataframe_analysisvariable = NULL;
    TCanvas canvas_analysisvariable = TCanvas("canvas_analysisvariable", "canvas_analysisvariable", 1920, 1080);
    ROOT::RDF::RResultPtr<TH1F> histogram_analysisvariable;
    TLegend legend_analysisvariable;
    vector<TGraph> legend_roc_sigvsnorm_graph;
    vector<TString> legend_roc_sigvsnorm_entries;

    // Loading dataframes which hold the variables needed to compute the efficiencies
    bool fileexist = TFile(data_holder->first.data()).IsZombie();
    bool treexist = (TTree *)(TFile(data_holder->first.data()).Get(data_holder->second.data()))->IsZombie();
    if (fileexist == false && treexist == false)
        dataframe_analysisvariable = new ROOT::RDataFrame(data_holder->second.data(), data_holder->first.data());

    // Looping over chosen variables in order to compute the efficiencies
    Float_t bins = var_to_be_analyzed->variable_bins, min = var_to_be_analyzed->variable_histmin, max = var_to_be_analyzed->variable_histmax;
    if (dataframe_analysisvariable->HasColumn(var_to_be_analyzed->variable_name) == false)
        return 0;

    histogram_analysisvariable = dataframe_analysisvariable->Filter(cut.data()).Fill(TH1F(TString::Format("%s_sig", var_to_be_analyzed->variable_prettyname), TString::Format("%s_sig", var_to_be_analyzed->variable_name),
                                                                       bins, min, max),
                                                                  {var_to_be_analyzed->variable_name});

    // Computing efficiencies and purities
    histogram_analysisvariable->Sumw2();
    histogram_analysisvariable->Scale(1 / histogram_analysisvariable->Integral());

    // Plotting

    THStack histcompare_stack = THStack("Histstack", "Histstack");

    canvas_analysisvariable.cd();
    canvas_analysisvariable.Clear();
    histogram_analysisvariable->SetStats(false);
    histogram_analysisvariable->SetFillStyle(1001);
    histogram_analysisvariable->SetMarkerStyle(kFullSquare);
    histogram_analysisvariable->SetMarkerSize(1);
    histogram_analysisvariable->SetTitle(TString::Format("%s occurencies", var_to_be_analyzed->variable_prettyname));
    histogram_analysisvariable->GetXaxis()->SetTitle(var_to_be_analyzed->Xlabel());
    histogram_analysisvariable->GetYaxis()->SetTitle("Normalized occurencies");
    histogram_analysisvariable->Draw();
    canvas_analysisvariable.Update();

    // Adjusting the legend
    TLegend *legend = var_to_be_analyzed->SetLegendPosAuto("TR", 1);
    legend->AddEntry(histogram_analysisvariable->GetName(), TString::Format("%s", var_to_be_analyzed->variable_prettyname), "PLC PMC");
    legend->SetTextSize(0.025);
    legend->Draw("SAME");
    if (var_to_be_analyzed->variable_logscale_flag == true)
    {
        histogram_analysisvariable->SetMinimum(1e-4);
        histogram_analysisvariable->SetMaximum(1e-1);
        gPad->SetLogy();
    }
    canvas_analysisvariable.Update();
    canvas_analysisvariable.Print(TString::Format("%s/%s_%s.png", var_to_be_analyzed->variable_histplotfolder, var_to_be_analyzed->variable_name, filelabel.data()));
    canvas_analysisvariable.Clear();
    legend->Clear();
    gPad->SetLogy(0);

    return 0;
}