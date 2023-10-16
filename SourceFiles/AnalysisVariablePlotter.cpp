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

int AnalysisVariable::VariableComparer(unordered_map<string, pair<string, string>> *map_data_holder,
                                       AnalysisVariable *var_to_be_analyzed,
                                       string cut,
                                       chrono::_V2::system_clock::time_point start,
                                       bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataframe_analysisvariable = NULL;
    TCanvas canvas_analysisvariable = TCanvas("canvas_analysisvariable", "canvas_analysisvariable", 1920, 1080);
    ROOT::RDF::RResultPtr<TH1F> histogram_analysisvariable;
    Float_t bins = 25, min = 1e15, max = -1e15;
    THStack histogram_varcomparer;
    TLegend *legend;
    if (var_to_be_analyzed->variable_legendauto)
        legend = var_to_be_analyzed->SetLegendPosAuto(1, map_data_holder->size());
    else
        legend = var_to_be_analyzed->SetLegendPosAuto();

    for (unordered_map<string, pair<string, string>>::iterator data_holder = map_data_holder->begin(); data_holder != map_data_holder->end(); data_holder++)
    {
        // Loading dataframes which hold the variables needed to compute the efficiencies

        TFile inputfile = TFile(data_holder->second.first.data());
        if (inputfile.IsZombie())
        {
            cout << "The file " << inputfile.GetName() << " doesn't exist" << endl;
            return 0;
        }        
        TTree *inputtree = (TTree *)(inputfile.Get(data_holder->second.second.data()));
        if (inputtree->IsZombie())
        {
            cout << "The tree " << inputtree->GetName() << " in the file " << inputfile.GetName() << " doesn't exist" << endl;
            return 0;
        }
        if (inputtree->GetMinimum(var_to_be_analyzed->variable_name) < min)
        {
            min = inputtree->GetMinimum(var_to_be_analyzed->variable_name);
            min = round(min*1e5)/1e5;
        }
        if (inputtree->GetMaximum(var_to_be_analyzed->variable_name) > max)
        {
            max = inputtree->GetMaximum(var_to_be_analyzed->variable_name);
            max = round(max*1e5)/1e5;
        }
        inputfile.Close();
    }

    for (unordered_map<string, pair<string, string>>::iterator data_holder = map_data_holder->begin(); data_holder != map_data_holder->end(); data_holder++)
    {
        dataframe_analysisvariable = new ROOT::RDataFrame(data_holder->second.second.data(), data_holder->second.first.data());
        ROOT::RDF::RNode dataframe_analysisvariablenode = ROOT::RDF::RNode(*dataframe_analysisvariable);

        // Looping over chosen variables in order to compute the efficiencies
        string filelabel = RootFileCreatorExtensionPathPurger(data_holder->second.first);
        if (var_to_be_analyzed->variable_binauto == false)
        {
            bins = var_to_be_analyzed->variable_bins; 
            min = var_to_be_analyzed->variable_histmin;
            max = var_to_be_analyzed->variable_histmax;
        }
        Float_t binswidth = (max - min) / (bins);

        if (dataframe_analysisvariable->HasColumn(var_to_be_analyzed->variable_name) == false)
        {
            cout << "The variable \"" << var_to_be_analyzed->variable_name  << "\" is not in the dataset" << endl;
            try 
            {
                dataframe_analysisvariablenode = dataframe_analysisvariable->Define(var_to_be_analyzed->variable_name, var_to_be_analyzed->variable_expression);
                cout << "The variable \"" << var_to_be_analyzed->variable_name  << "\" has been defined as " << var_to_be_analyzed->variable_expression << endl;
            }
            catch (...)
            {
                cout << "Not able to define the new column" << endl;
                return 0;
            }
        }    

        if (cut.compare("") != 0)
        {
            cout << "Cutting the events with the cut: " << endl;
            cout << cut << endl;
            histogram_analysisvariable = dataframe_analysisvariablenode.Filter(cut.data()).Fill(TH1F(TString::Format("%s%s", var_to_be_analyzed->variable_name, filelabel.data()), 
                                                                                             TString::Format("%s%s", var_to_be_analyzed->variable_name, filelabel.data()), 
                                                                                             bins, min, max),
                                                                                             {var_to_be_analyzed->variable_name});
        }
        else
        {
            histogram_analysisvariable = dataframe_analysisvariablenode.Fill(TH1F(TString::Format("%s%s", var_to_be_analyzed->variable_name, filelabel.data()), 
                                                                          TString::Format("%s%s", var_to_be_analyzed->variable_name, filelabel.data()),
                                                                          bins, min, max),
                                                                          {var_to_be_analyzed->variable_name});
        }

        cout << "Plotting the variable " << var_to_be_analyzed->variable_name << " in the file " << filelabel << endl;

        // Computing efficiencies and purities
        histogram_analysisvariable->Sumw2();
        if (var_to_be_analyzed->variable_plot_normalized_flag==true)
            histogram_analysisvariable->Scale(1 / histogram_analysisvariable->Integral());

        // Plotting
        canvas_analysisvariable.cd();
        if (map_data_holder->size() == 1) 
        {
            histogram_analysisvariable->SetStats(true);
            histogram_analysisvariable->Draw();
            gPad->Update();
        }
        else
            histogram_analysisvariable->SetStats(false);
        gPad->Clear();
        gPad->Modified();
        histogram_analysisvariable->SetFillStyle(1001);
        histogram_analysisvariable->SetMarkerStyle(kFullSquare);
        histogram_analysisvariable->SetMarkerSize(1);
        histogram_varcomparer.Add((TH1F*)histogram_analysisvariable.GetValue().Clone());
        histogram_varcomparer.Draw("HF PLC PMC NOSTACK");
        histogram_varcomparer.SetTitle(TString::Format("%s distribution comparison", var_to_be_analyzed->variable_prettyname));
        histogram_varcomparer.GetXaxis()->SetTitle(var_to_be_analyzed->Xlabel(var_to_be_analyzed->variable_prettyname, var_to_be_analyzed->variable_dimension));
        histogram_varcomparer.GetYaxis()->SetTitle(TString::Format("Events / #left(%.5f %s#right)", binswidth, var_to_be_analyzed->variable_dimension).Data());
        canvas_analysisvariable.Update();
        gPad->Modified();

        // Adjusting the legend
        Float_t padtextsize = var_to_be_analyzed->variable_padtextsize;
        legend->AddEntry(histogram_analysisvariable->GetName(), TString::Format("%s (%s)", var_to_be_analyzed->variable_prettyname, data_holder->first.data()), "PLC PMC");
        legend->SetTextSize(padtextsize);
        legend->Draw("SAME");
        if (var_to_be_analyzed->variable_logscale_flag==true)
        {
            cout << "y Log Scale Mode: ON" << endl;
            histogram_analysisvariable->SetMinimum(1e-4);
            histogram_analysisvariable->SetMaximum(1e-1);
            gPad->SetLogy();
        }
        canvas_analysisvariable.Update();
        if (map_data_holder->size() == 1)
        {    
            TPaveStats *pavestat = var_to_be_analyzed->SetStatAuto(histogram_analysisvariable.GetPtr(), legend);
            pavestat->SetOptStat(var_to_be_analyzed->variable_statoptions);
            pavestat->SetTextSize(padtextsize);
            pavestat->Draw("SAME");
        }
        canvas_analysisvariable.Update();
        histogram_analysisvariable->Delete();
    }
    canvas_analysisvariable.Print(TString::Format("%s/%s_DistroComparison.png", var_to_be_analyzed->variable_histplotfolder, var_to_be_analyzed->variable_name));
    canvas_analysisvariable.Clear();

    histogram_varcomparer.RecursiveRemove(histogram_varcomparer.GetHists());
    legend->Clear();
    gPad->SetLogy(0);

    return 0;
}

int AnalysisVariable::VariablePlotter(pair<string, string> *data_holder, AnalysisVariable *var_to_be_analyzed, string cut, chrono::_V2::system_clock::time_point start, bool debug)
{
    unordered_map<string, pair<string, string>> map_data_holder;
    map_data_holder["FILE"] = *data_holder;
    VariableComparer(&map_data_holder, var_to_be_analyzed, cut, chrono::system_clock::now(), false);
    return 0;
}

int AnalysisVariable::VariablePlotter2D(pair<string, string> *data_holder,
                                        AnalysisVariable *var_to_be_analyzed,
                                        AnalysisVariable *var_to_be_analyzed2,
                                        string cut,
                                        bool variable_logscale_zaxis_flag,
                                        chrono::_V2::system_clock::time_point start,
                                        bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    string filelabel = RootFileCreatorExtensionPathPurger(data_holder->first);
    ROOT::RDataFrame *dataframe_analysis = NULL;
    TCanvas canvas_analysisvariable = TCanvas("canvas_analysisvariable", "canvas_analysisvariable", 1920, 1080);
    ROOT::RDF::RResultPtr<TH2F> histogram_analysisvariable;
    TLegend legend_analysisvariable;

    // Loading dataframes which hold the variables needed to compute the efficiencies
    bool fileexist = TFile(data_holder->first.data()).IsZombie();
    bool treexist = (TTree *)(TFile(data_holder->first.data()).Get(data_holder->second.data()))->IsZombie();
    if (fileexist == false && treexist == false)
        dataframe_analysis = new ROOT::RDataFrame(data_holder->second.data(), data_holder->first.data());

    // Looping over chosen variables in order to compute the efficiencies
    Float_t bins = var_to_be_analyzed->variable_bins, min = var_to_be_analyzed->variable_histmin, max = var_to_be_analyzed->variable_histmax;
    Float_t bins2 = var_to_be_analyzed2->variable_bins, min2 = var_to_be_analyzed2->variable_histmin, max2 = var_to_be_analyzed2->variable_histmax;
    string x_varname = var_to_be_analyzed->variable_name;
    string x_varexpr = var_to_be_analyzed->variable_expression;    
    string y_varname = var_to_be_analyzed2->variable_name;
    string y_varexpr = var_to_be_analyzed2->variable_expression;
    bool dataframe_has_xvar = dataframe_analysis->HasColumn(x_varname);
    bool dataframe_has_yvar = dataframe_analysis->HasColumn(y_varname);
    auto dataframe_analysisvariable = std::make_unique<ROOT::RDF::RNode>(dataframe_analysis->Define("sqrt_mvacut_kbdt1", "sqrt(mvacut_kbdt1)"));
        
    if (dataframe_has_xvar==false)
    {
        dataframe_analysisvariable = std::make_unique<ROOT::RDF::RNode>(dataframe_analysisvariable->Define(x_varname, x_varexpr));
        cout << "Defining x variable" << endl;
        //return 0;
    }    
    if (dataframe_has_yvar==false)
    {
        dataframe_analysisvariable = std::make_unique<ROOT::RDF::RNode>(dataframe_analysisvariable->Define(y_varname, y_varexpr));
        cout << "Defining y variable" << endl;
    }
    
    Long_t histogram_analysisvariable_filteredentries = 0;
    if (cut.compare("") != 0)
    {
        histogram_analysisvariable_filteredentries = dataframe_analysisvariable->Filter(cut.data()).Count().GetValue(); 
        cout << "Entries of filtered dataset used for 2D histogram: " << histogram_analysisvariable_filteredentries << endl;
        histogram_analysisvariable = dataframe_analysisvariable->Filter(cut.data()).Fill(TH2F(TString::Format("%s", filelabel.data()), TString::Format("%s", filelabel.data()), bins, min, max, bins2, min2, max2), {x_varname, y_varname});
    }
    else
    {
        histogram_analysisvariable = dataframe_analysisvariable->Fill(TH2F(TString::Format("%s", filelabel.data()), TString::Format("%s", filelabel.data()),
                                                                                               bins, min, max, bins2, min2, max2),
                                                                                          {x_varname, y_varname});
    }
    // Computing efficiencies and purities
    histogram_analysisvariable->Sumw2();
    // histogram_analysisvariable->Scale(1 / histogram_analysisvariable->Integral());

    // Plotting

    canvas_analysisvariable.cd();
    canvas_analysisvariable.Clear();
    histogram_analysisvariable->SetStats(false);
    histogram_analysisvariable->SetFillStyle(1001);
    histogram_analysisvariable->SetMarkerStyle(kFullSquare);
    histogram_analysisvariable->SetMarkerSize(1);
    histogram_analysisvariable->SetTitle(TString::Format("%s occurencies", var_to_be_analyzed->variable_prettyname));
    histogram_analysisvariable->GetXaxis()->SetTitle(var_to_be_analyzed->Xlabel());
    histogram_analysisvariable->GetYaxis()->SetTitle(var_to_be_analyzed2->Xlabel());
    histogram_analysisvariable->Draw("COLZ");
    canvas_analysisvariable.Update();

    if (variable_logscale_zaxis_flag == true)
    {   
        cout << "Z Log Scale Mode: ON" << endl;
        gPad->SetLogz();
        canvas_analysisvariable.Update();
    }
    // Adjusting the legend
    /*    TLegend *legend = var_to_be_analyzed->SetLegendPosAuto();
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

        TPaveStats *pavestat = (TPaveStats *)(histogram_analysisvariable->GetListOfFunctions()->FindObject("stats"));
        pavestat->SetOptStat(111111);
        pavestat->SetX1NDC(legend->GetX1NDC());
        pavestat->SetX2NDC(legend->GetX2NDC());
        pavestat->SetY1NDC(legend->GetY1NDC() - 0.04 * 5);
        pavestat->SetY2NDC(legend->GetY1NDC());
        pavestat->Draw("SAME");
        canvas_analysisvariable.Update();
    */
    canvas_analysisvariable.Print(TString::Format("%s/%svs%s_%s.png", var_to_be_analyzed->variable_histplotfolder, x_varname.data(), y_varname.data(), filelabel.data()));
    canvas_analysisvariable.Clear();
    // legend->Clear();
    gPad->SetLogy(0);

    return 0;
}