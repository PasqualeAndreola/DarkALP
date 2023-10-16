/*!
 *  \file RootFileCreator.cpp
 *  \brief Source file for some functions used to create *.root files
 */

/*!
 *  \fn RootFileCreatorFilterer(pair<string, string> data_holder,
                            string output_file,
                            string cut,
                            chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                            bool outputprint = false,
                            bool debug = false) "";
 *  \brief Function used to produce *.root files choosing cuts or defining new columns
 *
 *  \param data_holder The first name of the pair represents the input file\n
 *                     The second name of the pair represents the name of the tree that holds the variables
 *  \param output_filename The name of the root file that will be produced (without the extension). If not specified, it will be the original filename+"_CREATED"
 *  \param cut Strings that express the cuts being used to produce the new *.root file
 *  \param output_path The name of the path where the output file will be written to (default is "OutputFiles/")
 *  \param start Starting time of the main function (default is the starting time of the function itself)
 *  \param outputprint Flag to decide to print the output to terminal
 *  \param debug Flag to acrivate debug mode
 */

#include "HeaderFiles/AnalysisVariableLambdaConstructor.h"
#include "HeaderFiles/RootFileCreator.h"
#include "HeaderFiles/PrintFuncInfo.h"

int RootFileCreatorDefine(pair<string, string> input_file_tree_name,
                          pair<string, string> output_file_tree_name,
                          unordered_map<string, vector<string>> anvar_tobedefined_variables,
                          string cut)
{
    // Enabling implicit multi-threading
    ROOT::EnableImplicitMT();

    // Building the dataframe and the node which will hold the new variables
    const char *input_filename = input_file_tree_name.first.data();
    const char *input_treename = input_file_tree_name.second.data();
    ROOT::RDataFrame input_dataframe = ROOT::RDataFrame(input_treename, input_filename);
    ROOT::RDF::RNode output_dataframe = input_dataframe;
    if (cut.compare("")!=0)
        output_dataframe = input_dataframe.Filter(cut);
        
    // Defining the new columns
    //if (output_dataframe.HasColumn("b_m_jpsiconstrained")==false)
    //    output_dataframe = input_dataframe.Define("b_m_jpsiconstrained", b_m_jpsiconstrained_3body, anvar_tobedefined_variables["b_m_jpsiconstrained"]);
    //if (output_dataframe.HasColumn("b_m_omegaconstrained")==false)
    //    output_dataframe = output_dataframe.Define("b_cons_omegaxpiz_m_best", b_m_omegaconstrained_3body, anvar_tobedefined_variables["b_m_omegaconstrained"]);
    //if (output_dataframe.HasColumn("x_m_pi0constrained")==false)
    //    output_dataframe = output_dataframe.Define("x_m_pi0constrained", x_m_pionstrained_3body, anvar_tobedefined_variables["x_m_pi0constrained"]);
    //if (output_dataframe.HasColumn("b_m_jpsiconstrained_k_pi0misidentifiedask")==false)
    //    output_dataframe = output_dataframe.Define("b_m_jpsiconstrained_k_pi0misidentifiedask", b_m_jpsiconstrained_k_pi0misidentifiedask, anvar_tobedefined_variables["b_m_jpsiconstrained_k_pi0misidentifiedask"]);
    if (output_dataframe.HasColumn("s_kpi")==false)
        output_dataframe = output_dataframe.Define("s_kpi", invariant_mass_couple, anvar_tobedefined_variables["s_kpi"]);
    if (output_dataframe.HasColumn("s_xpi")==false)
        output_dataframe = output_dataframe.Define("s_xpi", invariant_mass_couple, anvar_tobedefined_variables["s_xpi"]);
    if (output_dataframe.HasColumn("s_xk")==false)
        output_dataframe = output_dataframe.Define("s_xk", invariant_mass_couple, anvar_tobedefined_variables["s_xk"]);
    if (output_dataframe.HasColumn("s_kxpi")==false)
        output_dataframe = output_dataframe.Define("s_kxpi", invariant_mass_couple, anvar_tobedefined_variables["s_kxpi"]);
    if (output_dataframe.HasColumn("s_pixpi")==false)
        output_dataframe = output_dataframe.Define("s_pixpi", invariant_mass_couple, anvar_tobedefined_variables["s_pixpi"]);
        
    // Defining the variables that has to be written in the tree
    vector<string> vars_tobe_written = output_dataframe.GetColumnNames();

    // Modifying write options of the root data frame to overwrite other trees
    ROOT::RDF::RSnapshotOptions snapopt;
    snapopt.fMode = "RECREATE";
    snapopt.fOverwriteIfExists = "TRUE";

    // Writing the new node with the variables to the output file
    const char *output_filename = output_file_tree_name.first.data();
    const char *output_treename = output_file_tree_name.second.data();
    output_dataframe.Snapshot(output_treename, output_filename, vars_tobe_written, snapopt);

    ROOT::DisableImplicitMT();

    return 0;
}

int RootFileCreatorFilterer(pair<string, string> data_holder,
                            string output_filename,
                            string cut,
                            string output_path,
                            chrono::_V2::system_clock::time_point start,
                            bool outputprint,
                            bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *dataset = NULL;

    // Defining the name of the file and the name of the tree
    char *filename = data_holder.first.data();
    char *treename = data_holder.second.data();

    // Purging the input filename from the input path and its extension to collect a pure name
    string filename_clean = data_holder.first;
    size_t inputpath_name_position = filename_clean.rfind("/");
    if (inputpath_name_position != string::npos)
        filename_clean.erase(0, inputpath_name_position + 1);
    size_t extension_name_position = filename_clean.find(".root");
    if (extension_name_position != string::npos)
        filename_clean.erase(extension_name_position, string(".root").length());
    if (outputprint == false)
        cout << filename_clean.data() << endl;

    // Checking if the file and the tree do exist
    bool filezombie = TFile(filename).IsZombie();
    bool treezombie = TTree(filename, treename).IsZombie();

    // If the file and the tree do exist, extract a RDataFrame
    if (filezombie == false && treezombie == false)
        dataset = new ROOT::RDataFrame(treename, filename);

    // Filtering the new dataset using the cut given as a input to the function
    ROOT::RDF::RNode dataset_filtered = dataset->Filter(cut);

    // Looping over the columns that will be stored
    vector<string> columnnames_todataset;
    auto columnnames_fromdataset = dataset_filtered.GetColumnNames();
    for (auto &&colname : columnnames_fromdataset)
        columnnames_todataset.push_back(colname);

    // Modifying write options of the root data frame to overwrite other trees
    ROOT::RDF::RSnapshotOptions snapopt;
    snapopt.fMode = "UPDATE";
    snapopt.fOverwriteIfExists = "TRUE";

    // Write the new dataset to a new file. The tree name is the original one
    // Check if an output file name is specified, otherwise "_CREATED" is appended to the original one
    if (output_filename.compare("") == 0)
        dataset_filtered.Snapshot(treename, (output_path + filename_clean + "_CREATED.root").data(), dataset_filtered.GetColumnNames(), snapopt);
    else
        dataset_filtered.Snapshot(treename, (output_path + output_filename).data(), columnnames_todataset, snapopt);

    // Disabling implicit Multi-threading
    ROOT::DisableImplicitMT();

    return 0;
}

int RootFileCreatorFilterer(TChain *treechain_to_filter,
                            string output_filename,
                            string cut,
                            int nevents,
                            chrono::_V2::system_clock::time_point start,
                            bool outputprint,
                            bool debug)
{
    // Enabling implicit Multi-threading
    //ROOT::EnableImplicitMT();

    // Defining useful quantities
    const char *treename = treechain_to_filter->GetName();
    TObjArray *infiles = treechain_to_filter->GetListOfFiles();
    TIter next(infiles);
    TChainElement *chEl=0;
    TChain treechain(treename);
    while (( chEl=(TChainElement*)next() )) 
    {
        string filename(chEl->GetTitle());
        treechain.Add(filename.data());
    }

    // Filtering the new dataset using the cut given as a input to the function
    ROOT::RDataFrame *dataset = new ROOT::RDataFrame(treechain);
    dataset->Display({"x_cons_xpiz_m_best"});
    ROOT::RDF::RNode dataset_filtered = ROOT::RDF::RNode(*dataset);
    if (cut.compare("") != 0)
        dataset_filtered = dataset->Filter(cut);

    // Looping over the columns that will be stored
    vector<string> columnnames_todataset;
    auto columnnames_fromdataset = dataset_filtered.GetColumnNames();
    for (auto &&colname : columnnames_fromdataset)
    {
        bool forbiddenvar = colname.find("xpiz_cons_")==string::npos;
        forbiddenvar = forbiddenvar && colname.find("xpiz_tau")==string::npos;
        forbiddenvar = forbiddenvar && (colname.find("cons")==string::npos || (colname.find("cons")!=string::npos && colname.find("best")!=string::npos) );
        //forbiddenvar = forbiddenvar && colname.compare("xpiz_tau_chi2")!=0;
        if (forbiddenvar==true)
        {
            columnnames_todataset.push_back(colname);
        }
    }

    // Modifying write options of the root data frame to overwrite other trees
    ROOT::RDF::RSnapshotOptions snapopt;
    snapopt.fMode = "RECREATE";
    snapopt.fOverwriteIfExists = "TRUE";

    // Write the new dataset to a new file. The tree name is the original one
    // Check if an output file name is specified, otherwise "_CREATED" is appended to the original one
    if (output_filename.compare("") == 0)
        dataset_filtered.Snapshot(treename, (string(treename) + "CHAIN_CREATED.root").data(), dataset_filtered.GetColumnNames(), snapopt);
    else
        dataset_filtered.Range(0,nevents).Snapshot(treename, (output_filename).data(), columnnames_todataset, snapopt);

    // Disabling implicit Multi-threading
    ROOT::DisableImplicitMT();

    return 0;
}

/*!
 *  \fn RootFileCreatorTree(TTree *treetobewritten,
                            string output_filename,
                            string output_path,
                            chrono::_V2::system_clock::time_point start = chrono::_V2::system_clock::now(),
                            bool outputprint = false,
                            bool debug = false) "";
 *  \brief Function used to produce *.root files choosing cuts or defining new columns
 *
 *  \param treetobewritten Pointer to the tree that has to be written in the file\n
 *  \param output_filename The name of the root file that will be produced (without the extension). If not specified, it will be the original filename+"_CREATED"
 *  \param output_path The name of the path where the output file will be written to (default is "OutputFiles/")
 *  \param start Starting time of the main function (default is the starting time of the function itself)
 *  \param outputprint Flag to decide to print the output to terminal
 *  \param debug Flag to acrivate debug mode
 */

int RootFileCreatorTree(TTree *treetobewritten,
                        string output_filename,
                        string output_path,
                        vector<string> branchestobewritten,
                        chrono::_V2::system_clock::time_point start,
                        bool outputprint,
                        bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Purging the input filename from the tree and its extension to collect a pure name
    string filename_clean;
    if (output_filename.compare("") == 0)
        filename_clean = treetobewritten->GetCurrentFile()->GetName();
    else
        filename_clean = output_filename;

    RootFileCreatorExtensionPathPurger(&filename_clean);
    if (outputprint == false)
        cout << filename_clean.data() << endl;

    // Opening the output file and changing directory
    TFile *outputfile;

    // Write the tree to a new file. The tree name is the original one
    // Check if an output file name is specified, otherwise "_CREATED" is appended to the original one
    if (output_filename.compare("") == 0)
        outputfile = TFile::Open(string(output_path + filename_clean + "_CREATED.root").data(), "UPDATE");
    else
        outputfile = TFile::Open(string(output_path + filename_clean + ".root").data(), "UPDATE");

    outputfile->cd();
    cout << "Tree current file: " << gFile->GetName() << endl;
    TTree *treetobewritten_clone = treetobewritten->CloneTree();
    if (branchestobewritten.size()>0)
    {
        for(vector<string>::iterator branchit = branchestobewritten.begin(); branchit != branchestobewritten.end(); branchit++)
        {
            treetobewritten->GetBranch(branchit->data())->Write();
        }
    }
    else
    {
        treetobewritten_clone->Write("", TObject::kOverwrite);
    }
    outputfile->Close();

    // Disabling implicit Multi-threading
    ROOT::DisableImplicitMT();

    return 0;
}

int RootFileCreatorTreeBranches(TTree *treetobewritten,
                                vector<string> branchestobewritten,
                                string output_filename,
                                string output_path,
                                chrono::_V2::system_clock::time_point start,
                                bool outputprint,
                                bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Purging the input filename from the tree and its extension to collect a pure name
    string filename_clean;
    if (output_filename.compare("") == 0)
        filename_clean = treetobewritten->GetCurrentFile()->GetName();
    else
        filename_clean = output_filename;

    RootFileCreatorExtensionPathPurger(&filename_clean);
    if (outputprint == false)
        cout << filename_clean.data() << endl;

    // Opening the output file and changing directory
    TFile *outputfile;

    // Write the tree to a new file. The tree name is the original one
    // Check if an output file name is specified, otherwise "_CREATED" is appended to the original one
    if (output_filename.compare("") == 0)
        outputfile = TFile::Open(string(output_path + filename_clean + "_CREATED.root").data(), "UPDATE");
    else
        outputfile = TFile::Open(string(output_path + filename_clean + ".root").data(), "UPDATE");

    outputfile->cd();
    cout << "Tree current file: " << gFile->GetName() << endl;
    TTree *treetobewritten_clone = treetobewritten->CloneTree();
    treetobewritten_clone->Write();
    outputfile->Close();

    // Disabling implicit Multi-threading
    ROOT::DisableImplicitMT();

    return 0;
}
int RootFileCreatorExtensionPathPurger(string *stringtobepurged)
{
    size_t inputpath_name_position = stringtobepurged->rfind("/");
    if (inputpath_name_position != string::npos)
        stringtobepurged->erase(0, inputpath_name_position + 1);
    size_t extension_name_position = stringtobepurged->find(".root");
    if (extension_name_position != string::npos)
        stringtobepurged->erase(extension_name_position, string(".root").length());

    return 0;
}

string RootFileCreatorExtensionPathPurger(string stringtobepurged)
{
    string stringpurged = stringtobepurged;
    size_t inputpath_name_position = stringpurged.rfind("/");
    if (inputpath_name_position != string::npos)
        stringpurged.erase(0, inputpath_name_position + 1);
    size_t extension_name_position = stringpurged.find(".root");
    if (extension_name_position != string::npos)
        stringpurged.erase(extension_name_position, string(".root").length());

    return stringpurged;
}