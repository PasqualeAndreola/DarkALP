/*!
 *  \file FakeNNFactory.cpp
 *  \brief Source file for \ref FakeNNFactory function implementation
 */

/*!
 *  \fn int TreeReader(vector<MarkedNames> file2read, vector<MarkedNames> trees2read, vector<MarkedNames> var2read) "";
 *  \brief Function used to read selected variables stored in selected trees
 * 
 *  \param file2read name of the file which stores the trees
 *  \param trees2read vector of names of trees that the user wants to explore
 */

#include "HeaderFiles/MVACutOptimizerFactory.h"

int MVACutOptimizerFactoryMultiSigBkg(unordered_map<string, pair<string, string>> mva_cut_sig_factory_files,
                                      unordered_map<string, pair<string, string>> mva_cut_bkg_factory_files,
                                      string outputfilename,
                                      vector<string> features,
                                      vector<TMVAMethod> tmvamethods,
                                      chrono::_V2::system_clock::time_point start,
                                      bool debug_variableranking)
{
    // Initializing some variables 
    Int_t sig_nentries = 0, minsig_nentries = kMaxInt, bkg_nentries = 0;
    TFile *outputfile = TFile::Open(outputfilename.data(), "recreate");
    ROOT::RDataFrame *signal_dataframe = NULL, *background_dataframe = NULL;

    // Defining some variables
    vector<pair<string, string>> sigfactoryfiles, bkgfactoryfiles;
  
    // Enabling implicit Multi-threading
    //ROOT::EnableImplicitMT();
    std::cout << "Running with nthreads  = " << ROOT::GetThreadPoolSize() << std::endl;

    // Counting total entries of signal and background
    for (unordered_map<string, pair<string, string>>::iterator mvasigfileit = mva_cut_sig_factory_files.begin(); mvasigfileit != mva_cut_sig_factory_files.end(); mvasigfileit++)
    {
        TFile *input_signal_file = TFile::Open(mvasigfileit->second.first.data(), "read");
        TTree *signal = (TTree *)input_signal_file->Get(mvasigfileit->second.second.data());
        sig_nentries = signal->GetEntries();
        if (sig_nentries<=minsig_nentries)
          minsig_nentries = sig_nentries;
    }    
    for (unordered_map<string, pair<string, string>>::iterator mvabkgfileit = mva_cut_bkg_factory_files.begin(); mvabkgfileit != mva_cut_bkg_factory_files.end(); mvabkgfileit++)
    {
        TFile *input_bkg_file = TFile::Open(mvabkgfileit->second.first.data(), "read");
        TTree *bkg = (TTree *)input_bkg_file->Get(mvabkgfileit->second.second.data());
        bkg_nentries = bkg_nentries+bkg->GetEntries();
    }

    // Creating the signal and background tree chains. Number of entries chosen as described by the input maps    
    for (unordered_map<string, pair<string, string>>::iterator mvasigfileit = mva_cut_sig_factory_files.begin(); mvasigfileit != mva_cut_sig_factory_files.end(); mvasigfileit++)
    {
        string sig_filename = mvasigfileit->second.first.data();
        string sig_filetitle = RootFileCreatorExtensionPathPurger(mvasigfileit->second.first).data();

        // Printing on terminal the name of the file used as signal
        cout << "Signal file " << sig_filename;
        cout << " as " << stof(mvasigfileit->first)*100 << "% of the total signal events (";
        cout << sig_nentries << ")" << endl;

        // Creating dataframes from the signal
        signal_dataframe = new ROOT::RDataFrame(mvasigfileit->second.second.data(), mvasigfileit->second.first.data());
        string output_filename = (string(sig_filetitle) + string("_SIGFACTORY.root")).data();
        string output_test_filename = (string(sig_filetitle) + string("_SIGTEST.root")).data();
        string output_treename = mvasigfileit->second.second.data();
        vector<string> vars_tobe_written = signal_dataframe->GetColumnNames();

        ROOT::RDF::RSnapshotOptions snapopt;
        snapopt.fMode = "RECREATE";
        snapopt.fOverwriteIfExists = "TRUE";
        signal_dataframe->Range(0, round(stof(mvasigfileit->first)*minsig_nentries/2)).Snapshot(output_treename, output_filename, vars_tobe_written);
        signal_dataframe->Range(round(stof(mvasigfileit->first)*minsig_nentries/2), 0).Snapshot(output_treename, output_test_filename, vars_tobe_written);
        sigfactoryfiles.push_back(pair(output_filename, output_treename));
    }
    for (unordered_map<string, pair<string, string>>::iterator mvabkgfileit = mva_cut_bkg_factory_files.begin(); mvabkgfileit != mva_cut_bkg_factory_files.end(); mvabkgfileit++)
    {
        string bkg_filename = mvabkgfileit->second.first.data();
        string bkg_filetitle = RootFileCreatorExtensionPathPurger(mvabkgfileit->second.first).data();

        // Printing on terminal the name of the file used as background
        cout << "Background file " << bkg_filename;
        cout << " as " << stof(mvabkgfileit->first)*100 << "% of the total background events (";
        cout << bkg_nentries << ")" << endl;

        // Creating dataframes from the background
        background_dataframe = new ROOT::RDataFrame(mvabkgfileit->second.second.data(), mvabkgfileit->second.first.data());
        string output_filename = (string(bkg_filetitle) + string("_BKGFACTORY.root")).data();
        string output_test_filename = (string(bkg_filetitle) + string("_BKGTEST.root")).data();
        string output_treename = mvabkgfileit->second.second.data();
        vector<string> vars_tobe_written = background_dataframe->GetColumnNames();

        ROOT::RDF::RSnapshotOptions snapopt;
        snapopt.fMode = "RECREATE";
        snapopt.fOverwriteIfExists = "TRUE";
        background_dataframe->Range(0, round(stof(mvabkgfileit->first)*minsig_nentries/2*0.5)).Snapshot(output_treename, output_filename, vars_tobe_written);
        background_dataframe->Range(round(stof(mvabkgfileit->first)*minsig_nentries/2*0.5), 0).Snapshot(output_treename, output_test_filename, vars_tobe_written);
        bkgfactoryfiles.push_back(pair(output_filename, output_treename));
    }


    Float_t closure_test_threshold = 1;

    cout << "" << endl;
    cout << "Model" << endl;
    cout << "" << endl;

    vector<string> features_ct = features;
    features_ct.push_back("Target");

    // Setting up the MVA Factory
    Double_t signalweight = 1.0, backgroundweight = 1.0;
    TMVA::gConfig().GetIONames().fWeightFileDirPrefix = "TMVAResults";
    TMVA::gConfig().GetIONames().fWeightFileDir = "Weights";
    TMVA::DataLoader *loader = new TMVA::DataLoader("");
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputfile, "AnalysisType=Classification");
    for (vector<pair<string, string>>::iterator sigfileit=sigfactoryfiles.begin(); sigfileit!=sigfactoryfiles.end(); sigfileit++)
    {
        TFile *input_signal_file = TFile::Open(sigfileit->first.data(), "read");
        TTree *signal = (TTree *)input_signal_file->Get(sigfileit->second.data());
        loader->AddSignalTree(signal, signalweight);
    }
    for (vector<pair<string, string>>::iterator bkgfileit=sigfactoryfiles.begin(); bkgfileit!=sigfactoryfiles.end(); bkgfileit++)
    {
        TFile *input_bkg_file = TFile::Open(bkgfileit->first.data(), "read");
        TTree *bkg = (TTree *)input_bkg_file->Get(bkgfileit->second.data());
        loader->AddBackgroundTree(bkg, backgroundweight);
    }
    if (debug_variableranking == true)
    {
      vector<string> dataframenames = signal_dataframe->GetColumnNames();
      for (vector<string>::iterator mvavar = dataframenames.begin(); mvavar < dataframenames.end(); mvavar++)
      {
        if ( (signal_dataframe->HasColumn(mvavar->data()) == true) && ( (background_dataframe->HasColumn(mvavar->data())) == true) )
        {
            bool variablesnan = (mvavar->find("_ip_chi2_bvtx") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_m") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_e") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_q") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_mcor") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_x") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_y") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_z") != string::npos);            
            variablesnan = variablesnan || (mvavar->find("_px") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_py") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_pz") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_pid") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_hcal") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_l0") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_is_mu") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_pnn") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_prb_ghost") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_ecal") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_vid") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_best") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_tau") != string::npos);
            variablesnan = variablesnan || (mvavar->find("eventNumber") != string::npos);
            variablesnan = variablesnan || (mvavar->find("_MC") != string::npos);
            if (variablesnan==false) loader->AddVariable(mvavar->data(), mvavar->data(), "", 'F');
        }
      }
    }
    else
    {
      for (vector<string>::iterator mvavar = features.begin(); mvavar < features.end(); mvavar++)
      {
        if ( (signal_dataframe->HasColumn(mvavar->data()) == true) && ( (background_dataframe->HasColumn(mvavar->data())) == true) )
        {
            loader->AddVariable(mvavar->data(), mvavar->data(), "", 'F');
        }
      }
    }
    Float_t train_test_treshold = 0.8;
    long unsigned int ntrainsign = minsig_nentries*train_test_treshold/2; /*!< Number of events used to train signal identification in the sample */
    long unsigned int ntrainback = minsig_nentries*train_test_treshold/2; /*!< Number of events used to train background identification in the sample */
    long unsigned int ntestsign = minsig_nentries/2*(1-train_test_treshold);  /*!< Number of events used to test signal identification in the sample */
    long unsigned int ntestback = (minsig_nentries/2)*(1-train_test_treshold);  /*!< Number of events used to test background identification in the sample */
    TString dataString = TString::Format("nTrain_Signal=%lu", ntrainsign);
    dataString.Append(TString::Format(":nTrain_Background=%lu", ntrainback));
    dataString.Append(TString::Format(":nTest_Signal=%lu", ntestsign));
    dataString.Append(TString::Format(":nTest_Background=%lu", ntestback));
    dataString.Append(":SplitMode=Random:NormMode=NumEvents:!V");
    loader->PrepareTrainingAndTestTree("", "", dataString);
    TMVA::DataLoader *loader_with_variance_threshold = loader->VarTransform("VT(0.0)");
    //loader_with_variance_threshold->SetName("FakeReweighting");
    cout << loader_with_variance_threshold->GetName() << endl;

    for (vector<TMVAMethod>::iterator tmvamethit = tmvamethods.begin(); tmvamethit < tmvamethods.end(); tmvamethit++)
    {
      factory->BookMethod(loader_with_variance_threshold, tmvamethit->tmvamethodtype, tmvamethit->tmvamethodname, tmvamethit->tmvamethodconfiguration);
    }
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outputfile->Close();
    if(!gROOT->IsBatch()) TMVA::TMVAGui(outputfile->GetName());

    return 0;
}