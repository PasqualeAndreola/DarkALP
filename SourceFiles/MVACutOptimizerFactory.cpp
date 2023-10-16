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

int MVACutOptimizerFactory(unordered_map<string, pair<string, string>> mva_cut_factory_files,
                  vector<string> features,
                  vector<AnalysisVariable> *vartobefit,
                  vector<TMVAMethod> tmvamethods,
                  chrono::_V2::system_clock::time_point start,
                  bool debug_variableranking)
{
    // Initializing some variables 
    int num_threads = 0;   // use by default all threads
    // do enable MT running
    if (num_threads >= 0) {
        ROOT::EnableImplicitMT(num_threads);
        if (num_threads > 0) gSystem->Setenv("OMP_NUM_THREADS", TString::Format("%d",num_threads));
    }
    else
        gSystem->Setenv("OMP_NUM_THREADS", "1");
  
    TMVA::Config::Instance();
  
    std::cout << "Running with nthreads  = " << ROOT::GetThreadPoolSize() << std::endl;
    
    // Printing on terminal the name of the file used as signal/background
    cout << mva_cut_factory_files["Signal"].first.data() << endl;

    // Creating dataframes from the signal/background
    ROOT::RDataFrame* signal_dataframe = new ROOT::RDataFrame(mva_cut_factory_files["Signal"].second.data(), mva_cut_factory_files["Signal"].first.data());
    ROOT::RDataFrame* background_dataframe = new ROOT::RDataFrame(mva_cut_factory_files["Background"].second.data(), mva_cut_factory_files["Background"].first.data());

    Float_t closure_test_threshold = 1;

    cout << "" << endl;
    cout << "Model" << endl;
    cout << "" << endl;

    vector<string> features_ct = features;
    features_ct.push_back("Target");

    // Defining input/output files and trees
    TFile *input_signal_file = TFile::Open(mva_cut_factory_files["Signal"].first.data(), "read");
    TTree *signal = (TTree *)input_signal_file->Get(mva_cut_factory_files["Signal"].second.data());
    TFile *input_background_file = TFile::Open(mva_cut_factory_files["Background"].first.data(), "read");
    TTree *background = (TTree *)input_background_file->Get(mva_cut_factory_files["Background"].second.data());
    TFile *outputfile = TFile::Open(mva_cut_factory_files["Output"].first.data(), "recreate");
    
    // Setting up the MVA Factory
    Double_t signalweight = 1.0, backgroundweight = 1.0;
    TMVA::gConfig().GetIONames().fWeightFileDirPrefix = "TMVAResults";
    TMVA::gConfig().GetIONames().fWeightFileDir = "Weights";
    TMVA::DataLoader *loader = new TMVA::DataLoader("");
    TMVA::Factory *factory = new TMVA::Factory("TMVAClassification", outputfile, "AnalysisType=Classification");
    loader->AddSignalTree(signal, signalweight);
    loader->AddBackgroundTree(background, backgroundweight);
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
            variablesnan = variablesnan || (mvavar->find("iCand") != string::npos);
            variablesnan = variablesnan || (mvavar->find("k_ip_chi2") != string::npos);
            variablesnan = variablesnan || (mvavar->find("k_ip") != string::npos);
            variablesnan = variablesnan || (mvavar->find("k_p") != string::npos);
            variablesnan = variablesnan || (mvavar->find("pi_ip_chi2") != string::npos);
            variablesnan = variablesnan || (mvavar->find("pi_ip") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_fd_bvtx") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_fd_chi2") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_ip_chi2") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_ip") != string::npos);
            variablesnan = variablesnan || (mvavar->find("x_p") != string::npos);
            variablesnan = variablesnan || (mvavar->find("xpim_ip") != string::npos);
            variablesnan = variablesnan || (mvavar->find("xpim_p") != string::npos);
            variablesnan = variablesnan || (mvavar->find("xpip_ip") != string::npos);
            variablesnan = variablesnan || (mvavar->find("xpip_p") != string::npos);
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
    long unsigned int ntrainsign = signal->GetEntries()*train_test_treshold; /*!< Number of events used to train signal identification in the sample */
    long unsigned int ntrainback = background->GetEntries()*train_test_treshold/51; /*!< Number of events used to train background identification in the sample */
    long unsigned int ntestsign = signal->GetEntries()*(1-train_test_treshold);  /*!< Number of events used to test signal identification in the sample */
    long unsigned int ntestback = background->GetEntries()*(1-train_test_treshold)/51;  /*!< Number of events used to test background identification in the sample */
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
    input_signal_file->Close();
    input_background_file->Close();
    outputfile->Close();
    if(!gROOT->IsBatch()) TMVA::TMVAGui(outputfile->GetName());

    return 0;
}