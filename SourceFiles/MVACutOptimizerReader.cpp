/*!
 *  \file FakeNNReader.cpp
 *  \brief Source file for \ref FakeNNReader function implementation
 */

/*!
 *  \fn int TreeReader(vector<MarkedNames> file2read, vector<MarkedNames> trees2read, vector<MarkedNames> var2read) "";
 *  \brief Function used to read selected variables stored in selected trees
 *
 *  \param file2read name of the file which stores the trees
 *  \param trees2read vector of names of trees that the user wants to explore
 */

#include "HeaderFiles/MVACutOptimizerReader.h"
#include "HeaderFiles/RootFileCreator.h"
#include "HeaderFiles/TreeRecursiveSearch.h"

int MVACutOptimizerReader(string mva_cut_factory_filename,
                          pair<string, string> file_to_be_read,
                          vector<TMVAMethod> tmvamethods,
                          string outputfile_name,
                          bool testsamplefactoryefficency)
{
    // Looping over all the input files which should get a MVA cut variable

        string inputfilename = file_to_be_read.first;
        string inputfiletitle = RootFileCreatorExtensionPathPurger(file_to_be_read.first);
        string inputtreename = file_to_be_read.second;
        TFile *eventinputfile = TFile::Open(inputfilename.data(), "UPDATE");
        ROOT::RDataFrame *dataframe_inputfile = new ROOT::RDataFrame(inputtreename.data(), inputfilename.data());
        vector<string> features_inputfile = dataframe_inputfile->GetColumnNames();
        vector<string> mvareader_branches;

        TFile *mva_cut_factory_file = TFile::Open(mva_cut_factory_filename.data(), "write");
        if (mva_cut_factory_file->IsZombie() == true)
        {
            cout << "The MVA Cut Factory file does not exist\nAborting" << endl;
            return 0;
        }
        mva_cut_factory_file->cd();
        TDirectory *current_directory = gDirectory;
        vector<pair<string, string>> treeinthefactoryfile;
        TList tree_list_factoryfile;
        TreeRecursiveSearch(current_directory, &treeinthefactoryfile, &tree_list_factoryfile);

        // Creating a dataframe which holds only names of variables with a variance greater than the threshold set up before
        string mva_cut_factory_tree_signalname = treeinthefactoryfile[1].first + "/" + treeinthefactoryfile[1].second;
        if (testsamplefactoryefficency == true)    
        {    
            cout << mva_cut_factory_tree_signalname.data() << endl;
            ROOT::RDataFrame datafactoryframe = ROOT::RDataFrame(mva_cut_factory_tree_signalname.data(), mva_cut_factory_filename);
            cout << "Number of factory test signal events before the cut: " << datafactoryframe.Filter("classID==0").Count().GetValue() << endl;
            for (vector<TMVAMethod>::iterator tmvamethit = tmvamethods.begin(); tmvamethit < tmvamethods.end(); tmvamethit++)
            {            
                TString faccut = TString(("classID==0&&")+(tmvamethit->tmvamethodcut));
                cout << "Number of factory test signal events after the " << tmvamethit->tmvamethodname << " cut: " << datafactoryframe.Filter(faccut.Data()).Count().GetValue() << endl;
            }
        }
        ROOT::RDataFrame *dataframe_mva_factory = new ROOT::RDataFrame(mva_cut_factory_tree_signalname.data(), mva_cut_factory_filename);
        vector<string> features_reader = dataframe_mva_factory->GetColumnNames();
        vector<string> features_datareader;
        for (vector<string>::iterator inputdata_colname_iterator = features_inputfile.begin(); inputdata_colname_iterator != features_inputfile.end(); inputdata_colname_iterator++)
        {
            if (dataframe_mva_factory->HasColumn(inputdata_colname_iterator->data()) && dataframe_inputfile->HasColumn(inputdata_colname_iterator->data()))
            {
                features_datareader.push_back(inputdata_colname_iterator->data());
            }
        }
        vector<Double_t> Input_feature_var(features_datareader.size());
        unordered_map<string, Float_t> MVA_feature_var(features_datareader.size());

        int event_entries = dataframe_inputfile->Count().GetValue();

        for (vector<TMVAMethod>::iterator tmvamethit = tmvamethods.begin(); tmvamethit < tmvamethods.end(); tmvamethit++)
        {
            eventinputfile->cd();
            TTree *Event = (TTree *)eventinputfile->Get(inputtreename.data());
            for (vector<string>::iterator mvavar = features_datareader.begin(); mvavar < features_datareader.end(); mvavar++)
            {
                Event->SetBranchAddress(mvavar->data(), &(Input_feature_var[distance(features_datareader.begin(), mvavar)]));
            }
            string methoditname = string(tmvamethit->tmvamethodname);

            // Looking for any weight of the specified method in the folder containing the training results
            TString weightfilename = TString::Format("TMVAResults/%s/Weights/TMVAClassification_%s.weights.xml", "vt_transformed_dataset", methoditname.data());

            // Defining the reader that should do the MVA evaluation
            TMVA::Reader *reader = new TMVA::Reader(methoditname);
            for (vector<string>::iterator mvavar = features_datareader.begin(); mvavar != features_datareader.end(); mvavar++)
            {
                reader->AddVariable(mvavar->data(), &(MVA_feature_var[mvavar->data()]));
            }
            reader->BookMVA(methoditname, weightfilename);

            // Checking that input tree is read from the right file, because histograms are written to a different file
            eventinputfile->cd();
            //Event->SetName((string(Event->GetName()) + "_" + methoditname).data());
            Double_t mvacut_weight = 0, mva_cut_proba = 0, mva_cut_rarity = 0;
            char *branch_weight_name = (string("mvacut_") + methoditname).data();
            TBranch *Branch_Weight = Event->Branch(branch_weight_name, &mvacut_weight);
            mvareader_branches.push_back(string(branch_weight_name));
            char *branch_proba_name = (string("mvaproba_") + methoditname).data();
            TBranch *Branch_Proba = Event->Branch(branch_proba_name, &mva_cut_proba);
            mvareader_branches.push_back(string(branch_proba_name));            
            char *branch_rarity_name = (string("mvararity_") + methoditname).data();
            TBranch *Branch_Rarity = Event->Branch(branch_rarity_name, &mva_cut_rarity);
            mvareader_branches.push_back(string(branch_rarity_name));            
        
            for (int readerindex = 0; readerindex < event_entries; readerindex++)
            {
                Event->GetEntry(readerindex);
                for (vector<string>::iterator mvavar = features_datareader.begin(); mvavar < features_datareader.end(); mvavar++)
                {
                    MVA_feature_var[mvavar->data()] = (Double_t)Input_feature_var[distance(features_datareader.begin(), mvavar)];
                }
                mvacut_weight = reader->EvaluateMVA(methoditname);
                mva_cut_proba = reader->GetProba(methoditname);
                mva_cut_rarity = reader->GetRarity(methoditname);
                Branch_Weight->Fill();
                Branch_Proba->Fill();
                Branch_Rarity->Fill();
            }
            if (outputfile_name.compare("") == 0)
            {
                TFile *outputfile = TFile::Open(TString::Format("OutputFiles/%s_reader.root", inputfiletitle.data()), "RECREATE");
                outputfile->Close();
                RootFileCreatorTree(Event, inputfiletitle+string("_reader.root"));
            }
            else
            {
                TFile *outputfile = TFile::Open(outputfile_name.data(), "RECREATE");
                outputfile->Close();
                RootFileCreatorTree(Event, outputfile_name.data(), "OutputFiles/");
            }
            Branch_Weight->ResetAddress();
            Branch_Proba->ResetAddress();
            Branch_Rarity->ResetAddress();

            // Freeing memory allocated for the reader
            delete reader;
        }

    return 0;
}

int MVACutOptimizerReader(string mva_cut_factory_filename,
                          unordered_map<string, pair<string, string>> files_to_be_read,
                          vector<TMVAMethod> tmvamethods)
{
    // Looping over all the input files which should get a MVA cut variable
    for (unordered_map<string, pair<string, string>>::iterator fileit = files_to_be_read.begin(); fileit != files_to_be_read.end(); fileit++)
    {
        MVACutOptimizerReader(mva_cut_factory_filename, fileit->second, tmvamethods);
    }
    return 0;
}