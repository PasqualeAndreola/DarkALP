/*!
 *  \file TreeRecursiveSearch.cpp
 *  \brief Source file for \ref FakeNNReader function implementation
 */

/*!
 *  \fn int TreeReader(vector<MarkedNames> file2read, vector<MarkedNames> trees2read, vector<MarkedNames> var2read) "";
 *  \brief Function used to read selected variables stored in selected trees
 *
 *  \param file2read name of the file which stores the trees
 *  \param trees2read vector of names of trees that the user wants to explore
 */

#include "HeaderFiles/TreeRecursiveSearch.h"

int TreeRecursiveSearch(TDirectory *directory_to_explore, vector<pair<string, string>> *treeinthefile, TList *treelist)
{
    directory_to_explore->cd();
    TDirectory *current_directory = gDirectory;
    TIter nextkey(current_directory->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)nextkey()))
    {
        TObject *obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TDirectory::Class()))
        {
            TDirectory *new_directory = (TDirectory *)obj;
            cout << "Esploring the subfolder: " << new_directory->GetName() << endl << "\t";
            TreeRecursiveSearch(new_directory, treeinthefile, treelist);
        }
        if (obj->IsA()->InheritsFrom(TTree::Class()))
        {
            string treedirectory = string(directory_to_explore->GetName());
            string treename = string(obj->GetName());
            cout << "Tree " << treename << "found in directory " << treedirectory << endl;
            treeinthefile->push_back(pair(treedirectory, treename));
            treelist->Add((TTree *)obj);
        }
    }
    return 0;
};