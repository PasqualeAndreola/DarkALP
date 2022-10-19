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
#include "HeaderFiles/EfficienciesRatioAnalyzer.h"
#include "HeaderFiles/PrintFuncInfo.h"

int EfficienciesRatioAnalyzer(vector<AnalysisVariable> *vec_analysisvariables, ofstream& csv_efficiencies_ratios_file, unordered_map<string, pair<string, string>> data_holder)
{
    Double_t efficiencyratio;
    unordered_map<string, Double_t> SingleCutEfficiency_Dataset;

    csv_efficiencies_ratios_file.open("Efficiencies_Ratios_vs_Cut.csv");
    csv_efficiencies_ratios_file << "Cut, "
                                    "KPi3Pi_SIM_100ps, Ratio(KPi3Pi_SIM_100ps), KPi3Pi_SIM_1ps, Ratio(KPi3Pi_SIM_1ps), KPi3Pi_SIM_Omega, Ratio(KPi3Pi_SIM_Omega), KPi3Pi_2015_md_minCuts, ,"
                                    "KPiEtaPiPi_SIM_100ps, Ratio(KPiEtaPiPi_SIM_100ps), KPiEtaPiPi_SIM_1ps, Ratio(KPiEtaPiPi_SIM_1ps), KPiEtaPiPi_SIM_Etapr, Ratio(KPiEtaPiPi_SIM_Etapr), KPiEtaPiPi_2015_md_minCuts" 
                                << endl; 

    for (vector<AnalysisVariable>::iterator analvar_iterator = vec_analysisvariables->begin(); analvar_iterator != vec_analysisvariables->end(); analvar_iterator++)
    {
        string singlecut = "";
        if (analvar_iterator->variable_stringcut.compare("")!=0) 
        { 
        singlecut = analvar_iterator->variable_stringcut;
        }
        else
        {
        singlecut = analvar_iterator->variable_cut;
        }
        if ((analvar_iterator+2)==(vec_analysisvariables->end()))
        {
        csv_efficiencies_ratios_file << "Minimal cuts";
        }
        else if ((analvar_iterator+1)==(vec_analysisvariables->end()))
        {
        csv_efficiencies_ratios_file << "All cuts";
        }
        else
        {
        if (analvar_iterator->variable_cutchecker == false) csv_efficiencies_ratios_file << singlecut.data() << "(Minimal)";
        else csv_efficiencies_ratios_file << singlecut.data();
        }
        for (unordered_map<string, pair<string, string>>::iterator datahold_iterator = data_holder.begin(); datahold_iterator != data_holder.end(); datahold_iterator++)
        { 
        string datait_key = datahold_iterator->first;
        if (datait_key.find("KPiEtaPiPi") != string::npos)
        {
            SingleCutEfficiency_Dataset[datait_key] = EfficienciesAnalyzer(datahold_iterator->second, singlecut);
        }
        if (datait_key.find("KPi3Pi") != string::npos)
        {
            SingleCutEfficiency_Dataset[datait_key] = EfficienciesAnalyzer(datahold_iterator->second, singlecut);
        }
        }

        for (unordered_map<string, Double_t>::iterator eff_iterator = SingleCutEfficiency_Dataset.begin(); eff_iterator != SingleCutEfficiency_Dataset.end(); eff_iterator++)
        {
        if (eff_iterator->first.find("KPi3Pi_SIM") != string::npos)
        {
            efficiencyratio = (eff_iterator->second)/(SingleCutEfficiency_Dataset["KPi3Pi_2015_md_minCuts"]);
            csv_efficiencies_ratios_file << "," << eff_iterator->second << "," << efficiencyratio;
            printf("Simulation efficiency with the %s cut:\n", singlecut.data());
            printf("\033[1;36m%s\033[m efficiency = \033[1;32m%.4f\033[0m\n", eff_iterator->first.data(), eff_iterator->second);
            printf("Ratio of efficiencies with the %s cut:\n", singlecut.data());
            printf("\033[1;36m%s\033[m efficiency over \033[1;36mKPi3Pi_2015_md_minCuts\033[0m efficiency = \033[1;32m%.4f\033[0m\n\n", eff_iterator->first.data(), efficiencyratio);
        }
        }
        csv_efficiencies_ratios_file << "," << SingleCutEfficiency_Dataset["KPi3Pi_2015_md_minCuts"] << ","; 
        printf("Data efficiency with the %s cut:\n", singlecut.data());
        printf("\033[1;36mKPi3Pi_2015_md_minCuts\033[m efficiency = \033[1;32m%.4f\033[0m\n", SingleCutEfficiency_Dataset["KPi3Pi_2015_md_minCuts"]);
        for (unordered_map<string, Double_t>::iterator eff_iterator = SingleCutEfficiency_Dataset.begin(); eff_iterator != SingleCutEfficiency_Dataset.end(); eff_iterator++)
        {
        if (eff_iterator->first.find("KPiEtaPiPi_SIM") != string::npos)
        {
            efficiencyratio = (eff_iterator->second)/(SingleCutEfficiency_Dataset["KPiEtaPiPi_2015_md_minCuts"]);
            csv_efficiencies_ratios_file << "," << eff_iterator->second << "," << efficiencyratio;
            printf("Simulation efficiency with the %s cut:\n", singlecut.data());
            printf("\033[1;36m%s\033[m efficiency = \033[1;32m%.4f\033[0m\n", eff_iterator->first.data(), eff_iterator->second); 
            printf("Ratio of efficiencies with the %s cut:\n", singlecut.data());
            printf("\033[1;36m%s\033[m efficiency over \033[1;36mKKPiEtaPiPi_2015_md_minCuts\033[0m efficiency = \033[1;32m%.4f\033[0m\n\n", eff_iterator->first.data(), efficiencyratio); 
        }
        }
        csv_efficiencies_ratios_file << "," << SingleCutEfficiency_Dataset["KPiEtaPiPi_2015_md_minCuts"]; 
        printf("Data efficiency with the %s cut:\n", singlecut.data());
        printf("\033[1;36mKPiEtaPiPi_2015_md_minCuts\033[m efficiency = \033[1;32m%.4f\033[0m\n", SingleCutEfficiency_Dataset["KPiEtaPiPi_2015_md_minCuts"]);
        csv_efficiencies_ratios_file << endl;
    }
    csv_efficiencies_ratios_file.close();
    cout << endl;

    return 0;

}