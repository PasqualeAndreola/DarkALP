/*!
 *  \file TreeSkimmer.cpp
 *  \brief Source file for \ref TreeSkimmer function implementation
 */

/*!
 *  \fn int TreeSkimmer(TTree *treetoberead, vector<MarkedNames> trees2read, vector<MarkedNames> var2read) "";
 *  \brief Function used to read selected variables stored in selected trees
 *
 *  \param file2read name of the file which stores the trees
 *  \param trees2read vector of names of trees that the user wants to explore
 */

#include "HeaderFiles/TreeSkimmer.h"

int TreeSkimmer(string infile, string cuts)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();
    // cout << lengthbranchnames << endl;
    // cout << lengthpurgedbranchnames << endl;
    vector<string>* data2015magup = new vector<string>;
    string type, year, polarity, linex;
    ifstream fin(infile.data());
    std::getline(fin, type, ';');
    std::getline(fin, year, ';');
    std::getline(fin, polarity, ';');
    while (fin)
    {   //^^if you do not always have a dummy at the beginning of line
        std::getline(fin, linex, '\n');
        if (linex.compare("")!=0)
        {
            data2015magup->push_back(linex);
            cout << linex << endl;
        }
    }
    
    unordered_map<vector<string>*, tuple<string, string, string>> fileyearpolarity =
    {
                //{"00183579_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2016", "MagDown"}},
                //{"00183577_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2016", "MagUp"}},
                //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185824/0000/00185824_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2016", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185822/0000/00185822_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2016", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2015/DVNTUPLE.ROOT/00185826/0000/00185826_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2015", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2015/DVNTUPLE.ROOT/00185828/0000/00185828_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2015", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185822/0000/00185822_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2016", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185824/0000/00185824_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2016", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2017/DVNTUPLE.ROOT/00185817/0000/00185817_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2017", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2017/DVNTUPLE.ROOT/00185820/0000/00185820_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2017", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2018/DVNTUPLE.ROOT/00185809/0000/00185809_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2018", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2018/DVNTUPLE.ROOT/00185813/0000/00185813_00000001_1.dvntuple.root", {"MC_B2Kpiomega", "2018", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2015/DVNTUPLE.ROOT/00185842/0000/00185842_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2015", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2015/DVNTUPLE.ROOT/00185844/0000/00185844_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2015", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185838/0000/00185838_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2016", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/DVNTUPLE.ROOT/00185840/0000/00185840_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2016", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2017/DVNTUPLE.ROOT/00185834/0000/00185834_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2017", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2017/DVNTUPLE.ROOT/00185836/0000/00185836_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2017", "MagDown"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2018/DVNTUPLE.ROOT/00185830/0000/00185830_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2018", "MagUp"}},
            //{"root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2018/DVNTUPLE.ROOT/00185832/0000/00185832_00000001_1.dvntuple.root", {"MC_B2KpiJpsi", "2018", "MagDown"}},
        {data2015magup, {type, year, polarity}},
    };

    unordered_map<string, string> treestripping = /*! Remember to add the path (i.e. the name of the stripping plus an underscore)*/
        {
            {"B2KpiX2PiPiPi_Tuple/DecayTree", "B2KpiX2PiPiPi"},
            {"B2KpiX2PiPiPiM_Tuple/DecayTree", "B2KpiX2PiPiPiM"},
        };

    // List of particles stored in the input tuple file
    unordered_map<string, string> particletuplenames = /*! The first is the name in the tuple, the second is the alias*/
        {
            {"B0", "b"},
            {"xpim", "xpim"},
            {"xpip", "xpip"},
            {"xpiz", "xpiz"},
            {"Kplus", "k"},
            {"piminus", "pi"},
            {"KS0", "x"}};

    unordered_map<string, string> variabletuplenames = /*! The first is the name in the tuple, the second is the alias*/
        {
            {"P", "p"},
            {"PT", "pt"},
            {"PX", "px"},
            {"PY", "py"},
            {"PZ", "pz"},
            {"M", "m"},
            {"PE", "e"},
            {"ENDVERTEX_X", "x"},
            {"ENDVERTEX_Y", "y"},
            {"ENDVERTEX_Z", "z"},
            {"ENDVERTEX_CHI2", "vtx_chi2"},
            {"ENDVERTEX_NDOF", "vtx_ndof"},
            {"IP_TOPPV", "ip"},
            {"IPCHI2_TOPPV", "ip_chi2"},
            {"FD_TOPPV", "fd"},
            {"FDCHI2_TOPPV", "fd_chi2"},
            {"TAU", "tau"},
            {"CTAU", "ctau"},
            {"TAUERR", "tau_err"},
            {"TAUCHI2", "tau_chi2"},
            {"L0HadronDecision_Dec", "l0_hadron_dec"},
            {"L0HadronDecision_TIS", "l0_hadron_tis"},
            {"L0HadronDecision_TOS", "l0_hadron_tos"},
            {"DIRA_OWNPV", "dira"},
            {"MAXDOCA", "maxdoca"},
            {"MINDOCA", "mindoca"},
            {"DOCA12", "doca_kpi"},
            {"IP_OWNPV", "ip_best"},
            {"IPCHI2_OWNPV", "ip_chi2_best"},
            {"FD_OWNPV", "fd_best"},
            {"FDCHI2_OWNPV", "fd_chi2_best"},
            {"IP_ORIVX", "ip_bvtx"},
            {"IPCHI2_ORIVX", "ip_chi2_bvtx"},
            {"FD_ORIVX", "fd_bvtx"},
            {"FDCHI2_ORIVX", "fd_chi2_bvtx"},
            {"TRUEID", "truepid"},
            {"ID", "pid"},
            {"ProbNNe", "pnn_e"},
            {"ProbNNp", "pnn_p"},
            {"ProbNNk", "pnn_k"},
            {"ProbNNpi", "pnn_pi"},
            {"ProbNNmu", "pnn_mu"},
            {"ProbNNghost", "pnn_ghost"},
            {"TRACK_GhostProb", "prb_ghost"},
            {"CL", "cl"},
            {"MC_MOTHER_ID", "MC_mother_ID"},
            {"MC_GD_MOTHER_ID", "MC_GD_mother_ID"},
            {"MC_GD_GD_MOTHER_ID", "MC_GD_GD_mother_ID"},
            {"status", "fitstatus"}};

    unordered_map<string, string> particleconstrained =
        {
            {"xpiz", "xpiz"},
            {"xpizpimisidask", "xpiz_pimisidask"},
            {"jpsixpiz", "Jpsixpiz"},
            {"jpsixpizpimisidask", "Jpsixpiz_pimisidask"},
            {"omegaxpiz", "omegaxpiz"},
            {"omegaxpizpimisidask", "omegaxpiz_pimisidask"}};

    unordered_map<string, string> firststepdaughter =
        {
            {"omega_782", "omega"},
            {"J_psi_1S", "jpsi"},
            {"Kplus", "k"},
            {"piplus", "pi"},
            {"pi0", "xpiz"},
            {"KS0", "x"}};

    unordered_map<string, string> secondstepdaughter =
        {
            {"piplus_0", "xpim"},
            {"piplus", "xpip"}};

    // Here the RDataFrame is built using the tree from the input tuple file
    for (auto &&tuplefilenameit : fileyearpolarity)
    {
        string outputfilename = get<0>(tuplefilenameit.second) + "_" + get<1>(tuplefilenameit.second) + "_" + get<2>(tuplefilenameit.second);
        TChain datasetyearpolaritystripping; 
        for (auto &&tupletreenameit : treestripping)
        {
            // Defining the quantities that will be used
            ROOT::RDataFrame *dataset = NULL;
            vector<string> vars_tobe_written;

            //const char *tuplefilename = tuplefilenameit.first.data();
            const char *tupletreename = tupletreenameit.first.data();

            // Names of the file and tree used for the output
            char *outputtreename = "DecayTree";

            // Defining the chain
            TChain *inchain = new TChain(tupletreename);
            for (auto &&indatastring : *(tuplefilenameit.first))
                inchain->Add(indatastring.data());

            // Checking if the chain does exist
            bool chainzombie = inchain->IsZombie();

            // Checking if the file and the tree do exist
            //TFile* file = TFile::Open(tuplefilename);
            //cout << file->GetName() << endl;
            //TTree* tree = file->Get<TTree>(tupletreename);
            //bool filezombie = file->IsZombie();
            //bool treezombie = tree->IsZombie();

            // If the file and the tree do exist, extract a RDataFrame
            if (chainzombie == false)
            {
                dataset = new ROOT::RDataFrame(*inchain);
            }

            auto latestDF = std::make_unique<ROOT::RDF::RNode>(*dataset);
            for (auto &&particletuplenameit : particletuplenames)
            {
                for (auto &&variabletuplenameit : variabletuplenames)
                {
                    string invartupname = particletuplenameit.first + "_" + variabletuplenameit.first;
                    string outvartupname = particletuplenameit.second + "_" + variabletuplenameit.second;
                    if (dataset->HasColumn(invartupname) == true)
                    {
                        latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvartupname, invartupname.data()));
                        vars_tobe_written.push_back(outvartupname);
                        //cout << invartupname << endl
                        //     << outvartupname << endl
                        //     << endl;
                    }
                }
            }

            for (auto &&particleconstrainedit : particleconstrained)
            {
                string in_head = "B0_CONS_" + particleconstrainedit.first;
                for (auto &&firststepdaughterit : firststepdaughter)
                {
                    string in_first = in_head + "_" + firststepdaughterit.first;
                    for (auto &&secondstepdaughterit : secondstepdaughter)
                    {
                        string in_second = in_first + "_" + secondstepdaughterit.first;
                        for (auto &&variabletuplenameit : variabletuplenames)
                        {
                            string invar_head = in_head + "_" + variabletuplenameit.first;
                            string outvar_head = "b_cons_" + particleconstrainedit.second + "_" + variabletuplenameit.second;

                            string invar_first = in_first + "_" + variabletuplenameit.first;
                            string outvar_first = firststepdaughterit.second + "_cons_" + particleconstrainedit.second + "_" + variabletuplenameit.second;

                            string invar_second = in_second + "_" + variabletuplenameit.first;
                            string outvar_second = secondstepdaughterit.second + "_cons_" + particleconstrainedit.second + "_" + variabletuplenameit.second;

                            if (dataset->HasColumn(invar_head) == true && latestDF->HasColumn(outvar_head) == false)
                            {
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_head, invar_head.data()));
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_head + "_best", (invar_head + "[0]").data()));
                                vars_tobe_written.push_back(outvar_head);
                                vars_tobe_written.push_back((outvar_head + "_best").data());
                                //cout << invar_head << endl
                                //     << outvar_head << endl
                                //     << endl;
                            }
                            else if (dataset->HasColumn(invar_first) == true && latestDF->HasColumn(outvar_first) == false)
                            {
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_first, invar_first.data()));
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_first + "_best", (invar_first + "[0]").data()));
                                vars_tobe_written.push_back(outvar_first);
                                vars_tobe_written.push_back((outvar_first + "_best").data());
                                //cout << invar_first << endl
                                //     << outvar_first << endl
                                //     << endl;
                            }
                            else if (dataset->HasColumn(invar_second) == true && latestDF->HasColumn(outvar_second) == false)
                            {
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_second, invar_second.data()));
                                latestDF = std::make_unique<ROOT::RDF::RNode>(latestDF->Define(outvar_second + "_best", (invar_second + "[0]").data()));
                                vars_tobe_written.push_back(outvar_second);
                                vars_tobe_written.push_back((outvar_second + "_best").data());
                                //cout << invar_second << endl
                                //     << outvar_second << endl
                                //     << endl;
                            }
                        }
                    }
                }
            }

            latestDF = make_unique<ROOT::RDF::RNode>(latestDF->Define("Year", [tuplefilenameit]() -> int
                                                                      { return stoi(get<1>(tuplefilenameit.second)); }));
            vars_tobe_written.push_back("Year");
            latestDF = make_unique<ROOT::RDF::RNode>(latestDF->Define("Stripping", [tupletreenameit]() -> string
                                                                      { return tupletreenameit.second; }));
            latestDF = make_unique<ROOT::RDF::RNode>(latestDF->Define("Pi0Merged", [](string strip) -> int
                                                                      { if (strip.compare("B2KpiX2PiPiPi")==0) return 0; else return 1;}, {"Stripping"}));
            vars_tobe_written.push_back("Stripping");
            vars_tobe_written.push_back("Pi0Merged");
            vars_tobe_written.push_back("nCandidate");
            vars_tobe_written.push_back("eventNumber");
            vars_tobe_written.push_back("Polarity");

            // Modifying write options of the root data frame to overwrite other trees
            ROOT::RDF::RSnapshotOptions snapopt;
            snapopt.fMode = "RECREATE";
            snapopt.fOverwriteIfExists = "TRUE";

            // Write the new dataset to a new file. The tree name is the original one
            // Check if an output file name is specified, otherwise "_CREATED" is appended to the original one
            /*if (((string)outputfilename).compare("") == 0)
                dataset->Snapshot(outputtreename, (string(outputfilename)+ "_CREATED.root").data(), dataset->GetColumnNames(), snapopt);
            else*/
            latestDF->Filter(cuts).Snapshot(outputtreename, ("root://eosuser.cern.ch/eos/user/p/paandreo/RUN2DATA_QEE_BTOKSTARX_19MAY2023/"+outputfilename+"_"+tupletreenameit.second+ ".root").data(), vars_tobe_written, snapopt);
            
            delete dataset;
/*            TFile *input_signal_file = TFile::Open((outputfilename+"_"+tupletreenameit.second+ ".root").data(), "read");
            treetoberead = (TTree *)input_signal_file->Get(outputtreename);

            TIter nextbranch(treetoberead->GetListOfBranches());
            vector<string> *particlenames = new vector<string>;
            particlenames->push_back("b");
            particlenames->push_back("xpim");
            particlenames->push_back("xpip");
            particlenames->push_back("xpiz");
            particlenames->push_back("k");
            particlenames->push_back("pi");
            particlenames->push_back("x");
            TBranch *branch;
            int lengthbranchnames = 0, lengthpurgedbranchnames = 0;

            input_signal_file = TFile::Open("/home/pasqualeandreola/DarkALP/InputFiles/kpi3pi_2015_md.root", "read");
            treetoberead = (TTree *)input_signal_file->Get("T");

            nextbranch = treetoberead->GetListOfBranches();
            while ((branch = (TBranch *)nextbranch()))
            {
                string branchname = (string)branch->GetName();
                for (vector<string>::iterator partnameit = particlenames->begin(); partnameit != particlenames->end(); partnameit++)
                {
                    bool variablechecker = branchname.find(*partnameit + "_p") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_e") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_q") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_x") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_y") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_z") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_m") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_vtx") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_ip") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_fd") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_tau") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_l0_hadron_") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_dira") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_doca") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_hcal") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_type") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_is_mu") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_vid") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_cl") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_ctau") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_chi2") != string::npos;
                    variablechecker = variablechecker || branchname.find(*partnameit + "_ndof") != string::npos;
                    variablechecker = variablechecker && branchname.find("_electron_") == string::npos;
                    variablechecker = variablechecker && branchname.find("_muon_") == string::npos;
                    variablechecker = variablechecker && branchname.find("_dimuon_") == string::npos;
                    variablechecker = variablechecker && branchname.find("_photon_") == string::npos;
                    if (variablechecker == true)
                    {
                        TClass *branchclass = new TClass();
                        EDataType *branchetype = new EDataType();
                        branch->GetExpectedType(branchclass, *branchetype);
                        if (std::find(vars_tobe_written.begin(), vars_tobe_written.end(), branchname) == vars_tobe_written.end())
                        {
                            cout << branchname << endl;
                            cout << TDataType::GetDataType(*branchetype)->GetName() << endl
                                 << endl;
                        }
                    }
                }
                bool variablechecker = branchname.find("iCand") != string::npos;
                // variablechecker = variablechecker || branchname.find("dtf_")!=string::npos;
                if (variablechecker == true)
                {
                    TClass *branchclass = new TClass();
                    EDataType *branchetype = new EDataType();
                    branch->GetExpectedType(branchclass, *branchetype);
                    if (std::find(vars_tobe_written.begin(), vars_tobe_written.end(), branchname) == vars_tobe_written.end())
                    {
                        cout << branchname << endl;
                        cout << TDataType::GetDataType(*branchetype)->GetName() << endl
                             << endl;
                    }
                }
            }

            //datasetyearpolaritystripping->Add((outputfilename+"_"+tupletreenameit.second+ ".root/"+outputtreename).data()); 

        }*/}
        //ROOT::RDataFrame dataframeyearpolaritystripping(datasetyearpolaritystripping);
        //dataframeyearpolaritystripping.Snapshot("DecayTree", (outputfilename+".root").data());
    }
    return 0;
}