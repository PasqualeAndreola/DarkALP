/*!
 * \mainpage R-J-Psi
 *
 * \section intro_sec Introduction
 *
 * This project has been developed as a Master Thesis project used to search Lepton flavour
 * universality violation in CMS.
 */

/*!
 * \file main.cpp
 * \brief Main file
 */

#include "HeaderFiles/AnalysisVariable.h"
#include "HeaderFiles/EfficienciesAnalyzer.h"
#include "HeaderFiles/EfficienciesRatioAnalyzer.h"
#include "HeaderFiles/Fitter.h"
#include "HeaderFiles/MVACutOptimizerFactory.h"
#include "HeaderFiles/MVACutOptimizerReader.h"
#include "HeaderFiles/PrintFuncInfo.h"
#include "HeaderFiles/RootFileCreator.h"
#include "HeaderFiles/TMVAMethod.h"
#include <filesystem>
#include <fstream>
#include <ostream>
#include <RooFit.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooDataSet.h>
#include <RooChi2Var.h>
#include <RooChebychev.h>
#include <RooArgusBG.h>
#include <RooFFTConvPdf.h>
#include <RooAddPdf.h>
#include <TChain.h>
#include <TArrow.h>
#include <TPad.h>
#include <TF1.h>
#include <TMVA/Config.h>
#include <TMVA/Configurable.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Factory.h>
#include "TMVA/Reader.h"
#include <TMVA/TMVAGui.h>
#include <TMVA/RReader.hxx>
#include "Blue.h"
#include <TMatrixD.h>
#include <chrono>
#include <ctime>

// It can be useful to use these namespaces
using namespace std;
using namespace RooFit;

int main(int argc, char *argv[])
{
  gErrorIgnoreLevel = kWarning;
  int colors_palette[3];
  colors_palette[0] = kBlue;
  colors_palette[1] = kGreen + 1;
  colors_palette[2] = kRed;
  gStyle->SetPalette(3, colors_palette);
  bool cuts_as_analysis_variables = true;
  bool files_with_stripping_cuts = true;

  auto start = chrono::system_clock::now();
  auto end = chrono::system_clock::now();
  time_t start_time = chrono::system_clock::to_time_t(start);
  time_t end_time = chrono::system_clock::to_time_t(end);

  cout << endl
       << "Program compiled correctly." << endl
       << "Starting date: " << ctime(&start_time) << endl;

  unordered_map<string, pair<string, string>> data_holder;

  if (files_with_stripping_cuts == true)
  {
    data_holder["KPi3Pi_SIM_100ps"] = pair("InputFiles/KPi3Pi_SIM_100ps.root", "T");
    data_holder["KPi3Pi_SIM_1ps"] = pair("InputFiles/KPi3Pi_SIM_1ps.root", "T");
    data_holder["KPi3Pi_SIM_Omega"] = pair("InputFiles/KPi3Pi_SIM_Omega.root", "T");
    data_holder["KPiEtaPiPi_SIM_100ps"] = pair("InputFiles/KPiEtaPiPi_SIM_100ps.root", "T");
    data_holder["KPiEtaPiPi_SIM_1ps"] = pair("InputFiles/KPiEtaPiPi_SIM_1ps.root", "T");
    data_holder["KPiEtaPiPi_SIM_Etapr"] = pair("InputFiles/KPiEtaPiPi_SIM_Etapr.root", "T");
    data_holder["KPi3Pi_2015_md_minCuts"] = pair("InputFiles/KPi3Pi_2015_md_minCuts.root", "T");
    data_holder["KPiEtaPiPi_2015_md_minCuts"] = pair("InputFiles/KPiEtaPiPi_2015_md_minCuts.root", "T");
  }
  else
  {
    data_holder["KPi3Pi_SIM_100ps"] = pair("InputFiles/kpi3pi-sim-100ps.root", "T");
    data_holder["KPi3Pi_SIM_1ps"] = pair("InputFiles/kpi3pi-sim-1ps.root", "T");
    data_holder["KPi3Pi_SIM_Omega"] = pair("InputFiles/kpi3pi-sim-omega.root", "T");
    data_holder["KPiEtaPiPi_SIM_100ps"] = pair("InputFiles/kpietapipi-sim-100ps.root", "T");
    data_holder["KPiEtaPiPi_SIM_1ps"] = pair("InputFiles/kpietapipi-sim-1ps.root", "T");
    data_holder["KPiEtaPiPi_SIM_Etapr"] = pair("InputFiles/kpietapipi-sim-etapr.root", "T");
    data_holder["KPi3Pi_2015_md_minCuts"] = pair("InputFiles/kpi3pi_2015_md_minCuts.root", "T");
    data_holder["KPiEtaPiPi_2015_md_minCuts"] = pair("InputFiles/kpietapipi_2015_md_minCuts.root", "T");
  }

  /*MISSING X_DIRA AND XETA_CL AND XPIZ_CL*/
  string stripping_cuts = "&& xpip_prb_ghost<0.3 && xpip_p>3000 && xpip_ip_chi2_best>36 && xpip_pt>250 && xpip_pnn_pi>0.2"
                          "&& xpim_prb_ghost<0.3 && xpim_p>3000 && xpim_ip_chi2_best>36 && xpim_pt>250 && xpim_pnn_pi>0.2"
                          "&& x_m<5000 && x_pt>2000 && x_maxdoca<0.2 && x_fd_chi2>25 && x_vtx_chi2/x_vtx_ndof<10"
                          "&& k_pnn_k>0.1 && k_p>2000 && k_ip_chi2_best>9 && k_pt>250"
                          "&& pi_pnn_pi>0.2 && pi_p>2000 && pi_ip_chi2_best>9 && pi_pt>250"
                          "&& (b_m>4800 || b_m<5800) && b_pt>3000 && b_vtx_chi2/b_vtx_ndof<15 && b_tau>0.0002 && b_ip_chi2<10 && b_dira>0";
  string stripping_cuts_xeta = stripping_cuts + "&& xeta_p>2000 && (xeta_m>450||xeta_m<650) && xeta_pt>500";
  string stripping_cuts_xpiz = stripping_cuts + "&& xpiz_pt>500 ";

  string mincuts, cuts, cutspretty;
  unordered_map<string, string> mapofcuts;

  vector<AnalysisVariable> *vec_analysisvariables = new vector<AnalysisVariable>;
  AnalysisVariable B0_true_energy;
  B0_true_energy.variable_name = "b_true_e";
  B0_true_energy.variable_expression = "b_true_e";
  B0_true_energy.variable_prettyname = "E_{B^{0}}";
  B0_true_energy.variable_dimension = "GeV";
  B0_true_energy.SetVarCutOpt(100, 1e6);
  B0_true_energy.SetVarBinningOpt(15, 0, 1e6);
  B0_true_energy.variable_plot_flag = true;
  // vec_analysisvariables->push_back(B0_true_energy);

  AnalysisVariable B0_true_p;
  B0_true_p.variable_name = "b_true_p";
  B0_true_p.variable_expression = "b_true_p";
  B0_true_p.variable_prettyname = "p_{B^{0}}";
  B0_true_p.variable_dimension = "GeV";
  B0_true_p.SetVarCutOpt(100, 1e6);
  B0_true_p.SetVarBinningOpt(15, 0, 1e6);
  B0_true_p.variable_plot_flag = true;
  // vec_analysisvariables->push_back(B0_true_p);

  AnalysisVariable B0_mass;
  B0_mass.variable_name = "b_m";
  B0_mass.variable_expression = "b_m";
  B0_mass.variable_mvacutfeature = true;
  B0_mass.SetVarCutOpt(4800, 5800);
  vec_analysisvariables->push_back(B0_mass);

  AnalysisVariable B0_masscorr;
  B0_masscorr.variable_name = "b_mcor";
  B0_masscorr.variable_expression = "b_mcor";
  B0_masscorr.SetVarCutOpt(4000, 10000);
  vec_analysisvariables->push_back(B0_masscorr);

  AnalysisVariable B0_pt;
  B0_pt.variable_name = "b_pt";
  B0_pt.variable_expression = "b_pt";
  B0_pt.variable_mvacutfeature = true;
  B0_pt.SetVarCutOptMin(3000);
  vec_analysisvariables->push_back(B0_pt);

  AnalysisVariable B0_dira;
  B0_dira.SetVarNames("b_dira", "b_dira", "b_dira", "");
  B0_dira.SetVarCutOptMin(0.9995);
  B0_dira.variable_cutchecker = true;
  B0_dira.variable_mvacutfeature = true;
  vec_analysisvariables->push_back(B0_dira);

  AnalysisVariable B0_vtx_quality;
  B0_vtx_quality.SetVarNames("b_vtx_quality", "b_vtx_quality", "(b_vtx_chi2)/(b_vtx_ndof)");
  B0_vtx_quality.SetVarCutOptMax(6);
  B0_vtx_quality.variable_cutchecker = true;
  vec_analysisvariables->push_back(B0_vtx_quality);

  AnalysisVariable B0_doca_kpi;
  B0_doca_kpi.SetVarNames("b_doca_kpi", "b_doca_kpi", "b_doca_kpi");
  B0_doca_kpi.SetVarCutOptMax(0.2);
  B0_doca_kpi.variable_mvacutfeature = true;
  B0_doca_kpi.variable_cutchecker = true;
  vec_analysisvariables->push_back(B0_doca_kpi);

  AnalysisVariable B0_fd;
  B0_fd.SetVarNames("b_fd", "", "b_fd");
  B0_fd.SetVarCutOptMin(3);
  B0_fd.variable_cutchecker = true;
  B0_fd.variable_mvacutfeature = true;
  vec_analysisvariables->push_back(B0_fd);

  AnalysisVariable K_energy;
  K_energy.variable_name = "k_e";
  K_energy.variable_expression = "k_e";
  K_energy.variable_mvacutfeature = true;
  K_energy.SetVarCutOpt(4800, 5800);
  // vec_analysisvariables->push_back(K_energy);

  AnalysisVariable K_pt;
  K_pt.variable_name = "k_pt";
  K_pt.variable_expression = "k_pt";
  K_pt.variable_mvacutfeature = true;
  K_pt.SetVarCutOptMin(3000);
  // vec_analysisvariables->push_back(K_pt);

  AnalysisVariable Pi_energy;
  Pi_energy.variable_name = "pi_e";
  Pi_energy.variable_expression = "pi_e";
  Pi_energy.variable_mvacutfeature = true;
  Pi_energy.SetVarCutOpt(4800, 5800);
  // vec_analysisvariables->push_back(Pi_energy);

  AnalysisVariable Pi_pt;
  Pi_pt.variable_name = "pi_pt";
  Pi_pt.variable_expression = "pi_pt";
  Pi_pt.variable_mvacutfeature = true;
  Pi_pt.SetVarCutOptMin(3000);
  // vec_analysisvariables->push_back(Pi_pt);

  AnalysisVariable XPip_mass;
  XPip_mass.variable_name = "xpip_m";
  XPip_mass.variable_expression = "xpip_m";
  XPip_mass.variable_mvacutfeature = true;
  XPip_mass.SetVarCutOpt(4800, 5800);
  // vec_analysisvariables->push_back(XPip_mass);

  AnalysisVariable XPip_pt;
  XPip_pt.variable_name = "xpip_pt";
  XPip_pt.variable_expression = "xpip_pt";
  XPip_pt.variable_mvacutfeature = true;
  XPip_pt.SetVarCutOptMin(3000);
  // vec_analysisvariables->push_back(XPip_pt);

  AnalysisVariable XPim_mass;
  XPim_mass.variable_name = "xpim_m";
  XPim_mass.variable_expression = "xpim_m";
  XPim_mass.variable_mvacutfeature = true;
  XPim_mass.SetVarCutOpt(4800, 5800);
  // vec_analysisvariables->push_back(XPim_mass);

  AnalysisVariable XPim_pt;
  XPim_pt.variable_name = "xpim_pt";
  XPim_pt.variable_expression = "xpim_pt";
  XPim_pt.variable_mvacutfeature = true;
  XPim_pt.SetVarCutOptMin(3000);
  // vec_analysisvariables->push_back(XPim_pt);

  AnalysisVariable X_mass;
  X_mass.variable_name = "x_m";
  X_mass.variable_expression = "x_m";
  X_mass.variable_mvacutfeature = true;
  X_mass.SetVarCutOpt(4800, 5800);
  // vec_analysisvariables->push_back(X_mass);

  AnalysisVariable X_pt;
  X_pt.variable_name = "x_pt";
  X_pt.variable_expression = "x_pt";
  X_pt.variable_mvacutfeature = true;
  X_pt.SetVarCutOptMin(3000);
  // vec_analysisvariables->push_back(X_pt);

  AnalysisVariable X_dira;
  X_dira.SetVarNames("x_dira", "x_dira", "x_dira", "");
  X_dira.SetVarCutOptMin(0.9995);
  X_dira.variable_cutchecker = true;
  X_dira.variable_mvacutfeature = true;
  // vec_analysisvariables->push_back(X_dira);

  AnalysisVariable X_fd;
  X_fd.SetVarNames("x_fd", "", "x_fd");
  X_fd.SetVarCutOptMin(3);
  X_fd.variable_cutchecker = true;
  X_fd.variable_mvacutfeature = true;
  // vec_analysisvariables->push_back(X_fd);

  AnalysisVariable B0_fd_chi2;
  B0_fd_chi2.SetVarNames("b_fd_chi2", "b_fd_chi2", "b_fd_chi2");
  B0_fd_chi2.SetVarCutOptMin(500);
  B0_fd_chi2.variable_cutchecker = true;
  // vec_analysisvariables->push_back(B0_fd_chi2);

  AnalysisVariable X_fd_chi2;
  X_fd_chi2.SetVarNames("x_fd_chi2", "x_fd_chi2", "x_fd_chi2");
  X_fd_chi2.SetVarCutOptMin(500);
  X_fd_chi2.variable_cutchecker = true;
  // vec_analysisvariables->push_back(X_fd_chi2);

  AnalysisVariable bx_fd_chi2;
  bx_fd_chi2.SetVarNames("b_plus_x_fd_chi2", "", "b_fd_chi2+x_fd_chi2");
  bx_fd_chi2.SetVarCutOptMin(500);
  bx_fd_chi2.variable_cutchecker = true;
  vec_analysisvariables->push_back(bx_fd_chi2);

  AnalysisVariable X_vtx_quality;
  X_vtx_quality.SetVarNames("x_vtx_quality", "x_vtx_quality", "(x_vtx_chi2)/(x_vtx_ndof)");
  X_vtx_quality.SetVarCutOptMax(5);
  X_vtx_quality.variable_cutchecker = true;
  vec_analysisvariables->push_back(X_vtx_quality);

  AnalysisVariable XforwardofB0_or_XpromptfromB;
  XforwardofB0_or_XpromptfromB.variable_stringcut = "((x_z-b_z > -20) ||"
                                                    "( (b_maxdoca > 0.2) && (x_z-b_z > 0) ) ||"
                                                    "( (b_maxdoca < 1) && (x_z-b_z >20) ))";
  XforwardofB0_or_XpromptfromB.variable_cutchecker = true;
  vec_analysisvariables->push_back(XforwardofB0_or_XpromptfromB);

  AnalysisVariable kpixpipxpim_pt;
  kpixpipxpim_pt.variable_stringcut = "((k_pt > 100 )||"
                                      "(pi_pt > 100)||"
                                      "(xpip_pt > 100)||"
                                      "(xpim_pt > 100))";
  vec_analysisvariables->push_back(kpixpipxpim_pt);

  AnalysisVariable kpixpipxpim_p;
  kpixpipxpim_p.variable_stringcut = "((k_p > 1500 )||"
                                     "(pi_p > 1500)||"
                                     "(xpip_p > 1500)||"
                                     "(xpim_p > 1500))";
  vec_analysisvariables->push_back(kpixpipxpim_p);

  AnalysisVariable kpixpipxpim_ip_chi2_best;
  kpixpipxpim_ip_chi2_best.variable_stringcut = "(k_ip_chi2_best>10 || pi_ip_chi2_best>10 || xpip_ip_chi2_best>10 || xpim_ip_chi2_best>10)";
  kpixpipxpim_ip_chi2_best.variable_cutchecker = true;
  vec_analysisvariables->push_back(kpixpipxpim_ip_chi2_best);

  AnalysisVariable kpi_ip_chi2_bvtx;
  kpi_ip_chi2_bvtx.variable_stringcut = "(k_ip_chi2_bvtx < 10 || pi_ip_chi2_bvtx < 10)";
  kpi_ip_chi2_bvtx.variable_cutchecker = true;
  vec_analysisvariables->push_back(kpi_ip_chi2_bvtx);

  AnalysisVariable kpipxpipxpim_prb_ghost;
  kpipxpipxpim_prb_ghost.variable_stringcut = "(k_prb_ghost < 0.1 || pi_prb_ghost < 0.1 || xpip_prb_ghost < 0.1 || xpim_prb_ghost < 0.1)";
  kpipxpipxpim_prb_ghost.variable_cutchecker = true;
  vec_analysisvariables->push_back(kpipxpipxpim_prb_ghost);

  AnalysisVariable k_pnn_k;
  k_pnn_k.SetVarNames("k_pnn_k", "", "k_pnn_k");
  k_pnn_k.SetVarCutOptMin(0.3);
  k_pnn_k.variable_cutchecker = true;
  vec_analysisvariables->push_back(k_pnn_k);

  AnalysisVariable kpixpipxpim_sum_ip_chi2_best;
  kpixpipxpim_sum_ip_chi2_best.variable_stringcut = "(k_ip_chi2_best+pi_ip_chi2_best+xpip_ip_chi2_best+xpim_ip_chi2_best > 200.)";
  vec_analysisvariables->push_back(kpixpipxpim_sum_ip_chi2_best);

  vector<TMVAMethod> tmvamethods;
  {
    TMVAMethod kmlp1;
    kmlp1.tmvamethodname = "kmlp1";
    TString config = "CreateMVAPdfs:VarTransform=N:NCycles=100:HiddenLayers=N+2,N+1,N:NeuronType=sigmoid:EstimatorType=CE:";
    config += "TrainingMethod=BP:LearningRate=0.1:DecayRate=0.05:ConvergenceTests=50";
    kmlp1.tmvamethodconfiguration = config;
    kmlp1.tmvamethodtype = TMVA::Types::EMVA::kMLP;
    tmvamethods.push_back(kmlp1);
  }
  {
    TMVAMethod kmlp2;
    kmlp2.tmvamethodname = "kmlp2";
    TString config = "CreateMVAPdfs:VarTransform=N:NCycles=300:HiddenLayers=N+1:NeuronType=sigmoid:EstimatorType=CE:";
    config += "TrainingMethod=BP:LearningRate=0.05:DecayRate=0.05:CalculateErrors=True:ConvergenceTests=10";
    kmlp2.tmvamethodconfiguration = config;
    kmlp2.tmvamethodtype = TMVA::Types::EMVA::kMLP;
    tmvamethods.push_back(kmlp2);
  }
  {
    TMVAMethod kdl1;
    kdl1.tmvamethodname = "kdl1";
    TString config = "CreateMVAPdfs:!H:V";
    config += ":VarTransform=N";
    config += ":ErrorStrategy=CROSSENTROPY";
    config += ":WeightInitialization=XAVIERUNIFORM";
    config += ":Layout=TANH|100, TANH|50, TANH|10, LINEAR";
    config += ":TrainingStrategy=LearningRate=1e-2,Momentum=0.5, Repetitions=1,ConvergenceSteps=100,BatchSize=100,DropConfig=0.0+0.5+0.5+0.0";
    config += ",WeightDecay=0.001,Regularization=L2,TestRepetitions=15,Multithreading=True";
    kdl1.tmvamethodconfiguration = config;
    kdl1.tmvamethodtype = TMVA::Types::EMVA::kDL;
    tmvamethods.push_back(kdl1);
  }
  {
    TMVAMethod kbdt1;
    kbdt1.tmvamethodname = "kbdt1";
    TString config = "CreateMVAPdfs:Ntrees=100:MaxDepth=5:MinNodeSize=5%:nCuts=-1:BoostType=RealAdaBoost:UseRandomisedTrees=True";
    config += ":UseNvars=3:UsePoissonNvars=True";
    kbdt1.tmvamethodconfiguration = config;
    kbdt1.tmvamethodtype = TMVA::Types::EMVA::kBDT;
    tmvamethods.push_back(kbdt1);
  }

  vector<string> features;
  /*for (vector<AnalysisVariable>::iterator analvarit = vec_analysisvariables->begin(); analvarit!= vec_analysisvariables->end(); analvarit++)
  {
    if (analvarit->variable_mvacutfeature == true)
    {
        features.push_back(analvarit->variable_name);
    }
  }*/

  for (vector<AnalysisVariable>::iterator analvarcutiterator = vec_analysisvariables->begin(); analvarcutiterator != vec_analysisvariables->end(); analvarcutiterator++)
  {
    if (analvarcutiterator != vec_analysisvariables->begin())
      cuts.append("&&");
    if ((analvarcutiterator != vec_analysisvariables->begin()) && (analvarcutiterator - 1)->variable_cutchecker == false)
      mincuts.append("&&");
    if (analvarcutiterator->variable_stringcut.compare("") != 0)
    {
      cuts.append(analvarcutiterator->variable_stringcut);
      if (analvarcutiterator->variable_cutchecker == false)
      {
        mincuts.append(analvarcutiterator->variable_stringcut);
      }
    }
    else if ((analvarcutiterator->variable_cut_pidflag == true) && (analvarcutiterator->variable_cut_pid != 0))
    {
      cuts.append(analvarcutiterator->variable_cut);
      if (analvarcutiterator->variable_cutchecker == false)
      {
        mincuts.append(analvarcutiterator->variable_cut);
      }
    }
    else if ((analvarcutiterator->variable_cut_pidflag == false) && (analvarcutiterator->variable_cut_pid != 0))
    {
      cuts.append(analvarcutiterator->variable_cut);
      if (analvarcutiterator->variable_cutchecker == false)
      {
        mincuts.append(analvarcutiterator->variable_cut);
      }
    }
    else
    {
      cuts.append(analvarcutiterator->variable_cut);
      if (analvarcutiterator->variable_cutchecker == false)
      {
        mincuts.append(analvarcutiterator->variable_cut);
      }
    }
  }

  mapofcuts["MinCuts"] = mincuts;
  mapofcuts["Cuts"] = cuts;
  // EfficienciesAnalyzer(data_holder, mapofcuts);

  ElapsedTimeStamper(start);
  /*
  for(unordered_map<string, pair<string, string>>::iterator dataholder_iterator = data_holder.begin(); dataholder_iterator != data_holder.end(); dataholder_iterator++)
  {
    if (dataholder_iterator->first.find("KPiEtaPiPi")!=string::npos)
    {
      RootFileCreatorFilterer(dataholder_iterator->second, dataholder_iterator->first.data(), mincuts+stripping_cuts_xeta);
    }
    if (dataholder_iterator->first.find("KPi3Pi")!=string::npos)
    {
      RootFileCreatorFilterer(dataholder_iterator->second, dataholder_iterator->first.data(), mincuts+stripping_cuts_xpiz);
    }
  }
  */
  /*
  pair<string, string> signal_kpi3pi_data_holder = pair("InputFiles/kpi3pi-sim-jpsi.root", "T");
  RootFileCreatorFilterer(signal_kpi3pi_data_holder, string("KPi3Pi_SIM_JPsi"), cuts+stripping_cuts_xpiz);
  pair<string, string> background_data_holder = pair("InputFiles/kpietapipi-sim-jpsi.root", "T");
  RootFileCreatorFilterer(background_data_holder, string("KPiEtaPiPi_SIM_JPsi"), cuts+stripping_cuts_xeta);
  */
  ofstream csv_efficiencies_ratios_file;

  if (cuts_as_analysis_variables == true)
  {
    AnalysisVariable analvar_mincuts;
    analvar_mincuts.variable_stringcut = mincuts;
    vec_analysisvariables->push_back(analvar_mincuts);
    AnalysisVariable analvar_cuts;
    analvar_cuts.variable_stringcut = cuts;
    vec_analysisvariables->push_back(analvar_cuts);
  }

  // EfficienciesRatioAnalyzer(vec_analysisvariables, csv_efficiencies_ratios_file, data_holder);

  vector<string> *fake_factory_names = new vector<string>; /*!< 0 = MC signal path, 1-2= in/out file, 3-4= pass/fail trees, 5-6 pass/fail closure test trees*/

  // Initializing variables to set up the reader
  unordered_map<string, pair<string, string>> kpi3pi_mvacut_factory_files, kpi3pi_mvacut_reader_files;
  kpi3pi_mvacut_factory_files["Signal"] = pair("InputFiles/KPi3Pi_Signal.root", "T");
  kpi3pi_mvacut_factory_files["Background"] = pair("InputFiles/KPi3Pi_Data_Sidebands.root", "T");
  kpi3pi_mvacut_factory_files["Output"] = pair("OutputFiles/kpi3pi_MVACutFactory.root", "T");
  // MVACutOptimizerFactory(kpi3pi_mvacut_factory_files, features, vec_analysisvariables, tmvamethods, start, true);
  // kpi3pi_mvacut_reader_files["KPi3PiJPsi_Reader"] = pair("InputFiles/KPi3Pi_SIM_JPsi.root", "T");
  kpi3pi_mvacut_reader_files["KPi3PiJPsi_Reader"] = pair("OutputFiles/KPi3Pi_SIM_JPsi.root", "T");
  // MVACutOptimizerReader(kpi3pi_mvacut_factory_files["Output"].first, kpi3pi_mvacut_reader_files, tmvamethods);

  B0_masscorr.SetVarBinningOpt(15, 5000, 6000);
  B0_masscorr.SetVarNames("b_mcor", "B^{0}_{m,cor}", "b_mcor", "MeV");
  //B0_masscorr.VariablePlotter(&(kpi3pi_mvacut_reader_files["KPi3PiJPsi_Reader"]), &B0_masscorr, string("mvacut_kbdt1>0.5"));
  
  ROOT::RDataFrame dataframe_KPi3Pi_JPsi_SIM = ROOT::RDataFrame("T", "OutputFiles/KPi3Pi_SIM_JPsi.root");
  FitterMass(&dataframe_KPi3Pi_JPsi_SIM, &B0_masscorr);

  ElapsedTimeStamper(start);
  end = chrono::system_clock::now();
  end_time = chrono::system_clock::to_time_t(end);
  cout << "Ending date: " << ctime(&end_time) << endl;
  return 0;
}