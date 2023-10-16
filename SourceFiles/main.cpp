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
#include "HeaderFiles/AnalysisVariableLambdaConstructor.h"
#include "HeaderFiles/EfficienciesAnalyzer.h"
#include "HeaderFiles/EfficienciesRatioAnalyzer.h"
#include "HeaderFiles/Fitter.h"
#include "HeaderFiles/MVACutOptimizerFactory.h"
#include "HeaderFiles/MVACutOptimizerReader.h"
#include "HeaderFiles/PrintFuncInfo.h"
#include "HeaderFiles/RootFileCreator.h"
#include "HeaderFiles/TMVAMethod.h"
#include "HeaderFiles/TreeSkimmer.h"
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
#include <RooFormulaVar.h>
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
#include <TRatioPlot.h>
#include <chrono>
#include <ctime>

// It can be useful to use these namespaces
using namespace std;
using namespace RooFit;

int main(int argc, char *argv[])
{
/*   gErrorIgnoreLevel = kWarning;
  int colors_palette[6];
  colors_palette[0] = kBlue;
  colors_palette[1] = kGreen + 1;
  colors_palette[2] = kCyan;
  colors_palette[3] = kViolet;
  colors_palette[4] = kYellow;
  colors_palette[5] = kOrange;
  gStyle->SetPalette(6, colors_palette); */
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

  // Stripping line:
  // StrippingB2KpiX2PiPiPiDarkBosonLine
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


  unordered_map<string, vector<string>> var_tobe_defined;
  var_tobe_defined["b_m_jpsiconstrained"] = {
    "x_px", "x_py", "x_pz",
    "k_px", "k_py", "k_pz", "k_e",
    "pi_px", "pi_py", "pi_pz", "pi_e"
    };

  var_tobe_defined["b_m_omegaconstrained"] = {
    "x_px", "x_py", "x_pz",
    "k_px", "k_py", "k_pz", "k_e",
    "pi_px", "pi_py", "pi_pz", "pi_e"
    };

  var_tobe_defined["x_m_pi0constrained"] = {
    "xpiz_px", "xpiz_py", "xpiz_pz",
    "xpip_px", "xpip_py", "xpip_pz", "xpip_e",
    "xpim_px", "xpim_py", "xpim_pz", "xpim_e"
    };

  var_tobe_defined["b_m_jpsiconstrained_k_pi0misidentifiedask"] = {
    "x_px", "x_py", "x_pz",
    "k_px", "k_py", "k_pz", "k_e",
    "pi_px", "pi_py", "pi_pz", "pi_e"    
  };



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

  /*AnalysisVariable B0_masscorr;
  B0_masscorr.variable_name = "b_mcor";
  B0_masscorr.variable_expression = "b_mcor";
  B0_masscorr.SetVarCutOpt(4000, 10000);
  vec_analysisvariables->push_back(B0_masscorr);*/

  AnalysisVariable B0_PT;
  B0_PT.variable_name = "B0_PT";
  B0_PT.variable_expression = "B0_PT";
  B0_PT.variable_mvacutfeature = true;
  B0_PT.SetVarCutOptMin(700);
  vec_analysisvariables->push_back(B0_PT);

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
    kmlp1.tmvamethodcut = "kmlp1 > 0.5";
    TString config = "CreateMVAPdfs:VarTransform=N:NCycles=100:HiddenLayers=N+2,N+1,N:NeuronType=sigmoid:EstimatorType=CE:";
    config += "TrainingMethod=BP:LearningRate=0.1:DecayRate=0.05:ConvergenceTests=50";
    kmlp1.tmvamethodconfiguration = config;
    kmlp1.tmvamethodtype = TMVA::Types::EMVA::kMLP;
    tmvamethods.push_back(kmlp1);
  }
  {
    TMVAMethod kmlp2;
    kmlp2.tmvamethodname = "kmlp2";
    kmlp2.tmvamethodcut = "kmlp2 > 0.5";
    TString config = "CreateMVAPdfs:VarTransform=N:NCycles=300:HiddenLayers=N+1:NeuronType=sigmoid:EstimatorType=CE:";
    config += "TrainingMethod=BP:LearningRate=0.05:DecayRate=0.05:CalculateErrors=True:ConvergenceTests=10";
    kmlp2.tmvamethodconfiguration = config;
    kmlp2.tmvamethodtype = TMVA::Types::EMVA::kMLP;
    tmvamethods.push_back(kmlp2);
  }
  {
    TMVAMethod kdl1;
    kdl1.tmvamethodname = "kdl1";
    kdl1.tmvamethodcut = "kdl1 > 0.6";
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
    kbdt1.tmvamethodcut = "kbdt1 > 0.51";
    TString config = "CreateMVAPdfs:Ntrees=100:MaxDepth=5:MinNodeSize=5%:nCuts=-1:BoostType=RealAdaBoost:UseRandomisedTrees=True";
    config += ":UseNvars=3:UsePoissonNvars=True";
    kbdt1.tmvamethodconfiguration = config;
    kbdt1.tmvamethodtype = TMVA::Types::EMVA::kBDT;
    tmvamethods.push_back(kbdt1);
  }

  //vector<string> features;
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

  if (cuts_as_analysis_variables == true)
  {
    AnalysisVariable analvar_mincuts;
    analvar_mincuts.variable_stringcut = mincuts;
    vec_analysisvariables->push_back(analvar_mincuts);
    AnalysisVariable analvar_cuts;
    analvar_cuts.variable_stringcut = cuts;
    vec_analysisvariables->push_back(analvar_cuts);
  }

  ofstream csv_efficiencies_ratios_file;
  // EfficienciesRatioAnalyzer(vec_analysisvariables, csv_efficiencies_ratios_file, data_holder);

  // Initializing variables to set up the factory
  unordered_map<string, pair<string, string>> kpi3pi_mvacut_factory_files;
  kpi3pi_mvacut_factory_files["Signal"] = pair("InputFiles/KPi3Pi_SIM_Omega.root", "T");
  kpi3pi_mvacut_factory_files["Background"] = pair("InputFiles/KPi3Pi_Data_Sidebands.root", "T");
  kpi3pi_mvacut_factory_files["Output"] = pair("OutputFiles/KPi3Pi_PureOmega_MVA_Factory.root", "T");
  //MVACutOptimizerFactory(kpi3pi_mvacut_factory_files, features, vec_analysisvariables, tmvamethods, start, true);

  vector<string> features;
  unordered_map<string, pair<string, string>> kpi3pi_sig_mvacut_factory_files, kpi3pi_bkg_mvacut_factory_files;
  string mvacut_jpsiomegamixfactory_outputfile = "OutputFiles/KPi3Pi_OmegaJPsi_MVA_Factory.root";
  string mvacut_pureomegafactory_outputfile = "OutputFiles/KPi3Pi_PureOmega_MVA_Factory.root";
  unordered_map<string, pair<string, string>> kpi3pi_mvacut_reader_files;
  kpi3pi_sig_mvacut_factory_files["0.5"] = pair("omegaMC2016StrippingFiltered.root", "DecayTree");
  //kpi3pi_sig_mvacut_factory_files["OmegaREV"] = pair("OutputFiles/KPi3Pi_SIM_Omega_REV.root", "T");
  kpi3pi_sig_mvacut_factory_files["0.499"] = pair("JpsiMC2016StrippingFiltered.root", "DecayTree");
  //kpi3pi_sig_mvacut_factory_files["JPsiREV"] = pair("OutputFiles/KPi3Pi_SIM_JPsi_REV.root", "T");
  kpi3pi_bkg_mvacut_factory_files["1"] = pair("InputFiles/KPi3Pi_Data_Sidebands.root", "T");
  //MVACutOptimizerFactoryMultiSigBkg(kpi3pi_sig_mvacut_factory_files, kpi3pi_bkg_mvacut_factory_files, mvacut_jpsiomegamixfactory_outputfile, features, tmvamethods);
  //MVACutOptimizerReader(mvacut_jpsiomegamixfactory_outputfile, pair("JpsiMC2016StrippingFiltered_SIGTEST.root", "DecayTree"), tmvamethods, string("NonloSo.root"), true);
  TChain kpi3pi_fulldataset_treechain("T");
  bool kpi3pi_mc_jpsi = false,
       kpi3pi_mc_omega = false,
       kpi3pi_dataset_2015 = false, 
       kpi3pi_dataset_2016 = false,
       kpi3pi_dataset_2017 = false,
       kpi3pi_dataset_2018 = false;
  /*
  if (kpi3pi_mc_jpsi)
  {     
    pair<string, string> kpi3pi_sim_jpsi_filter = pair("OutputFiles/kpi3pi_sim_jpsi_filter.root", "T");
    pair<string, string> kpi3pi_sim_jpsi_rev = pair("OutputFiles/kpi3pi_sim_jpsi_REV.root", "T");
    pair<string, string> kpi3pi_sim_jpsi_mva_rev = pair("OutputFiles/kpi3pi_sim_jpsi_MVA_REV.root", "T");
    TChain kpi3pi_sim_jpsi_treechain("T");
    kpi3pi_sim_jpsi_treechain.Add("InputFiles/kpi3pi-sim-jpsi.root");

    RootFileCreatorFilterer(&kpi3pi_sim_jpsi_treechain, kpi3pi_sim_jpsi_filter.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_sim_jpsi_filter, kpi3pi_sim_jpsi_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_sim_jpsi_rev, tmvamethods, kpi3pi_sim_jpsi_mva_rev.first, true);
  }

  if (kpi3pi_mc_omega)
  {     
    pair<string, string> kpi3pi_sim_omega_filter = pair("OutputFiles/kpi3pi_sim_omega_filter.root", "T");
    pair<string, string> kpi3pi_sim_omega_rev = pair("OutputFiles/kpi3pi_sim_omega_REV.root", "T");
    pair<string, string> kpi3pi_sim_omega_mva_rev = pair("OutputFiles/kpi3pi_sim_omega_MVA_REV.root", "T");
    TChain kpi3pi_sim_omega_treechain("T");
    kpi3pi_sim_omega_treechain.Add("InputFiles/kpi3pi-sim-omega.root");

    RootFileCreatorFilterer(&kpi3pi_sim_omega_treechain, kpi3pi_sim_omega_filter.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_sim_omega_filter, kpi3pi_sim_omega_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_sim_omega_rev, tmvamethods, kpi3pi_sim_omega_mva_rev.first, true);
  }

  if (kpi3pi_dataset_2015)
  {     
    pair<string, string> kpi3pi_2015_mdu = pair("OutputFiles/kpi3pi_2015_mdu.root", "T");
    pair<string, string> kpi3pi_2015_mdu_rev = pair("OutputFiles/kpi3pi_2015_mdu_REV.root", "T");
    pair<string, string> kpi3pi_2015_mdu_mva_rev = pair("OutputFiles/kpi3pi_2015_mdu_MVA_REV.root", "T");
    TChain kpi3pi_2015_precuts_treechain("T");
    kpi3pi_2015_precuts_treechain.Add("InputFiles/kpi3pi_2015_md.root");
    kpi3pi_2015_precuts_treechain.Add("InputFiles/kpi3pi_2015_mu.root");

    RootFileCreatorFilterer(&kpi3pi_2015_precuts_treechain, kpi3pi_2015_mdu.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_2015_mdu, kpi3pi_2015_mdu_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_2015_mdu_rev, tmvamethods, kpi3pi_2015_mdu_mva_rev.first);
  }

  if (kpi3pi_dataset_2016)
  {    
    pair<string, string> kpi3pi_2016_mdu = pair("OutputFiles/kpi3pi_2016_mdu.root", "T");
    pair<string, string> kpi3pi_2016_mdu_rev = pair("OutputFiles/kpi3pi_2016_mdu_REV.root", "T");
    pair<string, string> kpi3pi_2016_mdu_mva_rev = pair("OutputFiles/kpi3pi_2016_mdu_MVA_REV.root", "T"); 
    TChain kpi3pi_2016_precuts_treechain("T");
    kpi3pi_2016_precuts_treechain.Add("InputFiles/kpi3pi_2016_md.root");
    kpi3pi_2016_precuts_treechain.Add("InputFiles/kpi3pi_2016_mu.root");
    
    RootFileCreatorFilterer(&kpi3pi_2016_precuts_treechain, kpi3pi_2016_mdu.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_2016_mdu, kpi3pi_2016_mdu_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_2016_mdu_rev, tmvamethods, kpi3pi_2016_mdu_mva_rev.first);
  }

  if (kpi3pi_dataset_2017)
  {    
    pair<string, string> kpi3pi_2017_mdu = pair("OutputFiles/kpi3pi_2017_mdu.root", "T");
    pair<string, string> kpi3pi_2017_mdu_rev = pair("OutputFiles/kpi3pi_2017_mdu_REV.root", "T");
    pair<string, string> kpi3pi_2017_mdu_mva_rev = pair("OutputFiles/kpi3pi_2017_mdu_MVA_REV.root", "T"); 
    TChain kpi3pi_2017_precuts_treechain("T");
    kpi3pi_2017_precuts_treechain.Add("InputFiles/kpi3pi_2017_md.root");
    kpi3pi_2017_precuts_treechain.Add("InputFiles/kpi3pi_2017_mu.root");
    
    RootFileCreatorFilterer(&kpi3pi_2017_precuts_treechain, kpi3pi_2017_mdu.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_2017_mdu, kpi3pi_2017_mdu_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_2017_mdu_rev, tmvamethods, kpi3pi_2017_mdu_mva_rev.first);
  }

  if (kpi3pi_dataset_2018)
  {    
    pair<string, string> kpi3pi_2018_mdu = pair("OutputFiles/kpi3pi_2018_mdu.root", "T");
    pair<string, string> kpi3pi_2018_mdu_rev = pair("OutputFiles/kpi3pi_2018_mdu_REV.root", "T");
    pair<string, string> kpi3pi_2018_mdu_mva_rev = pair("OutputFiles/kpi3pi_2018_mdu_MVA_REV.root", "T"); 
    TChain kpi3pi_2018_precuts_treechain("T");
    kpi3pi_2018_precuts_treechain.Add("InputFiles/kpi3pi_2018_md.root");
    kpi3pi_2018_precuts_treechain.Add("InputFiles/kpi3pi_2018_mu.root");
    
    RootFileCreatorFilterer(&kpi3pi_2018_precuts_treechain, kpi3pi_2018_mdu.first, cuts+stripping_cuts_xpiz);
    RootFileCreatorDefine(kpi3pi_2018_mdu, kpi3pi_2018_mdu_rev, var_tobe_defined);
    MVACutOptimizerReader(mvacut_pureomegafactory_outputfile, kpi3pi_2018_mdu_rev, tmvamethods, kpi3pi_2018_mdu_mva_rev.first);
  }
*/
  pair<string, string> kpi3pi_2015_mdu_mva_rev = pair("OutputFiles/kpi3pi_2015_mdu_MVA_REV.root", "T");
  pair<string, string> kpi3pi_2016_mdu_mva_rev = pair("OutputFiles/kpi3pi_2016_mdu_MVA_REV.root", "T");
  pair<string, string> kpi3pi_2017_mdu_mva_rev = pair("OutputFiles/kpi3pi_2017_mdu_MVA_REV.root", "T");
  pair<string, string> kpi3pi_2018_mdu_mva_rev = pair("OutputFiles/kpi3pi_2018_mdu_MVA_REV.root", "T");
  pair<string, string> kpi3pi_2015_2018_mvacuts_FullFactory_output = pair("OutputFiles/kpi3pi_2015_2018_bdtcut05_PureOmega50Variables.root", "T");
  pair<string, string> kpi3pi_2015_2018_mvacuts_PurgedFactory_output = pair("InputFiles/DATA_SignalRegion_PureOmega_MVA_SingleChannel_BDTCut.root", "DecayTree");
  pair<string, string> kKPi3Pi_SIM_JPsi_newvars = pair("InputFiles/KPi3Pi_SIM_JPsi_newvars.root", "T");
  pair<string, string> kpi3pi_sim_omega_mva_rev = pair("omegaMCStrippingFiltered.root", "DecayTree");
  pair<string, string> kpi3pi_sim_jpsi_mva_rev = pair("OutputFiles/kpi3pi_sim_jpsi_MVA_REV.root", "T");
  pair<string, string> kpi3pi_sim_jpsi_FullFactoryReader = pair("JpsiMCStrippingFiltered.root", "DecayTree");

  TChain kpi3pi_2015_2018_mvacuts_treechain("T");
  kpi3pi_2015_2018_mvacuts_treechain.Add("OutputFiles/kpi3pi_2015_mdu_MVA_REV.root");
  kpi3pi_2015_2018_mvacuts_treechain.Add("OutputFiles/kpi3pi_2016_mdu_MVA_REV.root");
  kpi3pi_2015_2018_mvacuts_treechain.Add("OutputFiles/kpi3pi_2017_mdu_MVA_REV.root");
  kpi3pi_2015_2018_mvacuts_treechain.Add("OutputFiles/kpi3pi_2018_mdu_MVA_REV.root");
  //RootFileCreatorFilterer(&kpi3pi_2015_2018_mvacuts_treechain, kpi3pi_2015_2018_mvacuts_output.first, "mvacut_kbdt1 > 0.5");

  AnalysisVariable B0_mass_jpsinotconstrained;
  B0_mass_jpsinotconstrained.SetVarBinningOpt(50, 5000, 5500);
  B0_mass_jpsinotconstrained.SetVarNames("b_m", "m_{NotConst}#left(#omegaK#pi#right)", "b_m", "#frac{MeV}{c^{2}}");  
  B0_mass_jpsinotconstrained.LegendPositionConstructor(1, 5, 0.22, 0.05);
  B0_mass_jpsinotconstrained.StatPositionConstructor(-1, 1, 0.24, 0.02);
  B0_mass_jpsinotconstrained.BoxHeaderPositionConstructor(3, 1, 0.22, 0.027);
  B0_mass_jpsinotconstrained.variable_padtextsize = 0.02;
  B0_mass_jpsinotconstrained.sampledata = "2015-2018";

  AnalysisVariable B0_mass_jpsiconstrained;
  B0_mass_jpsiconstrained.SetVarBinningOpt(50, 5000, 5500);
  B0_mass_jpsiconstrained.SetVarNames("b_cons_Jpsixpiz_m_best", "m_{Constrained}#left(J/#psiK#pi#right)", "b_m_jpsiconstrained", "#frac{MeV}{c^{2}}");  
  B0_mass_jpsiconstrained.LegendPositionConstructor(1, 5, 0.22, 0.05);
  B0_mass_jpsiconstrained.StatPositionConstructor(-1, 1, 0.24, 0.02);
  B0_mass_jpsiconstrained.BoxHeaderPositionConstructor(3, 1, 0.22, 0.027);
  B0_mass_jpsiconstrained.variable_padtextsize = 0.02;
  B0_mass_jpsiconstrained.sampledata = "2015-2018";
  //B0_mass_jpsiconstrained.VariablePlotter(&(kpi3pi_2015_2018_mvacuts_output), &B0_mass_jpsiconstrained, string("mvacut_kbdt1>0.5 && x_m_pi0constrained>3400 && x_m<3600"));
  //FitterMass(&kKPi3Pi_SIM_JPsi_newvars, &B0_mass_jpsiconstrained);

  AnalysisVariable jpsi_mass_pi0constrained;
  jpsi_mass_pi0constrained.SetVarBinningOpt(50, 2900, 3300);
  jpsi_mass_pi0constrained.SetVarNames("x_cons_xpiz_m_best", "m_{Constrained}#left(3#pi#right)", "x_m_pi0constrained", "#frac{MeV}{c^{2}}");  
  jpsi_mass_pi0constrained.LegendPositionConstructor(1, 5, 0.22, 0.05);
  jpsi_mass_pi0constrained.StatPositionConstructor(-1, 1, 0.22, 0.02);
  jpsi_mass_pi0constrained.BoxHeaderPositionConstructor(3, 1, 0.22, 0.027);
  jpsi_mass_pi0constrained.variable_padtextsize = 0.02;
  jpsi_mass_pi0constrained.sampledata = "2015-2018";
  //FitterMass(&kpi3pi_sim_jpsi_FullFactoryReader, &jpsi_mass_pi0constrained);
  
  //Fitter2DMass(&kpi3pi_2015_2018_mvacuts_PurgedFactory_output, &B0_mass_jpsiconstrained, &jpsi_mass_pi0constrained, "mvacut_kbdt1>0.5 && b_cons_xpiz_fitstatus_best==0 && b_cons_Jpsixpiz_fitstatus_best==0");

  AnalysisVariable B0_mass_omeganotconstrained;
  B0_mass_omeganotconstrained.SetVarBinningOpt(50, 5000, 5500);
  B0_mass_omeganotconstrained.SetVarNames("b_m", "m_{NotConst}#left(#omegaK#pi#right)", "b_m", "#frac{MeV}{c^{2}}");  
  B0_mass_omeganotconstrained.LegendPositionConstructor(1, 5, 0.22, 0.05);
  B0_mass_omeganotconstrained.StatPositionConstructor(-1, 1, 0.24, 0.02);
  B0_mass_omeganotconstrained.BoxHeaderPositionConstructor(3, 1, 0.22, 0.027);
  B0_mass_omeganotconstrained.variable_padtextsize = 0.02;
  B0_mass_omeganotconstrained.sampledata = "2015-2018";

  AnalysisVariable B0_mass_omegaconstrained;
  B0_mass_omegaconstrained.SetVarBinningOpt(50, 5000, 5500);
  B0_mass_omegaconstrained.SetVarNames("b_cons_omegaxpiz_m_best", "m_{Constrained}#left(#omegaK#pi#right)", "b_m_omegaconstrained", "#frac{MeV}{c^{2}}");  
  B0_mass_omegaconstrained.LegendPositionConstructor(1, 5, 0.22, 0.05);
  B0_mass_omegaconstrained.StatPositionConstructor(-1, 1, 0.24, 0.02);
  B0_mass_omegaconstrained.BoxHeaderPositionConstructor(3, 1, 0.22, 0.027);
  B0_mass_omegaconstrained.variable_padtextsize = 0.02;
  B0_mass_omegaconstrained.sampledata = "2015-2018";
  //B0_mass_omegaconstrained.VariablePlotter(&(kpi3pi_2015_2018_mvacuts_output), &B0_mass_omegaconstrained, string("mvacut_kbdt1>0.5 && x_m_pi0constrained>3400 && x_m<3600"));
  // FitterMass(&kpi3pi_sim_omega_mva_rev, &B0_mass_omegaconstrained);
  
  AnalysisVariable omega_mass;
  omega_mass.SetVarBinningOpt(50, 600, 1000);
  omega_mass.SetVarNames("x_m", "m_{Constrained}#left(3#pi#right)", "x_m_pi0constrained", "#frac{MeV}{c^{2}}");  
  omega_mass.LegendPositionConstructor(1, 5, 0.22, 0.05);
  omega_mass.StatPositionConstructor(-1, 1, 0.22, 0.02);
  omega_mass.BoxHeaderPositionConstructor(2, 1, 0.22, 0.027);
  omega_mass.variable_padtextsize = 0.02;
  omega_mass.sampledata = "2015-2018";
  //FitterMass(&kpi3pi_sim_omega_mva_rev, &omega_mass, "b_cons_xpiz_fitstatus_best==0 && b_cons_omegaxpiz_fitstatus_best==0");

  //B0_mass_omegaconstrained.VariablePlotter2D(&kpi3pi_sim_omega_mva_rev, &B0_mass_omeganotconstrained, &B0_mass_omegaconstrained);
  //Fitter2DMass(&kpi3pi_2015_2018_mvacuts_PurgedFactory_output, &B0_mass_omegaconstrained, &omega_mass, "mvacut_kbdt1>0.5  && b_cons_xpiz_fitstatus_best==0 && b_cons_omegaxpiz_fitstatus_best==0");
/*
  pair<string, string> kpi3pi_sim_jpsi = pair("InputFiles/kpi3pi-sim-jpsi.root", "T");
  pair<string, string> kpi3pi_sim_jpsi_rev = pair("OutputFiles/kpi3pi-sim-jpsi_REV.root", "T");
  X_mass_pi0constrained.SetVarBinningOpt(100, 3000, 3200);
  X_mass_pi0constrained.LegendPositionConstructor(4, 1, 0.22);
  X_mass_pi0constrained.StatPositionConstructor(-1, 1, 0.2);
  X_mass_pi0constrained.VariablePlotter(&(kpi3pi_sim_jpsi_rev), &X_mass_pi0constrained, "b_m_jpsiconstrained>5200&&b_m_jpsiconstrained<5500");  
  B0_mass_jpsiconstrained.SetVarBinningOpt(40, 4500, 6000);
  X_mass_pi0constrained.SetVarBinningOpt(40, 3000, 3200);
  X_mass_pi0constrained.VariablePlotter2D(&(kpi3pi_2015_2018_mvacuts_output), &B0_mass_jpsiconstrained, &X_mass_pi0constrained);
*/
/*
  ROOT::RDataFrame datafactoryframe = ROOT::RDataFrame("T", "OutputFiles/kpi3pi_sim_jpsi_MVA_REV.root");
  cout << datafactoryframe.Count().GetValue() << endl;
  for (vector<TMVAMethod>::iterator tmvamethit = tmvamethods.begin(); tmvamethit < tmvamethods.end(); tmvamethit++)
  {            
    cout << "Number of jpsi mc events after the " << tmvamethit->tmvamethodname << " cut: " << datafactoryframe.Filter(TString(("mvacut_")+(tmvamethit->tmvamethodcut)).Data()).Count().GetValue() << endl;
  }
*/

  //AnalysisVariable::VariablesComparer(&kpi3pi_files_distrostobecompared, vec_analvartobecompared);

  AnalysisVariable b_m_jpsiconstrained_k_pi0misidentifiedask;
  b_m_jpsiconstrained_k_pi0misidentifiedask.SetVarBinningOpt(50, 5000, 5800);
  b_m_jpsiconstrained_k_pi0misidentifiedask.SetVarNames("b_m_jpsiconstrained_k_pi0misidentifiedask", "m_{Constrained}#left(J/#psiK#pi^{K}_{MisId}#right)", "b_m_jpsiconstrained_k_pi0misidentifiedask", "#frac{MeV}{c^{2}}");  
  b_m_jpsiconstrained_k_pi0misidentifiedask.LegendPositionConstructor(1, 3, 0.22, 0.05);
  b_m_jpsiconstrained_k_pi0misidentifiedask.StatPositionConstructor(-1, 1, 0.22, 0.02);
  b_m_jpsiconstrained_k_pi0misidentifiedask.BoxHeaderPositionConstructor(2, 1, 0.22, 0.027);
  b_m_jpsiconstrained_k_pi0misidentifiedask.variable_padtextsize = 0.02;
  b_m_jpsiconstrained_k_pi0misidentifiedask.sampledata = "2015-2018";
  //B0_mass_jpsiconstrained.VariablePlotter(&(kKPi3Pi_SIM_JPsi_newvars), &b_m_jpsiconstrained_k_pi0misidentifiedask);
  //FitterMass(&kKPi3Pi_SIM_JPsi_newvars, &b_m_jpsiconstrained_k_pi0misidentifiedask);
  
  //TFile *input_signal_file = TFile::Open("/home/pasqualeandreola/DarkALP/InputFiles/kpi3pi_2015_md.root", "read");
  //TTree *signal = (TTree *)input_signal_file->Get("T");
  vector<string> infilestrings =
  {
    //"MC_B2KpiJpsi_2015_MagDown",
    //"MC_B2KpiJpsi_2015_MagUp",
    //"MC_B2KpiJpsi_2016_MagDown",
    //"MC_B2KpiJpsi_2016_MagUp",
    //"MC_B2KpiJpsi_2017_MagDown",
    //"MC_B2KpiJpsi_2017_MagUp",
    //"MC_B2KpiJpsi_2018_MagDown",
    //"MC_B2KpiJpsi_2018_MagUp",
    //"MC_B2Kpiomega_2015_MagDown",
    //"MC_B2Kpiomega_2015_MagUp",
    //"MC_B2Kpiomega_2016_MagDown",
    //"MC_B2Kpiomega_2016_MagUp",
    //"MC_B2Kpiomega_2017_MagDown",
    //"MC_B2Kpiomega_2017_MagUp",
    //"MC_B2Kpiomega_2018_MagDown",
    //"MC_B2Kpiomega_2018_MagUp"
    "DATA_2015_MagDown",
    //"DATA_2015_MagUp",
    //"DATA_2016_MagDown",
    //"DATA_2016_MagUp",
    //"DATA_2017_MagDown",
    //"DATA_2017_MagUp",
    //"DATA_2018_MagDown",
    //"DATA_2018_MagUp"
  };
  //for (auto &&infilestr : infilestrings)
  //  TreeSkimmer(infilestr, cuts);
  TChain kpi3pi_sim_omega_treechain("DecayTree");
  TChain kpi3pi_sim_jpsi_treechain("DecayTree");
  TChain kpi3pi_data_treechain("DecayTree");
  for (auto &&infilestr : infilestrings)
  {
    infilestr = "root://eoslhcb.cern.ch//eos/lhcb/user/p/paandreo/RUN2DATA_QEE_BTOKSTARX_19MAY2023/"+infilestr;
    if (infilestr.find("B2Kpiomega"))
    {
      kpi3pi_sim_omega_treechain.AddFile((infilestr+"_B2KpiX2PiPiPi.root").data());
      kpi3pi_sim_omega_treechain.AddFile((infilestr+"_B2KpiX2PiPiPiM.root").data());
    }
    if (infilestr.find("B2KpiJpsi")>=0)
    {
      kpi3pi_sim_jpsi_treechain.AddFile((infilestr+"_B2KpiX2PiPiPi.root").data());
      kpi3pi_sim_jpsi_treechain.AddFile((infilestr+"_B2KpiX2PiPiPiM.root").data());
    }
    if (infilestr.find("DATA"))
    {
      kpi3pi_data_treechain.AddFile((infilestr+"_B2KpiX2PiPiPi.root").data());
      kpi3pi_data_treechain.AddFile((infilestr+"_B2KpiX2PiPiPiM.root").data());
    }
  }
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagDown_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagDown_B2KpiX2PiPiPiM.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagUp_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagUp_B2KpiX2PiPiPiM.root");
  //TChain kpi3pi_sim_omega_treechainold("DecayTree");
  //kpi3pi_sim_omega_treechainold.Add("MC_B2Kpiomega_OLD_2016_MagDown_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechainold.Add("MC_B2Kpiomega_OLD_2016_MagDown_B2KpiX2PiPiPiM.root");
  //kpi3pi_sim_omega_treechainold.Add("MC_B2Kpiomega_OLD_2016_MagUp_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechainold.Add("MC_B2Kpiomega_OLD_2016_MagUp_B2KpiX2PiPiPiM.root");
  //TChain kpi3pi_sim_omega_treechain("DecayTree");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagDown_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagDown_B2KpiX2PiPiPiM.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagUp_B2KpiX2PiPiPi.root");
  //kpi3pi_sim_omega_treechain.Add("MC_B2Kpiomega_2016_MagUp_B2KpiX2PiPiPiM.root");
  //RootFileCreatorFilterer(&kpi3pi_sim_omega_treechain, "omegaMCtrippingFiltered.root", cuts+"&&x_truepid==223 && abs(x_MC_mother_ID)==511"
  //                                                                                             "&&abs(b_truepid)==511"
  //                                                                                             "&&abs(k_truepid)==321 && abs(k_MC_mother_ID)==313 && abs(k_MC_GD_mother_ID)==511"
  //                                                                                             "&&abs(pi_truepid)==211 && abs(k_MC_mother_ID)==313 && abs(pi_MC_GD_mother_ID)==511"
  //                                                                                             "&&abs(xpip_truepid)==211 && abs(xpip_MC_GD_mother_ID)==511 && abs(xpip_MC_mother_ID)==223"
  //                                                                                             "&&abs(xpim_truepid)==211 && abs(xpim_MC_GD_mother_ID)==511 && abs(xpim_MC_mother_ID)==223"
  //                                                                                             "&&abs(xpiz_truepid)==111 && abs(xpiz_MC_GD_mother_ID)==511 && abs(xpiz_MC_mother_ID)==223");
  //RootFileCreatorFilterer(&kpi3pi_sim_jpsi_treechain, "JpsiMCtrippingFiltered.root", cuts+"&&x_truepid==443 && abs(x_MC_mother_ID)==511"
  //                                                                                             "&&abs(b_truepid)==511"
  //                                                                                             "&&abs(k_truepid)==321 && abs(k_MC_mother_ID)==313 && abs(k_MC_GD_mother_ID)==511"
  //                                                                                             "&&abs(pi_truepid)==211 && abs(k_MC_mother_ID)==313 && abs(pi_MC_GD_mother_ID)==511"
  //                                                                                             "&&abs(xpip_truepid)==211 && abs(xpip_MC_GD_mother_ID)==511 && abs(xpip_MC_mother_ID)==443"
  //                                                                                             "&&abs(xpim_truepid)==211 && abs(xpim_MC_GD_mother_ID)==511 && abs(xpim_MC_mother_ID)==443"
  //                                                                                             "&&abs(xpiz_truepid)==111 && abs(xpiz_MC_GD_mother_ID)==511 && abs(xpiz_MC_mother_ID)==443");
  //TChain kpi3pi_sim_jpsi_old_treechain("T");
  //kpi3pi_sim_jpsi_old_treechain.Add("/home/pasqualeandreola/DarkALP/InputFiles/kpi3pi-sim-jpsi.root");
  //RootFileCreatorFilterer(&kpi3pi_data_treechain, "DATA_sidebands.root", "(b_m>4800 && b_m<5000) ||"
  //                                                                     "(b_m>5500 && b_m<5800)");
  
  //TChain kpi3pi_data_signalregion_mvacut_treechain("DecayTree");
  //kpi3pi_data_signalregion_mvacut_treechain.AddFile("root://eosuser.cern.ch//eos/user/p/paandreo/DATA_2018_MagUp_B2KpiX2PiPiPi.root");
  /*cuts = "(b_dira>0.999500)"
  "&&((b_vtx_chi2)/(b_vtx_ndof)<6.000000)"
  "&&(b_doca_kpi<0.200000)"
  "&&(b_fd>3.000000)"
  "&&(b_fd_chi2+x_fd_chi2>500.000000)"
  "&&((x_vtx_chi2)/(x_vtx_ndof)<4.500000)"
  "&&((x_z-b_z > -20) ||( (b_maxdoca > 0.2) && (x_z-b_z > 0) ) ||( (b_maxdoca < 1) && (x_z-b_z >20) ))"
  "&&((k_pt > 100 )||(pi_pt > 100)||(xpip_pt > 100)||(xpim_pt > 100))"
  "&&((k_p > 1500 )||(pi_p > 1500)||(xpip_p > 1500)||(xpim_p > 1500))"
  "&&(k_ip_chi2_best>10 || pi_ip_chi2_best>10 || xpip_ip_chi2_best>10 || xpim_ip_chi2_best>10)"
  "&&(k_ip_chi2_bvtx < 10 || pi_ip_chi2_bvtx < 10)"
  "&&(k_prb_ghost < 0.1 || pi_prb_ghost < 0.1 || xpip_prb_ghost < 0.1 || xpim_prb_ghost < 0.1)"
  "&&(k_pnn_k>0.300000)"
  "&&(b_m>4800.000000) && (b_m<5800.000000)"
  "&&(b_pt>3000.000000)"
  "&&(b_dira>0.999500)"
  "&&((b_vtx_chi2)/(b_vtx_ndof)<6.000000)"
  "&&(b_doca_kpi<0.200000)"
  "&&(b_fd>3.000000)"
  "&&(b_fd_chi2+x_fd_chi2>500.000000)"
  "&&((x_vtx_chi2)/(x_vtx_ndof)<4.500000)"
  "&&((x_z-b_z > -20) ||( (b_maxdoca > 0.2) && (x_z-b_z > 0) ) ||( (b_maxdoca < 1) && (x_z-b_z >20) ))"
  "&&((k_pt > 100 )||(pi_pt > 100)||(xpip_pt > 100)||(xpim_pt > 100))"
  "&&((k_p > 1500 )||(pi_p > 1500)||(xpip_p > 1500)||(xpim_p > 1500))"
  "&&(k_ip_chi2_best>10 || pi_ip_chi2_best>10 || xpip_ip_chi2_best>10 || xpim_ip_chi2_best>10)"
  "&&(k_ip_chi2_bvtx < 10 || pi_ip_chi2_bvtx < 10)"
  "&&(k_prb_ghost < 0.1 || pi_prb_ghost < 0.1 || xpip_prb_ghost < 0.1 || xpim_prb_ghost < 0.1)"
  "&&(k_pnn_k>0.300000)"
  "&&(k_ip_chi2_best+pi_ip_chi2_best+xpip_ip_chi2_best+xpim_ip_chi2_best > 200.)";
*/
  //RootFileCreatorFilterer(&kpi3pi_data_signalregion_mvacut_treechain, "DATA_signalregion_150kEVT.root", "", 1.5e5);
  
  // Checking if the file and the tree do exist
  int tuplebins = 100;
  ROOT::RDataFrame *newtuple=NULL, *oldtuple=NULL;
  oldtuple = new ROOT::RDataFrame("DecayTree", "DATA_SignalRegion_OmegaMVA_BDT_CUT.root");
  newtuple = new ROOT::RDataFrame("DecayTree", "DATA_SignalRegion_JpsiMVA_BDT_CUT_500kEVTS.root");
  auto oldtuple_def = std::make_unique<ROOT::RDF::RNode>(oldtuple->Range(0,5e5));
  //oldtuple_def->Snapshot("DecayTree", "DATA_SignalRegion_OmegaMVA_BDT_CUT_500kEVTS");
                                                                  //.Define("xvtxquality", "x_vtx_chi2/x_vtx_ndof")
                                                                  //.Define("b_vtxquality", "b_vtx_chi2/b_vtx_ndof"));
  auto newtuple_def = std::make_unique<ROOT::RDF::RNode>(newtuple->Range(0,5e5));
  //newtuple_def->Snapshot("DecayTree", "DATA_SignalRegion_JpsiMVA_BDT_CUT_500kEVTS");
                                                                  //.Define("xvtxquality", "x_vtx_chi2/x_vtx_ndof")
                                                                  //.Define("b_vtxquality", "b_vtx_chi2/b_vtx_ndof"));
  //auto newtuple_deffil = std::make_unique<ROOT::RDF::RNode>(newtuple->Filter("Stripping.compare(\"B2KpiX2PiPiPi\")==0")
  //                                                               .Define("xvtxquality", "x_vtx_chi2/x_vtx_ndof")
  //                                                                .Define("b_vtxquality", "b_vtx_chi2/b_vtx_ndof"));
  vec_analysisvariables->clear();
  AnalysisVariable xpipprbghost;
  xpipprbghost.SetVarNames("xpip_prb_ghost", "", "xpip_prb_ghost");
  xpipprbghost.SetVarBinningOpt(tuplebins,0,0.3);
  xpipprbghost.LegendPositionConstructor(1);
  //vec_analysisvariables->push_back(xpipprbghost);

  Double_t B0MASS = 5279.53, KMASS = 493.677, PIMASS = 139.57, PI0MASS = 134.9768;

  AnalysisVariable xpipp;
  xpipp.SetVarNames("xpip_p", "", "xpip_p");
  xpipp.SetVarBinningOpt(tuplebins,3000, 150000);
  //vec_analysisvariables->push_back(xpipp);
  AnalysisVariable xpippt;
  xpippt.SetVarNames("xpip_pt", "", "xpip_pt");
  xpippt.SetVarBinningOpt(tuplebins,250, 10000);
  //vec_analysisvariables->push_back(xpippt);
  AnalysisVariable xpipipchi2best;
  xpipipchi2best.SetVarNames("xpip_ip_chi2_best", "", "xpip_ip_chi2_best");
  xpipipchi2best.SetVarBinningOpt(tuplebins, 30, 10000);
  //vec_analysisvariables->push_back(xpipipchi2best);
  AnalysisVariable xpippnnpi;
  xpippnnpi.SetVarNames("xpip_pnn_pi", "", "xpip_pnn_pi");
  xpippnnpi.SetVarBinningOpt(tuplebins, 0.2, 1);
  xpippnnpi.LegendPositionConstructor(2);
  //vec_analysisvariables->push_back(xpippnnpi);
  AnalysisVariable xm;
  xm.SetVarNames("x_m", "", "x_m");
  xm.SetVarBinningOpt(tuplebins, 3000, 3300);
  //vec_analysisvariables->push_back(xm);
  AnalysisVariable xpt;
  xpt.SetVarNames("x_pt", "", "x_pt");
  xpt.SetVarBinningOpt(tuplebins, 2000, 25000);
  //vec_analysisvariables->push_back(xpt);
  AnalysisVariable xmaxdoca;
  xmaxdoca.SetVarNames("x_maxdoca", "", "x_maxdoca");
  xmaxdoca.SetVarBinningOpt(tuplebins, 0, 0.2);
  //vec_analysisvariables->push_back(xmaxdoca);
  AnalysisVariable xfdchi2;
  xfdchi2.SetVarNames("x_fd_chi2", "", "x_fd_chi2");
  xfdchi2.SetVarBinningOpt(tuplebins, 0, 15025);
  //vec_analysisvariables->push_back(xfdchi2);
  AnalysisVariable xvtxquality;
  xvtxquality.SetVarNames("xvtxquality", "", "xvtxquality");
  xvtxquality.SetVarBinningOpt(tuplebins, 0, 10);
  //vec_analysisvariables->push_back(xvtxquality);
  AnalysisVariable kp;
  kp.SetVarNames("k_p", "", "k_p");
  kp.SetVarBinningOpt(tuplebins,2000, 150000);
  //vec_analysisvariables->push_back(kp);
  AnalysisVariable kpt;
  kpt.SetVarNames("k_pt", "", "k_pt");
  kpt.SetVarBinningOpt(tuplebins,250, 25000);
  //vec_analysisvariables->push_back(kpt);
  AnalysisVariable kipchi2best;
  kipchi2best.SetVarNames("k_ip_chi2_best", "", "k_ip_chi2_best");
  kipchi2best.SetVarBinningOpt(tuplebins, 9, 10000);
  //vec_analysisvariables->push_back(kipchi2best);
  AnalysisVariable kpnnk;
  kpnnk.SetVarNames("k_pnn_k", "", "k_pnn_k");
  kpnnk.SetVarBinningOpt(tuplebins, 0.1, 1);
  kpnnk.LegendPositionConstructor(2);
  //vec_analysisvariables->push_back(kpnnk);
  AnalysisVariable pip;
  pip.SetVarNames("pi_p", "", "pi_p");
  pip.SetVarBinningOpt(tuplebins,2000, 150000);
  //vec_analysisvariables->push_back(pip);
  AnalysisVariable pipt;
  pipt.SetVarNames("pi_pt", "", "pi_pt");
  pipt.SetVarBinningOpt(tuplebins,250, 25000);
  //vec_analysisvariables->push_back(pipt);
  AnalysisVariable piipchi2best;
  piipchi2best.SetVarNames("pi_ip_chi2_best", "", "pi_ip_chi2_best");
  piipchi2best.SetVarBinningOpt(tuplebins, 9, 10000);
  //vec_analysisvariables->push_back(piipchi2best);
  AnalysisVariable pipnnpi;
  pipnnpi.SetVarNames("pi_pnn_pi", "", "pi_pnn_pi");
  pipnnpi.SetVarBinningOpt(tuplebins, 0.2, 1);
  pipnnpi.LegendPositionConstructor(2);
  //vec_analysisvariables->push_back(pipnnpi);
  AnalysisVariable b_m;
  b_m.SetVarNames("b_m", "m#left(B^{0}#right)", "b_m", "MeV");
  b_m.SetVarBinningOpt(tuplebins, 5000, 5500);
  //vec_analysisvariables->push_back(b_m);
  AnalysisVariable b_pt;
  b_pt.SetVarNames("b_pt", "", "b_pt");
  b_pt.SetVarBinningOpt(tuplebins, 3000, 30000);
  //vec_analysisvariables->push_back(b_pt);
  AnalysisVariable b_tau;
  b_tau.SetVarNames("b_tau", "", "b_tau");
  b_tau.SetVarBinningOpt(tuplebins, 0.0002, 0.02);
  //vec_analysisvariables->push_back(b_tau);
  AnalysisVariable b_ip_chi2;
  b_ip_chi2.SetVarNames("b_ip_chi2", "", "b_ip_chi2");
  b_ip_chi2.SetVarBinningOpt(tuplebins, 0, 10);
  //vec_analysisvariables->push_back(b_ip_chi2);
  AnalysisVariable b_dira;
  b_dira.SetVarNames("b_dira", "", "b_dira");
  b_dira.SetVarBinningOpt(tuplebins, 0.9997, 1);
  b_dira.LegendPositionConstructor(2);
  //vec_analysisvariables->push_back(b_dira);
  AnalysisVariable b_vtxquality;
  b_vtxquality.SetVarNames("b_vtxquality", "", "b_vtxquality");
  b_vtxquality.SetVarBinningOpt(tuplebins, 0, 15);
  //vec_analysisvariables->push_back(b_vtxquality);  
  AnalysisVariable b_cons_omegaxpiz_m_best;
  b_cons_omegaxpiz_m_best.SetVarNames("b_cons_xpiz_m_best", "m(K#pi3#pi)_{CONS#left(#omega,#pi^{0}#right)}", "b_cons_omegaxpiz_m_best", "MeV");
  b_cons_omegaxpiz_m_best.SetVarBinningOpt(20, 5e3, 5.5e3);
  vec_analysisvariables->push_back(b_cons_omegaxpiz_m_best);
  AnalysisVariable b_cons_jpsixpiz_m_best;
  b_cons_jpsixpiz_m_best.SetVarNames("b_cons_jpsixpiz_m_best", "m(K#pi3#pi)_{CONS#left(J/#Psi,#pi^{0}#right)}", "b_cons_Jpsixpiz_m_best", "MeV");
  b_cons_jpsixpiz_m_best.SetVarBinningOpt(20, 5e3, 5.5e3);
  vec_analysisvariables->push_back(b_cons_jpsixpiz_m_best);
  AnalysisVariable x_cons_xpiz_m_best_jpsi;
  x_cons_xpiz_m_best_jpsi.SetVarNames("x_cons_xpiz_m_best", "m(3#pi)_{CONS#left(#pi^{0}#right)}", "x_cons_xpiz_m_best_jpsi", "MeV");
  x_cons_xpiz_m_best_jpsi.SetVarBinningOpt(20, 2500, 3500);
  vec_analysisvariables->push_back(x_cons_xpiz_m_best_jpsi);
  AnalysisVariable x_cons_xpiz_m_best_omega;
  x_cons_xpiz_m_best_omega.SetVarNames("x_cons_xpiz_m_best", "m(3#pi)_{CONS#left(#pi^{0}#right)}", "x_cons_xpiz_m_best", "MeV");
  x_cons_xpiz_m_best_omega.SetVarBinningOpt(50, 600, 1100);
  vec_analysisvariables->push_back(x_cons_xpiz_m_best_omega);  
  AnalysisVariable x_cons_xpiz_m_best;
  x_cons_xpiz_m_best.SetVarNames("x_cons_xpiz_m_best", "m(3#pi)_{CONS#left(#pi^{0}#right)}", "x_cons_xpiz_m_best", "MeV");
  x_cons_xpiz_m_best.SetVarBinningOpt(tuplebins, 500, 1500);
  vec_analysisvariables->push_back(x_cons_xpiz_m_best);  
  AnalysisVariable s_kpi_fullrange;
  s_kpi_fullrange.SetVarNames("s_kpigev_fullrange", "s#left(K^{+},#pi^{-}#right)", "sqrt(s_kpi*1e-6)", "GeV");
  s_kpi_fullrange.SetVarBinningOpt(50, sqrt(1e-6*(KMASS+PIMASS)*(KMASS+PIMASS)), sqrt(3+1e-6*(B0MASS-(2*PIMASS+PI0MASS))*(B0MASS-(2*PIMASS+PI0MASS))));
  //vec_analysisvariables->push_back(s_kpi);  
  AnalysisVariable s_kpi;
  s_kpi.SetVarNames("s_kpigev_minirange", "s#left(K^{+},#pi^{-}#right)", "sqrt(s_kpi*1e-6)", "GeV");
  s_kpi.SetVarBinningOpt(50, sqrt(2.5+1e-6*(KMASS+PIMASS)*(KMASS+PIMASS)), sqrt(-20+1e-6*(B0MASS-(2*PIMASS+PI0MASS))*(B0MASS-(2*PIMASS+PI0MASS))));
  s_kpi.SetVarBinningOpt(50, 1.4, 1.8);
  vec_analysisvariables->push_back(s_kpi);  
  AnalysisVariable s_xpi_fullrange;
  s_xpi_fullrange.SetVarNames("s_xpigev_fullrange", "s#left(3#pi,#pi^{-}#right)", "s_xpi*1e-6", "GeV");
  s_xpi_fullrange.SetVarBinningOpt(50, 1e-6*((2*PIMASS+PI0MASS)+PIMASS)*((2*PIMASS+PI0MASS)+PIMASS), -20+1e-6*(B0MASS-KMASS)*(B0MASS-KMASS));
  vec_analysisvariables->push_back(s_xpi_fullrange);
  AnalysisVariable s_xpi;
  s_xpi.SetVarNames("s_xpigev", "s#left(3#pi,#pi^{-}#right)", "s_xpi*1e-6", "GeV");
  s_xpi.SetVarBinningOpt(100, -2+1e-6*((2*PIMASS+PI0MASS)+PIMASS)*((2*PIMASS+PI0MASS)+PIMASS), 5+1e-6*(B0MASS-KMASS)*(B0MASS-KMASS));
  vec_analysisvariables->push_back(s_xpi);    
  AnalysisVariable s_xk;
  s_xk.SetVarNames("s_xkgev", "s#left(3#pi,K^{+}#right)", "s_xk*1e-6", "GeV");
  s_xk.SetVarBinningOpt(100, -2+1e-6*((2*PIMASS+PI0MASS)+PIMASS)*((2*PIMASS+PI0MASS)+PIMASS), 3+1e-6*(B0MASS-PIMASS)*(B0MASS-PIMASS));
  vec_analysisvariables->push_back(s_xk);  
  AnalysisVariable s_kxpi;
  s_kxpi.SetVarNames("s_kxpigev", "s#left(K^{+}#pi^{-}_{x}#right)", "sqrt(s_kxpi*1e-6)", "GeV");
  s_kxpi.SetVarBinningOpt(50, 0, 4);
  vec_analysisvariables->push_back(s_kxpi);  
  AnalysisVariable s_pixpi;
  s_pixpi.SetVarNames("s_pixpigev", "s#left(#pi^{+}#pi^{-}_{x}#right)", "s_pixpi*1e-6", "GeV");
  s_pixpi.SetVarBinningOpt(50, 0, 16);
  vec_analysisvariables->push_back(s_pixpi);  
  AnalysisVariable mvacut_kbdt1;
  mvacut_kbdt1.SetVarNames("mvacut_kbdt1", "MVA Cut", "mvacut_kbdt1", "");
  mvacut_kbdt1.SetVarBinningOpt(100, 0.15, 0.7);
  vec_analysisvariables->push_back(mvacut_kbdt1);  
  
  /*
  for (vector<AnalysisVariable>::iterator analvar=vec_analysisvariables->begin(); analvar!=vec_analysisvariables->end(); analvar++)
  {
    string oldvarname = string(analvar->variable_name)+"old";
    string newvarname =(string(analvar->variable_name));
    Double_t minbin = analvar->variable_histmin, maxbin = analvar->variable_histmax;
    int oldentries = oldtuple_def->Count().GetValue(), newentries = newtuple_def->Count().GetValue();
    auto hist_old_var = oldtuple_def->Fill(TH1F(oldvarname.data(), oldvarname.data(), tuplebins, minbin, maxbin),{analvar->variable_expression});
    auto hist_new_var = newtuple_def->Fill(TH1F(newvarname.data(), newvarname.data(), tuplebins, minbin, maxbin),{analvar->variable_expression});
    hist_old_var->Sumw2(); hist_new_var->Sumw2();
    hist_old_var->Scale(1/hist_old_var->Integral()); hist_new_var->Scale(1/hist_new_var->Integral());
    hist_old_var->ClearUnderflowAndOverflow();
    hist_new_var->ClearUnderflowAndOverflow();
    hist_old_var->SetLineColor(kRed);
    hist_new_var->SetLineColor(kBlue);
    auto c1 = new TCanvas("c1", "A ratio example", 960, 1080);
    TRatioPlot* histratiocomparer = new TRatioPlot(hist_new_var.GetPtr(), hist_old_var.GetPtr(), "diff pois");
    gStyle->SetOptStat(0);
    c1->SetTicks(0, 1);
    histratiocomparer->Draw();
    gPad->Modified(); gPad->Update(); // make sure it’s really (re)drawn
    TPad *p = histratiocomparer->GetUpperPad();
    TLegend *l = analvar->SetLegendPosAuto(analvar->variable_padposition, 2); 
    l->SetTextSize(0.01);
    l->SetName(newvarname.data());
    l->AddEntry(hist_old_var->GetName(), string("PureOmega_MVA ("+to_string(oldentries)+"entries)").data());
    l->AddEntry(hist_new_var->GetName(), string("OmegaJpsi_MVA ("+to_string(newentries)+"entries)").data());
    l->Draw("SAME");
    p->Modified(); p->Update();
    c1->SaveAs(string("OutputFiles/PNGPlots/VariablesComparisonForDifferentProductions/"+string(newvarname)+".png").data());
    c1->SetCanvasSize(960,1080);
    c1->SaveAs(string("OutputFiles/PNGPlots/VariablesComparisonForDifferentProductionsHW/TwoStripping/"+string(newvarname)+".png").data());
  }
*/
/*
  ROOT::RDataFrame *deftuple = new ROOT::RDataFrame("DecayTree", "InputFiles/DATA_SignalRegion_JpsiMVA_BDT_20072023.root");

  vector<string> s_kpi_vars =    {"k_px", "k_py", "k_pz", "k_e",
                                  "pi_px", "pi_py", "pi_pz", "pi_e"
                                 };
  vector<string> s_xpi_vars =    {"x_px", "x_py", "x_pz", "x_e",
                                  "pi_px", "pi_py", "pi_pz", "pi_e"
                                 };
  vector<string> s_xk_vars =     {"x_px", "x_py", "x_pz", "x_e",
                                  "k_px", "k_py", "k_pz", "k_e"
                                 };
  vector<string> s_kxpi_vars =   {"k_px", "k_py", "k_pz", "k_e",
                                  "xpim_px", "xpim_py", "xpim_pz", "xpim_e"
                                 };
  vector<string> s_pixpi_vars =  {"pi_px", "pi_py", "pi_pz", "pi_e",
                                  "xpip_px", "xpip_py", "xpip_pz", "xpip_e"
                                 };
  
  auto deftuple_node = std::make_unique<ROOT::RDF::RNode>(deftuple->//Filter("mvacut_kbdt1>0.47")
                                                                   Define("s_kpi", invariant_mass_couple, s_kpi_vars)
                                                                   .Define("s_xpi", invariant_mass_couple, s_xpi_vars)
                                                                   .Define("s_xk", invariant_mass_couple, s_xk_vars)
                                                                   .Define("s_kxpi", invariant_mass_couple, s_kxpi_vars)
                                                                   .Define("s_pixpi", invariant_mass_couple, s_pixpi_vars)
                                                                   );
  ROOT::RDF::RSnapshotOptions snapopt;
  snapopt.fMode = "UPDATE";
  snapopt.fOverwriteIfExists = "TRUE";
  deftuple_node->Snapshot("DecayTree", "InputFiles/DATA_SignalRegion_JpsiMVA_BDT_20072023_OmegaCUT.root");
*/
  unordered_map<string, vector<string>> vardef;
  /*vardef["b_m_omegaconstrained"] = {
    "omega_PX", "omega_PY", "omega_PZ",
    "Kp_PX", "Kp_PY", "Kp_PZ", "Kp_E",
    "pim_PX", "pim_PY", "pim_PZ", "pim_E"
    };*/
  /*
  vardef["s_kpi"] = s_kpi_vars;
  vardef["s_xpi"] = s_xpi_vars;
  vardef["s_xk"]  = s_xk_vars;
  vardef["s_kxpi"] = s_kxpi_vars;
  vardef["s_pixpi"] = s_pixpi_vars;*/
  //pair<string, string> Bp2Kpipipomega_tree = pair("root://eosuser.cern.ch//eos/user/p/paandreo/B2Domega_tree.root", "DecayTree");
  //pair<string, string> B2K1omega_tree_filterrev = pair("B2Domega_tree.root", "DecayTree");
  pair<string, string> Bp2Kpipipomega_tree = pair("root://eosuser.cern.ch//eos/user/p/paandreo//DATA_SignalRegion_OmegaMVA_BDT_20072023_OmegaCUT.root", "DecayTree");
  pair<string, string> B2K1omega_tree_filterrev = pair("InputFiles/DATA_SignalRegion_OmegaMVA_BDT_20072023_OmegaCUT.root", "DecayTree");
  //RootFileCreatorDefine(Bp2Kpipipomega_tree, B2K1omega_tree_filterrev, vardef);
  unordered_map<string, pair<string, string>> part_reco_files;
  part_reco_files["K1omega"] = pair<string, string>("B2K1omega_tree.root", "DecayTree");
  part_reco_files["Kpigammaomega"] = pair<string, string>("B2Kpigammaomega_tree.root", "DecayTree");
  part_reco_files["Kpipiomega"] = pair<string, string>("B2Kpipiomega_tree.root", "DecayTree");
  part_reco_files["Dstar2Dgammaomega"] = pair<string, string>("B2Dstaromega2Dgamma3pi2Kpi3pi_tree.root", "DecayTree");
  part_reco_files["Dstar2Dpi0omega"] = pair<string, string>("B2Dstaromega2Dpiz3pi2Kpi3pi_tree.root", "DecayTree");
  part_reco_files["Domega"] = pair<string, string>("B2Domega_tree.root", "DecayTree");
  // part_reco_files["Kpipipomega"] = pair<string, string>("Bp2Kpipipomega_tree.root", "DecayTree");
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &s_kpi, &s_xpi);
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &s_kpi, &s_xk);
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &s_xk, &s_xpi);
  pair<string, string> *kpi3pi_darkbosonstrippingtuple = new pair<string, string>("/home/pasqualeandreola/Scaricati/B2KpiALP(500MeV,1ps)2pipipi_MagDown.root", "B2KpiX2PiPiPi_Tuple/DecayTree");
  // ROCCurve(*kpi3pi_darkbosonstrippingtuple, &B0_PT, "abs(B0_TRUEID)==511&&abs(KS0_TRUEID)==25&&abs(Kplus_TRUEID)==211&&abs(piminus_TRUEID)==321&&abs(xpiz_TRUEID)==111&&abs(xpim_TRUEID)==211&&abs(xpip_TRUEID)==211");
  pair<string, string> *kpi3pi_data_signalregion_omegamvacut = new pair<string, string>("InputFiles/DATA_SignalRegion_OmegaMVA_BDT_20072023_OmegaCUT.root", "DecayTree");
   pair<string, string> *kpi3pi_data_signalregion_jpsimvacut = new pair<string, string>("InputFiles/DATA_SignalRegion_JpsiMVA_BDT_20072023.root", "DecayTree");
  string varplot2d_omega_mincut = "(x_cons_xpiz_m_best<1100) & (x_cons_xpiz_m_best>600) & ((x_cons_xpiz_m_best<750) | (x_cons_xpiz_m_best>820)) & (b_cons_omegaxpiz_m_best<5500) & (b_cons_omegaxpiz_m_best>5000) & (Pi0Merged == 0) &";
  string varplot2d_jpsi_mincut = "(x_m<3500) & (x_m>2500) & (b_m<5500) & (b_m>5000) & (Pi0Merged == 0)";
  string varplot2d_kxpi_idcut = "( ((k_pid>0)&(xpim_pid<0)) | ((k_pid<0)&(xpim_pid>0)) ) &";
  string varplot2d_pixpi_idcut = "( ((pi_pid>0)&(xpip_pid<0)) | ((pi_pid<0)&(xpip_pid>0)) ) &";
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &x_cons_xpiz_m_best_omega, "(s_kpi<1.844*1.844*1e6 | s_kpi>1.884*1.884*1e6) && ((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6)) && (x_cons_xpiz_m_best>820 || x_cons_xpiz_m_best<740)");
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_kpi_fullrange, "s_kpi<1.844*1.844*1e6 | s_kpi>1.884*1.884*1e6");
  B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &x_cons_xpiz_m_best_omega, varplot2d_omega_mincut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_jpsimvacut, &b_cons_omegaxpiz_m_best, &x_cons_xpiz_m_best_jpsi, varplot2d_jpsi_mincut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_jpsimvacut, &mvacut_kbdt1, &x_cons_xpiz_m_best_jpsi, varplot2d_jpsi_mincut, true);
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_kxpi, varplot2d_omega_mincut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_pixpi, varplot2d_omega_mincut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &x_cons_xpiz_m_best_omega, &s_kxpi, varplot2d_omega_mincut+varplot2d_kxpi_idcut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  // B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &x_cons_xpiz_m_best_omega, &s_pixpi, varplot2d_omega_mincut+varplot2d_pixpi_idcut+"((s_kpi>0.6*0.6*1e6) & (s_kpi<1*1*1e6))");
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_xpi);
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_xpi_fullrange);
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_cons_omegaxpiz_m_best, &s_xk);
  //  //B0_mass_omegaconstrained.VariablePlotter2D(kpi3pi_data_signalregion_omegamvacut, &b_m, &s_kpi, "x_m>600 && x_m<1000");
  //  s_kpi_fullrange.StatPositionConstructor(-1, 1, 0.24, 0.02);
  b_cons_omegaxpiz_m_best.variable_plot_normalized_flag=true;
  string rapidcuts = "(B0_M>4.8) && (B0_M<5.8)"
                     "&& (B0_PT>3.000000)"
                     "&& (B0_FD>3.000000)"
                     "&& ((Kp_PT > 0.1 )||(pim_PT > 0.1)||(xpip_PT > 0.1)||(xpim_PT > 0.1))&&((Kp_P > 1.5 )||(pim_P > 1.5)||(xpip_P > 1.5)||(xpim_P > 1.5))";
  AnalysisVariable B0_M;
  B0_M.SetVarNames("B0_M", "m#left(B^{0}#right)", "B0_M", "MeV");
  B0_M.SetVarBinningOpt(tuplebins, 5, 5.5);
  B0_M.variable_plot_normalized_flag=true;
  AnalysisVariable xmfull;
  xmfull.SetVarNames("x_m", "m#left(3#pi#right)", "x_m", "MeV");
  xmfull.SetVarBinningOpt(tuplebins, 500, 3300);
  //xmfull.variable_statpadvertices(0,7,0.1);
  xmfull.variable_plot_normalized_flag=false;
  AnalysisVariable::VariablePlotter(kpi3pi_data_signalregion_jpsimvacut, &xmfull);
  AnalysisVariable::VariableComparer(&part_reco_files, &B0_M);
  //pair<string, string> *kpi3pi_data_signalregion_jpsimvacut = new pair<string, string>("DATA_SignalRegion_JpsiMVA_BDT_CUT.root", "DecayTree");
/*
  for (vector<AnalysisVariable>::iterator analvar=vec_analysisvariables->begin(); analvar!=vec_analysisvariables->end(); analvar++)
  {
    string oldvarname = string(analvar->variable_name)+"old";
    string newvarname =(string(analvar->variable_name));
    Double_t minbin = analvar->variable_histmin, maxbin = analvar->variable_histmax;
    int oldentries = oldtuple_def->Count().GetValue(), newentries = newtuple_deffil->Count().GetValue();
    auto hist_old_var = oldtuple_def->Fill(TH1F(oldvarname.data(), oldvarname.data(), tuplebins, minbin, maxbin),{analvar->variable_expression});
    auto hist_new_var = newtuple_deffil->Fill(TH1F(newvarname.data(), newvarname.data(), tuplebins, minbin, maxbin),{analvar->variable_expression});
    hist_old_var->Sumw2(); hist_new_var->Sumw2();
    hist_old_var->Scale(1/hist_old_var->Integral()); hist_new_var->Scale(1/hist_new_var->Integral());
    hist_old_var->ClearUnderflowAndOverflow();
    hist_new_var->ClearUnderflowAndOverflow();
    hist_old_var->SetLineColor(kRed);
    hist_new_var->SetLineColor(kBlue);
    auto c1 = new TCanvas("c1", "A ratio example", 1920, 1080);
    TRatioPlot* histratiocomparer = new TRatioPlot(hist_new_var.GetPtr(), hist_old_var.GetPtr(), "diff pois");
    gStyle->SetOptStat(0);
    c1->SetTicks(0, 1);
    histratiocomparer->Draw();
    gPad->Modified(); gPad->Update(); // make sure it’s really (re)drawn
    TPad *p = histratiocomparer->GetUpperPad();
    TLegend *l = analvar->SetLegendPosAuto(analvar->variable_padposition, 2); 
    l->SetName(newvarname.data());
    l->AddEntry(hist_old_var->GetName(), string("Old Tuple ("+to_string(oldentries)+"entries)").data());
    l->AddEntry(hist_new_var->GetName(), string("New Tuple ("+to_string(newentries)+"entries)").data());
    l->Draw("SAME");
    p->Modified(); p->Update();
    c1->SaveAs(string("OutputFiles/PNGPlots/VariablesComparisonForDifferentProductions/"+string(newvarname)+".png").data());
    c1->SetCanvasSize(960,1080);
    c1->SaveAs(string("OutputFiles/PNGPlots/VariablesComparisonForDifferentProductionsHW/OneStripping/"+string(newvarname)+".png").data());
    delete c1;
  }*/

  ElapsedTimeStamper(start);
  end = chrono::system_clock::now();
  end_time = chrono::system_clock::to_time_t(end);
  cout << "Ending date: " << ctime(&end_time) << endl;
  return 0;
}