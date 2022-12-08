#include "HeaderFiles/Fitter.h"

using namespace std;
using namespace RooFit;

int Fitter2DMass(pair<string, string> *input_file_tree,
               AnalysisVariable *analvar1_tobefit,
               AnalysisVariable *analvar2_tobefit,
               string selection_cut,
               string mc_selection_cut,
               bool debug,
               string closuretest)
{
    // Silencing output
    RooFit::PrintLevel(-1);
    RooMinimizer::PrintLevel(RooMinimizer::None);

    // Variable name
    const char *varfit1_name = analvar1_tobefit->variable_name;
    const char *varfit1_title = analvar1_tobefit->variable_prettyname;
    const char *varfit1_dimension = analvar1_tobefit->variable_dimension;
    const char *varfit1_sample_data = analvar1_tobefit->sampledata.data();
    const char *varfit2_name = analvar2_tobefit->variable_name;
    const char *varfit2_title = analvar2_tobefit->variable_prettyname;
    const char *varfit2_dimension = analvar2_tobefit->variable_dimension;
    const char *varfit2_sample_data = analvar2_tobefit->sampledata.data();

    // Fit variable parameters
    Double_t bins1 = analvar1_tobefit->variable_bins, minbin1 = analvar1_tobefit->variable_histmin, maxbin1 = analvar1_tobefit->variable_histmax, bin_width1 = (maxbin1 - minbin1) / bins1;
    Double_t bins2 = analvar2_tobefit->variable_bins, minbin2 = analvar2_tobefit->variable_histmin, maxbin2 = analvar2_tobefit->variable_histmax, bin_width2 = (maxbin2 - minbin2) / bins2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      Variables and PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector<RooRealVar *> roofitinputvar;
    RooRealVar varfit1 = RooRealVar(varfit1_name, varfit1_title, minbin1, maxbin1, varfit1_dimension);
    roofitinputvar.push_back(&varfit1);
    RooRealVar varfit2 = RooRealVar(varfit2_name, varfit2_title, minbin2, maxbin2, varfit2_dimension);
    roofitinputvar.push_back(&varfit2);
    RooRealVar mvacut_kbdt1 = RooRealVar("mvacut_kbdt1", "BDT1 MVA cut", 0, 1, "");
    roofitinputvar.push_back(&mvacut_kbdt1);
    RooRealVar x_m = RooRealVar("x_m", "Mass of the children cut", 0, 1e6, "");
    roofitinputvar.push_back(&x_m);
    RooRealVar b0_mass = RooRealVar("b0_mass", "b0_mass", 5279.34);
    RooFormulaVar jpsi_mbconstrained_jpsicorrected = RooFormulaVar("x_m_bconstraint", "@1-(@2-@3)", RooArgList(varfit2, b0_mass, varfit1));
    RooRealVar jpsimass = RooRealVar("jpsi_mass", "jpsi_mass", 3096.9);
    RooFormulaVar b_mjpsiconstrained_pi0corrected = RooFormulaVar("x_m_jpsiconstraint", "b_m_jpsiconstrained-(jpsi_mass-x_m_pi0constrained)", RooArgList(varfit1, jpsimass, varfit2));
    RooRealVar omegamass = RooRealVar("omega_mass", "omega_mass", 782.66);
    //RooFormulaVar b_momegaconstrained_pi0corrected = RooFormulaVar("x_m_omegaconstraint", "b_m_omegaconstrained-(omega_mass-x_m_pi0constrained)", RooArgList(varfit1, omegamass, varfit2));

    RooArgSet thevars = RooArgSet();
    for (vector<RooRealVar *>::iterator rooiter = roofitinputvar.begin(); rooiter != roofitinputvar.end(); rooiter++)
        thevars.add(**rooiter);

    const char *data_file_name = input_file_tree->first.data();
    const char *data_tree_name = input_file_tree->second.data();
    TTree *data_tree = (TTree *)(TFile::Open(data_file_name)->Get(data_tree_name));
    RooDataSet fulldata = RooDataSet("data", "data", data_tree->GetTree(), thevars, selection_cut.data());
    Int_t dataentries = fulldata.numEntries();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      varfit1 ranges
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*Double_t fit_range_lo = minbin1;
    Double_t mass_window_lo = (minbin1+maxbin1)/2 - 150;
    Double_t mass_window_hi = (minbin1+maxbin1)/2 + 150;
    Double_t fit_range_hi = maxbin1;*/
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // expo
    RooRealVar slope1 = RooRealVar("slope1", "slope1", -7.53684e-04, -1e5, 1e5);
    RooExponential bkg1 = RooExponential((string(varfit1_name)+"component_bkg_expo1").data(), "B^{0} Comb. BKG", varfit1, slope1);

    RooRealVar slope2 = RooRealVar("slope2", "slope2", -2.82889e-03, -1e5, 1e5);
    RooExponential bkg2 = RooExponential((string(varfit2_name)+"component_bkg_expo2").data(), "#omega Comb. BKG", varfit2, slope2);

    RooRealVar slope3 = RooRealVar("slope3", "slope3", -0.1, -1e6, 1e6);
    RooExponential bkg3 = RooExponential("bkg_expo3", "bkg_expo3", b_mjpsiconstrained_pi0corrected, slope3);

    RooRealVar b_omega_slope = RooRealVar("b_omega_slope", "B^{0}-#omega Slope", -0.1, -1e6, 1e6);
    RooExponential b_omega_expo = RooExponential("b_omega_bkg_expo", "B^{0}-#omega Expo", varfit1, b_omega_slope);

    // Gaussian
    RooRealVar gaussian_mean = RooRealVar("gaussian_mean", "gaussian_mean", 5279, minbin1, maxbin1);
    RooRealVar gaussian_width = RooRealVar("gaussian_width", "gaussian_width", 1.52423e+01, 1e-6, 200);
    RooGaussian signal_gauss = RooGaussian((string(varfit1_name)+"component_Signal_Gaus").data(), "B^{0} Gaussian", varfit1, gaussian_mean, gaussian_width);

    // Breit Wigner
    RooRealVar gaussian_mean2 = RooRealVar("gaussian_mean2", "gaussian_mean2", 3.09950e+03, minbin2, maxbin2);
    RooRealVar gaussian_width2 = RooRealVar("gaussian_width2", "gaussian_width2", 3.90478e+01, 1e-6, 200);
    RooGaussian signal_gauss2 = RooGaussian((string(varfit2_name)+"component_Signal_Gaus2").data(), "J/#psi Gaussian", varfit2, gaussian_mean2, gaussian_width2);

    // Breit Wigner
    RooRealVar breit_mean_omega = RooRealVar("breit_mean_omega", "#mu#left(m_{#omega}#right)", 782, minbin2, maxbin2);
    RooRealVar breit_width_omega = RooRealVar("breit_width_omega", "#sigma#left(m_{#omega}#right)", 25, 1e-6, 200);
    RooBreitWigner breit_wigner_omega = RooBreitWigner((string(varfit2_name)+"component_breit_wigner_omega").data(), "#omega Breit Wigner", varfit2, breit_mean_omega, breit_width_omega);

    // Voigtian
    RooRealVar voigtian_mean_omega = RooRealVar("voigtian_mean_omega", "#mu#left(m_{#omega}#right)", 782, minbin2, maxbin2);
    RooRealVar voigtian_width_omega = RooRealVar("voigtian_width_omega", "#sigma#left(m_{#omega}#right)", 25, 1e-6, 200);
    RooRealVar voigtian_sigma_omega = RooRealVar("voigtian_sigma_omega", "#sigma#left(m_{#omega}#right)", 25, 1e-6, 200);
    RooVoigtian voigtian_omega = RooVoigtian((string(varfit2_name)+"component_voigtian_omega").data(), "#omega Voigtian", varfit2, voigtian_mean_omega, voigtian_width_omega, voigtian_sigma_omega);

    // Breit Wigner
    RooRealVar gaussian_mean_omega = RooRealVar("gaussian_mean_omega", "#mu#left(m_{#omega}#right)", 782, minbin2, maxbin2);
    RooRealVar gaussian_width_omega = RooRealVar("gaussian_width_omega", "#sigma#left(m_{#omega}#right)", 25, 1e-6, 1e3);
    RooGaussian gaussian_omega = RooGaussian("gaussian_omega", "#omega Gaussian", varfit2, gaussian_mean_omega, gaussian_width_omega);

    // Gaussian
    RooRealVar gaussian_mean_b_jpsinonresonant = RooRealVar("gaussian_mean_b_jpsinonresonant", "gaussian_mean_b_jpsinonresonant", 3096, 1e3, 1e4);
    RooRealVar gaussian_width_b_jpsinonresonant = RooRealVar("gaussian_width_b_jpsinonresonant", "gaussian_width_b_jpsinonresonant", 250, 1e1, 1e3);
    RooGaussian gauss_b_jpsinonresonant = RooGaussian("gauss_b_jpsinonresonant", "gauss_b_jpsinonresonant", jpsi_mbconstrained_jpsicorrected, gaussian_mean_b_jpsinonresonant, gaussian_width_b_jpsinonresonant);
    
    // Gaussian
    RooRealVar gaussian_mean_jpsinonresonant = RooRealVar("gaussian_mean_jpsinonresonant", "gaussian_mean_jpsinonresonant", 5.28878e+03, 1e2, 1e4);
    RooRealVar gaussian_width2_jpsinonresonant = RooRealVar("gaussian_width2_jpsinonresonant", "gaussian_width2_jpsinonresonant", 5.86377e+01, 1e-6, 1e4);
    RooGaussian gauss_jpsinonresonant = RooGaussian("gauss_jpsinonresonant", "Signal gauss_jpsinonresonant", b_mjpsiconstrained_pi0corrected, gaussian_mean_jpsinonresonant, gaussian_width2_jpsinonresonant);

    // Gaussian
    RooRealVar gaussian_mean_omeganonresonant = RooRealVar("gaussian_mean_omeganonresonant", "gaussian_mean_omeganonresonant", 1000, 1e0, 1e5);
    RooRealVar gaussian_width2_omeganonresonant = RooRealVar("gaussian_width2_omeganonresonant", "gaussian_width2_omeganonresonant", 25, 1e-6, 1e5);
    RooGaussian gauss_omeganonresonant = RooGaussian("gauss_omeganonresonant", "#omega NonResonant Gaussian", varfit2, gaussian_mean_omeganonresonant, gaussian_width2_omeganonresonant);

    // B->Jpsi KPi crystal ball
    RooRealVar jpsikpi_mean = RooRealVar("jpsikpi_mean", "#mu", 5.27963e+03, minbin1, maxbin1);
    RooRealVar jpsikpi_sigma = RooRealVar("jpsikpi_sigma", "#sigma_{1}", 30, 0, 1e4);
    RooRealVar jpsikpi_n = RooRealVar("jpsikpi_n", "n_{1}", 0.794);
    RooRealVar jpsikpi_alpha = RooRealVar("jpsikpi_alpha", "#alpha_{1}", -1.129, -1e6, 1e6);
    RooCBShape jpsikpi_func = RooCBShape((string(varfit1_name)+"component_jpsikpi_func").data(), "component_jpsikpi_func", varfit1, jpsikpi_mean, jpsikpi_sigma, jpsikpi_alpha, jpsikpi_n);

    // B->Jpsi KPi second crystal ball
    RooRealVar jpsikpi_widthratio = RooRealVar("jpsikpi_widthratio", "#sigma Ratio", 1.493);
    RooFormulaVar jpsikpi_sigma2 = RooFormulaVar("jpsikpi_sigma2", "#sigma_{2}", "jpsikpi_widthratio*jpsikpi_sigma", RooArgList(jpsikpi_widthratio, jpsikpi_sigma));
    RooRealVar jpsikpi_n2 = RooRealVar("jpsikpi_n2", "n_{2}", 9.932);
    RooFormulaVar jpsikpi_alpha2 = RooFormulaVar("jpsikpi_alpha2", "#alpha_{2}", "-1*jpsikpi_alpha", RooArgList(jpsikpi_alpha));
    RooCBShape jpsikpi_func2 = RooCBShape((string(varfit1_name)+"component_jpsikpi_func2").data(), "jpsikpi_func2", varfit1, jpsikpi_mean, jpsikpi_sigma2, jpsikpi_alpha2, jpsikpi_n2);

    // bomegakpi KPi crystal ball
    RooRealVar bjpsikpi_mean = RooRealVar("bjpsikpi_mean", "#mu", 5280.784, minbin1, maxbin1);
    RooRealVar bjpsikpi_sigma = RooRealVar("bjpsikpi_sigma", "#sigma_{1}", 14.542, 0, maxbin1-minbin1);
    RooRealVar bjpsikpi_alpha = RooRealVar("bjpsikpi_alpha", "#alpha_{1}", 0.811);
    RooRealVar bjpsikpi_n = RooRealVar("bjpsikpi_n", "n_{1}", 3);
    RooCBShape bjpsikpi_func = RooCBShape((string(varfit1_name)+"component_bjpsikpi_func").data(), "#omega Crystal Ball", varfit1, bjpsikpi_mean, bjpsikpi_sigma, bjpsikpi_alpha, bjpsikpi_n);

    // bjpsikpi KPi second crystal ball
    RooRealVar bjpsikpi_widthratio = RooRealVar("bjpsikpi_widthratio", "#sigma Ratio", 1.034);
    RooFormulaVar bjpsikpi_sigma2 = RooFormulaVar("bjpsikpi_sigma2", "#sigma_{2}", "bjpsikpi_widthratio*bjpsikpi_sigma", RooArgList(bjpsikpi_widthratio, bjpsikpi_sigma));
    RooRealVar bjpsikpi_alpha2 = RooRealVar("bjpsikpi_alpha2", "#alpha_{2}", -0.887);
    RooRealVar bjpsikpi_n2 = RooRealVar("bjpsikpi_n2", "n_{2}", 3);
    RooCBShape bjpsikpi_func2 = RooCBShape((string(varfit1_name)+"component_bjpsikpi_func2").data(), "#omega Crystal Ball 2", varfit1, bjpsikpi_mean, bjpsikpi_sigma2, bjpsikpi_alpha2, bjpsikpi_n2);

    // JPsi->3Pi crystal ball
    RooRealVar jpsi_3pi_mean = RooRealVar("jpsi_3pi_mean", "#mu", 3096.9, minbin2, maxbin2);
    RooRealVar jpsi_3pi_sigma = RooRealVar("jpsi_3pi_sigma", "#sigma_{1}", 30, 0, maxbin2-minbin2);
    RooRealVar jpsi_3pi_n = RooRealVar("jpsi_3pi_n", "n_{1}", 3);
    RooRealVar jpsi_3pi_alpha = RooRealVar("jpsi_3pi_alpha", "#alpha_{1}", 0.027);
    RooCBShape jpsi_3pi_func = RooCBShape((string(varfit2_name)+"component_jpsi_3pi_func").data(), "Crystal Ball", varfit2, jpsi_3pi_mean, jpsi_3pi_sigma, jpsi_3pi_alpha, jpsi_3pi_n);

    // JPsi->3Pi second crystal ball
    RooRealVar jpsi_3pi_widthratio = RooRealVar("jpsi_3pi_widthratio", "#sigma Ratio", 0.794);
    RooFormulaVar jpsi_3pi_sigma2 = RooFormulaVar("jpsi_3pi_sigma2", "#sigma_{2}", "jpsi_3pi_widthratio*jpsi_3pi_sigma", RooArgList(jpsi_3pi_widthratio, jpsi_3pi_sigma));
    RooRealVar jpsi_3pi_n2 = RooRealVar("jpsi_3pi_n2", "n_{2}", 3);
    RooRealVar jpsi_3pi_alpha2 = RooRealVar("jpsi_3pi_alpha2", "#alpha_{2}", -0.848);
    RooCBShape jpsi_3pi_func2 = RooCBShape((string(varfit2_name)+"component_jpsi_3pi_func2").data(), "jpsi_3pi_func2", varfit2, jpsi_3pi_mean, jpsi_3pi_sigma2, jpsi_3pi_alpha2, jpsi_3pi_n2);

    // bomegakpi KPi crystal ball
    RooRealVar bomegakpi_mean = RooRealVar("bomegakpi_mean", "#mu", 5282, 0, 6000);
    RooRealVar bomegakpi_sigma = RooRealVar("bomegakpi_sigma", "#sigma_{1}", 54, 0, 1000);
    RooRealVar bomegakpi_alpha = RooRealVar("bomegakpi_alpha", "#alpha_{1}", 1.533);
    RooRealVar bomegakpi_n = RooRealVar("bomegakpi_n", "n_{1}", 3);
    RooCBShape bomegakpi_func = RooCBShape((string(varfit2_name)+"component_bomegakpi_func").data(), "#omega Crystal Ball", varfit1, bomegakpi_mean, bomegakpi_sigma, bomegakpi_alpha, bomegakpi_n);

    // bomegakpi KPi second crystal ball
    RooRealVar bomegakpi_widthratio = RooRealVar("bomegakpi_widthratio", "#sigma Ratio", 0.509);
    RooFormulaVar bomegakpi_sigma2 = RooFormulaVar("bomegakpi_sigma2", "#sigma_{2}", "bomegakpi_widthratio*bomegakpi_sigma", RooArgList(bomegakpi_widthratio, bomegakpi_sigma));
    RooRealVar bomegakpi_alpha2 = RooRealVar("bomegakpi_alpha2", "#alpha_{2}", -0.441);
    RooRealVar bomegakpi_n2 = RooRealVar("bomegakpi_n2", "n_{2}", 3);
    RooCBShape bomegakpi_func2 = RooCBShape((string(varfit2_name)+"component_bomegakpi_func2").data(), "#omega Crystal Ball 2", varfit1, bomegakpi_mean, bomegakpi_sigma2, bomegakpi_alpha2, bomegakpi_n2);

    // omega3pi KPi crystal ball
    RooRealVar omega3pi_mean = RooRealVar("omega3pi_mean", "#mu", 785, 0, 6000);
    RooRealVar omega3pi_sigma = RooRealVar("omega3pi_sigma", "#sigma_{1}", 44, 0, 1000);
    RooRealVar omega3pi_alpha = RooRealVar("omega3pi_alpha", "#alpha_{1}", 4.034);
    RooRealVar omega3pi_n = RooRealVar("omega3pi_n", "n_{1}", 3);
    RooCBShape omega3pi_func = RooCBShape((string(varfit2_name)+"component_omega3pi_func").data(), "#omega Crystal Ball", varfit2, omega3pi_mean, omega3pi_sigma, omega3pi_alpha, omega3pi_n);

    // omega3pi KPi second crystal ball
    RooRealVar omega3pi_widthratio = RooRealVar("omega3pi_widthratio", "#sigma Ratio", 0.3306);
    RooFormulaVar omega3pi_sigma2 = RooFormulaVar("omega3pi_sigma2", "#sigma_{2}", "omega3pi_widthratio*omega3pi_sigma", RooArgList(omega3pi_widthratio, omega3pi_sigma));
    RooRealVar omega3pi_alpha2 = RooRealVar("omega3pi_alpha2", "#alpha_{2}", -1.409);
    RooRealVar omega3pi_n2 = RooRealVar("omega3pi_n2", "n_{2}", 3);
    RooCBShape omega3pi_func2 = RooCBShape((string(varfit2_name)+"component_omega3pi_func2").data(), "#omega Crystal Ball 2", varfit2, omega3pi_mean, omega3pi_sigma2, omega3pi_alpha2, omega3pi_n2);

    // fractional yields
    // you need these and not absolute yields in combine
    // don"t fit with Extended!
    RooRealVar frac_sig1 = RooRealVar("frac_sig1", "f_{1}", 100, 0., dataentries);
    RooRealVar frac_sig2 = RooRealVar("frac_sig2", "f_{2}", 100, 0., dataentries);
    RooRealVar frac_jpsikpi1 = RooRealVar("frac_jpsikpi1", "frac_jpsikpi1", 1308.007);
    RooRealVar frac_jpsikpi2 = RooRealVar("frac_jpsikpi2", "frac_jpsikpi2", 500, 0., dataentries);
    RooRealVar frac_bomegakpi1 = RooRealVar("frac_frac_bomegakpi1", "frac_frac_bomegakpi1", 1000);
    RooRealVar frac_omega3pi1 = RooRealVar("frac_frac_omega3pi1", "frac_frac_omega3pi1", 1000);
    RooRealVar frac_omega3pi2 = RooRealVar("frac_frac_omega3pi2", "frac_frac_omega3pi2", 0.25);
    RooRealVar frac_bkg4 = RooRealVar("frac_bkg4", "frac_bkg4", 500, -1, dataentries);
    RooRealVar frac_bkg = RooRealVar("frac_bkg", "frac_bkg", 500, 0., dataentries);
    RooRealVar frac_bkg2 = RooRealVar("frac_bkg2", "frac_bkg2", 500, -1, dataentries);
    RooRealVar frac_bkg3 = RooRealVar("frac_bkg3", "frac_bkg3", 500, -1, dataentries);
    RooRealVar frac_bkg_jpsinonres = RooRealVar("frac_bkg_jpsinonres", "frac_bkg_jpsinonres", 500, 0., dataentries);
    // fixed to PDG (Jpsi K) / (Jpsi pi) value https://pdglive.lbl.gov/BranchingRatio.action?desig=14&parCode=S091
    Double_t frac_k_value = 0.079 / (1. + 0.079);
    RooRealVar frac_k = RooRealVar("frac_k", "frac_k", frac_k_value);

    RooAddPdf jpsikpi_twocb_mcshape = RooAddPdf((string(varfit1_name)+"component_jpsikpitwocbmc").data(), "Sum of the two B^{0} CB", RooArgList(jpsikpi_func, jpsikpi_func2), RooArgList(frac_jpsikpi1, frac_jpsikpi1));
    RooAddPdf omega3pi_twocb_mcshape = RooAddPdf((string(varfit2_name)+"component_omega3pitwocbmc").data(), "Sum of the two #omega CB", RooArgList(omega3pi_func, omega3pi_func2), RooArgList(0.25));
    RooAddPdf jpsi3pi_twocb_mcshape = RooAddPdf((string(varfit2_name)+"component_jpsi3pitwocbmc").data(), "Sum of the two #B JPsi3Pi", RooArgList(jpsi_3pi_func, jpsi_3pi_func2), RooArgList(0.061));
    RooAddPdf bomegakpi_twocb_mcshape = RooAddPdf((string(varfit1_name)+"component_Boomegakpitwocbmc").data(), "Sum of the two #B CB", RooArgList(bomegakpi_func, bomegakpi_func2), RooArgList(0.541));
    RooAddPdf bjpsikpi_twocb_mcshape = RooAddPdf((string(varfit1_name)+"component_Bjpsikpitwocbmc").data(), "Sum of the two #B CB", RooArgList(bjpsikpi_func, bjpsikpi_func2), RooArgList(0.446));

    // signal function
    RooProdPdf b_m_gaussian_cross_x_m_jpsigaussian = RooProdPdf(
        "b_m_gaussian_cross_x_m_jpsigaussian",
        "Gaussians",
        RooArgList(bjpsikpi_twocb_mcshape, jpsi3pi_twocb_mcshape));

    RooProdPdf b_m_gaussian_cross_x_m_jpsigaussian2 = RooProdPdf(
        "b_m_gaussian_cross_x_m_jpsigaussian2",
        "Gaussians2",
        RooArgList(jpsikpi_func2, signal_gauss2));

    RooProdPdf b_m_gaussian_cross_x_m_omegabreitwigner = RooProdPdf(
        "b_m_gaussian_cross_x_m_omegabreitwigner",
        "B^{0}Gaussian #times #omegaBreitWigner",
        RooArgList(bomegakpi_twocb_mcshape, omega3pi_twocb_mcshape));

    RooProdPdf b_m_gaussian_cross_x_m_omega3pi1 = RooProdPdf(
        "b_m_gaussian_cross_x_m_omega3pi1",
        "B^{0}Gaussian #times #omegaBreitWigner",
        RooArgList(signal_gauss, omega3pi_func));

    RooProdPdf b_m_gaussian_cross_x_m_omega3pi2 = RooProdPdf(
        "b_m_gaussian_cross_x_m_omega3pi2",
        "B^{0}Gaussian #times #omegaBreitWigner",
        RooArgList(signal_gauss, omega3pi_func2));

    RooProdPdf bkg_bm_cross_bkg_xm = RooProdPdf(
        "bkg_bm_cross_bkg_xm",
        "Prod_Bkgs",
        RooArgList(bkg1, bkg2));

    RooProdPdf b_jpsibkg_nonresonant = RooProdPdf(
        "b_jpsibkg_nonresonant",
        "b_jpsibkg_nonresonant",
        bjpsikpi_twocb_mcshape, bkg2);

    RooProdPdf jpsi_bbkg_nonresonant = RooProdPdf(
        "jpsi_bbkg_nonresonant",
        "jpsi_bbkg_nonresonant",
        gauss_b_jpsinonresonant, bkg1);

    RooProdPdf b_omegabkg_nonresonant = RooProdPdf(
        "b_omegabkg_nonresonant",
        "b_omegabkg_nonresonant",
        bomegakpi_twocb_mcshape, bkg2);

    RooProdPdf omega_bbkg_nonresonant = RooProdPdf(
        "omega_bbkg_nonresonant",
        "omega_bbkg_nonresonant",
        omega3pi_twocb_mcshape, bkg1);

    RooRealVar b_jpsikk_gauss_mean = RooRealVar("b_jpsikk_gauss_mean", "b_jpsikk_gauss_mean", 5398.557);
    RooRealVar b_jpsikk_gauss_width = RooRealVar("b_jpsikk_gauss_width", "b_jpsikk_gauss_width", 112.123);
    RooGaussian b_jpsikk_gauss = RooGaussian((string(varfit1_name)+"componentb_jpsikk_gauss").data(), "B^{c}#rightarrowJ/#psiKK (Gaussian)", varfit1, b_jpsikk_gauss_mean, b_jpsikk_gauss_width);

    RooRealVar b_jpsikk_cb_mean = RooRealVar("b_jpsikk_cb_mean", "b_jpsikk_cb_mean", 5332.403);
    RooRealVar b_jpsikk_cb_width = RooRealVar("b_jpsikk_cb_width", "b_jpsikk_cb_width", 23.454);
    RooRealVar b_jpsikk_cb_n = RooRealVar("b_jpsikk_cb_n", "b_jpsikk_cb_n", 5.540);
    RooRealVar b_jpsikk_cb_alpha = RooRealVar("b_jpsikk_cb_alpha", "b_jpsikk_cb_alpha", -0.454);
    RooCBShape b_jpsikk_cb = RooCBShape((string(varfit1_name)+"componentb_jpsikk_cb").data(), "B^{c}#rightarrowJ/#psiKK (CB)", varfit1, b_jpsikk_cb_mean, b_jpsikk_cb_width, b_jpsikk_cb_alpha, b_jpsikk_cb_n);

    RooAddPdf b_jpsikk_fitfunction = RooAddPdf(
        "b_jpsikk_fitfunction",
        "b_jpsikk_fitfunction",
        RooArgList(b_jpsikk_cb, b_jpsikk_gauss),
        RooArgList(RooConst(1034.205), RooConst(316.626))
    );

    RooAddPdf Total_Fit_function = RooAddPdf(
        "Total_Fit_function",
        "Fit function",
        RooArgList(b_m_gaussian_cross_x_m_jpsigaussian, bkg_bm_cross_bkg_xm, b_jpsibkg_nonresonant),
        RooArgList(frac_sig1, frac_bkg, frac_bkg2));
/*
    RooAddPdf Total_Fit_function = RooAddPdf(
        "Total_Fit_function",
        "Fit function",
        RooArgList(b_m_gaussian_cross_x_m_omegabreitwigner, bkg_bm_cross_bkg_xm, b_omegabkg_nonresonant, omega_bbkg_nonresonant),
        RooArgList(frac_sig1, frac_bkg, frac_bkg2, frac_bkg3));
*/
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // selection on data, plotting, fitting
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TH2F *fulldata_varfithist = new TH2F("", "", bins1, minbin1, maxbin1, bins2, minbin2, maxbin2);
    fulldata.fillHistogram(fulldata_varfithist, RooArgList(varfit1, varfit2));

    // plot
    TCanvas c1 = TCanvas("c1", "", 1920, 1080);
    c1.Draw();
    RooPlot *frame = varfit1.frame();
    frame->SetTitle(TString::Format("Fit of the %s variable", varfit1_title).Data());
    fulldata_varfithist->SetName(TString::Format("%s stats", varfit1_title));
    fulldata_varfithist->Draw();
    fulldata.plotOn(frame, Name("Data"), RooFit::Binning(bins1), RooFit::MarkerSize(1.5));
    gStyle->SetOptStat(analvar1_tobefit->variable_statoptions);

    // fit
    // results_data = fit_function.fitTo(fulldata, RooFit.Extended(True), RooFit.Save())
    RooFitResult *results_data = Total_Fit_function.fitTo(fulldata, RooFit::Save());

    Total_Fit_function.plotOn(frame, Name("Total_Fit_function"));
    Double_t chi2_datafit = frame->chiSquare();
    TIterator *fit_function_components = Total_Fit_function.getComponents()->createIterator();
    int i=1;
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit1_name)+"component");
        if (fitcomptoplot != string::npos)
        {
            Total_Fit_function.plotOn(frame, Name(fitcompit->GetName()), RooFit::Components(fitcompit->GetName()), 
                                    RooFit::LineStyle(kDashed), RooFit::LineColor(EColorPalette(i)));
            i++;
        }
    }
    frame->GetXaxis()->SetTitle(analvar1_tobefit->Xlabel().Data());
    frame->GetXaxis()->SetTitleSize(0.03);
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetYaxis()->SetTitleSize(0.03);
    frame->GetYaxis()->SetLabelSize(0.03);
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->Draw();
    c1.Update();

    Float_t padtextsize = analvar1_tobefit->variable_padtextsize;
    TLegend *leg1 = analvar1_tobefit->SetLegendPosAuto();
    /*leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(10);*/
    leg1->AddEntry(Total_Fit_function.GetName(), Total_Fit_function.GetTitle(), "L");
    fit_function_components = Total_Fit_function.getComponents()->createIterator();
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit1_name)+"component");
        if (fitcomptoplot != string::npos)
        {
         leg1->AddEntry(fitcompit->GetName(), fitcompit->GetTitle(), "L");
        }
    }
    leg1->AddEntry("Data", "Data in the MC sample", "EP");
    leg1->Draw("SAME");
    leg1->SetTextSize(padtextsize);
    c1.Update();

    TPaveStats *pvstat = analvar1_tobefit->SetStatAuto(fulldata_varfithist, leg1);
    pvstat->SetTextSize(padtextsize);
    pvstat->Draw("SAME");
    c1.Update();
    c1.Modified();

    // Printing parameters on the plot
    RooArgSet *params = Total_Fit_function.getParameters(varfit1);
    int params_size = params->getSize();
    Float_t pavetext_width = analvar1_tobefit->variable_textpadxlength, pavetext_entryheight = analvar1_tobefit->variable_textpadentryyheight;
    analvar1_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt = analvar1_tobefit->SetTextBoxHeaderAuto();
    pt->SetTextSize(padtextsize/2);
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->AddText(TString::Format("Fit function: %s", Total_Fit_function.getTitle().Data()))->SetTextAlign(22);    
    //pt->Draw("SAME");
    c1.Update();
    
    TPaveText *pt_left = analvar1_tobefit->SetTextBoxAuto(2, pt->GetX1NDC(), pt->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
    TPaveText *pt_right = analvar1_tobefit->SetTextBoxAuto(2, pt->GetX1NDC()+pavetext_width/2, pt->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
    pt_left->SetTextSize(0.018);
    pt_left->SetBorderSize(1);
    pt_left->SetFillColor(kWhite);
    pt_right->SetTextSize(0.018);
    pt_right->SetBorderSize(1);
    pt_right->SetFillColor(kWhite);
    pt_left->AddText("Parameter")->SetTextAlign(11);
    pt_left->AddLine();
    pt_right->AddText("Value")->SetTextAlign(31);
    pt_right->AddLine();
    TIterator *iter = params->createIterator();
    while(RooAbsArg* arg=(RooAbsArg*)iter->Next()) 
    { 
        pt_left->AddText(TString::Format("%s", arg->GetTitle()))->SetTextAlign(11); 
        pt_right->AddText(TString::Format("%.3f", ((RooRealVar*)arg)->getVal()))->SetTextAlign(31);
    }
    pt_left->AddText("#chi^{2}")->SetTextAlign(11);
    pt_right->AddText(TString::Format("%.3f", chi2_datafit))->SetTextAlign(31);
    //pt_left->Draw("SAME");
    //pt_right->Draw("SAME");
    c1.Update();
    TText* sample_data1 = new TText(0.1, 0.905, string("Dataset: "+string(varfit1_sample_data)).data());
    sample_data1->SetNDC(true);
    sample_data1->SetTextSize(0.025);
    sample_data1->SetTextAlign(10);  //align at top
    sample_data1->Draw("SAME");
    TText* sample_data1_2dprojection = new TText(0.805, 0.905, "2DFit Projection");
    sample_data1_2dprojection->SetBBoxX2(0.9);
    sample_data1_2dprojection->SetNDC(true);
    sample_data1_2dprojection->SetTextSize(0.025);
    sample_data1_2dprojection->SetTextAlign(10);  //align at top
    sample_data1_2dprojection->Draw("SAME");
    string outplot_filename = string(varfit1_name)+"_2DPROJECTION_"+varfit1_sample_data+"_FIT.png";
    c1.SaveAs(outplot_filename.data());

    // plot
    TCanvas c2 = TCanvas("c2", "", 1920, 1080);
    c2.Draw();
    RooPlot *frame2 = varfit2.frame();
    frame2->SetTitle(TString::Format("Fit of the %s variable", varfit2_title).Data());
    fulldata_varfithist->SetName(TString::Format("%s stats", varfit2_title));
    fulldata_varfithist->Draw();
    fulldata.plotOn(frame2, Name("Data"), RooFit::Binning(bins2), RooFit::MarkerSize(1.5));
    gStyle->SetOptStat(analvar2_tobefit->variable_statoptions);

    // fit
    // results_data = fit_function.fitTo(fulldata, RooFit.Extended(True), RooFit.Save())
    RooFitResult *results_data2 = Total_Fit_function.fitTo(fulldata, RooFit::Save());

    Total_Fit_function.plotOn(frame2, Name("Total_Fit_function"));
    Double_t chi2_datafit2 = frame2->chiSquare();

    fit_function_components = Total_Fit_function.getComponents()->createIterator();
    i = 1;
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit2_name)+"component_");
        if (fitcomptoplot != string::npos)
        {
            Total_Fit_function.plotOn(frame2, Name(fitcompit->GetName()), RooFit::Components(fitcompit->GetName()), 
                                       RooFit::LineStyle(kDashed), RooFit::LineColor(EColorPalette(i)));
            i++;
        }
    }
    frame2->GetXaxis()->SetTitle(analvar2_tobefit->Xlabel().Data());
    frame2->GetXaxis()->SetTitleSize(0.03);
    frame2->GetXaxis()->SetLabelSize(0.03);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetXaxis()->SetTitleOffset(1.3);
    frame2->Draw();
    c2.Update();

    Float_t padtextsize2 = analvar2_tobefit->variable_padtextsize;
    TLegend *leg2 = analvar2_tobefit->SetLegendPosAuto();
    /*leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetTextFont(10);*/
    leg2->AddEntry(Total_Fit_function.GetName(), Total_Fit_function.GetTitle(), "L");
    fit_function_components = Total_Fit_function.getComponents()->createIterator();
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit2_name)+"component_");
        if (fitcomptoplot != string::npos)
        {
            leg2->AddEntry(fitcompit->GetName(), fitcompit->GetTitle(), "L");
        }
    }
    leg2->AddEntry("Data", "Data in the MC sample", "EP");
    leg2->SetTextSize(padtextsize2);
    leg2->Draw("SAME");
    c2.Update();

    TPaveStats *pvstat2 = analvar2_tobefit->SetStatAuto(fulldata_varfithist, leg2);
    pvstat2->SetTextSize(padtextsize);
    pvstat2->Draw("SAME");
    c2.Update();
    c2.Modified();
    
    // Printing parameters on the plot
    RooArgSet *params2 = Total_Fit_function.getParameters(RooArgSet(varfit1, varfit2));
    int params_size2 = params2->getSize();
    Float_t pavetext_width2 = analvar2_tobefit->variable_textpadxlength, pavetext_entryheight2 = analvar2_tobefit->variable_textpadentryyheight;
    analvar2_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt2 = analvar2_tobefit->SetTextBoxHeaderAuto();
    pt2->SetTextSize(padtextsize/2);
    pt2->SetBorderSize(1);
    pt2->SetFillColor(kWhite);
    pt2->AddText(TString::Format("Fit function: %s", Total_Fit_function.getTitle().Data()))->SetTextAlign(22);    
    pt2->Draw("SAME");
    c2.Update();
    
    TPaveText *pt2_left = analvar2_tobefit->SetTextBoxAuto(2, pt2->GetX1NDC(), pt2->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
    TPaveText *pt2_right = analvar2_tobefit->SetTextBoxAuto(2, pt2->GetX1NDC()+pavetext_width/2, pt2->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
    pt2_left->SetTextSize(0.018);
    pt2_left->SetBorderSize(1);
    pt2_left->SetFillColor(kWhite);
    pt2_right->SetTextSize(0.018);
    pt2_right->SetBorderSize(1);
    pt2_right->SetFillColor(kWhite);
    pt2_left->AddText("Parameter")->SetTextAlign(11);
    pt2_left->AddLine();
    pt2_right->AddText("Value")->SetTextAlign(31);
    pt2_right->AddLine();
    TIterator *iter2 = params2->createIterator();
    while(RooAbsArg* arg=(RooAbsArg*)iter2->Next()) 
    { 
        pt2_left->AddText(TString::Format("%s", arg->GetTitle()))->SetTextAlign(11); 
        pt2_right->AddText(TString::Format("%f#pm%f", ((RooRealVar*)arg)->getVal(), ((RooRealVar*)arg)->getError()))->SetTextAlign(31);
    }
    pt2_left->AddText("#chi^{2}")->SetTextAlign(11);
    pt2_right->AddText(TString::Format("%.3f", chi2_datafit2))->SetTextAlign(31);
    pt2_left->Draw("SAME");
    pt2_right->Draw("SAME");
    c2.Modified();
    TText* sample_data2 = new TText(0.1, 0.905, string("Dataset: "+string(varfit2_sample_data)).data());
    sample_data2->SetNDC(true);
    sample_data2->SetTextSize(0.025);
    sample_data2->SetTextAlign(10);  //align at top
    sample_data2->Draw("SAME");
    TText* sample_data2_2dprojection = new TText(0.805, 0.905, "2DFit Projection");
    sample_data2_2dprojection->SetBBoxX2(0.9);
    sample_data2_2dprojection->SetNDC(true);
    sample_data2_2dprojection->SetTextSize(0.025);
    sample_data2_2dprojection->SetTextAlign(10);  //align at top
    sample_data2_2dprojection->Draw("SAME");
    string outplot_filename2 = string(varfit2_name)+"_2DPROJECTION_"+varfit2_sample_data+"_FIT.png";
    c2.SaveAs(outplot_filename2.data());

    // Printing parameters on the plot
    RooArgSet *params2d = Total_Fit_function.getParameters(RooArgList(varfit1, varfit2));

    TCanvas c2d = TCanvas("c2d", "", 1920, 1080);
    c2d.cd();
    TH2F *fit_2d_hist = new TH2F("fit_2d_hist", "fit_2d_hist", bins1, minbin1, maxbin1, bins2, minbin2, maxbin2);
    Total_Fit_function.fillHistogram(fit_2d_hist, RooArgList(varfit1, varfit2));    
    fulldata_varfithist->Draw("COLZ");
    fit_2d_hist->Draw("CONT2 SAME");
    gStyle->SetOptStat(false);
    string outplot_filename2d = string(varfit1_name)+"_"+string(varfit2_name)+"_"+varfit1_sample_data+"_FIT.png";
    c2d.SaveAs(outplot_filename2d.data());

    // 2D Pull plot
    TH2F* fulldata_varfithist_10bin = (TH2F*)fulldata_varfithist->Rebin2D(bins1/10, bins2/10);
    TH2F* pull_2D_histogram_10bins = (TH2F*)fit_2d_hist->Rebin2D(bins1/10, bins2/10);
	pull_2D_histogram_10bins->Scale(fulldata_varfithist_10bin->Integral()/pull_2D_histogram_10bins->Integral());
	for(int ibin=1; ibin<=10; ++ibin) {
		for(int jbin=1; jbin<=10; ++jbin) {
			double D = fulldata_varfithist_10bin->GetBinContent(ibin,jbin);
			double F = pull_2D_histogram_10bins->GetBinContent(ibin,jbin);
			double E = fulldata_varfithist_10bin->GetBinError(ibin,jbin);
			pull_2D_histogram_10bins->SetBinContent(ibin,jbin,(D-F)/E);
		}
	}

    Float_t target = 50.;
	int nBinAdapt = TMath::Nint(fulldata_varfithist_10bin->GetSumOfWeights()/target); //target 30 entries per bins
	if(nBinAdapt<16) nBinAdapt=16;//set a reasonable minimum
	AdaptBin binning("pullHist", nBinAdapt, minbin1, maxbin1, minbin2, maxbin2);
	std::vector<int> divisions_;
	TH2Poly* theHisto_ = new TH2Poly("name_", "", minbin1, maxbin1, minbin2, maxbin2);

	double x, y, n;
    int nBins(600);
	TH2* dataHist   = static_cast<TH2*>(fulldata.createHistogram("dataHist", varfit1, RooFit::Binning(nBins, minbin1, maxbin1),RooFit::YVar(varfit2,RooFit::Binning(nBins, minbin2, maxbin2))));
    TH2* modelHist  = static_cast<TH2*>(Total_Fit_function.createHistogram("fit_2d_hist", varfit1, RooFit::Binning(nBins, minbin1, maxbin1),RooFit::YVar(varfit2,RooFit::Binning(nBins, minbin2, maxbin2))));

    binning.loadDataFromHist(dataHist);
	TH2Poly* dataPoly = binning.getHisto("dataPoly");
	TH2Poly* modelPoly = binning.getHisto("modelPoly");

	//fill adaptive binning histograms
	for(uint ibin=1; ibin <= nBins; ++ibin) 
    {
		for(uint jbin=1; jbin <= nBins; ++jbin) 
        {
			x = dataHist->GetXaxis()->GetBinCenter(ibin);
			y = dataHist->GetYaxis()->GetBinCenter(jbin);
			n = dataHist->GetBinContent( ibin, jbin );
			dataPoly->Fill(x, y, n);
			x = modelHist->GetXaxis()->GetBinCenter(ibin);
			y = modelHist->GetYaxis()->GetBinCenter(jbin);
			n = modelHist->GetBinContent( ibin, jbin );
			modelPoly->Fill(x, y, n);
		}
	}
	for(int ibin=1; ibin<=dataPoly->GetNumberOfBins(); ++ibin) {
		double D = dataPoly->GetBinContent(ibin);
		double F = modelPoly->GetBinContent(ibin);
		double E = TMath::Sqrt(D);
		modelPoly->SetBinContent(ibin,(D-F)/E);
	}

    c2d.cd();
    c2d.SetLeftMargin(0.105);
    c2d.SetRightMargin(0.105);

	const Int_t NRGBs2 = 4;
	const Int_t NCont2 = 255;
	Double_t stops2[NRGBs2]  = { 0.00, 0.45, 0.55, 1.00};
	Double_t reds2[NRGBs2]   = { 0.00, 1.00, 1.00, 1.00};
	Double_t greens2[NRGBs2] = { 0.00, 1.00, 1.00, 0.00};
	Double_t blues2[NRGBs2]  = { 1.00, 1.00, 1.00, 0.00};
	TColor::CreateGradientColorTable(NRGBs2, stops2, reds2, greens2, blues2, NCont2);

    pull_2D_histogram_10bins->SetTitle("Pull Plot");
    pull_2D_histogram_10bins->GetXaxis()->SetTitle(analvar1_tobefit->Xlabel());
    pull_2D_histogram_10bins->GetXaxis()->SetTitleSize(0.025);
    pull_2D_histogram_10bins->GetXaxis()->SetTitleOffset(1.5);
    pull_2D_histogram_10bins->GetXaxis()->SetLabelSize(0.025);
    pull_2D_histogram_10bins->GetYaxis()->SetTitle(analvar2_tobefit->Xlabel());
    pull_2D_histogram_10bins->GetYaxis()->SetTitleSize(0.025);
    pull_2D_histogram_10bins->GetYaxis()->SetTitleOffset(1.25);
    pull_2D_histogram_10bins->GetYaxis()->SetLabelSize(0.025);
    pull_2D_histogram_10bins->GetZaxis()->SetTitle("#frac{Data-Model}{#sigma}");
    pull_2D_histogram_10bins->GetZaxis()->SetTitleSize(0.025);
    pull_2D_histogram_10bins->GetZaxis()->SetTitleOffset(0.8);
    pull_2D_histogram_10bins->GetZaxis()->SetLabelSize(0.025);
    pull_2D_histogram_10bins->Draw("colz");
    string outplot_filename2dpull = string(varfit1_name)+"_"+string(varfit2_name)+"_"+varfit1_sample_data+"_FITPULL.png";
    c2d.SaveAs(outplot_filename2dpull.data());

    cout << "Fit to data integral " << fulldata.numEntries() * (1 - frac_bkg.getVal()) << endl;
    cout << "Fit to data integral " << frac_sig1.getVal() << endl;
    cout << "Chi2 of the data fit to " << varfit1_name << " : " << chi2_datafit << endl;
    cout << "Chi2 of the data fit to " << varfit2_name << " : " << chi2_datafit2 << endl;

    c2d.cd();
    c2d.Clear();
    modelPoly->SetTitle(TString::Format("Pull Plot with Polygonal Bins (%.0f entries per bin)", target).Data());
    modelPoly->GetXaxis()->SetTitle(analvar1_tobefit->Xlabel());
    modelPoly->GetXaxis()->SetTitleSize(0.025);
    modelPoly->GetXaxis()->SetTitleOffset(1.5);
    modelPoly->GetXaxis()->SetLabelSize(0.025);
    modelPoly->GetYaxis()->SetTitle(analvar2_tobefit->Xlabel());
    modelPoly->GetYaxis()->SetTitleSize(0.025);
    modelPoly->GetYaxis()->SetTitleOffset(1.25);
    modelPoly->GetYaxis()->SetLabelSize(0.025);
    modelPoly->GetZaxis()->SetTitle("#frac{Data-Model}{#sigma}");
    modelPoly->GetZaxis()->SetTitleSize(0.025);
    modelPoly->GetZaxis()->SetTitleOffset(0.8);
    modelPoly->GetZaxis()->SetLabelSize(0.025);
	modelPoly->GetEntries();//TODO this fixes an issue where the histogram doesn't plot - I've given up trying to figure out why this works
	modelPoly->SetMinimum(-4.);
	modelPoly->SetMaximum( 4.);
	modelPoly->Draw("colz");
    c2d.SaveAs(("varbin_"+outplot_filename2dpull).data());

    return 0;
}