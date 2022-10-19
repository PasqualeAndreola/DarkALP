#include "HeaderFiles/Fitter.h"

using namespace std;
using namespace RooFit;

int FitterMass(ROOT::RDataFrame *input_dataframe,
               AnalysisVariable *analvar_tobefit,
               bool debug,
               string closuretest)
{
    // Silencing output
    RooFit::PrintLevel(-1);
    RooMinimizer::PrintLevel(RooMinimizer::None);

    // Variable name
    const char *varfit_name = analvar_tobefit->variable_name;
    const char *varfit_dimension = analvar_tobefit->variable_dimension;

    // Fit variable parameters
    Double_t bins = analvar_tobefit->variable_bins, minbin = analvar_tobefit->variable_histmin, maxbin = analvar_tobefit->variable_histmax, bin_width = (maxbin - minbin) / bins;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      Variables and PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector<RooRealVar*> roofitinputvar;
    RooRealVar mass = RooRealVar("b_mcor", "B corrected mass", 6.0, 6.6, "GeV");
    roofitinputvar.push_back(&mass);
    RooRealVar mvacut_kbdt1 = RooRealVar("mvacut_kbdt1", "BDT1 MVA cut", 0, 1, "");
    roofitinputvar.push_back(&mvacut_kbdt1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      mass ranges
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Double_t fit_range_lo = 4500;
    Double_t mass_window_lo = 5275 - 150;
    Double_t mass_window_hi = 6275 + 150;
    Double_t fit_range_hi = 6500;

    mass.setRange("left", fit_range_lo, mass_window_lo);
    mass.setRange("right", mass_window_hi, fit_range_hi);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //    1  argpar      -1.95759e+00   1.03831e+01   2.11846e-03  -1.97032e-01
    //    2  broad_width   5.62194e-02   5.57457e-03   7.32458e-05  -1.09202e+00
    //    3  frac_bkg     4.20044e-01   7.16860e-02   1.56468e-04  -1.60601e-01
    //    4  frac_pi      6.31013e-01   6.77992e-02   2.64847e-04   2.65121e-01
    //    5  frac_sig     2.67041e-01   2.28339e-01   5.99349e-04  -4.84672e-01
    //    6  maxM         6.20639e+00   2.25169e-01   8.23578e-04   7.09100e-01
    //    7  narrow_mean         6.26774e+00   8.02151e-03   7.24866e-05   1.18543e-01
    //    8  narrow_width   2.44845e-02   4.83913e-03   3.78671e-04  -5.35545e-01
    //    9  p1          -5.23507e-02   1.16627e-01   4.07071e-06  -5.23507e-04
    //   10  sg           1.14919e-02   1.00958e-02   1.07686e-03   2.99617e+00

    // combinatorial background poly
    //RooRealVar pol_c1 = RooRealVar("pol_c1", "coefficient of x^0 term", -5.23507e-02, -100, 100);
    // pol_c2 = RooRealVar("pol_c2", "coefficient of x^1 term", 0.6, -10, 10)
    // pol_c3 = RooRealVar("pol_c3", "coefficient of x^2 term", 0.5, -10, 10)
    //RooChebychev bkg = RooChebychev("bkg_pol", "1st order poly", mass, RooArgList(pol_c1));
    // bkg = RooChebychev("bkg_pol", "2nd order poly", mass, RooArgList(pol_c1, pol_c2))
    // bkg = RooChebychev("bkg_pol", "3rd order poly", mass, RooArgList(pol_c1, pol_c2, pol_c3))

    // expo
    RooRealVar slope = RooRealVar("slope", "slope", -0.001, -1e6, 1e6);
    RooExponential bkg = RooExponential("bkg_expo", "bkg_expo", mass, slope);

    // argus function, partially reconstructed decays
    RooRealVar argpar = RooRealVar("argpar", "argus shape parameter", -1.95759e+00, -20, 20);
    RooRealVar maxM = RooRealVar("maxM", "argus max m", 6.20639e+00, 6.0, 6.275); //6.2)
    RooArgusBG argus = RooArgusBG("argus", "Argus PDF", mass, maxM, argpar);

    // detector response function
    RooRealVar mg = RooRealVar("mg", "mg", 0);                         //, -0.01, 0.01)
    RooRealVar sg = RooRealVar("sg", "sg", 1.14919e-02, 0.0001, 0.03); //, 0.001,0.2)
    RooGaussian resGauss = RooGaussian("resGauss", "resGauss", mass, mg, sg);
    // construct convolution
    mass.setBins(10000, "fft");
    RooFFTConvPdf lxg = RooFFTConvPdf("lxg", "argus (X) gauss", mass, argus, resGauss);

    // B->Jpsi K crystal ball
    RooRealVar jpsik_mean = RooRealVar("jpsik_mean", "narrow_mean", 6.17, 6.10, 6.25);
    RooRealVar jpsik_sigma = RooRealVar("jpsik_sigma", "sigma", 0.03, 0.01, 0.1);
    RooRealVar jpsik_n     = RooRealVar("jpsik_n"    , "jpsik_n"    , 0.1 , 0.01,   3.  );
    RooRealVar jpsik_alpha = RooRealVar("jpsik_alpha", "jpsik_alpha", 2   ,  0.1,   4.  );
    RooCBShape jpsik_func = RooCBShape("jpsik_func", "jpsik_func", mass, jpsik_mean, jpsik_sigma, jpsik_alpha, jpsik_n);

    // signal narrow gaussian
    RooRealVar narrow_mean = RooRealVar("narrow_mean", "narrow_mean", 6.26774e+00, 6.1, 6.4);
    RooRealVar narrow_width = RooRealVar("narrow_width", "narrow_width", 2.44845e-02, 0., 0.1);
    RooGaussian narrow_gaus = RooGaussian("sig_narrow_gaus", "sig_narrow_gaus", mass, narrow_mean, narrow_width);

    // signal broad gaussian
    RooRealVar broad_mean = RooRealVar("broad_mean", "broad_mean", 6.27e+00, 6.1, 6.4);
    RooRealVar broad_width = RooRealVar("broad_width", "broad_width", 5.62194e-02, 0., 1.);
    RooGaussian broad_gaus = RooGaussian("sig_broad_gaus", "sig_broad_gaus", mass, broad_mean, broad_width);

    // absolute yields
    RooRealVar nsig = RooRealVar("signal_yield", "signal_yield", 800, 0., 10000.);
    RooRealVar nsig_narrow = RooRealVar("signal_yield_narrow", "signal_yield_narrow", 700, 0., 10000.);
    RooRealVar nsig_broad = RooRealVar("signal_yield_broad", "signal_yield_broad", 100, 0., 10000.);
    RooRealVar nbkgtot = RooRealVar("nbkgtot", "nbkgtot", 2000, 0., 10000.);
    RooRealVar nbkg = RooRealVar("nbkg", "nbkg", 7000, 0., 10000.);
    RooRealVar nPi = RooRealVar("nPi", "nPi", 1000, 0., 10000.);
    RooRealVar nK = RooRealVar("nK", "nK", 200, 0., 10000.);

    // fractional yields
    // you need these and not absolute yields in combine
    // don"t fit with Extended!
    RooRealVar frac_sig = RooRealVar("frac_sig", "frac_sig", 2.67041e-01, 0., 1.);
    RooRealVar frac_pi = RooRealVar("frac_pi", "frac_pi", 6.31013e-01, 0., 1.);
    RooRealVar frac_bkg = RooRealVar("frac_bkg", "frac_bkg", 4.20044e-01, 0., 1.);
    // fixed to PDG (Jpsi K) / (Jpsi pi) value https://pdglive.lbl.gov/BranchingRatio.action?desig=14&parCode=S091
    Double_t frac_k_value = 0.079 / (1. + 0.079);
    RooRealVar frac_k = RooRealVar("frac_k", "frac_k", frac_k_value);

    // signal function
    RooAddPdf signal_fit_function = RooAddPdf(
        "signal_fit_function",
        "signal_fit_function",
        RooArgList(narrow_gaus, broad_gaus),
        RooArgList(frac_sig));

    // signal Jpsi pi plus Jpsi K
    // RooAddPdf::pi_plus_k_fit_function[ frac_k * jpsik_func + [%] * signal_fit_function ]
    RooAddPdf pi_plus_k_fit_function = RooAddPdf(
        "pi_plus_k_fit_function",
        "pi_plus_k_fit_function",
        RooArgList(jpsik_func, signal_fit_function), // order matters for coefficients in next line https://www.nikhef.nl/~vcroft/SignalAndBackground-CompositeModels.html
        RooArgList(frac_k));

    // background function
    RooAddPdf bkg_fit_function = RooAddPdf(
        "bkg_fit_function",
        "bkg_fit_function",
        //     RooArgList(bkg, lxg, jpsik_func),
        //     RooArgList(frac_pi, frac_k)
        RooArgList(lxg, bkg),
        RooArgList(frac_pi));

    // total function
    RooAddPdf fit_function = RooAddPdf(
        "fit_function",
        "fit_function",
        RooArgList(bkg_fit_function, pi_plus_k_fit_function),
        RooArgList(frac_bkg));
/*
    // MC signal narrow gaussian
    RooRealVar mc_narrow_mean = RooRealVar("mc_narrow_mean", "mc_narrow_mean", 6.275, 5.5, 7.);
    RooRealVar mc_narrow_width = RooRealVar("mc_narrow_width", "mc_narrow_width", 0.038, 0., 1.);
    RooGaussian mc_narrow_gaus = RooGaussian("mc_sig_narrow_gaus", "mc_sig_narrow_gaus", mass, mc_narrow_mean, mc_narrow_width);

    // MC signal broad gaussian
    RooRealVar mc_broad_mean = RooRealVar("mc_broad_mean", "mc_broad_mean", 6.275, 5.5, 7.);
    RooRealVar mc_broad_width = RooRealVar("mc_broad_width", "mc_broad_width", 0.06, 0., 1.);
    RooGaussian mc_broad_gaus = RooGaussian("mc_sig_broad_gaus", "mc_sig_broad_gaus", mass, mc_broad_mean, mc_broad_width);

    RooRealVar mc_nsig = RooRealVar("mc_signal_yield", "mc_signal_yield", 800, 0, 100000);
    RooRealVar mc_nsig_narrow = RooRealVar("mc_signal_yield_narrow", "mc_signal_yield_narrow", 700, 0, 100000);
    RooRealVar mc_nsig_broad = RooRealVar("mc_signal_yield_broad", "mc_signal_yield_broad", 100, 0, 100000);

    // MC signal function
    RooAddPdf mc_signal_fitFunction = RooAddPdf(
        "mc_signal_fit_function",
        "mc_signal_fit_function",
        RooArgList(mc_narrow_gaus, mc_broad_gaus),
        RooArgList(mc_nsig_narrow, mc_nsig_broad));
*/
    RooArgSet thevars = RooArgSet();
    for (vector<RooRealVar*>::iterator rooiter=roofitinputvar.begin(); rooiter!=roofitinputvar.end(); rooiter++)
        thevars.add(**rooiter);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // selection on data, plotting, fitting
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // HLT_DoubleMu4_JpsiTrk_Displaced_v14
    const char *selection = "b_mcor > 4.6 && "
                            "b_mcor < 6.6 && "
                            "mvacut_kbdt1 > 0.5" ;

    // add gen matching
    const char *selection_mc = "x_true_pid==443";

    const char *data_files = "OutputFiles/KPi3Pi_SIM_JPsi.root";
    TTree *data_tree = (TTree *)(TFile::Open(data_files)->Get("T"));

    RooDataSet fulldata = RooDataSet("data", "data", data_tree->GetTree(), thevars, selection);

    // plot
    TCanvas c1 = TCanvas("c1", "", 1920, 1080);
    c1.Draw();
    RooPlot *frame = mass.frame();
    frame->SetTitle("");
    Int_t nbins = 80;
    fulldata.plotOn(frame, Name("Data"), RooFit::Binning(nbins), RooFit::MarkerSize(1.5));

    // fit
    // results_data = fit_function.fitTo(fulldata, RooFit.Extended(True), RooFit.Save())
    RooFitResult *results_data = fit_function.fitTo(fulldata, RooFit::Save());

    fit_function.plotOn(frame);
    Double_t chi2_datafit = frame->chiSquare();
    fit_function.plotOn(frame, Name("bkg_pol"), RooFit::Components("bkg_pol"), RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue));
    fit_function.plotOn(frame, Name("lxg"), RooFit::Components("lxg"), RooFit::LineStyle(kDashed), RooFit::LineColor(kOrange));
    fit_function.plotOn(frame, Name("signal_fit_function"), RooFit::Components("signal_fit_function"), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    fit_function.plotOn(frame, Name("jpsik_func"), RooFit::Components("jpsik_func"), RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen));

    frame->Draw();
    //CMS_lumi(c1, 4, 0, cmsText = "CMS", extraText = "   Preliminary", lumi_13TeV = "60 fb^{-1}");

    TLegend leg = TLegend(0.58, .65, .90, .90);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.035);
    // leg.SetNColumns(3);

    // RooFit
    leg.AddEntry("bkg_pol", "Combinatorial bkg");
    leg.AddEntry("lxg", "B_{c}#rightarrowJ/#Psi#pi + X", "L");
    leg.AddEntry("jpsik_func", "B_{c}#rightarrowJ/#PsiK", "L");
    leg.AddEntry("signal_fit_function", "B_{c}#rightarrowJ/#Psi#pi", "L");
    leg.AddEntry("Data", "Observed", "EP");
    leg.Draw("SAME");

    // gPad.SaveAs("sideband_fit.pdf")
    c1.SaveAs("KPi3Pi_JPsi_SIM_MassFit.png");



    RooArgSet *params = fit_function.getParameters(mass);
    params->writeToStream(cout, false);

    // Create a new empty workspace
    RooWorkspace *bchybridworkspace = new RooWorkspace("bchybridworkspace", "bchybridworkspace");

    // Import model and all its components into the workspace
    bchybridworkspace->import(fit_function);


    // Import data into the workspace
    bchybridworkspace->import(fulldata);

    // Print workspace contents
    bchybridworkspace->Print();

    // Save the workspace into a ROOT file
    bchybridworkspace->writeToFile("KPi3Pi_JPsiSim_Workspace.root");

    cout << "Fit to data integral " << fulldata.numEntries() * (1 - frac_bkg.getVal()) << endl;
    cout << "Fit to data integral " << frac_sig.getVal() << endl;
    cout << "Chi2 of the data fit: " <<  chi2_datafit << endl;
  
    return 0;
}