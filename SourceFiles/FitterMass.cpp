#include "HeaderFiles/Fitter.h"

using namespace std;
using namespace RooFit;

int FitterMass(pair<string, string> *input_file_tree,
               AnalysisVariable *analvar_tobefit,
               string selection_cut,
               string mc_selection_cut,
               bool debug,
               string closuretest)
{
    // Silencing output
    RooFit::PrintLevel(-1);
    RooMinimizer::PrintLevel(RooMinimizer::None);

    // Variable name
    const char *varfit_name = analvar_tobefit->variable_name;
    const char *varfit_title = analvar_tobefit->variable_prettyname;
    const char *varfit_dimension = analvar_tobefit->variable_dimension;
    const char *varfit_sample_data = analvar_tobefit->sampledata.data();
    // Fit variable parameters
    Double_t bins = analvar_tobefit->variable_bins, minbin = analvar_tobefit->variable_histmin, maxbin = analvar_tobefit->variable_histmax, bin_width = (maxbin - minbin) / bins;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      Variables and PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector<RooRealVar *> roofitinputvar;
    RooRealVar varfit = RooRealVar(varfit_name, varfit_title, minbin, maxbin, varfit_dimension);
    roofitinputvar.push_back(&varfit);
    RooRealVar mvacut_kbdt1 = RooRealVar("mvacut_kbdt1", "BDT1 MVA cut", 0, 1, "");
    roofitinputvar.push_back(&mvacut_kbdt1);
    RooRealVar x_m = RooRealVar("x_m", "Mass of the children cut", 0, 1e6, "");
    roofitinputvar.push_back(&x_m);
    RooRealVar x_m_pi0constrained = RooRealVar("x_m_pi0constrained", "Mass of the children constraining the pi0", 0, 1e6, "");
    roofitinputvar.push_back(&x_m_pi0constrained);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      varfit ranges
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Double_t fit_range_lo = minbin;
    Double_t mass_window_lo = (minbin+maxbin)/2 - 150;
    Double_t mass_window_hi = (minbin+maxbin)/2 + 150;
    Double_t fit_range_hi = maxbin;

    varfit.setRange("left", fit_range_lo, mass_window_lo);
    varfit.setRange("right", mass_window_hi, fit_range_hi);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // combinatorial background poly
    RooRealVar pol_c1 = RooRealVar("pol_c1", "coefficient of x^0 term", 100, -100, 1000);
    RooRealVar pol_c2 = RooRealVar("pol_c2", "coefficient of x^1 term", -5-02, -100, -0.);
    RooPolynomial polynomial_bkg = RooPolynomial("polybkg", "polybkg", varfit, RooArgList(pol_c1, pol_c2));
    // pol_c2 = RooRealVar("pol_c2", "coefficient of x^1 term", 0.6, -10, 10)
    // pol_c3 = RooRealVar("pol_c3", "coefficient of x^2 term", 0.5, -10, 10)
    // RooChebychev bkg = RooChebychev("bkg_pol", "1st order poly", varfit, RooArgList(pol_c1));
    // bkg = RooChebychev("bkg_pol", "2nd order poly", varfit, RooArgList(pol_c1, pol_c2))
    // bkg = RooChebychev("bkg_pol", "3rd order poly", varfit, RooArgList(pol_c1, pol_c2, pol_c3))

    // expo
    RooRealVar slope = RooRealVar("slope", "slope", -3e-3, -1e3, 1e3);
    RooExponential bkg = RooExponential("bkg_expo", "bkg_expo", varfit, slope);

    // B->Jpsi KPi crystal ball
    RooRealVar jpsikpi_mean = RooRealVar("jpsikpi_mean", "#mu", 5297, 0, 6000);
    RooRealVar jpsikpi_sigma = RooRealVar("jpsikpi_sigma", "#sigma_{1}", 30, 0, 1000);
    RooRealVar jpsikpi_n = RooRealVar("jpsikpi_n", "n_{1}", 50, 0., 1e6);
    RooRealVar jpsikpi_alpha = RooRealVar("jpsikpi_alpha", "#alpha_{1}", 2, -10, 10.);
    RooCBShape jpsikpi_func = RooCBShape("jpsikpi_func", "jpsikpi_func", varfit, jpsikpi_mean, jpsikpi_sigma, jpsikpi_alpha, jpsikpi_n);

    // B->Jpsi KPi second crystal ball
    RooRealVar jpsikpi_sigma2 = RooRealVar("jpsikpi_sigma2", "#sigma_{2}", 30, 0, 1000);
    RooRealVar jpsikpi_n2 = RooRealVar("jpsikpi_n2", "n_{2}", 1.77, 0., 500.);
    RooRealVar jpsikpi_alpha2 = RooRealVar("jpsikpi_alpha2", "#alpha_{2}", -2, -4, 4.);
    RooCBShape jpsikpi_func2 = RooCBShape("jpsikpi_func2", "jpsikpi_func2", varfit, jpsikpi_mean, jpsikpi_sigma2, jpsikpi_alpha2, jpsikpi_n2);

    // Gaussian
    RooRealVar gaussian_mean = RooRealVar("gaussian_mean", "gaussian_mean", 5279, 5100, 5300);
    RooRealVar gaussian_width = RooRealVar("gaussian_width", "gaussian_width", 25, 1e-6, 100);
    RooGaussian signal_gauss = RooGaussian("Signal_Gaus", "Signal Gaus", varfit, gaussian_mean, gaussian_width);

    // Gaussian
    RooRealVar gaussian_mean2 = RooRealVar("gaussian_mean2", "gaussian_mean2", 5450, 5300, 5600);
    RooRealVar gaussian_width2 = RooRealVar("gaussian_width2", "gaussian_width2", 25, 1e-6, 100);
    RooGaussian signal_gauss2 = RooGaussian("Signal_Gaus2", "Signal Gaus2", varfit, gaussian_mean2, gaussian_width2);

    // fractional yields
    // you need these and not absolute yields in combine
    // don"t fit with Extended!
    RooRealVar frac_sig = RooRealVar("frac_sig", "f_{1}", 0.15, 0., 1.);
    RooRealVar frac_pi = RooRealVar("frac_pi", "frac_pi", 6.31013e-01, 0., 1.);
    RooRealVar frac_bkg = RooRealVar("frac_bkg", "frac_bkg", 0.7, 0., 1.);
    // fixed to PDG (Jpsi K) / (Jpsi pi) value https://pdglive.lbl.gov/BranchingRatio.action?desig=14&parCode=S091
    Double_t frac_k_value = 0.079 / (1. + 0.079);
    RooRealVar frac_k = RooRealVar("frac_k", "frac_k", frac_k_value);

    // signal function
    RooAddPdf bkg_fit_function = RooAddPdf(
        "bkg_fit_function",
        "Expo",
        RooArgList(bkg),
        RooArgList(frac_bkg));

    RooAddPdf signal_fit_function = RooAddPdf(
        "signal_fit_function",
        "Gaussian+Expo",
        RooArgList(signal_gauss, bkg),
        RooArgList(frac_sig));

    // signal Jpsi pi plus Jpsi K
    // RooAddPdf::pi_plus_k_fit_function[ frac_k * jpsikpi_func + [%] * signal_fit_function ]
    /*RooAddPdf pi_plus_k_fit_function = RooAddPdf(
        "pi_plus_k_fit_function",
        "pi_plus_k_fit_function",
        RooArgList(jpsikpi_func, signal_fit_function), // order matters for coefficients in next line https://www.nikhef.nl/~vcroft/SignalAndBackground-CompositeModels.html
        RooArgList(frac_k));*/

    /*
        // MC signal narrow gaussian
        RooRealVar mc_narrow_mean = RooRealVar("mc_narrow_mean", "mc_narrow_mean", 6.275, 5.5, 7.);
        RooRealVar mc_narrow_width = RooRealVar("mc_narrow_width", "mc_narrow_width", 0.038, 0., 1.);
        RooGaussian mc_narrow_gaus = RooGaussian("mc_sig_narrow_gaus", "mc_sig_narrow_gaus", varfit, mc_narrow_mean, mc_narrow_width);

        // MC signal broad gaussian
        RooRealVar mc_broad_mean = RooRealVar("mc_broad_mean", "mc_broad_mean", 6.275, 5.5, 7.);
        RooRealVar mc_broad_width = RooRealVar("mc_broad_width", "mc_broad_width", 0.06, 0., 1.);
        RooGaussian mc_broad_gaus = RooGaussian("mc_sig_broad_gaus", "mc_sig_broad_gaus", varfit, mc_broad_mean, mc_broad_width);

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
    for (vector<RooRealVar *>::iterator rooiter = roofitinputvar.begin(); rooiter != roofitinputvar.end(); rooiter++)
        thevars.add(**rooiter);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // selection on data, plotting, fitting
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const char *data_file_name = input_file_tree->first.data();
    const char *data_tree_name = input_file_tree->second.data();
    TTree *data_tree = (TTree *)(TFile::Open(data_file_name)->Get(data_tree_name));
    RooDataSet fulldata = RooDataSet("data", "data", data_tree->GetTree(), thevars, selection_cut.data());
    TH1F *fulldata_varfithist = new TH1F("", "", bins, minbin, maxbin);
    fulldata.fillHistogram(fulldata_varfithist, varfit);

    // plot
    TCanvas c1 = TCanvas("c1", "", 1920, 1080);
    c1.Draw();
    RooPlot *frame = varfit.frame();
    frame->SetTitle(TString::Format("Fit of the %s variable", varfit_title).Data());
    Int_t nbins = bins;
    fulldata_varfithist->SetName(TString::Format("%s stats", varfit_title));
    fulldata_varfithist->Draw();
    fulldata.plotOn(frame, Name("Data"), RooFit::Binning(nbins), RooFit::MarkerSize(1.5));
    gStyle->SetOptStat(analvar_tobefit->variable_statoptions);

    // fit
    // results_data = fit_function.fitTo(fulldata, RooFit.Extended(True), RooFit.Save())
    RooFitResult *results_data = signal_fit_function.fitTo(fulldata, RooFit::Save());

    signal_fit_function.plotOn(frame, Name("signal_fit_function"));
    Double_t chi2_datafit = frame->chiSquare();
    signal_fit_function.plotOn(frame, Name("signal_gauss"), RooFit::Components(signal_gauss.GetName()), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
    frame->GetXaxis()->SetTitle(analvar_tobefit->Xlabel().Data());
    frame->GetXaxis()->SetTitleSize(0.03);
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetYaxis()->SetTitleSize(0.03);
    frame->GetYaxis()->SetLabelSize(0.03);
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->Draw();
    c1.Update();
    // CMS_lumi(c1, 4, 0, cmsText = "CMS", extraText = "   Preliminary", lumi_13TeV = "60 fb^{-1}");

    Float_t padtextsize = analvar_tobefit->variable_padtextsize;
    TLegend *leg = analvar_tobefit->SetLegendPosAuto();
    /*leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(10);*/
    leg->SetTextSize(padtextsize);
    leg->AddEntry("signal_fit_function", "Fit function", "L");
    leg->AddEntry("signal_gauss", "B^{0}#rightarrowJ/#PsiK#pi", "L");
    leg->AddEntry("Data", "Data in the MC sample", "EP");
    leg->Draw("SAME");
    c1.Update();

    TPaveStats *pvstat = analvar_tobefit->SetStatAuto(fulldata_varfithist, leg);
    pvstat->SetTextSize(padtextsize);
    pvstat->Draw("SAME");
    c1.Update();

    // Printing parameters on the plot
    RooArgSet *params = signal_fit_function.getParameters(varfit);
    params->writeToStream(cout, false);
    int params_size = params->getSize();
    Float_t pavetext_width = analvar_tobefit->variable_textpadxlength, pavetext_entryheight = analvar_tobefit->variable_textpadentryyheight;
    analvar_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt = analvar_tobefit->SetTextBoxHeaderAuto();
    pt->SetTextSize(padtextsize/2);
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->AddText(TString::Format("Fit function: %s", signal_fit_function.getTitle().Data()))->SetTextAlign(22);    
    pt->Draw("SAME");
    c1.Update();
    
    TPaveText *pt_left = analvar_tobefit->SetTextBoxAuto(2, pt->GetX1NDC(), pt->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
    TPaveText *pt_right = analvar_tobefit->SetTextBoxAuto(2, pt->GetX1NDC()+pavetext_width/2, pt->GetY1NDC(), params_size+2, pavetext_width/2, pavetext_entryheight);
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
    pt_left->Draw("SAME");
    pt_right->Draw("SAME");
    c1.Modified();
    TText* sample_data = new TText(0.1, 0.905, string("Dataset: "+string(varfit_sample_data)).data());
    sample_data->SetNDC(true);
    sample_data->SetTextSize(0.025);
    sample_data->SetTextAlign(10);  //align at top
    sample_data->Draw("SAME");
    string outplot_filename = string(varfit_name)+"_"+varfit_sample_data+"_FIT.png";
    c1.SaveAs(outplot_filename.data());
    // Create a new empty workspace
    RooWorkspace *bchybridworkspace = new RooWorkspace("bchybridworkspace", "bchybridworkspace");

    // Import model and all its components into the workspace
    bchybridworkspace->import(signal_fit_function);

    // Import data into the workspace
    bchybridworkspace->import(fulldata);

    // Print workspace contents
    bchybridworkspace->Print();

    // Save the workspace into a ROOT file
    bchybridworkspace->writeToFile("KPi3Pi_JPsiSim_Workspace.root");

    cout << "Fit to data integral " << fulldata.numEntries() * (1 - frac_bkg.getVal()) << endl;
    cout << "Fit to data integral " << frac_sig.getVal() << endl;
    cout << "Chi2 of the data fit: " << chi2_datafit << endl;

    return 0;
}