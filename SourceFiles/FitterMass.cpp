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
    //RooRealVar mvacut_kbdt1 = RooRealVar("mvacut_kbdt1", "BDT1 MVA cut", 0, 1, "");
    //roofitinputvar.push_back(&mvacut_kbdt1);
    RooRealVar x_m = RooRealVar("x_m", "Mass of the children cut", 0, 1e6, "");
    roofitinputvar.push_back(&x_m);
    //RooRealVar x_m_pi0constrained = RooRealVar("x_cons_xpiz_m_best", "Mass of the children constraining the pi0", 0, 1e6, "");
    //roofitinputvar.push_back(&x_m_pi0constrained);
    //RooRealVar mvacut_kbdt1 = RooRealVar("mvacut_kbdt1", "BDT1 MVA cut", 0, 1, "");
    //roofitinputvar.push_back(&mvacut_kbdt1);
    RooRealVar b_cons_xpiz_fitstatus_best = RooRealVar("b_cons_xpiz_fitstatus_best", "b_cons_xpiz_fitstatus_best", 0, 4, "");
    roofitinputvar.push_back(&b_cons_xpiz_fitstatus_best);
    RooRealVar b_cons_Jpsixpiz_fitstatus_best = RooRealVar("b_cons_Jpsixpiz_fitstatus_best", "b_cons_Jpsixpiz_fitstatus_best", 0, 4, "");
    roofitinputvar.push_back(&b_cons_Jpsixpiz_fitstatus_best);
    RooRealVar b_cons_omegaxpiz_fitstatus_best = RooRealVar("b_cons_omegaxpiz_fitstatus_best", "b_cons_omegaxpiz_fitstatus_best", 0, 4, "");
    roofitinputvar.push_back(&b_cons_omegaxpiz_fitstatus_best);
    RooArgSet thevars = RooArgSet();
    for (vector<RooRealVar *>::iterator rooiter = roofitinputvar.begin(); rooiter != roofitinputvar.end(); rooiter++)
        thevars.add(**rooiter);

    const char *data_file_name = input_file_tree->first.data();
    const char *data_tree_name = input_file_tree->second.data();
    TTree *data_tree = (TTree *)(TFile::Open(data_file_name)->Get(data_tree_name));
    RooDataSet fulldata = RooDataSet("data", "data", data_tree->GetTree(), thevars, selection_cut.data());
    Int_t dataentries = fulldata.numEntries();
    TH1F *fulldata_varfithist = new TH1F("", "", bins, minbin, maxbin);
    fulldata.fillHistogram(fulldata_varfithist, varfit);

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
    RooRealVar jpsikpi_mean = RooRealVar("jpsikpi_mean", "#mu", 5279, minbin, maxbin);
    RooRealVar jpsikpi_sigma = RooRealVar("jpsikpi_sigma", "#sigma_{1}", 30, 0, maxbin-minbin);
    RooRealVar jpsikpi_n = RooRealVar("jpsikpi_n", "n_{1}", 3);
    RooRealVar jpsikpi_alpha = RooRealVar("jpsikpi_alpha", "#alpha_{1}", 1, 0., 6.);
    RooCBShape jpsikpi_func = RooCBShape((string(varfit_name)+"component_jpsikpi_func").data(), "Crystal Ball", varfit, jpsikpi_mean, jpsikpi_sigma, jpsikpi_alpha, jpsikpi_n);

    // B->Jpsi KPi second crystal ball
    RooRealVar jpsikpi_widthratio = RooRealVar("jpsikpi_widthratio", "#sigma Ratio", 1.1, 0., 1e6);
    RooFormulaVar jpsikpi_sigma2 = RooFormulaVar("jpsikpi_sigma2", "#sigma_{2}", "jpsikpi_widthratio*jpsikpi_sigma", RooArgList(jpsikpi_widthratio, jpsikpi_sigma));
    RooRealVar jpsikpi_n2 = RooRealVar("jpsikpi_n2", "n_{2}", 3);
    RooRealVar jpsikpi_alpha2 = RooRealVar("jpsikpi_alpha2", "#alpha_{2}", -1., -6., -0.);
    RooCBShape jpsikpi_func2 = RooCBShape((string(varfit_name)+"component_jpsikpi_func2").data(), "Crystal Ball 2", varfit, jpsikpi_mean, jpsikpi_sigma2, jpsikpi_alpha2, jpsikpi_n2);

    // JPsi->3Pi crystal ball
    RooRealVar jpsi_3pi_mean = RooRealVar("jpsi_3pi_mean", "#mu", 3096.9, minbin, maxbin);
    RooRealVar jpsi_3pi_sigma = RooRealVar("jpsi_3pi_sigma", "#sigma_{1}", 30, 0, maxbin-minbin);
    RooRealVar jpsi_3pi_n = RooRealVar("jpsi_3pi_n", "n_{1}", 3);
    RooRealVar jpsi_3pi_alpha = RooRealVar("jpsi_3pi_alpha", "#alpha_{1}", 1, 0., 6.);
    RooCBShape jpsi_3pi_func = RooCBShape((string(varfit_name)+"component_jpsi_3pi_func").data(), "Crystal Ball", varfit, jpsi_3pi_mean, jpsi_3pi_sigma, jpsi_3pi_alpha, jpsi_3pi_n);

    // JPsi->3Pi second crystal ball
    RooRealVar jpsi_3pi_widthratio = RooRealVar("jpsi_3pi_widthratio", "#sigma Ratio", 1.1, 0., 1e6);
    RooFormulaVar jpsi_3pi_sigma2 = RooFormulaVar("jpsi_3pi_sigma2", "#sigma_{2}", "jpsi_3pi_widthratio*jpsi_3pi_sigma", RooArgList(jpsi_3pi_widthratio, jpsi_3pi_sigma));
    RooRealVar jpsi_3pi_n2 = RooRealVar("jpsi_3pi_n2", "n_{2}", 3);
    RooRealVar jpsi_3pi_alpha2 = RooRealVar("jpsi_3pi_alpha2", "#alpha_{2}", -1., -6., -0.);
    RooCBShape jpsi_3pi_func2 = RooCBShape((string(varfit_name)+"component_jpsi_3pi_func2").data(), "Crystal Ball 2", varfit, jpsi_3pi_mean, jpsi_3pi_sigma2, jpsi_3pi_alpha2, jpsi_3pi_n2);

    // omega3pi KPi crystal ball
    RooRealVar omega3pi_mean = RooRealVar("omega3pi_mean", "#mu", 783.5, minbin, maxbin);
    RooRealVar omega3pi_sigma = RooRealVar("omega3pi_sigma", "#sigma_{1}", 30, 0, 1000);
    RooRealVar omega3pi_n = RooRealVar("omega3pi_n", "n_{1}", 3);
    RooRealVar omega3pi_alpha = RooRealVar("omega3pi_alpha", "#alpha_{1}", 2.3, 0, 6);
    RooCBShape omega3pi_func = RooCBShape((string(varfit_name)+"component_omega3pi_func").data(), "Crystal Ball", varfit, omega3pi_mean, omega3pi_sigma, omega3pi_alpha, omega3pi_n);

    // omega3pi KPi second crystal ball
    RooRealVar omega3pi_widthratio = RooRealVar("omega3pi_widthratio", "#sigma Ratio", 0.32, 0.01, 0.99);
    RooFormulaVar omega3pi_sigma2 = RooFormulaVar("omega3pi_sigma2", "#sigma_{2}", "omega3pi_widthratio*omega3pi_sigma", RooArgList(omega3pi_widthratio, omega3pi_sigma));
    RooRealVar omega3pi_n2 = RooRealVar("omega3pi_n2", "n_{2}", 3);
    RooRealVar omega3pi_alpha2 = RooRealVar("omega3pi_alpha2", "#alpha_{2}", -1.3, -6., 0.);
    RooCBShape omega3pi_func2 = RooCBShape((string(varfit_name)+"component_omega3pi_func2").data(), "Crystal Ball 2", varfit, omega3pi_mean, omega3pi_sigma2, omega3pi_alpha2, omega3pi_n2);

    // Gaussian
    RooRealVar gaussian_mean = RooRealVar("gaussian_mean", "gaussian_mean", 5279, minbin, maxbin);
    RooRealVar gaussian_width = RooRealVar("gaussian_width", "gaussian_width", 25, 1e-6, 1e3);
    RooGaussian signal_gauss = RooGaussian((string(varfit_name)+"component_Signal_Gaus").data(), "Signal Gaus", varfit, gaussian_mean, gaussian_width);

    // Gaussian
    RooRealVar bomegakpi_gaussian_widthratio = RooRealVar("bomegakpi_gaussian_widthratio", "#sigma Ratio", 1.1, 0.2, 50);
    RooFormulaVar bomegakpi_gaussian_sigma2 = RooFormulaVar("bomegakpi_gaussian__sigma2", "#sigma_{2}", "bomegakpi_gaussian_widthratio*gaussian_width", RooArgList(bomegakpi_gaussian_widthratio, gaussian_width));
    RooGaussian signal_gauss2 = RooGaussian((string(varfit_name)+"component_Signal_Gaus2").data(), "Signal Gaus2", varfit, gaussian_mean, bomegakpi_gaussian_sigma2);

    // Breit Wigner
    RooRealVar breit_mean_omega = RooRealVar("breit_mean_omega", "#mu#left(m_{#omega}#right)", 5350, minbin, maxbin);
    RooRealVar breit_width_omega = RooRealVar("breit_width_omega", "#sigma#left(m_{#omega}#right)", 25, 1e-6, 200);
    RooBreitWigner breit_wigner_omega = RooBreitWigner((string(varfit_name)+"component_breit_wigner_omega").data(), "#omega Breit Wigner", varfit, breit_mean_omega , breit_width_omega);

    // fractional yields
    // you need these and not absolute yields in combine
    // don"t fit with Extended!
    RooRealVar frac_sig = RooRealVar("frac_sig", "f_{1}", 2e4, -0, dataentries);
    RooRealVar frac_sig2 = RooRealVar("frac_sig2", "f_{2}", 1e4, -0, dataentries);
    RooRealVar frac_pi = RooRealVar("frac_pi", "frac_pi", 6.31013e-01, 0., dataentries);
    RooRealVar frac_bkg = RooRealVar("frac_bkg", "frac_bkg", 0.7, -0, dataentries);
    // fixed to PDG (Jpsi K) / (Jpsi pi) value https://pdglive.lbl.gov/BranchingRatio.action?desig=14&parCode=S091
    Double_t frac_k_value = 0.079 / (1. + 0.079);
    RooRealVar frac_k = RooRealVar("frac_k", "frac_k", frac_k_value);
/*
    // signal function
    RooAddPdf Total_Fit_Function = RooAddPdf(
        "Total_Fit_Function",
        "Total_Fit_Function",
        RooArgList(jpsikpi_func, jpsikpi_func2),
        RooArgList(frac_sig, frac_sig));
*/
    RooAddPdf Total_Fit_Function = RooAddPdf(
        "Total_Fit_Function",
        "Sum of Two Crystal Ball",
        RooArgList(omega3pi_func, omega3pi_func2),
        RooArgList(frac_sig,frac_sig2));

    // signal Jpsi pi plus Jpsi K
    // RooAddPdf::pi_plus_k_fit_function[ frac_k * jpsikpi_func + [%] * Total_Fit_Function ]
    /*RooAddPdf pi_plus_k_fit_function = RooAddPdf(
        "pi_plus_k_fit_function",
        "pi_plus_k_fit_function",
        RooArgList(jpsikpi_func, Total_Fit_Function), // order matters for coefficients in next line https://www.nikhef.nl/~vcroft/SignalAndBackground-CompositeModels.html
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // selection on data, plotting, fitting
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    RooFitResult *results_data = Total_Fit_Function.fitTo(fulldata, RooFit::Save());

    Total_Fit_Function.plotOn(frame, Name("Total_Fit_Function"));
    Double_t chi2_datafit = frame->chiSquare();
    TIterator *fit_function_components = Total_Fit_Function.getComponents()->createIterator();
    int i=1;
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit_name)+"component");
        if (fitcomptoplot != string::npos)
        {
            Total_Fit_Function.plotOn(frame, Name(fitcompit->GetName()), RooFit::Components(fitcompit->GetName()), 
                                    RooFit::LineStyle(kDashed), RooFit::LineColor(EColorPalette(i)));
            i++;
        }
    }
    frame->GetXaxis()->SetTitle(analvar_tobefit->Xlabel().Data());
    frame->GetXaxis()->SetTitleSize(0.03);
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetYaxis()->SetTitleSize(0.03);
    frame->GetYaxis()->SetLabelSize(0.03);
    frame->GetXaxis()->SetTitleOffset(1.3);
    frame->Draw();
    c1.Update();

    Float_t padtextsize = analvar_tobefit->variable_padtextsize;
    TLegend *leg = analvar_tobefit->SetLegendPosAuto();
    /*leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(10);*/
    leg->AddEntry(Total_Fit_Function.GetName(), Total_Fit_Function.GetTitle(), "L");
    fit_function_components = Total_Fit_Function.getComponents()->createIterator();
    while (RooAbsArg* fitcompit=(RooAbsArg*)fit_function_components->Next())
    {
        size_t fitcomptoplot = string(fitcompit->GetName()).find(string(varfit_name)+"component");
        if (fitcomptoplot != string::npos)
        {
         leg->AddEntry(fitcompit->GetName(), fitcompit->GetTitle(), "L");
        }
    }
    leg->AddEntry("Data", "Data in the MC sample", "EP");
    leg->Draw("SAME");
    leg->SetTextSize(padtextsize);
    c1.Update();


    TPaveStats *pvstat = analvar_tobefit->SetStatAuto(fulldata_varfithist, leg);
    pvstat->SetTextSize(padtextsize);
    pvstat->Draw("SAME");
    c1.Update();

    // Printing parameters on the plot
    RooArgSet *params = Total_Fit_Function.getParameters(varfit);
    params->writeToStream(cout, false);
    int params_size = params->getSize();
    Float_t pavetext_width = analvar_tobefit->variable_textpadxlength, pavetext_entryheight = analvar_tobefit->variable_textpadentryyheight;
    analvar_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt = analvar_tobefit->SetTextBoxHeaderAuto();
    pt->SetTextSize(padtextsize/2);
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->AddText(TString::Format("Fit function: %s", Total_Fit_Function.getTitle().Data()))->SetTextAlign(22);    
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
    bchybridworkspace->import(Total_Fit_Function);

    // Import data into the workspace
    bchybridworkspace->import(fulldata);

    // Print workspace contents
    bchybridworkspace->Print();

    // Save the workspace into a ROOT file
    bchybridworkspace->writeToFile("KPi3Pi_JPsiSim_Workspace.root");

    cout << "Fit to data integral " << fulldata.numEntries() * (1 - frac_bkg.getVal()) << endl;
    cout << "Fit to data integral " << frac_sig.getVal() << endl;
    cout << "Chi2 of the data fit: " << chi2_datafit << endl;
    cout << "frac sig 1 " << frac_sig.getVal() << endl;
    cout << "frac sig 2 " << frac_sig2.getVal() << endl;
    return 0;
}