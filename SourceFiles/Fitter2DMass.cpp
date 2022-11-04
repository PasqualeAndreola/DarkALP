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
    RooRealVar x_m_pi0constrained = RooRealVar("x_m_pi0constrained", "Mass of the children constraining the pi0", 0, 1e6, "");
    roofitinputvar.push_back(&x_m_pi0constrained);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      varfit1 ranges
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Double_t fit_range_lo = minbin1;
    Double_t mass_window_lo = (minbin1+maxbin1)/2 - 150;
    Double_t mass_window_hi = (minbin1+maxbin1)/2 + 150;
    Double_t fit_range_hi = maxbin1;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //      PDFs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // combinatorial background poly
    RooRealVar pol_c1 = RooRealVar("pol_c1", "coefficient of x^0 term", 100, -100, 1000);
    RooRealVar pol_c2 = RooRealVar("pol_c2", "coefficient of x^1 term", -5-02, -100, -0.);
    RooPolynomial polynomial_bkg = RooPolynomial("polybkg", "polybkg", varfit1, RooArgList(pol_c1, pol_c2));
    RooRealVar pol2_c1 = RooRealVar("pol2_c1", "coefficient of x^0 term", 100, -100, 1000);
    RooRealVar pol2_c2 = RooRealVar("pol2_c2", "coefficient of x^1 term", -5-02, -100, -0.);
    RooPolynomial polynomial2_bkg = RooPolynomial("polybkg", "polybkg", varfit2, RooArgList(pol2_c1, pol2_c2));
    // pol_c2 = RooRealVar("pol_c2", "coefficient of x^1 term", 0.6, -10, 10)
    // pol_c3 = RooRealVar("pol_c3", "coefficient of x^2 term", 0.5, -10, 10)
    // RooChebychev bkg = RooChebychev("bkg_pol", "1st order poly", varfit1, RooArgList(pol_c1));
    // bkg = RooChebychev("bkg_pol", "2nd order poly", varfit1, RooArgList(pol_c1, pol_c2))
    // bkg = RooChebychev("bkg_pol", "3rd order poly", varfit1, RooArgList(pol_c1, pol_c2, pol_c3))

    // expo
    RooRealVar slope1 = RooRealVar("slope1", "slope1", -3e-3, -1e3, 1e3);
    RooExponential bkg1 = RooExponential("bkg_expo1", "bkg_expo1", varfit1, slope1);

    RooRealVar slope2 = RooRealVar("slope2", "slope2", -0.096, -1e2, 0.1);
    RooExponential bkg2 = RooExponential("bkg_expo2", "bkg_expo2", varfit2, slope2);

    // B->Jpsi KPi crystal ball
    RooRealVar jpsikpi_mean = RooRealVar("jpsikpi_mean", "#mu", 5297, 0, 6000);
    RooRealVar jpsikpi_sigma = RooRealVar("jpsikpi_sigma", "#sigma_{1}", 30, 0, 1000);
    RooRealVar jpsikpi_n = RooRealVar("jpsikpi_n", "n_{1}", 50, 0., 1e6);
    RooRealVar jpsikpi_alpha = RooRealVar("jpsikpi_alpha", "#alpha_{1}", 2, -10, 10.);
    RooCBShape jpsikpi_func = RooCBShape("jpsikpi_func", "jpsikpi_func", varfit1, jpsikpi_mean, jpsikpi_sigma, jpsikpi_alpha, jpsikpi_n);

    // B->Jpsi KPi second crystal ball
    RooRealVar jpsikpi_mean2 = RooRealVar("jpsikpi_mean2", "#mu_{2}", 3096, 3000, 3200);
    RooRealVar jpsikpi_sigma2 = RooRealVar("jpsikpi_sigma2", "#sigma_{2}", 30, 0, 1000);
    RooRealVar jpsikpi_n2 = RooRealVar("jpsikpi_n2", "n_{2}", 1.77, 0., 500.);
    RooRealVar jpsikpi_alpha2 = RooRealVar("jpsikpi_alpha2", "#alpha_{2}", -2, -4, 4.);
    RooCBShape jpsikpi_func2 = RooCBShape("jpsikpi_func2", "jpsikpi_func2", varfit2, jpsikpi_mean2, jpsikpi_sigma2, jpsikpi_alpha2, jpsikpi_n2);

    // Gaussian
    RooRealVar gaussian_mean = RooRealVar("gaussian_mean", "gaussian_mean", 5279, 5100, 5300);
    RooRealVar gaussian_width = RooRealVar("gaussian_width", "gaussian_width", 25, 1e-6, 100);
    RooGaussian signal_gauss = RooGaussian("Signal_Gaus", "Signal Gaus", varfit1, gaussian_mean, gaussian_width);

    // Gaussian
    RooRealVar gaussian_mean2 = RooRealVar("gaussian_mean2", "gaussian_mean2", 3096, 3000, 3200);
    RooRealVar gaussian_width2 = RooRealVar("gaussian_width2", "gaussian_width2", 25, 1e-6, 75);
    RooGaussian signal_gauss2 = RooGaussian("Signal_Gaus2", "Signal Gaus2", varfit2, gaussian_mean2, gaussian_width2);

    // fractional yields
    // you need these and not absolute yields in combine
    // don"t fit with Extended!
    RooRealVar frac_sig1 = RooRealVar("frac_sig1", "f_{1}", 0.15, 0., 1.);
    RooRealVar frac_sig2 = RooRealVar("frac_sig2", "f_{2}", 0.05, 0., 1.);
    RooRealVar frac_pi = RooRealVar("frac_pi", "frac_pi", 6.31013e-01, 0., 1.);
    RooRealVar frac_bkg = RooRealVar("frac_bkg", "frac_bkg", 0.7, 0., 1.);
    RooRealVar frac_bkg2 = RooRealVar("frac_bkg2", "frac_bkg2", 0.7, 0., 1.);
    // fixed to PDG (Jpsi K) / (Jpsi pi) value https://pdglive.lbl.gov/BranchingRatio.action?desig=14&parCode=S091
    Double_t frac_k_value = 0.079 / (1. + 0.079);
    RooRealVar frac_k = RooRealVar("frac_k", "frac_k", frac_k_value);

    // signal function
    RooAddPdf bkg_fit_function = RooAddPdf(
        "bkg_fit_function",
        "Expo",
        RooArgList(bkg1),
        RooArgList(frac_bkg));

    RooAddPdf signal_fit_function = RooAddPdf(
        "signal_fit_function",
        "Gaussian#left(B^{0}#right)+Gaussian#left(J/#psi#right)+Expo",
        RooArgList(signal_gauss, signal_gauss2, bkg_fit_function),
        RooArgList(frac_sig1, frac_sig2));

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
        RooGaussian mc_narrow_gaus = RooGaussian("mc_sig_narrow_gaus", "mc_sig_narrow_gaus", varfit1, mc_narrow_mean, mc_narrow_width);

        // MC signal broad gaussian
        RooRealVar mc_broad_mean = RooRealVar("mc_broad_mean", "mc_broad_mean", 6.275, 5.5, 7.);
        RooRealVar mc_broad_width = RooRealVar("mc_broad_width", "mc_broad_width", 0.06, 0., 1.);
        RooGaussian mc_broad_gaus = RooGaussian("mc_sig_broad_gaus", "mc_sig_broad_gaus", varfit1, mc_broad_mean, mc_broad_width);

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
    RooFitResult *results_data = signal_fit_function.fitTo(fulldata, RooFit::Save());

    signal_fit_function.plotOn(frame, Name("signal_fit_function"));
    Double_t chi2_datafit = frame->chiSquare();
    signal_fit_function.plotOn(frame, Name("signal_gauss"), RooFit::Components(signal_gauss.GetName()), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
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
    leg1->SetTextSize(padtextsize);
    leg1->AddEntry("signal_fit_function", "Fit function", "L");
    leg1->AddEntry("signal_gauss", "B^{0}#rightarrowJ/#PsiK#pi", "L");
    leg1->AddEntry("Data", "Data in the MC sample", "EP");
    leg1->Draw("SAME");
    c1.Update();

    TPaveStats *pvstat = analvar1_tobefit->SetStatAuto(fulldata_varfithist, leg1);
    pvstat->SetTextSize(padtextsize);
    pvstat->Draw("SAME");
    c1.Update();
    c1.Modified();

    // Printing parameters on the plot
    RooArgSet *params = signal_fit_function.getParameters(varfit1);
    params->writeToStream(cout, false);
    int params_size = params->getSize();
    Float_t pavetext_width = analvar1_tobefit->variable_textpadxlength, pavetext_entryheight = analvar1_tobefit->variable_textpadentryyheight;
    analvar1_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt = analvar1_tobefit->SetTextBoxHeaderAuto();
    pt->SetTextSize(padtextsize/2);
    pt->SetBorderSize(1);
    pt->SetFillColor(kWhite);
    pt->AddText(TString::Format("Fit function: %s", signal_fit_function.getTitle().Data()))->SetTextAlign(22);    
    pt->Draw("SAME");
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
    pt_left->Draw("SAME");
    pt_right->Draw("SAME");
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
    RooFitResult *results_data2 = signal_fit_function.fitTo(fulldata, RooFit::Save());

    signal_fit_function.plotOn(frame2, Name("signal_fit_function"));
    Double_t chi2_datafit2 = frame2->chiSquare();
    signal_fit_function.plotOn(frame2, Name("signal_gauss"), RooFit::Components(signal_gauss.GetName()), RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
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
    leg2->SetTextSize(padtextsize);
    leg2->AddEntry("signal_fit_function", "Fit function", "L");
    leg2->AddEntry("signal_gauss", "B^{0}#rightarrowJ/#PsiK#pi", "L");
    leg2->AddEntry("Data", "Data in the MC sample", "EP");
    leg2->Draw("SAME");
    c2.Update();

    TPaveStats *pvstat2 = analvar2_tobefit->SetStatAuto(fulldata_varfithist, leg2);
    pvstat2->SetTextSize(padtextsize);
    pvstat2->Draw("SAME");
    c2.Update();
    c2.Modified();
    
    // Printing parameters on the plot
    RooArgSet *params2 = signal_fit_function.getParameters(varfit2);
    params2->writeToStream(cout, false);
    int params_size2 = params2->getSize();
    Float_t pavetext_width2 = analvar2_tobefit->variable_textpadxlength, pavetext_entryheight2 = analvar2_tobefit->variable_textpadentryyheight;
    analvar2_tobefit->variable_headernumberofentries = params_size+2;
    TPaveText *pt2 = analvar2_tobefit->SetTextBoxHeaderAuto();
    pt2->SetTextSize(padtextsize/2);
    pt2->SetBorderSize(1);
    pt2->SetFillColor(kWhite);
    pt2->AddText(TString::Format("Fit function: %s", signal_fit_function.getTitle().Data()))->SetTextAlign(22);    
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
        pt2_right->AddText(TString::Format("%.3f", ((RooRealVar*)arg)->getVal()))->SetTextAlign(31);
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
    RooArgSet *params2d = signal_fit_function.getParameters(RooArgList(varfit1, varfit2));
    params2d->writeToStream(cout, false);

    TCanvas c2d = TCanvas("c2d", "", 1920, 1080);
    c2d.cd();
    //TH2* hh_pdf3 = signal_fit_function.createHistogram("llr", varfit1, Binning(bins1), YVar(varfit2, Binning(bins2)),Scaling(kFALSE)); 
    fulldata_varfithist->Draw("surf1");
    frame->Draw("SAME");
    string outplot_filename2d = string(varfit1_name)+"_"+string(varfit2_name)+"_"+varfit1_sample_data+"_FIT.png";
    c2d.SaveAs(outplot_filename2d.data());

    cout << "Fit to data integral " << fulldata.numEntries() * (1 - frac_bkg.getVal()) << endl;
    cout << "Fit to data integral " << frac_sig1.getVal() << endl;
    cout << "Chi2 of the data fit: " << chi2_datafit << endl;

    return 0;
}