/*!
 *  \file EfficienciesAnalyzer.cpp
 *  \brief Source file for \ref EfficienciesAnalyzer function implementation
 */

/*!
 *  \fn int EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                        vector<AnalysisVariable> *var_to_be_analyzed,
                        bool debug) "";
 *  \brief Function used to import selected variables from the files in input. Efficiencies of the cuts of the variables
 *         are computed and written in an output file defined by the user.
 * 
 *  \param data_holder Represents a list of pair of names: the first name of the pair represents the input file, 
 *                     whereas the second name of the pair represent the name of the tree that holds the variables.
 *  \param *var_to_be_analyzed Vector of instances of the \ref AnalysisVariable class.
 *  \param cuts String that expresses the cut being used in the analysis
 *  \param start Starting time of the main function (default is the starting time of the function)
 */
#include "HeaderFiles/VariablePlotter.h"
#include "HeaderFiles/PrintFuncInfo.h"

int EfficienciesAnalyzer(unordered_map<string, pair<string, string>> data_holder,
                        vector<AnalysisVariable> *var_to_be_analyzed,
                        string cuts,
                        chrono::_V2::system_clock::time_point start,
                        bool debug)
{
    // Enabling implicit Multi-threading
    ROOT::EnableImplicitMT();

    // Defining the quantities that will be used
    ROOT::RDataFrame *KPiEtaPiPi_SIM_100ps = NULL, *KPiEtaPiPi_SIM_1ps = NULL, *KPi3Pi_SIM_100ps = NULL, *KPi3Pi_SIM_1ps = NULL;
    ROOT::RDF::RResultPtr<TH1F> histcompare_norm, histcompare_sig, histcompare_signor, histcompare_bkg;
    TH1F histcompare_ratio;
    TCanvas canvas_comparison = TCanvas("canvas_comparison", "canvas_comparison", 1360, 768);
    TCanvas canvas_comparison_quadrati = TCanvas("canvas_comparison_quadrati", "canvas_comparison_quadrati", 1080, 1080);
    TCanvas canvas_comparison_ratio = TCanvas("canvas_comparison_ratio", "canvas_comparison_ratio", 1360, 768);

    TCanvas canvas_comparison_roc_sigvsnorm = TCanvas("canvas_comparison_roc_sigvsnorm", "canvas_comparison_roc_sigvsnorm", 1360, 768);
    TCanvas canvas_comparison_roc_sigvsbkg = TCanvas("canvas_comparison_roc_sigvsbkg", "canvas_comparison_roc_sigvsbkg", 1360, 768);
    TCanvas canvas_comparison_roc_normvsbkg = TCanvas("canvas_comparison_roc_normvsbkg", "canvas_comparison_roc_normvsbkg", 1360, 768);
    TCanvas canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg = TCanvas("canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg", "canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg", 1360, 768);
    TCanvas canvas_compraison_quadratiroc_sigvsnorm = TCanvas("canvas_compraison_quadratiroc_sigvsnorm", "canvas_compraison_quadratiroc_sigvsnorm", 1080, 1080);
    TCanvas canvas_compraison_quadratiroc_sigvsbkg = TCanvas("canvas_compraison_quadratiroc_sigvsbkg", "canvas_compraison_quadratiroc_sigvsbkg", 1080, 1080);
    TCanvas canvas_compraison_quadratiroc_normvsbkg = TCanvas("canvas_compraison_quadratiroc_normvsbkg", "canvas_compraison_quadratiroc_normvsbkg", 1080, 1080);
    TCanvas canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg = TCanvas("canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg", "canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg", 1080, 1080);
    TLegend legend_roc_sigvsnorm, legend_roc_sigvsbkg, legend_roc_normvsbkg, legend_roc_KPi3Pi_SIM_100psvsbkg;
    vector<TGraph> legend_roc_sigvsnorm_graph; vector<TString> legend_roc_sigvsnorm_entries;
    vector<TGraph> legend_roc_sigvsbkg_graph; vector<TString> legend_roc_sigvsbkg_entries;
    vector<TGraph> legend_roc_normvsbkg_graph; vector<TString> legend_roc_normvsbkg_entries;
    vector<TGraph> legend_roc_KPi3Pi_SIM_100psvsbkg_graph; vector<TString> legend_roc_KPi3Pi_SIM_100psvsbkg_entries;
    TMultiGraph roc_sig_norm, roc_sig_bkg, roc_norm_bkg, roc_KPi3Pi_SIM_100ps_bkg;

    TCanvas canvas_comparison_roc_inverted_sigvsnorm = TCanvas("canvas_comparison_roc_inverted_sigvsnorm", "canvas_comparison_roc_inverted_sigvsnorm", 1360, 768);
    TCanvas canvas_comparison_roc_inverted_sigvsbkg = TCanvas("canvas_comparison_roc_inverted_sigvsbkg", "canvas_comparison_roc_inverted_sigvsbkg", 1360, 768);
    TCanvas canvas_comparison_roc_inverted_normvsbkg = TCanvas("canvas_comparison_roc_inverted_normvsbkg", "canvas_comparison_roc_inverted_normvsbkg", 1360, 768);
    TCanvas canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg = TCanvas("canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg", "canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg", 1360, 768);
    TLegend legend_roc_inverted_sigvsnorm, legend_roc_inverted_sigvsbkg, legend_roc_inverted_normvsbkg, legend_roc_inverted_KPi3Pi_SIM_100psvsbkg;
    TMultiGraph roc_inverted_sig_norm, roc_inverted_sig_bkg, roc_inverted_norm_bkg, roc_inverted_KPi3Pi_SIM_100ps_bkg;

    TString roc_curve_title = TString();

    // Loading dataframes which hold the variables needed to compute the efficiencies
    bool fileexist = TFile(data_holder["KPiEtaPiPi_SIM_100ps"].first.data()).IsZombie();
    bool treeexist = TTree(data_holder["KPiEtaPiPi_SIM_100ps"].first.data(), data_holder["KPiEtaPiPi_SIM_100ps"].second.data()).IsZombie();
    if (fileexist == false && treeexist == false)
        KPiEtaPiPi_SIM_100ps = new ROOT::RDataFrame(data_holder["KPiEtaPiPi_SIM_100ps"].second.data(), data_holder["KPiEtaPiPi_SIM_100ps"].first.data());

    fileexist = TFile(data_holder["KPiEtaPiPi_SIM_1ps"].first.data()).IsZombie();
    treeexist = TTree(data_holder["KPiEtaPiPi_SIM_1ps"].first.data(), data_holder["KPiEtaPiPi_SIM_1ps"].second.data()).IsZombie();   
    if (fileexist == false && treeexist == false)
        KPiEtaPiPi_SIM_1ps = new ROOT::RDataFrame(data_holder["KPiEtaPiPi_SIM_1ps"].second.data(), data_holder["KPiEtaPiPi_SIM_1ps"].first.data());

    fileexist = TFile(data_holder["KPi3Pi_SIM_100ps"].first.data()).IsZombie();
    treeexist = TTree(data_holder["KPi3Pi_SIM_100ps"].first.data(), data_holder["KPi3Pi_SIM_100ps"].second.data()).IsZombie();     
    if (fileexist == false && treeexist == false)
        KPi3Pi_SIM_100ps = new ROOT::RDataFrame(data_holder["KPi3Pi_SIM_100ps"].second.data(), data_holder["KPi3Pi_SIM_100ps"].first.data());

    fileexist = TFile(data_holder["KPi3Pi_SIM_1ps"].first.data()).IsZombie();
    treeexist = TTree(data_holder["KPi3Pi_SIM_1ps"].first.data(), data_holder["KPi3Pi_SIM_1ps"].second.data()).IsZombie();  
    if (fileexist == false && treeexist == false)
        KPi3Pi_SIM_1ps = new ROOT::RDataFrame(data_holder["KPi3Pi_SIM_1ps"].second.data(), data_holder["KPi3Pi_SIM_1ps"].first.data());

    // Looping over chosen variables in order to compute the efficiencies
    for (vector<AnalysisVariable>::iterator analyzing = var_to_be_analyzed->begin(); analyzing != var_to_be_analyzed->end(); analyzing++)
    {
        Float_t bins = analyzing->variable_bins, min = analyzing->variable_histmin, max = analyzing->variable_histmax;
        if (KPiEtaPiPi_SIM_100ps->HasColumn(analyzing->variable_name) == false) continue;
        ElapsedTimeStamper(start);
        Int_t neventsprecut   = KPiEtaPiPi_SIM_100ps->Count().GetValue();
        ElapsedTimeStamper(start);
        Int_t neventsaftercut = KPiEtaPiPi_SIM_100ps->Filter(cuts).Count().GetValue();
        ElapsedTimeStamper(start);
        cout << "Number of events before the cuts: " << neventsprecut << endl;
        cout << "Number of events after the cuts: " << neventsaftercut << endl;
        histcompare_sig = KPiEtaPiPi_SIM_100ps->Fill(TH1F(TString::Format("%s_sig", analyzing->variable_prettyname), TString::Format("%s_sig", analyzing->variable_name),
                                            bins, min, max),
                                       {analyzing->variable_name});
                                       ElapsedTimeStamper(start);
        cout << "Number of events after the cuts: " << histcompare_sig->GetEntries() << endl;                                    
        if (KPiEtaPiPi_SIM_1ps->HasColumn(analyzing->variable_name) == false) continue;
        histcompare_norm = KPiEtaPiPi_SIM_1ps->Fill(TH1F(TString::Format("%s_norm", analyzing->variable_prettyname), TString::Format("%s_norm", analyzing->variable_name),
                                                    bins, min, max),
                                               {analyzing->variable_name});
        if (KPi3Pi_SIM_100ps->HasColumn(analyzing->variable_name) == false) continue;
        histcompare_signor = KPi3Pi_SIM_100ps->Fill(TH1F(TString::Format("%s_KPi3Pi_SIM_100ps", analyzing->variable_prettyname), TString::Format("%s_KPi3Pi_SIM_100ps", analyzing->variable_name),
                                                    bins, min, max),
                                            {analyzing->variable_name});
        if (KPi3Pi_SIM_1ps->HasColumn(analyzing->variable_name) == false) continue;            
        histcompare_bkg = KPi3Pi_SIM_1ps->Fill(TH1F(TString::Format("%s_bkg", analyzing->variable_prettyname), TString::Format("%s_bkg", analyzing->variable_name),
                                                bins, min, max),
                                           {analyzing->variable_name});

        // Computing efficiencies and purities
        histcompare_norm->Sumw2();
        histcompare_sig->Sumw2();
        histcompare_signor->Sumw2();
        histcompare_bkg->Sumw2();
        histcompare_norm->Scale(1 / histcompare_norm->Integral());
        histcompare_sig->Scale(1 / histcompare_sig->Integral());
        histcompare_signor->Scale(1 / histcompare_signor->Integral());
        histcompare_bkg->Scale(1 / histcompare_bkg->Integral());
        vector<Double_t> cuts, sig_eff, norm_eff, signor_eff, bkg_eff;
        vector<Double_t> purity_KPi3Pi_SIM_100ps, purity_KPi3Pi_SIM_100psbkg, purity_sigbkg, purity_normbkg;

        vector<Double_t> sig_eff_inverted, norm_eff_inverted, signor_eff_inverted, bkg_eff_inverted;
        vector<Double_t> purity_inverted_KPi3Pi_SIM_100ps, purity_inverted_KPi3Pi_SIM_100psbkg, purity_inverted_sigbkg, purity_inverted_normbkg;
        for (int i = 0; i < bins; i++)
        {
            cuts.push_back(min + i * (max - min) / bins);
            Double_t sig_integral = histcompare_sig->Integral(histcompare_sig->FindFixBin(cuts[i]), histcompare_sig->FindFixBin(max));
            Double_t norm_integral = histcompare_norm->Integral(histcompare_norm->FindFixBin(cuts[i]), histcompare_norm->FindFixBin(max));
            Double_t signor_integral = histcompare_signor->Integral(histcompare_signor->FindFixBin(cuts[i]), histcompare_signor->FindFixBin(max));
            Double_t bkg_integral = histcompare_bkg->Integral(histcompare_bkg->FindFixBin(cuts[i]), histcompare_bkg->FindFixBin(max));

            // Events at the right of the cut
            sig_eff.push_back(sig_integral);
            norm_eff.push_back(norm_integral);
            signor_eff.push_back(signor_integral);
            bkg_eff.push_back(bkg_integral);          
            purity_KPi3Pi_SIM_100ps.push_back(sig_integral / (sig_integral + 4.76 * norm_integral));
            purity_sigbkg.push_back(sig_integral / (sig_integral + bkg_integral));
            purity_normbkg.push_back(norm_integral / (norm_integral + bkg_integral));
            purity_KPi3Pi_SIM_100psbkg.push_back(signor_integral / (signor_integral + bkg_integral));

            // Events at the left of the cut
            sig_eff_inverted.push_back(1-sig_integral);
            norm_eff_inverted.push_back(1-norm_integral);
            signor_eff_inverted.push_back(1-signor_integral);
            bkg_eff_inverted.push_back(1-bkg_integral);  
            purity_inverted_KPi3Pi_SIM_100ps.push_back((1-sig_integral) / ((1-sig_integral) + 4.76 * (1-norm_integral)));
            purity_inverted_sigbkg.push_back((1-sig_integral) / ((1-sig_integral) + (1-bkg_integral)));
            purity_inverted_normbkg.push_back((1-norm_integral) / ((1-norm_integral) + (1-bkg_integral)));
            purity_inverted_KPi3Pi_SIM_100psbkg.push_back((1-signor_integral) / ((1-signor_integral) + (1-bkg_integral)));
        }

        // Plotting
        if (analyzing->variable_plot_flag == true)
        {
            THStack histcompare_stack = THStack("Histstack", "Histstack");

            canvas_comparison.cd();
            canvas_comparison.Clear();
            histcompare_signor->SetAxisRange(0, histcompare_signor->GetMaximum(), "Y");
            histcompare_signor->SetStats(false);
            histcompare_signor->SetFillStyle(1001);
            histcompare_signor->SetMarkerStyle(kFullSquare);
            histcompare_signor->SetMarkerSize(1);
            histcompare_stack.Add(histcompare_signor.GetPtr());
            canvas_comparison.Update();

            histcompare_bkg->SetAxisRange(0, histcompare_bkg->GetMaximum(), "Y");
            histcompare_bkg->SetStats(false);
            histcompare_bkg->SetFillStyle(1001);
            histcompare_bkg->SetMarkerStyle(kFullSquare);
            histcompare_bkg->SetMarkerSize(1);
            histcompare_stack.Add(histcompare_bkg.GetPtr());
            canvas_comparison.Update();

            histcompare_stack.SetTitle(TString::Format("%s pre and after cut", analyzing->variable_prettyname));
            histcompare_stack.Draw("NOSTACK PLC PMC");
            canvas_comparison.Update();
            histcompare_stack.GetXaxis()->SetTitle(analyzing->Xlabel());
            histcompare_stack.GetYaxis()->SetTitle("Normalized occurencies");

            //Adjusting the legend
            TLegend *legend = analyzing->SetLegendPosAuto("TR", 2);
            legend->AddEntry(histcompare_signor->GetName(), TString::Format("%s_Sig+Norm", analyzing->variable_prettyname), "PLC PMC");
            legend->AddEntry(histcompare_bkg->GetName(), TString::Format("%s_Bkg", analyzing->variable_prettyname), "PLC PMC");
            legend->SetTextSize(0.025);
            legend->Draw("SAME");
            if (analyzing->variable_logscale_flag == true) {histcompare_signor->SetMinimum(1e-4); histcompare_bkg->SetMinimum(1e-4);
                                                          histcompare_signor->SetMaximum(1e-1); histcompare_bkg->SetMaximum(1e-1); gPad->SetLogy();}
            canvas_comparison.Update();
            canvas_comparison.Print(TString::Format("%s/%s_distribution_KPi3Pi_SIM_100ps_bkg.png", analyzing->variable_histplotfolder, analyzing->variable_name));
            canvas_comparison.Clear();
            legend->Clear();
            gPad->SetLogy(0);
        }
        /*
        if (analyzing->flag_KPi3Pi_SIM_100pseffplot == true)
        {
            canvas_comparison.cd();
            canvas_comparison.Clear();
            TGraph eff_curve_sig = TGraph(sig_eff.size(), cuts.data(), sig_eff.data());
            eff_curve_sig.SetName("Signal efficiency");
            eff_curve_sig.SetTitle(TString::Format("analyzing %s: Signal - KPiEtaPiPi_SIM_1ps efficiencies", analyzing->analyzinglabel));
            eff_curve_sig.GetXaxis()->SetTitle(TString::Format("Cut Applied on %s[%s]", analyzing->analyzinglabel, analyzing->analyzingdimension));
            eff_curve_sig.GetYaxis()->SetTitle("Efficiency");
            eff_curve_sig.SetMarkerStyle(21);
            eff_curve_sig.SetMarkerColor(kBlue);
            eff_curve_sig.Draw("AP PLC");
            TGraph eff_curve_norm = TGraph(norm_eff.size(), cuts.data(), norm_eff.data());
            eff_curve_norm.SetName("KPiEtaPiPi_SIM_1ps Efficiency");
            eff_curve_norm.SetMarkerStyle(21);
            eff_curve_norm.SetMarkerColor(kRed);
            eff_curve_norm.Draw("P PLC SAME");

            //Adjusting the legend
            TLegend legend_roc_sigvsnorm = TLegend();
            legend_roc_sigvsnorm.SetHeader(analyzing->variable_name, "C");
            legend_roc_sigvsnorm.AddEntry(eff_curve_sig.GetName(), eff_curve_sig.GetName(), "PLC");
            legend_roc_sigvsnorm.AddEntry(eff_curve_norm.GetName(), eff_curve_norm.GetName(), "P PLC");
            legend_roc_sigvsnorm.Draw("same");
            gPad->Update();
            legend_roc_sigvsnorm.SetTextSize(0.019);
            legend_roc_sigvsnorm.SetX1NDC(0.75);
            legend_roc_sigvsnorm.SetX2NDC(0.9);
            legend_roc_sigvsnorm.SetY1NDC(0.78);
            legend_roc_sigvsnorm.SetY2NDC(0.9);
            gPad->Modified();
            canvas_comparison.Print(TString::Format("%s/%s_efficiency_sig_norm.png", analyzing->output_directory_png, analyzing->variable_name));
        }

        if (analyzing->flag_KPi3Pi_SIM_100psbkgplot == true)
        {
            THStack histcompare_stack = THStack("Histstack", "Histstack");

            canvas_comparison.cd();
            canvas_comparison.Clear();
            histcompare_norm->SetAxisRange(0, histcompare_norm->GetMaximum(), "Y");
            histcompare_norm->SetStats(false);
            histcompare_norm->SetFillStyle(1001);
            histcompare_norm->SetMarkerStyle(kFullSquare);
            histcompare_norm->SetMarkerSize(1);
            histcompare_stack.Add(histcompare_norm.GetPtr());
            canvas_comparison.Update();

            histcompare_sig->SetAxisRange(0, histcompare_sig->GetMaximum(), "Y");
            histcompare_sig->SetStats(false);
            histcompare_sig->SetFillStyle(1001);
            histcompare_sig->SetMarkerStyle(kFullSquare);
            histcompare_sig->SetMarkerSize(1);
            histcompare_stack.Add(histcompare_sig.GetPtr());
            canvas_comparison.Update();

            histcompare_bkg->SetAxisRange(0, histcompare_bkg->GetMaximum(), "Y");
            histcompare_bkg->SetStats(false);
            histcompare_bkg->SetFillStyle(1001);
            histcompare_bkg->SetMarkerStyle(kFullSquare);
            histcompare_bkg->SetMarkerSize(1);
            histcompare_stack.Add(histcompare_bkg.GetPtr());

            histcompare_stack.SetTitle(TString::Format("%s analyzing distibution: Signal, KPiEtaPiPi_SIM_1ps and KPi3Pi_SIM_1ps",analyzing->analyzinglabel));
            histcompare_stack.Draw("NOSTACK PLC PMC");
            canvas_comparison.Update();
            histcompare_stack.GetXaxis()->SetTitle(TString::Format("%s[%s]",analyzing->analyzinglabel, analyzing->analyzingdimension));
            histcompare_stack.GetYaxis()->SetTitle(analyzing->KPi3Pi_SIM_100psplot_ytitle);

            //Adjusting the legend
            TLegend *legend = new TLegend(0.1, 0.81, 0.31, 0.9);
            if (analyzing->distributionplotlegend != NULL) legend = (analyzing->distributionplotlegend);
            legend->AddEntry(histcompare_norm->GetName(), TString::Format("%s_Norm", analyzing->variable_prettyname), "PLC PMC");
            legend->AddEntry(histcompare_sig->GetName(), TString::Format("%s_Sig", analyzing->variable_prettyname), "PLC PMC");
            legend->AddEntry(histcompare_bkg->GetName(), TString::Format("%s_Bkg", analyzing->variable_prettyname), "PLC PMC");
            legend->SetTextSize(0.025);
            legend->Draw("SAME");
            if (analyzing->distributionlogscale != 0) {histcompare_sig->SetMinimum(1e-4); histcompare_norm->SetMinimum(1e-4); histcompare_bkg->SetMinimum(1e-4);
                                                          histcompare_sig->SetMaximum(1e-1); histcompare_norm->SetMaximum(1e-1); histcompare_bkg->SetMaximum(1e-1); gPad->SetLogy();}
            canvas_comparison.Update();
            canvas_comparison.Print(TString::Format("%s/%s_distribution_sig_norm_bkg.png", analyzing->output_directory_png, analyzing->variable_name));
            canvas_comparison.Clear();
            legend->Clear();
            gPad->SetLogy(0);

            canvas_comparison_quadrati.cd();
            canvas_comparison_quadrati.Clear();
            gPad->SetRightMargin(0.05);
            gPad->SetLeftMargin(0.15);
            histcompare_stack.Draw("NOSTACK PLC PMC");
            TLegend legend_quadrati = analyzing->SetLegendPosAuto(analyzing->legpos, 3, analyzing->legxlength);
            if (analyzing->distributionplotlegend != NULL) legend = (analyzing->distributionplotlegend);
            legend_quadrati.AddEntry(histcompare_norm->GetName(), TString::Format("%s_Norm", analyzing->variable_prettyname), "PLC PMC");
            legend_quadrati.AddEntry(histcompare_sig->GetName(), TString::Format("%s_Sig", analyzing->variable_prettyname), "PLC PMC");
            legend_quadrati.AddEntry(histcompare_bkg->GetName(), TString::Format("%s_Bkg", analyzing->variable_prettyname), "PLC PMC");
            legend_quadrati.SetTextSize(0.025);
            legend_quadrati.Draw("SAME");
            if (analyzing->distributionlogscale != 0) {histcompare_sig->SetMinimum(1e-4); histcompare_norm->SetMinimum(1e-4); histcompare_bkg->SetMinimum(1e-4);
                                                          histcompare_sig->SetMaximum(1e-1); histcompare_norm->SetMaximum(1e-1); histcompare_bkg->SetMaximum(1e-1); gPad->SetLogy();}
            canvas_comparison_quadrati.Update();
            canvas_comparison_quadrati.Print(TString::Format("%s/%s_distribution_sig_norm_bkg.png", TString::Format("%s_quadrati",analyzing->output_directory_png).Data(), analyzing->variable_name));
            canvas_comparison_quadrati.Clear();
            legend_quadrati.Clear();
            gPad->SetLogy(0);
        }

        if (analyzing->flag_KPi3Pi_SIM_100psbkgeffplot == true)
        {
            canvas_comparison.cd();
            canvas_comparison.Clear();
            TGraph eff_curve_signor = TGraph(signor_eff.size(), cuts.data(), signor_eff.data());
            eff_curve_signor.SetName("Sig+Norm efficiency");
            eff_curve_signor.SetTitle(TString::Format("analyzing %s: Signal+KPiEtaPiPi_SIM_1ps - KPi3Pi_SIM_1ps efficiencies", analyzing->analyzinglabel));
            eff_curve_signor.GetXaxis()->SetTitle(TString::Format("Cut Applied on %s[%s]", analyzing->analyzinglabel, analyzing->analyzingdimension));
            eff_curve_signor.GetYaxis()->SetTitle("Efficiency");
            eff_curve_signor.SetMarkerStyle(21);
            eff_curve_signor.SetMarkerColor(kBlue);
            eff_curve_signor.Draw("AP PLC");
            TGraph eff_curve_bkg = TGraph(bkg_eff.size(), cuts.data(), bkg_eff.data());
            eff_curve_bkg.SetName("KPi3Pi_SIM_1ps Efficiency");
            eff_curve_bkg.SetMarkerStyle(21);
            eff_curve_bkg.SetMarkerColor(kRed);
            eff_curve_bkg.Draw("P PLC SAME");

            //Adjusting the legend
            TLegend legend_roc_sigvsnorm = TLegend();
            legend_roc_sigvsnorm.SetHeader(analyzing->variable_name, "C");
            legend_roc_sigvsnorm.AddEntry(eff_curve_signor.GetName(), eff_curve_signor.GetName(), "PLC");
            legend_roc_sigvsnorm.AddEntry(eff_curve_bkg.GetName(), eff_curve_bkg.GetName(), "P PLC");
            legend_roc_sigvsnorm.Draw("same");
            gPad->Update();
            legend_roc_sigvsnorm.SetTextSize(0.019);
            legend_roc_sigvsnorm.SetX1NDC(0.75);
            legend_roc_sigvsnorm.SetX2NDC(0.9);
            legend_roc_sigvsnorm.SetY1NDC(0.78);
            legend_roc_sigvsnorm.SetY2NDC(0.9);
            gPad->Modified();
            canvas_comparison.Print(TString::Format("%s/%s_efficiency_KPi3Pi_SIM_100ps_bkg.png", analyzing->output_directory_png, analyzing->variable_name));
        }

        if (analyzing->flag_ratioplot == true)
        {
            canvas_comparison_ratio.cd();
            canvas_comparison_ratio.Clear();
            histcompare_ratio = (TH1F)(*histcompare_sig.GetPtr() / *histcompare_norm.GetPtr());

            histcompare_ratio.SetTitle(TString::Format("analyzing %s: Signal over KPiEtaPiPi_SIM_1ps", analyzing->analyzinglabel));
            histcompare_ratio.SetStats(kTRUE);
            histcompare_ratio.SetLineColorAlpha(kBlue, 1);
            histcompare_ratio.SetFillStyle(1001);
            histcompare_ratio.SetFillColor(kBlue);
            histcompare_ratio.SetFillColorAlpha(kBlue, 0.8);
            histcompare_ratio.SetMarkerStyle(kFullSquare);
            histcompare_ratio.SetMarkerColor(kBlue);
            histcompare_ratio.SetMarkerSize(1);
            histcompare_ratio.GetXaxis()->SetTitle(TString::Format("#frac{%s Sig}{%s KPiEtaPiPi_SIM_1ps}", analyzing->analyzinglabel, analyzing->analyzinglabel));
            histcompare_ratio.GetXaxis()->SetTitleOffset(1.7);
            histcompare_ratio.GetXaxis()->SetTitleSize(0.024);
            histcompare_ratio.GetYaxis()->SetTitle("Occurences");
            histcompare_ratio.SetAxisRange(min, max, "X");
            histcompare_ratio.SetAxisRange(0, histcompare_ratio.GetBinContent(histcompare_ratio.GetMaximumBin()), "Y");
            histcompare_ratio.Draw("HF P");
            canvas_comparison_ratio.Update();

            //Adjusting the stats
            TPaveStats *ratiostats = (TPaveStats *)histcompare_ratio.GetListOfFunctions()->FindObject("stats");
            ratiostats->SetName("Pisa Q^{2} distribution Ratio: Signal over KPiEtaPiPi_SIM_1ps");
            ratiostats->SetOptStat(112211);
            ratiostats->SetX1NDC(0.1);
            ratiostats->SetX2NDC(0.31);
            ratiostats->SetY1NDC(0.78);
            ratiostats->SetY2NDC(0.9);
            ratiostats->Draw("same");
            canvas_comparison_ratio.Print(TString::Format("%s/%s_ratio_sig_over_KPiEtaPiPi_SIM_1ps.png", analyzing->output_directory_png, analyzing->variable_name));

            // Q2 Signal over KPi3Pi_SIM_1ps
            canvas_comparison_ratio.cd();
            canvas_comparison_ratio.Clear();
            histcompare_ratio = (TH1F)(*histcompare_sig.GetPtr() / *histcompare_bkg.GetPtr());

            histcompare_ratio.SetTitle(TString::Format("analyzing %s: Signal over KPi3Pi_SIM_1ps", analyzing->analyzinglabel));
            histcompare_ratio.SetStats(kTRUE);
            histcompare_ratio.SetLineColorAlpha(kBlue, 1);
            histcompare_ratio.SetFillStyle(1001);
            histcompare_ratio.SetFillColor(kBlue);
            histcompare_ratio.SetFillColorAlpha(kBlue, 0.8);
            histcompare_ratio.SetMarkerStyle(kFullSquare);
            histcompare_ratio.SetMarkerColor(kBlue);
            histcompare_ratio.SetMarkerSize(1);
            histcompare_ratio.GetXaxis()->SetTitle(TString::Format("#frac{%s Sig}{%s KPi3Pi_SIM_1ps}", analyzing->analyzinglabel, analyzing->analyzinglabel));
            histcompare_ratio.GetXaxis()->SetTitleOffset(1.7);
            histcompare_ratio.GetXaxis()->SetTitleSize(0.024);  
            histcompare_ratio.GetYaxis()->SetTitle("Occurences");
            histcompare_ratio.SetAxisRange(min, max, "X");
            histcompare_ratio.SetAxisRange(0, histcompare_ratio.GetBinContent(histcompare_ratio.GetMaximumBin()), "Y");
            histcompare_ratio.Draw("HF P");
            canvas_comparison_ratio.Update();

            //Adjusting the stats
            ratiostats = (TPaveStats *)histcompare_ratio.GetListOfFunctions()->FindObject("stats");
            ratiostats->SetName("Pisa Q^{2} distribution Ratio: Signal over KPi3Pi_SIM_1ps");
            ratiostats->SetOptStat(112211);
            ratiostats->SetX1NDC(0.1);
            ratiostats->SetX2NDC(0.31);
            ratiostats->SetY1NDC(0.78);
            ratiostats->SetY2NDC(0.9);
            ratiostats->Draw("same");
            canvas_comparison_ratio.Print(TString::Format("%s/%s_ratio_sig_over_KPi3Pi_SIM_1ps.png", analyzing->output_directory_png, analyzing->variable_name));

            // Q2 KPiEtaPiPi_SIM_1ps over KPi3Pi_SIM_1ps
            canvas_comparison_ratio.cd();
            canvas_comparison_ratio.Clear();
            histcompare_ratio = (TH1F)(*histcompare_norm.GetPtr() / *histcompare_bkg.GetPtr());

            histcompare_ratio.SetTitle(TString::Format("analyzing %s: KPiEtaPiPi_SIM_1ps over KPi3Pi_SIM_1ps", analyzing->analyzinglabel));
            histcompare_ratio.SetStats(kTRUE);
            histcompare_ratio.SetLineColorAlpha(kBlue, 1);
            histcompare_ratio.SetFillStyle(1001);
            histcompare_ratio.SetFillColor(kBlue);
            histcompare_ratio.SetFillColorAlpha(kBlue, 0.8);
            histcompare_ratio.SetMarkerStyle(kFullSquare);
            histcompare_ratio.SetMarkerColor(kBlue);
            histcompare_ratio.SetMarkerSize(1);
            histcompare_ratio.GetXaxis()->SetTitle(TString::Format("#frac{%s Norm}{%s KPi3Pi_SIM_1ps}", analyzing->analyzinglabel, analyzing->analyzinglabel));
            histcompare_ratio.GetXaxis()->SetTitleOffset(1.7);
            histcompare_ratio.GetXaxis()->SetTitleSize(0.024);  
            histcompare_ratio.GetYaxis()->SetTitle("Occurences");
            histcompare_ratio.SetAxisRange(min, max, "X");
            histcompare_ratio.SetAxisRange(0, histcompare_ratio.GetBinContent(histcompare_ratio.GetMaximumBin()), "Y");
            histcompare_ratio.Draw("HF P");
            canvas_comparison_ratio.Update();

            //Adjusting the stats
            ratiostats = (TPaveStats *)histcompare_ratio.GetListOfFunctions()->FindObject("stats");
            ratiostats->SetName("Pisa Q^{2} distribution Ratio: KPiEtaPiPi_SIM_1ps over KPi3Pi_SIM_1ps");
            ratiostats->SetOptStat(112211);
            ratiostats->SetX1NDC(0.1);
            ratiostats->SetX2NDC(0.31);
            ratiostats->SetY1NDC(0.78);
            ratiostats->SetY2NDC(0.9);
            ratiostats->Draw("same");
            canvas_comparison_ratio.Print(TString::Format("%s/%s_ratio_norm_over_KPi3Pi_SIM_1ps.png", analyzing->output_directory_png, analyzing->variable_name));

            // Q2 Sig+Norm over KPi3Pi_SIM_1ps
            canvas_comparison_ratio.cd();
            canvas_comparison_ratio.Clear();
            histcompare_ratio = (TH1F)(*histcompare_signor.GetPtr() / *histcompare_bkg.GetPtr());

            histcompare_ratio.SetTitle(TString::Format("analyzing %s: Signal+KPiEtaPiPi_SIM_1ps over KPi3Pi_SIM_1ps", analyzing->analyzinglabel));
            histcompare_ratio.SetStats(kTRUE);
            histcompare_ratio.SetLineColorAlpha(kBlue, 1);
            histcompare_ratio.SetFillStyle(1001);
            histcompare_ratio.SetFillColor(kBlue);
            histcompare_ratio.SetFillColorAlpha(kBlue, 0.8);
            histcompare_ratio.SetMarkerStyle(kFullSquare);
            histcompare_ratio.SetMarkerColor(kBlue);
            histcompare_ratio.SetMarkerSize(1);
            histcompare_ratio.GetXaxis()->SetTitle(TString::Format("#frac{%s Sig+Norm}{%s KPi3Pi_SIM_1ps}", analyzing->analyzinglabel, analyzing->analyzinglabel));
            histcompare_ratio.GetXaxis()->SetTitleOffset(1.7);
            histcompare_ratio.GetXaxis()->SetTitleSize(0.024);
            histcompare_ratio.GetYaxis()->SetTitle("Occurences");
            histcompare_ratio.SetAxisRange(min, max, "X");
            histcompare_ratio.SetAxisRange(0, histcompare_ratio.GetBinContent(histcompare_ratio.GetMaximumBin()), "Y");
            histcompare_ratio.Draw("HF P");
            canvas_comparison_ratio.Update();

            //Adjusting the stats
            ratiostats = (TPaveStats *)histcompare_ratio.GetListOfFunctions()->FindObject("stats");
            ratiostats->SetName("Pisa Q^{2} distribution Ratio: Sig+Norm over KPi3Pi_SIM_1ps");
            ratiostats->SetOptStat(112211);
            ratiostats->SetX1NDC(0.1);
            ratiostats->SetX2NDC(0.31);
            ratiostats->SetY1NDC(0.78);
            ratiostats->SetY2NDC(0.9);
            ratiostats->Draw("same");
            canvas_comparison_ratio.Print(TString::Format("%s/%s_ratio_KPi3Pi_SIM_100ps_over_KPi3Pi_SIM_1ps.png", analyzing->output_directory_png, analyzing->variable_name));
        }
        
        if (analyzing->flag_rocplot == true)
        {   
            // Build string to print correct comparison title in ROC plots
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable) roc_curve_title.Append(TString::Format("%s", analyzing->analyzinglabel));
            else roc_curve_title.Append(TString::Format("%s,", analyzing->analyzinglabel));    
            
            // ROC curves: Signal with respect to KPiEtaPiPi_SIM_1ps
            canvas_comparison_roc_sigvsnorm.cd();
            TGraph roc_curve_sig_norm = TGraph(purity_KPi3Pi_SIM_100ps.size(), purity_KPi3Pi_SIM_100ps.data(), sig_eff.data());
            roc_curve_sig_norm.SetName(analyzing->variable_name);
            roc_curve_sig_norm.SetMarkerStyle(21);
            
            roc_sig_norm.Add((TGraph *)roc_curve_sig_norm.Clone());
            roc_sig_norm.Draw();

            //Adjusting the legend
            legend_roc_sigvsnorm.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_sigvsnorm.AddEntry(roc_curve_sig_norm.GetName(), analyzing->variable_name, "P");
            legend_roc_sigvsnorm.Draw();
            gPad->Update();
            legend_roc_sigvsnorm.SetTextSize(0.019);
            legend_roc_sigvsnorm.SetX1NDC(0.75);
            legend_roc_sigvsnorm.SetX2NDC(0.9);
            legend_roc_sigvsnorm.SetY1NDC(0.78);
            legend_roc_sigvsnorm.SetY2NDC(0.9);
            gPad->Modified();
            legend_roc_sigvsnorm_graph.push_back(roc_curve_sig_norm);
            legend_roc_sigvsnorm_entries.push_back(analyzing->variable_prettyname);
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_sig_norm.GetHistogram()->SetTitle(TString::Format("ROC (Sig vs Norm): %s comparison", roc_curve_title.Data()));
                roc_sig_norm.GetXaxis()->SetTitle("Purity := #frac{N(Sig)}{N(Sig)+N(Norm)}");
                roc_sig_norm.GetXaxis()->SetTitleOffset(1.5);
                roc_sig_norm.GetXaxis()->SetTitleSize(0.028);
                roc_sig_norm.GetXaxis()->SetLimits(0, 1);
                roc_sig_norm.GetYaxis()->SetTitle("Signal efficiency");
                roc_sig_norm.Draw("AP PLC PMC");
                legend_roc_sigvsnorm.Draw("same");
                canvas_comparison_roc_sigvsnorm.Print(TString::Format("%s/%s_roc_sig_vs_norm.png", analyzing->output_directory_png, analyzing->analyzingvariable));

                canvas_compraison_quadratiroc_sigvsnorm.cd();
                canvas_compraison_quadratiroc_sigvsnorm.Clear();
                gPad->SetRightMargin(0.05);
                gPad->SetLeftMargin(0.15);
                roc_sig_norm.Draw("AP PLC PMC");
                TLegend legend_quadrati = analyzing->SetLegendPosAuto("TR", 3, analyzing->legxlength);
                legend_quadrati.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
                int i = 0;
                for(auto graph=legend_roc_sigvsnorm_graph.begin(); graph<legend_roc_sigvsnorm_graph.end(); graph++)
                {
                    legend_quadrati.AddEntry(graph->GetName(), legend_roc_sigvsnorm_entries[i].Data(), "P");
                    i++;
                }
                legend_quadrati.SetTextSize(0.025);
                legend_quadrati.Draw("SAME");
                canvas_compraison_quadratiroc_sigvsnorm.Print(TString::Format("%s/%s_roc_sig_vs_norm.png", TString::Format("%s_quadrati",analyzing->output_directory_png).Data(), analyzing->analyzingvariable));

                roc_sig_norm.GetListOfGraphs()->Delete();
                legend_roc_sigvsnorm.Clear();    
                canvas_comparison_roc_sigvsnorm.Clear();
                legend_roc_sigvsnorm_graph.clear();
                legend_roc_sigvsnorm_entries.clear();
            }

            // ROC curves: Signal with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_sigvsbkg.cd();
            TGraph roc_curve_sig_bkg = TGraph(purity_sigbkg.size(), purity_sigbkg.data(), sig_eff.data());
            roc_curve_sig_bkg.SetName(analyzing->variable_name);
            roc_curve_sig_bkg.SetMarkerStyle(21);
            
            roc_sig_bkg.Add((TGraph *)roc_curve_sig_bkg.Clone());
            roc_sig_bkg.Draw();

            //Adjusting the legend
            legend_roc_sigvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_sigvsbkg.AddEntry(roc_curve_sig_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_sigvsbkg.Draw();
            gPad->Update();
            legend_roc_sigvsbkg.SetTextSize(0.019);
            legend_roc_sigvsbkg.SetX1NDC(0.75);
            legend_roc_sigvsbkg.SetX2NDC(0.9);
            legend_roc_sigvsbkg.SetY1NDC(0.78);
            legend_roc_sigvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            legend_roc_sigvsbkg_graph.push_back(roc_curve_sig_norm);
            legend_roc_sigvsbkg_entries.push_back(analyzing->variable_name);
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_sig_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Sig vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_sig_bkg.GetXaxis()->SetTitle("Purity := #frac{N(Sig)}{N(Sig)+N(Bkg)}");
                roc_sig_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_sig_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_sig_bkg.GetXaxis()->SetLimits(0, 1);
                roc_sig_bkg.GetYaxis()->SetTitle("Signal efficiency");
                roc_sig_bkg.Draw("AP PLC PMC");
                legend_roc_sigvsbkg.Draw("same");               
                canvas_comparison_roc_sigvsbkg.Print(TString::Format("%s/%s_roc_sig_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));

                canvas_compraison_quadratiroc_sigvsbkg.cd();
                canvas_compraison_quadratiroc_sigvsbkg.Clear();
                gPad->SetRightMargin(0.05);
                gPad->SetLeftMargin(0.15);
                roc_sig_bkg.Draw("AP PLC PMC");
                TLegend legend_quadrati = analyzing->SetLegendPosAuto("TR", 3, (analyzing->legxlength)-0.05);
                legend_quadrati.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
                int i = 0;
                for(auto graph=legend_roc_sigvsbkg_graph.begin(); graph<legend_roc_sigvsbkg_graph.end(); graph++)
                {
                    legend_quadrati.AddEntry(graph->GetName(), legend_roc_sigvsbkg_entries[i].Data(), "P");
                    i++;
                }
                legend_quadrati.SetTextSize(0.025);
                legend_quadrati.Draw("SAME");
                canvas_compraison_quadratiroc_sigvsbkg.Print(TString::Format("%s/%s_roc_sig_vs_bkg.png", TString::Format("%s_quadrati",analyzing->output_directory_png).Data(), analyzing->analyzingvariable));
                legend_roc_sigvsbkg.Clear();    
                legend_roc_sigvsbkg_graph.clear();
                legend_roc_sigvsbkg_entries.clear();
                
                roc_sig_bkg.GetListOfGraphs()->Delete();
                legend_roc_sigvsbkg.Clear();
                canvas_comparison_roc_sigvsbkg.Clear();
            } 

            // ROC curves: KPiEtaPiPi_SIM_1ps with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_normvsbkg.cd();
            TGraph roc_curve_norm_bkg = TGraph(purity_normbkg.size(), purity_normbkg.data(), norm_eff.data());
            roc_curve_norm_bkg.SetName(analyzing->variable_name);
            roc_curve_norm_bkg.SetMarkerStyle(21);
            
            roc_norm_bkg.Add((TGraph *)roc_curve_norm_bkg.Clone());
            roc_norm_bkg.Draw();

            //Adjusting the legend
            legend_roc_normvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_normvsbkg.AddEntry(roc_curve_norm_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_normvsbkg.Draw();
            gPad->Update();
            legend_roc_normvsbkg.SetTextSize(0.019);
            legend_roc_normvsbkg.SetX1NDC(0.75);
            legend_roc_normvsbkg.SetX2NDC(0.9);
            legend_roc_normvsbkg.SetY1NDC(0.78);
            legend_roc_normvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            legend_roc_normvsbkg_graph.push_back(roc_curve_norm_bkg);
            legend_roc_normvsbkg_entries.push_back(analyzing->variable_name);            
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_norm_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Norm vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_norm_bkg.GetXaxis()->SetTitle("Purity := #frac{N(Norm)}{N(Norm)+N(Bkg)}");
                roc_norm_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_norm_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_norm_bkg.GetXaxis()->SetLimits(0, 1);
                roc_norm_bkg.GetYaxis()->SetTitle("KPiEtaPiPi_SIM_1ps efficiency");
                roc_norm_bkg.Draw("AP PLC PMC");
                legend_roc_normvsbkg.Draw("same");
                canvas_comparison_roc_normvsbkg.Print(TString::Format("%s/%s_roc_norm_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));

                canvas_compraison_quadratiroc_normvsbkg.cd();
                canvas_compraison_quadratiroc_normvsbkg.Clear();
                gPad->SetRightMargin(0.05);
                gPad->SetLeftMargin(0.15);
                roc_norm_bkg.Draw("AP PLC PMC");
                TLegend legend_quadrati = analyzing->SetLegendPosAuto("TR", 3, (analyzing->legxlength)-0.05);
                legend_quadrati.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
                int i = 0;
                for(auto graph=legend_roc_normvsbkg_graph.begin(); graph<legend_roc_normvsbkg_graph.end(); graph++)
                {
                    legend_quadrati.AddEntry(graph->GetName(), legend_roc_normvsbkg_entries[i].Data(), "P");
                    i++;
                }
                legend_quadrati.SetTextSize(0.025);
                legend_quadrati.Draw("SAME");
                canvas_compraison_quadratiroc_normvsbkg.Print(TString::Format("%s/%s_roc_norm_vs_bkg.png", TString::Format("%s_quadrati",analyzing->output_directory_png).Data(), analyzing->analyzingvariable));
                legend_roc_normvsbkg.Clear();    
                legend_roc_normvsbkg_graph.clear();
                legend_roc_normvsbkg_entries.clear();

                roc_norm_bkg.GetListOfGraphs()->Delete();
                legend_roc_normvsbkg.Clear();
                canvas_comparison_roc_normvsbkg.Clear();
            }       

            // ROC curves: Signal-KPiEtaPiPi_SIM_1ps with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg.cd();
            TGraph roc_curve_KPi3Pi_SIM_100ps_bkg = TGraph(purity_KPi3Pi_SIM_100psbkg.size(), purity_KPi3Pi_SIM_100psbkg.data(), signor_eff.data());
            roc_curve_KPi3Pi_SIM_100ps_bkg.SetName(analyzing->variable_name);
            roc_curve_KPi3Pi_SIM_100ps_bkg.SetMarkerStyle(21);
            
            roc_KPi3Pi_SIM_100ps_bkg.Add((TGraph *)roc_curve_KPi3Pi_SIM_100ps_bkg.Clone());
            roc_KPi3Pi_SIM_100ps_bkg.Draw();


            //Adjusting the legend
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_KPi3Pi_SIM_100psvsbkg.AddEntry(roc_curve_KPi3Pi_SIM_100ps_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_KPi3Pi_SIM_100psvsbkg.Draw();
            gPad->Update();
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetTextSize(0.019);
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetX1NDC(0.75);
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetX2NDC(0.9);
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetY1NDC(0.78);
            legend_roc_KPi3Pi_SIM_100psvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            legend_roc_KPi3Pi_SIM_100psvsbkg_graph.push_back(roc_curve_KPi3Pi_SIM_100ps_bkg);
            legend_roc_KPi3Pi_SIM_100psvsbkg_entries.push_back(analyzing->variable_name);                
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {                
                roc_KPi3Pi_SIM_100ps_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Sig+Norm vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitle("Purity := #frac{N(Sign+Norm)}{N(Sign+Norm)+N(Bkg)}");
                roc_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetLimits(0, 1);
                roc_KPi3Pi_SIM_100ps_bkg.GetYaxis()->SetTitle("Signal+KPiEtaPiPi_SIM_1ps efficiency");
                roc_KPi3Pi_SIM_100ps_bkg.Draw("AP PLC PMC");
                legend_roc_KPi3Pi_SIM_100psvsbkg.Draw("same");
                canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg.Print(TString::Format("%s/%s_roc_KPi3Pi_SIM_100ps_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));

                canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg.cd();
                canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg.Clear();
                gPad->SetRightMargin(0.05);
                gPad->SetLeftMargin(0.15);
                roc_KPi3Pi_SIM_100ps_bkg.Draw("AP PLC PMC");
                TLegend legend_quadrati = analyzing->SetLegendPosAuto("TR", 3, (analyzing->legxlength)-0.05);
                legend_quadrati.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
                int i = 0;
                for(auto graph=legend_roc_KPi3Pi_SIM_100psvsbkg_graph.begin(); graph<legend_roc_KPi3Pi_SIM_100psvsbkg_graph.end(); graph++)
                {
                    legend_quadrati.AddEntry(graph->GetName(), legend_roc_KPi3Pi_SIM_100psvsbkg_entries[i].Data(), "P");
                    i++;
                }
                legend_quadrati.SetTextSize(0.025);
                legend_quadrati.Draw("SAME");
                canvas_compraison_quadratiroc_KPi3Pi_SIM_100psvsbkg.Print(TString::Format("%s/%s_roc_KPi3Pi_SIM_100ps_vs_bkg.png", TString::Format("%s_quadrati",analyzing->output_directory_png).Data(), analyzing->analyzingvariable));
                legend_roc_KPi3Pi_SIM_100psvsbkg.Clear();    
                legend_roc_KPi3Pi_SIM_100psvsbkg_graph.clear();
                legend_roc_KPi3Pi_SIM_100psvsbkg_entries.clear();

                roc_KPi3Pi_SIM_100ps_bkg.GetListOfGraphs()->Delete();
                legend_roc_KPi3Pi_SIM_100psvsbkg.Clear();
                canvas_comparison_roc_KPi3Pi_SIM_100psvsbkg.Clear();
            }

            // ROC curves inverted: Signal with respect to KPiEtaPiPi_SIM_1ps
            canvas_comparison_roc_inverted_sigvsnorm.cd();
            TGraph roc_inverted_curve_sig_norm = TGraph(purity_inverted_KPi3Pi_SIM_100ps.size(), purity_inverted_KPi3Pi_SIM_100ps.data(), sig_eff_inverted.data());
            roc_inverted_curve_sig_norm.SetName(analyzing->variable_name);
            roc_inverted_curve_sig_norm.SetMarkerStyle(21);
            
            roc_inverted_sig_norm.Add((TGraph *)roc_inverted_curve_sig_norm.Clone());
            roc_inverted_sig_norm.Draw();

            //Adjusting the legend
            legend_roc_inverted_sigvsnorm.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_inverted_sigvsnorm.AddEntry(roc_inverted_curve_sig_norm.GetName(), analyzing->variable_name, "P");
            legend_roc_inverted_sigvsnorm.Draw();
            gPad->Update();
            legend_roc_inverted_sigvsnorm.SetTextSize(0.019);
            legend_roc_inverted_sigvsnorm.SetX1NDC(0.75);
            legend_roc_inverted_sigvsnorm.SetX2NDC(0.9);
            legend_roc_inverted_sigvsnorm.SetY1NDC(0.78);
            legend_roc_inverted_sigvsnorm.SetY2NDC(0.9);
            gPad->Modified();
            
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_inverted_sig_norm.GetHistogram()->SetTitle(TString::Format("ROC (Sig vs Norm): %s comparison", roc_curve_title.Data()));
                roc_inverted_sig_norm.GetXaxis()->SetTitle("Inverted Purity := #frac{N(Sig)}{N(Sig)+N(Norm)}");
                roc_inverted_sig_norm.GetXaxis()->SetTitleOffset(1.5);
                roc_inverted_sig_norm.GetXaxis()->SetTitleSize(0.028);
                roc_inverted_sig_norm.GetXaxis()->SetLimits(0, 1);
                roc_inverted_sig_norm.GetYaxis()->SetTitle("Signal Inverted efficiency");
                roc_inverted_sig_norm.Draw("AP PLC PMC");
                legend_roc_inverted_sigvsnorm.Draw("same");
                canvas_comparison_roc_inverted_sigvsnorm.Print(TString::Format("%s/%s_roc_inverted_sig_vs_norm.png", analyzing->output_directory_png, analyzing->analyzingvariable));
                roc_inverted_sig_norm.GetListOfGraphs()->Delete();
                legend_roc_inverted_sigvsnorm.Clear();    
                canvas_comparison_roc_inverted_sigvsnorm.Clear();
            }

            // ROC curves inverted: Signal with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_inverted_sigvsbkg.cd();
            TGraph roc_inverted_curve_sig_bkg = TGraph(purity_inverted_sigbkg.size(), purity_inverted_sigbkg.data(), sig_eff_inverted.data());
            roc_inverted_curve_sig_bkg.SetName(analyzing->variable_name);
            roc_inverted_curve_sig_bkg.SetMarkerStyle(21);
            
            roc_inverted_sig_bkg.Add((TGraph *)roc_inverted_curve_sig_bkg.Clone());
            roc_inverted_sig_bkg.Draw();

            //Adjusting the legend
            legend_roc_inverted_sigvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_inverted_sigvsbkg.AddEntry(roc_inverted_curve_sig_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_inverted_sigvsbkg.Draw();
            gPad->Update();
            legend_roc_inverted_sigvsbkg.SetTextSize(0.019);
            legend_roc_inverted_sigvsbkg.SetX1NDC(0.75);
            legend_roc_inverted_sigvsbkg.SetX2NDC(0.9);
            legend_roc_inverted_sigvsbkg.SetY1NDC(0.78);
            legend_roc_inverted_sigvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_inverted_sig_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Sig vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_inverted_sig_bkg.GetXaxis()->SetTitle("Inverted Purity := #frac{N(Sig)}{N(Sig)+N(Bkg)}");
                roc_inverted_sig_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_inverted_sig_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_inverted_sig_bkg.GetXaxis()->SetLimits(0, 1);
                roc_inverted_sig_bkg.GetYaxis()->SetTitle("Signal Inverted efficiency");
                roc_inverted_sig_bkg.Draw("AP PLC PMC");
                legend_roc_inverted_sigvsbkg.Draw("same");               
                canvas_comparison_roc_inverted_sigvsbkg.Print(TString::Format("%s/%s_roc_inverted_sig_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));
                roc_inverted_sig_bkg.GetListOfGraphs()->Delete();
                legend_roc_inverted_sigvsbkg.Clear();
                canvas_comparison_roc_inverted_sigvsbkg.Clear();
            } 

            // ROC curves inverted: KPiEtaPiPi_SIM_1ps with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_inverted_normvsbkg.cd();
            TGraph roc_inverted_curve_norm_bkg = TGraph(purity_inverted_normbkg.size(), purity_inverted_normbkg.data(), norm_eff_inverted.data());
            roc_inverted_curve_norm_bkg.SetName(analyzing->variable_name);
            roc_inverted_curve_norm_bkg.SetMarkerStyle(21);
            
            roc_inverted_norm_bkg.Add((TGraph *)roc_inverted_curve_norm_bkg.Clone());
            roc_inverted_norm_bkg.Draw();

            //Adjusting the legend
            legend_roc_inverted_normvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_inverted_normvsbkg.AddEntry(roc_inverted_curve_norm_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_inverted_normvsbkg.Draw();
            gPad->Update();
            legend_roc_inverted_normvsbkg.SetTextSize(0.019);
            legend_roc_inverted_normvsbkg.SetX1NDC(0.75);
            legend_roc_inverted_normvsbkg.SetX2NDC(0.9);
            legend_roc_inverted_normvsbkg.SetY1NDC(0.78);
            legend_roc_inverted_normvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {
                roc_inverted_norm_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Norm vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_inverted_norm_bkg.GetXaxis()->SetTitle("Inverted Purity := #frac{N(Norm)}{N(Norm)+N(Bkg)}");
                roc_inverted_norm_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_inverted_norm_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_inverted_norm_bkg.GetXaxis()->SetLimits(0, 1);
                roc_inverted_norm_bkg.GetYaxis()->SetTitle("KPiEtaPiPi_SIM_1ps Inverted efficiency");
                roc_inverted_norm_bkg.Draw("AP PLC PMC");
                legend_roc_inverted_normvsbkg.Draw("same");
                canvas_comparison_roc_inverted_normvsbkg.Print(TString::Format("%s/%s_roc_inverted_norm_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));
                roc_inverted_norm_bkg.GetListOfGraphs()->Delete();
                legend_roc_inverted_normvsbkg.Clear();
                canvas_comparison_roc_inverted_normvsbkg.Clear();
            }       

            // ROC curves inverted: Signal-KPiEtaPiPi_SIM_1ps with respect to KPi3Pi_SIM_1ps
            canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg.cd();
            TGraph roc_inverted_curve_KPi3Pi_SIM_100ps_bkg = TGraph(purity_inverted_KPi3Pi_SIM_100psbkg.size(), purity_inverted_KPi3Pi_SIM_100psbkg.data(), signor_eff_inverted.data());
            roc_inverted_curve_KPi3Pi_SIM_100ps_bkg.SetName(analyzing->variable_name);
            roc_inverted_curve_KPi3Pi_SIM_100ps_bkg.SetMarkerStyle(21);
            
            roc_inverted_KPi3Pi_SIM_100ps_bkg.Add((TGraph *)roc_inverted_curve_KPi3Pi_SIM_100ps_bkg.Clone());
            roc_inverted_KPi3Pi_SIM_100ps_bkg.Draw();


            //Adjusting the legend
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetHeader(TString::Format("ROC curve for %s", analyzing->analyzingvariable), "C");
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.AddEntry(roc_inverted_curve_KPi3Pi_SIM_100ps_bkg.GetName(), analyzing->variable_name, "P");
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.Draw();
            gPad->Update();
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetTextSize(0.019);
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetX1NDC(0.75);
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetX2NDC(0.9);
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetY1NDC(0.78);
            legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.SetY2NDC(0.9);
            gPad->Modified();
            
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable)
            {                
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetHistogram()->SetTitle(TString::Format("ROC (Sig+Norm vs Bkg): %s comparison", roc_curve_title.Data()));
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitle("Inverted Purity := #frac{N(Sign+Norm)}{N(Sign+Norm)+N(Bkg)}");
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitleOffset(1.5);
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetTitleSize(0.028);
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetXaxis()->SetLimits(0, 1);
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetYaxis()->SetTitle("Signal+KPiEtaPiPi_SIM_1ps Inverted efficiency");
                roc_inverted_KPi3Pi_SIM_100ps_bkg.Draw("AP PLC PMC");
                legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.Draw("same");
                canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg.Print(TString::Format("%s/%s_roc_inverted_KPi3Pi_SIM_100ps_vs_bkg.png", analyzing->output_directory_png, analyzing->analyzingvariable));
                roc_inverted_KPi3Pi_SIM_100ps_bkg.GetListOfGraphs()->Delete();
                legend_roc_inverted_KPi3Pi_SIM_100psvsbkg.Clear();
                canvas_comparison_roc_inverted_KPi3Pi_SIM_100psvsbkg.Clear();
            }
            
            // Cleaning the memory of the roc title
            if (analyzing != tobediscriminated->end() && (analyzing+1)->analyzingvariable != analyzing->analyzingvariable) roc_curve_title.Clear();                    
        }*/
    }

    return 0;
}