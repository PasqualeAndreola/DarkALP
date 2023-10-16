/*!
 *  \file AnalysisVariable.h
 *  \brief Header file for \ref AnalysisVariable class.
 * 
 *  \class This class can be used to define a variable and its range which will be used throughout the analysis.
 *  
 */

#ifndef ANALYSISVARIABLE_H
#define ANALYSISVARIABLE_H

#include <TROOT.h>
#include <TH1.h>
#include <TLegend.h>
#include <TPad.h>
#include <iostream>
#include <chrono>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TStyle.h>

//It can be useful to use these namespaces
using namespace std;

class AnalysisVariable
{
public:
    char const * variable_name = "";            /*!< Expression by which the variable is stored in the input files */
    char const * variable_expression = "";      /*!< Expression by which the variable can be combined to form a new variable */
    char const * variable_prettyname = "";      /*!< Name that will be used as variable label (LaTex expression can be used) */
    char const * variable_dimension = "";       /*!< Unit of measurement of the variable */ 

    Double_t variable_cut_min = - 1e15;             /*!< Minimum of the variable: only values greater than this are considered*/
    Double_t variable_cut_max = 1e15;             /*!< Maximum of the variable: only values lower than this are considered*/
    Float_t variable_cut_pid = 0;              /*!< Particle identification number (please refer to PDG for further information) */
    bool variable_cut_pidflag = false;           /*!< If true, the particle has to be the one with this particle identification number. If false, the particle has to be different from the one with this particle identification*/
    bool variable_cutchecker = false;           /*!< If true, the event has to pass through the \ref variable_cut (if so, it is not a minimal cut)*/
    bool variable_mvacutfeature = false;           /*!< If true, the variable is used as a feature in the training of the MVA algorithm that optimizes the cut*/
    string variable_cut = "";                   /*!< String describing a linear cut between minimum and maximum or particle identification(has to be set with the constructor \ref SetVarCutOpt)*/
    string variable_stringcut = "";             /*!< String describing an additional cut (or a particularly complex one, if needed)*/
    Float_t variable_bins = 100;               /*!< Number of bins of the histogram that will represent the variable */
    Float_t variable_histmin = 0;              /*!< Minimum of bins of the histogram that will represent the variable */
    Float_t variable_histmax = 1;              /*!< Maximum of bins of the histogram that will represent the variable */
    bool variable_binauto = true;               /*!< Compute the minimum and the maximum of the histogram range that will be filled */

    int variable_padtextfont = 10;         /*!< Font of the text displayed in the various pads*/
    Float_t variable_padtextsize = 0.03;   /*!< Size of the text displayed in the various pads*/

    int     variable_padposition = 1;            /*!< Integer used to set the position of the legend in the canvas (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)*/
    Float_t variable_numberofentries = 1;        /*!< Number of entries in the pad (Default is 1) */
    Float_t variable_padxlength = 0.2;           /*!< Length of the  pad (Default is 0.2) */
    Float_t variable_padentryyheight = 0.04;     /*!< Height of each entry in the pad (Default is 0.04) */
    vector<Float_t> variable_padvertices = {0, 0, 0, 0}; /*!< Vector holding the vertices of the pad (x0, x1, y0, y1)*/
    bool variable_legendauto = true;            /*!< The legend is set up automatically in each function*/
    
    int     variable_headerposition = 1;          /*!< Integer used to set the position of the text header box in the canvas (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)*/
    Float_t variable_headernumberofentries = 1;   /*!< Number of entries in the header box (Default is 1) */
    Float_t variable_textpadxlength = 0.2;           /*!< Length of the text  pad (Default is 0.2) */
    Float_t variable_textpadentryyheight = 0.04;     /*!< Height of each entry in the text pad (Default is 0.04) */
    vector<Float_t> variable_textpadvertices = {0, 0, 0, 0}; /*!< Vector holding the vertices of the text pad (x0, x1, y0, y1)*/

    int     variable_statposition = 1;         /*!< Position of the stat pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)*/
    int     variable_statoptions = 1;          /*!< Stats that have to be printed in the stat pad (please refer to TPaveStats class)*/
    Float_t variable_statxlength = 0.2;        /*!< Width of the stat pad*/
    Float_t variable_statentryyheight = 0.04;  /*!< Height of each entry in the stat pad*/
    vector<Float_t> variable_statpadvertices = {0, 0, 0, 0}; /*!< Vector holding the vertices of the stat pad (x0, x1, y0, y1)*/

    char const * variable_histplotfolder = "OutputFiles/PNGPlots/Variables"; /*!< Folder which will store the plot of the histogram of the variable. Default is: "OutputFiles/PNGPlots/Variables/" */
    bool variable_plot_flag = false;            /*!< If true, \ref VariablePlotter will be used to plot an histogram of the variable in \ref variable_histplotfolder*/
    bool variable_plot_normalized_flag = false; /*!< If true, \ref VariablePlotter will be used to plot a normalized histogram of the variable in \ref variable_histplotfolder*/
    bool variable_logscale_flag = false;        /*!< If true, the histogram of the variable will be plot in a log scale */

    string filename = ""; /*!< Name of the file that holds the variable to be analyzed*/
    string treename = ""; /*!< Name of the tree that holds the variable to be analyzed*/
    string sampledata = "2015"; /*!< Data when the sample that holds the variable to be analyzed has been collected*/

    //! LegendPositionConstructor let the user construct the fields of the class used to plot the legend
    /*!
        \param legposition         Position of the legend pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param numberofentries     Number of entries that the legend pad will hold
        \param legxlength          Width of the legend pad
        \param legentryyheight     Height of each entry in the legend pad
    */
    int LegendPositionConstructor(int legposition, int numberofentries = 1, Float_t legxlength = 0.2, Float_t legentryyheight = 0.04)
    {
        variable_padposition = legposition;
        variable_numberofentries = numberofentries;
        variable_padxlength = legxlength;
        variable_padentryyheight = legentryyheight;
        variable_legendauto = false;
        
        return 0;
    }

    //! StatPositionConstructor let the user construct the fields of the class used to plot the stat pad
    /*!
        \param statposition        Position of the stat pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param statoptions         Stats that have to be printed in the stat pad (please refer to TPaveStats class)
        \param statxlength         Width of the stat pad
        \param statentryyheight    Height of each entry in the stat pad
    */
    int StatPositionConstructor(int statposition, int statoptions = 1, Float_t statxlength = 0.2, Float_t statentryyheight = 0.04)
    {
        variable_statposition = statposition;
        variable_statoptions = statoptions;
        variable_statxlength = statxlength;
        variable_statentryyheight = statentryyheight;
        
        return 0;
    }
    
    //! BoxHeaderPositionConstructor let the user construct the fields of the class used to plot the text header box
    /*!
        \param legposition         Position of the legend pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param numberofentries     Number of entries that the legend pad will hold
        \param legxlength          Width of the legend pad
        \param legentryyheight     Height of each entry in the legend pad
    */
    int BoxHeaderPositionConstructor(int headerboxposition, int headerboxnumberofentries = 1, Float_t headerboxxlength = 0.2, Float_t headerboxentryyheight = 0.04)
    {
        variable_headerposition = headerboxposition;
        variable_headernumberofentries = headerboxnumberofentries;
        variable_textpadxlength = headerboxxlength;
        variable_textpadentryyheight = headerboxentryyheight;
        
        return 0;
    }

    //! SetLegendPos returns the legend built from the four different vertices given by the user
    /*!
        \param x0 Left margin position
        \param x1 Right margin position
        \param y0 Bottom margin position
        \param y1 Top margin position
    */
    TLegend* SetLegendPos(Float_t x0, Float_t x1, Float_t y0, Float_t y1) 
    {
        variable_padvertices = {x0, x1, y0, y1};
        return new TLegend(x0, y0, x1, y1);
    };   

    //! SetLegendPosAuto returns the legend built from the different characteristics chosen by the user
    /*!
        \param variable_legendposition      Position of the legend pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param variable_numberofentries     Number of entries that the legend pad will hold
        \param variable_legxlength          Width of the legend pad
        \param variable_legentryyheight     Height of each entry in the legend pad
    */
    TLegend* SetLegendPosAuto(int variable_legendposition, int variable_numberofentries, Float_t variable_legxlength = 0.2, Float_t variable_legentryyheight = 0.04)
    {
        Float_t x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        switch (variable_legendposition)
        {
            case 1: 
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_legxlength;
                y1 = 0.9; y0 = y1 - variable_legentryyheight*variable_numberofentries;
                break;
            }
            case 2:
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_legxlength;
                y1 = 0.9; y0 = y1 - variable_legentryyheight*variable_numberofentries;
                break;
            }     
            case 3: 
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_legxlength;
                y0 = 0.1; y1 = y0 + variable_legentryyheight*variable_numberofentries;
                break;
            }
            case 4:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_legxlength;
                y0 = 0.1; y1 = y0 + variable_legentryyheight*variable_numberofentries;
                break;
            }
            case 5:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_legxlength;
                y1 = 0.5; y0 = y1 - variable_legentryyheight; 
                break;
            }
            default:
                break;        
        }

        variable_padvertices = {x0, x1, y0, y1};
        return new TLegend(x0, y0, x1, y1);
    };

    //! SetLegendPosAuto returns a Legend Pad with the different fields stored in the instance of the class
    TLegend* SetLegendPosAuto()
    {
        return SetLegendPosAuto(variable_padposition, variable_numberofentries, variable_padxlength, variable_padentryyheight);
    };

    //! SetStatPos returns the stat pad built from the four different vertices given by the user
    /*!
        \param x0 Left margin position
        \param x1 Right margin position
        \param y0 Bottom margin position
        \param y1 Top margin position
    */
    TPaveStats* SetStatPos(Float_t x0, Float_t x1, Float_t y0, Float_t y1) 
    {
        variable_statpadvertices = {x0, x1, y0, y1};
        return new TPaveStats(x0, y0, x1, y1);
    };   

    //! SetStatAuto returns the stat pad built from the different characteristics chosen by the user
    /*!
        \param variable_statposition         Position of the stat pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right; -1: below the legend)
        \param variable_statxlength          Width of the stat pad
        \param variable_statentryyheight     Height of each entry in the stat pad
    */
    TPaveStats* SetStatAuto(TH1* varfithist, TLegend *legend, int variable_statposition, Float_t variable_statxlength, Float_t variable_statentryyheight)
    {
        Float_t x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        TPaveStats *pavestat = (TPaveStats *)(varfithist->GetListOfFunctions()->FindObject("stats"));
        int statpad_numberofentries = pavestat->GetListOfLines()->GetSize();
        switch (variable_statposition)
        {
            case -1: 
            {
                if ((variable_padvertices[2] < 0.1 + variable_statentryyheight))
                {
                    variable_padvertices = {variable_padvertices[0], variable_padvertices[1],
                                            variable_padvertices[2]+variable_statentryyheight*statpad_numberofentries, variable_padvertices[3]+variable_statentryyheight*statpad_numberofentries};
                    legend->SetX1NDC(variable_padvertices[0]);
                    legend->SetX2NDC(variable_padvertices[1]); 
                    legend->SetY1NDC(variable_padvertices[2]); 
                    legend->SetY2NDC(variable_padvertices[3]);
                    legend->Draw("SAME");
                    gPad->Update();
                    gPad->Modified();
                }
                x0 = variable_padvertices[0];
                x1 = variable_padvertices[1];
                y0 = variable_padvertices[2]-variable_statentryyheight*statpad_numberofentries;
                y1 = variable_padvertices[2];
                break;
            }
            case 1: 
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_statxlength;
                y1 = 0.9; y0 = y1 - variable_statentryyheight*statpad_numberofentries;
                break;
            }
            case 2:
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_statxlength;
                y1 = 0.9; y0 = y1 - variable_statentryyheight*statpad_numberofentries;
                break;
            }     
            case 3: 
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_statxlength;
                y0 = 0.1; y1 = y0 + variable_statentryyheight*statpad_numberofentries;
                break;
            }
            case 4:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_statxlength;
                y0 = 0.1; y1 = y0 + variable_statentryyheight*statpad_numberofentries;
                break;
            }
            case 5:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_statxlength;
                y1 = 0.5; y0 = y1 - variable_statentryyheight; 
                break;
            }
            default:
                break;        
        }

        variable_statpadvertices = {x0, x1, y0, y1};
        pavestat->SetX1NDC(x0);
        pavestat->SetX2NDC(x1);
        pavestat->SetY1NDC(y0);
        pavestat->SetY2NDC(y1);

        return pavestat;
    };

    //! SetStatAuto returns a Stat Pad with the different fields stored in the instance of the class
    /*!
        \param varfithist Histogram which statistics have to be used
    */
    TPaveStats* SetStatAuto(TH1 *varfithist)
    {
        return SetStatAuto(varfithist, NULL, variable_statposition, variable_statxlength, variable_statentryyheight);
    };

    //! SetStatAuto returns a Stat Pad with the different fields stored in the instance of the class
    /*!
        \param varfithist Histogram which statistics have to be used
    */
    TPaveStats* SetStatAuto(TH1 *varfithist, TLegend *legend)
    {
        return SetStatAuto(varfithist, legend, variable_statposition, variable_statxlength, variable_statentryyheight);
    };

    //! SetTextBox returns the text pad built from the four different vertices given by the user
    /*!
        \param x0 Left margin position
        \param x1 Right margin position
        \param y0 Bottom margin position
        \param y1 Top margin position
    */
    TPaveText* SetTextBox(Float_t x0, Float_t x1, Float_t y0, Float_t y1) 
    {
        variable_textpadvertices = {x0, x1, y0, y1};
        return new TPaveText(x0, y0, x1, y1);
    };   
    
    //! SetTextBoxAuto returns the text pad built from the different characteristics chosen by the user
    /*!
        \param variable_textpadposition     Position of the legend pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param variable_numberofentries     Number of entries that the legend pad will hold
        \param variable_padxlength          Width of the legend pad
        \param variable_padentryyheight     Height of each entry in the legend pad
    */
    TPaveText* SetTextBoxAuto(int variable_textpadposition, int variable_numberofentries, Float_t variable_padxlength = 0.2, Float_t variable_padentryyheight = 0.04)
    {
        Float_t x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        switch (variable_textpadposition)
        {
            case 1: 
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y1 = 0.9; y0 = y1 - variable_padentryyheight*variable_numberofentries;
                break;
            }
            case 2:
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_padxlength;
                y1 = 0.9; y0 = y1 - variable_padentryyheight*variable_numberofentries;
                break;
            }     
            case 3: 
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_padxlength;
                y0 = 0.1; y1 = y0 + variable_padentryyheight*variable_numberofentries;
                break;
            }
            case 4:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y0 = 0.1; y1 = y0 + variable_padentryyheight*variable_numberofentries;
                break;
            }

            default:
                break;        
        }

        variable_textpadvertices = {x0, x1, y0, y1};
        return new TPaveText(x0, y0, x1, y1, "NDC");
    };

    //! SetTextBoxAuto returns the text pad built from the different characteristics chosen by the user
    /*!
        \param variable_textpadposition     Position of the text pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param x0                           Left margin of the text pad
        \param y1                           Top margin of the text pad
        \param variable_numberofentries     Number of entries that the text pad will hold
        \param variable_padxlength          Width of the text pad
        \param variable_padentryyheight     Height of each entry in the text pad
    */
    TPaveText* SetTextBoxAuto(int variable_textpadposition, Float_t x0, Float_t y1, int variable_numberofentries, Float_t variable_padxlength = 0.2, Float_t variable_padentryyheight = 0.04)
    {
        Float_t x1 = 0, y0 = 0;

        switch (variable_textpadposition)
        {
            case 1: 
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y1 = 0.9; y0 = y1 - variable_padentryyheight*variable_numberofentries;
                break;
            }
            case 2:
            {
                x1 = x0 + variable_padxlength;
                y0 = y1 - variable_padentryyheight*variable_numberofentries;
                break;
            }     
            case 3: 
            {
                x1 = x0 + variable_padxlength;
                y1 = y0 + variable_padentryyheight*variable_numberofentries;
                break;
            }
            case 4:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y0 = 0.1; y1 = y0 + variable_padentryyheight*variable_numberofentries;
                break;
            }

            default:
                break;        
        }

        variable_textpadvertices = {x0, x1, y0, y1};
        return new TPaveText(x0, y0, x1, y1, "NDC");
    };

    //! SetTextBoxAuto returns the text pad with the different fields stored in the instance of the class
    TPaveText* SetTextBoxAuto()
    {
        return SetTextBoxAuto(variable_padposition, int(variable_numberofentries), variable_textpadxlength, variable_textpadentryyheight);
    };

    //! SetTextBoxHeaderAuto returns a text box with one entry that lies on the top of another box with its entries
    /*!
        \param variable_boxheaderposition   Position of the text pad (1: top right -default; 2: top left; 3: bottom left; 4: bottom right)
        \param variable_numberofentries     Number of entries of the pad that will lie below the header
        \param variable_padxlength          Width of the text pad
        \param variable_padentryyheight     Height of each entry in the text pad below the header
    */
    TPaveText* SetTextBoxHeaderAuto(int variable_boxheaderposition, int variable_numberofentries, Float_t variable_padxlength = 0.2, Float_t variable_padentryyheight = 0.04)
    {
        Float_t x0 = 0, x1 = 0, y0 = 0, y1 = 0;
        switch (variable_boxheaderposition)
        {
            case -1: 
            {
                x0 = variable_statpadvertices[1]-variable_padxlength;
                x1 = variable_statpadvertices[1];
                y0 = variable_statpadvertices[2]-variable_padentryyheight;
                y1 = variable_statpadvertices[2];
                break;
            }
            case 2:
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_padxlength;
                y1 = 0.9; y0 = y1 - variable_padentryyheight;
                break;
            }
            case 3:
            {
                x0 = gPad->GetLeftMargin(); x1 = x0 + variable_padxlength;
                y0 = 0.1 + variable_padentryyheight*variable_numberofentries; y1 = y0 + variable_padentryyheight;
                break;
            }
            case 4:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y0 = 0.1 + variable_padentryyheight*variable_numberofentries; y1 = y0 + variable_padentryyheight;
                break;
            }
            case 5:
            {
                x1 = 1-gPad->GetRightMargin(); x0 = x1 - variable_padxlength;
                y1 = 0.5; y0 = y1 - variable_padentryyheight; 
                break;
            }

            default:
                break;
        }

        variable_textpadvertices = {x0, x1, y0, y1};
        return new TPaveText(x0, y0, x1, y1, "NDC");
    };

    //! SetTextBoxHeaderAuto returns a text box with one entry that lies on the top of another box with the different fields stored in the instance of the class
    TPaveText* SetTextBoxHeaderAuto()
    {
        return SetTextBoxHeaderAuto(variable_headerposition, int(variable_headernumberofentries), variable_textpadxlength, variable_textpadentryyheight);
    };

    //! SetVarBinningOpt is the default constructor used to set up the number of bins of the histogram that should represent the variable and the min/max of the histogram
    /*! 
    *   \param num_bins Number of the bins of the histogram that will represent the variable being analyzed
    *   \param hist_min Minimum of the histogram that will represent the variable being analyzed
    *   \param hist_max Maximum of the histogram that will represent the variable being analyzed
    */
    int SetVarBinningOpt(Float_t num_bins, Float_t hist_min, Float_t hist_max)
    {
        variable_bins = num_bins;
        variable_histmin = hist_min;
        variable_histmax = hist_max;
        variable_binauto = false;
        return 0;
    }

    //! SetVarCutOpt is the default constructor used to set up the cuts of the variable being analyzed. It returns a string that symbolizes the chosen cuts.
    /*! 
    *   \param cut_min Set the minimum threshold for the variable being analyzed
    *   \param cut_max Set the maximum threshold for the variable being analyzed
    */
    string SetVarCutOpt(Double_t cut_min = - 1e15, Double_t cut_max = 1e15)
    {
        if ( (cut_min != (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut_max = cut_max;
            variable_cut = TString::Format("(%s>%f) && (%s<%f)", variable_expression, variable_cut_min, variable_expression, variable_cut_max).Data();
        }
        else if ( (cut_min == (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_max = cut_max;
            variable_cut = TString::Format("(%s<%f)", variable_expression, variable_cut_max).Data();
        }
        else if ( (cut_min != (- 1e15 )) && (cut_max == (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut = TString::Format("(%s>%f)", variable_expression, variable_cut_min).Data();
        }
        return variable_cut;
    }

    //! SetVarCutOptMax is the default constructor used to set up the upper cuts of the variable being analyzed. It returns a string that symbolizes the chosen cuts.
    /*! 
    *   \param cut_max Set the maximum threshold for the variable being analyzed
    */
    string SetVarCutOptMax(Double_t cut_max = 1e15)
    {
        return SetVarCutOpt(-1e15, cut_max);
    }

    //! SetVarCutOptMin is the default constructor used to set up the bottom cuts of the variable being analyzed. It returns a string that symbolizes the chosen cuts.
    /*! 
    *   \param cut_max Set the minimum threshold for the variable being analyzed
    */
    string SetVarCutOptMin(Double_t cut_min = -1e15)
    {
        return SetVarCutOpt(cut_min, 1e15);
    }

    //! SetVarCutOptTruePId is the default constructor used to set up the cuts of the variable being analyzed. It returns a string that symbolizes the chosen cuts (particle has to be different from PId).
    /*! 
    *   \param cut_min Set the minimum threshold for the variable being analyzed
    *   \param cut_max Set the maximum threshold for the variable being analyzed
    *   \param cut_pid Set the particle identification code which the variable should be equal to
    */
    string SetVarCutOptTruePId(Double_t cut_min = - 1e15, Double_t cut_max = 1e15, Float_t cut_pid = 0)
    {
        if ( (cut_min != (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut_max = cut_max;
            variable_cut = TString::Format("(%s>%f) && (%s<%f)", variable_expression, variable_cut_min, variable_expression, variable_cut_max).Data();
        }
        else if ( (cut_min == (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_max = cut_max;
            variable_cut = TString::Format("(%s<%f)", variable_expression, variable_cut_max).Data();
        }
        else if ( (cut_min != (- 1e15 )) && (cut_max == (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut = TString::Format("(%s>%f)", variable_expression, variable_cut_min).Data();
        }
        if ( cut_pid != 0) 
        {
            variable_cut_pid = cut_pid; 
            variable_cut.append(TString::Format("&& (%s==%f)", variable_expression, variable_cut_pid).Data());
        }
        return variable_cut;
    }

    //! SetVarCutOptFalsePId is the default constructor used to set up the cuts of the variable being analyzed. It returns a string that symbolizes the chosen cuts (particle has to be different from PId).
    /*! 
    *   \param cut_min Set the minimum threshold for the variable being analyzed
    *   \param cut_max Set the maximum threshold for the variable being analyzed
    *   \param cut_pid Set the particle identification code which the variable should be equal to
    */
    string SetVarCutOptFalsePId(Double_t cut_min = - 1e15, Double_t cut_max = 1e15, Float_t cut_pid = 0)
    {
        if ( (cut_min != (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut_max = cut_max;
            variable_cut.append(TString::Format("(%s>%f) && (%s<%f)", variable_expression, variable_cut_min, variable_expression, variable_cut_max).Data());
        }
        else if ( (cut_min == (- 1e15 )) && (cut_max != (  1e15 )) ) 
        {
            variable_cut_max = cut_max;
            variable_cut.append(TString::Format("(%s<%f)", variable_expression, variable_cut_max).Data());
        }
        else if ( (cut_min != (- 1e15 )) && (cut_max == (  1e15 )) ) 
        {
            variable_cut_min = cut_min;
            variable_cut.append(TString::Format("(%s>%f)", variable_expression, variable_cut_min).Data());
        }
        if ( cut_pid != 0) 
        {
            variable_cut_pid = cut_pid; 
            variable_cut.append(TString::Format("&& (%s!=%f)", variable_expression, variable_cut_pid).Data());
        }
        return variable_cut;
    }

    //! SetVarNames is the default constructor used to set up the name (both plain and pretty in latex), the expression used to compute and the dimension of the variable
    /*! 
    *   \param var_name Set the name of the variable (\ref variable_name)
    *   \param var_prettyname Set the pretty name of the variable (\ref variable_prettyname)
    *   \param var_expression Set the expression used to compute the variable (\ref variable_expression)
    *   \param var_dimension Set the dimension of the variable going to be analyzed (\ref variable_dimension)
    */
    int SetVarNames(const char* var_name = "", const char* var_prettyname = "", const char* var_expression = "", const char* var_dimension = "")
    {
        variable_name = var_name;
        variable_prettyname = var_prettyname;
        variable_expression = var_expression;
        variable_dimension = var_dimension;
        return 0;
    }

    //! Xlabel returns the string formed by a variable name and a variable_dimension
    /*!
        \param name Name of the variable that will be plotted
        \param dimension Dimension of the variable that will be plotted
    */
    TString Xlabel(char const* name, char const* dimension) 
    {
        return TString::Format("%s#left[%s#right]", name, dimension);
    }; 

    //! Xlabel returns the string formed by \ref variable name and \ref variable_dimension
    /*!
        \param variable_prettyname Variable that will be plotted (\ref variable_prettyname)
        \param variable_dimension Dimension of the variable that will be plotted (\ref variable_dimension)
    */
    TString Xlabel() {return Xlabel(variable_prettyname, variable_dimension);};

    //! VariablePlotter is a function that plots the variables instantianed with this class
    static int VariableComparer(unordered_map<string, pair<string, string>> *vec_data_holder, AnalysisVariable *var_to_be_analyzed, string cut = "", chrono::_V2::system_clock::time_point start = chrono::system_clock::now(), bool debug = false);

    static int VariablesComparer(unordered_map<string, pair<string, string>> *vec_data_holder, vector<AnalysisVariable*> *vec_var_to_be_analyzed, string cut = "", chrono::_V2::system_clock::time_point start = chrono::system_clock::now(), bool debug = false)
    {
        for (vector<AnalysisVariable*>::iterator analvar = vec_var_to_be_analyzed->begin(); analvar != vec_var_to_be_analyzed->end(); analvar++)
        {
            AnalysisVariable::VariableComparer(vec_data_holder, *analvar);
        }
        
        return 0;
    };

    static int VariablePlotter(pair<string, string> *data_holder, AnalysisVariable *var_to_be_analyzed, string cut = "", chrono::_V2::system_clock::time_point start = chrono::system_clock::now(), bool debug = false);

    static int VariablePlotter2D(pair<string, string> *data_holder,AnalysisVariable *var_to_be_analyzed, AnalysisVariable *var_to_be_analyzed2, string cut = "", bool variable_logscale_zaxis_flag = false, chrono::_V2::system_clock::time_point start = chrono::system_clock::now(), bool debug = false);
};
#endif