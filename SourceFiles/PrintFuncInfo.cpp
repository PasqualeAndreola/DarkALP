/*!
 *  \file PrintFuncInfo.cpp
 *  \brief Source file for \ref PrintFuncInfo function implementation
 */

/*!
 *  \fn void PrintFuncInfo(vector<TString> string2print) "";
 *  \brief Function used to print a banner with function info
 * 
 *  This function is used to print a banner with some "*" characters.\n
 *  These "*" are arranged in a rectangular format.\n
 *  They contains information given by the user as a vector of TString.\n
 *  Global variable \c ncolumns can be defined by the user: the function will appropriately cut the strings.\n
 *  \param string2print vector of TString that the user wants to print
 */

#include "HeaderFiles/PrintFuncInfo.h"

/*These namespaces can be useful*/
using namespace std;

/*Global definition of number of columns of the banner*/
const int ncolumns = 55; /*!< Number of columns of "*" that will be printed */

/* Function printing on standard output std::cout */
void PrintFuncInfo(vector<TString> string2print)
{
    int nrows = string2print.size() * 2 + 3;
    int stringindex = 0;
    int printedcharactercounter = 0;    //Counter of printed characters to ensure all cut strings are considering their previous offset
    vector<int> partiallength;          //Length of the small strings remaining from cuts: even are left respect to the half
    vector<int> minchar, maxchar;       //Number of spaces from left/right column
    cout << endl;

    for (int isymb = 0; isymb < nrows; isymb++)
    {
        int flag = 0;
        for (int jsymb = 0; jsymb < ncolumns; jsymb++)
        {
            if (isymb == 0 || isymb == nrows - 1)
            {
                cout << "*";
            }
            else if (isymb % 2 == 0 && flag == 0)
            {
                /*Variables used to control strings longer than limits*/
                int length = string2print[stringindex].Length();
                int iterator = 1;          //Number of times that a string has been cut
                int maximum = length;           //Maximum length of cut string

                /*Controls how many times we have to crop a string*/
                while ( maximum  > (ncolumns-2))
                {
                    partiallength.clear();
                    minchar.clear();
                    maxchar.clear();
                    length = floor( string2print[stringindex].Length() / (iterator+1) );
                    maximum = 0;
                    iterator = iterator + 1;
                    if (string2print[stringindex].First(" ") > 0)
                    {
                        for (int ncut = 0; ncut < iterator; ncut++)
                        {
                            if (string2print[stringindex].Index(" ", 1,length+accumulate(partiallength.begin(),partiallength.end(),0), TString::ECaseCompare::kExact) > 0)
                            {
                                length = string2print[stringindex].Index(" ", 1,length+accumulate(partiallength.begin(),partiallength.end(),0), TString::ECaseCompare::kExact) - accumulate( partiallength.begin(), partiallength.end(), 0);
                            }
                            else
                            {
                                length = string2print[stringindex].Length() - accumulate( partiallength.begin(), partiallength.end(), 0);
                            }
                            partiallength.push_back(length);
                            maximum = *max_element(partiallength.begin(), partiallength.end());
                            minchar.push_back(round((ncolumns - length)/ 2));
                            maxchar.push_back(round((ncolumns + length)/ 2));
                        }                  
                    }
                }

                /*Print all strings*/
                if (string2print[stringindex].Length() > ncolumns)
                {                    
                    /*Print cut strings*/
                    printedcharactercounter = 0;
                    for (int i = 0; i < iterator; i++)
                    {
                        int pcc = 0;
                        for (jsymb = 0; jsymb <ncolumns; jsymb++)
                        {
                            if ( jsymb == 0 || jsymb == ncolumns-1)
                            {
                                cout << "*";                               
                            }
                            else if (minchar[i]<=jsymb && jsymb<=maxchar[i] && pcc<partiallength[i] && printedcharactercounter<string2print[stringindex].Length())
                            {
                                cout << string2print[stringindex][printedcharactercounter];
                                pcc++;
                                printedcharactercounter++;
                            }
                            else
                            {
                                cout << " ";
                            }
                        }
                        if (i < iterator - 1) cout << endl;
                    }
                    printedcharactercounter = 0;
                    partiallength.clear();
                    minchar.clear();
                    maxchar.clear();
                    stringindex++;
                    flag = 1;
                }
                /*Print uncut strings*/
                else
                {
                    printedcharactercounter = 0;
                    minchar.push_back(round((ncolumns - length) / 2));
                    maxchar.push_back(round((ncolumns + length) / 2));
                    for (jsymb=0; jsymb < ncolumns; jsymb++)
                    {
                        if (jsymb == 0 || jsymb == ncolumns -1)
                        {
                            cout << "*";
                        }
                        else if (minchar[0] <= jsymb && jsymb <= maxchar[0] && printedcharactercounter<string2print[stringindex].Length())
                        {
                            cout << string2print[stringindex][printedcharactercounter];
                            printedcharactercounter++;
                        }
                        else
                        {
                            cout << " ";
                        }
                    }
                    stringindex++;
                    minchar.clear();
                    maxchar.clear();
                    flag = 1;
                }
            }
            else if (jsymb == 0 || jsymb == ncolumns-1)
            {
                cout << "*";
            }
            else
            {
                cout << " ";
            }            
        }
        cout << endl;
    }
    cout << endl;
    return;
}

/*!
 *  \overload PrintFuncInfo(std::ostream& stream, vector<TString> string2print) "";
 *  \brief Function used to print a banner with function info
 * 
 *  This function is the overloaded version of \ref PrintFuncInfo(vector<TString> string2print)\n
 *  The user can set an ofstream to write the output to whatever he/she likes.\n
 *  The developer used this to make unit test.\n
 * 
 *  \param stream       Output where the starred banner has to be written
 *  \param string2print vector of TString that the user wants to print
 */

/* Function printing on different output*/
void PrintFuncInfo(std::ostream& stream, vector<TString> string2print)
{
    int nrows = string2print.size() * 2 + 3;
    int stringindex = 0;
    int printedcharactercounter = 0;    //Counter of printed characters to ensure all cut strings are considering their previous offset
    vector<int> partiallength;          //Length of the small strings remaining from cuts: even are left respect to the half
    vector<int> minchar, maxchar;       //Number of spaces from left/right column
    stream << endl;

    for (int isymb = 0; isymb < nrows; isymb++)
    {
        int flag = 0;
        for (int jsymb = 0; jsymb < ncolumns; jsymb++)
        {
            if (isymb == 0 || isymb == nrows - 1)
            {
                stream << "*";
            }
            else if (isymb % 2 == 0 && flag == 0)
            {
                /*Variables used to control strings longer than limits*/
                int length = string2print[stringindex].Length();
                int iterator = 1;          //Number of times that a string has been cut
                int maximum = length;           //Maximum length of cut string

                /*Controls how many times we have to crop a string*/
                while ( maximum  > (ncolumns-2))
                {
                    partiallength.clear();
                    minchar.clear();
                    maxchar.clear();
                    length = floor( string2print[stringindex].Length() / (iterator+1) );
                    maximum = 0;
                    iterator = iterator + 1;
                    if (string2print[stringindex].First(" ") > 0)
                    {
                        for (int ncut = 0; ncut < iterator; ncut++)
                        {
                            if (string2print[stringindex].Index(" ", 1,length+accumulate(partiallength.begin(),partiallength.end(),0), TString::ECaseCompare::kExact) > 0)
                            {
                                length = string2print[stringindex].Index(" ", 1,length+accumulate(partiallength.begin(),partiallength.end(),0), TString::ECaseCompare::kExact) - accumulate( partiallength.begin(), partiallength.end(), 0);
                            }
                            else
                            {
                                length = string2print[stringindex].Length() - accumulate( partiallength.begin(), partiallength.end(), 0);
                            }
                            partiallength.push_back(length);
                            maximum = *max_element(partiallength.begin(), partiallength.end());
                            minchar.push_back(round((ncolumns - length)/ 2));
                            maxchar.push_back(round((ncolumns + length)/ 2));
                        }                  
                    }
                }

                /*Print all strings*/
                if (string2print[stringindex].Length() > ncolumns)
                {                    
                    /*Print cut strings*/
                    printedcharactercounter = 0;
                    for (int i = 0; i < iterator; i++)
                    {
                        int pcc = 0;
                        for (jsymb = 0; jsymb <ncolumns; jsymb++)
                        {
                            if ( jsymb == 0 || jsymb == ncolumns-1)
                            {
                                stream << "*";                               
                            }
                            else if (minchar[i]<=jsymb && jsymb<=maxchar[i] && pcc<partiallength[i] && printedcharactercounter<string2print[stringindex].Length())
                            {
                                stream << string2print[stringindex][printedcharactercounter];
                                pcc++;
                                printedcharactercounter++;
                            }
                            else
                            {
                                stream << " ";
                            }
                        }
                        if (i < iterator - 1) stream << endl;
                    }
                    printedcharactercounter = 0;
                    partiallength.clear();
                    minchar.clear();
                    maxchar.clear();
                    stringindex++;
                    flag = 1;
                }
                /*Print uncut strings*/
                else
                {
                    printedcharactercounter = 0;
                    minchar.push_back(round((ncolumns - length) / 2));
                    maxchar.push_back(round((ncolumns + length) / 2));
                    for (jsymb=0; jsymb < ncolumns; jsymb++)
                    {
                        if (jsymb == 0 || jsymb == ncolumns -1)
                        {
                            stream << "*";
                        }
                        else if (minchar[0] <= jsymb && jsymb <= maxchar[0] && printedcharactercounter<string2print[stringindex].Length())
                        {
                            stream << string2print[stringindex][printedcharactercounter];
                            printedcharactercounter++;
                        }
                        else
                        {
                            stream << " ";
                        }
                    }
                    stringindex++;
                    minchar.clear();
                    maxchar.clear();
                    flag = 1;
                }
            }
            else if (jsymb == 0 || jsymb == ncolumns-1)
            {
                stream << "*";
            }
            else
            {
                stream << " ";
            }            
        }
        stream << endl;
    }
    stream << endl;
    return;
}

void ElapsedTimeStamper(chrono::_V2::system_clock::time_point start)    
{
    chrono::_V2::system_clock::time_point end = chrono::system_clock::now();
    chrono::_V2::system_clock::duration elapsed_seconds = end-start;
    cout << "Time elapsed from the starting date: " << elapsed_seconds.count()*1e-9 << "s" << endl;
    return;
}

