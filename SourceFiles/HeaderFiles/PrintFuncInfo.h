/*!
 *  \file PrintFuncInfo.h
 *  \brief Header file for \ref PrintFuncInfo functions and \ref ElapsedTimeStamper function.
 */

#ifndef PRINTFUNCINFO_H
#define PRINTFUNCINFO_H

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <TString.h>
#include <chrono>
#include <ctime>

/*These namespaces can be useful*/
using namespace std;

void PrintFuncInfo(vector<TString> string2print);
void PrintFuncInfo(std::ostream& stream, vector<TString> string2print);
void ElapsedTimeStamper(chrono::_V2::system_clock::time_point start = chrono::system_clock::now());

#endif