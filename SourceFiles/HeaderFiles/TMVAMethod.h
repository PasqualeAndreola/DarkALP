/*!
 *  \file TMVAMethod.h
 *  \brief Header file for \ref TMVAMethod class.
 * 
 *  \class This class is needed to associate a flag to each file
 */

#ifndef TMVAMETHOD_H
#define TMVAMETHOD_H

#include <string>
#include <TROOT.h>
#include <TString.h>
#include <TMVA/Types.h>

using namespace std;

class TMVAMethod
{
public:
    const char *tmvamethodname = "";
    TString tmvamethodconfiguration = "";
    TMVA::Types::EMVA tmvamethodtype;
};

#endif