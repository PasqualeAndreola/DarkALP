/*!
 *  \file Discriminants.h
 *  \brief Header file for some variables calculated with lambda function
 */

#ifndef ANALYSISVARIABLELAMBDACONSTRUCTOR_H
#define ANALYSISVARIABLELAMBDACONSTRUCTOR_H

#include <iostream>
#include <vector>
#include <TROOT.h>
#include <TLorentzVector.h>

/*These namespaces can be useful*/
using namespace std;

inline auto b_m_constrained_3body = [](Double_t jpsi_px, Double_t jpsi_py, Double_t jpsi_pz, Double_t jpsi_massPDG,
                                       Double_t kaon_px, Double_t kaon_py, Double_t kaon_pz, Double_t kaon_E,
                                       Double_t pion_px, Double_t pion_py, Double_t pion_pz, Double_t pion_E)
{
  TLorentzVector jpsiained;
  jpsiained.SetXYZM(jpsi_px, jpsi_py, jpsi_pz, jpsi_massPDG);
  TLorentzVector child_not_constrained_1(kaon_px, kaon_py, kaon_pz, kaon_E);
  TLorentzVector child_not_constrained_2(pion_px, pion_py, pion_pz, pion_E);

  return (jpsiained + child_not_constrained_1 + child_not_constrained_2).M();
};

inline auto b_m_jpsiconstrained_3body = [](Double_t jpsi_px, Double_t jpsi_py, Double_t jpsi_pz,
                                           Double_t kaon_px, Double_t kaon_py, Double_t kaon_pz, Double_t kaon_E,
                                           Double_t pion_px, Double_t pion_py, Double_t pion_pz, Double_t pion_E)
{
  Double_t JPSI_MASS = 3096.900;
  TLorentzVector jpsi_p4;
  jpsi_p4.SetXYZM(jpsi_px, jpsi_py, jpsi_pz, JPSI_MASS);
  TLorentzVector kaon_p4(kaon_px, kaon_py, kaon_pz, kaon_E);
  TLorentzVector pion_p4(pion_px, pion_py, pion_pz, pion_E);
  TLorentzVector b_fourmomentum_jpsiconstrained = jpsi_p4 + kaon_p4 + pion_p4;
  
  Double_t mass = b_fourmomentum_jpsiconstrained.M();

  return mass;
};

inline auto x_m_pionstrained_3body = [](Double_t pi0_px, Double_t pi0_py, Double_t pi0_pz,
                                           Double_t pionp_px, Double_t pionp_py, Double_t pionp_pz, Double_t pionp_E,
                                           Double_t pionm_px, Double_t pionm_py, Double_t pionm_pz, Double_t pionm_E)
{
  Double_t PI0_MASS = 134.9768;
  TLorentzVector pi0_p4;
  pi0_p4.SetXYZM(pi0_px, pi0_py, pi0_pz, PI0_MASS);
  TLorentzVector pionp_p4(pionp_px, pionp_py, pionp_pz, pionp_E);
  TLorentzVector pionm_p4(pionm_px, pionm_py, pionm_pz, pionm_E);
  TLorentzVector b_fourmomentum_jpsiconstrained = pi0_p4 + pionp_p4 + pionm_p4;
  
  Double_t mass = b_fourmomentum_jpsiconstrained.M();

  return mass;
};

#endif