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

inline auto b_m_jpsiconstrained_k_pi0misidentifiedask = [](Double_t jpsi_px, Double_t jpsi_py, Double_t jpsi_pz,
                                                           Double_t kaon_px, Double_t kaon_py, Double_t kaon_pz, Double_t kaon_E,
                                                           Double_t pion_px, Double_t pion_py, Double_t pion_pz, Double_t pion_E)
{
  Double_t JPSI_MASS = 3096.900;
  Double_t KAON_MASS = 493.677;
  Double_t PION_MASS = 139.570;

  TLorentzVector jpsi_p4, pionmisidentifiedask_p4;
  jpsi_p4.SetXYZM(jpsi_px, jpsi_py, jpsi_pz, JPSI_MASS);
  TLorentzVector kaon_p4(kaon_px, kaon_py, kaon_pz, kaon_E);
  pionmisidentifiedask_p4.SetXYZM(pion_px, pion_py, pion_pz, KAON_MASS-PION_MASS);
  TLorentzVector b_fourmomentum_jpsiconstrained = jpsi_p4 + kaon_p4 + pionmisidentifiedask_p4;
  
  Double_t mass = b_fourmomentum_jpsiconstrained.M();

  return mass;
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

inline auto b_m_omegaconstrained_3body = [](Double_t omega_px, Double_t omega_py, Double_t omega_pz,
                                           Double_t kaon_px, Double_t kaon_py, Double_t kaon_pz, Double_t kaon_E,
                                           Double_t pion_px, Double_t pion_py, Double_t pion_pz, Double_t pion_E)
{
  Double_t OMEGA_MASS = 782.66/1000;
  TLorentzVector omega_p4;
  omega_p4.SetXYZM(omega_px, omega_py, omega_pz, OMEGA_MASS);
  TLorentzVector kaon_p4(kaon_px, kaon_py, kaon_pz, kaon_E);
  TLorentzVector pion_p4(pion_px, pion_py, pion_pz, pion_E);
  TLorentzVector b_fourmomentum_omegaconstrained = omega_p4 + kaon_p4 + pion_p4;
  
  Double_t mass = b_fourmomentum_omegaconstrained.M()*1e3;

  return mass;
};

inline auto b_m_omegaconstrained_3body_masses = [](Double_t omega_px, Double_t omega_py, Double_t omega_pz,
                                                   Double_t kaon_px, Double_t kaon_py, Double_t kaon_pz, Double_t kaon_M,
                                                   Double_t pion_px, Double_t pion_py, Double_t pion_pz, Double_t pion_M)
{
  Double_t OMEGA_MASS = 782.66;
  TLorentzVector omega_p4, kaon_p4, pion_p4;
  omega_p4.SetXYZM(omega_px, omega_py, omega_pz, OMEGA_MASS);
  kaon_p4.SetXYZM(kaon_px, kaon_py, kaon_pz, kaon_M);
  pion_p4.SetXYZM(pion_px, pion_py, pion_pz, pion_M);
  TLorentzVector b_fourmomentum_omegaconstrained = omega_p4 + kaon_p4 + pion_p4;
  
  Double_t mass = b_fourmomentum_omegaconstrained.M();

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

inline auto invariant_mass_couple = [](Double_t particle1_px, Double_t particle1_py, Double_t particle1_pz, Double_t particle1_E,
                                       Double_t particle2_px, Double_t particle2_py, Double_t particle2_pz, Double_t particle2_E)
{
  TLorentzVector particle1_p4, particle2_p4;
  particle1_p4.SetPxPyPzE(particle1_px, particle1_py, particle1_pz, particle1_E); 
  particle2_p4.SetPxPyPzE(particle2_px, particle2_py, particle2_pz, particle2_E); 

  TLorentzVector couple12_p4 = particle1_p4 + particle2_p4;
  Double_t s12 = couple12_p4.Mag2();

  return s12;
};

#endif