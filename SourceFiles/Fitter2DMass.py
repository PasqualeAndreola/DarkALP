import mplhep
import numpy as np

import tensorflow as tf

import uproot
import pandas as pd

import zfit
from zfit import z

mplhep.style.use("LHCb2")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import hist

import calendar
import time

# Spaces (normalization ranges)

# J/Psi Mass Space
m_omega = zfit.Space("x_cons_xpiz_m_best",(600.,1100.), name=r"$m_{DTF}^{\pi^{0}}\left(\omega\right)\left[MeV/c^{2}\right]$")
m_jpsi = zfit.Space("x_cons_xpiz_m_best",(3000,3200.))
m_b = zfit.Space("b_cons_omegaxpiz_m_best", (5000, 5500),name=r"$m_{DTF}^{\pi^{0},\omega}\left(B^{0}\right)\left[MeV/c^{2}\right]$")
m_b_omega = m_b*m_omega
m_diag = zfit.Space("b_cons_omegaxpiz_m_best-x_cons_xpiz_m_best", (0, 6500))

n_bins = 50     #bestvalue 30
n_bins_2d = 20  #bestvalue 10
bins_2d_b_width = m_b.area()/n_bins_2d
bins_2d_omega_width = m_omega.area()/n_bins_2d
n_model_points_per_bin = 50
# Data

# Dataset
#file = uproot.open("DATAafterreader_signalregion_150kEVT.root")
file = uproot.open("InputFiles/DATA_SignalRegion_OmegaMVA_BDT_20072023_OmegaCUT.root")
#file = uproot.open("omegaMCStrippingFiltered.root")
#file = uproot.open("root://eosuser.cern.ch//eos/user/p/paandreo/DATA_signalregion.root")
tree = file["DecayTree"]
varstobeloaded = ["b_cons_omegaxpiz_m_best", "x_cons_xpiz_m_best"]
filter  = "(mvacut_kbdt1>0.4)"
filter += "& ((s_kpi<1.844*1.844*1e6) | (s_kpi>1.884*1.884*1e6))"
filter += "& ((s_kxpi<1.844*1.844*1e6) | (s_kxpi>1.884*1.884*1e6))"
filter += "& ((s_kpi>0.6**2*1e6) & (s_kpi<1**2*1e6))"
#filter += "& (Pi0Merged == 0)"
filter += "& (x_cons_xpiz_m_best<1100) & (x_cons_xpiz_m_best>600)"
filter += "& (b_cons_omegaxpiz_m_best<5500) & (b_cons_omegaxpiz_m_best>5000)"
#filter = ["x_cons_xpiz_m_best>600"]        
dataframe = tree.arrays(varstobeloaded, filter,  library='pd')

print(tree)
data_m_omega = zfit.Data.from_pandas(df=dataframe, obs=m_omega)
data_m_b0 = zfit.Data.from_pandas(df=dataframe, obs=m_b)
data_m_b_omega = zfit.Data.from_pandas(df=dataframe, obs=m_b_omega)
#data_mc_jpsi = zfit.Data.from_pandas(df=dataframe, obs=m_jpsi)
#data_mc_b    = zfit.Data.from_pandas(df=dataframe, obs=m_b)
#data_mc_bjpsi    = zfit.Data.from_pandas(df=dataframe, obs=m_b*m_jpsi)

# Parameters
# using a dict for the composite params
def param_product(params):
    return params["a"] * params["b"]

def param_sum(params):
    return params["a"] + params["b"]

# Shared signal parameters
yield_b_omega = zfit.param.Parameter(r"$n^{B^0-omega}$", 35000, 0, 1e5)
yield_bjpsi = zfit.param.Parameter(r"$n^{B^0-J/\psi}$", 1000, 0, )
yield_b = zfit.param.Parameter(r"$yield^{B^0}$", 80000, 0, )

# omega(782) parameters

# omega(782) signal parameters
# First Crystall Ball parameters
mu_omega_cb1 = zfit.Parameter(r"$\mu^{\omega}_{CB1}$", 783.77, 600, 1000)
sigma_omega_cb1 = zfit.Parameter(r"$\sigma^{\omega}_{CB1}$", 12, 0.1, 70)
alpha_omega_cb1 = zfit.param.ConstantParameter(r"$\alpha^{\omega}_{CB1}$", -1.25)
n_omega_cb1 = zfit.param.ConstantParameter(r"$n^{\omega}_{CB1}$", 3)
# Second Crystall Ball parameters
mu_omega_cb2 = mu_omega_cb1
sigmaratio_omega_cb2 = zfit.param.ConstantParameter(r"$\sigma_{CB2}/\sigma_{CB1}$", 2.43)
sigma_omega_cb2 = zfit.ComposedParameter(r"$\sigma^{\omega}_{CB2}$", param_product, params={"a": sigmaratio_omega_cb2, "b": sigma_omega_cb1})
alpha_omega_cb2 = zfit.param.ConstantParameter(r"$\alpha^{\omega}_{CB2}$", 1.73)
n_omega_cb2 = zfit.param.ConstantParameter(r"$n^{\omega}_{CB2}$", 3)

# Phi(1020) parameters
# First Crystall Ball parameters
deltam_phi_omega = zfit.param.ConstantParameter(r"$\Delta\left(m_{\phi}-m_{\phi}\right)$", 1019.461-782.66)
mu_phi_cb1 = zfit.ComposedParameter(r"$\mu^{\phi}_{CB1}$", param_sum, params={"a": deltam_phi_omega, "b": mu_omega_cb1})
ratiowidth_phi_omega = zfit.param.ConstantParameter(r"$\Gamma_{\phi}/\Gamma_{\omega}$", 4.249/8.68)
sigma_phi_cb1 = zfit.ComposedParameter(r"$\sigma^{\phi}_{CB1}$", param_product, params={"a": ratiowidth_phi_omega, "b": sigma_omega_cb1})
alpha_phi_cb1 = zfit.param.ConstantParameter(r"$\alpha^{\phi}_{CB1}$", -1.25)
n_phi_cb1 = zfit.param.ConstantParameter(r"$n^{\phi}_{CB1}$", 3)
# Second Crystall Ball parameters
mu_phi_cb2 = mu_phi_cb1
sigmaratio_phi_cb2 = zfit.param.ConstantParameter(r"$\sigma_{CB2,\phi}/\sigma_{CB1,\phi}$", 2.43)
sigma_phi_cb2 = zfit.ComposedParameter(r"$\sigma^{\phi}_{CB2}$", param_product, params={"a": sigmaratio_phi_cb2, "b": sigma_phi_cb1})
alpha_phi_cb2 = zfit.param.ConstantParameter(r"$\alpha^{\phi}_{CB2}$", 1.73)
n_phi_cb2 = zfit.param.ConstantParameter(r"$n^{\phi}_{CB2}$", 3)

# omega(782) background parameters
# Combinatorial background with Chebyshev polynomials
coeff_one_cheby = zfit.Parameter(r"$a_{1,Cheby}^{\omega}", 0.5, -100, 100)
coeff_two_cheby = zfit.Parameter(r"$a_{2,Cheby}^{\omega}", 0.5, -100, 100)

# Combinatorial background
lambda_omega_comb = zfit.Parameter(r"$\lambda^{\omega}_{COMB}$", -0.001, -10., 0.)
lambda_phi_comb = zfit.Parameter(r"$\lambda^{\phi}_{COMB}$", -0.13, -1, 0.)
# Part reco 
mu_omega_preco = zfit.Parameter(r"$\mu^{\omega}_{PReco}$", 685, 630, 750)
sigma_omega_preco = zfit.Parameter(r"$\sigma^{\omega}_{PReco}$", 45, 25, 65)
# Part reco High Mass
mu_omega_precohm = zfit.Parameter(r"$\mu^{\omega}_{PRecoHM}$", 950, 900, 1000)
sigma_omega_precohm = zfit.Parameter(r"$\sigma^{\omega}_{PRecoHM}$", 75, 25, 90)

# Yields
n_sig_omega_cb1 = zfit.Parameter(r"$n^{\omega}_{Sig-CB1}$", 1e3, 0, 1e7)
n_sig_omega_cb2 = zfit.Parameter(r"$n^{\omega}_{Sig-CB2}$", 1e3, 0, 1e7)
n_sig_phi_cb1 = zfit.Parameter(r"$n^{\phi}_{Sig-CB1}$", 1e3, 0, 1e7)
n_sig_phi_cb2 = zfit.Parameter(r"$n^{\phi}_{Sig-CB2}$", 1e3, 0, 1e7)
n_bkg_omega_comb = zfit.Parameter(r"$n^{\omega}_{Bkg-Comb}$", 2e3, 0, 1e7)
n_bkg_phi_comb = zfit.Parameter(r"$n^{\phi}_{Bkg-Comb}$", 2e3, 0, 1e7)
n_bkg_omega_preco = zfit.Parameter(r"$n^{\omega}_{Bkg-PReco}$", 2e3, 0, 1e7)
n_bkg_omega_precohm = zfit.Parameter(r"$n^{\omega}_{Bkg-PRecoHM}$", 2e3, 0, 1e7)

####### EXPERIMENTAL ########
#n_sig_omega_frac_cb2overcb1 = zfit.Parameter(r"$n^{\omega}_{Sig-CB2}/n^{\omega}_{Sig-CB1}", 3, 0, 100)
#n_sig_omega_cb2_asafrac = zfit.ComposedParameter(r"$n^{\omega}_{Sig-CB2}$", param_product, params={"a": n_sig_omega_frac_cb2overcb1, "b": n_sig_omega_cb1})
#############################

# J/Psi signal parameters
mu_jpsi = zfit.Parameter(r"$\mu^{J/\psi}$", 3096.9, 3096.9-100., 3096.9+100.)
sigma_jpsi = zfit.Parameter(r"$\sigma^{J/\psi}$", 14., 0., 100.)
alpha_jpsi = zfit.Parameter(r"$\alpha^{J/\psi}$", 1., 0, 6.)
n_jpsi = zfit.param.ConstantParameter(r"$n^{J/\psi}$", 3)
n_jpsi_sig = zfit.Parameter(r"$n_{J/\psi,1}^{Sig}$", 2e4, 0, 5e6)

mu2_jpsi = zfit.Parameter(r"$\mu_2^{J/\psi}$", 3096.9, 3096.9-100., 3096.9+100.)
sigmaratio_jpsi = zfit.Parameter(r"$\sigma_{Ratio}^{J/\psi}$", 3., 0., 100.)
alpha2_jpsi = zfit.Parameter(r"$\alpha_2^{J/\psi}$", -1., -6, 0.)
n2_jpsi = zfit.param.ConstantParameter(r"$n_2^{J/\psi}$", 3)
n2_jpsi_sig = zfit.Parameter(r"$n_{J/\psi,2}^{Sig}$", 2e4, 0, 5e6)

#B signal parameters
mu_b = zfit.Parameter(r"$\mu^{B^0}$", 5284, 5279-20., 5279+20.)
sigma_b = zfit.Parameter(r"$\sigma^{B^0}$", 20., 10., 60.)
alpha_b = zfit.param.ConstantParameter(r"$\alpha^{B^0}$", -1.27)
n_b = zfit.param.ConstantParameter(r"$n^{B^0}$", 3)
n_b_sig = zfit.Parameter(r"$n_{B^0,1}^{Sig}$", 4892, 0, 1e6)

sigmaratio_b = zfit.param.ConstantParameter(r"$\sigma_{Ratio}^{B^0}$", 2.64)
alpha2_b = zfit.param.ConstantParameter(r"$\alpha_2^{B^0}$", 1.93)
n2_b = zfit.param.ConstantParameter(r"$n_2^{B^0}$", 3)
n2_b_sig = zfit.Parameter(r"$n_{B^0,2}^{Sig}}$", 22106, 0, 1e6)

sigmaratio_jpsi_comp = zfit.ComposedParameter('sigmaratio_jpsi_comp', param_product, params={"a": sigmaratio_jpsi, "b": sigma_jpsi})
sigmaratio_b_comp = zfit.ComposedParameter('sigmaratio_b_comp', param_product, params={"a": sigmaratio_b, "b": sigma_b})

#########################################
#                 Models                #
#########################################

# Omega mass

# Signal - Sum of two CBs 
model_m_omega_cb1 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_omega_cb1, sigma=sigma_omega_cb1, alpha=alpha_omega_cb1, n=n_omega_cb1)
model_m_omega_cb1_extended = model_m_omega_cb1.create_extended(n_sig_omega_cb1, name='model_m_omega_cb1')
model_m_omega_cb2 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_omega_cb2, sigma=sigma_omega_cb2, alpha=alpha_omega_cb2, n=n_omega_cb2)
model_m_omega_cb2_extended = model_m_omega_cb2.create_extended(n_sig_omega_cb2, name='model_m_omega_cb2')
omegafrac = zfit.Parameter('omegafrac', 0.1, 0., 100)
model_m_omega_sig = zfit.pdf.SumPDF([model_m_omega_cb1, model_m_omega_cb2], omegafrac)
model_m_omega_sig_ext = zfit.pdf.SumPDF([model_m_omega_cb1_extended, model_m_omega_cb2_extended])
# BKG for Phi - Sum of two CBs
model_m_phi_cb1 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_phi_cb1, sigma=sigma_phi_cb1, alpha=alpha_phi_cb1, n=n_phi_cb1)
model_m_phi_cb1_extended = model_m_phi_cb1.create_extended(n_sig_phi_cb1, name='model_m_phi_cb1')
model_m_phi_cb2 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_phi_cb2, sigma=sigma_phi_cb2, alpha=alpha_phi_cb2, n=n_phi_cb2)
model_m_phi_cb2_extended = model_m_phi_cb2.create_extended(n_sig_phi_cb2, name='model_m_phi_cb2')
phifrac = zfit.Parameter('phifrac', 0.1, 0., 100)
model_m_phi_sig = zfit.pdf.SumPDF([model_m_phi_cb1, model_m_phi_cb2], phifrac)
model_m_phi_sig_ext = zfit.pdf.SumPDF([model_m_phi_cb1_extended, model_m_phi_cb2_extended])

# Background - Combinatorial
model_m_omega_cheby = zfit.pdf.Chebyshev(obs=m_omega, coeffs=[coeff_one_cheby, coeff_two_cheby])
model_m_omega_cheby_extended = model_m_omega_cheby.create_extended(n_bkg_omega_comb, name='model_m_omega_chebycomb')
model_m_omega_comb = zfit.pdf.Exponential(obs=m_omega, lam=lambda_omega_comb)
model_m_omega_comb_extended = model_m_omega_comb.create_extended(n_bkg_omega_comb, name='model_m_omega_comb')
model_m_phi_comb = zfit.pdf.Exponential(obs=m_omega, lam=lambda_phi_comb)
model_m_phi_comb_extended = model_m_phi_comb.create_extended(n_bkg_phi_comb, name='model_m_phi_comb')
model_m_omega_preco = zfit.pdf.Gauss(obs=m_omega, mu=mu_omega_preco, sigma=sigma_omega_preco)
model_m_omega_preco_extended = model_m_omega_preco.create_extended(n_bkg_omega_preco, name='model_m_omega_preco')
model_m_omega_precohm = zfit.pdf.Gauss(obs=m_omega, mu=mu_omega_precohm, sigma=sigma_omega_precohm)
model_m_omega_precohm_extended = model_m_omega_precohm.create_extended(n_bkg_omega_precohm, name='model_m_omega_precohm')
omegabkgfrac_preco = zfit.Parameter('omegabkgfrac_preco', 0.1, 0., 100)
model_m_omega_precotot = model_m_omega_bkg = zfit.pdf.SumPDF([model_m_omega_preco, model_m_omega_precohm], [omegabkgfrac_preco])
omegabkgfrac = zfit.Parameter('omegabkgfrac', 0.1, 0., 100)
model_m_omega_bkg = zfit.pdf.SumPDF([model_m_omega_comb, model_m_omega_preco], [omegabkgfrac])
#model_m_omega_bkg = zfit.pdf.SumPDF([model_m_omega_comb, model_m_omega_preco, model_m_omega_precohm], [omegabkgfrac, omegabkgfrac_preco]) # With pipiw background in high m_omega region
model_m_omega_bkg_ext = zfit.pdf.SumPDF([model_m_phi_sig_ext, model_m_omega_cheby_extended, model_m_omega_preco_extended])
#model_m_omega_bkg_ext = zfit.pdf.SumPDF([model_m_omega_comb_extended, model_m_omega_preco_extended, model_m_omega_precohm_extended])  # With pipiw background in high m_omega region

# Omega mass total fit function
model_m_omega = zfit.pdf.SumPDF([model_m_omega_sig_ext, model_m_omega_bkg_ext])

# J/Psi mass

# Signal - Sum of two CBs
model_m_jpsi_cb1 = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigma_jpsi, alpha=alpha_jpsi, n=n_jpsi)
model_m_jpsi = model_m_jpsi_cb1.create_extended(n_jpsi_sig, name='CrystalBall1')
model2_m_jpsi = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigmaratio_jpsi_comp, alpha=alpha2_jpsi, n=n2_jpsi)
model2_m_jpsi = model2_m_jpsi.create_extended(n2_jpsi_sig, name='CrystalBall2')
frac = zfit.Parameter('frac', 0.1, 0., 100)
model_m_jpsi = zfit.pdf.SumPDF([model_m_jpsi, model2_m_jpsi])

# J/psi mass combinatorial
lamba_jpsi_m_comb = zfit.Parameter(r"$\lambda^{J/\psi}_{Comb}$", -0.1, -1100, 1)
n_jpsi_m_comb = zfit.Parameter(r"$n_{J/\psi}^{Comb}$", 2e6, 0, 10e7)
model_jpsi_m_comb = zfit.pdf.Exponential(obs=m_jpsi, lam=lamba_jpsi_m_comb)
model_jpsi_m_comb = model_jpsi_m_comb.create_extended(n_jpsi_m_comb, name='jpsi_m_combbkg')
frac_jpsi_m_comb = zfit.Parameter('frac_jpsi_m_comb', 0.2, 0., 0.99)

# B mass

# B mass fit model
model_m_b0_cb1 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigma_b, alpha=alpha_b, n=n_b)
model_m_b0_cb1_extended = model_m_b0_cb1.create_extended(n_b_sig, name='B0_CrystalBall1B')
model_m_b0_cb2 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigmaratio_b_comp, alpha=alpha2_b, n=n2_b)
model_m_b0_cb2_extended = model_m_b0_cb2.create_extended(n2_b_sig, name='B0_CrystalBall2B')

model_m_b0_gauss1 = zfit.pdf.Gauss(mu=mu_b, sigma=sigma_b, obs=m_b)
model_m_b0_gauss1_extended = model_m_b0_gauss1.create_extended(n_b_sig, name='B0_Gauss1')
model_m_b0_gauss2 = zfit.pdf.Gauss(mu=mu_b, sigma=sigmaratio_b_comp, obs=m_b)
model_m_b0_gauss2_extended = model_m_b0_gauss2.create_extended(n2_b_sig, name='B0_Gauss2')

b_frac = zfit.Parameter('bfrac', 0.1, 0., 100)
model_m_b0_sig = zfit.pdf.SumPDF([model_m_b0_cb1, model_m_b0_cb2], b_frac)
model_m_b0_sig_ext = zfit.pdf.SumPDF([model_m_b0_cb1_extended, model_m_b0_cb2_extended])

# B mass combinatorial
lambda_b_m_comb = zfit.Parameter(r"$\lambda^{B^0}_{Comb}$", -0.0047, -1100, 1100)
n_b_m_comb = zfit.Parameter(r"$n^{B^0}_{Comb}$", 1.4e4, 0, 1e7)
model_b_m_comb = zfit.pdf.Exponential(obs=m_b, lam=lambda_b_m_comb)
model_b_m_comb_ext = model_b_m_comb.create_extended(n_b_m_comb, name='b_m_Combbkg')
frac_b_m_sig = zfit.Parameter('frac_b_m_sig', 0.2, 0., 0.999)
frac_b_m_comb = zfit.Parameter('frac_b_m_comb', 0.2, 0., 0.999)

# B mass partially reconstructed
kdegridpoints = 128
recofilter = "(B0_M>5) & (B0_M<5.5)"
xwide = np.linspace(5, 5.5, 5)
filepartreco = uproot.open("InputFiles/B2K1omega_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco = zfit.param.Parameter(r"$yield^{partreco}_{K(1)\omega}$", 2500, 0, )
kde_K1omega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
kde_K1omega_ext = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco, num_grid_points=kdegridpoints)

filepartreco = uproot.open("InputFiles/B2Kpigammaomega_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco_kpigammaomega = zfit.param.Parameter(r"$yield^{partreco}_{K\pi\gamma\omega}$", 2500, 0, )
kde_kpigammaomega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_kpigammaomega, num_grid_points=kdegridpoints)

filepartreco = uproot.open("InputFiles/B2Kpipiomega_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco_kpipiomega = zfit.param.Parameter(r"$yield^{partreco}_{K\pi\pi^{0}\omega}$", 2500, 0, )
kde_kpipiomega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_kpipiomega, num_grid_points=kdegridpoints)

filepartreco = uproot.open("InputFiles/B2Dstaromega2Dpiz3pi2Kpi3pi_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco_dstaromega = zfit.param.Parameter(r"$yield^{partreco}_{D^{*}_{\pi^{0}}\omega}$", 2500, 0, )
kde_dstaromega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_dstaromega, num_grid_points=kdegridpoints)

filepartreco = uproot.open("InputFiles/B2Dstaromega2Dgamma3pi2Kpi3pi_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco_dstargammaomega = zfit.param.Parameter(r"$yield^{partreco}_{D^{*}_{\gamma}\omega}$", 2500, 0, )
kde_dstargammaomega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
kde_dstargammaomega_ext = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_dstargammaomega, num_grid_points=kdegridpoints)

filepartreco = uproot.open("InputFiles/B2Domega_tree.root")
treepartreco = filepartreco["DecayTree"]
dataframepartreco = treepartreco.arrays(["b_cons_omegaxpiz_m_best"], recofilter, library='pd')
data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
yield_bpartreco_domega = zfit.param.Parameter(r"$yield^{partreco}_{D\omega}$", 2500, 0, )
kde_domega = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
kde_domega_ext = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_domega, num_grid_points=kdegridpoints)

frac_kde = zfit.Parameter('frac_kde', 0.2, 0., 0.999)
yield_kde = zfit.Parameter('yield_kde', 1e4, 0., 1e5)
#kde = zfit.pdf.SumPDF([kde_K1omega, kde_dstargammaomega, kde_dstaromega, kde_kpipiomega, kde_kpigammaomega, kde_domega])
kde = zfit.pdf.SumPDF([kde_K1omega, kde_dstargammaomega], frac_kde)
kde_ext = kde.create_extended(yield_kde, 'kde_ext')

mu_partiallyreco = zfit.Parameter(r"$\mu_{part}$", 5063, 5000, 5100)
sigma_partiallyreco = zfit.Parameter(r"$\sigma_{part}$", 95.5, 0, 300)
n_b_m_partiallyreco = zfit.Parameter(r"$n^{B^0}_{Part}$", 2.8e4, 0, 1e7)
model_b_m_partiallyreco = zfit.pdf.Gauss(mu=mu_partiallyreco, sigma=sigma_partiallyreco, obs=m_b)
model_b_m_partiallyreco_ext = model_b_m_partiallyreco.create_extended(n_b_m_partiallyreco, name='b_m_PartReco')

# Mass diagonal component
mu_diag = zfit.Parameter(r"$\mu_{diag}$", 3000, 600, 5500)
sigma_diag = zfit.Parameter(r"$\sigma_{diag}$", 95.5, 0, 1000)
n_b_m_diag = zfit.Parameter(r"$n^{B^0}_{diag}$", 2.8e4, 0, 1e7)
model_b_m_diag = zfit.pdf.Gauss(mu=mu_diag, sigma=sigma_diag, obs='b_cons_omegaxpiz_m_best-x_cons_xpiz_m_best')
model_b_m_diag_ext = model_b_m_diag.create_extended(n_b_m_diag, name='b_m_diag')

# 1D mass diagonal component
class DiagonalGaussianPDF(zfit.pdf.ZPDF):
    """ Diagonal component for the mass of the B and the omega

    The function is normalized over a finite range and therefore a PDF.

    Args:
        theta (`zfit.Parameter`): rotation angle for the Gaussian axis
        sigmax (`zfit.Parameter`): width along x axis
        sigmay (`zfit.Parameter`): width along y axis
        x0 (`zfit.Parameter`): x coordinate of the center
        y0 (`zfit.Parameter`): y coordinate of the center
        obs (`zfit.Space`):
        name (str):
        dtype (tf.DType):
    """

    _PARAMS = ['sigma', 'mean']
    _N_OBS = 2

    def _unnormalized_pdf(self, x):
        diagsigma = self.params['sigma']
        diagmean = self.params['mean']

        bmass, omegamass = z.unstack_x(x)

        pdf = (1 / (tf.sqrt(diagsigma))) * tf.exp( - ( tf.square( (bmass - omegamass)-diagmean ) / (2 * tf.square(diagsigma) ) ) )
        #pdf = bmass - omegamass

        return pdf
    
# 2D mass diagonal component
# Defining the custom model for the 2D Gaussian
class Diagonal2DGaussianPDF(zfit.pdf.ZPDF):
    """ Diagonal component for the mass of the B and the omega

    The function is normalized over a finite range and therefore a PDF.

    Args:
        theta (`zfit.Parameter`): rotation angle for the Gaussian axis
        sigmax (`zfit.Parameter`): width along x axis
        sigmay (`zfit.Parameter`): width along y axis
        x0 (`zfit.Parameter`): x coordinate of the center
        y0 (`zfit.Parameter`): y coordinate of the center
        obs (`zfit.Space`):
        name (str):
        dtype (tf.DType):
    """

    _PARAMS = ['theta', 'sigmax', 'sigmay', 'x0', 'y0']
    _N_OBS = 2

    def _unnormalized_pdf(self, x):
        theta = self.params['theta']
        sigmax = self.params['sigmax']
        sigmay = self.params['sigmay']
        x0 = self.params['x0']
        y0 = self.params['y0']

        bmass, omegamass = z.unstack_x(x)

        a = tf.square(tf.cos(theta)) / (2*tf.square(sigmax)) + tf.square(tf.sin(theta)) / (2*tf.square(sigmay))
        b = - tf.sin(2*theta) / (4*tf.square(sigmax)) + tf.sin(2*theta) / (4*tf.square(sigmay))
        c = tf.square(tf.sin(theta)) / (2*tf.square(sigmax)) + tf.square(tf.cos(theta)) / (2*tf.square(sigmay))

        pdf = 1 / tf.sqrt((sigmax+sigmay)) * tf.exp( - ( a*tf.square(bmass-x0) + 2*b*(bmass-x0)*(omegamass-y0) + c*tf.square(omegamass-y0) ) )

        return pdf

class Diagonal2DGaussianAmplitudePDF(zfit.pdf.ZPDF):
    """ Diagonal component for the mass of the B and the omega

    The function is normalized over a finite range and therefore a PDF.

    Args:
        theta (`zfit.Parameter`): rotation angle for the Gaussian axis
        sigmax (`zfit.Parameter`): width along x axis
        sigmay (`zfit.Parameter`): width along y axis
        x0 (`zfit.Parameter`): x coordinate of the center
        y0 (`zfit.Parameter`): y coordinate of the center
        obs (`zfit.Space`):
        name (str):
        dtype (tf.DType):
    """

    _PARAMS = ['A', 'theta', 'sigmax', 'sigmay', 'x0', 'y0']
    _N_OBS = 2

    def _unnormalized_pdf(self, x):
        amplitude = self.params['A']
        theta = self.params['theta']
        sigmax = self.params['sigmax']
        sigmay = self.params['sigmay']
        x0 = self.params['x0']
        y0 = self.params['y0']

        bmass, omegamass = z.unstack_x(x)

        a = tf.square(tf.cos(theta)) / (2*tf.square(sigmax)) + tf.square(tf.sin(theta)) / (2*tf.square(sigmay))
        b = - tf.sin(2*theta) / (4*tf.square(sigmax)) + tf.sin(2*theta) / (4*tf.square(sigmay))
        c = tf.square(tf.sin(theta)) / (2*tf.square(sigmax)) + tf.square(tf.cos(theta)) / (2*tf.square(sigmay))

        pdf = amplitude*tf.exp( - ( a*tf.square(bmass-x0) + 2*b*(bmass-x0)*(omegamass-y0) + c*tf.square(omegamass-y0) ) )

        return pdf
    
mass_diaggauss = zfit.Parameter(r"$m_{B^0}-m_{#omega}^{Diag}$", 3993, 0, 5000)
sigma_diaggauss = zfit.Parameter(r"$#sigma_{DiagGauss}$", 99, 0.1, 200.)
amplitude_diaggauss = zfit.Parameter(r"$A_{DiagGauss}$", 50, 0.1, 200)
theta_diaggauss = zfit.Parameter(r"$\theta_{DiagGauss}$", np.pi/4, -np.pi, np.pi)
sigma_b_diaggauss = zfit.Parameter(r"$DiagGauss-\sigma^{B^0}$", 100, 0.1, 500.)
sigma_omega_diaggauss = zfit.Parameter(r"$DiagGauss\sigma^{\omega}$", sigma_phi_cb1.value(), 0.1, 50.)
diag_gauss_pdf = DiagonalGaussianPDF(obs=m_b_omega, sigma=sigma_diaggauss, mean=mass_diaggauss)
#diag_gauss_pdf = Diagonal2DGaussianPDF(obs=m_b_omega, theta=theta_diaggauss, sigmax=sigma_b_diaggauss, sigmay=sigma_omega_diaggauss, x0=mu_b, y0=mu_phi_cb1)
#diag_gauss_pdf = Diagonal2DGaussianAmplitudePDF(obs=m_b_omega, A=amplitude_diaggauss, theta=theta_diaggauss, sigmax=sigma_b_diaggauss, sigmay=sigma_omega_diaggauss, x0=mu_b, y0=mu_phi_cb1)

''''def integral_func(limits, norm_range, params, model):

    b = params['b']
    c = params['c']

    lower, upper = limits.limit1d
    lower = z.convert_to_tensor(lower)  # the limits are now 1-D, for axis 1
    upper = z.convert_to_tensor(upper)

    # calculate the integral
    integral = cdf_poly(upper, b, c) - cdf_poly(lower, b, c)
    print("Integral called")
    return integral'''
#diag_gauss_pdf.register_analytic_integral(func=, limits=m_b_omega)

#1D model
model_m_b = zfit.pdf.SumPDF([model_m_b0_sig_ext, model_b_m_comb_ext, kde_ext])
#model_m_b = zfit.pdf.SumPDF([model_m_b0_sig, model_b_m_comb, kde], [frac_b_m_sig, frac_b_m_comb])
#model_m_b = model_m_b.create_extended(yield_b)

#2D model
#model_m_b_omega_sig = zfit.pdf.ProductPDF([model_m_b0_sig_ext, model_m_omega_sig_ext])
#model_m_b_omega_sig = model_m_b_omega_sig.create_extended(yield_b_omega)
#model_m_b_omega_bkg = zfit.pdf.ProductPDF([model_m_omega_comb_extended, model_b_m_comb_ext, model_b_m_partiallyreco])
#model_m_b_omega_bkg = model_m_b_omega_bkg.create_extended(yield_bjpsi)
n_m_b_omega_sig = zfit.Parameter(r"$n^{B^0-\omega}_{2DSIG}$", 0.24, 1e-4, 1)
n_m_b_omega_comb = zfit.Parameter(r"$n^{B^0-\omega}_{2DCOMB}$", 0.026, 1e-4, 1)
n_m_breco_omegasig = zfit.Parameter(r"$n^{B^0-\omega}_{2DRECO-SIG}$", 0.47, 1e-4, 1)
n_m_breco_omegacomb = zfit.Parameter(r"$n^{B^0-\omega}_{2DRECO-COMB}$", 0.1, 1e-4, 1)
n_m_bsig_omegacomb = zfit.Parameter(r"$n^{B^0-\omega}_{2DSIG-COMB}$", 0.01, 1e-4, 1)
n_m_bcomb_omegasig = zfit.Parameter(r"$n^{B^0-\omega}_{2DCOMB-SIG}$", 0.1, 1e-4, 1)
n_m_bcomb_diag = zfit.Parameter(r"$n^{B^0-\omega}_{Diag}$", 0.1, 1e-4, 1)
n_b_comb_per_omega_phi = zfit.Parameter(r"$n^{B^0-\phi}_{COMB-SIG}$", 0.1, 1e-4, 1)
model_m_b_omega = zfit.pdf.SumPDF([model_m_omega_sig*model_m_b0_sig,
                                   model_m_omega_cheby*model_b_m_comb,
                                   model_m_omega_sig*model_b_m_partiallyreco,
                                   model_m_omega_cheby*model_b_m_partiallyreco,
                                   model_m_omega_cheby*model_m_b0_sig,
                                   model_m_omega_sig*model_b_m_comb],
                                   [n_m_b_omega_sig, n_m_b_omega_comb, n_m_breco_omegasig, n_m_breco_omegacomb, n_m_bsig_omegacomb])
#model_m_b_omega = zfit.pdf.SumPDF([model_m_omega_sig*model_m_b0_sig,
#                                   model_m_omega_sig*model_b_m_partiallyreco,
#                                   model_m_omega_comb*model_b_m_partiallyreco,
#                                   model_m_omega_comb*model_m_b0_sig],
#                                   [n_m_b_omega_sig, n_m_breco_omegasig, n_m_breco_omegacomb])
""" model_m_b_omega = zfit.pdf.SumPDF([model_m_omega_sig*model_m_b0_sig,
                                   model_m_omega_comb*model_b_m_comb,
                                   model_m_omega_sig*model_b_m_partiallyreco,
                                   model_m_omega_comb*model_b_m_partiallyreco,
                                   model_m_omega_comb*model_m_b0_sig],
                                   [n_m_b_omega_sig, n_m_b_omega_comb, n_m_breco_omegasig, n_m_breco_omegacomb]) """
b_sig_per_omega_sig = zfit.pdf.ProductPDF([model_m_omega_sig, model_m_b0_sig], name = r"$B_{Sig}\cdot\omega_{Sig}$")
b_comb_per_omega_comb = zfit.pdf.ProductPDF([model_m_omega_cheby, model_b_m_comb], name = r"$B_{Comb}\cdot\omega_{Comb}$")
b_preco_per_omega_sig = zfit.pdf.ProductPDF([model_m_omega_sig, kde], name = r"$B_{PReco}\cdot\omega_{Sig}$")
b_preco_per_omega_comb = zfit.pdf.ProductPDF([model_m_omega_cheby, kde], name = r"$B_{PReco}\cdot\omega_{Comb}$")
b_sig_per_omega_comb = zfit.pdf.ProductPDF([model_m_omega_cheby, model_m_b0_sig], name = r"$B_{Sig}\cdot\omega_{Comb}$")
b_tot_per_omega_phi = zfit.pdf.ProductPDF([model_m_phi_sig, diag_gauss_pdf], name = r"$B_{Diag}\cdot\phi_{Sig}$")
b_sig_per_omega_preco = zfit.pdf.ProductPDF([model_m_omega_preco, model_m_b0_sig], name = r"$B_{Sig}\cdot\omega_{PReco}$")
b_sig_per_omega_precohm = zfit.pdf.ProductPDF([model_m_omega_precohm, model_m_b0_sig], name = r"$B_{Sig}\cdot\omega_{PRecoHM}$")
b_sig_per_omega_precotot = zfit.pdf.ProductPDF([model_m_omega_precotot, model_m_b0_sig], name = r"$B_{Sig}\cdot\omega_{PRecoTOT}$")
b_comb_per_omega_sig = zfit.pdf.ProductPDF([model_m_omega_sig, model_b_m_comb], name = r"$B_{Comb}\cdot\omega_{Sig}$")
b_comb_per_omega_phi = zfit.pdf.ProductPDF([model_m_phi_sig, model_b_m_comb], name = r"$B_{Comb}\cdot\phi{Sig}$")
b_comb_per_omega_preco = zfit.pdf.ProductPDF([model_m_omega_preco, model_b_m_comb], name = r"$B_{Comb}\cdot\omega_{PReco}$")
b_kde_per_omega_preco = zfit.pdf.ProductPDF([model_m_omega_preco, kde], name = r"$B_{KDE}\cdot\omega_{PReco}$")
b_per_omega_preco = zfit.pdf.ProductPDF([model_m_omega_preco, model_m_b], name = r"$B_{TOT}\cdot\omega_{PReco}$")
model_m_b_omega = zfit.pdf.SumPDF([b_sig_per_omega_sig,
                                   b_sig_per_omega_comb,
                                   b_comb_per_omega_comb,
                                   b_comb_per_omega_sig,
                                   b_preco_per_omega_comb,
                                   b_tot_per_omega_phi
                                   ],
                                   [n_m_b_omega_sig, n_m_b_omega_comb, n_m_breco_omegasig, n_m_breco_omegacomb, n_m_bcomb_diag])
""" model_m_b_omega = zfit.pdf.SumPDF([model_m_omega_sig_ext*model_m_b0_sig_ext,
                                   model_m_omega_comb_extended*model_b_m_comb_ext,
                                   model_m_omega_comb_extended*model_m_b0_sig_ext],
                                   [n_m_b_omega_sig, n_m_breco_omegacomb])  """
model_m_b_omega = model_m_b_omega.create_extended(yield_b_omega)

# Fit JPsi

## Stage 1: create an unbinned likelihood with the given PDF and dataset
#nll_m_jpsi = zfit.loss.ExtendedUnbinnedNLL(model=model_bm_jpsim, data=data_mc_bjpsi)
#
## Stage 2: instantiate a minimiser (in this case a basic minuit)
#minimizer = zfit.minimize.Minuit()
#
## Stage 3: minimise the given negative log-likelihood
#result = minimizer.minimize(nll_m_jpsi)
#
## Computing the errors
#param_errors = result.hesse()
#
## Information on all the parameters in the fit
#params = result.params
#print(params)

# Fit B

# Stage 1: create an unbinned likelihood with the given PDF and dataset
nll_m_omega = zfit.loss.ExtendedUnbinnedNLL(model=model_m_omega, data=data_m_omega)
nll_m_b = zfit.loss.ExtendedUnbinnedNLL(model=model_m_b, data=data_m_b0)

# Stage 2: instantiate a minimiser (in this case a basic minuit)
minimizer = zfit.minimize.Minuit()

# Stage 3: minimise the given negative log-likelihood
results_omega = minimizer.minimize(nll_m_omega)
current_GMT = time.gmtime()
time_stamp = calendar.timegm(current_GMT)
print("Current timestamp - Before Minuit:", time_stamp)
results_b = minimizer.minimize(nll_m_b)

current_GMT = time.gmtime()
time_stamp = calendar.timegm(current_GMT)
print("Current timestamp - After Minuit:", time_stamp)
# Stage 4: computing the errors
params_errors_omega = results_omega.hesse(method='approx')
current_GMT = time.gmtime()
time_stamp = calendar.timegm(current_GMT)
print("Current timestamp - Before Hesse:", time_stamp)
params_errors_b = results_b.hesse(method='approx')
current_GMT = time.gmtime()
time_stamp = calendar.timegm(current_GMT)
print("Current timestamp - After Hesse:", time_stamp)
# Saving Pars Infos
paramsomega = results_omega.params
paramsb = results_b.params

# Printing the results
print("Function minimum:", results_omega.fmin)
print("Converged:", results_omega.converged)
print("Full minimizer information:", results_omega)

print("Function minimum:", results_b.fmin)
print("Converged:", results_b.converged)
print("Full minimizer information:", results_b)

# Information on all the parameters in the fit
#print(paramsb[mu_b].keys)
# Printing information on specific parameters, e.g. mu
#print("mu_jpsi={}".format(paramsb[mu_b]['value']))

# Plotting

def plot_fit(dat: np.ndarray, basis: np.ndarray, model: np.ndarray, 
             obs: zfit.Space, params: zfit.minimizers.fitresult.ParamHolder, nbins : int=50, smodel: np.ndarray=None,
             drawstyle: str='default', zmodel: zfit.pdf.BasePDF=None, parlist: bool=True, pullplot : bool=True):
    """
    quick plotting function to visualise data and model. 
    Takes:
     - dat: (array) the data that are fitted
     - basis: (array) the points at which the model is evaluated
     - model: (array) the model that describes the data
     - obs: (zfit Space) the space in which the model lives
     - nbins: (int) the number of bins for the data histogram
     - smodel: (array) uncertainty on model (not needed)
     - drawstyle: (str) the drawstyle of plt.plot
     - zmodel: (BasePDF) for drawing submodels
    Returns:
     - None
    """
    # for normalising the pdf, scaled pdf = pdf * yield * area / bins
    limits = obs.limits 
    area = obs.area().numpy()
    area2d = m_b.area()*m_omega.area()
    # data in histogram over the full observable space
    histo = hist.Hist(hist.axis.Regular(nbins, *limits))
    histo.fill(dat)
    chi2_model_x_np = histo.axes.centers[0]
    chi2_model_np = zmodel.pdf(chi2_model_x_np).numpy() * (zmodel.get_yield().numpy()) *area/nbins
    chi2 = sum( ((histo.counts()-chi2_model_np))**2/chi2_model_np)
    print(r"$\chi^2=$", chi2)
    print("Normalized ",r"$\chi^2=$", chi2/len(params))

    # the figure with an errorbar for the data and a line for the model
    fig, (ax, pullsax) = plt.subplots(2,1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

    art_data = ax.errorbar(histo.axes.centers[0], histo.values(), 
                           xerr=histo.axes.widths[0]/2,
                           yerr=np.sqrt(histo.values()), fmt='.', 
                           label='Data', color='black', zorder=10)
    art_model = ax.plot(basis, model * area/nbins, color='darkturquoise', 
                        label='Model', zorder=8, drawstyle=drawstyle)[0]
    
    # if we have the uncertainty on the model we draw it as contour
    # and update the artist for the legend to reflect on the new model
    if smodel is not None:
        _art = ax.fill_between(basis, (model+smodel)*area/nbins, 
                               (model-smodel)*area/nbins, color='darkturquoise', 
                               alpha=0.5, zorder=-2)
        art_model = (art_model, _art)

    # define artists and labels for the legend
    artists = [art_data, art_model]
    labels = ['Data', 'Model']
    # if we want to plot the submodels of our model, we can iterate through
    # all of them and evaluate them at our basis. We will not bootstrap
    # all of their shape uncertainties though, this is just an illustration
    if hasattr(zmodel, 'get_models'):
        nmodels = len(zmodel.get_models())
        cmap = plt.get_cmap('turbo') # you can choose whatever you like. 
        norm = mpl.colors.Normalize(0, nmodels) # create a norm for the cmap
        pdfs = [(m.pdf(basis)*m.get_yield()).numpy()*area/nbins
                for m in zmodel.get_models()]
        names = [m.name.replace('_extended','') for m in zmodel.get_models()]
        labels.extend(names)
        for mdex, pdf in enumerate(pdfs):
            artists.append(ax.plot(basis, pdf, color=cmap(norm(mdex)), 
                                   linestyle='--', zorder=-1)[0])
        
    # legend and axis labels
    ax.legend(artists, labels, loc=1, 
              title='LHCb Preliminary\nData Entries: {0}'.format(len(dat)), title_fontsize=18, fontsize=18)
    ax.set_xlabel(obs.name)
    ax.set_ylabel('Counts [a.u.]')
    parliststring = 'List of parameters:'
    for ipar, par in enumerate(params):
        if ipar == 0:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['hesse']['error'])
        else:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['hesse']['error'])
        width = ax.get_xlim()[1] - ax.get_xlim()[0]
        height = ax.get_ylim()[1] - ax.get_ylim()[0]
    parliststring = parliststring+'\n'+r'$\chi^2={0:.2f}$'.format(chi2/len(params))
    props = dict(boxstyle='square', facecolor='white', alpha=1)

    # place a text box in upper left in axes coords
    if parlist==True:
        ax.text(0.0046, 0.493, parliststring, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        
    if pullplot == True:
        data_model_distance = histo.counts()-chi2_model_np
        pulls = data_model_distance/np.sqrt(histo.counts())
        pullserr = np.sqrt(histo.counts())
        # draw the +/- 2.0 horizontal lines
        [pullsax.axhline(y=i, color="red") for i in [-3, 3]]
        # draw the +/- 1.0 horizontal lines
        [pullsax.axhline(y=i, linestyle="dashdot", color="grey") for i in [-2, 2]]
        # draw a horizontal line at pull=0.0
        [pullsax.axhline(y=i, linestyle="dashed", color="black") for i in [0]]
        
        pullsax.scatter(chi2_model_x_np, pulls, s=20, color="black")
        pullsax.errorbar(chi2_model_x_np,
                               pulls,
                               color="black",
                               xerr=0,
                               yerr=pullserr/histo.counts(),
                               marker=".",
                               fmt="none",
                        ) 
        pullsax.set_ylabel('Pulls')
    return fig

def plot_projection_fit2d(dat: np.ndarray, basis: np.ndarray, model: np.ndarray, 
             obs: zfit.Space, obs2beintegrated : zfit.Space, params: zfit.minimizers.fitresult.ParamHolder, nbins : int=50, smodel: np.ndarray=None,
             drawstyle: str='default', zmodel: zfit.pdf.BasePDF=None, parlist: bool=True, pullplot : bool=True):
    """
    quick plotting function to visualise data and model. 
    Takes:
     - dat: (array) the data that are fitted
     - basis: (array) the points at which the model is evaluated
     - model: (array) the model that describes the data
     - obs: (zfit Space) the space in which the model lives
     - nbins: (int) the number of bins for the data histogram
     - smodel: (array) uncertainty on model (not needed)
     - drawstyle: (str) the drawstyle of plt.plot
     - zmodel: (BasePDF) for drawing submodels
    Returns:
     - None
    """
    # for normalising the pdf, scaled pdf = pdf * yield * area / bins
    limits = obs.limits 
    area = obs.area().numpy()
    obs_np = np.linspace(*obs.limit1d, nbins)

    # data in histogram over the full observable space
    histo = hist.Hist(hist.axis.Regular(nbins, *limits))
    histo.fill(dat)
    chi2_model_x_np = histo.axes.centers[0]
    chi2_model_np = zmodel.create_projection_pdf(limits=obs2beintegrated).pdf(obs_np).numpy()*(zmodel.get_yield().numpy())*area/n_bins
    chi2 = sum( ((histo.counts()-chi2_model_np))**2/chi2_model_np)
    print(r"$\chi^2=$", chi2)
    print("Normalized ",r"$\chi^2=$", chi2/len(params))

    # the figure with an errorbar for the data and a line for the model
    fig, (ax, pullsax) = plt.subplots(2,1, gridspec_kw={'height_ratios': [4, 1]}, sharex=True)

    art_data = ax.errorbar(histo.axes.centers[0], histo.values(), 
                           xerr=histo.axes.widths[0]/2,
                           yerr=np.sqrt(histo.values()), fmt='.', 
                           label='Data', color='black', zorder=10)
    art_model = ax.plot(basis, model * area/nbins, color='darkturquoise', 
                        label='Model', zorder=8, drawstyle=drawstyle)[0]
    
    # if we have the uncertainty on the model we draw it as contour
    # and update the artist for the legend to reflect on the new model
    if smodel is not None:
        _art = ax.fill_between(basis, (model+smodel)*area/nbins, 
                               (model-smodel)*area/nbins, color='darkturquoise', 
                               alpha=0.5, zorder=-2)
        art_model = (art_model, _art)

    # define artists and labels for the legend
    artists = [art_data, art_model]
    labels = ['Data', 'Model']
    # if we want to plot the submodels of our model, we can iterate through
    # all of them and evaluate them at our basis. We will not bootstrap
    # all of their shape uncertainties though, this is just an illustration
    if hasattr(zmodel, 'get_models'):
        pdfs = []
        for sumcomponent, sumfraction in zip(zmodel.pdfs, zmodel.params.values()):
            print(sumfraction)
            pdfs.append(sumcomponent.create_projection_pdf(limits=obs2beintegrated).pdf(obs_np)*(zmodel.get_yield())*area/n_bins*sumfraction)
        names = [m.name.replace('_extended','') for m in zmodel.get_models()]
        nmodels = len(pdfs)
        cmap = plt.get_cmap('turbo') # you can choose whatever you like. 
        norm = mpl.colors.Normalize(0, nmodels) # create a norm for the cmap
        labels.extend(names)
        for mdex, pdf in enumerate(pdfs):
            artists.append(ax.plot(obs_np, pdf, color=cmap(norm(mdex)), 
                                   linestyle='--', zorder=-1)[0])
            #plt.xlim([600, 1000])
        
    # legend and axis labels
    ax.legend(artists, labels, loc=1, 
              title='LHCb Preliminary\nData Entries: {0}'.format(len(dat)), title_fontsize=18, fontsize=18)
    ax.set_xlabel(obs.name)
    ax.set_ylabel('Counts [a.u.]')
    parliststring = 'List of parameters:'
    for ipar, par in enumerate(params):
        if ipar == 0:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['value'])
        else:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['value'])
        width = ax.get_xlim()[1] - ax.get_xlim()[0]
        height = ax.get_ylim()[1] - ax.get_ylim()[0]
    #parliststring = parliststring+'\n'+r'$\chi^2={0:.2f}$'.format(chi2/len(params))
    props = dict(boxstyle='square', facecolor='white', alpha=1)

    # place a text box in upper left in axes coords
    if parlist==True:
        ax.text(0.0046, 0.493, parliststring, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        
    if pullplot == True:
        data_model_distance = histo.counts()-chi2_model_np
        pulls = data_model_distance/np.sqrt(histo.counts())
        pullserr = np.sqrt(histo.counts())
        # draw the +/- 2.0 horizontal lines
        [pullsax.axhline(y=i, color="red") for i in [-3, 3]]
        # draw the +/- 1.0 horizontal lines
        [pullsax.axhline(y=i, linestyle="dashdot", color="grey") for i in [-2, 2]]
        # draw a horizontal line at pull=0.0
        [pullsax.axhline(y=i, linestyle="dashed", color="black") for i in [0]]
        
        pullsax.scatter(chi2_model_x_np, pulls, s=20, color="black")
        pullsax.errorbar(chi2_model_x_np,
                               pulls,
                               color="black",
                               xerr=0,
                               yerr=pullserr/histo.counts(),
                               marker=".",
                               fmt="none",
                        ) 
        pullsax.set_ylabel('Pulls')
    return fig

data_m_omega_np = np.sort(data_m_omega["x_cons_xpiz_m_best"].numpy())
data_m_b_np = np.sort(data_m_b0["b_cons_omegaxpiz_m_best"].numpy())
#data_mc_jpsi_np = np.sort(data_mc_jpsi["b_cons_Jpsixpiz_m_best"].numpy())
m_jpsi_area = m_jpsi.area().numpy()
m_b_area = m_b.area().numpy()
#model_x_np = np.linspace(*m_jpsi.limit1d, n_model_points_per_bin*n_bins)
#model_np = model_m_b.pdf(model_x_np).numpy() * (model_m_b.get_yield().numpy())
#fig = plot_fit(data_mc_jpsi_np, model_x_np, model_np, [m_b, m_jpsi], params, nbins=n_bins, zmodel=model_m_b)
#fig.savefig('test1.pdf')
#plt.show()

model_omega_np = np.linspace(*m_omega.limit1d, n_model_points_per_bin*n_bins)
model_np = model_m_omega.pdf(model_omega_np).numpy() * (model_m_omega.get_yield().numpy())
fig_omega = plot_fit(data_m_omega_np, model_omega_np, model_np, m_omega, paramsomega, nbins=n_bins, zmodel=model_m_omega, parlist=False)
fig_omega.savefig('omega1dfit.pdf')
fig_omega.savefig('omega1dfit.png')
fig_omega.savefig("OutputFiles/PNGPlots/FitResults/omega1dfit.png")
                  
model_b_np = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins)
model_np = model_m_b.pdf(model_b_np).numpy() * (model_m_b.get_yield().numpy())
fig_b = plot_fit(data_m_b_np, model_b_np, model_np, m_b, paramsb, nbins=n_bins, zmodel=model_m_b, parlist=False)
fig_b.savefig('b1dfit.pdf')
fig_b.savefig('b1dfit.png') 
fig_b.savefig("OutputFiles/PNGPlots/FitResults/b1dfit.png")

b_omega_constraints = []
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=frac_kde, observation=frac_kde.value(), uncertainty=frac_kde.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=phifrac, observation=phifrac.value(), uncertainty=phifrac.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=lambda_omega_comb, observation=lambda_omega_comb.value(), uncertainty=lambda_omega_comb.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=mu_omega_preco, observation=mu_omega_preco.value(), uncertainty=mu_omega_preco.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_omega_preco, observation=sigma_omega_preco.value(), uncertainty=sigma_omega_preco.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=mu_omega_precohm, observation=mu_omega_precohm.value(), uncertainty=mu_omega_precohm.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_omega_precohm, observation=sigma_omega_precohm.value(), uncertainty=sigma_omega_precohm.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_omega_cb1, observation=sigma_omega_cb1.value(), uncertainty=sigma_omega_cb1.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_omega_cb2, observation=sigma_omega_cb2.value(), uncertainty=sigma_omega_cb2.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=lambda_b_m_comb, observation=lambda_b_m_comb.value(), uncertainty=lambda_b_m_comb.value().numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=yield_bpartreco, observation=yield_bpartreco.value(), uncertainty=yield_bpartreco.numpy()*1e-2))
b_omega_constraints.append(zfit.constraint.GaussianConstraint(params=yield_bpartreco_dstargammaomega, observation=yield_bpartreco_dstargammaomega.value(), uncertainty=yield_bpartreco_dstargammaomega.numpy()*1e-2))
nll_m_b_omega = zfit.loss.ExtendedUnbinnedNLL(model=model_m_b_omega, data=data_m_b_omega, constraints=b_omega_constraints)
results_b_omega = minimizer.minimize(nll_m_b_omega)
paramsbomega = results_b_omega.params
params_errors_b_omega = results_b_omega.hesse()

print("Function minimum:", results_b_omega.fmin)
print("Converged:", results_b_omega.converged)
print("Full minimizer information:", results_b_omega)

model_omega_m_projection = model_m_b_omega.create_projection_pdf(limits=m_b)
model_b_omega_np_omegacounts = model_omega_m_projection.pdf(model_omega_np).numpy()*(yield_b_omega)
fig_omega_projection = plot_projection_fit2d(data_m_omega_np, model_omega_np, model_b_omega_np_omegacounts, m_omega, m_b, paramsbomega, nbins=n_bins, zmodel=model_m_b_omega, parlist=False)
fig_omega_projection.savefig('omega_2dfit_projection.pdf')
fig_omega_projection.savefig('OutputFiles/PNGPlots/FitResults/omega_2dfit_projection.png')

model_b_m_projection = model_m_b_omega.create_projection_pdf(limits=m_omega)
model_b_omega_np_bcounts = model_b_m_projection.pdf(model_b_np).numpy()*(yield_b_omega)
fig_b_projection = plot_projection_fit2d(data_m_b_np, model_b_np, model_b_omega_np_bcounts, m_b, m_omega, paramsbomega, nbins=n_bins, zmodel=model_m_b_omega, parlist=False)
fig_b_projection.savefig('b_2dfit_projection.pdf')
fig_b_projection.savefig('OutputFiles/PNGPlots/FitResults/b_2dfit_projection.png')

""" class AdaptBin:
    AdaptBin(TString name, int nbin, double xmin, double xmax, double ymin, double ymax)
        : name_(name), xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
        xs_(0), ys_(0), nentries_(0),
        loaded_(false), usingLists_(false), verbose_(false),
        theHisto_(0) {calcDivisions(nbin);}
    ~AdaptBin() {
        delete[] xs_;
        delete[] ys_;
        if(theHisto_) delete theHisto_;
    }
    bool loadDataFromTree(TString fname, TString tname, TString xname, TString yname);
    bool loadDataFromHist(TString fname, TString hname);
    bool loadDataFromHist(TH2* hist);
    bool addEntry(double x, double y);
    TH2Poly* getHisto(TString name="");
    void setVerbose(bool verbose=true) {verbose_ = verbose;}
	void calcDivisions(int ntarget)
	bool loadDataFromLists()
	void initHisto(double xmin, double xmax, double ymin, double ymax, uint iter=0)
	TString name_
	double xmin_, xmax_, ymin_, ymax_
	std::list<double> xList_, yList_
	double *xs_, *ys_
	std::vector<int> divisions_
	int nentries_
	bool loaded_, usingLists_, verbose_
	
 """
""" nentries_ = tree.size()
divisions_=[]
uproot.classes.
def initHisto(xmin : float, xmax : float, ymin: float , ymax : float , iter : int):

	#If it's the last iteration create the bin and return
	if(iter == divisions_.size()):
		x_new = np.zeros(5)
		y_new = np.zeros(5)
		x_new[0] = xmin; x_new[1] = xmin; x_new[2] = xmax; x_new[3] = xmax; x_new[4] = xmin
		y_new[0] = ymin; y_new[1] = ymax; y_new[2] = ymax; y_new[3] = ymin; y_new[4] = ymin
		theHisto_->AddBin(5, x_new, y_new)
		if(verbose_) std::cout << "INFO in AdaptBin::initHisto : Adding bin from (" << xmin << "," << ymin << ") to (" << xmax << "," << ymax << ")." << std::endl
		return
	

	#If not the last iteration then divide the bin
	n_divx=divisions_[iter]
	n_divy=divisions_[iter]

	#if(verbose_) std::cout << "INFO in AdaptBin::initHisto : Dividing bin from (" << xmin << "," << ymin << ") to (" << xmax << "," << ymax << ") into " << n_divx << " by " << n_divy << " subbins" << std::endl;

	xIn = np.zeros(nentries_)
	yIn = np.zeros(nentries_)
	xIndex = np.zeros(nentries_+2)
	yIndex = np.zeros(nentries_+2)

	int xCountIn = 0; 
	for(int i = 0; i<nentries_; ++i) {
		if ((xs_[i]<xmin)||(xs_[i]>xmax)||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
		xIn[xCountIn] = xs_[i];
		++xCountIn;
	}

	//find the delimitting x and y values for the sub bins
	double xLimits[n_divx + 1];
	double yLimits[n_divx][n_divy + 1];

	//first sort entries in x and divide bin into equally populated bins in x
	TMath::Sort(xCountIn, xIn, xIndex, false);

	xLimits[0] = xmin;
	xLimits[n_divx] = xmax;
	for (int nDivx = 0; nDivx < n_divx; ++nDivx){
		if (nDivx  < (n_divx-1)){
			if(xCountIn>0) {
				xLimits[nDivx+1] = xIn[xIndex[xCountIn*(nDivx+1)/n_divx]];
			} else {
				//if no entries then use equal bin widths
				xLimits[nDivx+1] = xmin + (xmax-xmin)*(nDivx+1)/n_divx;
			}
		}

		//for each bin in x divide into equally populated bins in y
		yLimits[nDivx][0] = ymin;
		yLimits[nDivx][n_divy] = ymax;    
		int yCountIn = 0;

		for(int i = 0; i<nentries_; ++i) {
			if ((xs_[i]<xmin)||(xs_[i]>xmax)||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
			if ((xs_[i]<xLimits[nDivx])||(xs_[i]>=xLimits[nDivx+1])||(ys_[i]<ymin)||(ys_[i]>ymax)) continue;
			yIn[yCountIn] = ys_[i];
			++yCountIn;
		}

		TMath::Sort(yCountIn, yIn, yIndex, false);

		for (int nDivy = 1; nDivy < n_divy; ++nDivy){
			if(yCountIn>0) {
				yLimits[nDivx][nDivy] = yIn[yIndex[yCountIn*nDivy/n_divy]];
			} else {
				//if no entries then use equal bin widths
				yLimits[nDivx][nDivy] = ymin + (ymax-ymin)*nDivy/n_divy;
			}
		}
	}

	delete[] xIn;
	delete[] yIn;
	delete[] xIndex;
	delete[] yIndex;

	//call for each sub bin
	for (int nDivx = 0; nDivx < n_divx; ++nDivx){
		for (int nDivy = 0; nDivy < n_divy; ++nDivy){
			initHisto(xLimits[nDivx], xLimits[nDivx + 1], yLimits[nDivx][nDivy], yLimits[nDivx][nDivy + 1],iter+1);
		}
	}
} """
fig, ax = plt.subplots(1,1)
model_b_2dprojection_values = model_m_omega.pdf(model_omega_np).numpy() * (model_m_omega.get_yield().numpy())
model_omega_2dprojection_values = model_m_b.pdf(model_b_np).numpy() * (model_m_b.get_yield().numpy())

model_omega_np2d = np.linspace(*m_omega.limit1d, n_model_points_per_bin*n_bins_2d)
model_b_np2d = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins_2d)
X, Y = np.meshgrid(model_omega_np2d, model_b_np2d)
m_b_omega_2dmodel_array = np.array([(yy, xx) for xx in model_b_np2d for yy in model_omega_np2d])
model_m_b_omega_values_reshaped = (model_m_b_omega.pdf(m_b_omega_2dmodel_array).numpy()*model_m_b_omega.get_yield().numpy()*(m_b_omega.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_model_points_per_bin*n_bins_2d,n_model_points_per_bin*n_bins_2d)
#model_m_b_omega_values = (model_m_b_omega.pdf(m_b_omega_2dmodel_array).numpy()*model_m_b_omega.get_yield().numpy()*(m_b_omega.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_bins_2d, n_bins_2d)

datahist, xdatabins, ydatabins, im = ax.hist2d(data_m_b0["b_cons_omegaxpiz_m_best"].numpy(), data_m_omega["x_cons_xpiz_m_best"].numpy(), (n_bins_2d, n_bins_2d))
ax.contour(Y, X, model_m_b_omega_values_reshaped, levels=10, cmap="autumn_r", linestyles="solid")
fig.savefig("OutputFiles/PNGPlots/FitResults/2dfit.png")

fig, ax = plt.subplots(1,1)
datahist, xdatabins, ydatabins, dataim = ax.hist2d(data_m_b0["b_cons_omegaxpiz_m_best"].numpy(), data_m_omega["x_cons_xpiz_m_best"].numpy(), (n_bins_2d, n_bins_2d))
omega_bincenters = np.linspace(*m_omega.limit1d, n_bins_2d)
b_bincenters = np.linspace(*m_b.limit1d, n_bins_2d)
pull2dmatrix = np.zeros((len(xdatabins)-1, len(ydatabins)-1))
b_omega_bincenters = np.array([(yy, xx) for xx in b_bincenters for yy in omega_bincenters])
model_b_omega_bincenters = (model_m_b_omega.pdf(b_omega_bincenters).numpy()*model_m_b_omega.get_yield().numpy()*(m_b_omega.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_bins_2d, n_bins_2d)
pull2d = (datahist-model_b_omega_bincenters)/np.sqrt(datahist)
for i in range(len(ydatabins)-1):
    for j in range(len(xdatabins)-1):
        #if np.sqrt(datahist.T[i,j]) == 0:
        #    print(datahist.T[i,j], model_b_omega_bincenters[i+j])
        #pull2d = (datahist.T[i,j]-model_b_omega_bincenters[i+j*(len(ydatabins)-1)])/np.sqrt(datahist.T[i,j])
        pull2dmatrix[i, j] = pull2d[i, j]
        #print(xdatabins[j], ydatabins[i], pull2d[i, j])
        ax.text(xdatabins[j]+bins_2d_b_width/2,ydatabins[i]+bins_2d_omega_width/2, r'${0:.2f}$'.format(pull2d[i,j]), color="w", ha="center", va="center", fontweight="bold")
fig.savefig("OutputFiles/PNGPlots/FitResults/2DDatawithtextpulls.png")
fig = plt.figure()          #create a canvas, tell matplotlib it's 3d
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot()
cmap = mpl.colormaps.get_cmap('jet') # Get desired colormap - you can change this!
xpos, ypos = np.meshgrid(xdatabins[:-1]+xdatabins[1:], ydatabins[:-1]+ydatabins[1:])

xpos = xpos.flatten()/2.
ypos = ypos.flatten()/2.
zpos = pull2dmatrix[:-1, :-1]

dx = xdatabins [1] - xdatabins [0]
dy = ydatabins [1] - ydatabins [0]
dz = pull2dmatrix.flatten()

cmap = mpl.colormaps.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(dz)   # get range of colorbars so we can normalize
min_height = np.min(dz)
# scale each z to [0,1], and get their rgb values
rgba = [cmap((k-min_height)/max_height) for k in dz] 

#ax.imshow(xpos, ypos, zpos, xdatabins[1]-xdatabins[0], ydatabins[1]-ydatabins[0], dz, zsort='average')
c = plt.imshow(zpos, vmin = min_height, vmax = max_height,
                 extent =[xdatabins.min(), xdatabins.max(), ydatabins.min(), ydatabins.max()])
plt.colorbar(c)
plt.savefig("OutputFiles/PNGPlots/FitResults/2DPulls_RegularBins.png")
#plt.colorbar()
#plt.hist2d(model_b_omega_np_counts, m_b, m_omega)
#plot_scaling_peakX = n_model_points_per_bin*(frac_yx_YY*n_yx_YY+frac_2x_YY*n_2x_YY+frac_3x_YY*n_3x_YY)/model_m_b_omega.sum()
#plt.show()