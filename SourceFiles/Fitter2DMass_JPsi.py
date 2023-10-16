import mplhep
import numpy as np

import tensorflow as tf

import uproot
import pandas as pd

import zfit
from zfit import z
from zfit import core

mplhep.style.use("LHCb2")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import cm
import hist

import calendar
import time

# Write the structure for the table for latex
with open("OutputFiles/LatexTables/PreliminaryFit/MVAScan/latex_table_nevts_uncert_vs_mvacut.tex", "w") as f:
    print('\\begin{table}[!h]', file=f)
    print('\centering', file=f)
    print('\\begin{adjustbox}{width=18cm, rotate=90}', file=f)
    print('\\begin{tabular}{*{1}{c|}*{16}{c}}', file=f)
    print('& \multicolumn{16}{c}{MVA BDT Cut}\\\\', file=f)

mvabdt_cut_list = []
n_b_sigevents_list = []
sigma_n_b_sigevents_list = []
n_jpsi_sigevents_list = []
sigma_n_jpsi_sigevents_list = []

for mvabdt_cut in np.linspace(0.40, 0.55, num=16):
    # Spaces (normalization ranges)

    # Mass Spaces
    m_jpsi = zfit.Space("x_cons_xpiz_m_best",(3000, 3200), name=r"$m_{DTF}^{\pi^{0}}\left(J/\Psi\right)\left[MeV/c^{2}\right]$")
    m_b = zfit.Space("b_cons_Jpsixpiz_m_best", (5280-150, 5280+150), name=r"$m_{DTF}^{\pi^{0},J/\Psi}\left(B^{0}\right)\left[MeV/c^{2}\right]$")
    # 2D Space for the 2D fit
    m_b_jpsi = m_b*m_jpsi

    n_bins = 50     #bestvalue 30
    n_bins_2d = 10  #bestvalue 10
    bins_2d_b_width = m_b.area()/n_bins_2d
    bins_2d_jpsi_width = m_jpsi.area()/n_bins_2d
    n_model_points_per_bin = 40

    # Data (data used as input and relative filters)

    # Dataset
    file = uproot.open("InputFiles/DATA_SignalRegion_JpsiMVA_BDT_20072023_OmegaCUT.root")
    #file = uproot.open("InputFiles/JpsiMCStrippingFiltered.root")
    tree = file["DecayTree"]
    varstobeloaded = ["b_cons_Jpsixpiz_m_best", "x_cons_xpiz_m_best"]
    # Filters
    filter  = f"(mvacut_kbdt1>{mvabdt_cut})"
    filter += "& ((s_kpi<1.844*1.844*1e6) | (s_kpi>1.884*1.884*1e6))"
    filter += "& ((s_kxpi<1.844*1.844*1e6) | (s_kxpi>1.884*1.884*1e6))"
    filter += "& ((s_kpi>0.6**2*1e6) & (s_kpi<1**2*1e6))"
    filter += "& (Pi0Merged == 0)"
    filter += f"& (x_cons_xpiz_m_best<{float(m_jpsi.upper[0])}) & (x_cons_xpiz_m_best>{float(m_jpsi.lower[0])})"
    filter += f"& (b_cons_Jpsixpiz_m_best<{float(m_b.upper[0])}) & (b_cons_Jpsixpiz_m_best>{float(m_b.lower[0])})"
    # Building the dataframe
    dataframe = tree.arrays(varstobeloaded, filter,  library='pd')
    # Building the variable's arrays
    data_m_jpsi = zfit.Data.from_pandas(df=dataframe, obs=m_jpsi)
    data_m_b0 = zfit.Data.from_pandas(df=dataframe, obs=m_b)
    data_m_b_jpsi = zfit.Data.from_pandas(df=dataframe, obs=m_b_jpsi)

    # Parameters

    # Defining a dict for the composite params
    def param_product(params):
        return params["a"] * params["b"]

    def param_sum(params):
        return params["a"] + params["b"]

    # Shared signal parameters
    zfit_parameters = []
    yield_b_jpsi = zfit.param.Parameter(r"$n^{B^0-jpsi}$"+f"{mvabdt_cut}", 35000, 0, 1e5)
    yield_bjpsi = zfit.param.Parameter(r"$n^{B^0-J/\psi}$"+f"{mvabdt_cut}", 1000, 0,)
    yield_b = zfit.param.Parameter(r"$yield^{B^0}$"+f"{mvabdt_cut}", 80000, 0, )
    
    # J/Psi signal parameters

    # First Crystall Ball parameters
    mu_jpsi = zfit.Parameter(r"$\mu^{J/\psi}$"+f"{mvabdt_cut}", 3091.7, float(m_jpsi.lower[0]), float(m_jpsi.upper[0]))
    sigma_jpsi = zfit.Parameter(r"$\sigma^{J/\psi}$"+f"{mvabdt_cut}", 46., 0., 110.)
    alpha_jpsi_cb1 = zfit.param.ConstantParameter(r"$\alpha^{J/\psi}$"+f"{mvabdt_cut}", 2.02201)
    n_jpsi_cb1 = zfit.param.ConstantParameter(r"$n^{J/\psi}$"+f"{mvabdt_cut}", 3)
    # Second Crystall Ball parameters
    sigmaratio_jpsi = zfit.param.ConstantParameter(r"sigma2/sigma1"+f"{mvabdt_cut}", 0.362191)
    sigma_jpsi_persigmaratio = zfit.ComposedParameter(r"$\sigma\times\sigma_{Ratio}"+f"{mvabdt_cut}", param_product, params={"a": sigma_jpsi, "b": sigmaratio_jpsi})
    alpha_jpsi_cb2 = zfit.param.ConstantParameter(r"$\alpha^{J/\psi}_{CB2}$"+f"{mvabdt_cut}", -0.954915)
    n_jpsi_cb2 = zfit.param.ConstantParameter(r"$n^{J/\psi}_{CB2}$"+f"{mvabdt_cut}", 3)

    # J/Psi background parameters

    # Combinatorial background with Chebyshev polynomials
    coeff_one_cheby = zfit.Parameter(r"$a_{1,Cheby}^{J/\Psi}"+f"{mvabdt_cut}", -0.4, -100, 100)
    coeff_two_cheby = zfit.Parameter(r"$a_{2,Cheby}^{J/\Psi}"+f"{mvabdt_cut}", -0.05, -100, 0)
    coeff_three_cheby = zfit.Parameter(r"$a_{3,Cheby}^{J/\Psi}"+f"{mvabdt_cut}", -0.5, -100, 100)
    # Combinatorial background with exponential
    lambda_jpsi_comb = zfit.Parameter(r"$\lambda^{\jpsi}_{COMB}$"+f"{mvabdt_cut}", -0.001, -10., 0.)

    # Yields
    n_b_events = zfit.Parameter(r"$n_{B^{0}}^{Tot}$"+f"{mvabdt_cut}", 1e4, 0, 5e6)
    n_jpsi_events = zfit.Parameter(r"$n_{J/\psi}^{Tot}$"+f"{mvabdt_cut}", 2e4, 0, 5e6)
    n2_jpsi_sig = zfit.Parameter(r"$n_{J/\psi,2}^{Sig}$"+f"{mvabdt_cut}", 2e4, 0, 5e6)
    n_sig_jpsi_cb1 = zfit.Parameter(r"$n^{\jpsi}_{Sig-CB1}$"+f"{mvabdt_cut}", 3735, 0,)
    n_sig_jpsi_cb2 = zfit.Parameter(r"$n^{\jpsi}_{Sig-CB2}$"+f"{mvabdt_cut}", 708, 0,)
    n_bkg_jpsi_comb = zfit.Parameter(r"$n^{\jpsi}_{Bkg-Comb}$"+f"{mvabdt_cut}", 45059, 0,)
    n_bkg_b_comb = zfit.Parameter(r"$n^{B}_{Bkg-Comb}$"+f"{mvabdt_cut}", 45697, 0, 1e7)

    #B signal parameters
    mu_b = zfit.Parameter(r"$\mu^{B^0}$"+f"{mvabdt_cut}", 5279.68, 5279-100., 5279+100.)
    sigma_b = zfit.Parameter(r"$\sigma^{B^0}$"+f"{mvabdt_cut}", 11.43, 0., 80.)
    alpha_b = zfit.param.ConstantParameter(r"$\alpha^{B^0}$"+f"{mvabdt_cut}", -1.86)
    n_b = zfit.param.ConstantParameter(r"$n^{B^0}$"+f"{mvabdt_cut}", 3)
    n_b_sig = zfit.Parameter(r"$n_{B^0,1}^{Sig}$"+f"{mvabdt_cut}", 3735, 0, 1e6)

    sigmaratio_b = zfit.param.ConstantParameter(r"$\sigma_{Ratio}^{B^0}$"+f"{mvabdt_cut}", 0.432)
    alpha2_b = zfit.param.ConstantParameter(r"$\alpha_2^{B^0}$"+f"{mvabdt_cut}", 1.64)
    n2_b = zfit.param.ConstantParameter(r"$n_2^{B^0}$"+f"{mvabdt_cut}", 3)
    n2_b_sig = zfit.Parameter(r"$n_{B^0,2}^{Sig}}$"+f"{mvabdt_cut}", 70, 0, 1e6)

    sigmaratio_jpsi_comp = zfit.ComposedParameter('sigmaratio_jpsi_comp'+f"{mvabdt_cut}", param_product, params={"a": sigmaratio_jpsi, "b": sigma_jpsi})
    sigmaratio_b_comp = zfit.ComposedParameter('sigmaratio_b_comp'+f"{mvabdt_cut}", param_product, params={"a": sigmaratio_b, "b": sigma_b})

    # B background parameters
    coeff_m_b_one_cheby = zfit.Parameter(r"$a_{1,Cheby}^{B^{0}}"+f"{mvabdt_cut}", -0.742762, -100, 100)
    coeff_m_b_two_cheby = zfit.Parameter(r"$a_{2,Cheby}^{B^{0}}"+f"{mvabdt_cut}", -0.183, -100, 100)
    coeff_m_b_three_cheby = zfit.Parameter(r"$a_{3,Cheby}^{B^{0}}"+f"{mvabdt_cut}", -0.5, -100, 100)
    #########################################
    #                 Models                #
    #########################################

    # jpsi mass

    # Signal - Sum of two CBs 
    model_m_jpsi_cb1 = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigma_jpsi, alpha=alpha_jpsi_cb1, n=n_jpsi_cb1)
    model_m_jpsi_cb1_extended = model_m_jpsi_cb1.create_extended(n_sig_jpsi_cb1, name='model_m_jpsi_cb1')
    model_m_jpsi_cb2 = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigma_jpsi_persigmaratio, alpha=alpha_jpsi_cb2, n=n_jpsi_cb2)
    model_m_jpsi_cb2_extended = model_m_jpsi_cb2.create_extended(n_sig_jpsi_cb2, name='model_m_jpsi_cb2')
    jpsisigfrac = zfit.Parameter('jpsisigfrac'+f"{mvabdt_cut}", 0.1, 0., 100)
    model_m_jpsi_sig = zfit.pdf.SumPDF([model_m_jpsi_cb1, model_m_jpsi_cb2], jpsisigfrac)
    model_m_jpsi_sig_ext = zfit.pdf.SumPDF([model_m_jpsi_cb1_extended, model_m_jpsi_cb2_extended])

    # Background - Combinatorial
    model_m_jpsi_cheby = zfit.pdf.Chebyshev(obs=m_jpsi, coeffs=[coeff_one_cheby, coeff_two_cheby])
    model_m_jpsi_cheby_extended = model_m_jpsi_cheby.create_extended(n_bkg_jpsi_comb, name='model_m_jpsi_chebycomb')
    model_m_jpsi_comb = zfit.pdf.Exponential(obs=m_jpsi, lam=lambda_jpsi_comb)
    model_m_jpsi_comb_extended = model_m_jpsi_comb.create_extended(n_bkg_jpsi_comb, name='model_m_jpsi_comb')
    #model_m_jpsi_bkg = zfit.pdf.SumPDF([model_m_jpsi_comb, model_m_jpsi_preco, model_m_jpsi_precohm], [jpsibkgfrac, jpsibkgfrac_preco]) # With pipiw background in high m_jpsi region
    model_m_jpsi_bkg_ext = zfit.pdf.SumPDF([model_m_jpsi_sig_ext, model_m_jpsi_cheby_extended])
    #model_m_jpsi_bkg_ext = zfit.pdf.SumPDF([model_m_jpsi_comb_extended, model_m_jpsi_preco_extended, model_m_jpsi_precohm_extended])  # With pipiw background in high m_jpsi region

    # jpsi mass total fit function
    jpsifrac = zfit.Parameter('jpsifrac'+f"{mvabdt_cut}", 0.1, 0., 100)
    model_m_jpsi = zfit.pdf.SumPDF([model_m_jpsi_sig, model_m_jpsi_cheby], [jpsifrac])
    model_m_jpsi = model_m_jpsi.create_extended(n_jpsi_events, name="n_jpsi_sig")

    # B mass

    # B mass fit model
    model_m_b0_cb1 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigma_b, alpha=alpha_b, n=n_b)
    model_m_b0_cb1_extended = model_m_b0_cb1.create_extended(n_b_sig, name='B0_CrystalBall1B')
    model_m_b0_cb2 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigmaratio_b_comp, alpha=alpha2_b, n=n2_b)
    model_m_b0_cb2_extended = model_m_b0_cb2.create_extended(n2_b_sig, name='B0_CrystalBall2B')

    b_sigfrac = zfit.Parameter('bsigfrac'+f"{mvabdt_cut}", 0.1, 0., 100)
    model_m_b0_sig = zfit.pdf.SumPDF([model_m_b0_cb1, model_m_b0_cb2], b_sigfrac)
    model_m_b0_sig_ext = zfit.pdf.SumPDF([model_m_b0_cb1_extended, model_m_b0_cb2_extended])

    # B mass combinatorial
    model_m_b_cheby = zfit.pdf.Chebyshev(obs=m_b, coeffs=[coeff_m_b_one_cheby, coeff_m_b_two_cheby])
    model_m_b_cheby_extended = model_m_b_cheby.create_extended(n_bkg_b_comb, name='model_m_b_chebycomb')
    lambda_b_m_comb = zfit.Parameter(r"$\lambda^{B^0}_{Comb}$"+f"{mvabdt_cut}", -0.0008, -1100, 0)
    n_b_m_comb = zfit.Parameter(r"$n^{B^0}_{Comb}$"+f"{mvabdt_cut}", 1.4e4, 0, 1e7)
    model_b_m_comb = zfit.pdf.Exponential(obs=m_b, lam=lambda_b_m_comb)
    model_b_m_comb_ext = model_b_m_comb.create_extended(n_b_m_comb, name='b_m_Combbkg')
    frac_b_m_sig = zfit.Parameter('frac_b_m_sig'+f"{mvabdt_cut}", 0.2, 0., 0.999)
    frac_b_m_comb = zfit.Parameter('frac_b_m_comb'+f"{mvabdt_cut}", 0.2, 0., 0.999)

    # B mass partially reconstructed
    kdegridpoints = 32
    """
    filepartreco = uproot.open("InputFiles/B2chic02Jpsi2pippimpi0_gamma_Kpi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["M_Kpi3pi_miss_chic0gamma"], library='pd')
    dataframepartreco["b_cons_Jpsixpiz_m_best"] = 1000*treepartreco.arrays(["M_Kpi3pi_miss_chic0gamma"], library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    kde_K1jpsi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)

    filepartreco = uproot.open("InputFiles/B2chic02pippimpipi_Kpi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["M_Kpi3pi_miss_chic0pion0"], library='pd')
    dataframepartreco["b_cons_Jpsixpiz_m_best"] = 1000*treepartreco.arrays(["M_Kpi3pi_miss_chic0pion0"], library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    kde_chic0_4pi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)

    filepartreco = uproot.open("InputFiles/B2chic12Jpsi2pippimpi0_gamma_Kpi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["M_Kpi3pi_miss_chic1gamma"], library='pd')
    dataframepartreco["b_cons_Jpsixpiz_m_best"] = 1000*treepartreco.arrays(["M_Kpi3pi_miss_chic1gamma"], library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    kde_chic1_jpsi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)

    filepartreco = uproot.open("InputFiles/B2chic12pippimpipi_Kpi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["M_Kpi3pi_miss_chic1pion0"], library='pd')
    dataframepartreco["b_cons_Jpsixpiz_m_best"] = 1000*treepartreco.arrays(["M_Kpi3pi_miss_chic1pion0"], library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    kde_chic1_4pi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
    """
    """
    filepartreco = uproot.open("InputFiles/B2Dstarjpsi2Dpiz3pi2Kpi3pi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["b_cons_Jpsixpiz_m_best"], recofilter, library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    yield_bpartreco_dstarjpsi = zfit.param.Parameter(r"$yield^{partreco}_{D^{*}_{\pi^{0}}\jpsi}$"+f"{mvabdt_cut}", 2500, 0, )
    kde_dstarjpsi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_dstarjpsi, num_grid_points=kdegridpoints)

    filepartreco = uproot.open("InputFiles/B2Dstarjpsi2Dgamma3pi2Kpi3pi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["b_cons_Jpsixpiz_m_best"], recofilter, library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    yield_bpartreco_dstargammajpsi = zfit.param.Parameter(r"$yield^{partreco}_{D^{*}_{\gamma}\jpsi}$"+f"{mvabdt_cut}", 2500, 0, )
    kde_dstargammajpsi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
    kde_dstargammajpsi_ext = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_dstargammajpsi, num_grid_points=kdegridpoints)

    filepartreco = uproot.open("InputFiles/B2Djpsi_tree.root")
    treepartreco = filepartreco["DecayTree"]
    dataframepartreco = treepartreco.arrays(["b_cons_Jpsixpiz_m_best"], recofilter, library='pd')
    data_rapidsim_b = zfit.Data.from_pandas(df=dataframepartreco, obs=m_b)
    yield_bpartreco_djpsi = zfit.param.Parameter(r"$yield^{partreco}_{D\jpsi}$"+f"{mvabdt_cut}", 2500, 0, )
    kde_djpsi = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, num_grid_points=kdegridpoints)
    kde_djpsi_ext = zfit.pdf.KDE1DimGrid(data=data_rapidsim_b, obs=m_b, padding=0.05, bandwidth=500/20, extended=yield_bpartreco_djpsi, num_grid_points=kdegridpoints)
    """
    frac_kde = zfit.Parameter('frac_kde'+f"{mvabdt_cut}", 0.05, 0., 0.999)
    frac_kde2 = zfit.Parameter('frac_kde2'+f"{mvabdt_cut}", 0.8, 0., 0.999)
    frac_kde3 = zfit.Parameter('frac_kde3'+f"{mvabdt_cut}", 0.05, 0., 0.999)
    yield_kde = zfit.Parameter('yield_kde'+f"{mvabdt_cut}", 100, 0., 1e5)
    #kde = zfit.pdf.SumPDF([kde_K1jpsi, kde_chic0_4pi, kde_chic1_jpsi, kde_chic1_4pi], [frac_kde, frac_kde2, frac_kde3])
    #kde_ext = kde.create_extended(yield_kde, 'kde_ext')

    #mu_partiallyreco = zfit.Parameter(r"$\mu_{part}$"+f"{mvabdt_cut}", 5063, 5000, 5100)
    #sigma_partiallyreco = zfit.Parameter(r"$\sigma_{part}$"+f"{mvabdt_cut}", 95.5, 0, 300)
    #n_b_m_partiallyreco = zfit.Parameter(r"$n^{B^0}_{Part}$"+f"{mvabdt_cut}", 2.8e4, 0, 1e7)
    #model_b_m_partiallyreco = zfit.pdf.Gauss(mu=mu_partiallyreco, sigma=sigma_partiallyreco, obs=m_b)
    #model_b_m_partiallyreco_ext = model_b_m_partiallyreco.create_extended(n_b_m_partiallyreco, name='b_m_PartReco')


    # Mass diagonal component
    mu_diag = zfit.Parameter(r"$\mu_{diag}$"+f"{mvabdt_cut}", 3000, 600, 5500)
    sigma_diag = zfit.Parameter(r"$\sigma_{diag}$"+f"{mvabdt_cut}", 95.5, 0, 1000)
    n_b_m_diag = zfit.Parameter(r"$n^{B^0}_{diag}$"+f"{mvabdt_cut}", 2.8e4, 0, 1e7)
    model_b_m_diag = zfit.pdf.Gauss(mu=mu_diag, sigma=sigma_diag, obs='b_cons_Jpsixpiz_m_best-x_cons_xpiz_m_best')
    model_b_m_diag_ext = model_b_m_diag.create_extended(n_b_m_diag, name='b_m_diag')

    # 1D mass diagonal component
    class DiagonalGaussianPDF(zfit.pdf.ZPDF):
        """ Diagonal component for the mass of the B and the jpsi

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

            bmass, jpsimass = z.unstack_x(x)

            pdf = (1 / (tf.sqrt(diagsigma))) * tf.exp( - ( tf.square( (bmass - jpsimass)-diagmean ) / (2 * tf.square(diagsigma) ) ) )
            #pdf = bmass - jpsimass

            return pdf
        
    mass_diaggauss = zfit.Parameter(r"$m_{B^0}-m_{#jpsi}^{Diag}$"+f"{mvabdt_cut}", 3993, 0, 5000)
    sigma_diaggauss = zfit.Parameter(r"$#sigma_{DiagGauss}$"+f"{mvabdt_cut}", 99, 0.1, 200.)
    diag_gauss_pdf = DiagonalGaussianPDF(obs=m_b_jpsi, sigma=sigma_diaggauss, mean=mass_diaggauss)

    #1D model
    b_frac = zfit.Parameter('bfrac'+f"{mvabdt_cut}", 0.1, 0., 100)
    model_m_b = zfit.pdf.SumPDF([model_m_b0_sig, model_m_b_cheby], b_frac)
    model_m_b = model_m_b.create_extended(n_b_events)

    #2D model
    #model_m_b_jpsi_sig = zfit.pdf.ProductPDF([model_m_b0_sig_ext, model_m_jpsi_sig_ext])
    #model_m_b_jpsi_sig = model_m_b_jpsi_sig.create_extended(yield_b_jpsi)
    #model_m_b_jpsi_bkg = zfit.pdf.ProductPDF([model_m_jpsi_comb_extended, model_b_m_comb_ext, model_b_m_partiallyreco])
    #model_m_b_jpsi_bkg = model_m_b_jpsi_bkg.create_extended(yield_bjpsi)
    n_m_b_jpsi_sig = zfit.Parameter(r"$n^{B^0-jpsi}_{2DSIG}$"+f"{mvabdt_cut}", 0.24, 1e-4, 1)
    n_m_b_jpsi_comb = zfit.Parameter(r"$n^{B^0-jpsi}_{2DCOMB}$"+f"{mvabdt_cut}", 0.026, 1e-4, 1)
    n_m_breco_jpsisig = zfit.Parameter(r"$n^{B^0-jpsi}_{2DRECO-SIG}$"+f"{mvabdt_cut}", 0.47, 1e-4, 1)
    n_m_breco_jpsicomb = zfit.Parameter(r"$n^{B^0-jpsi}_{2DRECO-COMB}$"+f"{mvabdt_cut}", 0.1, 1e-4, 1)
    n_m_bsig_jpsicomb = zfit.Parameter(r"$n^{B^0-jpsi}_{2DSIG-COMB}$"+f"{mvabdt_cut}", 0.01, 1e-4, 1)
    n_m_bcomb_jpsisig = zfit.Parameter(r"$n^{B^0-jpsi}_{2DCOMB-SIG}$"+f"{mvabdt_cut}", 0.1, 1e-4, 1)
    n_m_bcomb_diag = zfit.Parameter(r"$n^{B^0-jpsi}_{Diag}$"+f"{mvabdt_cut}", 0.1, 1e-4, 1)
    n_b_comb_per_jpsi_phi = zfit.Parameter(r"$n^{B^0-\phi}_{COMB-SIG}$"+f"{mvabdt_cut}", 0.1, 1e-4, 1)

    b_sig_per_jpsi_sig = zfit.pdf.ProductPDF([model_m_jpsi_sig, model_m_b0_sig], name = r"$B_{Sig}\cdot J\psi_{Sig}$"+f"{mvabdt_cut}")
    b_comb_per_jpsi_comb = zfit.pdf.ProductPDF([model_m_jpsi_cheby, model_m_b_cheby], name = r"$B_{Comb}\cdot J\psi_{Comb}$"+f"{mvabdt_cut}")
    b_sig_per_jpsi_comb = zfit.pdf.ProductPDF([model_m_jpsi_cheby, model_m_b0_sig], name = r"$B_{Sig}\cdot J\psi_{Comb}$"+f"{mvabdt_cut}")
    #b_tot_per_jpsi_phi = zfit.pdf.ProductPDF([model_m_phi_sig, diag_gauss_pdf], name = r"$B_{Diag}\cdot\phi_{Sig}$"+f"{mvabdt_cut}")
    b_comb_per_jpsi_sig = zfit.pdf.ProductPDF([model_m_jpsi_sig, model_m_b_cheby], name = r"$B_{Comb}\cdot J\psi_{Sig}$"+f"{mvabdt_cut}")
    #b_comb_per_jpsi_phi = zfit.pdf.ProductPDF([model_m_phi_sig, model_b_m_comb], name = r"$B_{Comb}\cdot\phi{Sig}$"+f"{mvabdt_cut}")
    model_m_b_jpsi = zfit.pdf.SumPDF([b_sig_per_jpsi_sig,
                                    b_sig_per_jpsi_comb,
                                    b_comb_per_jpsi_comb,
                                    b_comb_per_jpsi_sig,
                                    #b_preco_per_jpsi_comb,
                                    #b_tot_per_jpsi_phi
                                    ],
                                    [n_m_b_jpsi_sig, n_m_b_jpsi_comb, n_m_breco_jpsisig, n_m_bcomb_diag])
    model_m_b_jpsi = model_m_b_jpsi.create_extended(yield_bjpsi)

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
    nll_m_jpsi = zfit.loss.ExtendedUnbinnedNLL(model=model_m_jpsi, data=data_m_jpsi)
    nll_m_b = zfit.loss.ExtendedUnbinnedNLL(model=model_m_b, data=data_m_b0)

    # Stage 2: instantiate a minimiser (in this case a basic minuit)
    minimizer = zfit.minimize.Minuit()

    # Stage 3: minimise the given negative log-likelihood
    results_jpsi = minimizer.minimize(nll_m_jpsi)
    current_GMT = time.gmtime()
    time_stamp = calendar.timegm(current_GMT)
    print("Current timestamp - Before Minuit:", time_stamp)
    results_b = minimizer.minimize(nll_m_b)

    current_GMT = time.gmtime()
    time_stamp = calendar.timegm(current_GMT)
    print("Current timestamp - After Minuit:", time_stamp)
    # Stage 4: computing the errors
    params_errors_jpsi = results_jpsi.hesse()
    current_GMT = time.gmtime()
    time_stamp = calendar.timegm(current_GMT)
    print("Current timestamp - Before Hesse:", time_stamp)
    params_errors_b = results_b.hesse()
    current_GMT = time.gmtime()
    time_stamp = calendar.timegm(current_GMT)
    print("Current timestamp - After Hesse:", time_stamp)

    # Saving Pars Infos
    paramsjpsi = results_jpsi.params
    paramsb = results_b.params

    # Printing the results
    print("Function minimum:", results_jpsi.fmin)
    print("Converged:", results_jpsi.converged)
    print("Full minimizer information:", results_jpsi)

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
        area2d = m_b.area()*m_jpsi.area()
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
            pdfs = [#(m.pdf(basis)*m.get_).numpy()*area/nbins
                    zfit.run(m.pdf(basis, norm_range=obs)*zmodel.get_yield()*frac)
                    for m, frac in zip(zmodel.pdfs, zmodel.params.values())]
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
                parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['value'])
            else:
                parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'], params[par]['value'])
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

    data_m_jpsi_np = np.sort(data_m_jpsi["x_cons_xpiz_m_best"].numpy())
    data_m_b_np = np.sort(data_m_b0["b_cons_Jpsixpiz_m_best"].numpy())
    #data_mc_jpsi_np = np.sort(data_mc_jpsi["b_cons_Jpsixpiz_m_best"].numpy())
    m_jpsi_area = m_jpsi.area().numpy()
    m_b_area = m_b.area().numpy()
    #model_x_np = np.linspace(*m_jpsi.limit1d, n_model_points_per_bin*n_bins)
    #model_np = model_m_b.pdf(model_x_np).numpy() * (model_m_b.get_yield().numpy())
    #fig = plot_fit(data_mc_jpsi_np, model_x_np, model_np, [m_b, m_jpsi], params, nbins=n_bins, zmodel=model_m_b)
    #fig.savefig('test1.pdf')
    #plt.show()

    lower, upper = m_jpsi.limits
    model_jpsi_np = np.linspace(lower[-1][0], upper[0][0], num=1000)
    #model_jpsi_np = np.linspace(*m_jpsi.limit1d, n_model_points_per_bin*n_bins)
    model_np = model_m_jpsi.pdf(model_jpsi_np).numpy() * (model_m_jpsi.get_yield().numpy())
    fig_jpsi = plot_fit(data_m_jpsi_np, model_jpsi_np, model_np, m_jpsi, paramsjpsi, nbins=n_bins, zmodel=model_m_jpsi, parlist=False)
    fig_jpsi.savefig('jpsi1dfit.pdf')
    fig_jpsi.savefig('jpsi1dfit.png')
    fig_jpsi.savefig("OutputFiles/PNGPlots/FitResults/jpsi1dfit.png")
    fig_jpsi.set_size_inches(9,9)
    fig_jpsi.set_dpi(120)
    fig_jpsi.savefig(f"OutputFiles/PNGPlots/PreliminaryFit/MVAScan/jpsi1dfit_mva{mvabdt_cut:.2f}_jpsichannel.png")

    model_b_np = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins)
    model_np = model_m_b.pdf(model_b_np).numpy() * (model_m_b.get_yield().numpy())
    fig_b = plot_fit(data_m_b_np, model_b_np, model_np, m_b, paramsb, nbins=n_bins, zmodel=model_m_b, parlist=False)
    fig_b.savefig("OutputFiles/PNGPlots/FitResults/b1dfit_jpsichannel.png")
    fig_b.set_size_inches(9,9)
    fig_b.set_dpi(120)
    fig_b.savefig(f"OutputFiles/PNGPlots/PreliminaryFit/MVAScan/b1dfit_mva{mvabdt_cut:.2f}_jpsichannel.png")

    mvabdt_cut_list.append(mvabdt_cut)
    n_b_sigevents_list.append(n_b_events.value()*b_frac.value())
    sigma_n_b_sigevents_list.append(n_b_events.value()*b_frac.value()*(paramsb[n_b_events]['hesse']['error']/paramsb[n_b_events]['value']+paramsb[b_frac]['hesse']['error']/paramsb[b_frac]['value']))
    n_jpsi_sigevents_list.append(n_jpsi_events.value()*jpsifrac.value())
    sigma_n_jpsi_sigevents_list.append(n_jpsi_events.value()*jpsifrac.value()*(paramsjpsi[n_jpsi_events]['hesse']['error']/paramsjpsi[n_jpsi_events]['value']+paramsjpsi[jpsifrac]['hesse']['error']/paramsjpsi[jpsifrac]['value']))

with open("AnalysisNote/LatexTables/PreliminaryFit/MVAScan/latex_table_nevts_uncert_vs_mvacut.tex", "a") as f:
    for i, mvaval in enumerate(mvabdt_cut_list):
        print(f"& {mvaval:.3}", end='', file=f)
    print('\\\\ \\hline', file=f)
    print('$N_B$', end='', file=f)
    [print(f"& {nbevt:.3}", end='', file=f) for nbevt in n_b_sigevents_list]
    print('\\\\', file=f)
    print('$\\sigma_B$', end='', file=f)
    [print(f"& {sigmanb:.3}", end='', file=f) for sigmanb in sigma_n_b_sigevents_list]
    print('\\\\', file=f)
    print('$N_{J/\\psi}$', end='', file=f)
    [print(f"& {njpsievt:.3}", end='', file=f) for njpsievt in n_jpsi_sigevents_list]
    print('\\\\', file=f)
    print('$\\sigma_{J/\\psi}$', end='', file=f)
    [print(f"& {sigmanb:.3}", end='', file=f) for sigmanb in sigma_n_jpsi_sigevents_list]
    print('\\\\ \\hline', file=f)
    print('$FOM_{B}:=\\frac{N_B}{\\sigma_B}$', end='', file=f)
    for i in range(len(mvabdt_cut_list)):
        result=n_b_sigevents_list[i]/sigma_n_b_sigevents_list[i]
        print(f"& {result:.3}", end='', file=f) 
    print('\\\\', file=f)
    print('$FOM_{J/\\psi}:=\\frac{N_{J/\\psi}}{\\sigma_{J/\\psi}}$', end='', file=f)
    for i in range(len(mvabdt_cut_list)):
        result=n_jpsi_sigevents_list[i]/sigma_n_jpsi_sigevents_list[i]
        print(f"& {result:.3}", end='', file=f) 
    print('\\\\ \\hline', file=f)
    print('\\end{tabular}', file=f)
    print('\\end{adjustbox}', file=f)
    print('\\end{table}', file=f)
exit()

b_jpsi_constraints = []
#b_jpsi_constraints.append(zfit.constraint.GaussianConstraint(params=frac_kde, observation=frac_kde.value(), uncertainty=frac_kde.value().numpy()*1e-2))
#1b_jpsi_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_jpsi, observation=sigma_jpsi.value(), uncertainty=sigma_jpsi.value().numpy()*1e-2))
#1b_jpsi_constraints.append(zfit.constraint.GaussianConstraint(params=sigma_b, observation=sigma_b.value(), uncertainty=sigma_b.value().numpy()*1e-2))
nll_m_b_jpsi = zfit.loss.ExtendedUnbinnedNLL(model=model_m_b_jpsi, data=data_m_b_jpsi, constraints=b_jpsi_constraints)
results_b_jpsi = minimizer.minimize(nll_m_b_jpsi)
paramsbjpsi = results_b_jpsi.params
params_errors_b_jpsi = results_b_jpsi.hesse()

print("Function minimum:", results_b_jpsi.fmin)
print("Converged:", results_b_jpsi.converged)
print("Full minimizer information:", results_b_jpsi)

model_jpsi_m_projection = model_m_b_jpsi.create_projection_pdf(limits=m_b)
model_b_jpsi_np_jpsicounts = model_jpsi_m_projection.pdf(model_jpsi_np).numpy()*(yield_bjpsi)
fig_jpsi_projection = plot_projection_fit2d(data_m_jpsi_np, model_jpsi_np, model_b_jpsi_np_jpsicounts, m_jpsi, m_b, paramsbjpsi, nbins=n_bins, zmodel=model_m_b_jpsi, parlist=False)
fig_jpsi_projection.savefig('jpsi_2dfit_projection.pdf')
fig_jpsi_projection.savefig('OutputFiles/PNGPlots/FitResults/jpsi_2dfit_projection.png')

model_b_m_projection = model_m_b_jpsi.create_projection_pdf(limits=m_jpsi)
model_b_jpsi_np_bcounts = model_b_m_projection.pdf(model_b_np).numpy()*(yield_bjpsi)
fig_b_projection = plot_projection_fit2d(data_m_b_np, model_b_np, model_b_jpsi_np_bcounts, m_b, m_jpsi, paramsbjpsi, nbins=n_bins, zmodel=model_m_b_jpsi, parlist=False)
fig_b_projection.savefig('b_2dfit_projection.pdf')
fig_b_projection.savefig('OutputFiles/PNGPlots/FitResults/b_2djpsifit_projection.png')

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
model_b_2dprojection_values = model_m_jpsi.pdf(model_jpsi_np).numpy() * (model_m_jpsi.get_yield().numpy())
model_jpsi_2dprojection_values = model_m_b.pdf(model_b_np).numpy() * (model_m_b.get_yield().numpy())

model_jpsi_np2d = np.linspace(*m_jpsi.limit1d, n_model_points_per_bin*n_bins_2d)
model_b_np2d = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins_2d)
X, Y = np.meshgrid(model_jpsi_np2d, model_b_np2d)
m_b_jpsi_2dmodel_array = np.array([(yy, xx) for xx in model_b_np2d for yy in model_jpsi_np2d])
model_m_b_jpsi_values_reshaped = (model_m_b_jpsi.pdf(m_b_jpsi_2dmodel_array).numpy()*model_m_b_jpsi.get_yield().numpy()*(m_b_jpsi.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_model_points_per_bin*n_bins_2d,n_model_points_per_bin*n_bins_2d)
#model_m_b_jpsi_values = (model_m_b_jpsi.pdf(m_b_jpsi_2dmodel_array).numpy()*model_m_b_jpsi.get_yield().numpy()*(m_b_jpsi.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_bins_2d, n_bins_2d)

datahist, xdatabins, ydatabins, im = ax.hist2d(data_m_b0["b_cons_Jpsixpiz_m_best"].numpy(), data_m_jpsi["x_cons_xpiz_m_best"].numpy(), (n_bins_2d, n_bins_2d))
ax.contour(Y, X, model_m_b_jpsi_values_reshaped, levels=10, cmap="autumn_r", linestyles="solid")
fig.savefig("OutputFiles/PNGPlots/FitResults/2dfit_jpsi.png")

fig, ax = plt.subplots(1,1)
datahist, xdatabins, ydatabins, dataim = ax.hist2d(data_m_b0["b_cons_Jpsixpiz_m_best"].numpy(), data_m_jpsi["x_cons_xpiz_m_best"].numpy(), (n_bins_2d, n_bins_2d))
jpsi_bincenters = np.linspace(*m_jpsi.limit1d, n_bins_2d)
b_bincenters = np.linspace(*m_b.limit1d, n_bins_2d)
pull2dmatrix = np.zeros((len(xdatabins)-1, len(ydatabins)-1))
b_jpsi_bincenters = np.array([(yy, xx) for xx in b_bincenters for yy in jpsi_bincenters])
model_b_jpsi_bincenters = (model_m_b_jpsi.pdf(b_jpsi_bincenters).numpy()*model_m_b_jpsi.get_yield().numpy()*(m_b_jpsi.area().numpy())/n_bins_2d/n_bins_2d).reshape(n_bins_2d, n_bins_2d)
pull2d = (datahist-model_b_jpsi_bincenters)/np.sqrt(datahist)
for i in range(len(ydatabins)-1):
    for j in range(len(xdatabins)-1):
        #if np.sqrt(datahist.T[i,j]) == 0:
        #    print(datahist.T[i,j], model_b_jpsi_bincenters[i+j])
        #pull2d = (datahist.T[i,j]-model_b_jpsi_bincenters[i+j*(len(ydatabins)-1)])/np.sqrt(datahist.T[i,j])
        pull2dmatrix[i, j] = pull2d[i, j]
        #print(xdatabins[j], ydatabins[i], pull2d[i, j])
        ax.text(xdatabins[j]+bins_2d_b_width/2,ydatabins[i]+bins_2d_jpsi_width/2, r'${0:.2f}$'.format(pull2d[i,j]), color="w", ha="center", va="center", fontweight="bold")
fig.savefig("OutputFiles/PNGPlots/FitResults/2DDatawithtextpulls_jpsi.png")
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
plt.savefig("OutputFiles/PNGPlots/FitResults/2DPulls_RegularBins_jpsi.png")
#plt.colorbar()
#plt.hist2d(model_b_jpsi_np_counts, m_b, m_jpsi)
#plot_scaling_peakX = n_model_points_per_bin*(frac_yx_YY*n_yx_YY+frac_2x_YY*n_2x_YY+frac_3x_YY*n_3x_YY)/model_m_b_jpsi.sum()
#plt.show()