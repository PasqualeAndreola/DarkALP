import mplhep
import numpy as np

import uproot
import pandas as pd

import zfit

mplhep.style.use("LHCb2")
import matplotlib as mpl
import matplotlib.pyplot as plt
import hist

from apd import AnalysisData

#datasets = AnalysisData("qee", "btokstarx")
#eventtypes = ["11104498", "11134498", "90000000"]
#eventnames = ["MC_B2Kpiomega", "MC_B2KpiJpsi", "DATA"]
#datatypes  = ["2015", "2016", "2017", "2018"]
#polarities=["MagDown", "MagUp"]
#
#for evttype, evttypename in zip(eventtypes, eventnames):
#    for datatype in datatypes:
#        for polarityit in polarities:
#            with open(f'{evttypename}_{datatype}_{polarityit}', 'w') as outputfilesample:
#                print(f'{evttypename};{datatype};{polarityit};', file=outputfilesample)
#            mcsamples = datasets(eventtype=evttype, datatype=datatype, polarity=polarityit)
#            for index, mcsample in enumerate(mcsamples):
#                with open(f'{evttypename}_{datatype}_{polarityit}', 'a') as outputfilesample:
#                    print(f'{mcsample}', file=outputfilesample)

# Spaces (normalization ranges)

# J/Psi Mass Space
m_jpsi = zfit.Space("x_cons_xpiz_m_best",(600,1000))
m_omega = zfit.Space("x_cons_xpiz_m_best",(600,1000))
m_b = zfit.Space("b_cons_omegaxpiz_m_best",(5000,5600))
m_bpartreco = zfit.Space("b_cons_omegaxpiz_m_best",(5000, 5500))

# Parameters

# Shared signal parameters
# using a dict for the composite params
def param_ratio(params):
    return params["a"] * params["b"]

# J/Psi signal parameters
mu_jpsi = zfit.Parameter(r"$\mu^{J/\psi}$", 782.9, 782.9-100., 782.9+100.)
sigma_jpsi = zfit.Parameter(r"$\sigma^{J/\psi}$", 14., 0., 100.)
alpha_jpsi = zfit.Parameter(r"$\alpha^{J/\psi}$", 1., 0, 6.)
n_jpsi = zfit.param.ConstantParameter(r"$n^{J/\psi}$", 3)
n_jpsi_sig = zfit.Parameter(r"$n^{Sig}$", 200, 0, 50000)

mu2_jpsi = zfit.Parameter(r"$\mu_2^{J/\psi}$", 782.9, 782.9-100., 782.9+100.)
sigmaratio = zfit.Parameter(r"$\sigma_{Ratio}^{J/\psi}$", 0.5, 0., 100.)
sigmaratio2 = sigma_jpsi*sigmaratio
alpha2_jpsi = zfit.Parameter(r"$\alpha_2^{J/\psi}$", -1., -6, 0.)
n2_jpsi = zfit.param.ConstantParameter(r"$n_2^{J/\psi}$", 3)
n2_jpsi_sig = zfit.Parameter(r"$n_2^{Sig}$", 200, 0, 500000)
n_part_reco = zfit.Parameter(r"$n_{PartReco}$", 5000, 0, 500000)

# using a dict for the params
def mult_dict(params):
    return params["a"] / params["b"]

sigmaratio2 = zfit.ComposedParameter('mult_dict', mult_dict, params={"a": sigma_jpsi, "b": sigmaratio})
# Models

# J/Psi mass fit model
model_m_jpsi = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigma_jpsi, alpha=alpha_jpsi, n=n_jpsi)
model_m_jpsi = model_m_jpsi.create_extended(n_jpsi_sig, name='CrystalBall1')
model2_m_jpsi = zfit.pdf.CrystalBall(obs=m_jpsi, mu=mu_jpsi, sigma=sigmaratio2, alpha=alpha2_jpsi, n=n2_jpsi)
model2_m_jpsi = model2_m_jpsi.create_extended(n2_jpsi_sig, name='CrystalBall2')
frac = zfit.Parameter('frac', 0.1, 0., 0.99)
model_m_jpsi = zfit.pdf.SumPDF([model_m_jpsi, model2_m_jpsi])

# B mass fit model
#B signal parameters
mu_b = zfit.Parameter(r"$\mu^{B^0}$", 5279, 5279-60., 5279+60.)
sigma_b = zfit.Parameter(r"$\sigma^{B^0}$", 45., 30., 60.)
alpha_b = zfit.Parameter(r"$\alpha^{B^0}$", -1, -5, -0.1)
n_b = zfit.param.ConstantParameter(r"$n^{B^0}$", 3)
n_b_sig = zfit.Parameter(r"$n_{B^0,1}^{Sig}$", 15000, 0, 1e6)

sigmaratio_b = zfit.Parameter(r"$\sigma_{Ratio}^{B^0}$", 2, 1, 3.)
sigmaratio_b_comp = zfit.ComposedParameter('sigmaratio_b_comp', param_ratio, params={"a": sigmaratio_b, "b": sigma_b})
alpha2_b = zfit.Parameter(r"$\alpha_2^{B^0}$", 1, 0.1, 5)
n2_b = zfit.param.ConstantParameter(r"$n_2^{B^0}$", 3)
n2_b_sig = zfit.Parameter(r"$n_{B^0,2}^{Sig}}$", 15000, 0, 1e6)

# omega mass fit model
#omega signal parameters
mu_omega = zfit.Parameter(r"$\mu^{\omega}$", 782, 782-60., 782+60.)
sigma_omega = zfit.Parameter(r"$\sigma^{\omega}$", 15., 5., 60.)
alpha_omega = zfit.Parameter(r"$\alpha^{\omega}$", -1, -5, -0.1)
n_omega = zfit.param.ConstantParameter(r"$n^{\omega}$", 3)
n_omega_sig = zfit.Parameter(r"$n_{\omega,1}^{Sig}$", 15000, 0, 1e6)

sigmaratio_omega = zfit.Parameter(r"$\sigma_{Ratio}^{\omega}$", 2, 1, 3.)
alpha2_omega = zfit.Parameter(r"$\alpha_2^{\omega}$", 1, 0.1, 5)
n2_omega = zfit.param.ConstantParameter(r"$n_2^{\omega}$", 3)
n2_omega_sig = zfit.Parameter(r"$n_{\omega,2}^{Sig}}$", 15000, 0, 1e6)

model_m_b0_cb1 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigma_b, alpha=alpha_b, n=n_b)
model_m_b0_cb1_extended = model_m_b0_cb1.create_extended(n_b_sig, name='B0_CrystalBall1B')
model_m_b0_cb2 = zfit.pdf.CrystalBall(obs=m_b, mu=mu_b, sigma=sigmaratio_b_comp, alpha=alpha2_b, n=n2_b)
model_m_b0_cb2_extended = model_m_b0_cb2.create_extended(n2_b_sig, name='B0_CrystalBall2B')

sigmaratio_omega_comp = zfit.ComposedParameter('sigmaratio_omega_comp', param_ratio, params={"a": sigmaratio_omega, "b": sigma_omega})
model_m_omega_cb1 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_omega, sigma=sigma_omega, alpha=alpha_omega, n=n_omega)
model_m_omega_cb1_extended = model_m_omega_cb1.create_extended(n_omega_sig, name='omega_CrystalBall1B')
model_m_omega_cb2 = zfit.pdf.CrystalBall(obs=m_omega, mu=mu_omega, sigma=sigmaratio_omega_comp, alpha=alpha2_omega, n=n2_omega)
model_m_omega_cb2_extended = model_m_omega_cb2.create_extended(n2_omega_sig, name='omega_CrystalBall2B')

b_frac = zfit.Parameter('bfrac', 0.1, 0., 100)
model_m_b0_sig = zfit.pdf.SumPDF([model_m_b0_cb1, model_m_b0_cb2], b_frac)
model_m_b0_sig_ext = zfit.pdf.SumPDF([model_m_b0_cb1_extended, model_m_b0_cb2_extended])
model_m_omega_sig_ext = zfit.pdf.SumPDF([model_m_omega_cb1_extended, model_m_omega_cb2_extended])

# J/Psi dataset
file = uproot.open("omegaMCStrippingFiltered.root")
tree = file["DecayTree"]
filter = "(Pi0Merged == 0)"
filter += "& (x_cons_xpiz_m_best<1000) & (x_cons_xpiz_m_best>600)"
filter += "& (b_cons_omegaxpiz_m_best<5500) & (b_cons_omegaxpiz_m_best>5000)"
dataframe = tree.arrays(["b_cons_omegaxpiz_m_best", "x_cons_xpiz_m_best"], filter, library='pd')
print(tree)
mc_b_kpiomega = zfit.Data.from_pandas(df=dataframe, obs=m_b)
mc_omega_3pi = zfit.Data.from_pandas(df=dataframe, obs=m_omega)
#data_rapidsim_b = zfit.Data.from_pandas(df=dataframe, obs=m_bpartreco)

# Plotting

def plot_fit(dat: np.ndarray, basis: np.ndarray, model: np.ndarray, 
             obs: zfit.Space, nbins : int=50, smodel: np.ndarray=None,
             drawstyle: str='default', zmodel: zfit.pdf.BasePDF=None):
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

    # data in histogram over the full observable space
    histo = hist.Hist(hist.axis.Regular(nbins, *limits))
    histo.fill(dat)
    chi2_model_x_np = np.linspace(*obs.limit1d, n_bins)
    chi2_model_np = zmodel.pdf(chi2_model_x_np).numpy() * (zmodel.get_yield().numpy()) *area/nbins
    chi2 = sum( ((histo.counts()-chi2_model_np))**2/chi2_model_np)
    print(r"$\chi^2=$", chi2)
    print("Normalized ",r"$\chi^2=$", chi2/len(params))

    # the figure with an errorbar for the data and a line for the model
    fig, ax = plt.subplots()
    
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
        cmap = plt.get_cmap('autumn') # you can choose whatever you like. 
        norm = mpl.colors.Normalize(0, nmodels) # create a norm for the cmap
        pdfs = [(m.pdf(basis)*m.get_yield()).numpy()*area/nbins
                for m in zmodel.get_models()]
        names = [m.name.replace('_extended','') for m in zmodel.get_models()]
        labels.extend(names)
        for mdex, pdf in enumerate(pdfs):
            artists.append(ax.plot(basis, pdf, color=cmap(norm(mdex)), 
                                   linestyle='--', zorder=-1)[0])
        
    # legend and axis labels
    ax.legend(artists, labels, loc='best', 
              title='LHCb Preliminary\nData Entries: {0}'.format(len(dat)), title_fontsize=18, fontsize=18)
    ax.set_xlabel('Observable')
    ax.set_ylabel('Counts [a.u.]')
    parliststring = 'List of parameters:'
    for ipar, par in enumerate(params):
        if ipar == 0:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'],params[par]['hesse']['error'])
        else:
            parliststring = parliststring+'\n'+r'{0}={1:.2f}$\pm${2:.2f}'.format(par.name, params[par]['value'],params[par]['hesse']['error'])
        width = ax.get_xlim()[1] - ax.get_xlim()[0]
        height = ax.get_ylim()[1] - ax.get_ylim()[0]
    #parliststring = parliststring+'\n'+r'$\chi^2={0:.2f}$'.format(chi2/len(params))
    props = dict(boxstyle='square', facecolor='white', alpha=1)

    # place a text box in upper left in axes coords
    ax.text(0.0045, 0.9955, parliststring, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    return fig

# KDE Partially Reconstructed Decays
npoints = 1000
nbins = 50
x = np.linspace(*m_bpartreco.limit1d, npoints)
binwidth = (5500-5000)/nbins
xwide = np.linspace(*m_bpartreco.limit1d, 50)
""" kde = zfit.pdf.KDE1DimGrid(data_rapidsim_b,
                           obs=m_bpartreco, 
                           padding=0.05,
                           #bandwidth='silverman'
                           #, kernel, num_grid_points,
                           # binning_method, padding, weights, name
                           extended=False,
                           ) """
#kde = kde.create_extended(n_part_reco, name='KDE_PartReco')

""" data_rapidsim_b_np = np.sort(data_rapidsim_b["b_cons_omegaxpiz_m_best"].numpy())
model_x_np = np.linspace(*m_bpartreco.limit1d, 1000)
histo = hist.Hist(hist.axis.Regular(nbins, *m_bpartreco.limit1d))
histo.fill(data_rapidsim_b_np)
entries = histo.counts().sum()
print(entries)
# the figure with an errorbar for the data and a line for the model
fig, ax = plt.subplots()
art_data = ax.errorbar(histo.axes.centers[0], histo.values(), 
                       xerr=histo.axes.widths[0]/2,
                       yerr=np.sqrt(histo.values()), fmt='.', 
                       label='Data', color='black', zorder=10)
plt.plot(x, kde.pdf(x, m_bpartreco).numpy()*entries*binwidth, label='Grid KDE')
plt.legend()
plt.show() """
# Fit

# Stage 1: create an unbinned likelihood with the given PDF and dataset
nll_m_jpsi = zfit.loss.ExtendedUnbinnedNLL(model=model_m_b0_sig_ext, data=mc_b_kpiomega)
nll_m_omega = zfit.loss.ExtendedUnbinnedNLL(model=model_m_omega_sig_ext, data=mc_omega_3pi)

# Stage 2: instantiate a minimiser (in this case a basic minuit)
minimizer = zfit.minimize.Minuit()

# Stage 3: minimise the given negative log-likelihood
result = minimizer.minimize(nll_m_jpsi)
result_omega = minimizer.minimize(nll_m_omega)

# Computing the errors
param_errors = result.hesse(method='approx')
param_errors_omega = result_omega.hesse(method='approx')

# Printing the results
print (result)
print("Function minimum:", result.fmin)
print("Converged:", result.converged)
print("Full minimizer information:", result)

print("Full minimizer information:", result_omega)
# Information on all the parameters in the fit
params = result_omega.params
print(params)

data_mc_jpsi_np = np.sort(mc_b_kpiomega["b_cons_omegaxpiz_m_best"].numpy())
n_bins = 30
m_b_area = m_b.area().numpy()
n_model_points_per_bin = 30
model_x_np = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins)
model_np = model_m_b0_sig_ext.pdf(model_x_np).numpy() * (model_m_b0_sig_ext.get_yield().numpy())

data_mc_omega_np = np.sort(mc_omega_3pi["x_cons_xpiz_m_best"].numpy())
model_m_omega_np = np.linspace(*m_omega.limit1d, n_model_points_per_bin*n_bins)
model_counts_omega_np = model_m_omega_sig_ext.pdf(model_m_omega_np).numpy() * (model_m_omega_sig_ext.get_yield().numpy())
fig = plot_fit(data_mc_omega_np, model_m_omega_np, model_counts_omega_np, m_omega, nbins=n_bins, zmodel=model_m_omega_sig_ext)

data_mc_b_np = np.sort(mc_b_kpiomega["b_cons_omegaxpiz_m_best"].numpy())
model_m_b_np = np.linspace(*m_b.limit1d, n_model_points_per_bin*n_bins)
model_counts_b_np = model_m_b0_sig_ext.pdf(model_m_b_np).numpy() * (model_m_b0_sig_ext.get_yield().numpy())
fig = plot_fit(data_mc_b_np, model_m_b_np, model_counts_b_np, m_b, nbins=n_bins, zmodel=model_m_b0_sig_ext)
fig.savefig('test1.pdf')
fig.savefig('test1.png')
plt.show()
