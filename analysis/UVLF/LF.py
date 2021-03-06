import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.plt as fplt
import FLARE.photom as photom
from astropy.cosmology import Planck15 as cosmo


flares_dir = '/research/astro/flare'

binw = 0.5
bins = np.arange(-25,-16,binw)
binsc = bins[:-1]+binw/2

cbinw = 0.05
cbins = np.arange(-25,-17,cbinw)
cbinsc = cbins[:-1]+cbinw/2




for z in [5.0,6.0,7.0,8.0,9.0,10.0]:

    # -- initialise plot

    fig = plt.figure(figsize = (3.5,3.5))

    left  = 0.15
    bottom = 0.15
    height = 0.6
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    axH = fig.add_axes((left, bottom+height, width, 0.2))



    for sim, ls, alpha, c, label in zip(['ref','flares'],['--','-'], [1.0,1.0], ['0.5','k'],['EAGLE\ REF', 'FLARES']):


        if sim == 'flares':
            fl = flares.flares(f'{flares_dir}/simulations/FLARES/data/flares.hdf5', sim_type='FLARES')
            halo = fl.halos
            tag = fl.tags[fl.zeds.index(z)]

            L = fl.load_dataset('BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/FUV', arr_type='Galaxy') # stellar mass of galaxy

            df = pd.read_csv(f'{flares_dir}/modules/flares/weight_files/weights_grid.txt')
            weights = np.array(df['weights'])
            ws, x = np.array([]), np.array([])
            for ii in range(len(halo)):
                ws = np.append(ws, np.ones(np.shape(L[halo[ii]][tag]))*weights[ii])
                x = np.append(x, photom.lum_to_M(L[halo[ii]][tag]))

            # --- luminosity function

            N_weighted, edges = np.histogram(x, bins = bins, weights = ws)

            h = 0.6777
            vol = (4/3)*np.pi*(14/h)**3
            phi = N_weighted/(binw*vol)

        if sim == 'ref':

            eagle = flares.flares(f'{flares_dir}/simulations/FLARES/data/EAGLE_REF_sp_info.hdf5', sim_type='PERIODIC')

            tag_dict = {10.0: '002_z009p993', 9.0:'003_z008p988', 8.0:'004_z008p075', 7.0:'005_z007p050', 6.0:'006_z005p971', 5.0:'008_z005p037'}

            tag = tag_dict[z]

            L = eagle.load_dataset('BPASS_2.2.1/Chabrier300/Luminosity/DustModelI/FUV', arr_type='Galaxy') # stellar mass of galaxy

            x = photom.lum_to_M(L[tag])

            # --- luminosity function

            N, edges = np.histogram(x, bins = bins)

            vol = (100)**3
            phi = N/(binw*vol)


        ax.plot(binsc, np.log10(phi), ls = ls, lw = 1, c = c, alpha = alpha) # --- plot

        # --- cumulative number of galaxies

        N, edges = np.histogram(x, bins = cbins)


        Nc = np.cumsum(N)

        axH.plot(cbinsc, np.log10(Nc), ls = ls, lw = 1, c = c, alpha = alpha, label = rf'$\rm {label}$') # --- plot

        if sim == 'ref':
            Nc_ref = Nc
        if sim == 'flares':
            Nc_flares = Nc



    # --- surveys

    Surveys = {}

    # --- Euclid
    Surveys['Euclid'] = {}
    Surveys['Euclid']['Deep'] = {'area': 40*3600, 'Hlimit': 26.0, 'public': True}
    Surveys['Euclid']['Wide'] = {'area': 18000*3600, 'Hlimit': 24.0, 'public': True}

    # --- Roman
    # Surveys['Roman'] = {}
    # Surveys['Roman']['HLS'] = {'area': 2000*3600, 'Hlimit': 26.9, 'public': True}
    # Surveys['Roman']['SN-Wide'] = {'area': 18.5*3600, 'Hlimit': 28.7, 'public': True}
    # Surveys['Roman']['SN-Deep'] = {'area': 8.5*3600, 'Hlimit': 29.6, 'public': True}

    # --- Webb
    # Surveys['Webb'] = {}
    # Surveys['Webb']['JADES-Deep'] = {'area': 46., 'Hlimit': 30.7, 'public': False}
    # Surveys['Webb']['JADES-Medium'] = {'area': 144., 'Hlimit': 29.7, 'public': False}

    v1, v2 = cosmo.comoving_volume([z-0.5, z+0.5]).to('Mpc3')

    for Survey in Surveys.keys():

        for subSurvey, SS in Surveys[Survey].items():

            flux = photom.m_to_flux(SS['Hlimit'])
            lum = photom.flux_to_L(flux, cosmo, z)
            M = photom.lum_to_M(lum)

            low_lim = np.log10(1. / ((v2 - v1) * (SS['area']/(41253.*3600))).value)

            ax.fill_between([M, -30], [low_lim, low_lim], alpha=0.1, label = fr'$\rm {{\bf {Survey} }}/{subSurvey}$')

            axH.fill_between([M, -30], [0,0],[10,10], alpha=0.05)


            nc = np.interp(M, cbinsc, Nc_ref)
            if nc>0:
                axH.scatter(M, np.log10(nc), color='0.5', s=5)
                axH.text(M-0.1, np.log10(nc)+0.1, f'${nc:.0f}$', color='0.5', fontsize = 5)

            nc = np.interp(M, cbinsc, Nc_flares)
            if nc>0:
                axH.scatter(M, np.log10(nc), color='k', s=5)
                axH.text(M-0.1, np.log10(nc)+0.1, f'${nc:.0f}$', color='k', fontsize = 5)


    # --- plot parent volumes

    for l, ls in zip([6.4E3, 3.2E3, 0.3E3, 0.1E3],['-','--',':','-.']*2):

        ax.axhline(-3*np.log10(l), color = 'k', ls = ls, lw = 1, alpha = 0.3)

        if l>1E3:
            label = fr'$\rm ({l/1E3:.1f}\ Gpc)^3$'
        else:
            label = fr'$\rm ({int(l)}\ Mpc)^3$'

        ax.text(-23.5, -3*np.log10(l)+0.2, label, fontsize=7, color='0.3')


    ax.legend(loc='center left', fontsize=7)
    axH.legend(loc='upper right', fontsize=7)

    ax.set_xlim([-17.5, -25.])
    ax.set_ylim([-12., -2.1])

    axH.set_xlim([-17.5, -25.])
    axH.set_ylim([0.0, 5.])
    axH.set_xticks([])


    ax.set_xlabel(r'$\rm M_{FUV}$')
    ax.set_ylabel(r'$\rm\log_{10}(\phi/cMpc^{-3}\, mag^{-1})$')

    axH.set_ylabel(r'$\rm\log_{10}[N(>M)]$')


    ax.text(-23.5, -3, fr'$\rm z={z}$', fontsize = 12)

    fig.savefig(f'figures/UVLF_{z}.pdf')
    fig.clf()
