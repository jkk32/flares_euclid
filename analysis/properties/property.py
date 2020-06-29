import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.photom as photom
import FLARE.plt as fplt


# --- define colour scale

cmap = mpl.cm.plasma
norm = mpl.colors.Normalize(vmin=5., vmax=10.)


# --- used to define lower x-limit

plot_flux_limit = 2. # log10(nJy)
print(f'log10(flux_limit/nJy): {plot_flux_limit}')




# --- read in FLARS data

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos

# -- define data of interest
# -- CREATE A LOOP OVER PARAMETERS OF INTEREST

H = fl.load_dataset('H', arr_type='Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP') # M_sol/yr

df = pd.read_csv('/cosma/home/dp004/dc-wilk2/data/flare/modules/flares/weight_files/weights_grid.txt')
weights = np.array(df['weights'])


for property in ['SFR_inst_30']:

    Y = fl.load_dataset(property, arr_type='Galaxy') # SFR

    fig = plt.figure(figsize=(3,3))

    left  = 0.2
    bottom = 0.2
    width = 0.75
    height = 0.75

    ax = fig.add_axes((left, bottom, width, height))


    for i, tag in enumerate(fl.tags): # loop over snapshots (redshifts)

        z = fl.zeds[i]

        ws, x, y = np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(H[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(H[halo[ii]][tag]))
            y = np.append(y, np.log10(Y[halo[ii]][tag]))

        s = x>plot_flux_limit

        x = x[s]
        y = y[s]
        ws = ws[s]

        bins = np.arange(1.9, 4, 0.2) # log10(nJy)
        bincen = (bins[:-1]+bins[1:])/2.
        out = flares.binned_weighted_quantile(x,y,ws,bins,[0.84,0.50,0.16])

        N, bin_edges = np.histogram(x, bins=bins)

        Ns = N>10

        ax.plot(bincen, out[:,1], c=cmap(norm(z)), ls= ':')

        ax.plot(bincen[Ns], out[:,1][Ns], c=cmap(norm(z)), label = rf'$\rm z={int(z)}$')
        ax.fill_between(bincen[Ns], out[:,0][Ns], out[:,2][Ns], color=cmap(norm(z)), alpha = 0.2)



    ax.legend()

    ax.set_xlim([2., 4.])

    ax.axvline(np.log10(photom.m_to_flux(26.)), lw=3, c='k', alpha=0.3)

    ax.plot([2.0,4.0],[0.5, 2.5], lw=1, c='k', alpha=0.2)

    ax.set_xlabel(r'$\rm log_{10}(f_{H}/nJy)$')
    ax.set_ylabel(r'$\rm log_{10}(SFR/M_{\odot} yr^{-1})$')

    fig.savefig(f'figures/{property}.pdf')
    fig.clf()
