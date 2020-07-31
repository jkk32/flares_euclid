import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.photom as photom
import FLARE.plt as fplt


# -- plot of histogram of Euclid accessible FLARES galaxies compared to all galaxies



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


for property in ['SFR_inst_30', 'Mstar_30']:

    print('-'*10, property)

    Y = fl.load_dataset(property, arr_type='Galaxy') # SFR

    fig, axes = plt.subplots(2, len(fl.tags), figsize=(len(fl.tags), 2.5), sharey = False, sharex = True)

    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.9, wspace=0.0, hspace=0.0)

    for i, tag in enumerate(fl.tags): # loop over snapshots (redshifts)

        ax = axes[1, i]
        axC = axes[0, i]

        z = fl.zeds[i]

        ws, x, y = np.array([]), np.array([]), np.array([])
        for ii in range(len(halo)):
            ws = np.append(ws, np.ones(np.shape(H[halo[ii]][tag]))*weights[ii])
            x = np.append(x, np.log10(H[halo[ii]][tag]))
            y = np.append(y, np.log10(Y[halo[ii]][tag]))


        binw = 0.1

        if property == 'SFR_inst_30':
            br = [0, 3]
        if property == 'Mstar_30':
            y += 10. # convert from sub-find to sensible units
            br = [8., 11.]

        bins = np.arange(*br, binw)

        n, Bins, patches = ax.hist(y, weights = ws, bins=bins, color='k', histtype='stepfilled', alpha = 0.2, log = True, bottom = 0.0001)

        # print(np.min(n[n>0]), np.max(n))

        s = x>plot_flux_limit

        # hist, bin_edges = np.histogram(y[s], weights = ws[s], bins = bin_edges)

        ns, Bins, patches = ax.hist(y[s], weights = ws[s], bins=bins, color=cmap(norm(z)), histtype='step', log = True, bottom = 0.0001)

        axC.plot(bins[:-1]+binw/2, ns/n, color=cmap(norm(z)), lw = 1)

        # --- determine point at which Euclid is 95% Complete

        yp_ = np.linspace(*br, 100)

        C_cum = np.array([len(y[(y>yp)&s])/(len(y[y>yp])+1) for yp in yp_])

        if len(yp_[C_cum>0.90])>0:
            c_limit = yp_[C_cum>0.90][0]
            print(z, c_limit)
            axC.axvline(c_limit, c='k', lw=1, alpha=0.2, ls = ':')
            ax.axvline(c_limit, c='k', lw=1, alpha=0.2, ls = ':')


        if i>0:
            ax.set_yticklabels([])
            axC.set_yticklabels([])
        else:
            ax.set_ylabel(r'$\rm N_{w}$')
            axC.set_ylabel(r'$\rm Completeness$')

        ax.set_xlim(br)
        ax.set_ylim([0.0001, 20])

        axC.text(0.5, 1.05, rf'$\rm z={z:.0f}$', color=cmap(norm(z)), fontsize = 7, ha='center', transform=axC.transAxes)


    if property == 'SFR_inst_30':
        label = r'$\rm log_{10}(SFR/M_{\odot}yr^{-1})$'
    if property == 'Mstar_30':
        label = r'$\rm log_{10}(M^{\star}/M_{\odot})$'

    fig.text(0.5, 0.05, label, fontsize = 8, ha='center')


    # ax.set_xlabel(r'$\rm n$')
    # ax.set_xlabel(r'$\rm log_{10}(SFR/M_{\odot} yr^{-1})$')

    fig.savefig(f'figures/hist_{property.replace("/","_")}.pdf')
    fig.clf()
