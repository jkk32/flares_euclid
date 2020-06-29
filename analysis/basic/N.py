import numpy as np



import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

import FLARE.photom as photom
import FLARE.plt as fplt

flux_limit = photom.m_to_flux(26.) # nJy
print(f'flux_limit/nJy: {flux_limit}')

deltas = np.array([0.969639,0.918132,0.851838,0.849271,0.845644,0.842128,0.841291,0.83945,0.838891,0.832753,0.830465,0.829349,0.827842,0.824159,0.821425,0.820476,0.616236,0.616012,0.430745,0.430689,0.266515,0.266571,0.121315,0.121147,-0.007368,-0.007424,-0.121207,-0.121319,-0.222044,-0.222156,-0.311441,-0.311329,-0.066017,-0.066185,-0.00748,-0.007424,0.055076,0.054909,-0.47874,-0.433818])



cmap = mpl.cm.viridis
norm = mpl.colors.Normalize(vmin=np.min(deltas), vmax=np.max(deltas))


fig = plt.figure(figsize=(3,3))

left  = 0.2
bottom = 0.15
width = 0.75
height = 0.7

ax = fig.add_axes((left, bottom, width, height))


# ---- FLARES

print('-'*20, 'FLARES')

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
h = fl.load_dataset('H', arr_type='Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP')
resims = fl.halos

N_FLARES = []
N_FLARES_resim = {}
for resim in range(len(resims)):
    N_FLARES_resim[resim] = []

for tag in fl.tags:
    H = np.array([])
    for resim in range(len(resims)):
        h_sim = h[resims[resim]][tag]
        # print(ihalo, len(h_sim[h_sim[:]>flux_limit]))
        H = np.append(H, h_sim)
        N_FLARES_resim[resim].append(len(h_sim[h_sim>flux_limit]))

    N = len(H[H>flux_limit])
    N_FLARES.append(N)
    print(tag, N)

for resim in range(len(resims)):

    ax.plot(fl.zeds, np.log10(np.cumsum(N_FLARES_resim[resim])), alpha=0.5, lw=1, c=cmap(norm(deltas[resim])))

ax.plot(fl.zeds, np.log10(np.cumsum(N_FLARES)), c='k', label = 'FLARES')


# ---- EAGLE


print('-'*20, 'EAGLE')

eagle = flares.flares('../../data/EAGLE_REF_sp_info.hdf5', sim_type='PERIODIC')
h = eagle.load_dataset('H', arr_type='Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP')


N_EAGLE = []
for tag in eagle.tags:
    H = h[tag]
    N = len(H[H>flux_limit])
    N_EAGLE.append(N)
    print(tag, N)

ax.plot(eagle.zeds, np.log10(np.cumsum(N_EAGLE)), c='0.5', label = 'EAGLE', ls= '--')




# --- add colorbar

cbar_ax = fig.add_axes((left, bottom+height, width, 0.025))
cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')
cbar.set_label(r'$\rm \delta$')
cbar_ax.xaxis.set_ticks_position('top')
cbar_ax.xaxis.set_label_position('top')


ax.legend()

ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm log_{10}(N(>z))$')

fig.savefig(f'figures/N.pdf')
fig.clf()
