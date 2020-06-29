import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares

fl = flares.flares('../../data/flares.hdf5', sim_type='FLARES')
halo = fl.halos
tag = fl.tags[-1]  #This would be z=5

print(fl.tags)


volume = (4/3)*np.pi*(fl.radius**3)
quantiles = [0.84,0.50,0.16]

sfr = fl.load_dataset('SFR_inst_30', arr_type='Galaxy')

h = fl.load_dataset('H', arr_type='Galaxy/BPASS_2.2.1/Chabrier300/Flux/DustModelI/Euclid/NISP')

df = pd.read_csv('/cosma/home/dp004/dc-wilk2/data/flare/modules/flares/weight_files/weights_grid.txt')

weights = np.array(df['weights'])
ws, x, y = np.array([]), np.array([]), np.array([])
for ii in range(len(halo)):
    ws = np.append(ws, np.ones(np.shape(h[halo[ii]][tag]))*weights[ii])
    x = np.append(x, np.log10(h[halo[ii]][tag]))
    y = np.append(y, np.log10(sfr[halo[ii]][tag]))

bins = np.arange(0, 4, 0.5)
bincen = (bins[:-1]+bins[1:])/2.
out = flares.binned_weighted_quantile(x,y,ws,bins,quantiles)


fig = plt.figure(figsize=(3,3))

left  = 0.2
bottom = 0.2
width = 0.75
height = 0.75

ax = fig.add_axes((left, bottom, width, height))


ax.legend()

ax.set_xlabel(r'$\rm log_{10}(f_{H}/nJy)$')
ax.set_ylabel(r'$\rm log_{10}(N(>z))$')

fig.savefig(f'figures/SFR.pdf')
fig.clf()
