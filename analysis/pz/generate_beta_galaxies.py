

import numpy as np

import FLARE
cosmo = FLARE.default_cosmo()
import FLARE.SED.models as models
import FLARE.plt as fplt
import FLARE.photom as photo
import FLARE.surveys
import FLARE.filters
import h5py

h5 = h5py.File('../../data2/beta.hdf5', 'w')

filters = FLARE.surveys.surveys['Euclid'].fields['deep'].filters

lamz = np.arange(1000, 20000, 1) # observed frame wavelength
F = FLARE.filters.add_filters(filters, new_lam = lamz)

Ndim = 50


# --- First generate the fluxes for a range of galaxies in \beta and z.

redshifts = np.random.uniform(low=5, high=12, size=Ndim**2)
betas = np.random.uniform(low=-3, high=2, size=Ndim**2)
fluxes = {f: np.zeros(Ndim**2) for f in filters}

for i, (z, beta) in enumerate(zip(redshifts, betas)):
    print(i)
    lam = lamz/(1+z)

    m = models.beta(lam, beta, 1)

    m.get_fnu(cosmo, z)
    Fluxes = m.return_Fnu(F)

    for f in filters:
        fluxes[f][i] = Fluxes[f]


# --- Next repeat for many luminosities

h5.create_dataset('z', data = np.repeat(redshifts, Ndim))
h5.create_dataset('beta', data = np.repeat(betas, Ndim))
h5.create_dataset('log10_luminosity', data = np.random.uniform(low=28, high=31, size=Ndim**3))

for f in filters:
    data = 10**h5['log10_luminosity'][:]*np.repeat(fluxes[f], Ndim)
    h5.create_dataset(f, data = data)

h5.flush()
