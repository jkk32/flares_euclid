

import numpy as np
import FLARE.surveys
import h5py
import FLARE.obs.EAZY as eazy

input = h5py.File('../../data2/beta.hdf5', 'r') # contains only \beta, z, and fluxes


filters = FLARE.surveys.surveys['Euclid'].fields['deep'].filters
depths = FLARE.surveys.surveys['Euclid'].fields['deep'].depths


h5 = h5py.File('../../data2/beta_pz.hdf5', 'w')


N = len(input['z'][:])

print(N)

for f in filters:

    noise = depths[f]*np.random.randn(N)
    flux = input[f][:] + noise

    h5.create_dataset(f'obs/{f}/flux', data = flux)
    h5.create_dataset(f'obs/{f}/error', data = depths[f]*np.ones(N))


# ---  run EAZY

F = FLARE.filters.add_filters(filters)
hf_EAZY = eazy.eazy().run_new(h5, F, path = lambda f: f'obs/{f}/')
hf_EAZY.copy(f'EAZY', h5, name = f'EAZY/')
h5.flush()
