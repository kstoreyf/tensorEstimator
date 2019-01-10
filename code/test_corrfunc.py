import numpy as np
from os.path import dirname, abspath, join as pjoin
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_rp_pi_counts_to_wp

galaxy_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                     "../mocks/tests/data", "Mr19_mock_northonly.rdcz.ff")

# Read the supplied galaxies on a periodic box
RA, DEC, CZ = read_catalog(galaxy_catalog)
N = len(RA)
down = 100
RA = np.random.choice(RA, N/down)
DEC = np.random.choice(DEC, N/down)
CZ = np.random.choice(CZ, N/down)
N = len(RA)
print(N)

# Read the supplied randoms catalog
random_catalog=pjoin(dirname(abspath(Corrfunc.__file__)),
                     "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff")
rand_RA, rand_DEC, rand_CZ = read_catalog(random_catalog)
rand_N = len(rand_RA)
rand_RA = np.random.choice(rand_RA, rand_N/down)
rand_DEC = np.random.choice(rand_DEC, rand_N/down)
rand_CZ = np.random.choice(rand_CZ, rand_N/down)
rand_N = len(rand_RA)
print(rand_N)

# Setup the bins
nbins = 10
bins = np.linspace(0.1, 20.0, nbins + 1)
pimax = 40.0

cosmology = 1
nthreads = 2

# Auto pair counts in DD
autocorr=1
print('DD')
DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
                         RA, DEC, CZ)


# Cross pair counts in DR
autocorr=0
print('DR')
DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
                         RA, DEC, CZ,
                         RA2=rand_RA, DEC2=rand_DEC, CZ2=rand_CZ)


# Auto pairs counts in RR
autocorr=1
print('RR')
RR_counts = DDrppi_mocks(autocorr, cosmology, nthreads, pimax, bins,
                         rand_RA, rand_DEC, rand_CZ)

# All the pair counts are done, get the angular correlation function
wp = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                DD_counts, DR_counts,
                                DR_counts, RR_counts, nbins, pimax)

print wp
