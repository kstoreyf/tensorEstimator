#!/usr/bin/env python
import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM

from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.utils import compute_amps
from Corrfunc.utils import evaluate_xi

import utils
import corrfuncproj


def main():
    run_dr7_LRGs()


def run_dr7_LRGs():

    print "Running corrfunc estimator on LRGs"

    sample = 'Full'
    nproc = 1
    frac = 0.05
    K = 14
    rmin = 40
    rmax = 180

    datafn = '../data/DR7-{}.ascii'.format(sample)
    randfn = '../data/random-DR7-{}.ascii'.format(sample)

    print "Loading data..."
    data = pd.read_csv(datafn, index_col=0)
    rand = pd.read_csv(randfn, index_col=0)

    #saveto = None
    saveto = "../results/bao/xis_dr7_{}LRG_frac{}.npy".format(sample, frac)
    cosmo = LambdaCDM(H0=70, Om0=0.25,Ode0=0.75)

    #check if need to return or in place
    utils.write_comoving_dist(data, datafn, cosmo)
    utils.write_comoving_dist(rand, randfn, cosmo)


    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)

    data = data.sample(frac=frac)
    #frac *=0.5
    rand = rand.sample(frac=frac)
    #data = data[:int(frac*len(data.index))]
    #rand = rand[:int(frac*len(rand.index))]
    nd = len(data.index)
    nr = len(rand.index)
    print 'ndata=', nd
    print 'nrand=', nr

    weights_data = data['radial_weight']*data['fiber_coll_weight']
    weights_rand = rand['radial_weight']
    # weights_data = None
    # weights_rand = None

    mumax = 1.0 #max of cosine

    sbins = np.linspace(rmin, rmax, K + 1)
    print "bins:", sbins
    ss = []
    xis = []
    labels = []
    counts = []
    aa = []

    print "Running corrfunc..."
    start = time.time()

    dcm_col = 'dcm_Om0-{:.2f}'.format(cosmo.Om0)
    data_cz = data[dcm_col].values
    rand_cz = rand[dcm_col].values
    res = corrfuncproj.counts_smu(data['ra'].values, data['dec'].values,
                                  data_cz, rand['ra'].values, rand['dec'].values,
                                  rand_cz, sbins, mumax, cosmo, nproc=nproc,
                                  weights_data=weights_data, weights_rand=weights_rand,
                                  comoving=True)

    dd, dr, rr, qq, dd_orig, dr_orig, rr_orig = res
    nprojbins = len(sbins)-1
    # Note: dr twice because cross-correlations will be possible
    amps = compute_amps(nprojbins, nd, nd, nr, nr, dd, dr, dr, rr, qq)
    print 'Computed amplitudes'

    amps = np.array(amps)
    svals = np.array(0.5*(sbins[1:]+sbins[:-1]))
    nsvals = len(svals)
    sbins = np.array(sbins)
    nsbins = len(sbins)-1

    xi_proj = evaluate_xi(nprojbins, amps, nsvals, svals, nsbins, sbins)
    end = time.time()
    print 'Time for dr7 {} LRGs, nd={} nr={}: {}'.format(sample, nd, nr, end - start)

    xi_orig = convert_3d_counts_to_cf(nd, nd, nr, nr, dd_orig, dr_orig, dr_orig, rr_orig)


    print "s:", svals
    print "xi_orig", xi_orig
    print "xi_proj", xi_proj

    ss.append(svals)
    xis.append(xi_orig)
    labels.append("corrfunc orig")
    counts.append([dd_orig, dr_orig, rr_orig])
    aa.append(None)

    ss.append(svals)
    xis.append(xi_proj)
    labels.append("corrfunc projected")
    counts.append([dd, dr, rr, qq])
    aa.append(amps)

    if saveto:
        print "Saving to {}".format(saveto)
        np.save(saveto, [ss, xis, labels, counts, aa])


if __name__=="__main__":
    main()
