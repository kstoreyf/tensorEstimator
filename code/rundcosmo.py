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

    #sample = 'Full'
    #DR = 7
    sample = 'ns'
    DR = 3
    nproc = 24
    frac = 1
    randfrac = 1

    #proj_type = 'tophat'
    #proj_type = 'piecewise'
    proj_type = 'generalr'
    proj_label = proj_type+'_fast'
    #projfn = "/home/users/ksf293/Corrfunc/utils/dcosmos.dat"
    #projfn = None
    #projfn = "/home/users/ksf293/vectorizedEstimator/tables/dcosmos.dat" 
    projfn = "/home/users/ksf293/vectorizedEstimator/tables/dcosmos_norm.dat"
    #proj_type = 'dcosmo'
    #proj_label = proj_type+'_const'
    projtag = '_'+proj_label
    if randfrac==1:
      randfractag = ''
    else:
      randfractag = '_rand'+str(randfrac)
    
    saveto = "../results/dcosmo/xis_dr{}_{}LRG_frac{}{}{}.npy".format(DR, sample, frac, randfractag, projtag)

    K = 14
    smin = 40
    smax = 180
    #K=2
    #smin = 5
    #smax = 15
    sbins = np.linspace(smin, smax, K + 1)
    sbinsavg = np.array(0.5*(sbins[1:]+sbins[:-1]))
    print "bins:", sbins
    
    if proj_type=="tophat" or proj_type=="piecewise":
      nprojbins = len(sbins)-1
    elif proj_type=="powerlaw":
      nprojbins = 5
    elif proj_type=="dcosmo":
      nprojbins = 3 #FOR NOW
    elif proj_type=="generalr":
      nprojbins = 3;
    else:
      raise ValueError("Proj type {} not recognized".format(proj_type))
    print "nprojbins:", nprojbins

    datafn = '../data/DR{}-{}.ascii'.format(DR, sample)
    randfn = '../data/random-DR{}-{}.ascii'.format(DR, sample)

    print "Loading data..."
    data = pd.read_csv(datafn, index_col=0)
    rand = pd.read_csv(randfn, index_col=0)

    #saveto = None
    if DR==7:
      cosmo = LambdaCDM(H0=70, Om0=0.25,Ode0=0.75)
    elif DR==3:
      cosmo = LambdaCDM(H0=70, Om0=0.30,Ode0=0.70)
    else:
      raise ValueError("DR {} not recognized".format(DR))

    #check if need to return or in place
    utils.write_comoving_dist(data, datafn, cosmo)
    utils.write_comoving_dist(rand, randfn, cosmo)


    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)

    #data = data.sample(frac=frac)
    #rand = rand.sample(frac=frac*randfrac)
    step = int(1/frac)
    data = data[::step]
    rand = rand[::step]
    #data = data[:int(frac*len(data.index))]
    #rand = rand[:int(frac*len(rand.index))]
    nd = len(data.index)
    nr = len(rand.index)
    print 'ndata=', nd
    print 'nrand=', nr

    #weights_data = data['radial_weight']*data['fiber_coll_weight']
    #weights_rand = rand['radial_weight']
    weights_data = None
    weights_rand = None

    mumax = 1.0 #max of cosine


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
                                  comoving=True, proj_type=proj_type,
                                  nprojbins=nprojbins, projfn=projfn)
    print("Completed counts")
    
    dd, dr, rr, qq, dd_orig, dr_orig, rr_orig = res
    # Note: dr twice because cross-correlations will be possible
    amps = compute_amps(nprojbins, nd, nd, nr, nr, dd, dr, dr, rr, qq)
    print 'Computed amplitudes'

    amps = np.array(amps)
    svals = np.linspace(smin, smax, 300)
    #print "svals:",svals
    #svals = np.array(0.5*(sbins[1:]+sbins[:-1]))
    nsvals = len(svals)
    sbins = np.array(sbins)
    nsbins = len(sbins)-1

    xi_proj = evaluate_xi(nprojbins, amps, nsvals, svals, nsbins, sbins, proj_type, projfn=projfn)
    #print("xi:", xi_proj)
    end = time.time()
    print 'Time for dr7 {} LRGs, nd={} nr={}: {}'.format(sample, nd, nr, end - start)

    xi_orig = convert_3d_counts_to_cf(nd, nd, nr, nr, dd_orig, dr_orig, dr_orig, rr_orig)

    #print "savg:", sbinsavg
    #print "xi_orig", xi_orig
    #print "svals", svals
    #print "xi_proj", xi_proj

    ss.append(sbinsavg)
    xis.append(xi_orig)
    labels.append("orig")
    counts.append([dd_orig, dr_orig, rr_orig])
    aa.append(None)

    ss.append(svals)
    xis.append(xi_proj)
    labels.append(proj_label)
    counts.append([dd, dr, rr, qq])
    aa.append(amps)

    if saveto:
        print "Saving to {}".format(saveto)
        np.save(saveto, [ss, xis, labels, counts, aa])


if __name__=="__main__":
    main()
