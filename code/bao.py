#!/usr/bin/env python
import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM

import run
import plotter
import estimator_chunks


c_kms = 3.e5  # c in km/s
czmin = 0.02*c_kms
cz_lims = {7: [30900, 73500], 8: [19900, 47650], 9: [12600, 31900],
           10: [8050, 19250], 11: [5200, 12500], 12: [3200, 7850],
           21:[czmin, 73500], 20:[czmin, 59600], 19:[czmin, 47650], 18:[czmin, 39700],
           17:[czmin, 31900], 16:[czmin, 25450], 15:[czmin, 19250], 14:[czmin, 15750],
           13:[czmin, 12500]}
Mr_lims =  {7: [-23, -22], 8: [-22, -21], 9: [-21, -20],
           10: [-20, -19], 11: [-19, -18], 12: [-18, -17],
            21:-22, 20:-21.5, 19:-21, 18:-20.5, 17:-20, 16:-19.5,
            15:-19, 14:-18.5, 13:-18}
pimaxs = {7: 60, 8: 60, 9: 60, 10: 60, 11: 40, 12: 40,
          21: 60, 20: 60, 19:60, 18: 60, 17: 60, 16: 60, 15: 60, 14: 40, 13: 40}
labels_mr = {7: r'$-23<M_r<-22$', 8: r'$-22<M_r<-21$', 9: r'$-21<M_r<-20$',
            10: r'$-20<M_r<-19$', 11: r'$-19<M_r<-18$', 12: r'$-18<M_r<-17$',
            21: r'$M_r<-22$', 20: r'$M_r<-21.5$', 19: r'$M_r<-21$', 18: r'$M_r<-20.5$',
            17: r'$M_r<-20$', 16: r'$M_r<-19.5$', 15: r'$M_r<-19$', 14: r'$M_r<-18.5$',
            13: r'$M_r<-18$'}
colors = {7:'red', 8:'orange', 9:'green',
          10:'blue', 11:'cyan', 12:'magenta'}




def main():
    run_dr7_LRGs()


def run_dr7_LRGs():

    nproc = 24
    frac = 0.01
    #sample = 'Bright-no'
    #sample = 'Dim-no'
    print "Loading data..."
    sample = 'Full'
    datafn = '../data/DR7-{}.ascii'.format(sample)
    randfn = '../data/random-DR7-{}.ascii'.format(sample)
    data = pd.read_table(datafn, index_col=False, delim_whitespace=True, names=['ra', 'dec', 'z',
            'M_g', 'sector_completeness', 'n(z)*1e4', 'radial_weight', 'fiber_coll_weight',
            'fogtmain', 'ilss', 'icomb', 'sector'], dtype={'z':np.float64}, skiprows=1)
    rand = pd.read_table(randfn, index_col=False,  delim_whitespace=True, names=['ra', 'dec', 'z',
            'sector_completeness', 'n(z)*1e4', 'radial_weight', 'ilss', 'sector'], dtype={'z':np.float64},
             skiprows=1)

    #saveto = None
    saveto = "../results/bao/xis_dr7_{}LRG_frac{}_weights.npy".format(sample, frac)
    cosmo = LambdaCDM(H0=70, Om0=0.25,Ode0=0.75)

    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)
    #Sector completeness already cut to >0.6, not sure if still have to downsample randoms
    #and many have sector completness > 1!! ??
    # data = data[data['z']<0.36]
    # data = data[data['ra']>90][data['ra']<270] #NGC

    data = data.sample(frac=frac)
    rand = rand.sample(frac=frac)
    # data1 = data1[:int(frac*len(data1.index))]
    # rand1 = rand1[:int(frac*len(rand1.index))]
    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)

    print "Adding info..."
    data = run.add_info(data)
    rand = run.add_info(rand)

    weights_data = data['radial_weight']*data['fiber_coll_weight']
    weights_rand = rand['radial_weight']

    losmax = 1.0 #max of cosine
    #losmax = 40.0
    zspace = True
    if sample=='Bright-no':
        K = 21
        rmin = 60
        rmax = 200
    elif sample=='Full':
        K = 14
        rmin = 40
        rmax = 180
    elif sample=='Dim-no':
        K = 15
        rmin = 0.01
        rmax = 8.
    else:
        exit('ERROR')
    bins = np.linspace(rmin, rmax, K + 1)
    #bins = np.logspace(np.log10(rmin), np.log10(rmax), K + 1)
    ss = []
    xis = []
    labels = []

    print "Running corrfunc..."
    start = time.time()
    # or rp, wp
    s, xi = run.run_corrfunc(data, rand, data, rand, bins, losmax, cosmo,
                             weights_data=weights_data, weights_rand=weights_rand, zspace=zspace)
    ss.append(s)
    xis.append(xi)
    labels.append("corrfunc")
    end = time.time()
    print 'Time for dr7 {} LRGs, ndata={}: {}'.format(sample, len(data.index), end-start)

    wp = True
    basisfuncs = [estimator_chunks.tophat_xis]
    bin_arg = bins
    binwidth = (rmax-rmin)/float(K)
    pibinwidth = losmax
    vals = None
    print "Running estimator..."
    s_est, xi_est, a = run.run_chunks(data, rand, data, rand, losmax, rmin,
                                rmax, basisfuncs, K, cosmo, wp, bins,
                                vals, pibinwidth, zspace, nproc, bin_arg, binwidth)
    ss.append(s_est)
    xis.append(xi_est)
    labels.append("est tophat")

    #labels = ['dr7 {} LRGs'.format(sample)]

    if saveto:
        print "Saving to {}".format(saveto)
        np.save(saveto, [ss, xis, labels])
        #run.save_results(saveto, ss, xis, labels)
    #plotter.plot_xi_zspace(ss, xis, labels)


if __name__=="__main__":
    main()
