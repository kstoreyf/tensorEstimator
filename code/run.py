import numpy as np
import pandas as pd
from scipy import interpolate
from astropy.cosmology import LambdaCDM
import time

import treecorr
import corrfunc
from Corrfunc.utils import convert_3d_counts_to_cf

import pairs
import estimator
import plotter


funcs = {'tophat': estimator.tophat, 'tophat_orig': estimator.tophat,
         'top_z': estimator.top_z}

def main():

    nd = 10
    #nd = 31
    #nd = 102
    #nd = 307
    #nd = 3158
    #nd = 10015
    data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
    rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)
    data2fn = data1fn
    rand2fn = rand1fn

    print 'Running for n_data={}'.format(nd)

    K = 12
    #Separations should be given in Mpc/h
    pimax = 40 #Mpc/h
    rpmin = 0.1
    rpmax = 40. #Mpc/h
    bin_sep = np.log(rpmax / rpmin) / float(K) #natural log as treecorr expects


    #basisfuncs = [estimator.tophat, estimator.piecewise, estimator.gaussian, estimator.trig]
    #labels = ['tophat', 'piecewise', 'gaussian', 'trig']

    rpbins = np.logspace(np.log10(rpmin), np.log10(rpmax), K+1)
    logrpbins_avg = logbins_avg(rpbins)
    rpbins_avg = bins_logavg(rpbins)
    logwidth = log_width(rpbins)
    wp = True #vs 3d

    #basisfuncs = [estimator.tophat_fast]
    #bin_arg = np.log10(rpbins)
    #basisfuncs = [estimator.top_z]
    #basisfuncs = [estimator.gauss_z]
    basisfuncs = [estimator.tophat]
    labels = ['tophat']
    bin_arg = logrpbins_avg
    vals = None

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    data1, rand1 = load_data(data1fn, rand1fn)
    data2 = data1
    rand2 = rand1

    #vals = np.linspace(min(data1['z']), max(data1['z']), 8)
    #vals = [0.45, 0.5, 0.55, 0.6, 0.65]
    #labels = ["z={:.2f}".format(val) for val in vals]

    #vals = [0.55]
    #labels =['0.55']
    # vals = [None]
    # labels = ['top']
    #K*=3


    start = time.time()
    rps, wprps = run(data1, rand1, data2, rand2, pimax, rpmin, rpmax,
                bin_sep, basisfuncs, K, cosmo, wp, rpbins, vals, bin_arg, logwidth)

    end = time.time()
    print 'Time run: {:3f} s'.format(end-start)

    #labels += ['tophat_orig']
    #labels += ['treecorr']

    wprps = []
    rps = []
    labels = []
    start = time.time()
    # weights_data = np.full(len(data1.index), 1.0)
    # weights_rand = np.full(len(rand1.index), 1.0)
    # est_corrfunc, wprp_corrfunc = run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo,
    #                             weights_data=weights_data, weights_rand=weights_rand, nopi=True)
    rp_avg, est_corrfunc, wprp_corrfunc = run_corrfunc(data1, rand1, data2, rand2,
                                               rpbins, pimax, cosmo, pibinwidth=pimax)
    print wprp_corrfunc
    end = time.time()
    print 'Time corrfunc: {:3f} s'.format(end-start)
    #rps.append(rpbins_avg)
    rps.append(rp_avg)
    wprps.append(wprp_corrfunc)
    labels.append('corrfunc rp avg')

    rps.append(rpbins_avg)
    wprps.append(wprp_corrfunc)
    labels.append('corrfunc')

    # xi_tree = run_treecorr(data1, rand1, data1, rand1, rpmin, rpmax, bin_sep, pimax, wp)
    # rps.append(rpbins_avg)
    # wprps.append(xi_tree)
    # labels.append('treecorr')

    print bin_sep
    print rpbins

    print len(wprps)
    print len(rps)
    print labels
    colors = ['purple', 'red', 'orange', 'yellow', 'green', 'blue', 'cyan', 'magenta', 'grey']
    plotter.plot_wprp(rps, wprps, labels, colors=colors, wp_tocompare=None)



def load_data(datafn, randfn):
    print 'Loading data'
    data = pd.read_csv(datafn)
    rand = pd.read_csv(randfn)
    print 'Adding info to dataframes'
    data = add_info(data, zfile=None)
    rand = add_info(rand, zfile=None)
    return data, rand


def run(data1, rand1, data2, rand2, pimax, rmin, rmax, bin_sep, basisfuncs,
                   K, cosmo, wp, rpbins, vals, *args):

    print args
    # if wp, rmax means rpmax
    # TODO: now mine returns all up to max and treecorr's returns above rmin too - make consistent!
    #what would happen if I didn't include symmetric pairs for dd, rr?
    #have to do with 2 in dr maybe?
    #ooh poss bc wouldn't do if had 2 diff data catalogs...
    print 'Pairs'
    start = time.time()
    xi, d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs = pairs_treecorr(data1, rand1, data2, rand2,
                                                                rmin, rmax, bin_sep, pimax, wp)
    end = time.time()
    print "Time treecorr pairs:", end-start
    print 'Pairs (DD, DR, RD, RR):', len(d1d2pairs), len(d1r2pairs), len(d2r1pairs), len(r1r2pairs)

    print 'Est'
    if type(basisfuncs)!=list:
        basisfuncs = [basisfuncs]

    multi=True
    if multi:
        a = estimator.est_multi(d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs,
                          data1, rand1, data2, rand2, pimax, rmax, cosmo, basisfuncs, K, wp, *args)
    else:
        a = estimator.est(d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs,
                          data1, rand1, data2, rand2, pimax, rmax, cosmo, basisfuncs, K, wp, *args)

    # TODO: change for not projected? look at corrfunc

    # rps = []
    # wprps = []
    # for bb in range(len(basisfuncs)):
    #     for val in vals:
    #         rp, wprp = calc_wprp(a, basisfuncs[bb], K, rpbins, val, *args)
    #         rps.append(rp)
    #         wprps.append(wprp)
    print 'wprp!'
    rps, wprps = calc_wprp(a, basisfuncs, K, rpbins, vals, *args)

    # just tophat orig
    #rps_orig, wprps_orig = calc_wprp_orig(a, [basisfunc[0]], K, rpbins, *args)
    #print wprps
    # print wprps_orig
    #,
    # rps += rps_orig
    # wprps += wprps_orig

    # rps += rps_orig
    # wprps.append(xi*2)
    #test treecorr

    return rps, wprps


def bins_logavg(bins):
    logbins = np.log10(bins)
    logbins_avg =  0.5 * (logbins[1:] + logbins[:-1])
    avg_bins = 10**logbins_avg
    return avg_bins

def logbins_avg(bins):
    logbins = np.log10(bins)
    logbinsavg =  0.5 * (logbins[1:] + logbins[:-1])
    return logbinsavg

def log_width(bins):
    K = len(bins)-1
    return (np.log10(max(bins)) - np.log10(min(bins)))/float(K)

def calc_wprp(a, basisfunc, K, rpbins, vals, *args):

    #x = np.logspace(np.log10(min(rpbins)), np.log10(max(rpbins)), 500)
    x = bins_logavg(rpbins)
    #x = 0.5*np.array(rpbins[1:]+rpbins[:-1])

    rps = []
    wprps = []
    print 'WPRP'
    for bb in range(len(basisfunc)):

        if vals==None:
            bases = basisfunc[bb](None, None, None, None, x, *args)
            xi_rp = np.matmul(a[bb], bases)
            xi_rp = np.squeeze(np.asarray(xi_rp))
            rps.append(x)
            wprps.append(list(2 * xi_rp))
        else:
            for val in vals:
                bases = basisfunc[bb](None, None, None, None, x, *args, val=val)
                #xi_rp = np.zeros_like(x)
                #for k in range(len(bases)):
                #    xi_rp += a[bb][k]*ba.arrayses[k]
                xi_rp = np.matmul(a[bb], bases)
                xi_rp = np.squeeze(np.asarray(xi_rp))
                rps.append(x)
                wprps.append(list(2*xi_rp))
    return rps, wprps

def calc_wprp_orig(a, basisfunc, K, rpbins, *args):

    rpbins_avg = 0.5*np.array(rpbins[1:]+rpbins[:-1])
    xi_rp = np.zeros(K)
    for i in range(K):
        u = basisfunc(None, None, None, None, rpbins_avg[i], *args)
        xi_rp[i] = np.matmul(a, u)
    wprp = list(2*xi_rp)
    return rpbins_avg, wprp


def run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo, weights_data=None,
                 weights_rand=None, pibinwidth=1, zspace=False):
    print 'Running corrfunc'
    #can only do autocorrelations right now


    ndata = len(data1.index)
    nrand = len(rand1.index)

    if zspace:
        dd, dr, rr = corrfunc.counts(data1['ra'].values, data1['dec'].values, data1['z'].values,
                                     rand1['ra'].values, rand1['dec'].values, rand1['z'].values, rpbins, pimax,
                                     cosmo, weights_data=weights_data, weights_rand=weights_rand, comoving=True, zspace=True)
        xi = convert_3d_counts_to_cf(ndata, ndata, nrand, nrand, dd, dr, dr, rr)
        rp_avg = 0.5*(rpbins[1:]+rpbins[:-1])
        return rp_avg, xi
    else:
        dd, dr, rr, rp_avg = corrfunc.counts(data1['ra'].values, data1['dec'].values, data1['z'].values,
                                             rand1['ra'].values, rand1['dec'].values, rand1['z'].values, rpbins, pimax,
                                             cosmo, weights_data=weights_data, weights_rand=weights_rand, comoving=True)
        est_ls, wprp = corrfunc.calc_wprp(dd, dr, rr, len(data1), len(rand1), pibinwidth=pibinwidth)
        return rp_avg, est_ls, wprp




def run_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp):

    #TODO: make work for 2 data and 2 rand catalogs

    ra = data1['ra'].values
    dec = data1['dec'].values
    dist = data1['z'].apply(get_comoving_dist).values

    ra_rand = rand1['ra'].values
    dec_rand = rand1['dec'].values
    dist_rand = rand1['z'].apply(get_comoving_dist).values

    ndata = len(ra)
    nrand = len(ra_rand)

    print ndata, nrand
    idx = np.arange(ndata, dtype=long)
    idx_rand = np.arange(nrand, dtype=long)

    cat_data = treecorr.Catalog(ra=ra, dec=dec, r=dist, idx=idx, ra_units='deg', dec_units='deg')
    cat_rand = treecorr.Catalog(ra=ra_rand, dec=dec_rand, r=dist_rand, idx=idx_rand, ra_units='deg', dec_units='deg')

    if wp:
        metric = 'Rperp'
    else:
        metric = 'Euclidean'

    dd = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                res_size=ndata**2, min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0, num_threads=1)
    dd.process(cat_data, metric=metric)
    print 'DD pairs:', len(dd.idxpairs1)

    dr = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                res_size=ndata*nrand, min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0, num_threads=1)
    dr.process(cat_data, cat_rand, metric=metric)
    print 'DR pairs:', len(dr.idxpairs1)

    rd = dr

    rr_res = nrand**2
    if rr_res > 1e8:
        rr_res = 1e8
    rr = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                res_size=rr_res, min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0, num_threads=1)
    rr.process(cat_rand, metric=metric)
    print 'RR pairs:', len(rr.idxpairs1)

    xi, varxi = dd.calculateXi(rr, dr)
    print xi
    #TODO: Figure out factor of 2 issue
    # This gives diff answer - maybe because of fac of two issue (?)
    #xils = calc_ls(dd.npairs, dr.npairs, rr.npairs, ndata, nrand)
    #print xils
    #nd = float(len(ra))
    #nr = float(len(ra_rand))
    # These give same answer, and same as treecor calculateXi function
    #xi = (dd.npairs - 2*dr.npairs*(dd.tot/dr.tot) + rr.npairs*(dd.tot/rr.tot))/(rr.npairs*(dd.tot/rr.tot))
    #xi = (dd.npairs*(nr/nd)*(nr/nd) - dr.npairs*(nr/nd) + rr.npairs)/(rr.npairs)
    return xi, dd, dr, rd, rr


def run_treecorr_orig(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp):

    #TODO: make work for 2 data and 2 rand catalogs

    ra = data1['ra'].values
    dec = data1['dec'].values
    dist = data1['z'].apply(get_comoving_dist).values

    ra_rand = rand1['ra'].values
    dec_rand = rand1['dec'].values
    dist_rand = rand1['z'].apply(get_comoving_dist).values

    ndata = len(ra)
    nrand = len(ra_rand)

    cat_data = treecorr.Catalog(ra=ra, dec=dec, r=dist, ra_units='deg', dec_units='deg')
    cat_rand = treecorr.Catalog(ra=ra_rand, dec=dec_rand, r=dist_rand, ra_units='deg', dec_units='deg')

    if wp:
        metric = 'Rperp'
    else:
        metric = 'Euclidean'

    dd = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0)
    dd.process(cat_data, metric=metric)
    print 'DD npairs:', dd.npairs

    dr = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0)
    dr.process(cat_data, cat_rand, metric=metric)
    print 'DR npairs:', dr.npairs

    rd = dr

    rr = treecorr.NNCorrelation(min_sep=min_sep, max_sep=max_sep, bin_size=bin_size,
                                min_rpar=-pimax, max_rpar=pimax,
                                bin_slop=0, num_threads=1)
    rr.process(cat_rand, metric=metric)
    print 'RR npairs:', rr.npairs

    xi, varxi = dd.calculateXi(rr, dr)
    print xi
    return xi, dd, dr, rd, rr


def pairs_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp):

    xi, dd, dr, rd, rr = run_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)

    #double?
    d1d2pairs = zip(dd.idxpairs1, dd.idxpairs2, dd.dists) \
                + zip(dd.idxpairs2, dd.idxpairs1, dd.dists)
    d1r2pairs = zip(dr.idxpairs1, dr.idxpairs2, dr.dists)
                #+ zip(dr.idxpairs2, dr.idxpairs1, dr.dists)
    d2r1pairs = zip(rd.idxpairs1, rd.idxpairs2, rd.dists)
                #+ zip(rd.idxpairs2, rd.idxpairs1, rd.dists)
    r1r2pairs = zip(rr.idxpairs1, rr.idxpairs2, rr.dists) \
                + zip(rr.idxpairs2, rr.idxpairs1, rr.dists)

    return xi, d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs


def add_info(df, zfile=None):
    # Project onto unit sphere
    df['xproj'], df['yproj'], df['zproj'] = zip(*df.apply(ra_dec_to_unitxyz, axis=1))

    if zfile==None:
        zfile = '../tables/z_lookup_lcdm_H070_Om0.3_Ode0.7_4dec.csv'

    zdf = pd.read_csv(zfile)

    # Get comoving distances
    # interp_dcm = interpolate.interp1d(zdf['z_round'], zdf['dcm_mpc'])
    # interp_dcm_transverse = interpolate.interp1d(zdf['z_round'], zdf['dcm_transverse_mpc'])
    # df['dcm_mpc'] = df['z'].apply(interp_dcm)
    # df['dcm_transverse_mpc'] = df['z'].apply(interp_dcm_transverse)
    df['dcm_mpc'] = df['z'].apply(get_comoving_dist)
    df['dcm_transverse_mpc'] = df['dcm_mpc']

    #3d position in Mpc
    df['xpos'], df['ypos'], df['zpos'] = zip(*df.apply(unitxyz_to_xyz, axis=1))

    return df


def calc_rp(cat1, cat2, i, j):

    unitdist = np.sqrt((cat1['xproj'][i] - cat2['xproj'][j])**2
                     + (cat1['yproj'][i] - cat2['yproj'][j])**2
                     + (cat1['zproj'][i] - cat2['zproj'][j])**2)

    #making dcm_transverse average of two but idk what is correct
    dcm_transverse = 0.5*(cat1['dcm_transverse_mpc'][i]+cat2['dcm_transverse_mpc'][j])
    rp = unitdist * dcm_transverse * cosmo.h
    return rp

# Real-space distance in Mpc/h
def calc_r3d(cat1, cat2, i, j):

    dist = np.sqrt((cat1['xpos'][i] - cat2['xpos'][j])**2
                     + (cat1['ypos'][i] - cat2['ypos'][j])**2
                     + (cat1['zpos'][i] - cat2['zpos'][j])**2)
    dist *= cosmo.h
    return dist

# Assumes RA, dec in degrees
def ra_dec_to_unitxyz(row):
    ra = row['ra'] * np.pi / 180.
    dec = row['dec'] * np.pi / 180.
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    return x, y, z

def unitxyz_to_xyz(row):
    d = row['dcm_mpc']
    x = d*row['xproj']
    y = d*row['yproj']
    z = d*row['zproj']
    return x, y, z

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

def get_comoving_dist(z):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h

def calc_ls(dd_counts, dr_counts, rr_counts, ndata, nrand):
    fN = float(nrand)/float(ndata)
    return (fN*fN*dd_counts - 2*fN*dr_counts + rr_counts)/rr_counts

def save_results(fn, rps, wprps, labels):
    np.save(fn, [rps, wprps, labels])

def load_results(fn):
    rps, wprps, labels = np.load(fn)
    return rps, wprps, labels


if __name__=='__main__':
    main()
