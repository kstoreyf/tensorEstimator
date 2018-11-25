import numpy as np
import pandas as pd
from scipy import interpolate
from astropy.cosmology import LambdaCDM
import time

import pairs
import estimator
import plotter
import corrfunc


def main():

    print 'Running'
    nd = 31
    data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
    rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)
    data2fn = data1fn
    rand2fn = rand1fn

    K = 10
    pimax = 40 #Mpc/h
    rpmin = 0.1
    rpmax = 40 #Mpc/h
    #basisfunc = estimator.piecewise
    #basisfunc = estimator.tophat
    basisfunc = estimator.gaussian
    rpbins = np.logspace(np.log10(rpmin), np.log10(rpmax), K)
    #rpbins_fine = np.logspace(np.log10(0.1), np.log10(rpmax), K)
    wp = True

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

    start = time.time()
    rps, xi_rp, wprp = run(data1fn, rand1fn, data2fn, rand2fn, pimax, rpmax,
                           basisfunc, K, cosmo, wp, rpbins)
    end = time.time()

    est_corrfunc, wprp_corrfunc = run_corrfunc(data1fn, rand1fn, data2fn, rand2fn, rpbins, pimax)

    print 'Time: {:3f} s'.format(end-start)

    plotter.plot_wprp([rpbins, rpbins], [wprp, wprp_corrfunc], ['gaussian', 'corrfunc'])


def run(data1fn, rand1fn, data2fn, rand2fn, pimax, rpmax, basisfunc, K, cosmo, wp, *args):
    print 'Loading data'
    data1 = pd.read_csv(data1fn)
    rand1 = pd.read_csv(rand1fn)
    data2 = pd.read_csv(data2fn)
    rand2 = pd.read_csv(rand2fn)

    # should make so can take list
    print 'Adding info to dataframes'
    data1 = add_info(data1, zfile=None)
    rand1 = add_info(rand1, zfile=None)
    data2 = add_info(data2, zfile=None)
    rand2 = add_info(rand2, zfile=None)

    d1d2pairs, d1r1pairs, d2r2pairs, r1r2pairs = pairs.pairs(data1, rand1, data2, rand2,
                                                             pimax, rpmax, cosmo, wp)


    a = estimator.est(d1d2pairs, d1r1pairs, d2r2pairs, r1r2pairs,
                      data1, rand1, data2, rand2, pimax, rpmax, cosmo, basisfunc, K, wp, *args)

    # TODO: change for not projected? look at corrfunc
    rps, xi_rp, wprp = calc_wprp(a, basisfunc, K, *args)

    return rps, xi_rp, wprp


def calc_wprp(a, basisfunc, K, *args):
    rpbins = args[0]
    rps = 0.5*np.array(rpbins[1:]+rpbins[:-1])
    print rps

    xi_rp = np.zeros(K-1)
    for i in range(K-1):
        u = basisfunc(None, None, rps[i], *args)
        xi_rp[i] = np.matmul(a, u)[0][0]
        #xi_rp = np.array([[np.matmul(ael.T, u) for ael in arow] for arow in a])
    #wprp = 2*np.sum(xi_rp, axis=0)
    wprp = 2*xi_rp
    print a
    print wprp
    return rps, xi_rp, wprp


def run_corrfunc(data1fn, rand1fn, data2fn, rand2fn, rpbins, pimax):
    print 'Loading data'
    data1 = pd.read_csv(data1fn)
    rand1 = pd.read_csv(rand1fn)
    data2 = pd.read_csv(data2fn)
    rand2 = pd.read_csv(rand2fn)

    #can only do autocorrelations right now
    dd, dr, rr = corrfunc.counts(data1['ra'].values, data1['dec'].values, data1['z'].values,
                    rand1['ra'].values, rand1['dec'].values, rand1['z'].values,
                    rpbins, pimax, comoving=True)

    est_ls, wprp = corrfunc.calc_wprp_nopi(dd, dr, rr, len(data1), len(rand1))

    return est_ls, wprp

def add_info(df, zfile=None):
    # Project onto unit sphere
    df['xproj'], df['yproj'], df['zproj'] = zip(*df.apply(ra_dec_to_unitxyz, axis=1))

    if zfile==None:
        zfile = '../tables/z_lookup_lcdm_H070_Om0.3_Ode0.7_4dec.csv'

    zdf = pd.read_csv(zfile)

    # Get comoving distances
    interp_dcm = interpolate.interp1d(zdf['z_round'], zdf['dcm_mpc'])
    interp_dcm_transverse = interpolate.interp1d(zdf['z_round'], zdf['dcm_transverse_mpc'])
    df['dcm_mpc'] = df['z'].apply(interp_dcm)
    df['dcm_transverse_mpc'] = df['z'].apply(interp_dcm_transverse)

    return df



# Assumes RA, dec in degrees
def ra_dec_to_unitxyz(row):
    ra = row['ra'] * np.pi / 180.
    dec = row['dec'] * np.pi / 180.
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    return x, y, z


if __name__=='__main__':
    main()