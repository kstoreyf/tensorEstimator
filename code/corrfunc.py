import numpy as np
import pandas as pd
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import LambdaCDM



cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

def counts(ra_data, dec_data, z_data, ra_rand, dec_rand, z_rand, rpbins, pimax,
         weights_data=None, weights_rand=None, pibinwidth=1, comoving=False):

    assert(len(ra_data)==len(dec_data) and len(ra_data)==len(z_data))
    assert(len(ra_rand)==len(dec_rand) and len(ra_rand)==len(z_rand))

    ndata = len(ra_data)
    nrand = len(ra_rand)
    pibins = np.arange(0, pimax + pibinwidth, pibinwidth)

    if comoving:
        zdf = pd.DataFrame(z_data)
        z_data = zdf.apply(get_comoving_dist)[0].values
        rzdf = pd.DataFrame(z_rand)
        z_rand = rzdf.apply(get_comoving_dist)[0].values

    dd_res_corrfunc = DDrppi_mocks(1, 2, 0, pimax, rpbins, ra_data, dec_data, z_data, is_comoving_dist=comoving)
    dr_res_corrfunc = DDrppi_mocks(0, 2, 0, pimax, rpbins, ra_data, dec_data, z_data,
                                        RA2=ra_rand, DEC2=dec_rand, CZ2=z_rand, is_comoving_dist=comoving)
    rr_res_corrfunc = DDrppi_mocks(1, 2, 0, pimax, rpbins, ra_rand, dec_rand, z_rand, is_comoving_dist=comoving)

    dd_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
    dr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
    rr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))

    for m in range(len(pibins)-1):
        for n in range(len(rpbins)-1):
            idx = (len(pibins)-1) * n + m
            dd_rp_pi_corrfunc[m][n] = dd_res_corrfunc[idx][4]
            dr_rp_pi_corrfunc[m][n] = dr_res_corrfunc[idx][4]
            rr_rp_pi_corrfunc[m][n] = rr_res_corrfunc[idx][4]

    return dd_rp_pi_corrfunc, dr_rp_pi_corrfunc, rr_rp_pi_corrfunc


def calc_wprp_nopi(dd, dr, rr, ndata, nrand):
    dd = np.sum(dd, axis=0)
    dr = np.sum(dr, axis=0)
    rr = np.sum(rr, axis=0)

    est_ls = calc_ls(dd, dr, rr, ndata, nrand)
    #wprp = 2*np.sum(est_ls, axis=0)
    wprp = 2*est_ls

    return est_ls, wprp


def get_comoving_dist(z):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h

def calc_ls(dd_counts, dr_counts, rr_counts, ndata, nrand):
    fN = float(nrand)/float(ndata)
    return (fN*fN*dd_counts - 2*fN*dr_counts + rr_counts)/rr_counts
