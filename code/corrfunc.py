import numpy as np
import pandas as pd
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
from Corrfunc.utils import convert_rp_pi_counts_to_wp

from astropy.cosmology import LambdaCDM




def counts(ra_data, dec_data, z_data, ra_rand, dec_rand, z_rand, rpbins, pimax,
         cosmo, weights_data=None, weights_rand=None, comoving=False, zspace=False):

    assert(len(ra_data)==len(dec_data) and len(ra_data)==len(z_data))
    assert(len(ra_rand)==len(dec_rand) and len(ra_rand)==len(z_rand))

    ndata = len(ra_data)
    nrand = len(ra_rand)
    pibinwidth = 1
    pibins = np.arange(0, pimax + pibinwidth, pibinwidth)

    if comoving:
        zdf = pd.DataFrame(z_data)
        z_data = zdf.apply(get_comoving_dist, args=(cosmo,))[0].values
        rzdf = pd.DataFrame(z_rand)
        z_rand = rzdf.apply(get_comoving_dist, args=(cosmo,))[0].values

    if weights_data is None:
        weights_data = np.ones(ndata)
    if weights_rand is None:
        weights_rand = np.ones(nrand)

    cosmology = 1
    nthreads = 4
    verbose = False
    weight_type = 'pair_product'
    output_rpavg = True

    if zspace:
        mumax = pimax
        nmubins = 1
        print 'Computing DD pairs'
        dd_res_corrfunc = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                       weights1=weights_data, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type)
        print 'Computing DR pairs'
        dr_res_corrfunc = DDsmu_mocks(0, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                            RA2=ra_rand, DEC2=dec_rand, CZ2=z_rand, weights1=weights_data,
                                       weights2=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type)
        print 'Computing RR pairs'
        rr_res_corrfunc = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_rand, dec_rand, z_rand,
                                       weights1=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type)


        return dd_res_corrfunc, dr_res_corrfunc, rr_res_corrfunc

    else:

        print 'Computing DD pairs'
        dd_res_corrfunc = DDrppi_mocks(1, cosmology, nthreads, pimax, rpbins, ra_data, dec_data, z_data,
                                       weights1=weights_data, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, output_rpavg=output_rpavg)
        print 'Computing DR pairs'
        dr_res_corrfunc = DDrppi_mocks(0, cosmology, nthreads, pimax, rpbins, ra_data, dec_data, z_data,
                                            RA2=ra_rand, DEC2=dec_rand, CZ2=z_rand, weights1=weights_data,
                                       weights2=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, output_rpavg=output_rpavg)
        print 'Computing RR pairs'
        rr_res_corrfunc = DDrppi_mocks(1, cosmology, nthreads, pimax, rpbins, ra_rand, dec_rand, z_rand,
                                       weights1=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, output_rpavg=output_rpavg)

        dd_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
        dr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
        rr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))

        rp_avg = np.zeros(len(rpbins)-1)
        for m in range(len(pibins)-1):
            for n in range(len(rpbins)-1):
                idx = (len(pibins)-1) * n + m
                # = count * avg weight
                dd_rp_pi_corrfunc[m][n] = dd_res_corrfunc[idx][4]*dd_res_corrfunc[idx][5]
                dr_rp_pi_corrfunc[m][n] = dr_res_corrfunc[idx][4]*dr_res_corrfunc[idx][5]
                rr_rp_pi_corrfunc[m][n] = rr_res_corrfunc[idx][4]*rr_res_corrfunc[idx][5]
                rp_avg[n] += dd_res_corrfunc[idx][2] + dr_res_corrfunc[idx][2] + \
                             rr_res_corrfunc[idx][2]
        for n in range(len(rpbins) - 1):
            rp_avg[n] /= 3*(len(pibins)-1)

        wprp = convert_rp_pi_counts_to_wp(ndata, ndata, nrand, nrand, dd_res_corrfunc, dr_res_corrfunc,
                                 dr_res_corrfunc, rr_res_corrfunc, len(rpbins)-1, pimax, dpi=pibinwidth)

        return dd_rp_pi_corrfunc, dr_rp_pi_corrfunc, rr_rp_pi_corrfunc, rp_avg, wprp


def calc_wprp(dd, dr, rr, ndata, nrand, pibinwidth=1):

    assert type(pibinwidth) == int
    assert pibinwidth >= 1
    assert len(dd)%float(pibinwidth) == 0

    #reshape into different bin widths
    dd = dd.reshape(-1, pibinwidth, dd.shape[-1]).sum(axis=1)
    dr = dr.reshape(-1, pibinwidth, dr.shape[-1]).sum(axis=1)
    rr = rr.reshape(-1, pibinwidth, rr.shape[-1]).sum(axis=1)

    est_ls = calc_ls(dd, dr, rr, ndata, nrand)
    wprp = 2*pibinwidth*np.sum(est_ls, axis=0)
    #wprp = 2*np.sum(est_ls, axis=0)

    return est_ls, wprp


def get_comoving_dist(z, cosmo):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h

def calc_ls(dd_counts, dr_counts, rr_counts, ndata, nrand):
    fN = float(nrand)/float(ndata)
    return (fN*fN*dd_counts - 2*fN*dr_counts + rr_counts)/rr_counts
