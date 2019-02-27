import numpy as np
import pandas as pd
import time

from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from Corrfunc.mocks.DDsmu_mocks import DDsmu_mocks
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.utils import compute_amps
from Corrfunc.utils import evaluate_xi

#from Corrfunc.mocks.DDsmu_mocks import convert_3d_proj_counts_to_amplitude
#from Corrfunc._countpairs_mocks import countpairs_s_mu_mocks as s_mu_mocks
#from Corrfunc._countpairs_mocks import convert_3d_proj_counts_to_amplitude as compute_amplitude

from astropy.cosmology import LambdaCDM




def counts(ra_data, dec_data, z_data, ra_rand, dec_rand, z_rand, rpbins, pimax,
         cosmo, nproc=1, weights_data=None, weights_rand=None,
           comoving=False, zspace=False, proj=False):

    if proj:
        if not zspace:
            print "Projected corrfunc only available with smu (zspace==True)"
            exit(1)
        return counts_proj(ra_data, dec_data, z_data, ra_rand, dec_rand, z_rand, rpbins,
                    pimax, cosmo, nproc=nproc, weights_data=weights_data, weights_rand=weights_rand,
                    comoving=comoving)


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
    nthreads = nproc
    verbose = False
    weight_type = 'pair_product'
    output_rpavg = True
    isa = 'fallback'

    if zspace:
        mumax = pimax
        nmubins = 1

        print 'Computing DD pairs'
        start = time.time()
        dd_res_corrfunc, dd_proj, dd_projt = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                       weights1=weights_data, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, isa=isa)
        end = time.time()
        print "Time DD pairs:", end-start
        print "DD:",dd_proj

        print 'Computing DR pairs'
        start = time.time()
        dr_res_corrfunc, dr_proj, dr_projt = DDsmu_mocks(0, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                            RA2=ra_rand, DEC2=dec_rand, CZ2=z_rand, weights1=weights_data,
                                       weights2=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, isa=isa)
        end = time.time()
        print "Time DR pairs:", end-start
        print "DR:",dr_proj

        print 'Computing RR pairs'
        start = time.time()
        rr_res_corrfunc, rr_proj, rr_projt = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_rand, dec_rand, z_rand,
                                       weights1=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                       weight_type=weight_type, isa=isa)
        end = time.time()
        print "Time RR pairs:", end-start
        print "RR",rr_proj
        print "QQ",rr_projt

        # just for printing
        dd_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
        dr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
        rr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))

        for m in range(len(pibins)-1):
            for n in range(len(rpbins)-1):
                idx = (len(pibins)-1) * n + m
                # = count * avg weight
                dd_rp_pi_corrfunc[m][n] = dd_res_corrfunc[idx][4]*dd_res_corrfunc[idx][5]
                dr_rp_pi_corrfunc[m][n] = dr_res_corrfunc[idx][4]*dr_res_corrfunc[idx][5]
                rr_rp_pi_corrfunc[m][n] = rr_res_corrfunc[idx][4]*rr_res_corrfunc[idx][5]
        print "DD, DR, RR orig:"
        print dd_rp_pi_corrfunc, sum(dd_rp_pi_corrfunc.flatten())
        print dr_rp_pi_corrfunc, sum(dr_rp_pi_corrfunc.flatten())
        print rr_rp_pi_corrfunc, sum(rr_rp_pi_corrfunc.flatten())

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

        # wprp = convert_rp_pi_counts_to_wp(ndata, ndata, nrand, nrand, dd_res_corrfunc, dr_res_corrfunc,
        #                          dr_res_corrfunc, rr_res_corrfunc, len(rpbins)-1, pimax, dpi=pibinwidth)

        #np.save('../results/counts_bin20_frac0.1_randcz.npy',
        #       [dd_res_corrfunc, dr_res_corrfunc, rr_res_corrfunc])

        return dd_rp_pi_corrfunc, dr_rp_pi_corrfunc, rr_rp_pi_corrfunc, rp_avg


def counts_proj(ra_data, dec_data, z_data, ra_rand, dec_rand, z_rand, rpbins,
            losmax, cosmo, nproc=1, weights_data=None, weights_rand=None,
            comoving=False):

    assert(len(ra_data)==len(dec_data) and len(ra_data)==len(z_data))
    assert(len(ra_rand)==len(dec_rand) and len(ra_rand)==len(z_rand))

    nd1 = len(ra_data)
    nr1 = len(ra_rand)
    nd2 = nd1
    nr2 = nr1
    print "counts proj"
    print "comoving:",comoving
    if comoving:
        print "To comoving:"
        zdf = pd.DataFrame(z_data)
        print min(z_data), max(z_data)
        z_data = zdf.apply(get_comoving_dist, args=(cosmo,))[0].values
        print min(z_data), max(z_data)
        rzdf = pd.DataFrame(z_rand)
        print min(z_rand), max(z_rand)
        z_rand = rzdf.apply(get_comoving_dist, args=(cosmo,))[0].values
        print min(z_rand), max(z_rand)

    if weights_data is None:
        weights_data = np.ones(nd1)
    if weights_rand is None:
        weights_rand = np.ones(nr1)

    cosmology = 1
    nthreads = nproc
    verbose = False
    weight_type = 'pair_product'
    output_rpavg = True
    isa = 'fallback'

    mumax = losmax
    nmubins = 1

    print 'test'

    # dd_res_corrfunc, dd_proj, dd_projt = s_mu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
    #                               weights1=weights_data, is_comoving_dist=comoving, verbose=verbose,
    #                               weight_type=weight_type, isa=isa)


    print 'Computing DD pairs'
    start = time.time()
    dd_res_corrfunc, dd_proj, dd_projt = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                   weights1=weights_data, is_comoving_dist=comoving, verbose=verbose,
                                   weight_type=weight_type, isa=isa)
    end = time.time()
    print "Time DD pairs:", end - start
    print "DD:", dd_proj
    print dd_res_corrfunc

    print 'Computing DR pairs'
    start = time.time()
    dr_res_corrfunc, dr_proj, dr_projt = DDsmu_mocks(0, cosmology, nthreads, mumax, nmubins, rpbins, ra_data, dec_data, z_data,
                                        RA2=ra_rand, DEC2=dec_rand, CZ2=z_rand, weights1=weights_data,
                                   weights2=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                   weight_type=weight_type, isa=isa)
    end = time.time()
    print "Time DD pairs:", end - start
    print "DR:", dr_proj
    print dr_res_corrfunc

    print 'Computing RR pairs'
    start = time.time()
    rr_res_corrfunc, rr_proj, rr_projt = DDsmu_mocks(1, cosmology, nthreads, mumax, nmubins, rpbins, ra_rand, dec_rand, z_rand,
                                   weights1=weights_rand, is_comoving_dist=comoving, verbose=verbose,
                                   weight_type=weight_type, isa=isa)
    end = time.time()
    print "Time RR pairs:", end - start
    print "RR:", rr_proj
    print "QQ:", rr_projt
    print rr_res_corrfunc

    # a = []
    # dd = dd_proj * 1. / (nd1 * nd2)
    # dr = dr_proj * 1. / (nd1 * nr2)
    # #TODO: allow cross-correlations
    # rd = dr_proj * 1. / (nd2 * nr1)
    # rr = rr_proj * 1. / (nr1 * nr2)
    # qq = rr_projt * 1. / (nr1 * nr2)
    #a.append(calc_amplitudes(dd, dr, rd, rr, qq))

    nprojbins = len(rpbins)-1
    amps = compute_amps(nprojbins, nd1, nd2, nr1, nr2, dd_proj, dr_proj, dr_proj, rr_proj, rr_projt)
    print 'Computed amplitudes'
    print amps
    print amps
    amps = np.array(amps)
    svals = np.array(0.5*(rpbins[1:]+rpbins[:-1]))
    nsvals = len(svals)
    sbins = np.array(rpbins)
    nsbins = len(rpbins)-1
    print nprojbins
    print amps
    print nsvals
    print svals
    print nsbins
    print sbins

    xi = evaluate_xi(nprojbins, amps, nsvals, svals, nsbins, sbins)

    print "Evaluated xi"
    print xi

    return dd_res_corrfunc, dr_res_corrfunc, rr_res_corrfunc, svals, xi, amps



def split(a, n):
    k, m = divmod(len(a), n)
    return [(i * k + min(i, m), (i + 1) * k + min(i + 1, m)) for i in xrange(n)]
    # return list((a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)))


def calc_amplitudes(dd, dr, rd, rr, qq):
    print 'Calculating amplitudes'
    qqinv = np.matrix(qq).I
    numerator = dd - dr - rd + rr
    a = np.matmul(qqinv, numerator)
    return a.A1




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
