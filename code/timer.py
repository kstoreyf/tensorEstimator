import numpy as np
import time
from matplotlib import pyplot as plt
from astropy.cosmology import LambdaCDM
import pandas as pd

import run
import estimator
import pairs
import estimator_chunks
import pairgen
import pairgenz
import plotter

#globals
#ndata = [10, 31, 102, 307, 1012, 3158, 10015]
#ndata = [10, 31, 102, 307, 1012]
ndata = [10, 31, 102]
#colors = ['blue', 'orange', 'black', 'green']
#colors = ['red', 'orange', 'magenta']
#colors = ['magenta']
colors = ['red', 'orange', 'green', 'blue', 'purple']

def main():
    #time_pairs()
    #plot()
    plotz()

def plotz():

    nds = []
    nrs = []
    ts = []
    ls = []
    rs = []
    ws = []
    lps = []

    #ndmax = 1012
    ndmax = 3158
    nproc = 24
    fn = '../results/times/times_zshells_n{}_nproc{}.npy'.format(ndmax, nproc)
    saveto = 'plots_2019-02-15/times_zshells_n{}_nproc{}.png'.format(ndmax, nproc)

    ndatas, nrands, time_arrs, labels, rps, wprps = np.load(fn)

    print rps
    print len(rps)
    print len(rps[0])
    print wprps
    print len(wprps)
    print len(wprps[0])

    nds += list(ndatas)
    nrs += list(nrands)
    ts += list(time_arrs)
    ls += list(labels)

    # need bc messed up
    for i in range(len(rps)):
        rp = rps[i]
        wprp = wprps[i]
        if type(rp[0]) is np.float64:
            rp = 10 ** rp
            print rp
        else:
            rp = rp[0]
            wprp = wprp[0]
        rs.append(rp)
        ws.append(wprp)
    #rs += list(rps)
    #ws += list(wprps)
    lps += ['corrfunc', 'pairgen', 'pairgenz']*(len(rs)/3)

    #ls += ['chunks, nproc={}'.format(nprocs[i])]

    #plot_times(nrs, ts, ls, saveto=saveto)
    #plot_times(ndatas, time_arrs, labels)
    rs = rs[-3:]
    ws = ws[-3:]
    lps = lps[-3:]

    print len(rs)
    print len(ws)
    print len(lps)


    print rs
    print ws
    print lps

    plotter.plot_wprp(rs, ws, lps, wp_tocompare='corrfunc', colors=['red', 'green', 'orange'])


def plot():
    fns = ['../results/times_chunksonly_nproc4_n307.npy',
           '../results/times_chunksonly_nproc8_n307.npy',
           '../results/times_chunksonly_nproc16_n307.npy']
    nds = []
    ts = []
    ls = []
    nprocs = [4,8,12,16,24]
    for i in range(len(nprocs)):

        #fn = '../results/times_pairs_n{}.npy'.format(max(ndata))
        fn = '../results/times_chunksonly_nproc{}_n307.npy'.format(nprocs[i])
        ndatas, time_arrs, labels = np.load(fn)

        nds += list(ndatas)
        ts += list(time_arrs)
        #ls += list(labels)
        ls += ['chunks, nproc={}'.format(nprocs[i])]

    plot_times(nds, ts, ls)
    #plot_times(ndatas, time_arrs, labels)


def time_est():

    max_only = True

    times = np.zeros(len(ndata))

    K = 10
    pimax = 40. #Mpc/h

    rpmax = 40 #Mpc/h
    basisfunc = estimator.tophat
    rpbins = np.logspace(np.log10(0.1), np.log10(rpmax), K)
    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    wp = True

    for i in range(len(ndata)):

        nd = ndata[i]

        data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
        rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)

        start1 = time.time()
        rps, wprp = run.run(data1fn, rand1fn, data1fn, rand1fn, pimax, rpmax, basisfunc, K, cosmo, wp, rpbins)
        end1 = time.time()
        times[i] = end1 - start1

    time_arrs = [times]

    plot_times(time_arrs, ndata)


def time_pairs():

    times_cf = np.zeros(len(ndata))
    times_pg = np.zeros(len(ndata))
    times_pgz = np.zeros(len(ndata))

    rps = []
    wprps = []
    nrand = []

    nproc = 2

    K = 10
    pimax = 40. #Mpc/h
    pibinwidth = pimax

    min_sep = 0.1
    max_sep = 10. #Mpc/h
    basisfuncs = [estimator_chunks.tophat_robust]
    #bin_sep = np.log(rmax / rmin) / float(K)

    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = run.bins_logavg(rpbins)
    logrpbins_avg = run.logbins_avg(rpbins)
    logwidth = run.log_width(rpbins)

    bin_arg = np.log10(rpbins)

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    wp = True

    for i in range(len(ndata)):

        nd = ndata[i]
        print i, ndata

        data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
        rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)
        data2fn = data1fn
        rand2fn = rand1fn

        data1 = pd.read_csv(data1fn)
        rand1 = pd.read_csv(rand1fn)
        data2 = pd.read_csv(data2fn)
        rand2 = pd.read_csv(rand2fn)

        nrand.append(len(rand1))

        # should make so can take list
        print 'Adding info to dataframes'
        data1 = run.add_info(data1, zfile=None)
        rand1 = run.add_info(rand1, zfile=None)
        data2 = run.add_info(data2, zfile=None)
        rand2 = run.add_info(rand2, zfile=None)


        # start0 = time.time()
        # # run.run_treecorr(data1, rand1, data2, rand2, rmin, rmax, bin_sep, pimax, wp)
        # xi, d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs = run.pairs_treecorr(
        #     data1, rand1, data2, rand2, rmin, rmax, bin_sep, pimax, wp)
        # end0 = time.time()
        # print "Time treecorr pairs:", end0 - start0
        # times_tcp[i] = end0 - start0
        #
        # start1 = time.time()
        # #run.run_treecorr_orig(data1, rand1, data2, rand2, rmin, rmax, bin_sep, pimax, wp)
        # end1 = time.time()
        # print "Time treecorr:", end1 - start1
        # times_tc[i] = end1 - start1

        # start2 = time.time()
        # d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs = pairs.pairs(data1, rand1, data2, rand2,
        #                                                          rmax, cosmo, wp)
        # end2 = time.time()
        # print "Time pairs:", end2 - start2
        # times_kd[i] = end2 - start2

        start = time.time()
        rp, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax,
                                  cosmo, nproc=nproc, pibinwidth=int(pibinwidth))
        end = time.time()
        print "Time corrfunc:", end - start
        times_cf[i] = end - start
        rps.append(logrpbins_avg)
        wprps.append(wprp)

        vals = None

        start = time.time()
        ddgen = pairgen.PairGen(data1, data2, max_sep, cosmo, wp)
        drgen = pairgen.PairGen(data1, rand2, max_sep, cosmo, wp)
        rdgen = pairgen.PairGen(data2, rand1, max_sep, cosmo, wp)
        rrgen = pairgen.PairGen(rand1, rand2, max_sep, cosmo, wp)
        a = estimator_chunks.est(ddgen, drgen, rdgen, rrgen, pimax, max_sep,
                                    cosmo, basisfuncs, K, wp, nproc, bin_arg, logwidth)
        rp, wprp = run.calc_wprp(a,rpbins_avg, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)
        rps.append(rp)
        wprps.append(wprp)
        end = time.time()
        print "Time chunks:", end - start
        times_pg[i] = end - start

        start = time.time()
        ddgen = pairgenz.PairGen(data1, data2, max_sep, cosmo, wp, pimax)
        drgen = pairgenz.PairGen(data1, rand2, max_sep, cosmo, wp, pimax)
        rdgen = pairgenz.PairGen(data2, rand1, max_sep, cosmo, wp, pimax)
        rrgen = pairgenz.PairGen(rand1, rand2, max_sep, cosmo, wp, pimax)
        a = estimator_chunks.est(ddgen, drgen, rdgen, rrgen, pimax, max_sep, cosmo,
                                 basisfuncs, K, wp, nproc, bin_arg, logwidth)
        rp, wprp = run.calc_wprp(a,rpbins_avg, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)
        rps.append(rp)
        wprps.append(wprp)
        end = time.time()
        print "Time chunks:", end - start
        times_pgz[i] = end - start


    # time_arrs = [times_tc, times_kd]
    # labels = ['treecorr', 'kdtree']
    time_arrs = [times_cf, times_pg, times_pgz]
    ndatas = [ndata]*len(time_arrs)
    nrands = [nrand]*len(time_arrs)

    labels = ['corrfunc', 'pairgen', 'pairgen zshells']

    np.save('../results/times/times_zshells_n{}_nproc{}.npy'.format(max(ndata), nproc), [ndatas, nrands, time_arrs, labels, rps, wprps])

    #plot_times(ndata, time_arrs, labels)


def plot_times(ndatas, time_arrs, labels, saveto=None):

    plt.figure()

    for i in range(len(labels)):
        times = time_arrs[i]
        ndata = ndatas[i]

        logndata = np.log10(ndata)

        fit0 = np.polyfit(logndata, np.log10(times), 1)
        yy = np.array(logndata)*fit0[0]+fit0[1]
        plt.plot(logndata, yy, ls='--', color=colors[i])

        label = labels[i] + ': m={:.2f}'.format(fit0[0])
        plt.plot(logndata, np.log10(times), marker='o', label=label, ls='-', color=colors[i])

    plt.legend(loc='best')

    plt.xlabel(r'log(n$_{rand}$)')
    plt.ylabel('log(seconds)')

    if saveto:
        plt.savefig(saveto)
    plt.show()


if __name__=='__main__':
    main()