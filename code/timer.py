import numpy as np
import time
from matplotlib import pyplot as plt
from astropy.cosmology import LambdaCDM
import pandas as pd

import run
import estimator
import pairs

#globals
ndata = [10, 31, 102, 307, 1012, 3158, 10015]
#ndata = [10, 31, 102]
#ndata = [10, 31]
colors = ['blue', 'red', 'green']

def main():
    time_pairs()


def time_est():

    max_only = True

    times = np.zeros(len(ndata))

    K = 10
    pimax = 40 #Mpc/h
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

    max_only = True

    times_tc = np.zeros(len(ndata))
    times_kd = np.zeros(len(ndata))

    K = 10
    pimax = 40. #Mpc/h
    rmin = 0.1
    rmax = 10. #Mpc/h
    basisfunc = estimator.tophat
    bin_sep = np.log(rmax / rmin) / float(K)
    rpbins = np.logspace(np.log(rmin), np.log(rmax), K+1, base=np.e)
    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    wp = True

    for i in range(len(ndata)):

        nd = ndata[i]

        data1fn = '../../lss/mangler/samples/a0.6452_0001.v5_ngc_ifield_ndata{}.rdzw'.format(nd)
        rand1fn = '../../lss/mangler/samples/a0.6452_rand20x.dr12d_cmass_ngc_ifield_ndata{}.rdz'.format(nd)
        data2fn = data1fn
        rand2fn = rand1fn

        data1 = pd.read_csv(data1fn)
        rand1 = pd.read_csv(rand1fn)
        data2 = pd.read_csv(data2fn)
        rand2 = pd.read_csv(rand2fn)

        # should make so can take list
        print 'Adding info to dataframes'
        data1 = run.add_info(data1, zfile=None)
        rand1 = run.add_info(rand1, zfile=None)
        data2 = run.add_info(data2, zfile=None)
        rand2 = run.add_info(rand2, zfile=None)

        start1 = time.time()
        xi, d1d2pairs_tc, d1r2pairs_tc, d2r1pairs_tc, r1r2pairs_tc = run.pairs_treecorr(data1, rand1, data2, rand2,
                                                                                    rmin, rmax, bin_sep, pimax, wp)
        end1 = time.time()
        print "Time treecorr pairs:", end1 - start1
        times_tc[i] = end1 - start1

        # start2 = time.time()
        # d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs = pairs.pairs(data1, rand1, data2, rand2,
        #                                                          rmax, cosmo, wp)
        # end2 = time.time()
        # print "Time pairs:", end2 - start2
        # times_kd[i] = end2 - start2

    # time_arrs = [times_tc, times_kd]
    # labels = ['treecorr', 'kdtree']
    time_arrs = [times_tc]
    labels = ['treecorr']

    plot_times(time_arrs, labels, ndata)


def plot_times(time_arrs, labels, ndata):

    plt.figure()

    for i in range(len(labels)):
        times = time_arrs[i]

        logndata = np.log10(ndata)

        fit0 = np.polyfit(logndata, np.log10(times), 1)
        yy = np.array(logndata)*fit0[0]+fit0[1]
        plt.plot(logndata, yy, ls='--', color=colors[i])

        label = labels[i] + ': m={:.2f}'.format(fit0[0])
        plt.plot(logndata, np.log10(times), marker='o', label=label, ls='-', color=colors[i])

    plt.legend(loc='best')

    plt.xlabel(r'log(n$_{data}$)')
    plt.ylabel('log(seconds)')

    plt.show()


if __name__=='__main__':
    main()