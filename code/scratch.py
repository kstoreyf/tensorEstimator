import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from astropy.cosmology import LambdaCDM
import timeit
from time import time

import run


def main():
    scratch_args()
    #scratch_continuous()
    #piecewise_check()
    #piece_cf()
    #plot_bases()


def average_time(executable, *args, **kwargs):
    """Compute the average time over N runs"""
    N = 5
    t = 0
    for i in range(N):
        t0 = time()
        res = executable(*args, **kwargs)
        t1 = time()
        t += (t1 - t0)
    return res, t * 1. / N



def tophat(cat1, cat2, i, j, rp, logbins_avg, logwidth):
    logrp = np.log10(rp)
    # u = np.zeros(len(logbins_avg))
    # for i in range(len(u)):
    #     peak = logbins_avg[i]
    #     u[i] = top(logrp, peak, logwidth)
    u = np.array([top(logrp, peak, logwidth) for peak in logbins_avg])
    return u

def tophat_orig(rp, rpbins):
    u = np.zeros(len(rpbins)-1)
    pos = np.digitize(rp, rpbins)
    if pos>0 and pos<len(u):
        u[pos - 1] = 1
    return u

def top(x, peak, width):
    t = np.zeros_like(x)
    t[(x>peak-width/2.)&(x<peak+width/2.)] = 1
    return t

def piece(x, peak, width):
    p = np.array(1 - (1. / width) * abs(x - peak))
    p[(x<peak-width)^(x>peak+width)] = 0
    return p

def gauss(x, peak, width):
    sigma = width/(2.*np.sqrt(2*np.log(2.)))
    return np.array(np.exp(-(x - peak) ** 2 / (2. * sigma ** 2)))

def trig_u(x, N):
    bases = []
    logx = np.log10(x)
    for n in range(1, N+1):
        s = np.sin(n*logx)
        c = np.cos(n*logx)
        bases.append(s)
        bases.append(c)
    return bases

def trig(rp, logbins_avg, logwidth):
    N = len(logbins_avg)/2
    kbins = np.logspace(-1, 1, N)
    if np.isscalar(rp) and rp==0:
        return np.zeros(N*2)
    u = []
    logrp = np.log10(rp)
    L = max(logbins_avg)
    #for n in range(1, N+1):
    for i in range(len(kbins)):
        kmode = kbins[i]
        s = np.sin(2.*np.pi*kmode*logrp / L)
        c = np.cos(2.*np.pi*kmode*logrp / L)
        u.append(s)
        u.append(c)
    return np.array(u)

def gauss_u(x, bins):
    bins_avg = 0.5 * (bins[1:] + bins[:-1])
    logwidth = np.log10(bins_avg[1]) - np.log10(bins_avg[0])
    u = [gauss(np.log10(x), np.log10(peak), logwidth) for peak in bins_avg]
    return u

def tophat_u(x, bins):
    bins_avg = 0.5 * (bins[1:] + bins[:-1])
    logwidth = np.log10(bins_avg[1]) - np.log10(bins_avg[0])
    u = [top(np.log10(x), np.log10(peak), logwidth) for peak in bins_avg]
    return u

def piece_u(x, bins):
    bins_avg = 0.5 * (bins[1:] + bins[:-1])
    logwidth = np.log10(bins_avg[1]) - np.log10(bins_avg[0])
    u = [piece(np.log10(x), np.log10(peak), logwidth) for peak in bins_avg]
    return u


def tophat_arr(rp, rpbins):
    u = np.zeros(len(rpbins)-1)
    pos = np.digitize(rp, rpbins)

    if pos>0 and pos<len(u):
        u[pos - 1] = 1
    return u

def tophat_arr(rp, rpbins):
    u = np.zeros(len(rpbins)-1)
    pos = np.digitize(rp, rpbins)
    if pos>0 and pos<len(u):
        u[pos - 1] = 1
    return u

def gaussian(rp, bins):
    mu = rp
    pos = np.digitize(rp, bins)
    rploc = pos - 1
    sigma = bins[rploc] - bins[rploc - 1]
    return np.array(np.exp(-(bins - mu) ** 2 / (2. * sigma ** 2)))

def quadratic(rp, bins):
    pos = np.digitize(rp, bins)
    rploc = pos - 1
    a = bins[rploc] - bins[rploc - 1]
    u = -a*(bins-rp)**2 + 1
    u = [0 if uu<0 else uu for uu in u]
    return u



def piecewise(rp, rpbins, rpbins_fine):
    pos = np.digitize(rp, rpbins)
    rploc = pos-1
    u = 1-abs(rpbins_fine-rploc)
    u = [0 if uu<0 else uu for uu in u]
    return u

def piecewise_fine(rp, bins):
    #bins_avg = 0.5*(bins[1:]-bins[:-1])
    pos = np.digitize(rp, bins)
    rploc = pos-1
    width = bins[rploc]-bins[rploc-1]
    # print rpbins_fine
    # print rploc
    # print rpbins_fine-rploc
    u = 1-(1./(1*width))*abs(bins-rp)
    #u = 1 - abs(bins-rp)
    u = [0 if uu<0 else uu for uu in u]

    return u



def piecewise_exact(rp, bins):

    pos = np.digitize(rp, bins)
    rploc = pos-1

    u = 1 - abs(bins-rp)
    u = [0 if uu < 0 else uu for uu in u]

    bins_avg = 0.5 * (bins[1:] - bins[:-1])

    def piece(x, rploc):
        if rploc==0:
            rploc += 1
        elif rploc==len(bins_avg):
            rploc -=1
        if x<rp:
            width = bins_avg[rploc] - bins_avg[rploc+1]
        if x>=rp:
            width = bins_avg[rploc-1] - bins_avg[rploc]
        return 1-(1./(1*width))*abs(bins-rp)

    return piece(bins_avg, rploc, )




def get_rps(rpbins_fine, rpmax):
    bins = rpbins_fine[1:-1]
    rps = bins
    K = len(bins)
    div = 10
    for i in range(div-1):
        end = (i+1)*K/div
        rps = np.concatenate((rps, bins[:end]))
    return rps

#def piecewise_log(rploc, )

def piece_log_arr(x, peak, width):
    logx = np.log10(x)
    logpeak = np.log10(peak)
    tri = 1 - (1. / (width)) * abs(logx - logpeak)
    return [t if t>0 else 0 for t in tri]

def piece_log_bins(x, bins):
    bins_avg = 0.5 * (bins[1:] - bins[:-1])
    width = np.log10(bins_avg[1]) - np.log10(bins_avg[0])
    funcs = []
    for peak in bins_avg:
        funcs.append(piece_log_arr(x, peak, width))
    return funcs

def piecewise_log(x, bins):
    bins_avg = 0.5 * (bins[1:] - bins[:-1])
    logwidth = np.log10(bins_avg[1]) - np.log10(bins_avg[0])

    def piece_log(x, peak, width):
        logx = np.log10(x)
        logpeak = np.log10(peak)
        if logx > logpeak+logwidth or logx < logpeak-logwidth:
            return 0
        return 1 - (1. / logwidth) * abs(logx - logpeak)

    u = [piece_log(x, peak, logwidth) for peak in bins_avg]

    return u

def piecewise_check():
    K = 10
    rpmax = 10
    bins = np.logspace(np.log10(0.1), np.log10(rpmax), K)

    rp = 1
    u = piecewise_log(rp, bins)
    print u

    x = np.logspace(-2, 2, 1000)
    funcs = piece_log_bins(x, bins)
    for func in funcs:
        plt.semilogx(x, func, color='blue')
        plt.axvline(rp, color='red')

    # for i in range(len(bins_avg)-1):
    #     peak = bins_avg[i]
    #     width = 2*(np.log10(bins_avg[i+1]) - np.log10(bins_avg[i]))
    #     print width
    #     #width = np.log10(width)
    #     y = piece_log(x, peak, width)
    #     plt.semilogx(x, y)
    #     plt.axvline(peak, color='lightgrey')
    plt.show()


def plot_bases():
    K = 10
    rpbins = np.logspace(-1, 1, K+1)
    x = np.logspace(-2, 2, 1000)
    pbases = piece_u(x, rpbins)
    tbases = tophat_u(x, rpbins)
    gbases = gauss_u(x, rpbins)
    trbases = trig_u(x, K/2)

    for bb in range(len(pbases)):
        # plt.loglog(x, pbases[bb], color='red', ls='-')
        # plt.loglog(x, tbases[bb], color='blue', ls='-')
        # plt.loglog(x, gbases[bb], color='green', ls='-')
        #plt.loglog(x, trbases[bb], color='orange', ls='-')

        #plt.semilogx(x, pbases[bb], color='red', ls='-')
        #plt.semilogx(x, tbases[bb], color='blue', ls='-')
        plt.semilogx(x, gbases[bb], color='green', ls='-')
        #plt.semilogx(x, trbases[bb], color='orange', ls='-')

    plt.ylim(10**-1, 2*10**0)

    plt.show()


def piece_cf():
    K = 10
    rpbins = np.logspace(-1, 1, K+1)
    rpbins_avg = 0.5 * (rpbins[1:] + rpbins[:-1])

    N = 5
    up = np.zeros(K)
    ut = np.zeros(K)
    ug = np.zeros(K)
    utr = np.zeros(N*2)

    #rps = [0.05, 0.1, 1]
    rps = get_rps(rpbins, 1)

    plt.figure()
    for rp in rps:
        ut += tophat_u(rp, rpbins)
        up += piece_u(rp, rpbins)
        ug += gauss_u(rp, rpbins)
        utr += trig_u(rp, N)

    x = np.logspace(-2, 2, 1000)
    pbases = piece_u(x, rpbins)
    tbases = tophat_u(x, rpbins)
    gbases = gauss_u(x, rpbins)
    trbases = trig_u(x, N)

    pcf_sum = np.zeros_like(x)
    tcf_sum = np.zeros_like(x)
    gcf_sum = np.zeros_like(x)
    trcf_sum = np.zeros_like(x)

    for bb in range(len(pbases)):
        pbase = np.array(pbases[bb])
        pcf = pbase*up[bb]
        pcf_sum += pcf
        plt.semilogx(x, pcf, color='red', alpha=0.3)

        tbase = np.array(tbases[bb])
        tcf = tbase * ut[bb]
        tcf_sum += tcf
        plt.semilogx(x, tcf, color='blue', alpha=0.3)

        gbase = np.array(gbases[bb])
        gcf = gbase * ug[bb]
        gcf_sum += gcf
        plt.semilogx(x, gcf, color='green', alpha=0.3)

    for bb in range(len(trbases)):
        trbase = np.array(trbases[bb])
        trcf = trbase * utr[bb]
        trcf_sum += trcf
        plt.semilogx(x, trcf, color='orange', alpha=0.3)

    plt.loglog(x, pcf_sum, color='red', ls='-')
    plt.loglog(x, tcf_sum, color='blue', ls='-')
    plt.loglog(x, gcf_sum, color='green', ls='-')
    plt.loglog(x, trcf_sum, color='orange', ls='-')

    plt.ylim(10**-1, 10**2)

    plt.show()


def scratch_continuous():
    K = 10
    Kbig = 50
    rpmax = 10
    rpbins = np.logspace(np.log10(0.1), np.log10(rpmax), K)
    rpbins_fine = np.logspace(np.log10(0.1), np.log10(rpmax), Kbig)
    rpbins_avg = 0.5 * (rpbins[1:] - rpbins[:-1])

    # rpbins = np.linspace(0.1, rpmax, K)
    # rpbins_fine = np.linspace(0.1, rpmax, K*10)

    up = np.zeros(Kbig)
    ut = np.zeros(Kbig)
    ut_few = np.zeros(K-1)
    up_few = np.zeros(K-1)
    ug_few = np.zeros(K)
    uq_few = np.zeros(K)

    #for rp in [2]:
    #rps = get_rps(rpbins_fine, rpmax)
    #rps = rpbins_fine[:(K*10)/2]
    #rps = rpbins_fine
    #rps = rpbins_fine
    rps = [0.05, 0.1, 1]
    #print max(rps)
    #print max(rpbins_fine)
    #plt.figure()
    #plt.hist(rps)

    plt.figure()
    #for rp in [0.2, 0.3, 0.4, 1, 2, 6]:
    for rp in rps:
        #func = gaussian(rp, rpbins, rpbins_fine)
        #func = piecewise(rp, rpbins, rpbins_fine)
        #plt.loglog(rpbins, piecewise_fine(rp, rpbins), color='purple')
        #u += piecewise_fine(rp, rpbins_fine)
        ut_few += tophat(rp, rpbins)

        up_few += piecewise_log(rp, rpbins)
        #ug_few += gaussian(rp, rpbins)
        #plt.loglog(rpbins, gaussian(rp, rpbins), color='red')

        #plt.loglog(rpbins, quadratic(rp, rpbins), color='black')

        #uq_few += quadratic(rp, rpbins)

        #ut += tophat(rp, rpbins_fine)

    print up_few
    print ut_few

    x = np.logspace(-2, 2, 1000)
    pbases = piece_log_bins(x, rpbins)
    tbases = tophat(x, rpbins)
    pcf_sum = np.zeros_like(x)
    tcf_sum = np.zeros_like(x)
    for bb in range(len(pbases)):
        pbase = np.array(pbases[bb])
        pcf = pbase*up_few[bb]
        pcf_sum += pcf
        plt.semilogx(x, pcf, color='blue')

        tbase = np.array(tbases[bb])
        tcf = pbase * ut_few[bb]
        tcf_sum += tcf
        plt.semilogx(x, tcf, color='orange')

    #plt.loglog(rpbins_avg, up_few, color='cyan', ls='--', marker='o')
    plt.loglog(x, pcf_sum, color='cyan', ls='--')

    #plt.loglog(rpbins_avg, ut_few, color='limegreen', ls=':', marker='o')
    plt.loglog(x, tcf_sum, color='magenta', ls=':', marker=None)

    #plt.loglog(rpbins, ug_few, color='orange', ls='-.', marker='o')
    #plt.loglog(rpbins, uq_few, color='grey', ls='-', marker='o')

    # plt.loglog(rpbins_fine, u, color='blue', marker='o')
    # plt.loglog(rpbins_fine, ut, color='green', ls='--', marker='o')


    # plt.plot(rpbins_fine, u, color='blue')
    # plt.plot(rpbins, ut, color='green')

    plt.ylim(10**-1, 10**3)
    plt.show()

def check(arg1, kw1=None):
    print arg1
    print kw1

def argcheck(arg1, a1, a2):
    print arg1
    print a1
    print a2

def elsewhere(arg1, **kwargs):
    check(arg1, **kwargs)

def elsewhere2(arg1, *args):
    argcheck(arg1, *args)

def layer(arg1, *args):
    elsewhere2(arg1, *args)

#layer('yo', [1,2], 4.5)
def twoarrs():
    a = np.array([1,2])
    b = np.array([3,4])
    return a, b

def scratch_args():
    # print twoarrs()
    # a, b = 2.*twoarrs()
    # print a
    # print b
    elsewhere2('yo', 'a1', 'a2')

if __name__=='__main__':
    main()