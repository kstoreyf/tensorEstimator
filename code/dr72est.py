#!/usr/bin/env python
import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM
from matplotlib import pyplot as plt
import argparse

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
          21: 60, 20: 60, 19:60, 18: 60, 17: 60, 16: 60, 15: 60, 14: 40, 13: 40, 0: 60, 'bins':60}
labels_mr = {7: r'$-23<M_r<-22$', 8: r'$-22<M_r<-21$', 9: r'$-21<M_r<-20$',
            10: r'$-20<M_r<-19$', 11: r'$-19<M_r<-18$', 12: r'$-18<M_r<-17$',
            21: r'$M_r<-22$', 20: r'$M_r<-21.5$', 19: r'$M_r<-21$', 18: r'$M_r<-20.5$',
            17: r'$M_r<-20$', 16: r'$M_r<-19.5$', 15: r'$M_r<-19$', 14: r'$M_r<-18.5$',
            13: r'$M_r<-18$'}
colors = {7:'red', 8:'orange', 9:'green',
          10:'blue', 11:'cyan', 12:'magenta'}




def main(samplenum):
    #samplenums = [13,14,15,16,17,18,19,20,21]
    #samplenums = [13, 21]
    samplenums = [7,8,9,10,11,12]
    #samples_noczcut(samplenums)
    #samples_czcut(samplenums)

    #run_dr72()
    #check_samples()
    #run_dr7_LRGs()
    #fix_results()
    #read_dat(0)
    #plot_survey()
    #make_testsquare(samplenums)
    #make_testsquare([0])
    #write_randoms(samplenums)
    run_est(samplenum)
    #run_est_grid(samplenum)
    #eval_amps()
    #plot_bins()


def read_dat(samplenum):

    print samplenum
    fn = '../data/lss.dr72bright{}.dat'.format(samplenum)
    df = pd.read_csv(fn, header=None, delim_whitespace=True, names=['indx',
                    'sector', 'mregion', 'ra', 'dec', 'cz', 'fgotten', 'selection_fn'])
    df['z'] = df['cz']/c_kms
    print len(df.index)

    fn_photo = '../data/photoinfo.dr72bright{}.dat'.format(samplenum)
    df_photo = pd.read_csv(fn_photo, header=None, delim_whitespace=True, names=['indx',
                    'M_u', 'M_g', 'M_r', 'M_i', 'M_z', 'mu_{50}', 'r50/r90'])

    df = pd.merge(df, df_photo, on='indx') #This doesn't lose any data

    #df = df[df['cz']>cz_lims[samplenum][0]][df['cz']<cz_lims[samplenum][1]]
    print min(df['cz']), max(df['cz'])
    print min(df['z']), max(df['z'])

    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    curr_zdep = q0 * (1.0 + q1 * (df['z'] - qz0))
    df['M_rz'] = df['M_r'] + curr_zdep * (df['z'] - qz0)
    print min(df['M_rz']), max(df['M_rz'])

    fn_save = '../data/lss.dr72bright{}_photo.dat'.format(samplenum)
    df.to_csv(fn_save)


def write_randoms(samplenums):
    #tag = '_photo'
    tag = '_noczcut'
    for samplenum in samplenums:
        fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
        df_data = pd.read_csv(fn)
        df_rand = read_random(df_data)
        fn_save = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
        df_rand.to_csv(fn_save)


# fixed for bad assignment of Mrz
def read_random(df_data):
    fn = '../data/random-0.dr72bright.dat'
    # pretty sure fgotten is correct but couldn't find anywhere
    df_rand = pd.read_csv(fn, header=None, delim_whitespace=True, names=['ra',
                            'dec', 'sector', 'mregion', 'fgotten', 'min_mag?'])

    nrand = len(df_rand.index)

    df_rand['cz'] = df_data['cz'].sample(n=nrand, replace=True).values
    df_rand['z'] = df_rand['cz'] / c_kms

    # TODO: will need M_rz for random when do my estimator
    df_rand['M_r'] = df_data['M_r'].sample(n=nrand, replace=True).values
    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    curr_zdep = q0 * (1.0 + q1 * (df_rand['z'] - qz0))
    df_rand['M_rz'] = df_rand['M_r'] + curr_zdep * (df_rand['z'] - qz0)
    return df_rand


def make_testsquare(samplenums):

    for samplenum in samplenums:
        print samplenum
        fn_data = '../data/lss.dr72bright{}'.format(samplenum)
        #fn_rand = '../data/random-0.dr72bright'
        fn_rand = '../data/random-0.dr72bright{}'.format(samplenum)

        fns = [fn_data, fn_rand]
        #tags = ['_photo', '_photo']
        tags = ['_noczcut', '_noczcut']
        for i in range(len(tags)):
            df = pd.read_csv(fns[i]+tags[i]+'.dat')

            ### square 1k
            # ramin = 150
            # ramax = 152.5
            # decmin = 10
            # decmax = 15
            ### square 5k
            ramin = 150
            ramax = 155
            decmin = 10
            decmax = 25
            df = df[df['ra'] > ramin][df['ra'] < ramax][df['dec'] > decmin][df['dec'] < decmax]
            print len(df)
            print min(df['cz']), max(df['cz'])

            #plt.figure()
            #plt.scatter(df['ra'], df['dec'], s=1, color='red')

            fn_save = fns[i]+'_square5k_noczcut.dat'
            df.to_csv(fn_save)
        #plt.show()



def check_samples():
    #samplenums = [7, 8, 9, 10, 11, 12]
    samplenums = [0]
    #tags = []
    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75)

    for samplenum in samplenums:

        #fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        fn = '../data/lss.dr72bright{}_photo.dat'.format(samplenum)
        data = pd.read_csv(fn)
        #print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])
        print 'ndata=', len(data.index)
        #Mr_h = data['M_r'] + 5*np.log10(cosmo.h)
        curr_zdep = q0 * (1.0 + q1 * (data['z'] - qz0))
        Mr_z = data['M_r'] + curr_zdep*(data['z']-qz0)
        print min(data['M_r']), max(data['M_r']), np.mean(data['M_r']), "|", \
            min(Mr_z), max(Mr_z), np.mean(Mr_z)
        # so the given data['M_r'] has already been evolution corrected i think?




def run_est(samplenum):

    nproc = 2
    frac = 1

    tag = '_square1k'
    basistag = 'toprobust'
    basisfuncs = [estimator_chunks.tophat_robust]
    Kfac = 1
    vals = None

    # tag = '_square5k'
    # basistag = 'topMrz0linsumnoabs'
    # basisfuncs = [estimator_chunks.top_Mrz]
    # Kfac = 2
    # vals = [-18.5, -19.5, -20.5, -21.5, -22.5]
         
    #tag = '_czcut'
    #basistag = 'top'
    #basisfuncs = [estimator_chunks.tophat]
    #Kfac = 1
    #vals = None

    if samplenum=='bins':
      bintag = '_bins'
    else:
      bintag = '_bin{}'.format(samplenum)
    fractag = '_frac'.format(frac)
    #saveto = '../results/amps/amps{}{}{}{}.npy'.format(bintag, tag, basistag, fractag)
    saveto = None
    #Separations should be given in Mpc/h
    min_sep = 0.13
    max_sep = 40. #Mpc/h
    #max_sep = 1. #Mpc/h
    bin_size = 0.2
    K = int((np.log10(max_sep) - np.log10(min_sep))/bin_size)
    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = run.bins_logavg(rpbins)
    logrpbins_avg = run.logbins_avg(rpbins)
    logwidth = run.log_width(rpbins)
    pimax = pimaxs[int(samplenum)]
    pibinwidth = int(pimax)

    bin_arg = np.log10(rpbins)
    #bin_arg = logrpbins_avg

    wp = True
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75)


    ###############
    print 'Load samples'
    fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
    data1 = pd.read_csv(fn)
    fn_rand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
    rand1 = pd.read_csv(fn_rand)
    data1 = run.add_info(data1)
    rand1 = run.add_info(rand1)

    data1 = data1.sample(frac=frac)
    #if samplenum>0:
    #  frac /= 10.
    rand1 = rand1.sample(frac=frac)
    #data1 = data1[:int(frac*len(data1.index))]
    #rand1 = rand1[:int(frac*len(rand1.index))]
    print len(data1), len(rand1)
    data2 = data1
    rand2 = rand1
    ##################

    #vals = [np.mean(data1['M_rz'])]
    #vals = np.linspace(min(data1['M_rz']), max(data1['M_rz']), 6)
    #vals = [-18, -18.5, -19, -19.5, -20]

    rps = []
    wprps = []
    # weights_data = data1['fgotten']
    # weights_rand = rand1['fgotten']
    # weights_data = None
    # weights_rand = None
    # print 'Corrfuncing'
    # rpcf, wprpcf = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo,
    #                                 weights_data=weights_data, weights_rand=weights_rand,
    #                                 pibinwidth=pibinwidth) #pimax because my est isnt binning in pi
    # rps.append(rpbins_avg)
    # wprps.append(wprpcf)

    K *= Kfac
    print 'Running vecest'
    #rp, wprp = run.run(data1, rand1, data2, rand2, pimax, min_sep, max_sep,
    #                  bin_size, basisfuncs, K, cosmo, wp, rpbins, vals, pibinwidth, bin_arg, logwidth)
    rpsest, wprpsest, a = run.run_chunks(data1, rand1, data2, rand2, pimax, min_sep,
                                max_sep, basisfuncs, K, cosmo, wp, rpbins,
                                vals, pibinwidth, nproc, bin_arg, logwidth)
    rps += [rpbins_avg]*len(wprpsest)
    wprps += wprpsest

    if saveto:
        print 'Saved to {}'.format(saveto)
        np.save(saveto, [a, K, rpbins, pibinwidth, bin_arg, logwidth])


    print rps
    print wprps
    #labels = ['corrfunc'] + ['bright0 square1k, M_rz={}'.format(val) for val in vals]

    #plotter.plot_wprp(rps, wprps, labels)

def run_est_grid(samplenum):

    nproc = 1
    frac = 0.01

    tag = '_square1k'
    #basistag = 'gridMrz0lin'
    #basisfuncs = [estimator_chunks.grid_Mrz]
    #basistag = 'matchdimmest'
    basistag = '_matchbrighter'
    basisfuncs = [estimator_chunks.match_bins]
    Kfac = 1
    vals = [-18.5, -19.5, -20.5, -21.5, -22.5]
    #savetag = '_binwidthdec'
    savetag = ''
    saveto = '../results/amps/amps{}{}{}_frac{}{}.npy'.format(samplenum, tag, basistag, frac, savetag)

    #Separations should be given in Mpc/h
    min_sep = 0.13
    max_sep = 40. #Mpc/h
    K = 100
    #bin_size = 0.2
    #K = int((np.log10(max_sep) - np.log10(min_sep))/bin_size)
    print 'K:', K
    print 'Kfac:', Kfac
    print 'frac:', frac

    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    print 'rpbins:',rpbins
    rpbins_avg = run.bins_logavg(rpbins)
    logrpbins_avg = run.logbins_avg(rpbins)
    logwidth = run.log_width(rpbins)
    pimax = pimaxs[samplenum]
    pibinwidth = int(pimax)

    bin_arg = np.log10(rpbins)

    wp = True
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75)

    ###############
    print 'Load samples'
    fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
    data1 = pd.read_csv(fn)
    fn_rand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
    rand1 = pd.read_csv(fn_rand)
    data1 = run.add_info(data1)
    rand1 = run.add_info(rand1)

    # data1 = data1.sample(frac=frac)
    # rand1 = rand1.sample(frac=frac)
    data1 = data1[:int(frac*len(data1.index))]
    rand1 = rand1[:int(frac*len(rand1.index))]
    print len(data1), len(rand1)
    data2 = data1
    rand2 = rand1
    ##################

    rps = []
    wprps = []

    K *= Kfac
    print 'Running vecest'
    rpsest, wprpsest, a = run.run_chunks(data1, rand1, data2, rand2, pimax, min_sep, max_sep, basisfuncs, K, cosmo, wp, rpbins, vals, pibinwidth, nproc, bin_arg, logwidth)
    rps += [rpbins_avg] * len(wprpsest)
    wprps += wprpsest

    print 'Saving to {}'.format(saveto)
    np.save(saveto, [a, K, rpbins, pibinwidth, bin_arg, logwidth])

    print 'Rp, wprps, a'
    print rps
    print wprps
    print a



def eval_amps():

    # evalat = '_many'
    # tag = '_square5k'
    # basistag = '_topMrz0linsumnoabs'
    # #basistag = '_toplogMrz0linmean'
    # #basistag = 'topMrzallmean'
    # basisfuncs = [estimator_chunks.top_Mrz]
    # fractag = '_frac1'
    # #extratag = '_valx0.5'
    # #extratag = '_doublecheck'
    # savetag = '
    # '
    # evalat = '_midpoints'
    # samplenum = 'bins'
    # tag = '_square5k'
    # fractag = '_frac1'
    # basistag = '_grid'
    # basisfuncs = [estimator_chunks.grid_Mrz]
    # extratag = ''
    # savetag = extratag

    samplenum = 11
    evalat = ''
    #tag = '_bin{}_czcut'.format(samplenum)
    tag = '_square5k'
    basistag = '_toprobust'
    basisfuncs = [estimator_chunks.tophat_robust]
    #fractag = '_frac1rand0.1'
    fractag = '_frac'
    savetag = ''
    extratag = ''

    # samplenum = 7
    # evalat = ''
    # tag = '_bin{}_square1k'.format(samplenum)
    # samplenum = ''
    # basistag = '_top'
    # basisfuncs = [estimator_chunks.tophat]
    # fractag = '_frac1'
    #extratag = ''
    # savetag = ''

    # evalat = '_midpoints'
    # tag = '_square5k'
    # fractag = '_frac1'
    # #basistag = '_matchbins'
    # #basistag = '_matchbinsandrp2'
    # #basistag = '_matchbrighter'
    # #basistag = '_matchdimmest'
    # basistag = '_matchbins'
    # basisfuncs = [estimator_chunks.match_bins]
    # #extratag = '_binwidthdec'
    # extratag = ''
    # savetag = extratag
    # samplenum = 0

    if samplenum=='bins':
      bintag = '_bins'
    else:
      bintag = '_bin{}'.format(samplenum)

    fn = '../results/amps/amps{}{}{}{}{}.npy'.format(bintag, tag, basistag, fractag, extratag)
    saveto = 'plots_2019-02-06/wprp{}{}{}{}{}.png'.format(bintag, tag, basistag, evalat, savetag)
    savetowp = '../results/wprps/wprp{}{}{}{}{}.png'.format(bintag, tag, basistag, evalat, savetag)
    #basisfuncs = [estimator_chunks.tophat]
    #basisfuncs = [estimator_chunks.top_rand]

    if evalat=='_midpoints':
        vals = [-18.5, -19.5, -20.5, -21.5, -22.5]

        #vals = [-18.3, -18.75, -19.5, -21]
        #vals = [-18.5, -20, -21.5]


        labels = ['M_rz={} '.format(val) for val in vals]
    elif evalat=='_many':
        vals = np.arange(-18.5, -22.6, -0.1)
        labelarr = [-18.5, -19.5, -20.5, -21.5, -22.5]
        labels = ['M_rz={} '.format(val) if round(val,1) in labelarr else None for val in vals]
    elif evalat=='':
        vals = None
        labels = ['bin{}'.format(samplenum)]



    a, K, rpbins, pibinwidth, bin_arg, logwidth = np.load(fn)
    x = run.bins_logavg(rpbins)
    print ':::a:::'
    print a

    rps, wprps = run.calc_wprp(a, x, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)

    print len(rps)
    print len(labels)
    #print rps
    #print wprps

    np.save(savetowp ,[rps, wprps, labels])
    #fn = 'plots_2019-01-29/wprp_Mrzlin_square5k_midpoints.png'
    #fn = 'plots_2019-01-29/wprp_Mrzall_square5k_-22.5to-18.5.png'
    #plotter.plot_wprp(rps, wprps, labels, saveto=saveto)


def eval_amps_grid():

    evalat = '_many'
    tag = '_square1k'
    basistag = '_topMrz0lingrid'
    basisfuncs = [estimator_chunks.grid_Mrz]
    fractag = '_frac1'
    extratag = ''

    fn = '../results/amps/amps{}{}{}.npy'.format(tag, basistag, fractag)
    saveto = 'plots_2019-02-06/wprp{}{}{}{}.png'.format(tag, basistag, evalat, extratag)


    if evalat=='_midpoints':
        vals = [-18.5, -19.5, -20.5, -21.5, -22.5]
        labels = ['M_rz={} '.format(val) for val in vals]
    elif evalat=='_many':
        vals = np.arange(-18.5, -22.6, -0.1)
        labelarr = [-18.5, -19.5, -20.5, -21.5, -22.5]
        labels = ['M_rz={} '.format(val) if round(val,1) in labelarr else None for val in vals]
    elif evalat=='':
        vals = None
        labels = ['bin{}'.format(samplenum)]

    a, K, rpbins, pibinwidth, bin_arg, logwidth = np.load(fn)
    x = run.bins_logavg(rpbins)
    print ':::a:::'
    print a

    rps, wprps = run.calc_wprp(a, x, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)

    print len(rps)
    print len(labels)
    print rps
    print wprps

    np.save('../results/wprps/wprp{}{}{}{}'.format(tag, basistag, fractag, evalat),
            [rps, wprps, labels])
    #fn = 'plots_2019-01-29/wprp_Mrzlin_square5k_midpoints.png'
    #fn = 'plots_2019-01-29/wprp_Mrzall_square5k_-22.5to-18.5.png'
    plotter.plot_wprp(rps, wprps, labels, saveto=saveto)



def plot_bins():

    #samplenums = [7,8,9,10,11]
    #samplenums = [11, 10, 9, 8, 7]
    samplenums = [11,10,9,8,7]
    #samplenums = [7]
    rps = []
    wprps = []
    labels = []
    frac = 1
    tag = '_square5k_noczcut'
    #tag = '_square5k'
    for samplenum in samplenums:
        vals = None
        fn = '../results/amps/amps_bin{}{}_top_frac{}.npy'.format(samplenum, tag, frac)
        basisfuncs = [estimator_chunks.tophat]

        a, K, rpbins, pibinwidth, bin_arg, logwidth = np.load(fn)
        x = run.bins_logavg(rpbins)

        rp, wprp = run.calc_wprp(a, x, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)
        rps += rp
        wprps += wprp
        labels.append('bin {}, {}'.format(samplenum, labels_mr[samplenum]))
    fn = 'plots_2019-01-29/wprp_bins{}.png'.format(tag)
    plotter.plot_wprp(rps, wprps, labels, saveto=fn)


def run_sample_corrfunc(samplenum, tag, min_sep, max_sep, rpbins, pimax, wp, cosmo,
                        frac=1, bin_size=None, pibinwidth=2):

    fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
    data1 = pd.read_csv(fn)
    fn_rand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
    rand1 = pd.read_csv(fn_rand)
    #rand1 = get_random(data1)

    print 'Sample {}'.format(samplenum)

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data1 = data1.sample(frac=frac)
    rand1 = rand1.sample(frac=frac)
    # data1 = data1[:int(frac*len(data1.index))]
    # rand1 = rand1[:int(frac*len(rand1.index))]

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data2 = data1
    rand2 = rand1

    start = time.time()
    weights_data = data1['fgotten']
    weights_rand = rand1['fgotten']
    rp_avg, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo,
                                    weights_data=weights_data, weights_rand=weights_rand,
                                    pibinwidth=pibinwidth)

    end = time.time()
    print 'Time for sample {}, ndata={}, nrand={}: {}'.format(
        samplenum, len(data1.index), len(rand1.index), end-start)
    return rp_avg, wprp


def get_random(df_data):
    fn = '../data/random-0.dr72bright.dat'
    # pretty sure fgotten is correct but couldn't find anywhere
    df_rand = pd.read_csv(fn, header=None, delim_whitespace=True, names=['ra',
            'dec', 'sector', 'mregion', 'fgotten', 'min_mag?'])

    nrand = len(df_rand.index)

    df_rand['cz'] = df_data['cz'].sample(n=nrand, replace=True).values
    df_rand['z'] = df_rand['cz'] / c_kms

    # TODO: will need M_rz for random when do my estimator
    # df_rand['M_r'] = np.random.random(nrand)*(Mr_lims[samplenum][1]-Mr_lims[samplenum][0]) \
    #            + Mr_lims[samplenum][0]
    return df_rand


def samples_noczcut(samplenums):

    for samplenum in samplenums:
        print samplenum, labels_mr[samplenum]
        fn = '../data/lss.dr72bright{}.dat'.format(samplenum)
        df = pd.read_csv(fn, header=None, delim_whitespace=True, names=['indx',
                        'sector', 'mregion', 'ra', 'dec', 'cz', 'fgotten', 'selection_fn'])
        df['z'] = df['cz']/c_kms
        print len(df.index)

        fn_photo = '../data/photoinfo.dr72bright{}.dat'.format(samplenum)
        df_photo = pd.read_csv(fn_photo, header=None, delim_whitespace=True, names=['indx',
                        'M_u', 'M_g', 'M_r', 'M_i', 'M_z', 'mu_{50}', 'r50/r90'])

        df = pd.merge(df, df_photo, on='indx') #This doesn't lose any data

        print min(df['cz']), max(df['cz'])
        q0 = 2.0
        q1 = -1.0
        qz0 = 0.1
        curr_zdep = q0 * (1.0 + q1 * (df['z'] - qz0))
        df['M_rz'] = df['M_r'] + curr_zdep * (df['z'] - qz0)
        print min(df['M_rz']), max(df['M_rz'])

        print len(df.index)

        fn_save = '../data/lss.dr72bright{}_noczcut.dat'.format(samplenum)
        df.to_csv(fn_save)



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("samplenum", help="number of sample to process")
    args = parser.parse_args()
    return args


if __name__=="__main__":
    arguments = parse_arguments()
    main(arguments.samplenum)
