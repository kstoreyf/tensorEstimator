#!/usr/bin/env python
import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM
from matplotlib import pyplot as plt

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
    #samplenums = [13,14,15,16,17,18,19,20,21]
    #samplenums = [13, 21]
    #samples_czcut(samplenums, threshold=True)
    #run_dr72()
    #check_samples()
    #run_dr7_LRGs()
    #fix_results()
    #read_dat(0)
    #plot_survey()
    #make_testsquare([0])
    #write_random()
    run_est()
    #eval_amps()



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


def write_random():
    samplenum = 0
    tag = '_photo'
    fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
    df_data = pd.read_csv(fn)
    df_rand = read_random(df_data)
    fn_save = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
    df_rand.to_csv(fn_save)


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
        fn_data = '../data/lss.dr72bright{}'.format(samplenum)
        #fn_rand = '../data/random-0.dr72bright'
        fn_rand = '../data/random-0.dr72bright{}'.format(samplenum)

        fns = [fn_data, fn_rand]
        tags = ['_photo', '_photo']
        for i in range(len(tags)):
            df = pd.read_csv(fns[i]+tags[i]+'.dat')

            ramin = 150
            ramax = 152.5
            decmin = 10
            decmax = 15
            df = df[df['ra'] > ramin][df['ra'] < ramax][df['dec'] > decmin][df['dec'] < decmax]
            print len(df)

            fn_save = fns[i]+'_square1k.dat'
            df.to_csv(fn_save)


def plot_survey():

    samplenums = [0, 21]
    tags = ['_photo', '_czcut']

    for i in range(len(samplenums)):
        fn = '../data/lss.dr72bright{}{}.dat'.format(samplenums[i], tags[i])
        data = pd.read_csv(fn)

        frac = 0.01

        #data = data.sample(1000).reset_index(drop=True)
        #data = data[data['ra']<160][data['ra']>150][data['dec']>10][data['dec']<15]
        data = data[data['ra']<152.5][data['ra']>150][data['dec']>10][data['dec']<15]

        print len(data)
        #data = data.sample(1000).reset_index(drop=True)

        plt.figure()
        plt.scatter(data['ra'], data['dec'], s=1, color='red')
    plt.show()


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




def run_est():
    #Separations should be given in Mpc/h
    min_sep = 0.13
    max_sep = 40. #Mpc/h
    bin_size = 0.2
    K = int((np.log10(max_sep) - np.log10(min_sep))/bin_size)
    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = run.bins_logavg(rpbins)
    logrpbins_avg = run.logbins_avg(rpbins)
    logwidth = run.log_width(rpbins)
    pimax = 60.
    pibinwidth = int(pimax)

    wp = True
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75)

    samplenum = 0
    tag = '_square1k'

    ###############
    print 'Load samples'
    fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
    data1 = pd.read_csv(fn)
    fn_rand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
    rand1 = pd.read_csv(fn_rand)
    data1 = run.add_info(data1)
    rand1 = run.add_info(rand1)

    frac = 1
    #data1 = data1.sample(frac=frac)
    #rand1 = rand1.sample(frac=frac)
    data1 = data1[:int(frac*len(data1.index))]
    rand1 = rand1[:int(frac*len(rand1.index))]
    print len(data1), len(rand1)

    data2 = data1
    rand2 = rand1
    ##################

    #vals = [np.mean(data1['M_rz'])]
    #vals = np.linspace(min(data1['M_rz']), max(data1['M_rz']), 6)
    #vals = [-18, -18.5, -19, -19.5, -20]
    vals = [-18,-19,-20,-21,-22,-23]
    print vals
    bin_arg = logrpbins_avg

    rps = []
    wprps = []
    # weights_data = data1['fgotten']
    # weights_rand = rand1['fgotten']
    weights_data = None
    weights_rand = None
    print 'Corrfuncing'
    rpcf, wprpcf = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo,
                                    weights_data=weights_data, weights_rand=weights_rand,
                                    pibinwidth=pibinwidth) #pimax because my est isnt binning in pi
    rps.append(rpbins_avg)
    wprps.append(wprpcf)

    #basisfuncs = [estimator_chunks.tophat]

    #basisfuncs = [estimator_chunks.top_Mrz]
    basisfuncs = [estimator_chunks.top_Mrz_0lin]
    K *= 2
    print 'Running vecest'
    #rp, wprp = run.run(data1, rand1, data2, rand2, pimax, min_sep, max_sep,
    #                  bin_size, basisfuncs, K, cosmo, wp, rpbins, vals, pibinwidth, bin_arg, logwidth)
    nproc = 16
    rpsest, wprpsest, a = run.run_chunks(data1, rand1, data2, rand2, pimax, min_sep,
                                max_sep, bin_size, basisfuncs, K, cosmo, wp, rpbins,
                                vals, pibinwidth, nproc, bin_arg, logwidth)
    rps += [rpbins_avg]*len(wprpsest)
    wprps += wprpsest

    np.save('../results/amps/amps_square1k_topMrz0lin.npy', [a, K, rpbins, pibinwidth, bin_arg, logwidth])

    print rps
    print wprps
    labels = ['corrfunc'] + ['bright0 square1k, M_rz={}'.format(val) for val in vals]

    #plotter.plot_wprp(rps, wprps, labels)


def eval_amps():
    vals = [-18,-19,-20,-21,-22,-23]

    fn = '../results/amps/amps_square1k_topMrzlin.npy'
    basisfuncs = [estimator_chunks.top_Mrz_lin]
    a, K, rpbins, pibinwidth, bin_arg, logwidth = np.load(fn)
    x = run.bins_logavg(rpbins)

    rps, wprps = run.calc_wprp(a, x, basisfuncs, K, rpbins, vals, pibinwidth, bin_arg, logwidth)
    labels = ['bright0 square1k, M_rz={}'.format(val) for val in vals]
    plotter.plot_wprp(rps, wprps, labels)


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






if __name__=="__main__":
    main()
