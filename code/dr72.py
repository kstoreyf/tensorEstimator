#!/usr/bin/env python
import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM

import run
import plotter
import estimator


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
    run_dr72()
    #check_samples()
    #run_dr7_LRGs()
    #fix_results()

def run_dr72():

    #Separations should be given in Mpc/h
    min_sep = 0.13
    max_sep = 40. #Mpc/h
    #K = 12
    #bin_size_treecorr = np.log(max_sep / min_sep) / float(K)
    bin_size = 0.2

    wp = True
    print 'Running'
    run_bins(min_sep, max_sep, bin_size, wp)
    #run_together(min_sep, max_sep, bin_size, K, pimax, wp)


def run_dr7_LRGs():

    #sample = 'Full'
    #sample = 'Bright-no'
    sample = 'Dim-no'
    datafn = '../data/DR7-{}.ascii'.format(sample)
    randfn = '../data/random-DR7-{}.ascii'.format(sample)
    data = pd.read_table(datafn, index_col=False, delim_whitespace=True, names=['ra', 'dec', 'z',
            'M_g', 'sector_completeness', 'n(z)*1e4', 'radial_weight', 'fiber_coll_weight',
            'fogtmain', 'ilss', 'icomb', 'sector'], dtype={'z':np.float64}, skiprows=1)
    rand = pd.read_table(randfn, index_col=False,  delim_whitespace=True, names=['ra', 'dec', 'z',
            'sector_completeness', 'n(z)*1e4', 'radial_weight', 'ilss', 'sector'], dtype={'z':np.float64},
             skiprows=1)

    frac = 1
    #saveto = None
    saveto = "../results/wp_dr7_{}LRG_frac{}_weights.npy".format(sample, frac)
    cosmo = LambdaCDM(H0=70, Om0=0.25,Ode0=0.75)

    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)
    #Sector completeness already cut to >0.6, not sure if still have to downsample randoms
    #and many have sector completness > 1!! ??
    # data = data[data['z']<0.36]
    # data = data[data['ra']>90][data['ra']<270] #NGC

    print len(data.index)

    data = data.sample(frac=frac)
    rand = rand.sample(frac=frac)
    # data1 = data1[:int(frac*len(data1.index))]
    # rand1 = rand1[:int(frac*len(rand1.index))]
    print 'ndata=', len(data.index)
    print 'nrand=', len(rand.index)

    weights_data = data['radial_weight']*data['fiber_coll_weight']
    weights_rand = rand['radial_weight']

    #losmax = 1.0
    losmax = 40.0
    zspace = False
    if sample=='Bright-no':
        K = 21
        rmin = 60
        rmax = 200
    elif sample=='Full':
        K = 14
        rmin = 40
        rmax = 180
    elif sample=='Dim-no':
        K = 15
        rmin = 0.01
        rmax = 8.
    else:
        exit('ERROR')
    #bins = np.linspace(rmin, rmax, K + 1)
    bins = np.logspace(np.log10(rmin), np.log10(rmax), K + 1)

    start = time.time()
    # or rp, wp
    s, xi = run.run_corrfunc(data, rand, data, rand, bins, losmax, cosmo,
                             weights_data=weights_data, weights_rand=weights_rand, zspace=zspace)
    end = time.time()
    print 'Time for dr7 {} LRGs, ndata={}: {}'.format(sample, len(data.index), end-start)

    ss = [s]
    xis = [xi]
    labels = ['dr7 {} LRGs'.format(sample)]
    if saveto:
        run.save_results(saveto, ss, xis, labels)
    #plotter.plot_xi_zspace(ss, xis, labels)



def check_samples():
    samplenums = [7, 8, 9, 10, 11, 12]
    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75)

    for samplenum in samplenums:
        fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        data = pd.read_csv(fn)
        print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])
        print 'ndata=', len(data.index)
        #Mr_h = data['M_r'] + 5*np.log10(cosmo.h)
        curr_zdep = q0 * (1.0 + q1 * (data['z'] - qz0))
        Mr_z = data['M_r'] + curr_zdep*(data['z']-qz0)
        print min(data['M_r']), max(data['M_r']), np.mean(data['M_r']), "|", \
            min(Mr_z), max(Mr_z), np.mean(Mr_z)
        # so the given data['M_r'] has already been evolution corrected i think?


def run_bins(min_sep, max_sep, bin_size, wp, saveto=None):
    #samplenums = [7, 8, 9, 10, 11, 12]
    #samplenums = [14, 16, 18, 19, 20, 21]
    samplenums = [10]

    frac = 1
    saveto = "../results/dr72_bin{}_frac{}_pibinmaxweights1rand0.1n1.npy".format(samplenums[0], frac)
    print saveto

    cosmo = LambdaCDM(H0=70, Om0=0.25, Ob0=0.045, Ode0=0.75) 
    K = (np.log10(max_sep) - np.log10(min_sep))/bin_size
    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = run.bins_logavg(rpbins)


    rps = []
    #pibinwidth = 1 #Mpc/h
    wprps = []
    labels = []
    cols = []
    for samplenum in samplenums:
        if samplenum in labels_mr:
            labels.append(labels_mr[samplenum])
        else:
            labels.append(samplenum)
        pimax = pimaxs[samplenum]
        pibinwidth = pimax
        rp_avg, wprp = run_sample(samplenum, min_sep, max_sep, rpbins, pimax, wp, cosmo,
                                  frac=frac, pibinwidth=pibinwidth)
        print samplenum
        print wprp
        rps.append(rp_avg)
        wprps.append(wprp)
        if samplenum in colors:
            cols.append(colors[samplenum])

    if saveto:
        run.save_results(saveto, rps, wprps, labels)

    #plotter.plot_wprp(rps, wprps, labels)
    #plotter.plot_wprp(rps, wprps, labels, colors=cols)


def run_together(min_sep, max_sep, bin_size, K, pimax, wp):
    #samplenums = [7, 8, 9, 10, 11, 12]
    samplenums = [8, 9, 10]

    data1, rand1 = combine_samples(samplenums)
    data2 = data1.copy()
    rand2 = rand1.copy()
    print 'ndata={}, nrand={}'.format(len(data1.index), len(rand1.index))

    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = run.bins_logavg(rpbins)
    logwidth = run.log_width(rpbins)

    #basisfuncs = [estimator.top_z]
    #basisfuncs = [estimator.top_Mr]
    basisfuncs = [estimator.gauss_Mr]
    #basisfuncs = [estimator.tophat]

    K *= 3
    #vals = [0.1, 0.15, 0.2, 0.25]
    vals = [-22.5, -21.5, -20.5, -19.5, -18.5]
    #vals = [-21., -20.5, -20, -19.5]

    #labels = ["M_r={:.2f}".format(val) for val in vals]
    labels = ['corrfunc all']
    #labels = [0.15]
    #vals = [None]
    #labels = ['top']
    cols = ['purple', 'red', 'orange', 'green', 'blue', 'cyan', 'magenta', 'grey']


    print 'Run'
    start = time.time()
    est_ls, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax)
    rps = [rpbins_avg]
    wprps = [wprp]
    # rps, wprps = run.run(data1, rand1, data2, rand2, pimax, min_sep, max_sep, bin_size, basisfuncs,
    #     K, cosmo, wp, rpbins, vals, logrpbins_avg, logwidth)
    end = time.time()

    print 'Time for all, ndata={}: {}'.format(len(data1.index), end-start)

    plotter.plot_wprp(rps, wprps, labels, colors=cols)



def combine_samples(samplenums):
    datas = []
    rands = []
    for samplenum in samplenums:
        fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        data = pd.read_csv(fn, index_col=0)
        rand = get_random(samplenum)

        datas.append(data)
        rands.append(rand)
    datadf = pd.concat(datas, ignore_index=True)
    randdf = pd.concat(rands, ignore_index=True)

    frac = 0.2
    datadf = datadf.sample(frac=frac).reset_index(drop=True)
    randdf = randdf.sample(frac=frac).reset_index(drop=True)
    print 'Nums'
    print len(datadf.index)
    print len(randdf.index)
    #datadf = run.add_info(datadf, zfile=None)
    #randdf = run.add_info(randdf, zfile=None)

    return datadf, randdf

def combine_bins(samplenums):
    datas = []
    rands = []
    tag = '_square5k_noczcut'
    for samplenum in samplenums:
        print samplenum
        fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
        data = pd.read_csv(fn, index_col=0)
        #rand = get_random(samplenum)
        fnrand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)
        rand = pd.read_csv(fnrand, index_col=0)
        print len(data)
        print len(rand)
        datas.append(data)
        rands.append(rand)
    datadf = pd.concat(datas, ignore_index=True)
    randdf = pd.concat(rands, ignore_index=True)

    frac = 1
    datadf = datadf.sample(frac=frac).reset_index(drop=True)
    randdf = randdf.sample(frac=frac).reset_index(drop=True)

    datadf['M_rz'] = calc_Mrz(datadf['M_r'], datadf['z'])

    print 'Nums'
    print len(datadf.index)
    print len(randdf.index)
    #datadf = run.add_info(datadf, zfile=None)
    #randdf = run.add_info(randdf, zfile=None)

    fn_save = '../data/lss.dr72bright{}{}.dat'.format('bins', tag)
    fnrand_save = '../data/random-0.dr72bright{}{}.dat'.format('bins', tag)

    datadf.to_csv(fn_save)
    randdf.to_csv(fnrand_save)


def run_sample(samplenum, min_sep, max_sep, rpbins, pimax, wp, cosmo, frac=1, bin_size=None, pibinwidth=2):

    fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
    data1 = pd.read_csv(fn)
    rand1 = get_random(data1)

    print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data1 = data1.sample(frac=frac)
    print 'subsampling rands'
    frac /= 10.
    rand1 = rand1.sample(frac=frac)
    # data1 = data1[:int(frac*len(data1.index))]
    # rand1 = rand1[:int(frac*len(rand1.index))]

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data2 = data1
    rand2 = rand1

    start = time.time()
    #xi, dd, dr, rd, rr = run.run_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    #xi, dd, dr, rd, rr = run.run_treecorr_orig(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    #weights_data = data1['fgotten']
    #weights_rand = rand1['fgotten']
    weights_data = None
    weights_rand = None
    rp_avg, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, cosmo,
                                    weights_data=weights_data, weights_rand=weights_rand,
                                    pibinwidth=pibinwidth)

    end = time.time()
    print 'Time for sample {}, ndata={}: {}'.format(
        samplenum, len(data1.index), end-start)
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


def samples_czcut(samplenums):

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

        df = df[df['cz']>cz_lims[samplenum][0]][df['cz']<cz_lims[samplenum][1]]
        print min(df['cz']), max(df['cz'])
        q0 = 2.0
        q1 = -1.0
        qz0 = 0.1
        curr_zdep = q0 * (1.0 + q1 * (df['z'] - qz0))
        df['M_rz'] = df['M_r'] + curr_zdep * (df['z'] - qz0)
        print min(df['M_rz']), max(df['M_rz'])

        print len(df.index)

        fn_save = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        df.to_csv(fn_save)



def fix_results():
    samplenums = [7,8,9,10,11,12]
    for samplenum in samplenums:
        fn = '../results/dr72_bin{}_all.npy'.format(samplenum)
        rps, wprps, labels = run.load_results(fn)
        wprp_fixed = wprps[0]*2
        wprps_fixed = [wprp_fixed]
        #run.save_results(fn, rps, wprps_fixed, labels)


if __name__=="__main__":
    main()
