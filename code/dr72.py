import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM

import run
import plotter
import estimator


cz_lims = {7: [30900, 73500], 8: [19900, 47650], 9: [12600, 31900],
           10: [8050, 19250], 11: [5200, 12500], 12: [3200, 7850]}
Mr_lims =  {7: [-23, -22], 8: [-22, -21], 9: [-21, -20],
           10: [-20, -19], 11: [-19, -18], 12: [-18, -17]}
pimaxs = {7: 60, 8: 60, 9: 60, 10: 60, 11: 40, 12: 40}
labels_mr = {7: r'$-23<M_r<-22$', 8: r'$-22<M_r<-21$', 9: r'$-21<M_r<-20$',
            10: r'$-20<M_r<-19$', 11: r'$-19<M_r<-18$', 12: r'$-18<M_r<-17$'}
colors = {7:'red', 8:'orange', 9:'green',
          10:'blue', 11:'cyan', 12:'magenta'}
c_kms = 3.e5  # c in km/s

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)


def main():
    #samples_czcut()
    run_dr72()
    #check_samples()

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


def check_samples():
    samplenums = [7, 8, 9, 10, 11, 12]
    print cosmo.h
    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    for samplenum in samplenums:
        fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        data = pd.read_csv(fn)
        print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])
        print 'ndata=', len(data.index)
        Mr_h = data['M_r'] + 5*np.log10(cosmo.h)
        curr_zdep = q0 * (1.0 + q1 * (data['z'] - qz0))
        Mr_z = data['M_r'] + curr_zdep*(data['z']-qz0)
        print min(data['M_r']), max(data['M_r']), np.mean(data['M_r']), "|", \
            min(Mr_z), max(Mr_z), np.mean(Mr_z)
        # so the given data['M_r'] has already been evolution corrected i think?
        #min(Mr_h), max(Mr_h), np.mean(Mr_h)

def run_bins(min_sep, max_sep, bin_size, wp):
    #samplenums = [7, 8, 9, 10, 11, 12]
    #samplenums = [9,10,11]
    samplenums = [7]

    #samplenums = []

    K = (np.log10(max_sep) - np.log10(min_sep))/bin_size
    print K
    rpbins = np.logspace(np.log(min_sep), np.log(max_sep), K+1, base=np.e)
    print rpbins
    rpbins_avg = 0.5 * (rpbins[1:] + rpbins[:-1])
    rps = [rpbins_avg]*len(samplenums)
    pibinwidth = 2 #Mpc/h
    wprps = []
    labels = []
    cols = []
    for samplenum in samplenums:
        if samplenum in labels_mr:
            labels.append(labels_mr[samplenum])
        else:
            labels.append(samplenum)
        pimax = pimaxs[samplenum]
        xi = run_sample(samplenum, min_sep, max_sep, rpbins, pimax, wp, pibinwidth=pibinwidth)
        wprps.append(xi) #*2?
        cols.append(colors[samplenum])

    #plotter.plot_wprp(rps, wprps, labels, colors=cols)


def run_together(min_sep, max_sep, bin_size, K, pimax, wp):
    #samplenums = [7, 8, 9, 10, 11, 12]
    samplenums = [8, 9, 10]

    data1, rand1 = combine_samples(samplenums)
    data2 = data1.copy()
    rand2 = rand1.copy()
    print 'ndata={}, nrand={}'.format(len(data1.index), len(rand1.index))

    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)
    rpbins_avg = 0.5 * (rpbins[1:] + rpbins[:-1])
    logrpbins_avg = np.log10(rpbins_avg)
    logwidth = np.log10(rpbins_avg[1]) - np.log10(rpbins_avg[0])

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


def run_sample(samplenum, min_sep, max_sep, rpbins, pimax, wp, bin_size=None, pibinwidth=pibinwidth):

    fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
    data1 = pd.read_csv(fn)
    rand1 = get_random(samplenum)

    print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    frac = 0.08
    data1 = data1.sample(frac=frac)
    rand1 = rand1.sample(frac=frac)

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data2 = data1
    rand2 = rand1

    start = time.time()
    #xi, dd, dr, rd, rr = run.run_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    #xi, dd, dr, rd, rr = run.run_treecorr_orig(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    est_ls, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax, pibinwidth=pibinwidth)

    end = time.time()
    print 'Time for sample {}, ndata={}: {}'.format(
        samplenum, len(data1.index), end-start)
    return wprp


def get_random(samplenum):
    fn = '../data/random-0.dr72bright.dat'
    df_rand = pd.read_csv(fn, header=None, delim_whitespace=True, names=['ra',
            'dec', 'sector', 'mregion', '?', 'min_mag?'])

    nrand = len(df_rand.index)
    df_rand['cz'] = np.random.random(nrand)*(cz_lims[samplenum][1]-cz_lims[samplenum][0]) \
               + cz_lims[samplenum][0]
    df_rand['z'] = df_rand['cz'] / c_kms

    df_rand['M_r'] = np.random.random(nrand)*(Mr_lims[samplenum][1]-Mr_lims[samplenum][0]) \
               + Mr_lims[samplenum][0]
    #TODO: adjust for log10(h) thing
    #df_rand['M_r'] += np.log10 ...
    return df_rand

def samples_czcut():

    for samplenum in cz_lims:

        fn = '../data/lss.dr72bright{}.dat'.format(samplenum)
        df = pd.read_csv(fn, header=None, delim_whitespace=True, names=['indx',
                        'sector', 'mregion', 'ra', 'dec', 'cz', 'fgotten', 'selection_fn'])
        df['z'] = df['cz']/c_kms

        fn_photo = '../data/photoinfo.dr72bright{}.dat'.format(samplenum)
        df_photo = pd.read_csv(fn_photo, header=None, delim_whitespace=True, names=['indx',
                        'M_u', 'M_g', 'M_r', 'M_i', 'M_z', 'mu_{50}', 'r50/r90'])

        df = pd.merge(df, df_photo, on='indx')
        print len(df.index)

        df = df[df['cz']>cz_lims[samplenum][0]][df['cz']<cz_lims[samplenum][1]]
        print len(df.index)

        fn_save = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
        df.to_csv(fn_save)



if __name__=="__main__":
    main()