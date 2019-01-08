import pandas as pd
import numpy as np
import time

import run
import plotter


cz_lims = {7: [30900, 73500], 8: [19900, 47650], 9: [12600, 31900],
           10: [8050, 19250], 11: [5200, 12500], 12: [3200, 7850]}
labels_mr = {7: r'$-23<M_r<-22$', 8: r'$-22<M_r<-21$', 9: r'$-21<M_r<-20$',
            10: r'$-20<M_r<-19$', 11: r'$-19<M_r<-18$', 12: r'$-18<M_r<-17$'}
colors = {7:'red', 8:'orange', 9:'yellow',
          10:'green', 11:'blue', 12:'purple'}
c_kms = 3.e5  # c in km/s


def main():
    #samples_czcut()
    run_dr72()

def run_dr72():

    #Separations should be given in Mpc/h
    pimax = 40. #Mpc/h
    min_sep = 0.1
    max_sep = 40. #Mpc/h
    K = 12
    bin_size = np.log(max_sep / min_sep) / float(K)

    wp = True
    print 'Running'
    run_bins(min_sep, max_sep, bin_size, K, pimax, wp)


def run_bins(min_sep, max_sep, bin_size, K, pimax, wp):
    samplenums = [7, 8, 9, 10, 11, 12]
    #samplenums = [7]

    rpbins = np.logspace(np.log(min_sep), np.log(max_sep), K+1, base=np.e)
    rpbins_avg = 0.5 * (rpbins[1:] + rpbins[:-1])
    rps = [rpbins_avg]*len(samplenums)
    wprps = []
    labels = []
    cols = []
    for samplenum in samplenums:
        if samplenum in labels_mr:
            labels.append(labels_mr[samplenum])
        else:
            labels.append(samplenum)
        xi = run_sample(samplenum, min_sep, max_sep, bin_size, rpbins, pimax, wp)
        wprps.append(xi) #*2?
        cols.append(colors[samplenum])

    plotter.plot_wprp(rps, wprps, labels, colors=cols)

def run_sample(samplenum, min_sep, max_sep, bin_size, rpbins, pimax, wp):

    fn = '../data/lss.dr72bright{}_czcut.dat'.format(samplenum)
    data1 = pd.read_csv(fn)
    rand1 = get_random(samplenum)

    print 'Sample {}, {}'.format(samplenum, labels_mr[samplenum])

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    frac = 0.2
    data1 = data1.sample(frac=frac)
    rand1 = rand1.sample(frac=frac)

    print 'ndata=', len(data1.index)
    print 'nrand=', len(rand1.index)

    data2 = data1
    rand2 = rand1

    start = time.time()
    #xi, dd, dr, rd, rr = run.run_treecorr(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    #xi, dd, dr, rd, rr = run.run_treecorr_orig(data1, rand1, data2, rand2, min_sep, max_sep, bin_size, pimax, wp)
    est_ls, wprp = run.run_corrfunc(data1, rand1, data2, rand2, rpbins, pimax)

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