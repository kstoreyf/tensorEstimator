import pandas as pd
import numpy as np
import time
from astropy.cosmology import LambdaCDM


# globals
c_kms = 3.e5  # c in km/s
czmin = 0.02*c_kms
cz_lims = {7: [30900, 73500], 8: [19900, 47650], 9: [12600, 31900],
           10: [8050, 19250], 11: [5200, 12500], 12: [3200, 7850],
           21:[czmin, 73500], 20:[czmin, 59600], 19:[czmin, 47650], 18:[czmin, 39700],
           17:[czmin, 31900], 16:[czmin, 25450], 15:[czmin, 19250], 14:[czmin, 15750],
           13:[czmin, 12500]}


def main():
    samplenums = [7,8,9,10,11]
    tag = '_square5k'
    #combine_bins(samplenums, tag)
    write_samples()


def write_samples():
    #samplenums = [7,8,9,10,11]
    samplenums = [0]
    #tags = ['_czcut']
    squaretag = '_square5k'
    for samplenum in samplenums:
        merge_photo(samplenum)
        #write_sample(samplenum, tags)
        data, rand = make_testsquare(samplenum, squaretag)
        rand = write_random(data, squaretag, str(samplenum)+squaretag)

        print 'Data:'
        print len(data)
        print min(data['M_rz']), max(data['M_rz'])
        print min(data['cz']), max(data['cz'])

        print 'Rand:'
        print len(rand)
        print min(rand['M_rz']), max(rand['M_rz'])
        print min(rand['cz']), max(rand['cz'])





def czcut(samplenum, tags):


    fn = '../data/lss.dr72bright{}.dat'.format(samplenum)
    df = pd.read_csv(fn)

    if '_czcut' in tags:
        df = df[df['cz'] > cz_lims[samplenum][0]][df['cz'] < cz_lims[samplenum][1]]


def merge_photo(samplenum):

    print samplenum
    # Original data
    fn = '../data/lss.dr72bright{}.dat'.format(samplenum)
    df = pd.read_csv(fn, header=None, delim_whitespace=True, names=['indx',
                                                                    'sector', 'mregion', 'ra', 'dec', 'cz',
                                                                    'fgotten', 'selection_fn'])
    df['z'] = df['cz'] / c_kms

    fn_photo = '../data/photoinfo.dr72bright{}.dat'.format(samplenum)
    df_photo = pd.read_csv(fn_photo, header=None, delim_whitespace=True, names=['indx',
                                                                                'M_u', 'M_g', 'M_r', 'M_i',
                                                                                'M_z', 'mu_{50}', 'r50/r90'])

    df = pd.merge(df, df_photo, on='indx')  # This doesn't lose any data

    df['M_rz'] = calc_Mrz(df['M_r'], df['z'])

    print len(df.index)
    print min(df['cz']), max(df['cz'])
    print min(df['M_rz']), max(df['M_rz'])

    fn_save = '../data/lss.dr72bright{}_photo.dat'.format(samplenum)
    df.to_csv(fn_save)




def write_random(df_data, loadtag, savetag):
    fn = '../data/random-0.dr72bright{}.dat'.format(loadtag)
    # pretty sure fgotten is correct but couldn't find anywhere
    # df_rand = pd.read_csv(fn, header=None, delim_whitespace=True, names=['ra',
    #         'dec', 'sector', 'mregion', 'fgotten', 'min_mag?'])
    df_rand = pd.read_csv(fn)

    nrand = len(df_rand.index)

    df_rand['cz'] = df_data['cz'].sample(n=nrand, replace=True).values
    df_rand['z'] = df_rand['cz'] / c_kms

    df_rand['M_rz'] = df_data['M_rz'].sample(n=nrand, replace=True).values

    print len(df_rand.index)
    print min(df_rand['cz']), max(df_rand['cz'])
    print min(df_rand['M_rz']), max(df_rand['M_rz'])

    saveto = '../data/random-0.dr72bright{}.dat'.format(savetag)
    df_rand.to_csv(saveto)
    return df_rand


def make_testsquare(samplenum, squaretag):

    print samplenum
    fn_data = '../data/lss.dr72bright{}'.format(samplenum)
    fn_rand = '../data/random-0.dr72bright'

    if squaretag == '_square1k':
        ramin = 150
        ramax = 152.5
        decmin = 10
        decmax = 15
    elif squaretag == '_square5k':
        ramin = 150
        ramax = 155
        decmin = 10
        decmax = 25

    df_data = pd.read_csv(fn_data+'_photo.dat')
    df_rand = pd.read_csv(fn_rand+'.dat', header=None, delim_whitespace=True, names=['ra',
                'dec', 'sector', 'mregion', 'fgotten', 'min_mag?'])
    df_data = df_data[df_data['ra'] > ramin][df_data['ra'] < ramax][df_data['dec'] > decmin][df_data['dec'] < decmax]
    df_rand = df_rand[df_rand['ra'] > ramin][df_rand['ra'] < ramax][df_rand['dec'] > decmin][df_rand['dec'] < decmax]

    print len(df_data), len(df_rand)

    saveto_data = fn_data+squaretag+'.dat'
    saveto_rand = fn_rand+squaretag+'.dat'

    df_data.to_csv(saveto_data)
    df_rand.to_csv(saveto_rand)

    return df_data, df_rand

def combine_bins(samplenums, tag):
    datas = []
    rands = []
    for samplenum in samplenums:
        print samplenum

        fn = '../data/lss.dr72bright{}{}.dat'.format(samplenum, tag)
        fnrand = '../data/random-0.dr72bright{}{}.dat'.format(samplenum, tag)

        data = pd.read_csv(fn, index_col=0)
        rand = pd.read_csv(fnrand, index_col=0)
        print len(data), len(rand)
        datas.append(data)
        rands.append(rand)
    datadf = pd.concat(datas, ignore_index=True)
    randdf = pd.concat(rands, ignore_index=True)

    frac = 1
    datadf = datadf.sample(frac=frac).reset_index(drop=True)
    randdf = randdf.sample(frac=frac).reset_index(drop=True)

    #datadf['M_rz'] = calc_Mrz(datadf['M_r'], datadf['z'])

    print 'Nums'
    print len(datadf.index)
    print len(randdf.index)
    #datadf = run.add_info(datadf, zfile=None)
    #randdf = run.add_info(randdf, zfile=None)

    fn_save = '../data/lss.dr72bright{}{}.dat'.format('bins', tag)
    fnrand_save = '../data/random-0.dr72bright{}{}.dat'.format('bins', tag)

    datadf.to_csv(fn_save)
    randdf.to_csv(fnrand_save)



def calc_Mrz(Mr, z):
    q0 = 2.0
    q1 = -1.0
    qz0 = 0.1
    curr_zdep = q0 * (1.0 + q1 * (z - qz0))
    Mrz = Mr + curr_zdep * (z - qz0)
    return Mrz


if __name__=='__main__':
    main()