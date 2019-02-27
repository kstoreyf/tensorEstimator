import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from astropy.cosmology import LambdaCDM
from time import time




def main():
    check_kdtrees()
    #check_kdtrees_periodic()


def average_time(executable, *args, **kwargs):
    """Compute the average time over N runs"""
    N = 1
    t = 0
    for i in range(N):
        t0 = time()
        res = executable(*args, **kwargs)
        t1 = time()
        t += (t1 - t0)
    return res, t * 1. / N




def check_kdtrees_periodic():
    from scipy.spatial import KDTree
    from scipy.spatial import cKDTree

    ### Load dataframe ###
    fn = '../../emulator/CMASS_BIAS/COS_2000HOD/galaxy_mock/cosmo_0_HOD_0_test_0.mock'
    rand = pd.read_table(fn, index_col=False, delim_whitespace=True, names=['x','y','z','vx','vy','vz','Mh'])
    print len(rand.index)

    print min(rand['x']), min(rand['y']), min(rand['z'])
    print max(rand['x']), max(rand['y']), max(rand['z'])

    print 'Downsampling'
    #frac = 1
    #rand = rand.sample(frac=frac)
    #print len(rand.index)

    rmax = 40.
    #boxsize = 1050.000001
    boxsize = None
    print 'box:',boxsize
    ### Build trees ###
    points = np.array([rand['x'], rand['y'], rand['z']])

    print 'Building ckdtree'
    ckd = cKDTree(list(points.T), boxsize=boxsize)
    #ckd, t = average_time(cKDTree, list(points.T))
    #print 'time: {:.3f}s'.format(t)

    ### Query ###
    n = 1
    for i in range(n):

        ipoint = points[:, i]

        print 'Querying ckdtree'
        print ipoint
        ntree = len(ckd.data)
        # dists, locs = ckd.query(ipoint, k=ntree, distance_upper_bound=rmax_tree)
        kw = {'k': ntree, 'distance_upper_bound': rmax}
        res, t = average_time(ckd.query, ipoint, **kw)
        print 'time: {:.8f}ms'.format(t * 1.e3)
        dists, locs = res
        if locs[-1] == ntree:
            imax = next(index for index, value in enumerate(locs) if value == ntree)
            locs = locs[:imax]
            dists = dists[:imax]
            print len(dists)
        #print dists
        #print locs
        print imax


def check_kdtrees():
    from scipy.spatial import KDTree
    from scipy.spatial import cKDTree

    ### Load dataframe ###
    print 'Loading data'
    sample = 'Full'
    randfn = '../data/random-DR7-{}.ascii'.format(sample)
    rand = pd.read_table(randfn, index_col=False, delim_whitespace=True, names=['ra', 'dec', 'z',
                                                                                'sector_completeness', 'n(z)*1e4',
                                                                                'radial_weight', 'ilss', 'sector'],
                         dtype={'z': np.float64},
                         skiprows=1)
    print len(rand.index)

    print 'Downsampling'
    frac = 1
    rand = rand.sample(frac=frac)
    print len(rand.index)

    print 'Adding info'
    add_info(rand)
    print min(rand['dcm_mpc']), max(['dcm_mpc'])
    pimax = 60.
    print (max(['dcm_mpc'])-min(rand['dcm_mpc']))/(pimax*2)

        ### Set up ###
    print 'Setting up'
    rmax = 40.
    cosmo = LambdaCDM(H0=70, Om0=0.25, Ode0=0.75)

    ### Build trees ###
    points = np.array([rand['xproj'], rand['yproj'], rand['zproj']])
    print 'Building kdtree'
    kd = KDTree(list(points.T))
    #kd, t = average_time(KDTree, list(points.T))
    #print 'time: {:.3f}s'.format(t)

    print 'Building ckdtree'
    ckd = cKDTree(list(points.T))
    #ckd, t = average_time(cKDTree, list(points.T))
    #print 'time: {:.3f}s'.format(t)

    ### Query ###
    n = 5
    for i in range(n):

        rmax_tree = rmax / (rand['dcm_transverse_mpc'].values[i] * cosmo.h)  # turn bin into unit dist
        ipoint = points[:, i]

        print 'Querying kdtree'
        ntree = len(kd.data)
        # dists, locs = kd.query(ipoint, k=ntree, distance_upper_bound=rmax_tree)
        kw = {'k': ntree, 'distance_upper_bound': rmax_tree}
        res, t = average_time(kd.query, ipoint, **kw)
        print 'time: {:.8f}ms'.format(t * 1.e3)

        print 'Querying ckdtree'
        ntree = len(ckd.data)
        # dists, locs = ckd.query(ipoint, k=ntree, distance_upper_bound=rmax_tree)
        kw = {'k': ntree, 'distance_upper_bound': rmax_tree}
        res, t = average_time(ckd.query, ipoint, **kw)
        print 'time: {:.8f}ms'.format(t * 1.e3)


def add_info(df):
    df['xproj'], df['yproj'], df['zproj'] = zip(*df.apply(ra_dec_to_unitxyz, axis=1))

    df['dcm_mpc'] = df['z'].apply(get_comoving_dist)
    df['dcm_transverse_mpc'] = df['dcm_mpc']


def ra_dec_to_unitxyz(row):
    ra = row['ra'] * np.pi / 180.
    dec = row['dec'] * np.pi / 180.
    x = np.cos(ra) * np.cos(dec)
    y = np.sin(ra) * np.cos(dec)
    z = np.sin(dec)
    return x, y, z


cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

def get_comoving_dist(z):
    comov = cosmo.comoving_distance(z)
    return comov.value*cosmo.h

if __name__=='__main__':
    main()