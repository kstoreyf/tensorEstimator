import numpy as np
from scipy.spatial import KDTree
from scipy.spatial import cKDTree
import random
import time


def main():


    tbuilds = []
    tqueries = []
    #treesizes = np.logspace(1e1, 1e7)
    treesizes = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]

    for treesize in treesizes:
        tbuild, tquery = timekdtree(treesize, periodic=True)
        tbuilds.append(tbuild)
        tqueries.append(tquery)

    fn = '../../clust/results/times/timekd_center_scipyperiodic_together.txt'
    np.savetxt(fn, np.array([treesizes, tbuilds, tqueries]).T, fmt='%d', delimiter=',')



def timekdtree(treesize, periodic=False):

    reps = 1000
    ndens = 0.00047
    L = (treesize/ndens)**(1./3.)
    r = L/10.

    if periodic:
        boxsize = L
    else:
        boxsize = None

    #corner
    # px = r/3.
    # py = r/3.
    # pz = r/3.
    # center
    px = L/2.
    py = L/2.
    pz = L/2.
    print treesize, L, r, px


    ### Build trees ###

    x = [random.random()*L for _ in range(int(treesize))]
    y = [random.random()*L for _ in range(int(treesize))]
    z = [random.random()*L for _ in range(int(treesize))]


    points = np.array([x,y,z])

    print 'Building ckdtree'
    start = time.time()
    ckd = cKDTree(list(points.T), boxsize=boxsize)
    end = time.time()
    tbuild = (end-start)*1000
    print 'Build time: {:.3f}ms'.format(tbuild)

    ### Query ###
    print 'Querying ckdtree'
    ntree = len(ckd.data)
    start = time.time()
    for i in range(reps):
        dists, locs = ckd.query([px,py,pz], k=ntree, distance_upper_bound=r)
    end = time.time()
    tquery = (end-start)*1000
    if locs[-1] == ntree:
        imax = next(index for index, value in enumerate(locs) if value == ntree)
    print 'found {} time: {:.8f}ms'.format(imax, tquery)

    return tbuild, tquery



if __name__=="__main__":
    main()