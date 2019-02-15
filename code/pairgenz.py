import numpy as np
from scipy.spatial import KDTree
from scipy.spatial import cKDTree



class PairGen():


    def __init__(self, cat1, cat2, rmax, cosmo, wp, pimax, chunksize=20):
        assert len(cat2) >= len(cat1), 'The second catalog should be the larger one'

        cat1['idx'] = range(len(cat1))
        cat2['idx'] = range(len(cat2))
        self.cat1 = cat1
        self.cat2 = cat2
        self.rmax = rmax
        self.cosmo = cosmo
        self.wp = wp
        self.pimax = pimax
        self.chunksize = chunksize
        self.construct_trees_zshells()
        self.num_pairs = 0



    def get_pairchunk(self, cat1start):

        #print self.cat1loc
        #print self.chunksize

        # if self.cat1loc >= min(self.cat1stop, len(self.cat1)):
        #     return []

        cat1loc = cat1start
        pairchunk = []
        while len(pairchunk) < self.chunksize:
            pairs = self.get_neighbors(cat1loc)
            # Eliminate self-pairs with zero separation
            pairs = [p for p in pairs if p[2] > 0]
            #print len(pairs)
            pairchunk += pairs
            cat1loc += 1
            #print len(pairchunk)
            #print
        return pairchunk, cat1loc


    def get_neighbors(self, cat1loc):

        i = cat1loc
        zshell = self.cat1['zshell'].values[i]
        tree = self.trees[zshell]
        treecat = self.treecats[zshell]

        # this might be possible if there's a big chunk of non-overlap bw catalogs
        # (cross-correlations??) #check
        if tree is None:
            print 'No tree!'
            return []

        ipoint = np.array([self.cat1['xproj'].values[i], self.cat1['yproj'].values[i], self.cat1['zproj'].values[i]])
        if not self.wp:
            ipoint *= self.cat1['dcm_mpc'].values[i]  # turn projected into real space
            rmax_tree = self.rmax
        if self.wp:
            # here the given rmax rmax is really rpmax
            rmax_tree = self.rmax / (self.cat1['dcm_transverse_mpc'].values[i] * self.cosmo.h)  # turn bin into unit dist

        ntree = len(tree.data)
        dists, locs = tree.query(ipoint, k=ntree, distance_upper_bound=rmax_tree)

        # Query returns infinities when <k neighbors are found, cut these
        if locs[-1] == ntree:
            imax = next(index for index, value in enumerate(locs) if value == ntree)
            locs = locs[:imax]
            dists = dists[:imax]

        if self.wp:
            dists = np.array([dists[k] * 0.5 * (self.cat1['dcm_transverse_mpc'].values[i]
                                                + treecat['dcm_transverse_mpc'].values[locs[k]]) for k in range(len(locs))])

        # if wp: dists are transverse distances in mpc/h
        # if not wp: dists are real space 3d dists in mpc/h
        # TODO: don't know why don't need this??
        #dists *= self.cosmo.h

        js = treecat['idx'].values
        pairs = [(i, js[locs[k]], dists[k]) for k in range(len(locs))]
        #pairs = [(i, locs[k], dists[k]) for k in range(len(locs))]
        #Eliminate self-pairs (though they shouldn't pass through later if rmin > 0
        pairs = [p for p in pairs if p[2] > 0]
        #self.cat1loc += 1
        #print i, len(pairs)
        self.num_pairs += len(pairs)
        return pairs



    def construct_trees_zshells(self):

        zedges = []

        dcmh1 = self.cat1['dcm_mpc'] * self.cosmo.h
        dcmh2 = self.cat2['dcm_mpc'] * self.cosmo.h
        dcmh_min = min(min(dcmh1), min(dcmh2))
        dcmh_max = max(max(dcmh1), max(dcmh2))

        zmin = np.floor(dcmh_min)
        dshell = 4. * float(self.pimax)
        while zmin < dcmh_max:
            zedges.append((zmin, zmin + dshell))
            #probs a nicer way to do this in loop, but without then
            #sometimes adds an extra unnecesarry shell
            if zmin + dshell > dcmh_max:
                    break
            zmin += dshell/2.

        self.zedges = zedges

        self.cat1['zshell'] = self.cat1['dcm_mpc'].apply(self.pick_zshell,
                                                    args=(self.cosmo.h, zedges, dshell))

        # this splits up cat1 by which shell they belong to - but don't think
        # it's needed unless you're parallelizing on the shell groups rather than
        # my current way
        #groups1 = self.cat1.groupby(self.cat1['zshell'])
        #dfs1 = []
        #for name, group in groups1:
        #    dfs1.append(group.reset_index())

        print 'Constructing KDTrees...'
        self.trees = []
        self.treecats = []
        for i in range(len(zedges)):
            zedge = zedges[i]
            shell = self.cat2.loc[(zedge[0] <= dcmh2) & (dcmh2 < zedge[1])]
            points = np.array([shell['xproj'], shell['yproj'], shell['zproj']])
            if not self.wp:
                points *= shell['dcm_mpc']

            if len(shell) > 1:
                tree = cKDTree(list(points.T))
                self.trees.append(tree)
                self.treecats.append(shell)
            else:
                self.trees.append(None)
                self.treecats.append(None)
        print 'N_trees:', len(self.trees)
        print 'pimax:', self.pimax, 'dshell:', dshell
        print 'dcms - min:', dcmh_min, 'max:', dcmh_max, 'diff:', dcmh_max - dcmh_min
        print 'zedges (dcm)', zedges
        print 'Datapoints in trees:', [len(t.data) for t in self.trees]


    # def a smarter way to do this but
    def pick_zshell(self, dcm, h, zedges, dshell):
        dcmh = dcm * h
        for i in range(len(zedges)):
            zshell = zedges[i]
            if zshell[0] <= dcmh < zshell[1] and abs(dcmh - zshell[0]) > dshell / 4. \
                    and abs(dcmh - zshell[1]) > dshell / 4.:
                #print dcmh, zshell[0], zshell[1]
                return i
        # edge cases
        if zedges[0][0] <= dcmh < zedges[0][1]:
            #print dcmh, zedges[0][0], zedges[0][1]
            return 0
        last = len(zedges) - 1
        if zedges[last][0] <= dcmh < zedges[last][1]:
            #print dcmh, zedges[last][0], zedges[last][1]
            return last
        print 'Error: not sorted into redshift shell'
        return None