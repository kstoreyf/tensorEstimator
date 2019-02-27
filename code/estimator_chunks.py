import numpy as np
import multiprocessing as mp
import time

import estimator


c_kms = 3.e5
czmin = 0.02*c_kms
cz_lims = {7: [30900, 73500], 8: [19900, 47650], 9: [12600, 31900],
           10: [8050, 19250], 11: [5200, 12500], 12: [3200, 7850],
           21:[czmin, 73500], 20:[czmin, 59600], 19:[czmin, 47650], 18:[czmin, 39700], 17:[czmin, 31900], 16:[czmin, 25450], 15:[czmin, 19250], 14:[czmin, 15750], 13:[czmin, 12500]}

def est(ddpg, drpg, rdpg, rrpg, pimax, rmax, cosmo, basisfunc, K, wp, nproc, *args):

    nd1 = len(ddpg.cat1)
    nr1 = len(rrpg.cat1)
    nd2 = len(ddpg.cat2)
    nr2 = len(rrpg.cat2)

    if nproc==0:
        locs = None
        print 'Calculating dd vector'
        dds = project_pairs(ddpg, locs, pimax, rmax, cosmo, basisfunc, K, wp, False, None, None, *args)
        print dds, sum(np.array(dds).flatten())
        #print sdfdf
        print 'Calculating dr vector'
        drs = project_pairs(drpg, locs, pimax, rmax, cosmo, basisfunc, K, wp, False, None, None, *args)
        print drs, sum(np.array(drs).flatten())

        print 'Calculating rd vector'
        # should check if need to compute both aka if d1==d2 (drs and rds should be same bc using symmetric defs of dcm)
        rds = project_pairs(rdpg, locs, pimax, rmax, cosmo, basisfunc, K, wp, False, None, None, *args)

        print 'Calculating rr vector'
        rrs, qqs = project_pairs(rrpg, locs, pimax, rmax, cosmo, basisfunc, K, wp, True, None, None, *args)
        print rrs, sum(np.array(rrs).flatten())

    else:
        pair_arrs = [[ddpg, False],
                     [drpg, False],
                     [rdpg, False],
                     [rrpg, True]]

        print mp.cpu_count(), nproc

        bnum = len(basisfunc)
        count_arrs = [[] for _ in range(len(basisfunc))]
        for pg, tensor in pair_arrs:
            count_arr = mp.Array('d', np.zeros(K*bnum))

            if tensor:
                count_tensor_arr = mp.Array('d', np.zeros(K*K*bnum))
            else:
                count_tensor_arr = None

            start = time.time()
            subsize = int(np.ceil(float(len(pg.cat1))/float(nproc)))

            #loc_arr = [(nn*subsize, min((nn+1)*subsize, len(pg.cat1))) for nn in range(nproc)]
            loc_arr = split(range(len(pg.cat1)), nproc)
            print loc_arr

            pargs = [[pg, locs, pimax, rmax, cosmo, basisfunc, K, wp, tensor]+[count_arr, count_tensor_arr]+list(args) for locs in loc_arr]
            processes = [mp.Process(target=project_pairs, args=parg) for parg in pargs]
            print 'Starting'
            for p in processes:
                p.start()
            for p in processes:
                p.join()

            for bb in range(bnum):
                print count_arr[bb*K:(bb+1)*K]
                count_arrs[bb].append(np.array(count_arr[bb*K:(bb+1)*K]))
                if tensor:
                    count_arrs[bb].append(np.array(count_tensor_arr[bb*K*K:(bb+1)*K*K]).reshape(K,K))

            end = time.time()
            print 'TIME:', end-start

    a = []
    for bb in range(len(basisfunc)):
        if nproc>0:
            dds, drs, rds, rrs, qqs = zip(*count_arrs)
        print dds[bb]
        print drs[bb]
        print rrs[bb]
        print qqs[bb]
        dd = dds[bb] * 1./(nd1*nd2)
        dr = drs[bb] * 1./(nd1*nr2)
        rd = rds[bb] * 1./(nd2*nr1)
        rr = rrs[bb] * 1./(nr1*nr2)
        qq = qqs[bb] * 1./(nr1*nr2)
        a.append(calc_amplitudes(dd, dr, rd, rr, qq))
    print 'Computed amplitudes'

    return a


def split(a, n):
    k, m = divmod(len(a), n)
    return [(i * k + min(i, m), (i + 1) * k + min(i + 1, m)) for i in xrange(n)]
    #return list((a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n)))


def calc_amplitudes(dd, dr, rd, rr, qq):
    print 'Calculating amplitudes'
    qqinv = np.matrix(qq).I
    numerator = dd - dr - rd + rr
    a = np.matmul(qqinv, numerator)
    return a.A1


def project_pairs(pg, locs, pimax, rmax, cosmo, basisfunc, K, wp, tensor, count_arr, count_tensor_arr, *args):

    bnum = len(basisfunc)
    counts = [np.zeros(K) for _ in range(bnum)]
    counts_tensor = [np.zeros((K, K)) for _ in range(bnum)]
    outarr = np.empty((K, K))

    dcm1 = pg.cat1['dcm_mpc'].values
    dcm2 = pg.cat2['dcm_mpc'].values

    h = cosmo.h

    if not locs:
        locs = (0, len(pg.cat1))

    startloc, endloc = locs
    loc = startloc

    #TODO: why do i need this??
    pimax *= h

    while loc<endloc:
        pairs = pg.get_neighbors(loc)
        for i,j,r in pairs:

            if wp:
                # pi = h * abs(dcm1[i] - dcm2[j])
                pi = get_pi(dcm1, dcm2, int(i), int(j), h)
                #TODO: implement binning in pi
                #if pi>pimax:
                #    continue
                #r = rp_unit2mpch(dcm1_transverse, dcm2_transverse, i, j, r, h)
                # dcm_transverse_avg = 0.5 * (dcm1_transverse[i] + dcm2_transverse[j])
                # r *= dcm_transverse_avg * h
            else:
                pi = None

            # Now this is rp if wp or just r if 3d
            if r<=rmax:
                for bb in range(bnum):
                    #print (i, j, r)
                    u_ij = basisfunc[bb](pg.cat1, pg.cat2, i, j, r, pi, *args)

                    counts[bb] += u_ij

                    if tensor:
                        if u_ij.any() != 0:
                            ein = np.einsum('i,j', u_ij, u_ij, out=outarr)
                            counts_tensor[bb] += ein

        #pairs = pg.get_neighbors(loc)
        loc += 1

    if tensor:
        if count_arr:
            for bb in range(bnum):
                count_arr[bb*K:(bb+1)*K] += counts[bb]
                ct = counts_tensor[bb].reshape(K*K)
                count_tensor_arr[bb*K*K:(bb+1)*K*K] += ct
        else:
            return counts, counts_tensor
    else:
        if count_arr:
            for bb in range(bnum):
                count_arr[bb*K:(bb+1)*K] += counts[bb]
        else:
            return counts


def get_pi(dcm1, dcm2, i, j, h):
    pi = h * abs(dcm1[i] - dcm2[j])
    return pi

# Goes from unit rp to physical rp
def rp_unit2mpch(dcm1_t, dcm2_t, i, j, r, h):
    dcm_transverse_avg = 0.5 * (dcm1_t[i] + dcm2_t[j])
    r *= dcm_transverse_avg * h
    return r

def get_obj(cat, idx):
    return cat.iloc[idx]

# def calc_rp(cat1, cat2, i, j, cosmo):
#
#     unitdist = np.sqrt((cat1['xproj'][i] - cat2['xproj'][j])**2
#                      + (cat1['yproj'][i] - cat2['yproj'][j])**2
#                      + (cat1['zproj'][i] - cat2['zproj'][j])**2)
#
#     #making dcm_transverse average of two but idk what is correct
#     dcm_transverse = 0.5*(cat1['dcm_transverse_mpc'][i]+cat2['dcm_transverse_mpc'][j])
#     rp = unitdist * dcm_transverse * cosmo.h
#     return rp

def top(x, peak, width):
    if np.isscalar(x):
        if (x>peak-width/2.)&(x<peak+width/2.):
            return 1
        return 0
    t = np.zeros(len(x))
    t[(x>peak-width/2.)&(x<peak+width/2.)] = 1
    return t

def piece(x, peak, width):
    p = np.array(1 - (1. / width) * abs(x - peak))
    p[(x<peak-width)^(x>peak+width)] = 0
    return p

def gauss(x, peak, width):
    sigma = width/(2.*np.sqrt(2*np.log(2.)))
    return np.array(np.exp(-(x - peak) ** 2 / (2. * sigma ** 2)))

def tophat(cat1, cat2, i, j, rp, logbins_avg, logwidth):
    logrp = np.log10(rp)
    # u = np.zeros(len(logbins_avg))
    # for i in range(len(u)):
    #     peak = logbins_avg[i]
    #     u[i] = top(logrp, peak, logwidth)
    u = np.array([top(logrp, peak, logwidth) for peak in logbins_avg])
    return u

def tophat_robust(cat1, cat2, i, j, rp, logbins, logwidth):

    rp = np.log10(rp)
    ins_rp = -1
    for nn in range(len(logbins)-1):
        if logbins[nn] <= rp and rp < logbins[nn+1]:
            ins_rp = nn
            break
    u = np.zeros(len(logbins)-1)
    if ins_rp>=0 and ins_rp<len(logbins)-1:
        u[ins_rp] = 1
    return u

def tophat_fast(cat1, cat2, i, j, rp, logbins, logwidth):
    #print logbins
    #print rp
    ins = int((np.log10(rp) - min(logbins)) / logwidth)
    #print ins
    #print
    u = np.zeros(len(logbins)-1)
    if ins>=0 and ins<len(logbins)-1:
        u[ins] = 1
    return u

def tophat_robust(cat1, cat2, i, j, rp, logbins, logwidth):
    
    rp = np.log10(rp)
    ins_rp = -1
    for nn in range(len(logbins)-1):
        if logbins[nn] <= rp and rp < logbins[nn+1]:
            ins_rp = nn
            break 
    u = np.zeros(len(logbins)-1)
    if ins_rp>=0 and ins_rp<len(logbins)-1:
        u[ins_rp] = 1
    return u


def tophat_xis(cat1, cat2, i, j, rp, pi, bins, width):
    if pi:
        s = np.sqrt(rp**2 + pi**2)
    else:
        s = rp
    #if s <= max(bins):
    #    print rp, pi, s
    ins_s = -1
    for nn in range(len(bins) - 1):
        if bins[nn] <= s and s < bins[nn + 1]:
            ins_s = nn
            break
    u = np.zeros(len(bins) - 1)
    if ins_s >= 0 and ins_s < len(bins) - 1:
        u[ins_s] = 1
    return u

def tophat_xis_pow(cat1, cat2, i, j, rp, pi, bins, width):
    if pi:
        s = np.sqrt(rp**2 + pi**2)
    else:
        s = rp

    u = []
    s0 = 50
    gamma = 2
    powlaw = (s/s0)**-gamma
    u.append(powlaw)

    amp = 0.1
    sigma = 5
    mean = 107
    bump = amp*np.exp(-(s-mean)**2/(2*sigma**2))
    u.append(bump)

    return np.array(u)

#def quad(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None)

def top_z(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None):
    logrp = np.log10(rp)
    u = np.array([top(logrp, peak, logwidth) for peak in logbins_avg])
    if not val:
        val = 0.5*(cat1['z'][i] + cat2['z'][j])
    u = np.concatenate((u, u*val, u*val**2))
    return u

def top_Mrz(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None):
    logrp = np.log10(rp)
    u = np.array([top(logrp, peak, logwidth) for peak in logbins_avg])
    Mrsun = 4.46
    if val:
        val = val
        #val = 10**(0.4*(Mrsun-val))
    #   weight = 1.0
    else:
        m1 = cat1['M_rz'].values[i]
        m2 = cat2['M_rz'].values[j]
        l1 = 10**(0.4*(Mrsun-m1))
        l2 = 10**(0.4*(Mrsun-m2))
        #val = 0.5*(cat1['M_rz'][i] + cat2['M_rz'][j])
        #diff = abs(cat1['M_rz'][i]-cat2['M_rz'][j])
        #weight = (abs(val) - diff)/abs(val)
        lval = np.mean([l1,l2])
        mval = Mrsun - 2.5*np.log10(lval)
        val = mval
        #if diff > 1:
        #    weight = 0.0
        #else:
        #    weight = 1.0
    #u = np.concatenate((u, u*abs(val), u*abs(val)**2))
    u = np.concatenate((u, u*val))
    #u = u*abs(val)
    return u

# be sure inputting bin edges not centers/avgs
def grid_Mrz(cat1, cat2, i, j, rp, logbins, logwidth, val=None):

    Mrzbins = np.linspace(-23, -18, 6)
    nrpbins = len(logbins)-1
    nMrzbins = len(Mrzbins)-1
    # Calc Mrz val
    Mrsun = 4.46
    if not val:
        m1 = cat1['M_rz'].values[i]
        m2 = cat2['M_rz'].values[j]
        l1 = 10**(0.4*(Mrsun-m1))
        l2 = 10**(0.4*(Mrsun-m2))
        lval = np.sum([l1,l2])
        mval = Mrsun - 2.5*np.log10(lval)
        val = mval

    nrpbins = len(logbins)-1
    nMrzbins = len(Mrzbins)-1
    # Mrz bin
    ins_Mrz = -1
    for mm in range(nMrzbins):
      if Mrzbins[mm] <= val and val < Mrzbins[mm+1]:
        ins_Mrz = mm
        break  

    # rp bin
    #ins_rp = int((np.log10(rp) - min(logbins)) / logwidth)
    rp = np.log10(rp)
    ins_rp = -1
    for nn in range(nrpbins):
        if logbins[nn] <= rp and rp < logbins[nn+1]:
            ins_rp = nn
            break 

    Kval = nrpbins*nMrzbins
    u = np.zeros(Kval)

    if ins_rp>=0 and ins_rp<nrpbins and ins_Mrz>=0 and ins_Mrz<nMrzbins:
        u[ins_Mrz*nrpbins + ins_rp] = 1

    return u

# be sure inputting bin edges not centers/avgs
def match_bins(cat1, cat2, i, j, rp, logbins, logwidth, val=None):

    Mrzbins = np.linspace(-23, -18, 6)
    #Mrzbins = np.array([-19,-18]) #11
    ##Mrzbins = np.array([-20,-19]) #10
    #Mrzbins = np.array([-21,-20]) #9
    #Mrzbins = np.array([-22,-21]) #8
    Mrzbins = np.array([-23,-22]) #7
    #Mrzbins = np.array([-19,-18])
    #Mrzbins = np.array([-22,-21])
    #Mrzbins = np.array([-23, -22])

    nrpbins = len(logbins)-1
    nMrzbins = len(Mrzbins)-1
    #samplenum = 11
    #print Mrzbins
    if val:
        ins_Mrz = -1
        for mm in range(nMrzbins):
          if Mrzbins[mm] <= val and val < Mrzbins[mm+1]:
            ins_Mrz = mm
            break
    else:
        ins_Mrz = -1
        m1 = cat1['M_rz'].values[i]
        m2 = cat2['M_rz'].values[j]
        for mm in range(nMrzbins):
            if Mrzbins[mm] <= m1 and m1 < Mrzbins[mm + 1] \
                    and Mrzbins[mm] <= m2 and m2 < Mrzbins[mm + 1]:
                #if cz_lims[samplenum][0]...
                ins_Mrz = mm
                break

    # rp bin
    #ins_rp = int((np.log10(rp) - min(logbins)) / logwidth)
    rp = np.log10(rp)
    ins_rp = -1
    for nn in range(nrpbins):
        if logbins[nn] <= rp and rp < logbins[nn+1]:
            ins_rp = nn
            break

    Kval = nrpbins*nMrzbins
    u = np.zeros(Kval)

    if ins_rp>=0 and ins_rp<nrpbins and ins_Mrz>=0 and ins_Mrz<nMrzbins:
        u[ins_Mrz*nrpbins + ins_rp] = 1

    return u


def top_rand(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None):
    logrp = np.log10(rp)
    u = np.array([top(logrp, peak, logwidth) for peak in logbins_avg])
    if not val:
        nr1 = np.random.rand() * 17 - 30
        nr2 = np.random.rand() * 17 - 30
        val = 0.5 * (nr1 + nr2)
        #val = 0.5*(np.random.rand() + np.random.rand())
    #u = u*abs(val)
    u = np.concatenate((u, u*abs(val), u*abs(val)**2))
    return u

def gauss_Mr(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None):
    logrp = np.log10(rp)
    u = np.array([gauss(logrp, peak, logwidth) for peak in logbins_avg])
    if not val:
        val = 0.5*(cat1['M_r'][i] + cat2['M_r'][j])
    u = np.concatenate((u, u*val, u*val**2))
    return u

def gauss_z(cat1, cat2, i, j, rp, logbins_avg, logwidth, val=None):
    logrp = np.log10(rp)
    u = np.array([gauss(logrp, peak, logwidth) for peak in logbins_avg])
    if not val:
        val = 0.5*(cat1['z'][i] + cat2['z'][j])
    u = np.concatenate((u, u*val, u*val**2))
    return u

def piecewise(cat1, cat2, i, j, rp, logbins_avg, logwidth):
    logrp = np.log10(rp)
    u = np.array([piece(logrp, peak, logwidth) for peak in logbins_avg])
    return u

def gaussian(cat1, cat2, i, j, rp, logbins_avg, logwidth):
    logrp = np.log10(rp)
    u = np.array([gauss(logrp, peak, logwidth) for peak in logbins_avg])
    return u

def trig(cat1, cat2, i, j, rp, logbins_avg, logwidth):
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
