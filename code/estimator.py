import numpy as np




def est(d1d2pairs, d1r2pairs, d2r1pairs, r1r2pairs, data1, rand1, data2, rand2,
        pimax, rpmax, cosmo, basisfunc, K, wp, *args):

    nd1 = len(data1)
    nr1 = len(rand1)
    nd2 = len(data2)
    nr2 = len(rand2)

    print 'Calculating dd vector'
    dd = 1./(nd1*nd2) * project_pairs(d1d2pairs, data1, data2, pimax, rpmax,
                                      cosmo, basisfunc, K, wp, False, *args)
    # is d1 r1 right?
    dr = 1./(nd1*nr2) * project_pairs(d1r2pairs, data1, rand2, pimax, rpmax,
                                      cosmo, basisfunc, K, wp, False, *args)
    # could check if need to compute both aka if d1==d2 (tho will be slightly diff even w same catalog)
    rd = 1./(nd2*nr1) * project_pairs(d2r1pairs, data2, rand1, pimax, rpmax,
                                      cosmo, basisfunc, K, wp, False, *args)
    print 'Calculated dr, rd vector'

    rr, qq = project_pairs(r1r2pairs, rand1, rand2, pimax, rpmax,
                                      cosmo, basisfunc, K, wp, True, *args)
    rr = 1./(nr1*nr2) * rr
    qq = 1./(nr1*nr2) * qq
    print 'Calculated rr vector'

    a = calc_amplitudes(dd, dr, rd, rr, qq)

    return a


def calc_amplitudes(dd, dr, rd, rr, qq):
    print qq
    qqinv = np.matrix(qq).I
    numerator = dd - dr - rd + rr
    a = np.matmul(qqinv, numerator)
    return a


def project_pairs(pairs, cat1, cat2, pimax, rpmax, cosmo, basisfunc, K, wp, tensor, *args):

    counts = np.zeros(K)
    counts_tensor = np.zeros((K, K))
    outarr = np.empty((K, K))

    dcm1 = cat1['dcm_mpc'].values
    dcm2 = cat2['dcm_mpc'].values
    dcm1_transverse = cat1['dcm_transverse_mpc'].values
    dcm2_transverse = cat2['dcm_transverse_mpc'].values
    h = cosmo.h

    for i,j,r in pairs:

        if wp:
            pi = cosmo.h * abs(dcm1[i] - dcm2[j])
            if wp and pi>pimax:
                continue
            dcm_transverse_avg = 0.5 * (dcm1_transverse[i] + dcm2_transverse[j])
            r *= dcm_transverse_avg * h

        if r<=rpmax:
            u_ij = basisfunc(cat1.iloc[i], cat2.iloc[j], r, *args)
            counts += u_ij

            if tensor:
                if u_ij.any() != 0:
                    ein = np.einsum('i,j', u_ij, u_ij, out=outarr)
                    counts_tensor += ein

    if tensor:
        return counts, counts_tensor
    else:
        return counts


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

def calc_rp(obj1, obj2, dist, cosmo):

    #making dcm_transverse average of two but idk what is correct
    dcm_transverse = 0.5*(obj1['dcm_transverse_mpc']+obj2['dcm_transverse_mpc'])
    rp = dist * dcm_transverse * cosmo.h
    return rp

def tophat(obj1, obj2, rp, rpbins):
    u = np.zeros(len(rpbins))
    pos = np.digitize(rp, rpbins)
    u[pos - 1] = 1
    return u


def gaussian(obj1, obj2, rp, bins):
    mu = rp
    pos = np.digitize(rp, bins)
    rploc = pos - 1
    sigma = bins[rploc] - bins[rploc - 1]
    return np.array(np.exp(-(bins - mu) ** 2 / (2. * sigma ** 2)))

# def piecewise(obj1, obj2, rp, rpbins, rpbins_fine):
#     pos = np.digitize(rp, rpbins)
#     rploc = pos-1
#     u = 1-abs(rpbins_fine-rploc)
#     u = [0 if uu<0 else uu for uu in u]
#     return u

def piecewise(obj1, obj2, rp, rpbins_fine):
    pos = np.digitize(rp, rpbins_fine)
    rploc = pos - 1
    width = rpbins_fine[rploc] - rpbins_fine[rploc - 1]
    u = 1 - (1. / (5 * width)) * abs(rpbins_fine - rp)
    u = np.array([0 if uu < 0 else uu for uu in u])
    return u

# def tophat(obj1, obj2, rp, rpbins):
#     print rpbins
#     u = np.zeros(len(rpbins))
#     pos = np.digitize(rp, rpbins)
#     u[pos-1] = 1
#     return u