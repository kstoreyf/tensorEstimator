import numpy as np
from matplotlib import pyplot as plt



def main():
    check()


def check():
    dd_res_corrfunc, dr_res_corrfunc, rr_res_corrfunc = np.load(
        '../results/counts_bin20_frac0.1_randcz.npy')

    min_sep = 0.13
    max_sep = 40. #Mpc/h
    bin_size = 0.2
    pimax = 40

    pibinwidth = 1
    pibins = np.arange(0, pimax + pibinwidth, pibinwidth)

    K = (np.log10(max_sep) - np.log10(min_sep))/bin_size
    rpbins = np.logspace(np.log10(min_sep), np.log10(max_sep), K+1)

    dd_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
    dr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))
    rr_rp_pi_corrfunc = np.zeros((len(pibins) - 1, len(rpbins) - 1))

    rp_avg = np.zeros(len(rpbins) - 1)
    for m in range(len(pibins) - 1):
        for n in range(len(rpbins) - 1):
            idx = (len(pibins) - 1) * n + m
            # = count * avg weight
            dd_rp_pi_corrfunc[m][n] = dd_res_corrfunc[idx][4] * dd_res_corrfunc[idx][5]
            dr_rp_pi_corrfunc[m][n] = dr_res_corrfunc[idx][4] * dr_res_corrfunc[idx][5]
            rr_rp_pi_corrfunc[m][n] = rr_res_corrfunc[idx][4] * rr_res_corrfunc[idx][5]
            rp_avg[n] += dd_res_corrfunc[idx][2] + dr_res_corrfunc[idx][2] + \
                         rr_res_corrfunc[idx][2]
    for n in range(len(rpbins) - 1):
        rp_avg[n] /= 3 * (len(pibins) - 1)

    dd = np.sum(dd_rp_pi_corrfunc, axis=0)
    dr = np.sum(dr_rp_pi_corrfunc, axis=0)
    rr = np.sum(rr_rp_pi_corrfunc, axis=0)

    plot(dd, dr, rr, rp_avg)


def plot(dd, dr, rr, rp):
    fig = plt.figure()
    ax0 = fig.gca()
    rp = [0.17, 0.27, 0.42, 0.67, 1.1, 1.7, 2.7, 4.2, 6.7, 10.6, 16.9, 26.8]

    ax0.loglog(rp, dd, label='DD', marker='o', color='blue')
    ax0.loglog(rp, dr, label='DR', marker='o', color='green')
    ax0.loglog(rp, rr, label='RR', marker='o', color='black')
    plt.legend(loc='best')
    plt.show()

if __name__=="__main__":
    main()