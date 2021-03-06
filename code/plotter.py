from matplotlib import pyplot as plt
import numpy as np

import run


color_dict = {'corrfunc':'grey', 'tophat':'blue', 'piecewise':'red',
          'tophat_orig':'cyan', 'gaussian':'green', 'trig':'magenta',
          'treecorr': 'orange', 'top_quad':'purple'}

color_list = ['magenta', 'red', 'black', 'limegreen', 'blue', 'cyan']
#color_list = ['magenta', 'red', 'limegreen', 'blue', 'cyan']

def main():
    #plot_dr72bins()
    #plot()
    #plot_bins()
    plot_bao()

def plot():
    #fn = "../results/dr72_brightLRG_frac0.1.npy"
    # fn = "../results/dr7_FullLRG_frac1_weights.npy"
    # ss, xis, labels = run.load_results(fn)
    # plot_xi_zspace(ss, xis, labels)

    sample = 'Dim-no'
    frac = 0.001
    #fn = "../results/wp_dr7_{}LRG_frac{}_weights.npy".format(sample, frac)
    fn = "../results/dr72_bin{}_frac{}_randcz.npy".format(20, frac)
    rps, wps, labels = run.load_results(fn)
    print rps
    print wps
    rps_paper = [[0.17, 0.27, 0.42, 0.67, 1.1, 1.7, 2.7, 4.2, 6.7, 10.6, 16.9, 26.8]]
    plot_wprp(rps_paper, wps, labels)


def plot_bao():
    #fn = "../results/bao/xis_dr7_FullLRG_frac1_rand0.5.npy"
    #fn = "../results/bao/xis_dr7_FullLRG_frac0.1_copy2.npy"
    fn = "../results/bao/xis_dr7_FullLRG_frac0.005_corrfunc.npy"

    ss, xis, aa, labels = np.load(fn)

    print ss
    print xis

    #ss = [ss[0], ss[1][0]]
    #xis = [xis[0], xis[1][0]]
    #print ss
    #print xis
    plot_xi_zspace(ss, xis, list(labels), xi_tocompare='corrfunc orig')
    #
    # sample = 'Dim-no'
    # frac = 0.1
    # #fn = "../results/wp_dr7_{}LRG_frac{}_weights.npy".format(sample, frac)
    # fn = "../results/dr72_bin{}_frac{}_randcz.npy".format(20, frac)
    # rps, wps, labels = run.load_results(fn)
    # print rps
    # print wps
    # rps_paper = [[0.17, 0.27, 0.42, 0.67, 1.1, 1.7, 2.7, 4.2, 6.7, 10.6, 16.9, 26.8]]
    # plot_wprp(rps_paper, wps, labels)





def plot_dr72bins():
    rps = []
    wprps = []
    labels = []
    #samplenums = [7,8,9,10,11,12]
    #samplenums = [7,8,9,10,11]
    samplenums = [21,20,19,18,16,14]
    #fns = ["../results/dr72_bin11_all.npy", "../results/dr72_bin11_pi2.npy"]
         #  "../results/dr72_bin11_pi1.npy", "../results/dr72_bin11_pi1corrfunc.npy"]
    #fns = ["../results/dr72_bin{}_rmin0.1_rmax50.npy".format(samplenum) for samplenum in samplenums]
    #fns = ['../results/dr72_bin11_all.npy', '../results/dr72_bin11_lsno2.npy']
    fns = ["../results/dr72_bin{}.npy".format(samplenum) for samplenum in samplenums]
    #fns = ['../results/dr72_bin7_all.npy', '../results/dr72_bin7_weights1.npy']
    #fns = ['../results/dr72_bin11_all.npy', '../results/dr72_bin11_rmin0.1.npy',
    #       '../results/dr72_bin11_rmin0.1_rmax50.npy', '../results/dr72_bin11_H0-100_Om0-0.3.npy']

    for fn in fns:

        #fn = '../results/dr72_bin8_frac0.05_try4.npy'
        rp, wprp, label = run.load_results(fn)
        print wprp
        #rps.append(rp[0])
        rps_paper = [0.17,0.27,0.42,0.67,1.1,1.7,2.7,4.2,6.7,10.6,16.9,26.8]
        rps_paper_all = [0.17,0.27,0.42,0.67,1.1,1.7,2.7,4.2,6.7,10.6,16.9,26.8,42.3]

        print len(rp[0])
        print len(rps_paper)
        if len(rp[0])==len(rps_paper):
            rps.append(rps_paper)
            print 'here'
        elif len(rp[0])==len(rps_paper_all):
            rps.append(rps_paper_all)
            print 'there'
        wprps.append(wprp[0])
        labels.append(label[0])
    plot_wprp(rps, wprps, labels)



def plot_bins():
    #samplenums = [7,8,9,10,11]
    samplenums = [11,10,9,8,7]
    #samplenums = [7]
    rps = []
    wps = []
    labels = []
    for samplenum in samplenums:
        fn = '../results/wprps/wprp_bin{}_square5k_toprobust.npy'.format(samplenum)
        hi = np.load(fn)
        print hi
        # print '0',hi[0]
        # print '1',hi[1]
        # print
        # rp = hi[2]
        # wprp = hi[0][0]
        # print rp
        # print wprp
        rp, wprp, label = np.load(fn)
        rps_paper = [0.17, 0.27, 0.42, 0.67, 1.1, 1.7, 2.7, 4.2, 6.7, 10.6, 16.9, 26.8]
        rps.append(rps_paper)
        wps.append(wprp[0])
        labels.append('bin'+str(samplenum))
    plot_wprp(rps, wps, labels)


def plot_wprp(rps, wprps, labels, colors=None, wp_tocompare=None, saveto=None):
    # if np.array(rps).ndim==1:
    #     rps = [rps]
    #     wprps = [wprps]
    #     labels = [labels]
    assert len(rps)==len(wprps) and len(wprps)==len(labels)
    color_idx = np.linspace(0, 1, len(wprps))
    if wp_tocompare:
        compidx = labels.index(wp_tocompare)
        fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
    else:
        fig = plt.figure()
        ax0 = fig.gca()

    for i in range(len(rps)):
        rp = rps[i]
        wprp = np.array(wprps[i])
        label = labels[i]
        if label in color_dict:
            color = color_dict[label]
        elif colors:
            color = colors[i]
        else:
            #color = color_list[i]
            color = plt.cm.rainbow(color_idx[i])
        # if label=='tophat':
        #     plt.step(rps, [0]+wprp, marker=None, label=basisname)
        #     plt.xscale('log')
        #     plt.yscale('log')
        # else:

        ax0.loglog(rp, wprp, label=label, color=color, marker='o', markersize=4, ls='-')
        #ax0.semilogx(rp, wprp, label=label, color=color, marker='o')

        plt.xlabel(r'$r_p$ (Mpc/h)')
        plt.xlim(0.1, 40.)
        ax0.set_ylim(1, 2000)
        ax0.set_ylabel(r'$w_p$($r_p$)')

        ax0.legend(loc='best')

        if wp_tocompare:
            # if wp.all() != wp_tocompare.all():
            wpcomp = wprps[compidx]

            if len(wprp)==len(wpcomp):
                #ax1.semilogx(rp, (np.log10(wprp)-np.log10(wpcomp)) / np.log10(wpcomp), color=colors[label])
                ax1.semilogx(rp, wprp/wpcomp, color=color)
                ax1.set_ylabel(r'$w_p$/$w_{{p,\mathrm{{{0}}}}}$'.format(wp_tocompare))

    if saveto:
        plt.savefig(saveto)
    plt.show()


def plot_xi_zspace(ss, xis, labels, colors=None, xi_tocompare=None, saveto=None):

    # if np.array(ss).ndim==1:
    #     ss = [ss]
    #     xis = [xis]
    #     labels = [labels]

    if xi_tocompare:
        compidx = labels.index(xi_tocompare)
        fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig = plt.figure()
        ax0 = fig.gca()

    for i in range(len(ss)):
        s = ss[i]
        xi = np.array(xis[i])
        label = labels[i]

        if label in color_dict:
            color = color_dict[label]
        elif colors:
            color = colors[i]
        else:
            color = color_list[i]

        ax0.plot(s, xi, label=label, color=color, marker='o', ls='-')

        #plt.xlabel(r'$s$ (Mpc/h)')
        #ax0.set_ylim(-0.005, 0.025)
        #plt.xlim(60, 210)
        ax0.set_ylabel(r'$\xi(s)$')

        ax0.legend(loc='best')

        if xi_tocompare:
            # if wp.all() != wp_tocompare.all():
            xicomp = xis[compidx]

            if len(xi)==len(xicomp):
                #ax1.semilogx(rp, (np.log10(wprp)-np.log10(wpcomp)) / np.log10(wpcomp), color=colors[label])
                ax1.semilogx(s, xi/xicomp, color=color)
                ax1.set_ylabel(r'$w_p$/$w_{{p,\mathrm{{{0}}}}}$'.format(xi_tocompare))

    if saveto:
        plt.savefig(saveto)
    plt.show()

if __name__=="__main__":
    main()