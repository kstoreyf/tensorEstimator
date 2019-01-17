from matplotlib import pyplot as plt
import numpy as np

import run


color_dict = {'corrfunc':'grey', 'tophat':'blue', 'piecewise':'red',
          'tophat_orig':'cyan', 'gaussian':'green', 'trig':'magenta',
          'treecorr': 'orange', 'top_quad':'purple'}

color_list = ['red', 'orange', 'green', 'blue', 'cyan', 'magenta']

def main():
    plot_dr72bins()
    #plot()

def plot():
    #fn = "../results/dr72_brightLRG_frac0.1.npy"
    fn = "../results/dr7_FullLRG_frac1_weights.npy"
    ss, xis, labels = run.load_results(fn)
    print ss
    print xis
    print labels
    plot_xi_zspace(ss, xis, labels)


def plot_dr72bins():
    rps = []
    wprps = []
    labels = []
    #samplenums = [7,8,9,10,11,12]
    #samplenums = [7,8,9,10,11]
    fns = ["../results/dr72_bin11_frac0.1_pi1.npy", "../results/dr72_bin11_frac0.1_pi2.npy",
           "../results/dr72_bin11_frac0.1_pi2bad.npy"]#], "../results/dr72_bin11_frac0.1_pi1corrfunc.npy",]
    #fns = ["../results/dr72_bin{}_rmin0.1_rmax50.npy".format(samplenum) for samplenum in samplenums]
    #fns = ['../results/dr72_bin11_all.npy', '../results/dr72_bin11_lsno2.npy']
    #fns = ["../results/dr72_bin{}_all.npy".format(samplenum) for samplenum in samplenums]
    #fns = ['../results/dr72_bin7_all.npy', '../results/dr72_bin7_weights1.npy']
    #fns = ['../results/dr72_bin11_all.npy', '../results/dr72_bin11_rmin0.1.npy',
    #       '../results/dr72_bin11_rmin0.1_rmax50.npy', '../results/dr72_bin11_H0-100_Om0-0.3.npy']

    for fn in fns:

        #fn = '../results/dr72_bin8_frac0.05_try4.npy'
        rp, wprp, label = run.load_results(fn)
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




def plot_wprp(rps, wprps, labels, colors=None, wp_tocompare=None):
    if np.array(rps).ndim==1:
        rps = [rps]
        wprps = [wprps]
        labels = [labels]

    if wp_tocompare:
        compidx = labels.index(wp_tocompare)
        fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
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
            color = color_list[i]
        # if label=='tophat':
        #     plt.step(rps, [0]+wprp, marker=None, label=basisname)
        #     plt.xscale('log')
        #     plt.yscale('log')
        # else:
        ax0.loglog(rp, wprp, label=label, color=color, marker='o')
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

    plt.show()


def plot_xi_zspace(ss, xis, labels, colors=None, xi_tocompare=None):

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

        ax0.plot(s, xi, label=label, color=color, marker='o')

        plt.xlabel(r'$s$ (Mpc/h)')
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

    plt.show()

if __name__=="__main__":
    main()