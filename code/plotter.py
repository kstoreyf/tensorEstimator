from matplotlib import pyplot as plt
import numpy as np




cols = {'corrfunc':'grey', 'tophat':'blue', 'piecewise':'red',
          'tophat_orig':'cyan', 'gaussian':'green', 'trig':'magenta',
          'treecorr': 'yellow'}


def plot_wprp(rps, wprps, labels, colors=None, wp_tocompare=None):

    if wp_tocompare:
        compidx = labels.index(wp_tocompare)

    fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    for i in range(len(labels)):
        print rps
        print i
        print rps[i]
        rp = rps[i]
        wprp = np.array(wprps[i])
        label = labels[i]
        if colors:
            color = colors[i]
        else:
            colors = cols[label]
        # if label=='tophat':
        #     plt.step(rps, [0]+wprp, marker=None, label=basisname)
        #     plt.xscale('log')
        #     plt.yscale('log')
        # else:
        ax0.loglog(rp, wprp, label=label, color=color, marker='o')

        plt.xlabel(r'$r_p$ (Mpc/h)')
        ax0.set_ylabel(r'$w_p$($r_p$)')
        ax1.set_ylabel(r'$w_p$/$w_{{p,\mathrm{{{0}}}}}$'.format(wp_tocompare))

        ax0.legend(loc='best')

        if wp_tocompare:
            # if wp.all() != wp_tocompare.all():
            wpcomp = wprps[compidx]

            if len(wprp)==len(wpcomp):
                #ax1.semilogx(rp, (np.log10(wprp)-np.log10(wpcomp)) / np.log10(wpcomp), color=colors[label])
                ax1.semilogx(rp, wprp/wpcomp, color=color)

    plt.show()