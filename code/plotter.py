from matplotlib import pyplot as plt
import numpy as np

import run


color_dict = {'corrfunc':'grey', 'tophat':'blue', 'piecewise':'red',
          'tophat_orig':'cyan', 'gaussian':'green', 'trig':'magenta',
          'treecorr': 'orange', 'top_quad':'purple'}

color_list = ['red', 'orange', 'green', 'blue', 'cyan', 'magenta']

def main():
    rps = []
    wprps = []
    labels = []
    #samplenums = [7,8,9,10,11,12]
    #fns = ["../results/dr72_bin{}_all.npy".format(samplenum) for samplenum in samplenums]
    fns = ['../results/dr72_bin7_all.npy', '../results/dr72_bin7_weights1.npy']
    for fn in fns:

        #fn = '../results/dr72_bin8_frac0.05_try4.npy'
        rp, wprp, label = run.load_results(fn)
        rps.append(rp[0])
        wprps.append(wprp[0])
        labels.append(label[0])
    plot_wprp(rps, wprps, labels)

def plot_wprp(rps, wprps, labels, colors=None, wp_tocompare=None):

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


if __name__=="__main__":
    main()