from make_statistics import *
import numpy as N
import matplotlib.pyplot as plt

def haloStatistics(pathtobin, boxsize, ng, saveto=None):
    """ Compute statistics for dark matter from binary """

    dNb = MASHalos(pathtobin, boxsize, ng)

    kbins, PkNb = fastPk(d=dNb, boxsize=boxsize, kb=1)

    #Pkmup = fastPk(d=dmup, boxsize=boxsize)

    Xk = fastXk(d1=dNb, d2=dmup, boxsize=boxsize, kb=0)

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(
        4.5 * 2., 4.5 * 2), gridspec_kw={'hspace': 0.0})

    for i in range(2):
        ax[i].grid(ls='--')
        ax[i].tick_params(direction='in', length=6, width=2, colors='k', grid_color='k',
                          grid_alpha=0.5, labelsize='x-large', which='both', top=True, right=True)
        ax[i].set_xlim([kbins[0], kbins[-1]])
        ax[i].patch.set_edgecolor('black')
        ax[i].patch.set_linewidth('1')

    ax[0].set_ylim([0.75, 1.05])
    ax[1].set_ylim([0.1, 1.15])

    ax[1].set_xlabel('$k \\ [h/Mpc]$', fontsize='xx-large')
    ax[0].set_ylabel('$X(k)$', fontsize='xx-large')
    ax[1].set_ylabel('$P_{approx}(k)/P_{Nbody}(k)$', fontsize='xx-large')

    ax[0].semilogx(kbins, Xk, lw=2.5, label='current')

    ax[1].loglog(kbins, PkNb, lw=2.5, label='current')

    fig.subplots_adjust(right=0.8)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles, labels, ncol=1, fontsize='xx-large',
                 loc='center', bbox_to_anchor=(1.15, 0.))
    if saveto is not None:
        plt.savefig(saveto, dpi=150)
    plt.show()

    return

if __name__ == "__main__":
    haloStatistics('/home/federico/Quijote/snapdir_004/', 1000., 32, saveto=None)
