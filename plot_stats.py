import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_statistics import *
import numpy as N

boxsize = 1000.
ng = 512

def ComputeMatterStats():
    global ng, boxsize

    print("Quijote")
    snaps_quijote = '/home/wp3i/Quijote/10000/snapdir_004/snap_004'
    posQ = get_pos(snaps_quijote, norma=1e+03)
    dQ = MAS(posQ, boxsize, ng//2)
    kbins, pkQ = fastPk(d=dQ, boxsize=boxsize, kb=1)
    N.save('data/pkQ',pkQ)
    N.save('data/kbins',kbins)
 
    print("2lpt")
    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/2lpt/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    pk = fastPk(d=d, boxsize=boxsize, kb=0)
    xk = fastXk(d1=dQ, d2=d, boxsize=boxsize, kb = 0)
    N.save('data/pk_2lpt',pk)
    N.save('data/xk_2lpt',xk)

    print("alpt")
    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/alpt/z0.0sigmaalpt4.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    pk = fastPk(d=d, boxsize=boxsize, kb=0)
    xk = fastXk(d1=dQ, d2=d, boxsize=boxsize, kb = 0)
    N.save('data/pk_alpt',pk)
    N.save('data/xk_alpt',xk)

    print("Rockstar + muscleups")
    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/Rockstar/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    pk = fastPk(d=d, boxsize=boxsize, kb=0)
    xk = fastXk(d1=dQ, d2=d, boxsize=boxsize, kb = 0)
    N.save('data/pk_rock_mup',pk)
    N.save('data/xk_rock_mup',xk)

    return

def PlotMatterSpectra():
    kbins = N.load('data/kbins.npy')
    xk_rck = N.load('data/xk_rock_mup.npy')
    xk_2lpt = N.load('data/xk_2lpt.npy')
    xk_alpt = N.load('data/xk_alpt.npy')
    pkQ = N.load('data/pkQ.npy')
    pk_rck = N.load('data/pk_rock_mup.npy')
    pk_2lpt = N.load('data/pk_2lpt.npy')
    pk_alpt = N.load('data/pk_alpt.npy')
 
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

    ax[0].semilogx(kbins, xk_rck, lw=2.5, label='muscleups (rockstar)')
    ax[0].semilogx(kbins, xk_2lpt, lw=2.5, label='2lpt')
    ax[0].semilogx(kbins, xk_alpt, lw=2.5, label='alpt')

    ax[1].loglog(kbins, pk_rck/pkQ, lw=2.5, label='muscleups (rockstar)')
    ax[1].loglog(kbins, pk_2lpt/pkQ, lw=2.5, label='2lpt')
    ax[1].loglog(kbins, pk_alpt/pkQ, lw=2.5, label='alpt')

    fig.subplots_adjust(right=0.8)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles, labels, ncol=1, fontsize='xx-large',
                 loc='center', bbox_to_anchor=(1.15, 0.))
    plt.savefig('images/spectra.png', dpi=150)

    return

def plot_dens():
    """ Eulerian density snapshot """
    NG = ng//2
    extent = 0.0, boxsize, 0.0, boxsize
    snaps_quijote = '/home/wp3i/Quijote/10000/snapdir_004/snap_004'
    pos = get_pos(snaps_quijote, norma=1e+03)
    d = MAS(pos, boxsize, NG)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
            grid_color='k', grid_alpha=0.5,labelsize='x-large')
    im = ax.imshow(N.log(2.+d[NG//2,:,:]),extent=extent)
    plt.colorbar(im)
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('Quijote density', fontsize='x-large')
    plt.savefig('images/density_quijote.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/2lpt/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
            grid_color='k', grid_alpha=0.5,labelsize='x-large')
    im = ax.imshow(N.log(2.+d[NG//2,:,:]),extent=extent)
    plt.colorbar(im)
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('2LPT density', fontsize='x-large')
    plt.savefig('images/density_2lpt.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/alpt/z0.0sigmaalpt4.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
	    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    im = ax.imshow(N.log(2.+d[NG//2,:,:]),extent=extent)
    plt.colorbar(im)
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('ALPT density', fontsize='x-large')
    plt.savefig('images/density_alpt.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/Rockstar/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    d = MAS(pos, boxsize, ng//2)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
	    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    im = ax.imshow(N.log(2.+d[NG//2,:,:]),extent=extent)
    plt.colorbar(im)
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('rockstar + muscleups density', fontsize='x-large')
    plt.savefig('images/density_rockstar.png',dpi=150)

    return


def plot_snaps():

    NG = ng
    snaps_quijote = '/home/wp3i/Quijote/10000/snapdir_004/snap_004'
    pos = get_pos(snaps_quijote, norma=1e+03)
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
            grid_color='k', grid_alpha=0.5,labelsize='x-large')
    ax.scatter(pos[0][:,:,NG//2].flatten(),pos[1][:,:,NG//2].flatten(),s=.2,lw=0.,color='black')
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('Quijote', fontsize='x-large')
    plt.savefig('images/lagslice_quijote.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/2lpt/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
            grid_color='k', grid_alpha=0.5,labelsize='x-large')
    ax.scatter(pos[0][:,:,NG//2].flatten(),pos[1][:,:,NG//2].flatten(),s=.2,lw=0.,color='black')
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('2LPT', fontsize='x-large')
    plt.savefig('images/lagslice_2lpt.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/alpt/z0.0sigmaalpt4.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
	    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    ax.scatter(pos[0][:,:,NG//2].flatten(),pos[1][:,:,NG//2].flatten(),s=.2,lw=0.,color='black')
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('ALPT', fontsize='x-large')
    plt.savefig('images/lagslice_alpt.png',dpi=150)


    sim = '/home/wp3i/Quijote/muscleups/sims/bx1000.0_ng512_z0.0_Om0.30/Rockstar/z0.0__0.dat'
    pos = get_pos(sim, norma=1.0)
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
	    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    ax.scatter(pos[0][:,:,NG//2].flatten(),pos[1][:,:,NG//2].flatten(),s=.2,lw=0.,color='black')
    ax.grid(ls='--',color='blue')
    ax.set_xlim([0,boxsize])
    ax.set_ylim([0,boxsize])
    ax.set_title('rockstar + muscleups', fontsize='x-large')
    plt.savefig('images/lagslice_rockstar.png',dpi=150)
    return


#def haloStatistics(pathtobin, boxsize, ng, saveto=None):
#    """ Compute statistics for dark matter from binary """
#
#    dNb = MASHalos(pathtobin, boxsize, ng)
#
#    kbins, PkNb = fastPk(d=dNb, boxsize=boxsize, kb=1)
#
#    #Pkmup = fastPk(d=dmup, boxsize=boxsize)
#
#    Xk = fastXk(d1=dNb, d2=dmup, boxsize=boxsize, kb=0)
#
#    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(
#        4.5 * 2., 4.5 * 2), gridspec_kw={'hspace': 0.0})
#
#    for i in range(2):
#        ax[i].grid(ls='--')
#        ax[i].tick_params(direction='in', length=6, width=2, colors='k', grid_color='k',
#                          grid_alpha=0.5, labelsize='x-large', which='both', top=True, right=True)
#        ax[i].set_xlim([kbins[0], kbins[-1]])
#        ax[i].patch.set_edgecolor('black')
#        ax[i].patch.set_linewidth('1')
#
#    ax[0].set_ylim([0.75, 1.05])
#    ax[1].set_ylim([0.1, 1.15])
#
#    ax[1].set_xlabel('$k \\ [h/Mpc]$', fontsize='xx-large')
#    ax[0].set_ylabel('$X(k)$', fontsize='xx-large')
#    ax[1].set_ylabel('$P_{approx}(k)/P_{Nbody}(k)$', fontsize='xx-large')
#
#    ax[0].semilogx(kbins, Xk, lw=2.5, label='current')
#
#    ax[1].loglog(kbins, PkNb, lw=2.5, label='current')
#
#    fig.subplots_adjust(right=0.8)
#    handles, labels = ax[0].get_legend_handles_labels()
#    ax[0].legend(handles, labels, ncol=1, fontsize='xx-large',
#                 loc='center', bbox_to_anchor=(1.15, 0.))
#    if saveto is not None:
#        plt.savefig(saveto, dpi=150)
#    plt.show()
#
#    return

if __name__ == "__main__":
    #ComputeMatterStats()
    #PlotMatterSpectra()
    #plot_dens()
    plot_snaps()
    #haloStatistics('/home/federico/Quijote/snapdir_004/', 1000., 32, saveto=None)
