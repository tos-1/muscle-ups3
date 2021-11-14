""" Utility functions for analyzing results from muscle-ups """
import numpy as N
import struct
import pynbody
import os
import matplotlib.pyplot as plt
from ext._mas import mas
from ext._masHalos import masHalos

def getkgrid(boxsize, ng, what=None):
    '''
    It returns a meshgrid of kx,ky ,kz and of modulus k
    '''
    kmin = 2*N.pi/N.float(boxsize)
    thirdim = ng//2+1
    sh = (ng, ng, thirdim)
    kx, ky, kz = N.mgrid[0:sh[0], 0:sh[1], 0:sh[2]].astype(N.float64)

    kx[N.where(kx > ng//2)] -= ng
    ky[N.where(ky > ng//2)] -= ng
    kz[N.where(kz > ng//2)] -= ng

    kx *= kmin
    ky *= kmin
    kz *= kmin

    k = N.sqrt(kx**2+ky**2+kz**2)

    if what == None:
        return k
    elif what == 'all':
        return kx, ky, kz, k


def get_pos(pathtogadget="/home/federico/Documenti/PhD/LSS/HMUSCLE/sims/bx256_ng256_z50.0_Om0.30/Gadget/IC_2lpt_z50_008", flat=0):

    snap = pynbody.load(pathtogadget)
    #boxsize = snap.header.BoxSize
    ng = snap.header.npart[1]
    ng = int(round(ng**(1/3)))
    pos = snap['pos']
    iord = snap['iord']-1
    indices = N.argsort(iord)
    pos = pos[indices, :]

    if flat == 1:
        pass

    else:
        pos = N.reshape(pos, (ng, ng, ng, 3))
        pos = N.rollaxis(pos, -1)

    return pos


def computeStatistics(pathtobin=None, saveto=None):
    """ Compute statistics from binary """
    # if not os.path.exists(pathtobin):
    #        raise ValueError("Select correctly the path to binary")
    snap = pynbody.load(pathtobin)
    boxsize = snap.header.BoxSize
    ng = snap.header.npart[1]
    ng = int(round(ng**(1/3)))
    pos = snap['pos']
    iord = snap['iord']-1
    indices = N.argsort(iord)
    pos = pos[indices, :]
    pos = N.reshape(pos, (ng, ng, ng, 3))
    pos = N.rollaxis(pos, -1)

    folder = '/home/federico/Documenti/PhD/LSS/HMUSCLE/'

    kbins = N.load(folder+'data_plots/kbins.npy')
    PkNb = N.load(folder+'data_plots/PkNb.npy')

    Xk_moses = N.load(folder+'data_plots/Xk_hmsc_hm_opt.npy')
    Pk_moses = N.load(folder+'data_plots/Pk_hmsc_hm_opt.npy')
    Xk_alpt = N.load(folder+'data_plots/Xk_alpt.npy')
    Pk_alpt = N.load(folder+'data_plots/Pk_alpt.npy')

    Pk = fastPk(pos=pos, boxsize=256., kb=0)
    posNb = get_pos(pathtogadget=folder +
                    'sims/bx256_ng256_z50.0_Om0.30/Gadget/IC_2lpt_z50_008')
    Xk = fastXk(pos1=pos, pos2=posNb, boxsize=256., kb=0)

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(
        4.5*2., 4.5*2), gridspec_kw={'hspace': 0.0})

    for i in range(2):
        ax[i].grid(ls='--')
        ax[i].tick_params(direction='in', length=6, width=2, colors='k', grid_color='k',
                          grid_alpha=0.5, labelsize='x-large', which='both', top=True, right=True)
        ax[i].set_xlim([kbins[0], kbins[-1]])
        ax[i].patch.set_edgecolor('black')
        ax[i].patch.set_linewidth('1')

    ax[0].set_ylim([0.75, 1.05])
    ax[1].set_ylim([0.1, 1.15])

    ax[1].set_xlabel('$k \ [h/Mpc]$', fontsize='xx-large')
    ax[0].set_ylabel('$X(k)$', fontsize='xx-large')
    ax[1].set_ylabel('$P_{approx}(k)/P_{Nbody}(k)$', fontsize='xx-large')

    ax[0].semilogx(kbins, Xk, lw=2.5, label='current')
    ax[0].semilogx(kbins, Xk_moses, lw=2.5, label='muscleups')
    ax[0].semilogx(kbins, Xk_alpt, lw=2.5, label='alpt')

    ax[1].semilogx(kbins, Pk/PkNb, lw=2.5, label='current')
    ax[1].semilogx(kbins, Pk_moses/PkNb, lw=2.5, label='muscleups')
    ax[1].semilogx(kbins, Pk_alpt/PkNb, lw=2.5, label='alpt')

    fig.subplots_adjust(right=0.8)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles, labels, ncol=1, fontsize='xx-large',
                 loc='center', bbox_to_anchor=(1.15, 0.))
    if saveto is not None:
        plt.savefig(saveto, dpi=150)
    plt.show()

    return


def haloStatistics(pathtobin,ng,saveto=None):
    """ Compute statistics from binary """
    # if not os.path.exists(pathtobin):
    #        raise ValueError("Select correctly the path to binary")

    dmup = MAShalos(pathtobin, 256., ng) # muscle-ups halo field
    snap = pynbody.load(pathtobin)
    boxsize = snap.header.BoxSize
 
    dNb  = rockstarHalos(boxsize,ng)

    kbins, PkNb = fastPk(d=dNb, boxsize=boxsize, kb=1)

    Pkmup = fastPk(d=dmup, boxsize=boxsize)

    Xk = fastXk(d1=dNb, d2=dmup, boxsize=256., kb=0)

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=False, figsize=(
        4.5*2., 4.5*2), gridspec_kw={'hspace': 0.0})

    for i in range(2):
        ax[i].grid(ls='--')
        ax[i].tick_params(direction='in', length=6, width=2, colors='k', grid_color='k',
                          grid_alpha=0.5, labelsize='x-large', which='both', top=True, right=True)
        ax[i].set_xlim([kbins[0], kbins[-1]])
        ax[i].patch.set_edgecolor('black')
        ax[i].patch.set_linewidth('1')

    ax[0].set_ylim([0.75, 1.05])
    ax[1].set_ylim([0.1, 1.15])

    ax[1].set_xlabel('$k \ [h/Mpc]$', fontsize='xx-large')
    ax[0].set_ylabel('$X(k)$', fontsize='xx-large')
    ax[1].set_ylabel('$P_{approx}(k)/P_{Nbody}(k)$', fontsize='xx-large')

    ax[0].semilogx(kbins, Xk, lw=2.5, label='current')

    ax[1].semilogx(kbins, Pkmup/PkNb, lw=2.5, label='current')

    fig.subplots_adjust(right=0.8)
    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(handles, labels, ncol=1, fontsize='xx-large',
                 loc='center', bbox_to_anchor=(1.15, 0.))
    if saveto is not None:
        plt.savefig(saveto, dpi=150)
    plt.show()

    return


def MAS(pos, boxsize, ng):
    """ Mass assignment scheme. It uses CiC algorhitm as explained in Leclerque's thesis. It is calling an external C module, _mas,
        which simply assigns mass to density on grid.
    Inputs:
        pos: tuple (3,ng,ng,ng) of positions
        boxsize: side of the box of the sim
        np: number of grid resolution 
    """

    np = N.shape(pos[0])[0]
    sh = (ng, ng, ng,)
    dens = N.zeros(sh, dtype=N.float32).flatten()
    pos = N.asarray(pos).flatten().astype(dtype=N.float32)
    pos[pos == 256.] = 1e-06
    mas(ng, np, boxsize, pos, dens)
    dens = dens.reshape(sh)

    return dens


def MAShalos(path_to_binary, boxsize, ng):
    """ Mass assignment scheme for haloes. It uses CiC algorhitm as 
	explained in Leclerque's thesis. It is calling an external C module, _mas,
        which simply assigns mass to density on grid.
    Inputs:
        boxsize: side of the box of the sim
        ng: number of grid resolution 
    """
    path_to_catalogue = path_to_binary[:-4]+'_hc.dat'
    filename = os.path.basename(path_to_catalogue)
    path = os.path.dirname(path_to_catalogue)+'/'
    
    with open(path_to_catalogue, mode='rb') as file: # b is important -> binary
          fileContent = file.read()
    Nhalos = struct.unpack('i',fileContent[:4])[0]
    
    sizeofall=4#remove Nhalos
    size=N.zeros(Nhalos,dtype=N.int32)
    
    #extract the sizes first
    j=0
    while j<Nhalos:
        size[j] = struct.unpack('i',fileContent[sizeofall:sizeofall+4])[0]
        sizeofall = size[j]+sizeofall
        j+=1
    
    csize = N.append(0,size.cumsum())
    
    # size, haloID, haloCID, np, Mass, rv, cn, posHC[3], halopartID[]
    hc = [struct.unpack('iiiidddddd'+'i'*((size[j]-4*4-8*6)//4),fileContent[4+csize[j]:4+csize[j+1]]) for j in range(Nhalos)]
    #masses = N.asarray([_hc[4] for _hc in hc])
    Nparticles = N.asarray([_hc[3] for _hc in hc]).astype(dtype=N.int32)
    posHC  = N.asarray([_hc[7:10] for _hc in hc]) # position of center of mass
    posHCx, posHCy, posHCz = posHC[:,0], posHC[:,1], posHC[:,2]
    posHC = N.concatenate((posHCx, posHCy, posHCz),axis=0)
    #print(posHC.shape)

    sh = (ng, ng, ng,)
    densH = N.zeros(sh, dtype=N.float64).flatten()
    posHC = N.asarray(posHC).astype(dtype=N.float64)
    #print(posHC[0],posHC[1],posHC[2])
    #print(posHC[0],posHC[Nhalos],posHC[2*Nhalos]) # coincides with this
    posHC[posHC == 256.] = 1e-06
    masHalos( ng, Nhalos, boxsize, posHC, densH, Nparticles)
    densH = densH.reshape(sh)

    return densH


def fastXk(dk1=None, dk2=None, d1=None, d2=None, pos1=None, pos2=None, boxsize=None, gadgetfile=None, nkbins=None, kb=0):
    '''
    cross correlation between fields on the grid
    '''

    if N.any(d1 == None) and N.any(dk1 == None) and N.any(pos1 == None):
        print('provide only one among pos, d or dk')
        assert 0

    if N.any(d2 == None) and N.any(dk2 == None) and N.any(pos2 == None):
        print('provide only one among pos, d or dk')
        assert 0

    boxsize = float(boxsize)

    ############
    if N.any(pos1 != None):
        sh = N.shape(pos1[0])
        ng1 = sh[0]
        d1 = MAS(pos1, boxsize, ng1)
        dk1 = N.fft.rfftn(d1)
    elif N.any(d1 != None):
        sh = N.shape(d1)
        ng1 = sh[0]
        dk1 = N.fft.rfftn(d1)
    else:
        ng1 = N.shape(dk1)[0]
        sh = (ng1, ng1, ng1/2+1)

    if gadgetfile is not None:
        print('reading binary')
        pos2 = get_pos(gadgetfile, flat=0)

    ############
    if N.any(pos2 != None):
        sh = N.shape(pos2[0])
        ng2 = sh[0]
        d2 = MAS(pos2, boxsize, ng2)
        dk2 = N.fft.rfftn(d2)
    elif N.any(d2 != None):
        sh = N.shape(d2)
        ng2 = sh[0]
        dk2 = N.fft.rfftn(d2)
    else:
        ng2 = N.shape(dk2)[0]
        sh = (ng2, ng2, ng2//2+1)

    ############
    try:
        len(N.shape(dk1)) == len(N.shape(dk2))
        ng = N.shape(dk1)[0]
        sh = len(N.shape(dk1))
    except:
        raise ValueError('the input field must have the same size')

    kgrid = getkgrid(boxsize, ng)
    norm = (boxsize/float(ng)**2.)**sh

    if N.any(pos1 != None):
        kny = N.pi*ng/boxsize
        dk1 = dk1/N.sqrt(1.-N.sin(N.pi * kgrid/2./kny)**2. * 2./3.)

    if N.any(pos2 != None):
        kny = N.pi*ng/boxsize
        dk2 = dk2/N.sqrt(1.-N.sin(N.pi * kgrid/2./kny)**2. * 2./3.)

    if nkbins == None:
        nkbins = ng//2

    kgrid = kgrid.flatten()
    Xk = (N.conjugate(dk1) * dk2 + N.conjugate(dk2) * dk1).real.flatten()
    Xk *= 0.5
    Xk = Xk[kgrid > 0.]
    knz = kgrid[kgrid > 0.]

    delta = (ng*N.pi/boxsize - 2*N.pi/boxsize)/float(nkbins+1)/2.
    kbin = N.linspace(2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta, nkbins+1)
    xs = N.histogram(N.log10(knz), N.log10(kbin), weights=Xk, range=(
        2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta))[0]
    counts = N.histogram(N.log10(knz), N.log10(kbin), range=(
        2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta))[0]
    binvals = kbin[0:nkbins]+N.diff(kbin)/2.
    binvals = binvals[counts > 0]
    xs = xs[counts > 0]
    counts = counts[counts > 0]
    xs = xs/counts

    # if pk1==None:
    pk1 = fastPk(dk=dk1, boxsize=boxsize)

    # if pk2==None:
    pk2 = fastPk(dk=dk2, boxsize=boxsize)

    xs /= N.sqrt(pk1 * pk2)
    xs *= norm

    if kb == 0:
        return xs
    else:
        return binvals, xs


def fastPk(d=None, dk=None, pos=None, boxsize=None, nkbins=None, kb=0):
    '''
    Compute Pk for a field on the grid.
    Inputs::
        dk: fft of field. It can be a 2d field as well
        boxsize: of the simulation
        nkbins: number of bins in k-modulus space
        kb: boolean to decide whether to return (kb,ps)
        linear: boolean to decide whether the bin width is linear
    '''

    if N.any(d == None) and N.any(dk == None) and N.any(pos == None):
        print('provide only one among pos, d or dk')
        assert 0

    boxsize = float(boxsize)

    if N.any(pos != None):
        sh = N.shape(pos[0])
        ng = sh[0]
        d = MAS(pos, boxsize, ng)
        dk = N.fft.rfftn(d)

    elif N.any(d != None):
        sh = N.shape(d)
        ng = sh[0]
        dk = N.fft.rfftn(d)
    else:
        ng = N.shape(dk)[0]
        sh = (ng, ng, ng//2+1)

    #ng = N.shape(dk)[0]
    sh = len(N.shape(dk))  # dimensions
    kgrid = getkgrid(boxsize, ng)
    norm = (boxsize/float(ng)**2.)**sh

    if N.any(pos != None):
        kny = N.pi*ng/boxsize
        dk = dk/N.sqrt(1.-N.sin(N.pi * kgrid/2./kny)**2. * 2./3.)

    kgrid = kgrid.flatten()
    dk2 = abs(dk.flatten()) ** 2.
    dk2 = dk2[kgrid > 0.]
    knz = kgrid[kgrid > 0.]

    if nkbins == None:
        nkbins = ng//2

    delta = (ng*N.pi/boxsize - 2*N.pi/boxsize)/float(nkbins+1)/2.
    kbin = N.linspace(2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta, nkbins+1)
    ps = N.histogram(N.log10(knz), N.log10(kbin), weights=dk2,
                     range=(2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta))[0]
    counts = N.histogram(N.log10(knz), N.log10(kbin), range=(
        2*N.pi/boxsize-delta, ng*N.pi/boxsize-delta))[0]
    binvals = kbin[0:nkbins]+N.diff(kbin)/2.
    binvals = binvals[counts > 0]
    ps = ps[counts > 0]
    counts = counts[counts > 0]
    ps = ps/counts

    if kb == 0:
        return ps * norm
    else:
        return binvals, ps * norm

def rockstarHalos(boxsize,ng):
    folder = '/home/federico/Documenti/PhD/LSS/HMUSCLE/sims/bx256_ng256_z50.0_Om0.30/Gadget/rockstar_halos/'
    from BGC2 import read_bgc2
    halos = N.empty((1,14))
    for num in range(8):
        input_file = folder+'halos_0.'+str(num)+'.bgc2'
        _ , _halos, _ = read_bgc2( input_file )
        halos = N.concatenate((halos,_halos),axis=0)
    halos = halos[1:] # remove empty
    
    #halo-id, parent-id, nparticles, _, rvir, mass, pos[3], vel[3]
    
    #halocat = N.empty(nhalos, dtype=[('Position', ('f8', 3)), ('Velocity', ('f8', 3)), ('Np', 'i4')])
    Nparticles = N.asarray(halos[:,2],dtype=N.int32)
    posHC = N.asarray(halos[:,6:9],dtype=N.float64)
    del halos
    import gc
    gc.collect()

    posHC = posHC[Nparticles>20,:]
    Nparticles = Nparticles[Nparticles>20]

    Nhalos = len(Nparticles)
    print('total number of haloes',Nhalos)
    posHCx, posHCy, posHCz = posHC[:,0], posHC[:,1], posHC[:,2]
    posHC = N.concatenate((posHCx, posHCy, posHCz),axis=0)
    #print(posHC.shape)

    sh = (ng, ng, ng,)
    densH = N.zeros(sh, dtype=N.float64).flatten()
    posHC = N.asarray(posHC).astype(dtype=N.float64)
    #print(posHC[0],posHC[1],posHC[2])
    #print(posHC[0],posHC[Nhalos],posHC[2*Nhalos]) # coincides with this
    posHC[posHC == 256.] = 1e-06
    masHalos( ng, Nhalos, boxsize, posHC, densH, Nparticles)
    densH = densH.reshape(sh)

    return densH


def generateHMSC(ng=256, boxsize=256., get_pk=False, info='', window='gauss', hm='PS', sigmaalpt=3.0, redshift=0.0):
    """ generate realizations with HaloMUSCLE """
    HMSC = lpt.lpt(scheme='hmuscle', z_pk=50., redshift=redshift, ng=ng, boxsize=boxsize, get_pk=get_pk,
                   return_pos=False, aug=True, extra_info=info, window=window, hm=hm, sigmaalpt=sigmaalpt)
    HMSC.generate()
    return


def computeHMF(path_to_hc, saveto=None):
    """ return the halo mass function from the halo catalogue """
    boxsize = 256.
    ng = 256
    nph = 0
    z = 0.

    from colossus.cosmology import cosmology
    params = {'flat': True, 'H0': 70., 'Om0': 0.3,
              'Ob0': 0.046, 'sigma8': 0.8, 'ns': 0.96}
    colossus = cosmology.setCosmology('myCosmo', params)
    from colossus.lss import mass_function

    # matched HMF
    with open('hmf_at_z0.0_Om0.30.dat', mode='rb') as file:
        fileContent = file.read()
        Nbinedg = struct.unpack('i', fileContent[:4])[0]
        cumsize = 4  # cumulative size (bytes)
        logM = N.zeros(Nbinedg, dtype=N.float64)
        t08 = N.zeros(Nbinedg, dtype=N.float64)
        j = 0
        while j < Nbinedg:
            logM[j] = struct.unpack('d', fileContent[cumsize:cumsize+8])[0]
            cumsize += 8
            j += 1
        while j < Nbinedg*2:
            t08[j-Nbinedg] = struct.unpack('d', fileContent[cumsize:cumsize+8])[0]
            cumsize += 8
            j += 1

    with open(path_to_hc[:-4]+'_hc.dat', mode='rb') as file:
        fileContent = file.read()
        Nhalos = struct.unpack('i', fileContent[:4])[0]
        sizeofall = 4  # remove Nhalos
        size = N.zeros(Nhalos, dtype=N.int32)
        j = 0
        while j < Nhalos:
            size[j] = struct.unpack('i', fileContent[sizeofall:sizeofall+4])[0]
            sizeofall = size[j]+sizeofall
            j += 1

    csize = N.append(0, size.cumsum())
    hc = [struct.unpack('iiiidddddd'+'i'*((size[j]-4*4-6*8)//4),
                        fileContent[4+csize[j]:4+csize[j+1]]) for j in range(Nhalos)]
    np = [hc[i][2] for i in range(Nhalos)]
    masses = N.asarray([hc[i][2] for i in range(Nhalos)])
    masses = masses.astype(N.float32)
    # pygadgetreader.readheader('sims/bx256_ng256_z50.0_Om0.30/Gadget/IC_2lpt_z50_008','O0')
    O0 = 0.3
    mass_res = 2.77536627e+11 * O0 * (boxsize/float(ng))**3.
    masses *= mass_res

    logbins = N.linspace(12., 15., int((15.-12.)/0.3))
    nM, bindg = N.histogram(N.log10(masses), bins=logbins)
    dlogM = bindg[1]-bindg[0]
    nM = nM/boxsize**3
    nM = nM/N.log(10**dlogM)
    binmean = (bindg[1:]+bindg[:-1])/2.
    binmean = binmean[nM > 0]
    nM = nM[nM > 0]
    nM = nM[10**binmean/mass_res >= nph]
    binmean = binmean[10**binmean/mass_res >= nph]
    hmf = N.vstack((binmean, nM))

    t08 = N.interp(hmf[0], logM, t08)
    #t08 = mass_function.massFunction( 10**hmf[0], z, mdef = '340m', model = 'tinker08', q_out = 'dndlnM')

    fig, ax = plt.subplots(figsize=(6*1.62, 6*1.3), sharey=True)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size='35%', pad=0.)
    ax.figure.add_axes(ax2)

    ax.plot(hmf[0], hmf[1], lw=2.5)  # ,label='top-hat',lw=2.5)
    ax2.plot(hmf[0], hmf[1]/t08-1., lw=2.5)
    ax.semilogy(hmf[0], t08, label='Tinker+08', lw=2.5, ls='--', color='black')
    ax2.set_ylim([-1., 1.])
    ax2.set_yticks([-1., -0.5, 0., 0.5, 1.])
    ax.legend(fontsize='xx-large', frameon=False)
    ax.set_ylabel("$d n/d \ln{M} \ [Mpc/h]^3$", fontsize='xx-large')
    ax2.set_ylabel("$\Delta$", fontsize='xx-large')
    ax2.set_xlabel("$\log_{10}{(m/(M_\odot / h))}$", fontsize='xx-large')
    ax2.hlines(0., hmf[0, 0], hmf[0, -1], colors='black', linestyles='--')
    for axx in [ax, ax2]:
        axx.set_xlim([hmf[0, 0], hmf[0, -1]])
        axx.grid(ls='--')
        axx.tick_params(direction='in', length=6, width=2, colors='k', grid_color='k',
                        grid_alpha=0.5, labelsize='x-large', which='both', top=True, right=True)
        axx.patch.set_edgecolor('black')
        axx.patch.set_linewidth('1')
    if saveto is not None:
        plt.savefig(saveto, dpi=150, bbox_inches='tight')
    plt.show()
    return


def lagslice(pathtobin=None,saveto='sims/bx256_ng256_z0.0_Om0.30/images/'):
    """ Eulerian positions snapshot, color coding halo particles """

    if not os.path.exists(pathtobin):
            raise ValueError("Select correctly the path to binary")

    snap = pynbody.load(pathtobin)
    boxsize = snap.header.BoxSize
    ng = snap.header.npart[1]
    ng = int( round( ng**(1./3.) ) )
    pos = snap['pos']
    iord = snap['iord']-1
    indices = N.argsort(iord)
    pos = pos[indices,:]
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)

    # import halonum
    with open( pathtobin[:-4]+'_hnum.dat', mode='rb') as file:
        fileContent = file.read()
        halonum = N.zeros(ng**3,dtype=N.int32)
        sizeofall=0
        j=0
        while j<ng**3:
            halonum[j] = struct.unpack('i',fileContent[sizeofall:sizeofall+4])[0]
            sizeofall +=4
            j+=1

    haloIDh = N.where(N.logical_and(halonum!=-1,halonum!=ng**3)==1)[0]

    kk = haloIDh%(ng)
    jj = ((haloIDh-kk)//ng)%ng
    ii = (haloIDh-jj*ng-kk)//ng//ng
    slicek = N.where( kk==ng//2 )

    all_id = halonum[ii[slicek]*ng*ng+jj[slicek]*ng+kk[slicek]]

    vals = N.linspace(0,1,256*4)
    N.random.shuffle(vals)
    cmap = plt.cm.colors.ListedColormap(plt.cm.jet(vals))

    fig,ax = plt.subplots(figsize=(7,7))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    ax.scatter(pos[0][:,:,ng//2].flatten(),pos[1][:,:,ng//2].flatten(),s=.2,lw=0.,color='black')
    ax.scatter(pos[0][ii[slicek],jj[slicek],kk[slicek]].flatten(),pos[1][ii[slicek],jj[slicek],kk[slicek]].flatten(),s=1.,lw=0.,c=all_id,cmap=cmap)
    ax.set_xlim([0,ng])
    ax.set_ylim([0,ng])
    if saveto is not None:
        plt.savefig(saveto,dpi=150,bbox_inches='tight')
    plt.show()
    return

def muscleups_particles(ng=256,boxsize=256., path_to_hc=None,saveto=None):
    """ plot lagrangian slice of particles color coded by haloID """

    with open( path_to_hc[:-4]+'_hnum.dat', mode='rb') as file:
        fileContent = file.read()
        halonum = N.zeros(ng**3,dtype=N.int32)
        sizeofall=0
        j=0
        while j<ng**3: #j<ng*(ng+16)**2:
            halonum[j] = struct.unpack('i',fileContent[sizeofall:sizeofall+4])[0]
            sizeofall +=4
            j+=1

    vals = N.linspace(0,1,256*4)
    r = N.random.RandomState(1234)
    r.shuffle(vals)
    cmap = plt.cm.colors.ListedColormap(plt.cm.jet(vals))
    fig,ax = plt.subplots(figsize=(10,10))
    ax.set_xlabel('x [Mpc/h]',fontsize='xx-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='xx-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',
    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    extent = 0, ng, 0,ng
    halonum = halonum.reshape(ng,ng,ng)
    #fig.tight_layout(rect=[0., 0., 0.97, 0.9])
    ax.tick_params(direction='in', length=6, width=2, colors='k',
    grid_color='k', grid_alpha=0.5,labelsize='x-large')
    masked_array = N.ma.masked_where( N.logical_or(halonum==ng**3, halonum==-1), halonum )
    cmap.set_bad(color='white')
    ax.set_title('serial protohaloes', fontsize='x-large')
    im = ax.imshow(masked_array[:,:,ng//2], cmap=cmap,origin='lower')

    if saveto is not None:
        plt.savefig(saveto,bbox_inches='tight',dpi=150)
    else:
        plt.show()
    return


def dispfield(pathtobin=None,saveto='sims/bx256_ng256_z0.0_Om0.30/images/'):
    """ imshow of scalar displacement field """
    if not os.path.exists(pathtobin):
            raise ValueError("Select correctly the path to binary")

    snap = pynbody.load(pathtobin)
    boxsize = snap.header.BoxSize
    ng = snap.header.npart[1]
    ng = int( round( ng**(1./3.) ) )
    pos = snap['pos']
    iord = snap['iord']-1
    indices = N.argsort(iord)
    pos = pos[indices,:]
    pos = N.reshape(pos,(ng,ng,ng,3))
    pos = N.rollaxis(pos,-1)

    psi = utils.get_psi(pos=pos,boxsize=boxsize)
    #psi[psi>6.]=8.
    #psi[psi<-6.]=-8.
    fig,ax = plt.subplots(nrows=1,ncols=1,figsize=(6,6))
    #ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.set_xlabel('x [Mpc/h]',fontsize='x-large')
    ax.set_ylabel('y [Mpc/h]',fontsize='x-large')
    ax.tick_params(direction='in', length=6, width=2, colors='k',grid_color='k', grid_alpha=0.5,labelsize='large')
    ax.grid(color='blue', linestyle='-', linewidth=0.4)
    ax.set_xlim([0,ng])
    ax.set_ylim([0,ng])
    im = ax.imshow(psi[:,:,ng/2],origin='lower')

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "7%", pad="-0.1%")
    cb=plt.colorbar(im, cax=cax,)
    cb.ax.tick_params(direction='in', length=6, width=2, colors='k',
    grid_color='k', grid_alpha=0.5,labelsize='xx-large')
    if saveto is not None:
        plt.savefig(saveto,dpi=150,bbox_inches='tight')
    plt.show()

    return



if __name__ == "__main__":
    #dh = MAShalos('sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat', 256., 32)

    #dh = rockstarHalos(256.,32)

    #haloStatistics('sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat',20,saveto=None)

    #computeStatistics( pathtobin='sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat', saveto=None)

    #computeHMF(path_to_hc='sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat', saveto=None)

    #lagslice(pathtobin='sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat', saveto=None)

    #muscleups_particles(ng=256,boxsize=256., path_to_hc='sims/bx256.0_ng256_z0.0_Om0.30/muscleups/z0.0__0.dat',saveto=None)
