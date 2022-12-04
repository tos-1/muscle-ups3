""" Utility functions for analyzing results from muscle-ups """
import numpy as N
import os
from ext._mas import mas
from ext._masHalos import masHalos
from utils.BGC2 import read_bgc2
import readgadget
import gc

def get_pos(snapshot, norma=1e+03, binic=None):
    ''' get the positions of particles from binaries
	Inputs::
	  snapshot: do not specify format, only root of the name, e.g. snap_004 for Quijote
	  norma: the units of Gadget positions, either Mpc (norma=1), or kpc (norma=1e+03)
          binic: flag of Quijote simulations
    '''
    # read the positions and IDs of the ICs
    print('normalization of Gadget units=',norma)
    ptype = [1]    #1(CDM) 2(neutrinos)
    pos = readgadget.read_block(snapshot, "POS ", ptype)/norma #Mpc/h
    IDs = readgadget.read_block(snapshot, "ID  ", ptype)-1     #IDs begin from 0
    indexes = N.argsort(IDs)
    pos = pos[indexes]

    if binic is not None:
        ng = 512
        boxsize = 1000.
        pos_ICs = readgadget.read_block(
            binic, "POS ", ptype) / 1e3  # Mpc/h
        IDs_ICs = readgadget.read_block(
            binic, "ID  ", ptype) - 1  # IDs begin from 0
        indexes = N.argsort(IDs_ICs)
        pos_ICs = pos_ICs[indexes]

        grid_index = ( N.round( (pos_ICs / boxsize) * ng,
                      decimals=0 ) ).astype(N.int32)
        del pos_ICs, IDs_ICs
        gc.collect()
        
        grid_index[N.where(grid_index == ng)] = 0
        grid_index = grid_index[:, 0] * ng**2 + \
            grid_index[:, 1] * ng + grid_index[:, 2]
        
        indexes = N.argsort(grid_index)
        pos = pos[indexes]
    return pos

def get_pos_ics(binic, ng=512, boxsize=1000., norma=1e+03):
    ''' get the initial positions of particles from binaries
	Inputs::
	  snapshot: do not specify format, only root of the name, e.g. snap_004 for Quijote
	  norma: the units of Gadget positions, either Mpc (norma=1), or kpc (norma=1e+03)
          binic: flag of Quijote simulations
    '''
    # read the positions and IDs of the ICs
    print('normalization of Gadget units=',norma)
    ptype = [1]    #1(CDM) 2(neutrinos)
    pos_ICs = readgadget.read_block(
        binic, "POS ", ptype) / 1e3  # Mpc/h
    IDs_ICs = readgadget.read_block(
        binic, "ID  ", ptype) - 1  # IDs begin from 0
    indexes = N.argsort(IDs_ICs)
    pos_ICs = pos_ICs[indexes]

    grid_index = ( N.round( (pos_ICs / boxsize) * ng,
                  decimals=0 ) ).astype(N.int32)
    gc.collect()
    
    grid_index[N.where(grid_index == ng)] = 0
    grid_index = grid_index[:, 0] * ng**2 + \
        grid_index[:, 1] * ng + grid_index[:, 2]

    indexes = N.argsort(grid_index)
    pos_ICs = pos_ICs[indexes]
    return pos_ICs

def get_psi( disp, boxsize=1000.):

    disp[N.where(disp> boxsize/2.0)] -= boxsize
    disp[N.where(disp<-boxsize/2.0)] += boxsize

    ng = int(N.cbrt(N.shape(disp)[0]))
    cellsize = boxsize / ng
    disp = N.reshape(disp,(ng,ng,ng,3))
    disp = N.rollaxis(disp,-1)

    # partial derivatives
    for i in range(3):
        disp[i] = (2. / 3. * (N.roll(disp[i], -1, i) - N.roll(disp[i], 1, i)) - 1. /
                   12. * (N.roll(disp[i], -2, i) - N.roll(disp[i], 2, i))) / cellsize
    # divergence 
    return disp[0] + disp[1] + disp[2]


def getkgrid(boxsize, ng, what=None):
    '''
    It returns a meshgrid of kx,ky ,kz and of modulus k
    '''
    kmin = 2 * N.pi / N.float(boxsize)
    thirdim = ng // 2 + 1
    sh = (ng, ng, thirdim)
    kx, ky, kz = N.mgrid[0:sh[0], 0:sh[1], 0:sh[2]].astype(N.float32)

    kx[N.where(kx > ng // 2)] -= ng
    ky[N.where(ky > ng // 2)] -= ng
    kz[N.where(kz > ng // 2)] -= ng

    kx *= kmin
    ky *= kmin
    kz *= kmin

    k = N.sqrt(kx**2 + ky**2 + kz**2)

    if what == None:
        return k
    elif what == 'all':
        return kx, ky, kz, k
    elif what == 'cart':
        return kx, ky, kz


def fastXk(dk1=None, dk2=None, d1=None, d2=None, pos1=None,
           pos2=None, boxsize=None, nkbins=None, kb=0):
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
        sh = (ng1, ng1, ng1 // 2 + 1)

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
        sh = (ng2, ng2, ng2 // 2 + 1)

    ############
    try:
        len(N.shape(dk1)) == len(N.shape(dk2))
        ng = N.shape(dk1)[0]
        sh = len(N.shape(dk1))
    except:
        raise ValueError('the input field must have the same size')

    kgrid = getkgrid(boxsize, ng)
    norm = (boxsize / float(ng)**2.)**sh

    if nkbins == None:
        nkbins = ng // 2

    kgrid = kgrid.flatten()
    Xk = (N.conjugate(dk1) * dk2 + N.conjugate(dk2) * dk1).real.flatten()
    Xk *= 0.5
    Xk = Xk[kgrid > 0.]
    knz = kgrid[kgrid > 0.]

    kbin = N.linspace(
    2 * N.pi / boxsize,
    ng * N.pi / boxsize,
     nkbins + 1)
    xs = N.histogram(N.log10(knz), N.log10(kbin), weights=Xk, range=(
        2 * N.pi / boxsize, ng * N.pi / boxsize))[0]
    counts = N.histogram(N.log10(knz), N.log10(kbin), range=(
        2 * N.pi / boxsize, ng * N.pi / boxsize ))[0]
    binvals = kbin[0:nkbins] + N.diff(kbin) / 2.
    binvals = binvals[counts > 0]
    xs = xs[counts > 0]
    counts = counts[counts > 0]
    xs = xs / counts

    pk1 = fastPk(dk=dk1, boxsize=boxsize)

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
        del pos
        gc.collect()

    elif N.any(d != None):
        sh = N.shape(d)
        ng = sh[0]
        dk = N.fft.rfftn(d)
    else:
        ng = N.shape(dk)[0]
        sh = (ng, ng, ng // 2 + 1)

    sh = len(N.shape(dk))  # dimensions
    kx, ky, kz, kgrid = getkgrid(boxsize, ng, what='all')
    norm = (boxsize / float(ng)**2.)**sh

    # basic deconvolution
    kny = N.pi * ng / boxsize
    dk = dk / (N.sinc(kx/2./kny)*N.sinc(ky/2./kny)*N.sinc(kz/2./kny))**2.
    del kx, ky, kz
    gc.collect()

    kgrid = kgrid.flatten()
    dk2 = abs(dk.flatten()) ** 2.
    dk2 = dk2[kgrid > 0.]
    knz = kgrid[kgrid > 0.]

    if nkbins == None:
        nkbins = ng // 2

    kbin = N.linspace(
    2 * N.pi / boxsize,
    ng * N.pi / boxsize,
     nkbins + 1)
    ps = N.histogram(N.log10(knz), N.log10(kbin), weights=dk2,
                     range=(2 * N.pi / boxsize, ng * N.pi / boxsize ))[0]
    counts = N.histogram(N.log10(knz), N.log10(kbin), range=(
        2 * N.pi / boxsize, ng * N.pi / boxsize))[0]
    binvals = kbin[0:nkbins] + N.diff(kbin) / 2.
    binvals = binvals[counts > 0]
    ps = ps[counts > 0]
    counts = counts[counts > 0]
    ps = ps / counts

    if kb == 0:
        return ps * norm
    else:
        return binvals, ps * norm

def MAS(pos, boxsize, ng):
    """ Mass assignment scheme. It is calling an external C module, _mas,
        which simply assigns mass to density on grid.
    Inputs:
        pos: array (ng*ng*ng,3)/(3,ng,ng,ng) of the positions, output of get_pos
        boxsize: side of the box of the sim
        ng: number of grid resolution
    """
    pos = N.asarray(pos)
    if pos.ndim==4:
        np = N.shape(pos)[-1]
        pos = pos.flatten().astype(N.float32)
    else:
        np = int(N.cbrt(N.shape(pos)[0]))
        pos = N.concatenate((pos[:,0], pos[:,1], pos[:,2]), axis=0).astype(N.float32)
    sh = (ng, ng, ng,)
    dens = N.zeros(sh, dtype=N.float32).flatten()
    pos[pos == boxsize] = 1.0e-06
    mas(ng, np, boxsize, pos, dens)
    dens = dens.reshape(sh)
    return dens

def MASHalos( Nparticles, posHC, boxsize, ng, minpart=20 ):
    """ It computes halo overdensity with Cloud in Cell MAS 
	Inputs::
          Npart: array or list containing number of particles per halo, or mass of each halo
          posHC: (len(Npart),3) array or list of positions of center of mass of haloes in the catalogue
          boxsize: of simulation
          ng: meshgrid size
          minpart: min number of particles
        Returns::
          halo over-density field of shape (ng, ng, ng)
    """

    Nparticles = N.float32(Nparticles)
    print(N.shape(Nparticles))
    posHC = posHC[Nparticles >= minpart, :]
    Nparticles = Nparticles[Nparticles > minpart]

    Nhalos = len(Nparticles)
    print('total number of haloes', Nhalos)
    posHCx, posHCy, posHCz = posHC[:, 0], posHC[:, 1], posHC[:, 2]
    posHC = N.concatenate((posHCx, posHCy, posHCz), axis=0)

    sh = (ng, ng, ng,)
    densH = N.zeros(sh, dtype=N.float32).flatten()
    posHC = N.asarray(posHC).astype(dtype=N.float32)
    posHC = posHC%boxsize
    masHalos(ng, Nhalos, boxsize, posHC, densH, Nparticles)
    densH /= densH.mean()
    densH -= 1.0
    densH = densH.reshape(sh)

    return densH
