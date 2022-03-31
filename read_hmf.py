import readgadget
import numpy as N
import matplotlib.pyplot as plt
import struct

def read_hmf(pathtobin,bins=20):
    ''' This function reads the binary file containing the particle information
        about haloes, and builds the halo mass function '''

    head = readgadget.header(pathtobin)
    boxsize = head.boxsize
    omega_m = head.omega_m
    npart   = head.npart[1]
    ng = int(N.cbrt(npart))
    mass_res = 2.77536627e+11 * omega_m * (boxsize/float(ng))**3.

    with open( pathtobin[:-4]+'_hc.dat', mode='rb') as file:
        fileContent = file.read()
        Nhalos = struct.unpack('i',fileContent[:4])[0]
        sizeofall=4 #remove Nhalos
        size=N.zeros(Nhalos,dtype=N.int32)
        j=0
        while j<Nhalos:
            size[j] = struct.unpack('i',fileContent[sizeofall:sizeofall+4])[0]
            sizeofall = size[j]+sizeofall
            j+=1

    csize = N.append(0,size.cumsum())
    hc = [struct.unpack('iiiidddddd'+'i'*((size[j]-4*4-6*8)//4),fileContent[4+csize[j]:4+csize[j+1]]) for j in range(Nhalos)]
    np =  [hc[i][2] for i in range(Nhalos)]
    masses = N.asarray([ hc[i][2] for i in range(Nhalos) ])
    masses = masses.astype(N.float32)
    masses *= mass_res
    logMmin = N.log10(20*mass_res)
    logM = N.linspace(logMmin, 15., bins+1).astype(N.float64)
    nM, bindg = N.histogram(N.log10(masses),bins=logM)            # compute histogram of halo masses
    dlogM = bindg[1]-bindg[0]
    nM = nM/boxsize**3
    nM = nM/N.log(10**dlogM)
    binmean = (bindg[1:]+bindg[:-1])/2.                               # right edge of mass bin
    binmean = binmean[nM>0]
    nM = nM[nM>0]
    nM = nM[10**binmean/mass_res>=20]
    binmean = binmean[10**binmean/mass_res>=20]
    hmf = N.vstack((binmean,nM))

    return hmf

if __name__ == "__main__":
    hmf = read_hmf('sims/bx200.0_ng300_z0.0_Om0.30/muscleups/z0.0__0.dat')
    plt.semilogy(hmf[0],hmf[1])
    plt.show()
