"""
    Bunch of functions needed to write the outputs and initial conditions in gadget format.
"""
import numpy as N
import os
import shutil

def writedir(sigmaalpt,
             extra_info,
             scheme='2lpt',
             smallscheme=None,
             redshift=50.,
             boxsize=256.,
             ngrid=256,
             hubble=.7, Omega0=0.3, makeic=False, paramfile=None):
    '''
    it returns the path and it creates the output directory where to save the binaries
    '''

    # name of the output directory
    sims = os.getcwd() + '/sims/'
    if not os.path.exists(sims):
        os.mkdir(sims)
        print('creating the folder ', sims)
    uscr = '_'

    # folder with the specifics of the simulation
    folder = 'bx' + str(boxsize) + uscr + 'ng' + str(ngrid) + uscr + \
        'z' + str(redshift) + uscr + 'Om' + '{0:.2f}'.format(Omega0) + '/'

    folder_sims = sims + folder

    if smallscheme is not None:
        folder_scheme = folder_sims + 'alpt' + '/'

    else:
        folder_scheme = folder_sims + scheme + '/'

    # output folder scheme, in alpt save sigmaalpt
    if (smallscheme == 'sc'):
        fileroot_out = extra_info + 'z' + \
            str(redshift) + 'sigmaalpt' + '{0:.1f}'.format(sigmaalpt) + '__0'

    else:
        fileroot_out = extra_info + 'z' + str(redshift) + '__0'

    if not os.path.exists(folder_sims):
        os.mkdir(folder_sims)
        print('creating the folder ', folder_sims)

    if not os.path.exists(folder_scheme):
        os.mkdir(folder_scheme)
        print('creating the folder ', folder_scheme)

    # case of Initial Conditions
    if makeic == True:
        fileroot_ic = 'IC_' + scheme + '_z' + str(redshift)
        if os.path.exists(folder_sims + fileroot_ic + '.dat'):
            print('these initial conditions already existed!')
            print('not overwriting them')
            return [], []
        else:
            print('initial conditions at ', folder_sims + fileroot_ic + '.dat')
            return folder_sims, fileroot_ic

    # normal run, do not overwrite binaries
    else:
        p = 0
        while os.path.exists(folder_scheme + fileroot_out + '.dat'):
            p += 1
            sep = '__'
            fileroot_out = fileroot_out.split(sep, 1)[0]
            fileroot_out = fileroot_out + sep + str(p)
        if paramfile is not None:
            shutil.copyfile(paramfile, folder_scheme + 'params_' + fileroot_out + '.txt')
        return folder_scheme, fileroot_out


def writeparam(
        path_sims,
        fileroot,
        scheme='2lpt',
        redshift=69,
        boxsize=256.,
        ngrid=256,
        hubble=.7,
        ombh2=0.045 *
        0.73**2,
        Omega0=0.3,
        Softening=None,
        SofteningPhys=None):
    ''' It writes the .param file to give to gadget '''

    # I fix here the redshift of output, in scale factor form
    redshifts = N.array([1 / 3., 1 / 2.2, 1 / 1.7, 1 / 1.5, 1 / 1.3]).T

    if Softening is None:
        Softening = boxsize / ngrid / 35.
    if SofteningPhys is None:
        SofteningPhys = boxsize / ngrid / 35.

    # do not really need these, I use the redshifts
    TimeBetSnapshot = N.sqrt(2.)
    TimeMax = 1.
    TimeOfFirstSnapshot = 1. / 16.

    F = open(path_sims + fileroot + '.param', 'w')
    F.write('InitCondFile  	' + path_sims + fileroot + '.dat\n')

    F.write('OutputDir          ' + path_sims + fileroot + '.OUT/\n')
    if not os.path.exists(path_sims + fileroot + '.OUT/'):
        os.mkdir(path_sims + fileroot + '.OUT')

    F.write('EnergyFile         energy.txt\n')
    F.write('InfoFile           info.txt\n')
    F.write('TimingsFile        timings.txt\n')
    F.write('CpuFile            cpu.txt\n')

    F.write('RestartFile        restart\n')
    F.write('SnapshotFileBase   ' + fileroot + '\n')

    # write the output redshifts
    N.savetxt(
        path_sims +
        'zsnaps.txt',
        redshifts,
        fmt='%.1f',
        delimiter=' ',
        newline='\n',
        header='',
        footer='',
        comments='#',
        encoding=None)
    N.savetxt(
        path_sims +
        fileroot +
        '.OUT/' +
        'zsnaps.txt',
        redshifts,
        fmt='%.1f',
        delimiter=' ',
        newline='\n',
        header='',
        footer='',
        comments='#',
        encoding=None)

    F.write(
        'OutputListFilename         ' +
        path_sims +
        fileroot +
        '.OUT/zsnaps.txt\n')

    F.write('% CPU time -limit\n')

    F.write('TimeLimitCPU      345600  % = 96 hours\n')
    F.write('ResubmitOn        0\n')
    F.write('ResubmitCommand   my-scriptfile  \n')
    F.write('\n')
    F.write('\n')
    F.write('% Code options\n')
    F.write('\n')
    F.write('\n')
    F.write('ICFormat                 1\n')
    F.write('SnapFormat               1\n')
    F.write('ComovingIntegrationOn    1\n')
    F.write('\n')
    F.write('TypeOfTimestepCriterion  0\n')

    F.write('OutputListOn             0\n')
    F.write('PeriodicBoundariesOn     1\n')
    F.write('\n')
    F.write('%  Characteristics of run\n')
    F.write('\n')
    F.write('TimeBegin           %g\n' % (1. / (1. + redshift)))
    F.write('TimeMax             %g\n' % TimeMax)
    F.write('\n')
    F.write('Omega0                %g\n' % Omega0)
    F.write('OmegaLambda           %g\n' % (1. - Omega0))
    F.write('OmegaBaryon           %g\n' % (ombh2 / hubble**2))
    F.write('HubbleParam           %g\n' % hubble)
    F.write('BoxSize               %g\n' % boxsize)
    F.write('\n')
    F.write('% Output frequency\n')
    F.write('TimeBetSnapshot       ' + str(TimeBetSnapshot) + '\n')
    F.write('TimeOfFirstSnapshot   ' + str(TimeOfFirstSnapshot) + '\n')
    F.write('\n')
    F.write('CpuTimeBetRestartFile     36000.0    ; here in seconds\n')
    F.write('TimeBetStatistics         0.05\n')
    F.write('\n')
    F.write('NumFilesPerSnapshot       1\n')
    F.write('NumFilesWrittenInParallel 1\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('% Accuracy of time integration\n')
    F.write('\n')
    F.write('ErrTolIntAccuracy      0.025 \n')
    F.write('\n')
    F.write('MaxRMSDisplacementFac  0.2\n')
    F.write('\n')
    F.write('CourantFac             0.15     \n')
    F.write('\n')
    F.write('MaxSizeTimestep       0.03\n')
    F.write('MinSizeTimestep       0.0\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('\n')
    F.write('% Tree algorithm, force accuracy, domain update frequency\n')
    F.write('\n')
    F.write('ErrTolTheta            0.5            \n')
    F.write('TypeOfOpeningCriterion 1\n')
    F.write('ErrTolForceAcc         0.005\n')
    F.write('\n')
    F.write('\n')
    F.write('TreeDomainUpdateFrequency    0.1\n')
    F.write('\n')
    F.write('\n')
    F.write('%  Further parameters of SPH\n')
    F.write('\n')
    F.write('DesNumNgb              60\n')
    F.write('MaxNumNgbDeviation     2\n')
    F.write('ArtBulkViscConst       0.8\n')
    F.write('InitGasTemp            0        % always ignored if set to 0 \n')
    F.write('MinGasTemp             50.0    \n')
    F.write('\n')
    F.write('\n')
    F.write('% Memory allocation\n')
    F.write('\n')
    F.write('PartAllocFactor       1.25\n')
    F.write('TreeAllocFactor       0.8\n')
    F.write('BufferSize            100          % in MByte\n')
    F.write('\n')
    F.write('\n')
    F.write('% System of units\n')
    F.write('\n')
    F.write('UnitLength_in_cm         3.085678e24        ;  1.0 Mpc/h \n')
    F.write('UnitMass_in_g            1.989e43           ;  1.0e10 solar masses/h \n')
    F.write('UnitVelocity_in_cm_per_s 1e5                ;  1 km/sec \n')
    F.write('GravityConstantInternal  0\n')
    F.write(' \n')
    F.write('\n')
    F.write('% Softening lengths\n')
    F.write('\n')
    F.write('MinGasHsmlFractional 0.25\n')
    F.write('\n')
    F.write('SofteningGas       0\n')
    F.write('SofteningHalo      ' + str(Softening) + '\n')
    F.write('SofteningDisk      0\n')
    F.write('SofteningBulge     0           \n')
    F.write('SofteningStars     0\n')
    F.write('SofteningBndry     0\n')
    F.write('\n')
    F.write('SofteningGasMaxPhys       0.0\n')
    F.write('SofteningHaloMaxPhys      ' + str(SofteningPhys) + '\n')
    F.write('SofteningDiskMaxPhys      0\n')
    F.write('SofteningBulgeMaxPhys     0           \n')
    F.write('SofteningStarsMaxPhys     0\n')
    F.write('SofteningBndryMaxPhys     0\n')
    F.close()

    return 1


def writegadget(
        pos,
        vel,
        redshift,
        boxsize,
        OmegaM,
        OmegaL,
        HubbleParam,
        path_sims,
        fileroot,
        id=None):

    pos = N.array(pos)
    vel = N.array(vel)

    pos = N.asfortranarray(
        [pos[0].flatten(), pos[1].flatten(), pos[2].flatten()])
    vel = N.asfortranarray(
        [vel[0].flatten(), vel[1].flatten(), vel[2].flatten()])

    ng = pos.shape[1]

    npartarr = N.array([0, ng, 0, 0, 0, 0]).astype(N.uint32)
    # in 10^10 Msolar/h
    mass = N.array([0., 27.7536627 * (boxsize**3 / float(ng))
                    * OmegaM, 0., 0., 0., 0.]).astype(N.float64)
    time = 1. / (1. + redshift)  # scale factor
    F = open(path_sims + fileroot + '.dat', 'wb')
    F.write(N.array(256).astype(N.int32))  # f77dummy
    F.write(npartarr)  # 6-member array
    F.write(mass)  # 6-member array
    F.write(N.array([time]).astype(N.float64))
    F.write(N.array([redshift]).astype(N.float64))
    F.write(N.array([0, 0]).astype(N.int32))  # flagsfr, flag_feedback
    F.write(N.array(npartarr).astype(N.int32))
    F.write(N.array([0, 1]).astype(N.int32))  # FlagCooling,NumFiles
    F.write(N.array([boxsize, OmegaM, OmegaL, HubbleParam]).astype(N.float64))
    # F.write(N.array([0,0]).astype(N.int32)) #FlagCooling,NumFiles
    bytesleft = 256 - 6 * 4 - 6 * 8 - 8 - 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8  # -2*4
    la = N.zeros(bytesleft // 4, dtype=N.int32)
    F.write(la)
    F.write(N.array(256).astype(N.int32))  # f77dummy

    #pos, vel, id
    F.write(N.array([12 * ng]).astype(N.int32))  # f77dummy
    F.write((pos.astype(N.float32)).T)
    F.write(N.array([12 * ng]).astype(N.int32))  # f77dummy
    F.write(N.array([12 * ng]).astype(N.int32))  # f77dummy
    F.write(((vel / N.sqrt(time)).astype(N.float32)).T)
    # note 1/sqrt(a)
    F.write(N.array([12 * ng]).astype(N.int32))  # f77dummy

    if id is None:
        id = 1 + N.arange(ng)  # Note infuriating "1"!
        id = N.array(id).astype(N.uint32)

    F.write(N.array([4 * ng]).astype(N.int32))  # f77dummy
    F.write(id)  # id array
    F.write(N.array([4 * ng]).astype(N.int32))  # f77dummy

    F.close()

    return
