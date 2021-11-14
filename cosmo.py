'''
    Cosmological module to compute linear theory related quantities.
    The transfer functions of Eisenstein and Hu power spectrum were plagiarized from COLOSSUS.
    There is the possibility to use Class
'''

import numpy as N
import math as M
import os
from scipy.integrate import romberg, quad
from scipy.optimize import newton
import scipy.interpolate
from colossus.cosmology import cosmology
from colossus.lss import mass_function
try:
    from classy import Class
except ImportError:
    pass


class PowerSpectrum(object):
    ''' Parent class to compute the linear power spectrum. You just need to specify
        the cosmological parameters'''

    def __init__(self, h=0.7, omega_b=0.0225, Omega_cdm=0.25, ns=0.96,
                 sigma8=0.8, num_sigma_bins=100, log_mass_lim=(11, 15), redshift=0.):
        """
        Initialize cosmological parameters
        """
        self.h = h
        self.omega_b = omega_b
        self.omega_0 = (omega_b + Omega_cdm * h**2)
        self.Omega_0 = omega_b / h / h + Omega_cdm
        self.ns = ns
        self.sigma8 = sigma8

        # need to change this if the box of sim
        # is huge or vert high res
        self.log10kmin = -3.0
        self.log10kmax = N.log10(10.)
        self.ksteps = 1000

        # halo model part
        self.redshift = redshift
        self.num_sigma_bins = num_sigma_bins
        self.log_mass_max = log_mass_lim[1]
        self.log_mass_min = log_mass_lim[0]
        self.rho_c0 = 2.77536627e+11  # in units of h^2 Msolar/Mpc^3
        rho0 = self.Omega_0 * self.rho_c0  # in units of h^2 Msolar/Mpc^3
        # units of Msolar/h
        self.Marray = N.logspace(
            self.log_mass_min, self.log_mass_max, self.num_sigma_bins)
        # Comoving radius is in units of h^-1 Mpc
        self.Rarray = N.power(3. * self.Marray / 4. / N.pi / rho0, 1. / 3.)

    def E(self, z):
        """
        Omega as a function of redshift
        """
        return (self.Omega_0 * (1. + z)**3 + 1. - self.Omega_0)**0.5

    def Om0z(self, z):
        """
        Omega as a function of redshift
        """
        return self.Omega_0 * (1. + z)**3 / self.E(z)**2.

    def OmLz(self, z):
        """
        Lambda as a function of redshift
        """
        return (1 - self.Omega_0) / self.E(z)**2.

    def d1(self, z):
        """
        Growth function, as in Bernardeau01. At early times d1~a
        """
        return 5 * self.Om0z(z) / (1.0 + z) / 2.0 / (self.Om0z(z)**(4. / 7.) -
                                                     self.OmLz(z) + (1.0 + self.Om0z(z) / 2.) * (1. + self.OmLz(z) / 70.))

    def Delta_v(self, z):
        """
        Spherical overdensity ratio at virialization
        """
        x = self.Om0z(z) - 1.
        return (18 * N.pi**2. + 82. * x - 39. * x**2.) / self.Om0z(z)

    def SC_delta_c(self, z):
        return 1.686

    def tablehmf(self, z, boxsize, ng, Dm=0.05):
        ''' write target hmf in a binary '''

        params = {'flat': True, 'H0': self.h * 100., 'Om0': self.Omega_0,
                  'Ob0': self.omega_b / self.h**2., 'sigma8': self.sigma8, 'ns': self.ns}

        colossus = cosmology.setCosmology('myCosmo', params)

        filen = 'hmf_at_z' + str(z) + '_' + 'Om' + \
            '{0:.2f}'.format(self.Omega_0) + '.dat'
        path = os.getcwd()
        filename = path + '/' + filen
        print(filename)

        # min halo (depends on res)
        logMmin = N.log10(self.rho_c0 * self.Omega_0 * (boxsize / ng)**3.)
        print('logMmin=', logMmin)

        # max halo (depends on boxsize)
        logM = N.linspace(10., 17., 300).astype(N.float64)
        DlogM = logM[1] - logM[0]
        DlnM = N.log(10**DlogM)
        t08 = mass_function.massFunction(
            10**logM, z, mdef='340m', model='tinker08', q_out='dndlnM').astype(N.float64)
        t08 *= boxsize**3 * DlnM
        log10t08 = N.log10(t08)
        intp = scipy.interpolate.InterpolatedUnivariateSpline(
            logM, log10t08, k=1)

        def tmp(logmass): return 10**intp(logmass) - 1.0
        logMmax = newton(func=tmp, x0=11., tol=1e-03, maxiter=100)
        print('logMmax=', logMmax)

        # write lookup table
        Nm = int(round((logMmax - logMmin) / Dm))
        logM = N.linspace(logMmin, logMmax, Nm).astype(N.float64)
        t08 = mass_function.massFunction(
            10**logM, z, mdef='340m', model='tinker08', q_out='dndlnM').astype(N.float64)

        fileobj = open(filename, mode='wb')
        Nm = N.asarray([Nm], dtype=N.int32).tofile(fileobj)
        logM.tofile(fileobj)
        t08.tofile(fileobj)
        fileobj.close

        return filen


class EisHu(PowerSpectrum):
    '''
    Numerical Calculation of the Power Spectrum, using Eisenstein&Hu 1998b fitting formulae
    (26) and (28)-(31). It includes growth factor given by Lahav 1991 fit, to rescale it.
    Output k is in units of h/Mpc, but be careful, in EHU paper k is in 1/Mpc
    '''

    def __init__(self, h=0.7, omega_b=0.0225,
                 Omega_cdm=0.25, ns=0.96, sigma8=0.8):
        """
        Initialize cosmological parameters
        """
        PowerSpectrum.__init__(
            self,
            h=0.7,
            omega_b=0.0225,
            Omega_cdm=0.25,
            ns=0.96,
            sigma8=0.8)
        self.delH = 1.94e-5 * self.Omega_0**(-0.785 - 0.05 * N.log(
            self.Omega_0)) * N.exp(-0.95 * (self.ns - 1) - 0.169 * (self.ns - 1)**2)

        self.sigma_squared_of_8_raw = self.sigma_squared_of_R(R=8.)
        self.Normalization = self.sigma8**2. / self.sigma_squared_of_8_raw
        self.sigmaarray = [
            (N.sqrt(self.Normalization * self.sigma_squared_of_R(R))) for R in self.Rarray]

    def T0(self, k):
        """
        transfer function, no wiggle case
        (COLOSSUS)
        """
        # Eq. 31
        alfa_T = 1 - 0.328 * N.log(431. * self.omega_0) * self.omega_b / self.omega_0 \
            + 0.38 * N.log(22.3 * self.omega_0) * \
            (self.omega_b / self.omega_0)**2.

        # Eq. 26
        s = 44.5 * N.log(9.83 / self.omega_0) / \
            N.sqrt(1 + 10 * self.omega_b**0.75)

        # Eq. 30
        GammaEffk = self.Omega_0 * self.h * \
            (alfa_T + (1 - alfa_T) / (1 + (0.43 * k * s * self.h)**4.))

        # Eq. 28
        q = k * (2.725 / 2.7)**2 / GammaEffk

        # Eq. 29
        L0 = N.log(2 * M.e + 1.8 * q)
        C0 = 14.2 + 731. / (1 + 62.5 * q)

        return L0 / (L0 + C0 * q**2)

    def T(self, k):
        '''
        Transfer function of baryons and cold dark matter
        (COLOSSUS)
        '''
        # Define shorter expressions
        Tcmb0 = 2.725
        Om0 = self.Omega_0
        om0h2 = self.omega_0
        omc = self.Omega_0 - self.omega_b / self.h**2.
        ombom0 = self.omega_b / self.omega_0
        h = self.h
        h2 = self.h**2
        om0h2 = self.Omega_0 * h2
        ombh2 = self.omega_b
        theta2p7 = Tcmb0 / 2.7
        theta2p72 = theta2p7**2
        theta2p74 = theta2p72**2

        # Convert kh from h/Mpc to 1/Mpc
        kh = k * h

        # Equation 2
        zeq = 2.50e4 * om0h2 / theta2p74

        # Equation 3
        keq = 7.46e-2 * om0h2 / theta2p72

        # Equation 4
        b1d = 0.313 * om0h2**-0.419 * (1.0 + 0.607 * om0h2**0.674)
        b2d = 0.238 * om0h2**0.223
        zd = 1291.0 * om0h2**0.251 / \
            (1.0 + 0.659 * om0h2**0.828) * (1.0 + b1d * ombh2**b2d)

        # Equation 5
        Rd = 31.5 * ombh2 / theta2p74 / (zd / 1e3)
        Req = 31.5 * ombh2 / theta2p74 / (zeq / 1e3)

        # Equation 6
        s = 2.0 / 3.0 / keq * N.sqrt(6.0 / Req) * N.log((N.sqrt(1.0 + Rd) +
                                                         N.sqrt(Rd + Req)) / (1.0 + N.sqrt(Req)))

        # Equation 7
        ksilk = 1.6 * ombh2**0.52 * om0h2**0.73 * (1.0 + (10.4 * om0h2)**-0.95)

        # Equation 10
        q = kh / 13.41 / keq

        # Equation 11
        a1 = (46.9 * om0h2)**0.670 * (1.0 + (32.1 * om0h2)**-0.532)
        a2 = (12.0 * om0h2)**0.424 * (1.0 + (45.0 * om0h2)**-0.582)
        ac = a1**(-ombom0) * a2**(-ombom0**3)

        # Equation 12
        b1 = 0.944 / (1.0 + (458.0 * om0h2)**-0.708)
        b2 = (0.395 * om0h2)**-0.0266
        bc = 1.0 / (1.0 + b1 * ((omc / Om0)**b2 - 1.0))

        # Equation 15
        y = (1.0 + zeq) / (1.0 + zd)
        Gy = y * (-6.0 * N.sqrt(1.0 + y) + (2.0 + 3.0 * y)
                  * N.log((N.sqrt(1.0 + y) + 1.0) / (N.sqrt(1.0 + y) - 1.0)))

        # Equation 14
        ab = 2.07 * keq * s * (1.0 + Rd)**(-3.0 / 4.0) * Gy

        # Get CDM part of transfer function

        # Equation 18
        f = 1.0 / (1.0 + (kh * s / 5.4)**4)

        # Equation 20
        C = 14.2 / ac + 386.0 / (1.0 + 69.9 * q**1.08)

       # Equation 19
        T0t = N.log(N.e + 1.8 * bc * q) / \
            (N.log(N.e + 1.8 * bc * q) + C * q * q)

        # Equation 17
        C1bc = 14.2 + 386.0 / (1.0 + 69.9 * q**1.08)
        T0t1bc = N.log(N.e + 1.8 * bc * q) / \
            (N.log(N.e + 1.8 * bc * q) + C1bc * q * q)
        Tc = f * T0t1bc + (1.0 - f) * T0t

        # Get baryon part of transfer function

        # Equation 24
        bb = 0.5 + ombom0 + (3.0 - 2.0 * ombom0) * \
            N.sqrt((17.2 * om0h2) * (17.2 * om0h2) + 1.0)

        # Equation 23
        bnode = 8.41 * om0h2**0.435

        # Equation 22
        st = s / (1.0 + (bnode / kh / s) * (bnode / kh / s)
                  * (bnode / kh / s))**(1.0 / 3.0)

        # Equation 21
        C11 = 14.2 + 386.0 / (1.0 + 69.9 * q**1.08)
        T0t11 = N.log(N.e + 1.8 * q) / (N.log(N.e + 1.8 * q) + C11 * q * q)
        Tb = (T0t11 / (1.0 + (kh * s / 5.2)**2) + ab / (1.0 + (bb / kh / s)**3) * N.exp(-(kh / ksilk)**1.4)) \
            * N.sin(kh * st) / (kh * st)

        # Total transfer function
        Tk = ombom0 * Tb + omc / Om0 * Tc

        return Tk

    def Pk_prim(self, k):
        """
        unnormalized primordial power spectrum
        """
        return k**self.ns

    def pk_ehu(self, k, z):
        """
        linear power spectrum.
        """
        return self.delH**2. * \
            self.Pk_prim(k) * self.T(k)**2 * (self.d1(z) / self.d1(0))**2.

    def sigma_squared_of_R(self, R=8.):
        """
           calculates sigma^2(R). Integrates the sigma_squared_integrand
           parameter from 0 to infinity. R is in Mpc/h
        """
        result = quad(self.sigma_squared_integrand, 0., N.Inf,
                      args=(R,), epsrel=1e-16, epsabs=1e-16)[0]
        result = result / (2.0 * M.pi**2)
        return result

    def sigma_squared_integrand(self, k, R):
        """
        integrand for integral to get sigma^2(R).
        Parameters: R is Mpc/h units
        """
        f = k**2 * self.pk_ehu(k, 0.) * self.WofK(R * k)**2.0
        return f

    def WofK(self, x):
        """
        returns W(k*R), which is the fourier transform of the top-hat function.
        """
        return 3.0 * (N.sin(x) - x * N.cos(x)) / x**3

    def pk_lin(self, k, z):
        """
        CDM linear power spectrum with sigma8 accounted for.
        It is in units of Mpc^3/h^3
        """
        return self.sigma8**2. / \
            self.sigma_squared_of_8_raw * self.pk_ehu(k, z)

    def table(self, z):
        """
        write the power spectrum and k on a table
        """
        import os
        file = 'Pk_at_z' + \
            '{0:.1f}'.format(z) + '_' + 'Om' + \
            '{0:.2f}'.format(self.Omega_0) + '.txt'
        k = N.logspace(self.log10kmin, self.log10kmax, self.ksteps)
        pk = self.pk_lin(k, z)
        kPk = N.vstack((k, pk)).T
        path = os.getcwd()
        N.savetxt(path + '/' + file, kPk, fmt=['%.4f', '%.4f'],
                  header='k   Pk', comments='# ')
        return 0

    def sigmaMz(self, mass, z):
        """
        define the interpolant for sigma(M,z). M is in Msolar/h units
        """
        self.intp = scipy.interpolate.InterpolatedUnivariateSpline(
            self.Marray, self.sigmaarray, k=1)
        return self.d1(z) / self.d1(0) * self.intp(mass)

    def Mzcrit(self, z):
        """
        function whose root is Mcrit in Msolar/h units
        """
        def tmp(M, z): return self.sigmaMz(M, z) - self.SC_delta_c(z)
        return newton(func=tmp, x0=1e+10, args=(z,), tol=1e+7, maxiter=100)


class PSClass(PowerSpectrum):
    '''Class to initialize CLASS instance.
       In addition to cosmological paramters, you need to specify the redshift at which Class
       computes the power spectrum (z_pk), and the max k (Pk_max_h0Mpc)
    '''

    def __init__(self, h=0.7, omega_b=0.0225, Omega_cdm=0.25,
                 ns=0.96, sigma8=0.8, z_pk=50.):
        PowerSpectrum.__init__(
            self,
            h=h,
            omega_b=omega_b,
            Omega_cdm=Omega_cdm,
            ns=ns,
            sigma8=sigma8)
        Pk_max_h0Mpc = 10**self.log10kmax
        self.class_params = {
            'output': 'mPk',
            # 'non linear': '',
            'P_k_max_h/Mpc': Pk_max_h0Mpc,
            'z_pk': z_pk,
            'A_s': 2.3e-9,
            'n_s': ns,
            'h': h,
            'omega_b': omega_b,
            'Omega_cdm': Omega_cdm,
        }

        # instance of CLASS
        self.C = self.Compute()
        Rphys = self.Rarray / h  # CLASS wants physical length
        self.sigmaarray = [   # already normalized
            (N.sqrt(self.C.sigma(r, self.redshift))) for r in Rphys]

    def Compute(self):
        '''initialize the cosmological module for the parameters'''
        c = Class()
        c.set(self.class_params)
        c.compute()
        sig8 = c.sigma8()
        A_s = c.pars['A_s']
        c.struct_cleanup()  # does not clean the input class_params, cosmo.empty() does that
        c.set(A_s=A_s * (self.sigma8 * 1. / sig8)**2)
        # recompute P_k
        c.compute()
        return c

    def pk_lin(self, k, z):
        """
        Class linear power spectrum, to account for Mpc^3/h^3 units,
        and to deal with arrays
        """
        pk = N.array([self.C.pk_lin(x, z) for x in k * self.h]) * self.h**3.

        return pk

    def table(self, z):
        """
        write the power spectrum and k on a table
        """
        import os
        file = 'Pk_at_z' + \
            '{0:.1f}'.format(z) + '_' + 'Om' + \
            '{0:.2f}'.format(self.Omega_0) + '.txt'
        k = N.logspace(self.log10kmin, self.log10kmax, self.ksteps)
        pk = self.pk_lin(k, z)
        kPk = N.vstack((k, pk)).T
        path = os.getcwd()
        N.savetxt(path + '/' + file, kPk, fmt=['%.4f', '%.4f'],
                  header='k   Pk', comments='# ')
        return 0

    def sigmaMz(self, mass, z):
        """
        define the interpolant for sigma(M,z). M is in Msolar/h units
        """
        self.intp = scipy.interpolate.InterpolatedUnivariateSpline(
            self.Marray, self.sigmaarray, k=1)
        return self.d1(z) / self.d1(0) * self.intp(mass)

    def Mzcrit(self, z):
        """
        function whose root is Mcrit in Msolar/h units
        """
        def tmp(M, z): return self.sigmaMz(M, z) - self.SC_delta_c(z)
        return newton(func=tmp, x0=1e+10, args=(z,), tol=1e+7, maxiter=100)
