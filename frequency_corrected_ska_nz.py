#!/usr/bin/python
"""
Calculate dn/dz for SKA HI galaxy redshift surveys, using dn/dz curves for 
given flux thresholds (from Mario Santos) and flux scalings with redshift for 
specific arrays.
"""
import numpy as np
import pylab as P
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import sys
from fit_Euclid import dvdz

DEBUG_PLOT = True #True # Whether to plot fitting functions or not

NU_LINE = 1.420 # HI emxission line freq. in GHz
FULLSKY = (4.*np.pi * (180./np.pi)**2.) # deg^2 in the full sky
NGAL_MIN = 1e3 # Min. no. of galaxies to tolerate in a redshift bin
CBM = 1. #np.sqrt(1.57) # Correction factor due to effective beam for MID/MK (OBSOLETE)
CTH = 0.5 # Correction factor due to taking 5 sigma (not 10 sigma) cuts for SKA1
SBIG = 500. # Flux rms to extrapolate dn/dz out to (constrains behaviour at large Srms)
C = 3e5 # Speed of light, km/s



# Planck-only best-fit parameters, from Table 2 of Planck 2013 XVI.
cosmo = {
    'omega_M_0':        0.316,
    'omega_lambda_0':   0.684,
    'omega_b_0':        0.049,
    'omega_HI_0':       6.50e-4,
    'N_eff':            3.046,
    'h':                0.67,
    'ns':               0.962,
    'sigma_8':          0.834,
    'gamma':            0.55,
    'w0':               -1.,
    'wa':               0.,
    'fNL':              0.,
    'mnu':              0.,
    'k_piv':            0.05, # n_s
    'aperp':            1.,
    'apar':             1.,
    'bHI0':             0.702,
    'A':                1.,
    'sigma_nl':         7.,
    'beta_1':           0.,         # Scale-dependent bias (k^1 term coeff. [Mpc])
    'beta_2':           0.          # Scale-dependent bias (k^2 term coeff. [Mpc^2])
}




def background_evolution_splines(cosmo, zmax=10., nsamples=500):
    """
    Get interpolation functions for background functions of redshift:
      * H(z), Hubble rate in km/s/Mpc
      * r(z), comoving distance in Mpc
      * D(z), linear growth factor
      * f(z), linear growth rate
    """
    _z = np.linspace(0., zmax, nsamples)
    a = 1. / (1. + _z)
    H0 = (100.*cosmo['h']); w0 = cosmo['w0']; wa = cosmo['wa']
    om = cosmo['omega_M_0']; ol = cosmo['omega_lambda_0']
    ok = 1. - om - ol
    
    # Sample Hubble rate H(z) and comoving dist. r(z) at discrete points
    omegaDE = ol * np.exp(3.*wa*(a - 1.)) / a**(3.*(1. + w0 + wa))
    E = np.sqrt( om * a**(-3.) + ok * a**(-2.) + omegaDE )
    _H = H0 * E
    
    r_c = np.concatenate( ([0.], scipy.integrate.cumtrapz(1./E, _z)) )
    if ok > 0.:
        _r = C/(H0*np.sqrt(ok)) * np.sinh(r_c * np.sqrt(ok))
    elif ok < 0.:
        _r = C/(H0*np.sqrt(-ok)) * np.sin(r_c * np.sqrt(-ok))
    else:
        _r = (C/H0) * r_c
    
    # Integrate linear growth rate to find linear growth factor, D(z)
    # N.B. D(z=0) = 1.
    a = 1. / (1. + _z)
    Oma = cosmo['omega_M_0'] * (1.+_z)**3. * (100.*cosmo['h']/_H)**2.
    _f = Oma**cosmo['gamma']
    _D = np.concatenate( ([0.,], scipy.integrate.cumtrapz(_f, np.log(a))) )
    _D = np.exp(_D)
    
    # Construct interpolating functions and return
    r = scipy.interpolate.interp1d(_z, _r, kind='linear', bounds_error=False)
    H = scipy.interpolate.interp1d(_z, _H, kind='linear', bounds_error=False)
    D = scipy.interpolate.interp1d(_z, _D, kind='linear', bounds_error=False)
    f = scipy.interpolate.interp1d(_z, _f, kind='linear', bounds_error=False)
    return H, r, D, f

def extend_with_linear_interp(xnew, x, y):
    """
    Extend an array using a linear interpolation from the last two points.
    """
    dx = x[-1] - x[-2]
    dy = y[-1] - y[-2]
    ynew = y[-1] + dy * (xnew - x[-1]) / dx
    y = np.concatenate((y, [ynew,]))
    return y

def n_bin(zmin, zmax, dndz, bias=None):
    """
    Number density of galaxies in a given z bin (assumes full sky). Also 
    returns volume of bin. dndz argument expects an interpolation fn. in units 
    of deg^-2.
    """
    _z = np.linspace(zmin, zmax, 500)
    vol = 4.*np.pi*C * scipy.integrate.simps(r(_z)**2. / H(_z), _z)
    N_bin = FULLSKY * scipy.integrate.simps(dndz(_z), _z)
    nz = N_bin / vol
    
    # Calculate mean bias (weighted by number density)
    if bias is not None:
        b = scipy.integrate.simps(bias(_z)*dndz(_z), _z) / (N_bin / FULLSKY)
        #b = bias(0.5*(zmin+zmax))
        return nz, vol, b
    return nz, vol

def redshift_bins(dz=0.1, Nbins=None):
    """
    Calculate redshift bins.
    """
    zmin = NU_LINE*1e3 / numax[ID] - 1.
    zmax = NU_LINE*1e3 / numin[ID] - 1.
    print zmin, zmax, 100000000
    if zmin < 0.: zmin = 0.
    if Nbins is not None:
        zbins = np.linspace(zmin, zmax, Nbins+1)
    else:
        Nbins = np.floor((zmax - zmin) / dz)
        zbins = np.linspace(zmin, zmin + dz*Nbins, Nbins+1)
        if zmax - np.max(zbins) > 0.04:
            zbins = np.concatenate((zbins, [zmax,]))
    return zbins

def flux_redshift(z):
    """
    Flux rms as a function of redshift.
    """
    z = np.atleast_1d(z)
    nu = NU_LINE / (1. + z)
    if nucrit[ID] is not None:
        Sz = fluxrms[ID] * Scorr[ID] * np.ones(nu.size)
        idxs = np.where(nu*1e3 > nucrit[ID])
        Sz[idxs] *= (nu[idxs]*1e3 / nucrit[ID])
    else:
        Sz = fluxrms[ID] * Scorr[ID]
        Sz = nu * Sz if not Sconst[ID] else Sz * np.ones(nu.size)
    return Sz

def k_max(z,k_NL):
    kmax = k_NL*(1.+z)**(2./(2.+cosmo['ns']))
    return kmax

if __name__ == '__main__': 



	# Specify which experiment to calculate for
	try:
    		print 'try 0,  3,  5, 9 for SKA1 and  14, 15, 17 for SKA2 or try 10 for SKA-SUR1'
    		ID =int(raw_input('input a number from 0 to 10 to identify the experiment:'))
	except:
    		print "Expects 1 argument: int(experiment_id)"
    		sys.exit()
	#ax1 = P.subplot(211)
	#ax2 = P.subplot(212)
	#for ID in [4, 5, 9, 10]:

	# Flux rms limits at 1GHz for various configurations
	name = ['0n', 'SKA1MID-B1','SKA1MID_B2','SKA1MID_B2_Opt', 'SKA1MID_B2_Real', 'SKA1MID_B2_Pess', 'MEERKAT-B1', 'MEERKAT-B2', 'MID+MK-B1',
        'MID_MK_B2_Real', 'SKA1SUR-B1', 'SKA1SUR-B2', 'ASKAP', 'SUR+ASKAP', 'SKA2_Opt', 'SKA2_Real','SKA2_Pess' ,'SKA2']
	#fluxrms = [251., 149., 555., 598., 197., 122., 174., 192., 645., 139., 5.] # Old
	fluxrms = [0.0, 315., 187., 70,  100., 200.,  696., 750., 247., 152., 174., 192., 645., 179., 3., 7.3, 23.,  5.4]
	numin =   [500, 350., 950., 950., 950., 950.,  580., 900., 580., 950., 350., 650., 700., 700.,500., 500., 500., 500.]
	numax = [1200, 1050., 1760., 1760.,1760., 1760.,1015., 1670., 1015., 1670., 900., 1670., 1800., 1670., 1200., 1200., 1200., 1200.]
	nucrit= [None, None, None, None, None, None, None, None, None, None, 710., 1300., 1250., 1300., None, None, None, None]
	Sarea = [30e3, 5e3, 5e3, 5e3,5e3 ,5e3 , 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 5e3, 30e3, 30e3, 30e3, 30e3]
	Sconst= [True, False, False, False, False, False, False, False, False, False, False, False, False, False,True,True ,True,True]
	Scorr = [1, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CBM*CTH, CTH, CTH, CTH, CTH, 1.,1., 1. ,1.]

	# Define fitting coefficients from Mario's note (HI_specs.pdf)
	Srms = np.array([0., 1., 3., 5., 6., 7.3, 10., 23., 40., 70., 100., 150.,])# 200.,])
        c1 = [6.124, 6.556, 6.532, 6.551, 6.577, 6.555, 6.443, 6.02, 5.74, 5.62, 5.63, 5.48,]# 5.00]
	c2 = [1.717, 2.0185, 1.932, 1.932, 1.946, 1.923, 1.831, 1.43, 1.22, 1.11, 1.41, 1.33,]# 1.04]
	c3 = [0.789, 3.810, 5.225, 6.225,6.685, 7.078, 7.585, 9.03, 10.58, 13.03, 15.49, 16.62,]# 17.52]
	#c4 = [0.8695, 0.5863, 0.4780, 0.5884, 0.5908, 0.5088, 0.4489, 0.5751, 0.5125, 0.6193, 0.6212, 1., 1.]#, 1.]
	#c5 = [0.2338, 0.6410, 0.9181, 0.8076, 0.8455, 1.0222, 1.2069, 0.9993, 1.1842, 1.0179, 1.0759, 0., 0.]#, 0.]
	c4 =  [0.587, 0.497, 0.530, 0.550, 0.547, 0.562, 0.593, 0.607, 0.628, 0.609, 0.605,0.637, 1.0]
	c5 =  [0.358, 0.720, 0.781, 0.801, 0.829, 0.823, 0.807, 0.852, 0.844, 0.929, 1.086,0.965, 0.0]
	c1 = np.array(c1); c2 = np.array(c2); c3 = np.array(c3)
	c4 = np.array(c4); c5 = np.array(c5)
	Smax = np.max(Srms)    
	# Extrapolate fitting functions to high flux rms
	c1 = extend_with_linear_interp(SBIG, Srms, c1)
	c2 = np.concatenate((c2, [1.,])) # Asymptote to linear fn. of redshift
	c3 = extend_with_linear_interp(SBIG, Srms, c3)
	Srms = np.concatenate((Srms, [SBIG,]))

	# Construct grid of dn/dz (deg^-2) as a function of flux rms and redshift and 
	# then construct 2D interpolator
	z = np.linspace(0., 4., 400)
	nu = NU_LINE / (1. + z)
	_dndz = np.array([10.**c1[j] * z**c2[j] * np.exp(-c3[j]*z) for j in range(Srms.size)])
	_bias = np.array([c4[j] * np.exp(c5[j]*z) for j in range(Srms.size)])

	#print _dndz
	dndz = scipy.interpolate.RectBivariateSpline(Srms, z, _dndz, kx=1, ky=1)
	bias = scipy.interpolate.RectBivariateSpline(Srms, z, _bias, kx=1, ky=1)

	#print dndz
	# Construct dndz(z) interpolation fn. for the sensitivity of actual experiment
	fsky = Sarea[ID] / FULLSKY
	Sz = flux_redshift(z)
	dndz_expt = scipy.interpolate.interp1d(z, dndz.ev(Sz, z))
	bias_expt = scipy.interpolate.interp1d(z, bias.ev(Sz, z))
	print 'dndz_expt', dndz_expt

	# Fit function to dn/dz [deg^-2]
	_z = np.linspace(0., 1., 100)
	dndz_vals = dndz_expt(_z)
	bias_vals = bias_expt(_z)
	p0 = [100.*np.max(dndz_vals), 2., 10.]
	def lsq(params):
    		A, c2, c3 = params
    		model = A * _z**c2 * np.exp(-c3*_z)
    		return model - dndz_vals
	p = scipy.optimize.leastsq(lsq, p0)[0]

	# Fit function to bias
	p0 = [np.max(bias_vals), 0.5]
	def lsq(params):
    		c4, c5 = params
    		model = c4 * np.exp(c5*_z)
    		return model - bias_vals
	pb = scipy.optimize.leastsq(lsq, p0)[0]

	# Print best-fit coefficients
	print "-"*30
	print name[ID]
	print "-"*30
	print "Fitting coeffs."
	print "c1: %6.4f" % np.log10(p[0])
	print "c2: %6.4f" % p[1]
	print "c3: %6.4f" % p[2]
	print "c4: %6.4f" % pb[0]
	print "c5: %6.4f" % pb[1]
	print name[ID], '&', " & ".join(["%6.4f" % n for n in [ np.log10(p[0]), p[1], p[2], pb[0], pb[1]]])

	# Calculate cosmo. functions
	H, r, D, f = background_evolution_splines(cosmo)

	# Calculate binned number densities
	zbins = redshift_bins(dz=0.1)
	zc = np.array([0.5*(zbins[i] + zbins[i+1]) for i in range(zbins.size-1)])
	nz, vol, b = np.array( [n_bin(zbins[i], zbins[i+1], dndz_expt, bias_expt) 
                        for i in range(zbins.size-1)] ).T
	vol *= fsky
	k_NL = 0.2
	h= 0.67
	_kmax = k_max(zc,k_NL)# Output survey info
	print "-"*30
	print "zc      n Mpc^-3 h^3    bias     kmax    Vsur Mpc^3/h^3     dvdz      Srms"
	for i in range(zc.size):
    		#Szz = fluxrms[ID] * Scorr[ID]
    		#Szz = NU_LINE/(1.+zc[i]) * Szz if not Sconst[ID] else Szz
    		Szz = flux_redshift(zc[i])

    		print "%2.1f     %3.3e    %6.3f    %6.3f       %5.3e    %5.3e          %6.2f" % \
    		(zc[i],dndz_expt(zc[i]), b[i], _kmax[i], vol[i]*(h**3), dvdz(zc[i]), Szz),

    		#(zc[i], zbins[i], zbins[i+1], nz[i], b[i], vol[i]/1e9, nz[i]*vol[i], Szz),
    		if (nz[i]*vol[i]) < NGAL_MIN: print "*",
    		if Szz > Smax: print "#",
    		print ""
		#print len(zc) , len(b), len(vol), len(_kmax)
	data = np.concatenate((np.reshape(zc,(len(zc),1)),np.reshape(dndz_expt(zc),(len(zc),1)),np.reshape(b,(len(zc),1)),np.reshape(_kmax,(len(zc),1)),np.reshape(vol*(h**3),(len(zc),1)),np.reshape(dvdz(zc),(len(zc),1))),axis=1)
	np.savetxt('FotranCodes/'+name[ID]+'_Corrected_nz_Srms_2.txt',data)
	#data2 = np.concatenate((np.reshape(zc,(len(zc),1)),np.reshape(dndz_expt(zc),(len(zc),1)),np.reshape(b,(len(zc),1)) ,np.reshape(dvdz(zc),(len(zc),1))),axis=1)
	#np.savetxt('inputs/'+name[ID]+'_Corrected_nz_Srms_dndz_2.txt',data2)

	print "-"*30
	print "Ntot: %3.3e" % np.sum(nz * vol)
	print "fsky: %3.3f" % fsky
	_zmin = (NU_LINE*1e3 / numax[ID] - 1.)
	print "zmin: %3.3f" % (_zmin if _zmin >= 0. else 0.)
	print "zmax: %3.3f" % (NU_LINE*1e3 / numin[ID] - 1.)
	print "Srms const: %s" % Sconst[ID]
	print "Exp name", name[ID]
	print "-"*30
	print "\n"

	# Comparison plot of dndz, bias, and fitting function
	if DEBUG_PLOT:
    		P.subplot(211)
    		P.plot(_z, (dndz_expt(_z)),'r--', linewidth=2.0,label=name[ID])
    		P.plot(_z, (p[0] * _z**p[1] * np.exp(-p[2]*_z)),'b--',linewidth=2.0, label='fitting function')
    		P.ylabel('n(z)')
    		P.legend(loc='upper right')

    		P.subplot(212)
    		P.plot(_z, bias_expt(_z), 'r--', linewidth= 2.0)
    		P.plot(_z, pb[0] * np.exp(pb[1]*_z), 'b--', linewidth=2.0)
    		P.ylabel('b(z)')
    		#P.legend(loc='upper left')
    		P.savefig('outputs/'+name[ID]+'_dndz_also_bias_and_fitting_function_2.pdf')
    		P.show()
