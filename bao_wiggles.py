#!/usr/bin/python
"""
Separate P(k) from CAMB into a smooth reference power spectrum and a BAO 
wiggles function, using the method in Bull et al. (2014), arXiv:1405.1452.
  -- Phil Bull (2014) <philbull@gmail.com>
"""
import numpy as np
import pylab as P
import scipy.interpolate

C = 3e5 # Speed of light, km/s
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



def spline_pk_nobao(k_in, pk_in, kref=[2.15e-2, 4.5e-1]):
    """
    Construct a smooth power spectrum with BAOs removed, and a corresponding 
    BAO template function, by using a two-stage splining process.
    """
    # Get interpolating function for input P(k) in log-log space
    _interp_pk = scipy.interpolate.interp1d( np.log(k_in), np.log(pk_in), 
                                        kind='quadratic', bounds_error=False )
    interp_pk = lambda x: np.exp(_interp_pk(np.log(x)))
    
    # Spline all (log-log) points except those in user-defined "wiggle region",
    # and then get derivatives of result
    idxs = np.where(np.logical_or(k_in <= kref[0], k_in >= kref[1]))
    _pk_smooth = scipy.interpolate.UnivariateSpline( np.log(k_in[idxs]), 
                                                    np.log(pk_in[idxs]), k=3, s=0 )
    pk_smooth = lambda x: np.exp(_pk_smooth(np.log(x)))

    # Construct "wiggles" function using spline as a reference, then spline it 
    # and find its 2nd derivative
    fwiggle = scipy.interpolate.UnivariateSpline(k_in, pk_in / pk_smooth(k_in), k=3, s=0)
    derivs = np.array([fwiggle.derivatives(_k) for _k in k_in]).T
    d2 = scipy.interpolate.UnivariateSpline(k_in, derivs[2], k=3, s=1.0) #s=1.
    # (s=1 to get sensible smoothing)
    
    # Find maxima and minima of the gradient (zeros of 2nd deriv.), then put a
    # low-order spline through zeros to subtract smooth trend from wiggles fn.
    wzeros = d2.roots()
    wzeros = wzeros[np.where(np.logical_and(wzeros >= kref[0], wzeros <= kref[1]))]
    wzeros = np.concatenate((wzeros, [kref[1],]))
    wtrend = scipy.interpolate.UnivariateSpline(wzeros, fwiggle(wzeros), k=3, s=0)
    
    # Construct smooth "no-bao" function by summing the original splined function and 
    # the wiggles trend function
    idxs = np.where(np.logical_and(k_in > kref[0], k_in < kref[1]))
    pk_nobao = pk_smooth(k_in)
    pk_nobao[idxs] *= wtrend(k_in[idxs])
    fk = (pk_in - pk_nobao)/pk_nobao
    
    # Construct interpolating functions
    ipk = scipy.interpolate.interp1d( k_in, pk_nobao, kind='linear',
                                      bounds_error=False, fill_value=0. )
    ifk = scipy.interpolate.interp1d( k_in, fk, kind='linear',
                                      bounds_error=False, fill_value=0. )
    return ipk, ifk

# Calculate cosmo. functions
H, r, D, f = background_evolution_splines(cosmo)
# Load P(k) data from CAMB (IMPORTANT: P(k) should have lots of k samples!)
k_full, pk_full = np.genfromtxt("inputs/linear_matterpower_1.dat").T
# Construct interpolation functions for smooth reference P(k), and BAO wiggles 
# function, f_BAO(k)
h= 0.67

kref=np.linspace(2.15e-2/h, 4.5e-1/h, 10)
pk_ref, fbao = spline_pk_nobao(k_full, pk_full, kref=[2.15e-2/h, 4.5e-1/h])
#print kref
print 'pk(0.2)=', pk_ref(0.2)*(1+fbao(0.2))
print D(0)
test_sin = k_full**2
data2 = np.concatenate((np.reshape(k_full,(len(k_full),1)),np.reshape(test_sin,(len(k_full),1))),axis=1)
data = np.concatenate((np.reshape(k_full,(len(k_full),1)),np.reshape(fbao(k_full),(len(k_full),1)),np.reshape(pk_ref(k_full),(len(k_full),1))),axis=1)
#print data[:,0], data[:,1]
#np.savetxt('bao_wiggles_test_sin.txt',data2)#fmt='%1.20f')
np.savetxt('bao_wiggles_powerspectrum.txt',data)#fmt='%1.20f')
#np.savetxt('bao_wiggles_powerspectrum.txt',(k_full, pk_ref(k_full), fbao(k_full)))
# Plot wiggles function
fig = P.figure(figsize=(8,8), dpi=100)
P.subplot(211)
P.plot(k_full, fbao(k_full))
P.xlim(7e-3, 1e+0)
P.ylim(-0.15, 0.15)
P.xscale('log')
#P.yscale('log')
# Plot full P(k), P_ref(k), and P_ref(k) * (1 + f_BAO(k))
P.subplot(212)
P.plot(k_full, pk_full, 'k-', lw=1.5, label = 'P(k)') # P(k) from CAMB
#P.plot(k_full, pk_ref(k_full), 'b-', lw=1., label= 'P(k) ref') # P_ref(k)
P.plot(k_full, pk_ref(k_full) * (1. + fbao(k_full)), 'y--', lw=1.5, label = 'P(k) ref from bao function') # P_ref(k)
P.legend(loc='upper right')
P.yscale('log')
P.xscale('log')
P.xlim((7e-3, 1e+0))
P.ylim((9e2, 3e4))
P.savefig('test_wiggles_func.eps')
P.show()
