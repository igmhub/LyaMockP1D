import numpy as np

# code to make mock Lya spectra following McDonald et al. (2006)
# copied from c++ code in Cosmology/LNLyaF

def power_amplitude(z):
  """adds redshift evolution to the 1D power spectrum"""
  return 58.6*pow((1+z)/4.0,-2.82)

def tau_amplitude(z):
  """adds redshift evolution to the mean optical depth (tau)"""
  return 0.374*pow((1+z)/4.0,5.10)

def power_kms(z,k_kms,dv_kms=10.0,simple=False):
  """at a fixed redshift z, compute P1D at different wavenumbers k (in s/km)"""
  # option for debugging (return constant power)
  if simple: return np.ones_like(k_kms)*100.0
  # power used to make mocks in from McDonald et al. (2006)
  A = power_amplitude(z)
  k1 = 0.001
  n = 0.7
  R1 = 5.0
  # compute term without smoothing
  P = A * (1.0+pow(0.01/k1,n)) / (1.0+pow(k_kms/k1,n))
  # smooth with Gaussian and top hat
  kdv = np.fmax(k_kms*dv_kms,0.000001)
  P *= np.exp(-pow(k_kms*R1,2)) * pow(np.sin(kdv/2)/(kdv/2),2)
  return P

def get_density(z_c,var_delta,z,delta):
  """transform Gaussian field delta to lognormal density, at each z"""
  tau_pl=2.0
  # relative amplitude 
  rel_amp = power_amplitude(z)/power_amplitude(z_c)
  return np.exp(tau_pl*(delta*np.sqrt(rel_amp)-0.5*var_delta*rel_amp))

def get_tau(z,density):
  """transform lognormal density to optical depth, at each z"""
  A = tau_amplitude(z)
  return A*density

def get_flux(z_c,tau):
  """transform optical depth to transmitted flux fraction"""
  return np.exp(-tau)

def get_redshifts(z_c,N2,dv_kms):
  """get redshift for each cell in the array (centered at z_c)"""
  N = np.power(2,N2)
  L_kms = N * dv_kms
  c_kms = 2.998e5
  if (L_kms > 4 * c_kms):
    print('Array is too long, approximations break down.')
    exit()
  # get indices
  i = range(N)
  z = (1+z_c)*pow(1-(i-N/2+1)*dv_kms/2.0/c_kms,-2)-1
  return z

def get_gaussian_field(z_c, N2, dv_kms, seed=666, simple=False):
  """generate Gaussian field at redshift z_c"""
  # length of array
  N = np.power(2,N2)
  # number of Fourier modes
  NF=int(N/2+1)
  # setup random number generator
  gen = np.random.RandomState(seed)
  # generate random Fourier modes
  modes = np.empty(NF,dtype=complex)
  modes[:].real = gen.normal(size=NF)
  modes[1:-1].imag = gen.normal(size=NF-2)
  # get frequencies (wavenumbers in units of s/km)
  k_kms = np.fft.rfftfreq(N)*2*np.pi/dv_kms
  # get power evaluated at each k
  P_kms = power_kms(z_c,k_kms,simple)
  # normalize to desired power
  modes[-1:1].real *= np.sqrt(P_kms[-1:1])
  modes[-1:1].imag = 0
  modes[1:-1].real *= np.sqrt(0.5*P_kms[1:-1])
  modes[1:-1].imag *= np.sqrt(0.5*P_kms[1:-1])
  # inverse FFT to get (normalized) delta field
  delta = np.fft.irfft(modes) * np.sqrt(N/dv_kms)
  # compute also expected variance, will be used in lognormal transform
  dk_kms = 2*np.pi/(N*dv_kms)
  var_delta=np.sum(P_kms)*dk_kms/np.pi
  # Nyquist frecuency is counted twice in variance, and it should not be
  var_delta *= NF/(NF+1)
  return delta, var_delta

def get_lya_skewer(z_c=3.0, N2=15, dv_kms=10.0, seed=666):
  """directly compute flux skewer, using functions above"""
  z = get_redshifts(z_c,N2,dv_kms)
  delta, var_delta = get_gaussian_field(z_c,N2,dv_kms,seed)
  #var_delta = np.var(delta)
  density = get_density(z_c,var_delta,z,delta)
  tau = get_tau(z,density)
  flux = get_flux(z_c,tau)
  wave = 1215.67*(1+z)
  return wave, flux


