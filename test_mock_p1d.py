import numpy as np
import matplotlib.pyplot as plt
import lya_mock_p1d as mock

# number of cells (power of two)
N2 = 15
# cell width (in km/s)
dv_kms=10
# whether to use white noise power or not
white_noise=False
# random seed
seed=555

# generate a mock maker
mock_maker=mock.MockMaker(N2,dv_kms,seed,white_noise)

# get redshift for each cell
z = mock_maker.get_redshifts()

# get Gaussian field
delta, var_delta = mock_maker.get_gaussian_field()
print('mean delta =', np.mean(delta))
print('var delta =', np.var(delta))
print('expected var delta =',var_delta)
plt.plot(z,delta)
plt.xlabel('z')
plt.ylabel('Gaussian field')
plt.show()

# from Gaussian field to lognormal density
density = mock_maker.get_density(var_delta,z,delta)
print('mean density =', np.mean(density))
print('var density =', np.var(density))
plt.semilogy(z,density)
plt.xlabel('z')
plt.ylabel('density')
plt.show()

# from lognormal density to optical depth
tau = mock.get_tau(z,density)
print('mean tau =', np.mean(tau))
print('var tau =', np.var(tau))
plt.semilogy(z,tau)
plt.xlabel('z')
plt.ylabel('optical depth')
plt.show()

# from optical depth to flux
flux = np.exp(-tau)
print('mean flux =', np.mean(flux))
print('var flux =', np.var(flux))
plt.plot(z,flux)
plt.xlabel('z')
plt.ylabel('transmitted flux fraction')
plt.xlim(3.0,3.1)
plt.show()

# generate multiple mock spectra
for i in range(3):
    # get everything on one shot
    wave, flux = mock_maker.get_lya_skewer()
    plt.plot(wave,flux)
    plt.xlabel('wave')
    plt.ylabel('transmitted flux fraction')
    plt.xlim(4000,4200)
    plt.show()

