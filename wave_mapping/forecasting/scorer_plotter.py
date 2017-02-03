#!/usr/bin/env python
import numpy
import scipy

import matplotlib.pyplot as plt

import parsers.gsd
import meteorology.sounding

sounding_getter = parsers.gsd.GSDParser()

sounding = sounding_getter.from_RUC_soundings(
    station='kunv',
    month=2,
    day=4,
    hour=15)

g = 9.806

z = numpy.linspace(sounding._z[0], sounding._z[-1], 100)
N = numpy.sqrt(g / sounding.theta(z)[0] * sounding.theta_gradient(z)[0])

U = sounding.wind_speed(z)[0]
dUdz = numpy.linalg.norm(sounding.wind_gradient(z)[0], axis=1)
dU2dz2 = numpy.gradient(dUdz, z)

l = numpy.sqrt(
    numpy.power(N/U, 2.0) - dU2dz2 / U)

plt.figure()
plt.plot(l, z)
plt.grid()
plt.ylim(0, 5000.0)
plt.xlabel('scorer parameter (1/m)')
plt.ylabel('altitude (m)')
plt.title('Scorer Parameter Profile')

#roughly figured from google earth, edge of plateau to top of bald eagle
lambda_ridge = 7.5e3

plt.figure()
plt.plot(2 * numpy.pi / l, z)
plt.plot([lambda_ridge, lambda_ridge], [z[0], z[-1]], 'k')
plt.grid()
plt.ylim(0, 5000.0)
plt.xlabel('wavelength (m)')
plt.ylabel('altitude (m)')
plt.title('Wavelength Profile')

plt.figure()
plt.plot(U, z)
plt.grid()
plt.ylim(0, 5000.0)
plt.xlabel('wind magnitude (m/s)')
plt.ylabel('altitude (m)')
plt.title('Wind Profile')

plt.figure()
plt.plot(N, z)
plt.grid()
plt.ylim(0, 5000.0)
plt.xlabel('B-V frequency (1/s)')
plt.ylabel('altitude (m)')
plt.title('Brunt-Vaisala Frequency Profile')

plt.show()
