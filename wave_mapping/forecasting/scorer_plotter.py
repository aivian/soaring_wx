#!/usr/bin/env python
import pdb
import numpy
import scipy
import datetime

import matplotlib.pyplot as plt

import geodesy.conversions

import parsers.gsd
import parsers.bufrgruven
import meteorology.sounding

import parse_soundings

#sounding_getter = parsers.gsd.GSDParser()
sounding_getter = parsers.bufrgruven.BufrGruvenParser()

year = 2013
month = 11
day = 19
hour = 19

#sounding = sounding_getter.from_RUC_soundings(
#    station='kunv',
#    #lat=34.7,
#    #lon=-95.1,
#    month=1,
#    day=10,
#    hour=01,
#    model='Op40')

datestamp = '{}{:02d}{:02d}'.format(year, month, day)
data_file = '/opt/bufrgruven/metdat/ascii/{}00_nam.prof.abe'.format(
    datestamp)
print('loading: {}'.format(data_file))
soundings = sounding_getter.parse_ascii(data_file)

t_desired = datetime.datetime(year, month, day, hour, 0, 0, 0)
gps_t_desired = geodesy.conversions.datetime_to_gps(t_desired)
best_idx = 0
delta_min = numpy.inf
for idx, sounding in enumerate(soundings):
    delta = numpy.abs(sounding.time - gps_t_desired)
    if delta < delta_min:
        best_idx = idx
        delta_min = delta
sounding = soundings[best_idx]

g = 9.806

z = numpy.linspace(sounding._z[0], sounding._z[-1], 100)
N = numpy.sqrt(g / sounding.theta(z)[0] * sounding.theta_gradient(z)[0])

U = sounding.wind_speed(z)[0]
dUdz = numpy.linalg.norm(sounding.wind_gradient(z)[0], axis=1)
dU2dz2 = numpy.gradient(dUdz, z)

l = numpy.sqrt(
    numpy.power(N/U, 2.0) - dU2dz2 / U)
#numpy.save('l', l)
#numpy.save('z', z)

plt.figure()
plt.scatter(l, z)
plt.grid()
plt.ylim(0, 10000.0)
plt.xlim(0.0, 0.01)
plt.xlabel('scorer parameter (1/m)')
plt.ylabel('altitude (m)')
plt.title('Scorer Parameter Profile')

#roughly figured from google earth, edge of plateau to top of bald eagle
lambda_ridge = 7.5e3

plt.figure()
plt.plot(2 * numpy.pi / l, z)
#plt.plot([lambda_ridge, lambda_ridge], [z[0], z[-1]], 'k')
plt.grid()
plt.ylim(0, 10000.0)
plt.xlabel('wavelength (m)')
plt.ylabel('altitude (m)')
plt.title('Wavelength Profile')

plt.figure()
plt.plot(U, z)
plt.grid()
plt.ylim(0, 10000.0)
plt.xlabel('wind magnitude (m/s)')
plt.ylabel('altitude (m)')
plt.title('Wind Profile')

plt.figure()
plt.plot(N, z)
plt.grid()
plt.ylim(0, 10000.0)
plt.xlabel('B-V frequency (1/s)')
plt.ylabel('altitude (m)')
plt.title('Brunt-Vaisala Frequency Profile')

plt.show()
