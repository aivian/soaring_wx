#!/usr/bin/env python
import pdb
import os
import numpy
import scipy
import datetime
import subprocess

import matplotlib.pyplot as plt

import geodesy.conversions

import parsers.gsd
import parsers.bufrgruven
import meteorology.sounding

import parse_soundings

gsd = False

ruc_sounding_getter = parsers.gsd.GSDParser()
sounding_getter = parsers.bufrgruven.BufrGruvenParser()

load_stamp = {
    'year': 2018,
    'month': 11,
    'day': 06,
    'station': 'msv',
    'cycle': 12,
}
plot_stamp = {
    'year': 2018,
    'month': 11,
    'day': 07,
    'hour': 17,
    }

plot_lim = (0, 8000)

datestamp = '{}{:02d}{:02d}'.format(
    load_stamp['year'], load_stamp['month'], load_stamp['day'])

if not gsd:
    data_file = '/opt/bufrgruven/metdat/ascii/{}{:02d}_nam.prof.{}'.format(
        datestamp, load_stamp['cycle'], load_stamp['station'])

    if not os.path.isfile(data_file):
        call = '/opt/bufrgruven/run_bufr.sh k{} {} {}'.format(
            load_stamp['station'], datestamp, load_stamp['cycle'])
        subprocess.call(call, shell=True)
    print('loading: {}'.format(data_file))
    soundings = sounding_getter.parse_ascii(data_file)

    t_desired = datetime.datetime(
        plot_stamp['year'],
        plot_stamp['month'],
        plot_stamp['day'],
        plot_stamp['hour'],
        0, 0, 0)
    gps_t_desired = geodesy.conversions.datetime_to_gps(t_desired)
    best_idx = 0
    delta_min = numpy.inf
    for idx, sounding in enumerate(soundings):
        delta = numpy.abs(sounding.time - gps_t_desired)
        if delta < delta_min:
            best_idx = idx
            delta_min = delta
    sounding = soundings[best_idx]
else:
    sounding = ruc_sounding_getter.from_RUC_soundings(
        station='k{}'.format(load_stamp['station']),
    #    #lat=34.7,
    #    #lon=-95.1,
        month=plot_stamp['month'],
        day=plot_stamp['day'],
        hour=plot_stamp['hour'],
        model='Op40')

g = 9.806

z = numpy.linspace(sounding._z[0], sounding._z[-1], 100)
N = numpy.sqrt(g / sounding.theta(z)[0] * sounding.theta_gradient(z)[0])

wind = sounding.wind(z)[0]
base_wind = numpy.mean(wind[numpy.where(z<2000)[0],:], axis=0)
base_direction = base_wind / numpy.linalg.norm(base_wind)
wind_gradient = sounding.wind_gradient(z)[0]
#dUdz = numpy.linalg.norm(sounding.wind_gradient(z)[0], axis=1)
U = wind.dot(base_direction)
dUdz = wind_gradient.dot(base_direction)
dU2dz2 = numpy.gradient(dUdz, z)

theta = sounding.theta(z=z)[0]
dthetadz = sounding.theta_gradient(z)[0]

l = numpy.sqrt(
    numpy.power(N/U, 2.0) - dU2dz2 / U)
#numpy.save('l', l)
#numpy.save('z', z)

plt.figure()
plt.scatter(l, z)
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlim(0.0, 0.002)
plt.xlabel('scorer parameter (1/m)')
plt.ylabel('altitude (m)')
plt.title('Scorer Parameter Profile')

#roughly figured from google earth, edge of plateau to top of bald eagle
lambda_ridge = 7.5e3

f_lambda = plt.figure(figsize=(3,4))
plt.plot(2 * numpy.pi / l / 1000.0, z, linewidth=2)
#plt.plot([lambda_ridge, lambda_ridge], [z[0], z[-1]], 'k')
plt.xlim(0, 30.0)
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlabel('wavelength (km)')
plt.ylabel('altitude (m)')
plt.title('Wavelength Profile')
plt.tight_layout()

plt.figure()
plt.plot(dthetadz * 1000.0, z)
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlabel('lapse rate (K/km)')
plt.ylabel('altitude (m)')
plt.title('Lapse Rate Profile')

f_wind = plt.figure(figsize=(5,5))
plt.plot(U, z, '-r', label='magnitude', linewidth=2)
plt.plot(wind[:,0], z, '--g', label='u', linewidth=2)
plt.plot(wind[:,1], z, '--b', label='v', linewidth=2)
#plt.legend()
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlabel('wind magnitude (m/s)')
plt.ylabel('altitude (m)')
plt.title('Wind Profile')
plt.tight_layout()

plt.figure()
plt.scatter(N, z)
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlabel('B-V frequency (1/s)')
plt.ylabel('altitude (m)')
plt.title('Brunt-Vaisala Frequency Profile')

f_theta = plt.figure(figsize=(4,5))
plt.plot(theta, z, linewidth=2)
plt.grid()
plt.ylim(plot_lim[0], plot_lim[1])
plt.xlim(260, 360)
plt.xlabel(r'$\theta$ (K)')
plt.ylabel('altitude (m)')
plt.title('Potential Temperature Profile')
plt.tight_layout()

f_lambda.savefig('lambda.png', format='png', dpi=300)
f_theta.savefig('theta.png', format='png', dpi=300)
f_wind.savefig('wind.png', format='png', dpi=300)
plt.show()
