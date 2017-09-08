#!/usr/bin/env python
import os
import datetime

import numpy
import scipy

import netCDF4

import parsers.igc
import geodesy.conversions

# Right now this is just to work out how to parse a flight for idv. It is set
# up to parse my flight from 2017-02-17

# for track file format see
# https://www.unidata.ucar.edu/support/help/MailArchives/idv/msg00523.html
# for netcdf usage see
# https://unidata.github.io/netcdf4-python/#section1

info = ('795v2vy4.igc', 'Tony-2017-09-05')
# info = (igc_file_path, flight_name)

flight = parsers.igc.IGC(os.path.expanduser(info[0]))

flight_track = netCDF4.Dataset(
    '{}_{}.nc'.format(info[1], flight.start_date), 'w', format='NETCDF4')

# Create the dimensions
flight_track.createDimension('time', None)

base_time = geodesy.conversions.datetime_to_gps(datetime.datetime(
    flight.start_year, flight.start_month, flight.start_day, 0, 0, 0, 0))
time = flight.time
lla = flight.lla()

interp_time = numpy.arange(time[0], time[-1], 10.0)
lla_interpolant = scipy.interpolate.interp1d(time, lla, axis=0)
interp_lla = lla_interpolant(interp_time)

times = flight_track.createVariable('time', 'f8', ('time',))
times.units = 'seconds since {} 00:00:00.0'.format(flight.start_date)
times[:] = interp_time - base_time

latitude = flight_track.createVariable('latitude', 'f8', ('time',))
latitude.units = 'degrees_N'
latitude[:] = numpy.rad2deg(interp_lla[:,0])

longitude = flight_track.createVariable('longitude', 'f8', ('time',))
longitude.units = 'degrees_E'
longitude[:] = numpy.rad2deg(interp_lla[:,1])

altitude = flight_track.createVariable('altitude', 'f8', ('time',))
altitude.units = 'm'
altitude[:] = interp_lla[:,2]

flight_track.close()

