#!/usr/bin/env python
import os
import sys
import datetime

import numpy
import scipy

import netCDF4

import parsers.igc
import geodesy.conversions

# for track file format see
# https://www.unidata.ucar.edu/support/help/MailArchives/idv/msg00523.html
# for netcdf usage see
# https://unidata.github.io/netcdf4-python/#section1

def parse_flight(igc_path, flight_name):
    """Parse an IGC file into a netCDF for IDV

    Arguments:
        igc_path: path to the igc file to parse
        flight_name: name of the output file, will have date appended

    Returns:
        no returns
    """
    flight = parsers.igc.IGC(os.path.expanduser(igc_path))

    flight_track = netCDF4.Dataset(
        '{}_{}.nc'.format(
            flight_name, flight.start_date), 'w', format='NETCDF4')

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

if __name__ == '__main__':
    """The parse function can be run as a standalone script

    Call it with the input igc file and output file name as arguments
    Ex:
        ./parse_flight some_file.igc some_nc_file

    It will parse the igc into timestamped lat/lon/alt points and output
    a netcdf file some_nc_file_yyyy-mm-dd.nc

    where yyy-mm-dd is the datestamp of the beginning of the igc file
    """
    assert len(sys.argv) >= 3,\
        'must specify both an igc file and an output name'
    assert os.path.exists(os.path.expanduser(sys.argv[1])),\
        'specified input file not found'
    assert os.path.splitext(sys.argv[1])[1].lower() == '.igc',\
        'input file must be an IGC file'

    parse_flight(sys.argv[1], sys.argv[2])
