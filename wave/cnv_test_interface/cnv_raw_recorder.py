"""
A script to record raw serial data to file
John Bird 30 Dec 2014
before running make sure you set the comm port to reflect your system
"""

import serial as serial
import time

comm_port = '/dev/ttyUSB0'
baud_rate = 9600
cnv_serial = serial.Serial(comm_port, baud_rate)

logfile = open('raw_cnv_data.bcnv', 'w')
buf = []

while True:
    buf.extend(cnv_serial.read(cnv_serial.inWaiting()))
    if len(buf) > 100:
        logfile.write(buf)
    time.sleep(0.01)
