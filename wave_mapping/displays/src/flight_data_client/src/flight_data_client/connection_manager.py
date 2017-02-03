import sys
from PyQt4 import QtGui, QtCore
import time
import copy

import rospy

#from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as pgl
import numpy as np
from pyqtgraph.dockarea import *

# import protocol
from communications.AVIa_tcp import AVIaTCPParser as AVIaParser
from communications.socket_connections import pysock

# import serial packets
from communications.AVIa_packet import AVIaPacket
from wave_mapping_packets.AircraftStatePacket import AircraftStatePacket
from wave_mapping_packets.SplineModelPacket import SplineModelPacket

# import widgets
import graphics.plot_widgets as plot_widgets
from graphics.glider3d import glider3d

# more stuff
from geometry.quaternion import Quaternion
from geometry.bSpline import TensorProductSpline

class SocketServicer(QtCore.QObject): 
    """ class to service the comms from the data server
    
    really, this should exist in a separate node, but if I put it
    here it'll be easier to build into a standalone application
    
    qt signals are implemented as class variables...because that was
    how I could get them to work... I need to learn more about this...
    """
    newPkt = QtCore.pyqtSignal(AVIaPacket)
    newIAS = QtCore.pyqtSignal(float, list)
    newWind = QtCore.pyqtSignal(float, np.ndarray, np.ndarray)
    newEuler = QtCore.pyqtSignal(float, np.ndarray)
    newSpline = QtCore.pyqtSignal(np.ndarray)
    newMap = QtCore.pyqtSignal(float, np.ndarray)

    def __init__(self, addr=None):
        """ constructor
        
        calls super-class constructor (QObject) then sets up connection
        to the flight data server. At the moment this can't recover
        well from disconnections, but I'll take care of that later
        """
        super(SocketServicer, self).__init__()
        
        self.addr = 'www.cathartes.org'
        self.addr = '10.0.0.12'
        if addr is not None:
            self.addr = addr
        self.connection = self.connect_to_server()
        self.connection._conn = self.connection.sock
        self.parser = AVIaParser(self.connection)
        self.newPkt.connect(self.parse_packet)
    
    
    def connect_to_server(self):
        server_socket = pysock(
            self.addr,
            5006,
            1024)
        server_socket.start_send()
        print 'connected'
        return server_socket
        
    def service_comms(self):
        """ service the comms
        
        look for new packets in the buffer and parse them out if there
        are any
        
        No Args
        No Returns, but will emit a newPkt signal for each new packet
        """
        new_packets = self.parser.read_buffer()
        if new_packets:
            for recv_packet in self.parser.packet_buffer:
                self.newPkt.emit(recv_packet)
        self.parser.packet_buffer = []
        
    def parse_packet(self, pkt):
        if pkt.packet_type() == 0x10:
            state_pkt = AircraftStatePacket(pkt)
            now = time.clock()
            self.newIAS.emit(now, [state_pkt.v_ias()])
            
            airV_b = np.zeros((3,))
            airV_b[0] = state_pkt.v_ias()
            
            euler = state_pkt.euler()
            q_bi = Quaternion()
            q_bi.from_euler(euler)
            
            airV_i = q_bi.rot(airV_b,False)
            iV_i = np.asarray(state_pkt.v_gps())
            
            wind = airV_i - iV_i
            lla = state_pkt.ac_lla()
            
            self.newWind.emit(now, wind, np.asarray(lla))
            self.newEuler.emit(now, np.asarray(euler))
            
        if pkt.packet_type() == 0x11:
            spline_pkt = SplineModelPacket(pkt)            
            now = time.clock()
            c = np.array(spline_pkt.coefficients())
            self.newSpline.emit(c)
        
