from wave_mapping_packets.AircraftStatePacket import AircraftStatePacket
from wave_mapping_packets.SplineModelPacket import SplineModelPacket
from wave_mapping_msgs.msg import AircraftState
from wave_mapping_msgs.msg import SplineModel

from communications import AVIa_tcp
from communications import socket_connections

import rospy
import numpy as np
import socket
import errno
import sys

class FlightDataSharingCallbacks:
    def __init__( self, fdata_server ):
        self.fdata_server = fdata_server
        
    def aircraft_state( self, msg ):
        # callback on aircraft state received
        # only continue for every 3rd msg received (because there are a lot of aircraft state messages)
        if ( ( msg.header.seq % 3 ) == 0 ):
            # declare an avia aircraft state packet
            pkt = AircraftStatePacket()
            # include aircraft location
            pkt.ac_lla( [ msg.ac_lla.x, msg.ac_lla.y, msg.ac_lla.z ] )
            # include inertial velocities
            pkt.v_gps( [ msg.v_gps.x, msg.v_gps.y, msg.v_gps.z ] )
            # include euler state
            pkt.euler( [ msg.euler.x, msg.euler.y, msg.euler.z ] )
            # include airspeed
            pkt.v_ias( msg.v_ias )
            # write full packet out on tcp
            self.fdata_server.sock.send(pkt.fulPkt())
            
class DataTerminalConnection:
    def __init__(self):
        self.client = None
        self.connected = False
        
    def prepare_connection(self):
        self.client = socket_connections.pysock(
            '',
            5006,
            1024)
        self.client.start_listen()
        
    def wait_for_client(self):
        self.client.start_recv()
        rospy.loginfo(self.client.remote_addr)
        return self.client

    
    def state_sender(self, state_msg):
        """ send a state message to the client
        
        takes a state message and packs it up to send to the client
        
        Args:
            state_msg: AircraftState ROS message to send
            
        Returns:
            no returns
        """
        state_pkt = AircraftStatePacket()
    
        state_pkt.ac_lla([state_msg.ac_lla.x, state_msg.ac_lla.y, state_msg.ac_lla.z])
        state_pkt.v_gps([state_msg.v_gps.x, state_msg.v_gps.y, state_msg.v_gps.z])
        state_pkt.euler([state_msg.euler.x, state_msg.euler.y, state_msg.euler.z])
        state_pkt.v_ias(state_msg.v_ias)
        
        self._send_to_client(state_pkt)
        
    def spline_sender(self, spline_msg):
        """ send a spline message to client
        
        takes a spline message and packs it up to send to client
        
        Args:
            spline_msg: SplineModel ROS mesage
            
        Returns:
            no returns
        """
        spline_pkt = SplineModelPacket()
        spline_pkt.coefficients(spline_msg.c)
        np.save('/home/bird/test.npy', spline_pkt.payload)
        print('sending spline')
        self._send_to_client(spline_pkt)
        
    def _send_to_client(self, pkt):
        """ sends a packet to a client
        
        this does the actual sending and deals with closed connections
        and whatnont
        """
        if self.connected:
            try:
                self.client._conn.send(pkt.fulPkt())
            except socket.error, e:
                self.connected = False
