#import AVIaParser as avia_iface
#import AVIaPacket as avia
import AircraftStatePacket
from soaring_msgs.msg import AircraftState

class FlightDataSharingCallbacks:
    def __init__( self, fdata_server ):
        self.fdata_server = fdata_server
        
    def aircraft_state( self, msg ):
        # callback on aircraft state received
        # only continue for every 3rd msg received (because there are a lot of aircraft state messages)
        if ( ( msg.header.seq % 3 ) == 0 ):
            # declare an avia aircraft state packet
            pkt = AircraftStatePacket.AircraftStatePacket()
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
