from struct import *
from communications.AVIa_packet import AVIaPacket

class AircraftStatePacket(AVIaPacket):

    def __init__(self,pkt='NaN'):
        if pkt == 'NaN':
            super(AircraftStatePacket,self)
            self.pkt_type = 0x10
            self.payload = ['0']*40
        else:
            super(AircraftStatePacket,self)
            self.pktType(pkt.pktType())
            self.payload = pkt.payload
            
    def ac_lla(self, setVal='NaN'):
        if (setVal == 'NaN'):
            lla = []
            for i in range(0,3):
                lla.append(unpack('f',self.payload[0+4*i:4+4*i])[0])
            return lla
        else:
            # expecting setVal to be a list of three floats [lat, lon, alt]
            self.payload[0:4] = pack('f',setVal[0])
            self.payload[4:8] = pack('f',setVal[1])
            self.payload[8:12] = pack('f',setVal[2])
    
    def euler(self, setVal='NaN'):
        if (setVal == 'NaN'):
            ea = []
            for i in range(0,3):
                ea.append(unpack('f',self.payload[12+4*i:16+4*i])[0])
            return ea
        else:
            # expecting setVal to be a list of three floats [phi, theta, psi]
            self.payload[12:16] = pack('f',setVal[0])
            self.payload[16:20] = pack('f',setVal[1])
            self.payload[20:24] = pack('f',setVal[2])
    
    def v_gps(self, setVal='NaN'):
        if (setVal == 'NaN'):
            v = []
            for i in range(0,3):
                v.append(unpack('f',self.payload[24+4*i:28+4*i])[0])
            return v
        else:
            # expecting setVal to be a list of three floats [phi, theta, psi]
            self.payload[24:28] = pack('f',setVal[0])
            self.payload[28:32] = pack('f',setVal[1])
            self.payload[32:36] = pack('f',setVal[2])
            
    def v_ias(self, setVal='NaN'):
        if (setVal == 'NaN'):
            v_ias = unpack('f',self.payload[36:40])[0]
            return v_ias
        else:
            self.payload[36:40] = pack('f',setVal)
