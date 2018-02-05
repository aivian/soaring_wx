from struct import *
from communications.AVIa_packet import AVIaPacket

class SplineModelPacket(AVIaPacket):

    def __init__(self, pkt=None):
        if pkt is not None:
            super(SplineModelPacket,self)
            self.pktType(pkt.packet_type())
            self.payload = pkt.payload
        else:
            super(SplineModelPacket, self)
            self.pkt_type = 0x11
            self.payload = ['0']*1296
            
    def coefficients(self, setVal=None):
        """ set/retrieve spline model coefficients
        
        Args:
            setVal: optionally, numpy array of coefficients
        
        Returns:
            coefficients: if setVal is not specified then returns the 
                contents of the payload interpreted as for an 18x18
                model stored as single precision floats
        """
        if setVal is not None:
            byte_count = 0
            for i in setVal:
                self.payload[byte_count:byte_count+4] = pack('f', i)
                byte_count += 4
        else:
            return unpack('f'*324, self.payload[0:1296])
