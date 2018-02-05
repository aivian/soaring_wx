from sbxc_ground_control import AircraftStatePacket:    
from soaring_msgs.msg import AircraftState
 
def parse_state(recv_packet, pub_state, time_stamp):
	state_packet = AircraftStatePacket.AircraftStatePacket(recv_packet)
	
	msg = AircraftState()
	msg.header.stamp = time_stamp
	
	euler = state_packet.euler()
	msg.euler.x = euler[0]
	msg.euler.y = euler[1]
	msg.euler.z = euler[2]
	
	msg.v_ias = state_packet.v_ias()
	
	v_gps = state_packet.v_gps()
	msg.v_gps.x = v_gps[0]
	msg.v_gps.y = v_gps[1]
	msg.v_gps.z = v_gps[2]
	
	lla = state_packet.ac_lla()
	msg.ac_lla.x = ac_lla[0]
	msg.ac_lla.y = ac_lla[1]
	msg.ac_lla.z = ac_lla[2]
	
	pub_state.publish(msg)	
	
