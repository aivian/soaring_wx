# CenteringParams.msg
# message contains data about the MacCready ring state

std_msgs/Header  header
float32 setting		# absolute ring setting (expected climb rate, m/s)
float32 threshold	# approximate threshold wind speed to achieve climb rate (m/s)