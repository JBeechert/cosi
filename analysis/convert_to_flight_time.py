import sys
import gzip

# Took the times from the 10 det .tra.gz files. Almost identical to 9 det times
# Example of equation:
#  slope = (flight_final - flight_inital) / (activ_final - activ_inital)
#  flight_initial = slope*activ_inital + b
#  flight_time = slope*activ_time + b
def convert_to_flight_time(tra_time,particle):

	# for activation .tra.gz files
	if (particle == "alpha"):
		flight_time = 50.74949*tra_time + 1463443196.1237
	if (particle == "neutron"):
		flight_time = 50.749809*tra_time + 1463443201.2555
	if (particle == "positron"):
		flight_time = 50.7499*tra_time + 1463443175.6877
	if (particle == "proton"):
		flight_time = 50.74944*tra_time + 1463443213.4404

	# for Ling .tra file
	if (particle == "ling"):
		flight_time = 46.688*tra_time + 1463443214.70457

	return flight_time

def write_flight_time_tra(tra_file,flight_time_tra_file,particle):
	# Ling files isn't compressed
	if (particle == "ling"):
		tra = open(tra_file,"r")
		flight_time_tra = open(flight_time_tra_file,"w")

	else:	
		tra = gzip.open(tra_file,"rt")
		flight_time_tra = gzip.open(flight_time_tra_file,"wt")

	for line in tra:
		line = line.rstrip()
		if not line.startswith("TI"):
			print(line,file=flight_time_tra)
		if line.startswith("TI"):
			time = float(line.strip(" ")[3:])
			flight_time = convert_to_flight_time(time,particle)
			print("TI "+str(flight_time),file=flight_time_tra)

if __name__ == '__main__':

	if len(sys.argv) < 4:
		print("")
		print("Usage: python convert_to_flight_time.py tra_file flight_time_tra_file particle")
		print("")
		print("- Converts the times in a simulated .tra.gz file to the times in .tra.gz files from the 2016 flight.")
		print("- Anyone can adapt this script to their needs by editing 'particles' and \n the equations that convert the times inside the code.")
		print("")
		print("tra_file: required. Path of the .tra.gz file whose times are read in for conversion to flight times")
		print("flight_time_tra_file: required. Path of a .tra.gz file to be written with times converted to flight times")
		print("particle: required. Type of activation particle or Ling in tra_file: alpha, neutron, positron, proton, or ling")
		print("")
		sys.exit(0)
	
	tra_file = sys.argv[1]
	flight_time_tra_file = sys.argv[2]
	particle = sys.argv[3]

	write_flight_time_tra(tra_file,flight_time_tra_file,particle)
