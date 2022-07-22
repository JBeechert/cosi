# Split an orientation file in to given number of files.


import numpy as np

nSplits = 30

times,aspect1,aspect2,aspect3,aspect4 = np.genfromtxt("OP_160612-14.10secgap.ori",skip_header = 3, skip_footer = 1, unpack=True, usecols=(1,2,3,4,5))

length = len(times)/nSplits

endTimes = []

for i in range(nSplits):
	fout = open("OP_160612-14.10secgap.split"+str(i)+".ori","w")
	fout.write("\nType OrientationsGalactic\n\n")

	index = i*length

	if i == nSplits-1:
		length += len(times)%nSplits

	lastTime = 0
	for j in range(length):
		fout.write("OG %f %f %f %f %f\n"%(times[index+j],aspect1[index+j],aspect2[index+j],aspect3[index+j],aspect4[index+j]))
		lastTime = times[index+j] #intentionally overwritten every time

	endTimes.append(lastTime+1)

	fout.write("EN\n\n")
	fout.close()





