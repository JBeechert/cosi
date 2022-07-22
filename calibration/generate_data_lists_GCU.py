# Adapted by J. Beechert from a script by C. Sleator

import os
import sys
#import pickle
import glob
import numpy as np
from datetime import datetime


# Take start/end dates/times as arguments in the function
def generate_data_list(startdatetime,stopdatetime,directory,source,length):
	startdate,starttime = startdatetime.strip().split("_")
	stopdate,stoptime = stopdatetime.strip().split("_") 	

	path_to_dat = "/mnt/DATA_C/"
	potential_files = sorted(glob.glob(path_to_dat+"COSIc"+startdate+"*dat"))
        
        if startdate != stopdate:
                dates = [date for date in range(int(startdate)+1,int(stopdate)+1)]
                for date in dates:
		    #potential_files.extend(sorted(glob.glob(path_to_dat+"COSIc"+stopdate+"*dat")))
                    potential_files.extend(sorted(glob.glob(path_to_dat+"COSIc"+str(date)+"*dat")))

	startdate = int(startdate)
	stopdate = int(stopdate)
	starttime = int(starttime)
	stoptime = int(stoptime)
        
	files_for_obs = []
	for f in potential_files:
                         
		f = f.split("/")[len(f.split("/"))-1]
	
		time = f.split("__")[1].strip(".dat")
		time = int(time)
		date = f.split("__")[0].strip("COSIc")
		date = int(date)

	
		if startdate == stopdate:
			if time >= starttime and time <= stoptime:
				files_for_obs.append(path_to_dat+f)
		else:
			if (date == startdate and time >= starttime) or (date == stopdate and time <= stoptime):
				files_for_obs.append(path_to_dat+f)      
                        
                        if (date > startdate and date < stopdate):
                            files_for_obs.append(path_to_dat+f)
        
	if length != None:
		file_list_short = []
		total_tdiff = 0
		length_seconds = 60*length
		file_list_short.append(files_for_obs[0])
		for i in range(1,len(files_for_obs)):
			fprev = file_list_short[i-1].split("/")[len(file_list_short[i-1].split("/"))-1]
			f = files_for_obs[i].split("/")[len(files_for_obs[i].split("/"))-1]
	
			date_prev = fprev.split("__")[0].strip("COSIc")
			time_prev = fprev.split("__")[1].strip(".dat")
			date = f.split("__")[0].strip("COSIc")
			time = f.split("__")[1].strip(".dat")
	
			time_prev_dt = datetime.strptime(date_prev+time_prev,"%y%m%d%H%M%S")
			time_dt = datetime.strptime(date+time,"%y%m%d%H%M%S")
	
			tdiff = (time_dt-time_prev_dt).seconds

			file_list_short.append(path_to_dat+f)
			total_tdiff += tdiff

			if total_tdiff >= length_seconds:
				break

		files_for_obs = file_list_short
	
		print("time in file list:",total_tdiff/60.,"minutes")
	
	if directory[len(directory)-1] != "/":
		directory += "/"
	
	#fout_name = directory+calib_dict["Run"+run]["source"]+"_Run"+run
	fout_name = directory+"COSI_GCUcalib-"+startdatetime+"-"+stopdatetime+"-"+source	
	if length != None:
		fout_name += "_"+str(length)+"mins"
	fout_name += ".dat"
	
	fout = open(fout_name,"w")
        for f in files_for_obs:
		fout.write(f+'\n')
	
	fout.close()
		
        return fout_name

if __name__ == '__main__':

	if len(sys.argv) < 5:
		print("")
		print("Usage: python generate_data_lists_GCU.py Startdatetime Stopdatetime Directory Source Time")
		print("")
		#print "Run#: required, self explanatory"
		print("Startdatetime (UTC): required, beginning of desired observation. YYMMDD_HHMMSS")
		print("Stopdatetime (UTC): required, end of desired observation. YYMMDD_HHMMSS")
                print("Directory (string): where you want the file list to go")
                print("Source (string): Atomic symbol-mass number")
		print("Time: *optional argument*. Time in minutes that you want the observation to be. Default: entire calibration run")
		print("")
		sys.exit(0)
	
	startdatetime = sys.argv[1]
	stopdatetime = sys.argv[2]
        directory = sys.argv[3]
	source = sys.argv[4]
	length = None
	if len(sys.argv) == 6:
		length = int(sys.argv[5])

	generate_data_list(startdatetime,stopdatetime,directory,source,length)



