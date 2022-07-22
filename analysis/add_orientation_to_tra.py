import gzip
import numpy as np
import math
import sys


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def write_ori_to_tra(ori_file,tra_file,ori_tra_file,tra_gz=True):
    #ori_path = "/Users/jacqueline/AllData_160517-160702_506-516keV.ori"
    time,Xthetalat,Xphilong,Zthetalat,Zphilong = np.genfromtxt(ori_file,skip_header=3,skip_footer=1,usecols=(1,2,3,4,5),unpack=True,delimiter=" ",dtype=str)
    time = np.array([float(x) for x in time])
    Xthetalat = np.array([float(x) for x in Xthetalat])
    Xphilong = np.array([float(x) for x in Xphilong])
    Zthetalat = np.array([float(x) for x in Zthetalat])
    Zphilong = np.array([float(x) for x in Zphilong])
    if (tra_gz == False):
        tra = open(tra_file,"rt")
        ori_tra = open(ori_tra_file,"wt") 
    else:
        tra = gzip.open(tra_file,"rt")
        ori_tra = gzip.open(ori_tra_file,"wt")
    # tra = gzip.open("ActivationStep3_positron_46days_10dets.p1_whole_flight_times.tra.gz","rt")
    # ori_tra = gzip.open("ActivationStep3_positron_46days_10dets.p1_whole_flight_times_ori.tra.gz","wt")

    for line in tra:
        line = line.rstrip()
        if not line.startswith("TI"):
        	print(line,file=ori_tra)
        if line.startswith("TI"):
                print(line,file=ori_tra)
                input_time = float(line.strip(" ")[3:])
                idx = find_nearest(time,input_time)
                
                closest_Xlat = Xthetalat[idx]
                closest_Xlon = Xphilong[idx]
                closest_Zlat = Zthetalat[idx]
                closest_Zlon = Zphilong[idx]
                print("GX",closest_Xlon,closest_Xlat,file=ori_tra)
                print("GZ",closest_Zlon,closest_Zlat,file=ori_tra)

if __name__ == '__main__':

    if len(sys.argv) < 5:
        print("")
        print("Usage: python add_orientation_to_tra.py ori_file tra_file ori_tra_file")
        print("")
        print("- Finds the time in a given .ori file which is closest to the time in a given .tra.gz file event. \n Then, it writes the X and Z orientations from the .ori file at the closest time as 'GX' and 'GZ' lines in a new .tra.gz file.")
        print("")
        print("ori_file: required. Path of the .ori file with the orientations you need. Times in tra_file are compared to times in this .ori file.")
        print("tra_file: required. Path of the .tra.gz file without GX, GZ orientations")
        print("ori_tra_file: required. Path of a .tra.gz file to be written with GX, GZ orientations")
        print("tra_gz: optional. Default tra_gz=True to indicate input and output .tra.gz files. Set to False for a .tra file")
        print("")
        print("Note: .ori file is OG time Xthetalat Xphilong Zthetalat Zphilong \n .tra file is GX Xlong Xlat, GZ Zlong Zlat")
        print("")
        sys.exit(0)
    
    ori_file = sys.argv[1]
    tra_file = sys.argv[2]
    ori_tra_file = sys.argv[3]
    tra_gz = sys.argv[4]

    write_ori_to_tra(ori_file,tra_file,ori_tra_file,tra_gz)

