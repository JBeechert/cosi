import numpy as np
import matplotlib.pyplot as plt

def parse_peaks(ecal_fname):
  f = open(ecal_fname,"r")
  lines = f.readlines()

  fwhm_energy_ls = []
  fwhm_percent_ls = []

  for line in lines:
    line = line.split(' ')
    if line[0] == "CP":
      num_peaks = int(line[6])
      det = int(line[2])

      for peak in range(num_peaks):
        line_adc = float(line[7+3*peak])
        line_energy = float(line[8+3*peak])
        fwhm_adc = 0.0
        fwhm_energy = float(line[9+3*peak])
        fwhm_percent = 100*fwhm_energy/line_energy
        if line_energy == 1836.06: # optionally select one energy to histogram 
        #if line[4]=='p':         # optionally select one side to histogram
            #if det==11:
          fwhm_energy_ls.append(fwhm_energy)
          fwhm_percent_ls.append(fwhm_percent)

  f.close()
  return fwhm_energy_ls,fwhm_percent_ls

# Calculate mean, standard deviation
def mean_stddev(list):
  mean = sum(list) / len(list)  
  variance = sum([((x - mean) ** 2) for x in list]) / len(list) 
  stddev = variance ** 0.5  
  error = stddev / len(list)**0.5
  return mean,stddev,error

# Main program make histograms
if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(
      description="Create histogram of all FWHMs [keV,%] from a melinator .ecal file.", 
      epilog="Plot 0: FWHMs [keV]. Plot 1: FWHMs [%], where % = 100*FWHM[kev]/line energy[kev]")  
  parser.add_argument('ecal', default=None,help='.ecal file with FWHMs of interest.')
  parser.add_argument('--ecal2',default=None,help='Optionally provide a second .ecal file')
  args = parser.parse_args()

  ecal_file = args.ecal
  fwhm_energy_ls,fwhm_percent_ls = parse_peaks(ecal_file)

  fwhm_energy_mean = mean_stddev(fwhm_energy_ls)[0];fwhm_energy_stddev = mean_stddev(fwhm_energy_ls)[1];fwhm_energy_error = mean_stddev(fwhm_energy_ls)[2]
  fwhm_percent_mean = mean_stddev(fwhm_percent_ls)[0];fwhm_percent_stddev = mean_stddev(fwhm_percent_ls)[1];fwhm_percent_error = mean_stddev(fwhm_percent_ls)[2] 

  # Histogram of FWHMs [keV]
  if args.ecal2 is not None:
    fwhm_energy_ls2 = parse_peaks(args.ecal2)[0]
    fwhm_energy_mean2 = mean_stddev(fwhm_energy_ls2)[0];fwhm_energy_stddev2 = mean_stddev(fwhm_energy_ls2)[1]
    fwhm_energy_error2 = mean_stddev(fwhm_energy_ls2)[2]
    plt.hist(fwhm_energy_ls,bins=200,histtype='step',label=ecal_file,density=True)
    plt.hist(fwhm_energy_ls2,bins=200,histtype='step',label=args.ecal2,density=True)
    plt.xlabel('FWHM [keV]');plt.title('DC strips, 661.657 keV: '+ecal_file+' and '+args.ecal2)
    #plt.xlim(0,max(np.max(fwhm_energy_ls),np.max(fwhm_energy_ls2))+1)
    #plt.xlim(0,1)
    plt.legend(loc='best')
    plt.text(5,0.5,ecal_file+': mean={} +/- {} keV'.format(round(fwhm_energy_mean,3),round(fwhm_energy_error,3))) 
    plt.text(5,0.4,args.ecal2+': mean={} +/- {} keV'.format(round(fwhm_energy_mean2,3),round(fwhm_energy_error2,3)))   
    plt.show()

  else:
    plt.hist(fwhm_energy_ls,bins=50,histtype='step',density=True)
    plt.xlabel('FWHM [keV]');plt.title(ecal_file)
    plt.xlim(0,15)
    plt.text(5,4000,ecal_file+': mean={} +/- {} keV'.format(round(fwhm_energy_mean,3),round(fwhm_energy_error,3)))    
    plt.show()


  # Histogram of FWHMs [%]
  if args.ecal2 is not None:
    fwhm_percent_ls2 = parse_peaks(args.ecal2)[1]
    fwhm_percent_mean2 = mean_stddev(fwhm_percent_ls2)[0];fwhm_percent_stddev2 = mean_stddev(fwhm_percent_ls2)[1]
    fwhm_percent_error2 = mean_stddev(fwhm_percent_ls2)[2]
    plt.hist(fwhm_percent_ls,bins=200,histtype='step',label=ecal_file,density=True)
    plt.hist(fwhm_percent_ls2,bins=200,histtype='step',label=args.ecal2,density=True)
    plt.xlabel('FWHM [%]');plt.title('DC strips, 661.657 keV: '+ecal_file+' and '+args.ecal2)
    #plt.xlim(0,max(np.max(fwhm_percent_ls),np.max(fwhm_percent_ls2))+0.5)
    #plt.xlim(0.2,0.8)
    plt.legend(loc='best')
    plt.text(0.45,5,ecal_file+': mean={} +/- {} %'.format(round(fwhm_percent_mean,3),round(fwhm_percent_error,3))) 
    plt.text(0.45,4,args.ecal2+': mean={} +/- {} %'.format(round(fwhm_percent_mean2,3),round(fwhm_percent_error2,3)))   
    plt.show()
    print("2020:",round(fwhm_percent_mean,3),round(fwhm_percent_error,3))
    print("2016:",round(fwhm_percent_mean2,3),round(fwhm_percent_error2,3))
  else:
    plt.hist(fwhm_percent_ls,bins=50,histtype='step',density=True)
    plt.xlabel('FWHM [%]');plt.title(ecal_file)
    #plt.xlim(0,max(np.max(fwhm_energy_ls),np.max(fwhm_energy_ls2))+1)
    #plt.xlim(0,8)
    plt.text(1.6,0.5,ecal_file+': mean={} +/- {} %'.format(round(fwhm_percent_mean,3),round(fwhm_percent_error,3)))    
    plt.show()

 
