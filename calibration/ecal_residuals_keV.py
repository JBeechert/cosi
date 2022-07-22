import numpy as np
import matplotlib.pyplot as plt

def d_energyADC(file):
    file_object = open(str(file),"r")
    lines = file_object.readlines()
    d = {}
    
    for line in lines:
        if line.startswith("CP"):
            if int(line.split()[6]) > 0:
                det = line.split()[2]; strip = line.split()[3]; side = line.split()[4] 
                key = det+" "+strip+" "+side
                ADC = np.array([float(i) for i in line.split()[7::3]])
                energy = np.array([float(i) for i in line.split()[8::3]])
                FWHM = np.array([float(i) for i in line.split()[9::3]])

                d[key] = [ADC,energy,FWHM]
    return d

def d_model(file):
    file_object = open(str(file),"r")
    lines = file_object.readlines()
    model = []; d = {}
    
    for line in lines:
        if line.startswith("CP"):
            det = line.split()[2]; strip = line.split()[3]; side = line.split()[4] 
            key = det+" "+strip+" "+side
        
        if line.startswith("CM"):
            if key in line:
                model = [float(x) for x in line.split()[6:]]
                max_num_coeffs = 5
                if len(model) < max_num_coeffs:
                    model.extend([0]*(max_num_coeffs - len(model)))

                if len(model) > max_num_coeffs:
                    print("You're using a calibration model with more than 5 coefficients. Code isn't updated for that yet. Sorry.")
                d[key] = model
    return d

def ADC_to_energy(x,model,modeltype):
    energy = 0
    if modeltype == "poly":
        energy = [model[0] + model[1]*i + model[2]*i**2 + model[3]*i**3 + model[4]*i**4 for i in x]
    if modeltype == "poly1inv1":
        energy = [model[0] + model[1]*i + (1. / (model[2] + model[3]*i)) for i in x]
    if modeltype == "poly1inv1zero":
        energy = [model[0] + model[1]*i + (1. / (model[2]*i)) for i in x]
    if modeltype == "poly2inv1":
        energy = [model[0] + model[1]*i + model[2]*i**2 + (1. / (model[3] + model[4]*i)) for i in x]
    if modeltype == "poly2inv1zero":
        energy = [model[0] + model[1]*i + model[2]*i**2 + (1. / (model[3]*i)) for i in x]
    if modeltype == "poly1ln1":
        energy = [model[0] + model[1]*i + model[2]*np.log(i) for i in x]
    if modeltype == "poly2ln1":
        energy = [model[0] + model[1]*i + model[2]*i**2 + model[3]*np.log(i) for i in x]
    if modeltype == "poly1exp1":
        energy = [model[0] + model[1]*i + model[2]*np.exp((-i/model[3])) for i in x]
    if modeltype == "poly1exp2":
        energy = [model[0] + model[1]*i + model[2]*np.exp((-i/model[3])**2) for i in x]
    if modeltype == "poly1exp3":
        energy = [model[0] + model[1]*i + model[2]*np.exp((-i/model[3])**3) for i in x]

    return energy

def d_calib_result(ecal,modeltype,photopeak=None):
    energyADC_dict = d_energyADC(ecal)
    model_dict = d_model(ecal)
    all_residuals = []

    d = {}

    # Each key is a strip
    for key in energyADC_dict:
        
        # List of identified, true energies on the strip
        true_energy = energyADC_dict[key][1]
        
        # ADCs at which these true energies were located 
        ADC = energyADC_dict[key][0]
        
        # Model's conversion of these ADCs to energy
        model_energy = ADC_to_energy(ADC,model_dict[key],modeltype)

        # Calculate residuals from all photopeaks
        if photopeak == None:
            residual = model_energy - true_energy
            d[key] = residual

            all_residuals.extend(residual)
 
        # Or, optionally, select residual at a certain energy
        #  Yes, this code is getting sloppy, leave me alone
        if photopeak != None:
            for energy in true_energy:
               if energy == photopeak:
                  idx = np.where(photopeak==true_energy)[0][0]
                  residual = model_energy[idx] - energy
                  d[key] = residual

                  all_residuals.append(residual)

    return d,all_residuals


if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description="Check predicted energies from a melinator .ecal file.")  
  parser.add_argument('ecal', default=None,help='.ecal file of interest.')
  args = parser.parse_args()

  ecal_file = args.ecal

  # Photopeak options are :
    # None for all photopeaks
    # Or, choose 59.541 80.997 122.061 136.474 276.398 302.853 356.017 383.851 510.99 661.657 898.042 1173.24 1274.58 1332.5 1836.06
  peak = None
  d_residuals,residuals = d_calib_result(ecal_file,"poly",photopeak=None)

  # Remove outliers (residuals > y keV)
  all_residuals = []
  for i in residuals:
  	if abs(i) < 2:
  		all_residuals.append(i)

  # Histogram all residuals
  mean_all_residuals = np.mean(all_residuals)
  stddev_all_residuals = np.std(all_residuals)
  plt.hist(all_residuals,bins=100)
  plt.title("All residuals, all strips"+"\n"+ecal_file+"\n"+"Mean of all residuals: {} keV".format(round(mean_all_residuals,3))+", std.dev. of all residuals: {} keV".format(round(stddev_all_residuals,3)))
  plt.xlabel("Model energy - true energy [keV]")
  plt.show()


  # # Histogram the mean residual of all peaks on each strip (879 mean residuals)
  # mean_residuals = []

  # for key in d_residuals.keys():
  #    # If you want all strips
  #    mean_residuals.append(np.mean(d_residuals[key]))
    
  #    # If you want individual detector
  #    # if key[0:2] == "7 ":  
  #    #    mean_residuals.append(np.mean(d_residuals[key]))

  # # mean of the mean residual on each strip
  # mean_residual = np.mean(mean_residuals)

  # plt.hist(mean_residuals,bins=100)
  # plt.title("Mean of the mean residual on each strip: {}".format(round(mean_residual,3))+"keV"+"\n"+ecal_file)
  # #plt.title("Mean residual on each strip, {} keV: {}".format(peak,round(mean_residual,3))+"keV"+"\n"+ecal_file)
  # plt.xlabel("Mean residual per strip (879 values here): Model energy - true energy [keV]")
  # plt.show()

