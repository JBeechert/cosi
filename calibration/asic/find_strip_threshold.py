# 20/09/25:
# Used Cs-137 .hdf.roa files corresponding to the ASIC 0 1us and 2us peaking time and 18.4mV/fC gain settings.
# For each, use find_strip_threshold.py to generate histograms of each channel and visually identify the threshold in ADC.
# Each channel's histogram is saved in the "_plots" folder corresponding to its peaking time and gain settings.
# Enter each channel's threshold in ADC into the code and convert to keV, using a specified .ecal file.
# The threshold results are in google drive folder ASIC Thresholds/NRL1v2B ASIC 0 Thresholds

import numpy as np
import matplotlib.pyplot as plt

# Parse .ecal, return dictionary of {strip:[polynomial coefficients]}
def parse_ecal(ecal_fname):
  f = open(ecal_fname,"r")
  lines = f.readlines()
  dict = {}

  for line in lines:
    line = line.split(' ')
    if line[0] == "CM":
      strip = str(line[2]) + str(line[4]) + str(line[3])
      params = line[6:-1]
      params = [float(i) for i in params]
      if strip in dict:
        dict[strip].append(params)
      else:
        dict[strip] = params

  f.close()
  return dict

# Parse .roa, return dictionary of {strip:[all ADC values]}
def parse_roa(roa_fname):
  f = open(roa_fname,"r")
  lines = f.readlines()  
  dict = {}

  for line in lines:
    line = line.split(' ')
    if line[0] == "UH":
      strip = str(line[1]) + str(line[2]) + line[3]
      ADC = float(line[4])
      
      if strip in dict:
        dict[strip].append(ADC)
      else:
        dict[strip] = [ADC]
  
  f.close()
  return dict


# Convert from ADC to keV, return dictionary of {strip: threshold in keV} given
#  ecal_dict = {strip:[polynomial coefficients} and threshold_dict {strip:threshold identified from histogram}
def polynomial(p):
    return lambda x: sum(a*x**i for i, a in enumerate(p))

def ADC_to_keV(ecal_dict,threshold_dict):
  keV_dict = {}
  for key,value in ecal_dict.items():
    for key2,value2 in threshold_dict.items():
      if key==key2:
        value = [float(i) for i in value]
        poly = polynomial(value)
        keV = poly(float(value2))
        keV_dict[key] = keV

  return keV_dict


# Make the histograms
#roa_fname = "asic_gse_200922234020.hdf.roa"
roa_fname = "asic_gse_200922205740.hdf.roa"

#roa_dict = parse_roa(roa_fname)

#for key, value in roa_dict.items():
#  plt.hist(value,bins=np.arange(2**14),histtype='step')
#  plt.title("ASIC0 1us 18.4mV/fC Cs-137 \n" + roa_fname + " strip: "+ key)
#  plt.xlabel("ADC");plt.ylabel("counts")
#  plt.show()

# Enter the threshold [ADC] visually identified from each histogram above
threshold_dict = {'0p13':1800.95,'0p15':1853.69,'0p17':1951.7,'0p23':1843.29,'0p25':1828.14,'0p27':1756.13,'0p29':1913.6,'0p31':1861.18,'0p33':1800.96,'0p35':1889.44,'0p37':1883.5,'0p39':1793.33,'0p41':1801.9,'0p43':1749.38,'0p45':1851.4,'0p47':1818.61,'0p49':1753.53, \
                  '0p51':1876.65,'0p53':1839.85,'0p55':1877.5,'0p57':1850.91,'0p59':1857.6,'0p61':1818.89,'0p63':1825.23,'0p65':1805.63,'0p67':1793.04,'0p69':1831.42,'0p71':1777.69,'0p73':1814.27,'0p75':1840.56}

# Do the conversion to keV given the thresholds in ADC
#ecal_dict = parse_ecal("am_bar_cs_1us_18.4mVfC_200923_NRL1v2B_ASIC0_poly3.ecal")
ecal_dict = parse_ecal("am_bar_cs_2us_18.4mVfC_200923_NRL1v2B_ASIC0_poly3.ecal")

keV_dict = ADC_to_keV(ecal_dict,threshold_dict)
print(keV_dict)
