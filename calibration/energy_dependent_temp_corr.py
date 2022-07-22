# Attempt an energy-dependent temperature correction
# Fit slope vs. offset and offset vs. energy plots with a piecewise function.

import numpy as np
import matplotlib.pyplot as plt

# 2020 method: used Carolyn's Cs-137 and Am-241 data
_2016_59keV_2020method_slope_param = -4.0E-4
_2016_59keV_2020method_slope_param_err = 3.4E-5
_2016_59keV_2020method_offset_param = 1.010
_2016_59keV_2020method_offset_param_err = 6.911E-4

_2016_662keV_2020method_slope_param = -3.5E-4
_2016_662keV_2020method_slope_param_err = 2.1E-5
_2016_662keV_2020method_offset_param = 1.011
_2016_662keV_2020method_offset_param_err = 4.322E-4

# 2020 method: checked all peaks except 59 keV
energies = [80.997,122.061,136.474,276.398,302.9,356.0,383.85,511.0,661.7,898.0,1173.2,1274.5,1332.5,1836.0]
_2020_slope_params = [-4.6E-4,-2.08E-4,-1.53E-4,-3.5E-4,-3.3E-4,-3.1E-4,-3.0E-4,-2.01E-4,-1.01E-4,-1.41E-4,-1.04E-4,-5.94E-5,-1.02E-4,-1.36E-4]
_2020_slope_params_err = [8.3E-5,1.92E-4,1.90E-4,3.0E-5,3.0E-5,2.8E-5,2.9E-5,2.44E-5,1.93E-5,9.97E-6,8.88E-6,1.2E-5,1.03E-5,1.33E-5]
_2020_offset_params = [1.011,1.004,1.003,1.0097,1.009,1.009,1.008,1.005,1.002,1.004,1.003,1.001,1.003,1.004]
_2020_offset_params_err = [2.46E-3,6.07E-3,6.02E-3,8.715E-4,8.788E-4,8.097E-4,8.381E-4,7.857E-4,5.763E-4,2.875E-4,2.74E-4,3.892E-4,3.18E-4,3.832E-4]

# 2016 method: took the average of all slope and offset parameters in /nuclearizer/resource/calibration/COSI16/Wanaka/TempCalibration_041118.txt
# See 2016_temp_mean.py
#  Since this was a strip-by-strip calibration, it's an average over all strips
_2016_avg_slope = -3.4E-4
_2016_avg_slope_stddev = 4.9E-5
_2016_avg_offset = 1.012
_2016_avg_offset_stddev = 0.0012


# Plot slope parameters
plt.errorbar(energies,_2020_slope_params,yerr=_2020_slope_params_err,color="k",fmt=".",label="2020 data and method")
plt.axhline(y=np.mean(_2020_slope_params),color="k",linestyle="--",label="2020 average across all peaks (shading = stddev)")
plt.fill_between(np.arange(0,2000,1),np.mean(_2020_slope_params) - np.std(_2020_slope_params),np.mean(_2020_slope_params) + np.std(_2020_slope_params),color="k",alpha=0.2)


plt.errorbar(59.541,_2016_59keV_2020method_slope_param,yerr=_2016_59keV_2020method_slope_param_err,color="r",fmt=".",label="2016 Am-241 data, 2020 method")
plt.errorbar(661.657,_2016_662keV_2020method_slope_param,yerr=_2016_662keV_2020method_slope_param_err,color="r",fmt=".",label="2016 Cs-137 data, 2020 method")
plt.axhline(y=_2016_avg_slope,color="r",linestyle="--",label="2016 average across all strips, 2016 method (shading = stddev)")
plt.fill_between(np.arange(0,2000,1),_2016_avg_slope - _2016_avg_slope_stddev,_2016_avg_slope + _2016_avg_slope_stddev,color="r",alpha=0.2)

plt.xlabel("Photopeak energy [keV]")
plt.ylabel("Slope parameter [1/C]")
plt.xticks([59.541]+energies,rotation=75)


# Fit piecewise linear fit (excluding Co-57 outliers)
import pwlf
energies_fit = [80.997,276.398,302.9,356.0,383.85,511.0,661.7,898.0,1173.2,1274.5,1332.5,1836.0]
_2020_slope_params_fit = [-4.6E-4,-3.5E-4,-3.3E-4,-3.1E-4,-3.0E-4,-2.01E-4,-1.01E-4,-1.41E-4,-1.04E-4,-5.94E-5,-1.02E-4,-1.36E-4]
my_pwlf = pwlf.PiecewiseLinFit(energies_fit,_2020_slope_params_fit)
res = my_pwlf.fit(2)
print(res)
xHat = np.linspace(min(energies_fit),max(energies_fit),num=10000)
yHat = my_pwlf.predict(xHat)
plt.plot(xHat,yHat,"-",label="pwise fit")
plt.legend(loc="best")
plt.xlim(0,2000)
plt.show()

from sympy import Symbol
from sympy.utilities import lambdify
x = Symbol('x')


def get_symbolic_eqn(pwlf_, segment_number):
    if pwlf_.degree < 1:
        raise ValueError('Degree must be at least 1')
    if segment_number < 1 or segment_number > pwlf_.n_segments:
        raise ValueError('segment_number not possible')
    # assemble degree = 1 first
    for line in range(segment_number):
        if line == 0:
            my_eqn = pwlf_.beta[0] + (pwlf_.beta[1])*(x-pwlf_.fit_breaks[0])
        else:
            my_eqn += (pwlf_.beta[line+1])*(x-pwlf_.fit_breaks[line])
    # assemble all other degrees
    if pwlf_.degree > 1:
        for k in range(2, pwlf_.degree + 1):
            for line in range(segment_number):
                beta_index = pwlf_.n_segments*(k-1) + line + 1
                my_eqn += (pwlf_.beta[beta_index])*(x-pwlf_.fit_breaks[line])**k
    return my_eqn.simplify()


eqn_list = []
f_list = []
for i in range(my_pwlf.n_segments):
    eqn_list.append(get_symbolic_eqn(my_pwlf, i + 1))
    print('Equation number: ', i + 1)
    print(eqn_list[-1])
    f_list.append(lambdify(x, eqn_list[-1]))




# Plot offset parameters
plt.errorbar(energies,_2020_offset_params,yerr=_2020_offset_params_err,color="k",fmt=".",label="2020 data and method")
plt.axhline(y=np.mean(_2020_offset_params),color="k",linestyle="--",label="2020 average across all peaks (shading = stddev)")
plt.fill_between(np.arange(0,2000,1),np.mean(_2020_offset_params) - np.std(_2020_offset_params),np.mean(_2020_offset_params) + np.std(_2020_offset_params),color="k",alpha=0.2)


plt.errorbar(59.541,_2016_59keV_2020method_offset_param,yerr=_2016_59keV_2020method_offset_param_err,color="r",fmt=".",label="2016 Am-241 data, 2020 method")
plt.errorbar(661.657,_2016_662keV_2020method_offset_param,yerr=_2016_662keV_2020method_offset_param_err,color="r",fmt=".",label="2016 Cs-137 data, 2020 method")
plt.axhline(y=_2016_avg_offset,color="r",linestyle="--",label="2016 average across all strips, 2016 method (shading = stddev)")
plt.fill_between(np.arange(0,2000,1),_2016_avg_offset - _2016_avg_offset_stddev,_2016_avg_offset + _2016_avg_offset_stddev,color="r",alpha=0.2)

plt.xlabel("Photopeak energy [keV]")
plt.ylabel("Offset parameter")
plt.xticks([59.541]+energies,rotation=75)


# Fit piecewise linear fit (excluding Co-57 outliers)
_2020_offset_params_fit = [1.011,1.0097,1.009,1.009,1.008,1.005,1.002,1.004,1.003,1.001,1.003,1.004]
my_pwlf_2 = pwlf.PiecewiseLinFit(energies_fit,_2020_offset_params_fit)
res_2 = my_pwlf_2.fit(2)
print(res_2)
xHat_2 = np.linspace(min(energies_fit),max(energies_fit),num=10000)
yHat_2 = my_pwlf_2.predict(xHat_2)
plt.plot(xHat_2,yHat_2,"-",label="pwise fit")
plt.legend(loc="best")
plt.xlim(0,2000)
plt.show()

eqn_list_2 = []
f_list_2 = []
for i in range(my_pwlf_2.n_segments):
    eqn_list_2.append(get_symbolic_eqn(my_pwlf_2, i + 1))
    print('Equation number: ', i + 1)
    print(eqn_list_2[-1])
    f_list_2.append(lambdify(x, eqn_list_2[-1]))
