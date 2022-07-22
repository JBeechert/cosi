# Plot results from Al-26 simulations over the COMPTEL map.
# These results (alpha, beta, significance, flux, etc.) were saved to a text file in the jupyter notebook which ran
#  this analysis. That text file is read in here and plotted.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

for num in [1,2,5,10,20,30,40,50]:
	path = "/volumes/eos/users/jacqueline/analysis/Al26/COMPTEL/choose_random_COMPTEL_runs/"
	inc_num,alpha,alpha_stds,beta,beta_stds,significance,Y,X,flux,flux_stds = np.genfromtxt(path+"{}COMPTEL/results_{}COMPTELruns_1flightactivation_1flightLinginc2.txt".format(num,num),unpack=True)

	num_itrs = len(alpha)

	significance_mean = np.mean(significance)
	significance_std = np.std(significance)
	flux_mean = np.mean(flux)
	flux_std = np.std(flux)

	# Plot flux vs. significance
#	print(num)
#	print(flux,significance,"\n")
	plt.scatter(flux,significance,label=r"$n$ = {}".format(num))

plt.title("$n$ flights COMPTEL, 1 flight activation, 1 day Ling \n {} entries, each a random selection of $n$ COMPTEL simulations \n 1 Ling simulation: inc2".format(num_itrs),fontsize=12)
plt.xlabel(r"Simulated flux measured in ($\ell$ $\leq$ 65$^{\circ}$, $b$ $\leq$ 45$^{\circ}$) [10$^{-4}$ ph cm$^{-2}$ s$^{-1}$]",fontsize=12)
plt.ylabel("Significance [$\sigma$ above background]",fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.show()
