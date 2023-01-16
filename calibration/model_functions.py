import numpy as np
import matplotlib.pyplot as plt
import emcee
import time

print(f'Need emcee version 3\nCurrent version: {emcee.__version__}')


def mcmc(init, ndim, fit_func, ADC, energy, nwalkers=10, iters=2000):
	# mcmc fit

	# intial guess for fit parameters
	init = init

	# vary the initial guess 
	init_var = init*1e-4

	# dimensionality of the problem and defining multiple starting points for fit
	ndim, nwalkers = ndim, 10
	pos = [init + np.random.randn(ndim)*init_var for i in range(nwalkers)]

	# taking time
	start = time.time()

	# setting up the sampler, i.e. the thing that holds the information about the fit
	# and performs the iterative steps; otherwise similar to calling statement from above
	sampler = emcee.EnsembleSampler(nwalkers,
									ndim,
									fit_func,
									args = (ADC,
											energy))

	# do the actual Monte Carlo sampling for 2000 iterations
	iters = 2000
	_ = sampler.run_mcmc(pos, iters, progress=True)

	# taking time again
	end = time.time()

	# extract samples
	samples = sampler.get_chain()
	samplesf = sampler.flatchain
	
	ttime = end - start
	print("Processing took {0:.1f} seconds".format(ttime))

	return sampler, samples, samplesf    


def parameter_table(samples, ndim, burnin):
	# output here
	# see explanation below for what is actually happening

	print('Results:')

	spec_params = np.zeros((ndim, 7))

	# formatting the table
	row_format ='{:>10}' * 8

	# first table row
	print(row_format.format(*['Parameter', 'mean', 'std', '0.15', '15.85', '50.00', '84.15', '99.85']))

	for i in range(ndim):
		mean_val   = np.mean(samples[burnin:, :, i])
		std_val    = np.std(samples[burnin:, :, i])
		median_val = np.median(samples[burnin:, :, i])
		ub1_val    = np.percentile(samples[burnin:, :, i], 50+68.3/2)
		lb1_val    = np.percentile(samples[burnin:, :, i], 50-68.3/2)
		ub3_val    = np.percentile(samples[burnin:, :, i], 50+99.73/2)
		lb3_val    = np.percentile(samples[burnin:, :, i], 50-99.73/2)
		spec_params[i,:] = [mean_val, std_val, lb3_val, lb1_val, median_val, ub1_val, ub3_val]

		print(row_format.format(str(i)+':',
								str('%1.2e' % mean_val),
								str('%1.2e' % std_val),
								str('%1.2e' % lb3_val),
								str('%1.2e' % lb1_val),
								str('%1.2e' % median_val),
								str('%1.2e' % ub1_val),
								str('%1.2e' % ub3_val)))


def evaluate_model(energy, nwalkers, iters, samplesf, ndim):

	# evaluate models according to sampling

	# we calculate only the last 100 samples from each walker,
	# because it will just take forever otherwise
	n_use = 100
	n_plot_samples = nwalkers*n_use

	y_models = np.zeros((len(energy), n_plot_samples))

	last_x_samples = iters - n_use

	for i in range(nwalkers*last_x_samples, nwalkers*last_x_samples + n_plot_samples):

		if ndim == 2:
			y_models[:, i-nwalkers*last_x_samples] = samplesf[i,0] + samplesf[i,1]*energy

		elif ndim == 3:
			y_models[:, i-nwalkers*last_x_samples] = samplesf[i,0] + samplesf[i,1]*energy + samplesf[i,2]*energy**2

		elif ndim == 4:
			y_models[:, i-nwalkers*last_x_samples] = samplesf[i,0] + samplesf[i,1]*energy + samplesf[i,2]*energy**2 + samplesf[i,3]*energy**3

		elif ndim == 5:
			y_models[:, i-nwalkers*last_x_samples] = samplesf[i,0] + samplesf[i,1]*energy + samplesf[i,2]*energy**2 + samplesf[i,3]*energy**3 + samplesf[i,4]*energy**4

		else:
			print("Function only evaluates linear, quadratic, cubic, and quartic models. Check ndim or update the function.")
		
	return y_models


def plot_model(y_models, ADC, energy, model_label, color):

	# data
	plt.plot(ADC, energy, 'k.', label='Data')

	# looping over some self-defined confidence levels, where
	# 0 is the median (one central line), 68.3 is 1sigma, 95.4 is 2sigma
	for level in [0, 68.3, 95.4]:
		# total model
		plt.fill_between(ADC,
						 np.percentile(y_models, 50-level/2, axis=1),
						 np.percentile(y_models, 50+level/2, axis=1),
						 color=color, alpha=0.3, label=model_label)

	# and a useful legend
	plt.legend(['Data', model_label])
	plt.xlabel('ADC')
	plt.ylabel('Energy [keV]');


def plot_residual_per_deviation(y_models, ADC, energy, model_label, color):

	# median model and stddev
	median_model = np.percentile(y_models, 50, axis=1)
	median_model_std = np.sqrt(median_model)

	# data
	plt.step(ADC, (energy-np.percentile(y_models, 50, axis=1))/median_model_std, 'k', where='mid')

	# looping over some self-defined confidence levels, where
	# 0 is the median (one central line), 68.3 is 1sigma, 95.4 is 2sigma
	for level in [0, 68.3, 95.4]:
		# total model
		plt.fill_between(ADC,
						 (energy-np.percentile(y_models, 50-level/2, axis=1))/median_model_std,
						 (energy-np.percentile(y_models, 50+level/2, axis=1))/median_model_std,
						 color=color, 
						 alpha=0.2, 
						 label=model_label, 
						 step='mid')

	# and a useful legend
	plt.legend(['Median residual', model_label], fontsize=14)
	plt.xlabel('ADC', fontsize=14)
	plt.ylabel(r'Residuals [$\sigma$]', fontsize=14)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.axhline(y=0, color='k', linestyle='--');

