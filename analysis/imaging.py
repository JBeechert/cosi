import numpy as np
import matplotlib.pyplot as plt
from tqdm.autonotebook import tqdm

from COSIpy import *
from COSIpy_tools import zenaziGrid, one_func, cashstat
import response

import pickle
import pywt
import pystan

import warnings
warnings.filterwarnings('ignore')


class Imaging:
    """
    A class to conduct imaging analyses with Richardson-Lucy (RL) deconvolution.

    Attributes
    ----------
    

    Methods
    -------

    """
    # conversion from degrees to radians
    deg2rad = np.pi/180    


    def __init__(self, pixel_size, delta_map_scaling=2000, init_map_flux=1.0, fitted_bg=np.array([0.2])):
        """
        Constructs the necessary objects for the imaging class.
        Sets up the pixels in the sky and initializes the variables and initial image for RL deconvolution.

        Parameters
        ----------
        pixel_size : float, required. One-dimensional size of pixels in Galactic longitude/latitude [deg]. 
                                      Typically 5. or 6. degrees.
        delta_map_scaling : int, optional. Scaling factor for the 'delta map', allowing it to be stronger 
                                           than the original RL algorithm suggests (Knoedlseder+1997). Default = 2000
        init_map_flux : float, optional. Flux of the initial, isotropic map. Default = 1.0. 
                                         Consider smaller values (0.1, 0.01, etc.) for weak or BG-dominated sources.
        fitted_bg : np.ndarray, optional. Inital guess of bg/total ratio. Default = np.array([0.2]).
                                          For data without background, i.e. source only, the suggested value is np.array([1e-6])
        """
        print(f"\nInitializing the imaging class with \n\
            \t pixel size = {pixel_size} deg \n\
            \t delta map scaling parameter = {delta_map_scaling} \n\
            \t initial map flux = {init_map_flux} \n\
            \t inital guess of bg/total ratio = {fitted_bg}")

        ############## Set up the sky 
        self.pixel_size = pixel_size

        # number of pixels in l and b
        self.n_l = int(360/self.pixel_size)
        self.n_b = int(180/self.pixel_size)

        # l and b pixel edges
        self.l_arrg = np.linspace(-180, 180, self.n_l+1)
        self.b_arrg = np.linspace(-90, 90, self.n_b+1)

        # making a grid
        self.L_ARRg, self.B_ARRg = np.meshgrid(self.l_arrg, self.b_arrg)

        # choosing the centre points as representative
        self.l_arr = self.l_arrg[0:-1] + self.pixel_size/2
        self.b_arr = self.b_arrg[0:-1] + self.pixel_size/2
        self.L_ARR, self.B_ARR = np.meshgrid(self.l_arr, self.b_arr)

        # define solid angle for each pixel for normalisations later
        self.domega = (self.pixel_size*self.deg2rad)*(np.sin(np.deg2rad(self.B_ARR + self.pixel_size/2)) 
            - np.sin(np.deg2rad(self.B_ARR - self.pixel_size/2)))

        # calculate the zeniths and azimuths on that grid for all times
        self.zensgrid, self.azisgrid = zenaziGrid(pointing1.ypoins[:,0], pointing1.ypoins[:,1],
                               pointing1.xpoins[:,0], pointing1.xpoins[:,1],
                               pointing1.zpoins[:,0], pointing1.zpoins[:,1],
                               self.L_ARR.ravel(), self.B_ARR.ravel())

        # reshape for next routines ...
        self.zenith = self.zensgrid.reshape(self.n_b, self.n_l, len(pointing1.xpoins))
        self.azimuth = self.azisgrid.reshape(self.n_b, self.n_l, len(pointing1.xpoins))


        ########### Initialize variables for RL

        # Scaling factor for the 'Delta map', allowing it to be stronger than the original RL algorithm suggests (Knoedlseder+1997)
        self.delta_map_scaling = delta_map_scaling

        # initial guess for bg parameter (e.g. ratio of Ling BG level to total spectrum)
        self.fitted_bg = fitted_bg

        ###### Define initial image for RL (isotropic flat, small amplitude A0)
        A0 = init_map_flux
        shape = np.ones(self.L_ARR.shape)
        norm = np.sum(shape*(self.pixel_size*np.pi/180)*(np.sin(np.deg2rad(self.B_ARR + self.pixel_size/2)) - np.sin(np.deg2rad(self.B_ARR - self.pixel_size/2))))
        self.current_image = A0*shape/norm


    def next_iteration(self, map_init, expectation_init, sky_response_scaled, expo_map, background1, analysis1, show_its=False, iteration=None):
        """
        Apply one iteration of the RL deconvolution algorithm to the provided image.

        Parameters:
             map_init (np.ndarray, required): The sky response from COSIpy's response class
             expectation_init (np.ndarray, required): Initial expectation for the image (initial isotropic map convolved with the response)
             sky_response_scaled (np.ndarray, required): Scaled sky response. Shape: (# energy bins, n_b, n_l, FISBEL)
             expo_map (np.ndarray, required): Exposure map. Shape: (n_b, n_l)
             background1 (COSIpy.BG object, required): Contains the background information returned by the COSIpy BG class
             analysis1 (COSIpy.COSIpy object, required): Contains the data/simulation to be analyzed
             show_its (bool, optional): Display plots of the image and delta image of the iteration. Default: False
             iteration (int, optional): Number of this iteration. Default: None. If show_its is True, iteration should not be None.

        Returns:
             [map_new, expectation_new, map_likelihood, intermediate_lp, acc_par, bg_pars] (list): 
                         Image, expectation image, map likelihood, fit likelihood, acceleration parameter, 
                        and fitted background parameter of the iterated iteration
        """
        # as zeroth iteration, copy initial map to become the 'old map' (see below)
        map_old = map_init #self.current_image

        # set old expectation (in data space bins) to new expectation (convolved image)
        expectation_old = expectation_init


        # cf. Knoedlseder+1997 what the values denominator etc are
        # this is the response R summed over the CDS and the time bins
        denominator = expo_map

        ###### now that we have an expectation for the image, look at the BG
        # define bad exposure
        bad_expo = np.where(expo_map/self.domega < 0)

        map_old[bad_expo[0], bad_expo[1]] = 0

        # check for each pixel to be finite
        map_old[np.where(np.isnan(map_old) == True)] = 0

        # make new background for the next iteration
        bg_cuts, idx_arr, Ncuts = background1.bg_cuts, background1.idx_arr, background1.Ncuts

        d2h = background1.bg_model.shape[0]

        background_model = background1.bg_model_reduced[ebin]

        dataset = analysis1.dataset.binned_data[:, ebin, :, :].reshape(d2h, 30*1145)[:, nonzero_idx]

        # temporary background model
        tmp_model_bg = np.zeros((d2h, background_model.shape[1]))#background1.bg_model_reduced[ebin].shape[1]))

        for g in range(d2h):
            tmp_model_bg[g, :] = background_model[g, :]*self.fitted_bg[idx_arr-1][g]
                
        # expectation (in data space) is the image (expectation_old) plus the background (tmp_model_bg)
        expectation_tot_old = expectation_old + tmp_model_bg 

        # calculate likelihood of current total expectation
        map_likelihood = cashstat(dataset.ravel(), expectation_tot_old.ravel())
        
        # calculate numerator of RL algorithm
        numerator = 0
        print('Calculating Delta image')
        for i in tqdm(range(d2h)):
            for j in range(dataset.shape[1]):
                numerator += (dataset[i, j]/expectation_tot_old[i, j]-1)*sky_response_scaled[i, :, :, j]


        # calculate delta map (denominator scaled by fourth root to avoid exposure edge effects. Can change this exponent to 0, 0.5, etc.)
        delta_map_tot_old = (numerator/denominator)*map_old*(denominator)**0.25
        ######

        # check again for finite values and zero our bad exposure regions
        nan_idx = np.where(np.isnan(delta_map_tot_old) == 1)
        delta_map_tot_old[nan_idx[0], nan_idx[1]] = 0
        delta_map_tot_old[bad_expo[0], bad_expo[1]] = 0


        # plot the image and delta image
        if show_its:
            plt.figure(figsize=(16,6))
            plt.subplot(121)
            plt.pcolormesh(self.L_ARRg, self.B_ARRg, np.roll(map_old, axis=1, shift=0))
            plt.colorbar()
            plt.title(f'expectation old, iteration {iteration}')
            plt.subplot(122)
            plt.pcolormesh(self.L_ARRg, self.B_ARRg, np.roll(delta_map_tot_old, axis=1, shift=0))
            plt.colorbar()
            plt.title(f'delta map, iteration {iteration}')
            plt.show()

        # convolve delta image
        print('Convolving Delta image')
        conv_delta_map_tot = 0
        for i in range(self.n_b):
            for j in range(self.n_l):
                conv_delta_map_tot += sky_response_scaled[:, i, j, :]*delta_map_tot_old[i, j]
        
        # find maximum acceleration parameter to multiply delta image with
        # so that the total image is still positive everywhere
        print('Finding maximum acceleration parameter')
        try:
            len_arr = []
            for i in range(0,10000):
                len_arr.append(len(np.where((map_old + delta_map_tot_old*i/self.delta_map_scaling) < 0)[0]))
            len_arr = np.array(len_arr)
            afl = np.max(np.where(len_arr == 0)[0])
            print(f'Maximum acceleration parameter found: {afl/self.delta_map_scaling}')

             
            # fit delta map and current map to speed up RL algorithm
            print('Fitting delta-map in addition to old map')#, iteration {its}')
            # dictionary for data set and prior
            # note that here the value for N should be your response CDS dimension
            # should be last dimension of the scaled response thing (change to your value)
            data_multimap = dict(N = dataset.shape[1],
                         Nh = d2h,
                         Ncuts = Ncuts,
                         Nsky = 2,
                         acceleration_factor_limit = afl*0.95,
                         bg_cuts = bg_cuts,
                         bg_idx_arr = idx_arr,
                         y = dataset.ravel().astype(int),
                         bg_model = tmp_model_bg,
                         conv_sky = np.concatenate([[expectation_old], [conv_delta_map_tot/self.delta_map_scaling]]),
                         mu_flux = np.array([1, afl/2]),
                         sigma_flux = np.array([1e-2, afl]),
                         mu_Abg = self.fitted_bg,
                         sigma_Abg = self.fitted_bg)

            # initial values for fit (somewhat sensitive here with COSI data)
            init = {}
            init['flux'] = np.array([1., afl/2.])
            init['Abg'] = np.repeat(self.fitted_bg, Ncuts)#np.repeat(0.99,Ncuts)
            # fit: might take some time but it shouldn't be more than a minute
            op2D = model_multimap.optimizing(data=data_multimap, init=init, as_vector=False, verbose=True, tol_rel_grad=1e3, tol_obj=1e-20)

            # save values
            print('Saving new map, and fitted parameters')
            intermediate_lp = op2D['value']
            acc_par = op2D['par']['flux'][1]
            bg_pars = op2D['par']['Abg']
      
            # make new map as old map plus scaled delta map
            map_new = map_old + op2D['par']['flux'][1]*delta_map_tot_old/self.delta_map_scaling

            # same with expectation (data space)
            expectation_new = expectation_old + op2D['par']['flux'][1]*conv_delta_map_tot/self.delta_map_scaling
    
        except:
            # if the fit failed...
            # this shouldn't happen too often (or at all)
            print('############## Fit failed! proceeding without acceleration ##############')
            map_new = map_old + delta_map_tot_old
            expectation_new = expectation_old + conv_delta_map_tot
            #### break if the fit fails?
            raise ValueError("Fit failed! Please adjust parameters and try again")

        
        # check finite values again
        map_new[bad_expo[0], bad_expo[1]] = 0
        map_new[np.where(np.isnan(map_new) == True)] = 0


        return [map_new, expectation_new, map_likelihood, intermediate_lp, acc_par, bg_pars]


    def get_image_response_from_pixelhit_general(self, Response, cut, analysis1, pointing1, altitude_correction=False, al=None):
        """
        Get Compton response from hit pixels for each zenith/azimuth.

        Parameters:
             Response (np.ndarray, required): The sky response from COSIpy's response class
             cut (float, required): Angle from COSI's zenith [deg] at which to truncate the response 
             analysis1 (COSIpy.COSIpy object, required): Contains the data/simulation to be analyzed 
             pointing1 (COSIpy.COSIpy object, required): Contains the pointing information from COSIpy pointing class
             altitude_correction (bool, optional): Use interpolated transmission probability, normalised to 33 km and 500 keV,
                                                   to modify number of expected photons as a function of altitude and zenith 
                                                   angle of cdxervation. Default: False
             al (np.ndarray, optional): Altitude values according to dt from construct_pointings(); 
                                        Used if altitude_correction is set to True. Default: None

        Returns:
             image_response (np.ndarray): Scaled sky response. Shape: (# energy bins, n_b, n_l, FISBEL)
        """
        # Number of hours in cdxervation
        n_hours = analysis1.dataset.times.n_ph

        dt = pointing1.dtpoins

        # check which pixel (index) was hit on regular grid
        # zenith: Zenith positions of all points of predefined sky grid with respect to the instrument (in deg)
        # azimuth: Azimuth positions of all points of predefined sky grid with respect to the instrument (in deg)
        hit_pixel_zi = np.floor(self.zenith/self.pixel_size)
        hit_pixel_ai = np.floor(self.azimuth/self.pixel_size)

        # and which pixel centre
        hit_pixel_z = (hit_pixel_zi + 0.5)*self.pixel_size
        hit_pixel_a = (hit_pixel_ai + 0.5)*self.pixel_size

        # check which zeniths are beyond threshold
        bad_idx = np.where(hit_pixel_z > cut)#self.deg2rad)

        # set hit pixels to output array
        za_idx = np.array([hit_pixel_zi, hit_pixel_ai]).astype(int)

        nz = self.zenith.shape[2]

        # take care of regular grid by applying weighting with latitude
        weights = ((self.pixel_size*np.pi/180)*(np.sin(np.deg2rad(self.B_ARR + self.pixel_size/2)) - np.sin(np.deg2rad(self.B_ARR - self.pixel_size/2)))).repeat(nz).reshape(self.n_b, self.n_l, nz)
        weights[bad_idx] = 0

        
        # check for negative weights and indices and remove
        weights[za_idx[0, :] < 0] = 0.
        weights[za_idx[1, :] < 0] = 0.
        za_idx[0, za_idx[0, :] < 0] = 0.
        za_idx[1, za_idx[1, :] < 0] = 0.
        
        
        if altitude_correction == True:
            altitude_response = return_altitude_response()
        else:
            altitude_response = one_func

        # get responses at pixels       
        image_response = np.zeros((n_hours, self.n_b, self.n_l, Response.shape[2]))

        for c in tqdm(range(n_hours)):
            cdx = np.where((pointing1.cdtpoins > analysis1.dataset.times.times_min[analysis1.dataset.times.n_ph_dx[c]]) &
                           (pointing1.cdtpoins <= analysis1.dataset.times.times_max[analysis1.dataset.times.n_ph_dx[c]]))[0]
            
            # this calculation is basically a look-up of the response entries. In general, weighting (integration) with the true shape can be introduced, however with a lot more computation time (Simpson's rule in 2D ...)
            #altitude_weights = altitude_response(self.zenith[:,:,cdx].ravel(),al[cdx])[np.argsort(np.argsort(self.zenith[:,:,cdx].ravel()))].reshape(n_lat,n_lon)
            image_response[c, :, :, :] += np.sum(Response[za_idx[0, :, :, cdx], za_idx[1, :, :, cdx], :]*np.einsum('klij->iklj', weights[:, :, cdx, None])*dt[cdx, None, None, None], axis=0)#*altitude_weights[:,:,None]

        return image_response


    def calculate_exposure_map(self, sky_response_scaled):
        """
        Calculate the exposure map.

        Parameters:
             sky_response_scaled (np.ndarray, required): Scaled sky response. Shape: (# energy bins, n_b, n_l, FISBEL)

        Returns: 
             expo_map (np.ndarray): Exposure map. Shape: (n_b, n_l)
        """
        print('\nCalculating the exposure map')
        expo_map = np.zeros((self.n_b, self.n_l))

        for i in range(sky_response_scaled.shape[0]):
            expo_map += np.sum(sky_response_scaled[i, :, :, :], axis=2)

        return expo_map


    def set_bad_exposure(self, image, threshold=1e4):
        """
        Mask regions of bad exposure in a given array containing image data.

        Parameters:
             image (np.ndarray, required): One image to analyze. Shape (n_b, n_l)
             threshold (float, optional): Bad exposure criterion defined as expo_map/domega <= threshold. Default: 1e4

        Returns: 
             image_nan (np.ndarray): The given image array with pixels of bad exposure set to 'nan'.
        """
        image_nan = np.copy(image)
        bad_expo = np.where(expo_map/self.domega <= threshold)
        image_nan[bad_expo[0], bad_expo[1]] = np.nan

        return image_nan


    def plot_exposure_map(self, expo_map):
        """
        Plot the exposure map.

        Parameters:
             expo_map (np.ndarray, required): The exposure map to plot. Shape: (n_b, n_l)

        Returns: 
             Plot of the exposure map.
        """
        plt.subplot(projection='aitoff')

        # exposure map weighted by pixel size
        emap = np.roll(expo_map/self.domega, axis=1, shift=0)

        # exposure map, un-weighted
        #emap = np.roll(expo_map, axis=1, shift=0)

        # exposure map * domega = effective area?
        #emap = np.roll(expo_map*self.domega, axis=1, shift=0)

        p = plt.pcolormesh(self.L_ARRg*self.deg2rad, self.B_ARRg*self.deg2rad, emap)
        plt.contour(self.L_ARR*self.deg2rad, self.B_ARR*self.deg2rad, emap, colors='black')
        plt.colorbar(p, orientation='horizontal')
        plt.show()


    def get_total_image_flux(self, image):
        """
        Get total flux contained in an image.

        Parameters:
             image (np.ndarray, required): One image to analyze. Shape (n_b, n_l)

        Returns: 
             total_flux (float): Total flux in the given image
        """
        total_flux = np.sum(image*self.domega)

        return total_flux


    def get_pixel_flux(self, image, l, b, num_neighbors=1):
        """
        Find the flux at a given (l, b) coordinate of interest in the supplied image.
        Also find the flux in the neighboring pixels.
        
        Parameters:
            l (float, required): Galactic longitude of interest [deg]
            b (float, required): Galactic latitude of interest [deg]
            num_neighbors (int, optional): Number of neighboring pixels to examine in each direction.
                                           Default: 1, i.e. 1 pixels on either side of the pixel of 
                                                    interest for a total of 4 neighboring pixels
        
        Returns:
            pixels (dictionary): Dictionary of pixel boundaries and the flux contained therein. 
                                 key: l: (l0, l1), b: (b0, b1)
                                 value: flux in the specified pixel
            Primary pixel of interest is listed first, followed by neighboring pixels.
        """
        
        for i in range(1, len(self.l_arrg)):
            l0 = self.l_arrg[i-1]
            l1 = self.l_arrg[i]

            for j in range(1, len(self.b_arrg)):
                b0 = self.b_arrg[j-1]
                b1 = self.b_arrg[j]

                if (l0 < l <= l1) and (b0 < b <= b1):

                    pixels = {}
                    
                    pixels[f'l: ({l0}, {l1}), b: ({b0}, {b1})'] = image[j-1, i-1]
                    
                    for n in range(1, num_neighbors+1):
                        
                        # move up n
                        pixels[f'l: ({self.l_arrg[i-1]}, {self.l_arrg[i]}), b: ({self.b_arrg[j-1+n]}, {self.b_arrg[j+n]})'] = image[j-1+n, i-1]
                        
                        # move down n
                        pixels[f'l: ({self.l_arrg[i-1]}, {self.l_arrg[i]}), b: ({self.b_arrg[j-1-n]}, {self.b_arrg[j-n]})'] = image[j-1-n, i-1]
                        
                        # move right n
                        pixels[f'l: ({self.l_arrg[i-1+n]}, {self.l_arrg[i+n]}), b: ({self.b_arrg[j-1]}, {self.b_arrg[j]})'] = image[j-1, i-1+n]
                        
                        # move left n
                        pixels[f'l: ({self.l_arrg[i-1-n]}, {self.l_arrg[i-n]}), b: ({self.b_arrg[j-1]}, {self.b_arrg[j]})'] = image[j-1, i-1-n]
                        
                        # move right n, up n
                        pixels[f'l: ({self.l_arrg[i-1+n]}, {self.l_arrg[i+n]}), b: ({self.b_arrg[j-1+n]}, {self.b_arrg[j+n]})'] = image[j-1+n, i-1+n]
                        
                        # move left n, down n
                        pixels[f'l: ({self.l_arrg[i-1-n]}, {self.l_arrg[i-n]}), b: ({self.b_arrg[j-1-n]}, {self.b_arrg[j-n]})'] = image[j-1-n, i-1-n]

                        # move right n, down n
                        pixels[f'l: ({self.l_arrg[i-1+n]}, {self.l_arrg[i+n]}), b: ({self.b_arrg[j-1-n]}, {self.b_arrg[j-n]})'] = image[j-1-n, i-1+n]
                        
                        # move left n, up n
                        pixels[f'l: ({self.l_arrg[i-1-n]}, {self.l_arrg[i-n]}), b: ({self.b_arrg[j-1+n]}, {self.b_arrg[j+n]})'] = image[j-1+n, i-1-n]
                    
        
        return pixels


    def plot_image(self, image, iteration=None, mask_bad_expo=True, save_fig=False):
        """
        Plot an image, e.g. from the RL deconvolution.

        Parameters:
             image (np.ndarray, required): One image to analyze. Shape (n_b, n_l)
             iteration (int, optional): Iteration number of the given image. Default: None
             mask_bad_expo (bool, optional): Mask pixels with bad exposure as gray. Default: True
             save_fig (bool, optional): Save the figure as RL_image.pdf. Default: False

        Returns: 
             Image of specified array.
             If save_fig=True, also return RL_image.pdf
        """
        # bad exposures will be gray
        cmap = plt.get_cmap('viridis')
        cmap.set_bad('lightgray')

        if mask_bad_expo:
            # select here which pixels should be gray
            image = imaging.set_bad_exposure(image)

        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection':'aitoff'}, nrows=1, ncols=1)

        ax.set_xticks(np.array([-120, -60, 0, 60, 120])*self.deg2rad)
        ax.tick_params(axis='x', colors='orange')
        ax.set_xticklabels([r'$-120^{\circ}$'+'\n',
                            r'$-60^{\circ}$'+'\n',
                            r'$0^{\circ}$'+'\n',
                            r'$60^{\circ}$'+'\n',
                            r'$120^{\circ}$'+'\n'])
        ax.set_yticks(np.array([-60, -30, 0, 30, 60])*self.deg2rad)
        ax.tick_params(axis='y', colors='orange')

        img = ax.pcolormesh(self.L_ARRg*self.deg2rad, self.B_ARRg*self.deg2rad, image, cmap=plt.cm.viridis)

        cbar = fig.colorbar(img, orientation='horizontal')
        cbar.ax.set_xlabel(r'Flux [UNITS??? 10$^{-2}$ ph cm$^{-2}$ s$^{-1}$]')
            
        ax.grid()

        if iteration is not None:
            plt.title(f'Iteration {iteration}')

        plt.show()

        if save_fig:
             plt.savefig(f'RL_image.pdf', bbox_inches='tight')


    def plot_bg_pars(self, iterations, bg_pars):
        """
        Plot the background parameter of each RL iteration.

        Parameters:
             iterations (int, required): Number of iterations in the RL deconvolution

        Returns:
             Plot of background parameter vs. iteration number of the RL deconvolution
        """
        plt.figure(figsize=(14, 6))
        plt.plot(range(iterations-1), bg_pars, '.-')#self.bg_pars, '.-')
        plt.xlabel('Iteration')
        plt.ylabel('BG params')
        plt.show()


    def plot_flux(self, iterations, images):
        """
        Plot the flux of each RL iteration.

        Parameters:
             iterations (int, required): Number of iterations in the RL deconvolution
             images (np.ndarray, required): Array of images from the RL deconvolution. Shape: (n_b, n_l, iterations)

        Returns:
             Plot of flux vs. iteration number of the RL deconvolution
        """
        map_fluxes = np.zeros(iterations)

        for i in range(iterations):
            map_fluxes[i] = np.sum(images[:, :, i]*self.domega)
        
        plt.figure(figsize=(14, 6))
        plt.plot(range(iterations), map_fluxes, 'o-')
        plt.xlabel('Iteration')
        plt.ylabel('Flux [units?]')
        plt.show()


    def plot_likelihood(self, iterations, likelihood):
        """
        Plot the likelihood of each RL iteration.

        Parameters:
             iterations (int, required): Number of iterations in the RL deconvolution

        Returns:
             Plot of likelihood vs. iteration number of the RL deconvolution
        """
        plt.figure(figsize=(14, 6))
        plt.plot(range(iterations-1), likelihood, '.-')# self.intermediate_lp, '.-')
        plt.xlabel('Iteration')
        plt.ylabel('Likelihood')
        plt.show()


    # # this isn't right I don't think...
    # def plot_image_lightcurve(self, image, analysis1_binned):
    #     conv_map = 0
    #     for i in range(self.n_b):
    #         for j in range(self.n_l):
    #             if np.isfinite(image[i, j]) == True:
    #                 conv_map += sky_response_scaled[:, i, j, :]*image[i, j]

    #     plt.plot(np.sum(conv_map, axis=1), label='Image counts')
    #     plt.plot(analysis1_binned, label='Data counts')
    #     plt.xlabel('Time')
    #     plt.ylabel('Expected counts')
    #     plt.show()


    def create_gif_of_iterations(self, images):
        """
        Create gif of deconvolution.

        Parameters:
             images (list, required): List of images to combine as a GIF

        Returns: 
             RL_image_iterations.gif, GIF of images
        """
        from matplotlib import animation

        cmap = plt.get_cmap('viridis') 

        # bad exposures will be gray
        cmap.set_bad('lightgray')

        # plotting for GIF
        fig, ax = plt.subplots(figsize=(10, 7), subplot_kw={'projection':'aitoff'}, nrows=1, ncols=1)

        ax.set_xticks(np.array([-120, -60, 0, 60, 120])*self.deg2rad)
        ax.tick_params(axis='x', colors='orange')
        ax.set_xticklabels([r'$-120^{\circ}$'+'\n',
                            r'$-60^{\circ}$'+'\n',
                            r'$0^{\circ}$'+'\n',
                            r'$60^{\circ}$'+'\n',
                            r'$120^{\circ}$'+'\n'])
        ax.set_yticks(np.array([-60, -30, 0, 30, 60])*self.deg2rad)
        ax.tick_params(axis='y', colors='orange')

        # ims is a list of lists, each row is a list of artists to draw in the
        # current frame; here we are just animating one artist, the image, in each frame
        ims = []
        for i in range(images.shape[-1]):
            ttl = plt.text(0.5, 1.01, r'RL iteration {0:1.0f}'.format(i), horizontalalignment='center', 
                           verticalalignment='bottom', transform=ax.transAxes)
            
            image = images[:, :, i]#np.roll(images[:, :, i], axis=1, shift=0) 
            img = ax.pcolormesh(self.L_ARRg*self.deg2rad, self.B_ARRg*self.deg2rad, image, cmap=plt.cm.viridis)    
            ims.append([img, ttl])

        cbar = fig.colorbar(img, orientation='horizontal')
        #cbar.ax.set_xlabel(r'Flux [10$^{-2}$ ph cm$^{-2}$ s$^{-1}$]')
            
        ax.grid()         
            
        ani = animation.ArtistAnimation(fig, ims, interval=350, blit=True, repeat_delay=0)
        ani.save('RL_image_iterations.gif')



###########################################


# Read the data
data_dir = '/volumes/eos/users/jacqueline/analysis/Al26/'
#tra = 'RL/point_source/fixed_ori/l-3_b-3/point_l-3_b-3_1809keV.inc1.id1.Mimrec_v3.tra'
#tra = 'RL/point_source/fixed_ori/Gaussian_2D/FarFieldGauss/X.tra'
#tra = 'RL/DataChallenge/Al26/Al26.10xFlux.inc1.id1.extracted.Mimrec1809_v3.tra'
tra = 'RL/DataChallenge/DC_Al2610xFlux_Ling_Mimrec1809_v3.tra'

analysis1 = COSIpy(data_dir, tra)
analysis1.read_COSI_DataSet()


# Define time binning
Delta_T = 1800 # seconds
analysis1.dataset.time_binning_tags(time_bin_size=Delta_T)


# Define energy binning
energy_bin_edges = np.array([1803, 1817])
ebin = 0 # energy bin of interest (only 1 bin here)


# Bin the data
pixel_size = 6.

analysis1.dataset.init_binning(energy_bin_edges=energy_bin_edges, pixel_size=pixel_size)
print('\nBinning the data')
analysis1.dataset.get_binned_data()
print(f'\nShape of binned data (time, energy, phi (Compt. scatt. angle), FISBEL): \n \t{analysis1.dataset.binned_data.shape}')


# Define pointings
pointing1 = Pointing(dataset=analysis1.dataset,)


# Initialize the class now that you have the pointing object.
# Optionally specify delta_map_scaling=2000, init_map_flux=1.0, fitted_bg=np.array([0.2])
imaging = Imaging(
    pixel_size=pixel_size,
    delta_map_scaling=2000, 
    init_map_flux=0.01, 
    #fitted_bg=np.array([1E-6])) # 1E-6 for no BG in the dataset
    fitted_bg=np.array([0.73]))


# Define the background
print('\nDefining the background model')
background1 = BG(dataset=analysis1.dataset,mode='sim 6deg despina')
nonzero_idx = background1.calc_this[ebin]


# Load the response
rsp = response.SkyResponse(filename='/volumes/eos/users/jacqueline/analysis/Al26/RL/Isotropic.1809keV.coriv1v2v3saviov1.binnedimaging.imagingresponse.npz',
                           pixel_size=pixel_size)


# Truncate the response [cut] degrees off-axis (COSI FOV ~ 60 deg)
cut = 60.


# Get the sky response
sky_response_CDS = rsp.rsp.response_grid_normed_efinal.reshape(
    imaging.n_b,
    imaging.n_l,
    analysis1.dataset.phis.n_phi_bins*\
    analysis1.dataset.fisbels.n_fisbel_bins, 1)[:, :, nonzero_idx, ebin]

sky_response_scaled = imaging.get_image_response_from_pixelhit_general(
    Response=sky_response_CDS,
    analysis1 = analysis1,
    pointing1 = pointing1,
    cut=cut,
    altitude_correction=False,
    al=np.ones(len(pointing1.dtpoins)))


# Calculate the exposure map
expo_map = imaging.calculate_exposure_map(sky_response_scaled=sky_response_scaled)


# Plot exposure map
imaging.plot_exposure_map(expo_map=expo_map)


# Compile the model
try:
    #read COSImodefit.pkl (if already compiled)
    model_multimap = pickle.load(open('fit_COSI_conved_2D_new2_multimap_indivBG.pkl', 'rb'))

except:
    print('\nModel not yet compiled, doing that now (might take a while).')
    ## compile model (if not yet compiled):
    model_multimap = pystan.StanModel('fit_COSI_conved_2D_multimap_indivBG.stan')

    
    # for line response, save it to the file 'filename.pkl' for later use
    with open('/volumes/eos/users/jacqueline/analysis/Al26/RL/fit_COSI_conved_2D_new2_multimap_indivBG.pkl', 'wb') as f:
        pickle.dump(model_multimap, f)


# Run imaging
iterations = 10

print('\nStarting the RL image deconvolution:')
print(f'\t pixel size = {pixel_size} deg\n\
    \t Number of iterations = {iterations}\n\
    \t Response cut = {cut} deg\n')


# Initial image
image = imaging.current_image

# Convolve map with the response
expectation_init = 0
print('Convolving with response (init expectation), iteration 0')
for i in range(imaging.n_b):
    for j in range(imaging.n_l):
        expectation_init += sky_response_scaled[:, i, j, :]*image[i, j]


# Iterate over images
RL_image_array = np.zeros((imaging.n_b, imaging.n_l, iterations))
RL_image_array[:, :, 0] = image

# likelihood of maps (vs. initial i.e. basically only background)
map_likelihood = []

# fit likelihoods, ie fit quality
fit_likelihood = []

# acceleration parameters
acc_par = []

# fitted bg parameters
bg_pars = []

for i in range(1, iterations):
    next_ = imaging.next_iteration(
        map_init = RL_image_array[:, :, i-1],
        expectation_init = expectation_init,
        sky_response_scaled = sky_response_scaled,
        expo_map = expo_map,
        background1 = background1,
        analysis1 = analysis1,
        show_its=True,
        iteration=i)

    new_image = next_[0]
    new_expectation = next_[1]

    map_likelihood.append(next_[2])
    fit_likelihood.append(next_[3])
    acc_par.append(next_[4])
    bg_pars.append(next_[5])

    RL_image_array[:, :, i] = new_image
    expectation_init = new_expectation


# # Save the array of images as a .npy file for later use
# (could also save the likelihood and other parameter lists)
np.save("RL_images", RL_image_array)



###### Load the saved array of images
RL_image_array = np.load("RL_images.npy")


# Plot an image of interest
idx = iterations-1
image = RL_image_array[:, :, idx]
imaging.plot_image(image, iteration=idx)


# Get total image flux
image_flux = imaging.get_total_image_flux(image)
print(f'\nImage flux: {image_flux} UNITS???')


# Get flux in pixel of interest (and neighboring pixels)
l0 = -3.
b0 = -3.
pixel_flux = imaging.get_pixel_flux(image=image, l=l0, b=b0, num_neighbors=1)
print(f'Flux at (l, b) = ({l0}, {b0}) deg [flux UNITS????]:')
for k, v in sorted(pixel_flux.items(), key=lambda x: x[1], reverse=True):
    print(k, v)


# Plot BG parameters
imaging.plot_bg_pars(iterations=iterations, bg_pars=bg_pars)


# Plot flux
imaging.plot_flux(iterations=iterations, images=RL_image_array)


# Plot map likelihood
imaging.plot_likelihood(iterations=iterations, likelihood=map_likelihood) 

# Plot fit likelihood
imaging.plot_likelihood(iterations=iterations, likelihood=fit_likelihood) 


# Make GIF of RL iterations
print('\nSaving gif of RL iterations...')
imaging.create_gif_of_iterations(RL_image_array)
print('...done')