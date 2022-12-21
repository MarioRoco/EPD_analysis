import datetime as dt
from pathlib import Path
from astropy import units
import numpy as np
import matplotlib.pyplot as plt

from SIS_auxfuncs import * # import all functions and variables of SIS_auxfuncs

#input that needs to be adapted by user:
#infile_path_SISdata="/media/marior/Seagate Expansion Drive/MASTER/TFM/SOARdata/EPD/SIS/"
infile_path_SISdata="/media/marior/47BF-82B6/aa_PAPER_SolarOrbiter/SOARdata/EPD/SIS/"





class SIS:
    
    def clarificarions(self):
        """
        This function returns the next clarifications:
        
            - "elem": ID of one ion specie. Type: string. It can be: 'H', 'He3', 'He4', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca' or 'Fe'.
            
            - "channels": the number of energy channels you want. It can be an array, a tuple or a 1-dimensional array. Each item should be an integer, the lowest value should be 0 and the largest depends on the ion species, for example for hydrogen it's 10. Example: "channels = [0,1,2,4,7,11]". 
            
            - "sampling_factor": it defines the time resolution you want to use, for example in a time series. The new time resolution in the data is sampling_factor multiplied by the original time resolution (the original time resolution is the time resolution of the data from Solar Orbiter Archive). For example, if the original time resolution is 30 seconds and you define sampling_factor=120, the new time resolution is 30*120 seconds, i.e. 1 hour. 
            
            - "limits": initial and final dates and hours where you want to restrict the time series, or calculation of spectrum, or another analysis or plot... It should be a list or tuple with 2 items, each item is another list or tuple with 4 (or 5) sub-items which are in this order: year, month, day, hour (and minute). For example: "limits=[(2022,3,5,12), (2022,4,17,22,47)]" if you want to cut from 5th March 2022 at 12:00 am to 17th April at 22:47 hrs. The dates should belong to the time interval described by the input variables "first_date" and "last_date". 
            
            - "selected dates": when we talk about "selected dates" in the comments, we mean: from "first_date" to "last_date", or from "first_date_quiet" to "last_date_quiet" in case of quiet time. 
            
            - "varget" is a function of "cdflib" ("cdflib" explained below) that returns the data inside a variable of the file that you load. An example of code:
                import cdflib 
                file_name = cdflib.CDF('/path/to/cdf_file.cdf')
                variable1 = file_name.varget('He4_Flux') # we load the variable "He4_Flux" that is in the file "cdf_file.cdf". We can do this also: "variable1 = file_name['He4_Flux']"
                
            - "cdflib": "cdflib is a python module to read/write CDF (Common Data Format .cdf) files without needing to install the CDF NASA library." (https://pypi.org/project/cdflib/)
            
            - "flux" and "intensity": in this code when we talk about fluxes and about uncertainties, they are the same. 
                
            - Nils said that in spectrum I should sum the counts (rates * time step) or intensities (fluxes) during the selected interval of time, but in time series I should do the average of counts (rates * time step) or intensities (fluxes) for each time step (time resolution). 
            
            - class, object and module. When you do this "from SIS_class import SIS", you are importing the class "SIS" from the module "SIS_class". And when you do "sis = SIS(direction='sun', elements=['H','He3','He4','O','Fe'], first_date=(2022,3,5), last_date=(2022,3,10), first_date_quiet=(2020,11,7), last_date_quiet=(2020,11,12), show_loaded_files='no')", you are creating the object "sis" of the class "SIS". 
            
        """
    
    def clean_data(self, rate, rate_uncertainty, flux, flux_uncertainty, Elements):
        """
        We put a mask in values smaller than zero (for rates, fluxes and their uncertainties)
        """
        
        rate_masked = {} # They are dictionaries whose items corresponds to one ion species. 
        rate_uncertainty_masked = {}
        flux_masked = {}
        flux_uncertainty_masked = {}
        for elem in Elements:
            rate_masked[elem] = np.ma.masked_array(rate[elem], rate[elem]<0)
            rate_uncertainty_masked[elem] = np.ma.masked_array(rate_uncertainty[elem], rate_uncertainty[elem]<0)
            flux_masked[elem] = np.ma.masked_array(flux[elem], flux[elem]<0)
            flux_uncertainty_masked[elem] = np.ma.masked_array(flux_uncertainty[elem], flux_uncertainty[elem]<0)

        return [rate_masked, rate_uncertainty_masked, flux_masked, flux_uncertainty_masked]
        
    
    def __init__(self, direction, elements, first_date, last_date, first_date_quiet, last_date_quiet, show_loaded_files):
        """
        "The special method __init__ is the Python constructor." (https://www.udacity.com/blog/2021/11/__init__-in-python-an-overview.html)
        
        
        Input variables:
            
            - "direction": the sensor EPD/SIS has two different directions: sunward or anti-sunward; the input variable "direction" should be a string: 'sun' (if you want the data of the "sunward" direction) or 'asun' or 'antisun' (if if you want the data of the "anti-sunward" direction). 
            
            - "elements": is a list or a tuple whose items are strings corresponding to the ions that the sensor SOLO/EPD/SIS con measure, so: 'H', 'He3', 'He4', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca' and/or 'Fe'. 
            
            - "first_date", "last_date", "first_date_quiet" and "last_date_quiet": correspond to: first day of the data you want to load, last day of the data you want to load, first day of the data corresponding to quiet time, and last day of the data corresponding to quiet time, respectively. Each one should be a tuple of three items in the format (year, month, day). For instance (2022,3,10) which corresponds to the date March 10, 2022. 
            
            - "show_loaded_files": it should be the string 'yes' if you want the class to show you if each file has been loaded correctly or not. Otherwise this variable should be something different to the string 'yes'. 
        
        
        Output variables:
            
            - "self.first_date_qt", "self.last_date_qt": they are the input variables "first_date_quiet" and "last_date_quiet", respectively. 
            - "self.elements": it is the input variable "elements". 
            - "self.original_time_resolution": (list, seconds) list that contains the time resolution of the original files. Each item of the list is the time resolution of the first "epoch" item of each day. For the files 'solo_L2_epd-sis-a-rates-medium_' and 'solo_L2_epd-sis-b-rates-medium_' the original time resolution is always 30 seconds (I've checked it). 
            - "self.cdf_files": the files of the selected period of time and that exist in the corresponding path or folder. 
            - "self.n_days": (integer) number of selected days (but not counting missing files). self.n_days=len(self.cdf_files) 
            - "self.Epoch": (list, dates) list with the epochs of all the selected days. 
            - "self.Epoch_sec": (list, seconds) list with the epochs of all the selected days, but in seconds instead of dates format. 
            - "self.Flux": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) each item of this dictionary corresponds to one ion especie and contains a 2-dimensional array with the fluxes (or intensities). In these 2-dimensional arrays the columns represent the energy channels and the rows represent the epochs. 
            - "self.Flux_unc": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) same as "self.Flux" but instead of fluxes are uncertainties of the fluxes. 
            - "self.Rate": (dictionary, [counts/second]) same as "self.Flux" but instead of fluxes are rates. 
            - "self.Rate_unc": (dictionary, [counts/second]) same as "self.Flux" but instead of fluxes are uncertainties of the rates. 
            - "self.Enuc_low": (array, [MeV/nucleon]) lower boundary of the energy channels. 
            - "self.Enuc_width": (array, [MeV/nucleon]) energy width of the energy channels. 
            - "self.Enuc_up": (array, [MeV/nucleon]) upper boundary of the energy channels (Up = Low + Width). 
            - "self.Enuc": (array, [MeV/nucleon]) sqrt(Low^2 * Up^2)
            - "self.original_time_resolution_quiet", "self.cdf_files_quiet", "self.n_days_quiet", "self.Epoch_quiet", "self.Flux_quiet", "self.Flux_unc_quiet", "self.Rate_quiet", "self.Rate_unc_quiet", "self.Enuc_low_quiet", "self.Enuc_width_quiet", "self.Enuc_up_quiet", "self.Enuc_quiet": they are similar to their above respective variables (but without "_quiet") but corresponding to the quiet time. 
            
            
        The code I usually use to execute this class SIS is:
            import SIS_class as SIS #so: "import name_of_this_file" , or "import name_of_this_file as another_name"
            sis = SIS(direction, elements, first_date, last_date, first_date_quiet, last_date_quiet, show_loaded_files)
        """
        
        self.first_date_qt, self.last_date_qt = first_date_quiet, last_date_quiet
        self.elements = elements
        date1, date2 = [first_date, first_date_quiet], [last_date, last_date_quiet]
        
        for noquiet_quiet in [0,1]: #0 corresponds to period of selected days and 1 corresponds to quiet time
            
            # load the files of the selected dates (you can read what means "selected dates" using the function above "clarifications"). 
            ## select the direction
            if direction == 'sun':
                inf_corname = 'solo_L2_epd-sis-a-rates-medium_'
            elif direction == 'asun' or direction == 'anti-sun':
                inf_corname = 'solo_L2_epd-sis-b-rates-medium_'
            else: print("direction should be 'sun', 'asun' or 'antisun'")
            ## load the files of all the selected days. 
            Files_of_days = load_datefiles(date1[noquiet_quiet], date2[noquiet_quiet], infile_path=infile_path_SISdata, infile_corname=inf_corname, show_loaded_files=show_loaded_files) # "load_datefiles" is a function of the file "SIS_auxfuncs.py" that returns the existing files of the selected dates. It returns a list which contains "[files, n_days, first_datetime]"
            self_cdf_files = Files_of_days[0] # the files of the selected period of time and that exist in the corresponding path or folder. 
            self_n_days = len(self_cdf_files) # number of days
            
            
            # load the data: epoch, flux, rate, flux uncertainties, rate uncertainties
            i=0
            Epoch, Epoch_sec, Flux, Flux_unc, Rate, Rate_unc, original_time_resolution=[],[],{},{},{},{},[] # empty lists and dictionaries
            for elem in self.elements:
                Flux[elem], Flux_unc[elem], Rate[elem], Rate_unc[elem] = [], [], [], [] # Filling the above dictionaries with pairs of selected elements (of the variable self.elements) and empty lists. e.g.: {'H': [] , 'He3': [] , ...}
                
            while i<self_n_days: # we go day by day. 
                original_time_resolution.append(self_cdf_files[i].varget("DELTA_EPOCH")[0]) # (list, seconds). filling the list that contains the time resolution of the original files. Each item of the list is the time resolution of the first "epoch" item of each day. Units: [seconds]. 

                Epoch_nanosec = self_cdf_files[i].varget("EPOCH") # (array, nanoseconds). array whose items are the epochs of the day corresponding to the index "i". Units: [nanoseconds] 
                Epoch_seconds = epoch_to_unixseconds(Epoch_nanosec) # (array, seconds). same as Epoch_nanosec but in seconds
                Epoch_date = unixseconds_to_date(Epoch_seconds) # (list, dates). same as  Epoch_nanosec but in date format (datetime.datetime(year, month, day,...), however this is a list, not an array. 
                Epoch.append(Epoch_date) # (list, dates). each item of this list is an array of dates corresponding to one day. 
                Epoch_sec.append(Epoch_seconds) # (list, seconds). same as "Epoch" but in seconds instead of date format. 

                for elem in elements: # we go elem by elem, in the day corresponding to index "i".
                    Rate_1day = self_cdf_files[i].varget("%s_Rate"%(elem)) # (2-dimensional array, [counts/second]). array with the rates of the day corresponding to the index "i" and of the ion_specie "elem". In this 2-dimensional array the columns represent the energy channels and the rows represent the epochs.  
                    #Rate_1day = Rate_1day * original_time_resolution[i] # Before I multiplied Rate_1day by the original time resolution, thus I had the total number of counts during each time step. Later in the time series I sum all the counts during the resampling. However it is better to do the average, thus it is incorrect to multiply by the original time resolution. Anyway I maintain this line "Rate_1day = Rate_1day * original_time_resolution[i]" (but with a "#") just in case I calculate the sum of counts in the future. If I change this, I should also change all the comments in which I wrote "[counts/second]" to "[counts/original_time_resolution]"). 
                    Flux_1day = self_cdf_files[i].varget("%s_Flux"%(elem)) # (2-dimensional array, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]). same as "Rate_1day" but with intensities. 
                    Rate[elem].append(Rate_1day) # (list, lists of 2-dimensional arrays of rates [counts/second]). A list in which each item is a 2-dimensional array of the rates of one day (each item is one "Rate_1day" (see above)). 
                    Rate_unc[elem].append(np.sqrt(Rate_1day) / np.sqrt(original_time_resolution[i])) # (list, lists of 2-dimensional arrays). same as Rate[elem] but with the uncertainties of the rates
                    Flux[elem].append(Flux_1day) # (list, lists of 2-dimensional arrays). same as Rate[elem] but with intensities
                    #Flux_unc_nan = Flux_1day/np.sqrt(Rate_1day) #To compare with below way to calculate uncertainties. Yes, I have compare them and they are similar, so the uncertainties are well calculated. 
                    Flux_unc_nan = self_cdf_files[i].varget("%s_UNCERTAINTY"%(elem)) # (2-dimensional array, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$])
                    Flux_unc_nan[Rate_1day == 0] = 0 # (2-dimensional array, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$])
                    Flux_unc[elem].append(Flux_unc_nan) # (list, lists of 2-dimensional arrays). same as Rate[elem] but with the uncertainties of the intensities
                i=i+1
                
                
            self_Epoch = np.concatenate(Epoch) # (list, dates). we concatenate the above array "Epoch" to obtain a list with the epochs of all the selected days. 
            self_Epoch_sec = np.concatenate(Epoch_sec) # (list, seconds). same as self_Epoch, but in seconds instead of dates. 
            self_Flux, self_Flux_unc, self_Rate, self_Rate_unc = {},{},{},{} # we create these empty dictionaries. 
            for elem in elements:
                # we fill each item of the dictionaries. Each item corresponds to one ion especie. 
                self_Flux[elem] = np.vstack(np.array(Flux[elem])) # (list, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$])
                self_Flux_unc[elem] = np.vstack(np.array(Flux_unc[elem])/2) # (list, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$])
                self_Rate[elem] = np.vstack(np.array(Rate[elem])) # (list, [counts/second])
                self_Rate_unc[elem] = np.vstack(np.array(Rate_unc[elem])/2)  # (list, [counts/second])
                
            
            self_Enuc_low, self_Enuc_width, self_Enuc_up, self_Enuc = {},{},{},{} # we create these empty dictionaries
            for elem in elements:
                # we fill each item of the dictionaries. Each item corresponds to one ion especie. 
                self_Enuc_low[elem] = self_cdf_files[0].varget("%s_Bins_Low_Energy"%(elem)) # (array, [MeV/nucleon]). Lower boundary of the energy channels. 
                self_Enuc_width[elem] = self_cdf_files[0].varget("%s_Bins_Width"%(elem)) # (array, [MeV/nucleon]). Energy width of the energy channels. 
                self_Enuc_up[elem] = self_Enuc_low[elem] + self_Enuc_width[elem]  # (array, [MeV/nucleon]). Upper boundary of the energy channels (Up = Low + Width). 
                self_Enuc[elem] = np.sqrt(self_Enuc_low[elem] * self_Enuc_up[elem])  # (array, [MeV/nucleon]). sqrt(Low^2 * Up^2)
                
            #Put a mask in false values of rate, flux and their uncertainties. For it we use the function "clean_data" which is above in this file. 
            self_Rate_masked, self_Rate_unc_masked, self_Flux_masked, self_Flux_unc_masked = self.clean_data(self_Rate, self_Rate_unc, self_Flux, self_Flux_unc, elements)
            
            
            #Separate data from the selected days and from quiet time
            if noquiet_quiet == 0: # selected days
                self.original_time_resolution = original_time_resolution
                self.cdf_files, self.n_days, self.Epoch, self.Epoch_sec = self_cdf_files, self_n_days, self_Epoch, self_Epoch_sec
                self.Flux, self.Flux_unc, self.Rate, self.Rate_unc = self_Flux_masked, self_Flux_unc_masked, self_Rate_masked, self_Rate_unc_masked
                self.Enuc_low, self.Enuc_width, self.Enuc_up, self.Enuc = self_Enuc_low, self_Enuc_width, self_Enuc_up, self_Enuc
            
            elif noquiet_quiet == 1: # quiet time
                self.original_time_resolution_quiet = original_time_resolution
                self.cdf_files_quiet, self.n_days_quiet, self.Epoch_quiet = self_cdf_files, self_n_days, self_Epoch
                self.Flux_quiet, self.Flux_unc_quiet, self.Rate_quiet, self.Rate_unc_quiet = self_Flux_masked, self_Flux_unc_masked, self_Rate_masked, self_Rate_unc_masked
                self.Enuc_low_quiet, self.Enuc_width_quiet, self.Enuc_up_quiet, self.Enuc_quiet = self_Enuc_low, self_Enuc_width, self_Enuc_up, self_Enuc
        
        
    def timeseries_fluxes_rates(self, sampling_factor, limits):
        """
        - Input variables:
            
            - "sampling_factor": it is the variable that modifies the time resolution from the original one to the time resolution you want in the time series. Example: If the original time resolution is 30 seconds and you want the time series to have a time step of 1 hour, sampling_factor should be 1h/30sec = 60min/0.5min = 120. 
            
            - "limits": it is a list or a tuple with 2 items. Each item is a tuple or a list that represent a date in format (yyyy,mm,dd,hh) or (yyyy,mm,dd,hh,min). The second date shoud be larger than the first one because these dates are the limits where you want to cut the time series. Example: "limits=((2022,2,12), (2022,2,6))"
            
        - Output variables: 
        
            - self.ind1: index of the first date where you want to cut.
            
            - self.ind2: index of the second date where you want to cut. 
            
            - self.Epoch_ts: (1-dimensional array, dates) Resampled epoch. The format of each item is date as datetime.datetime
            
            - self.Epoch_ts_sec: (1-dimensional array, seconds) Resample Epoch but in seconds
            
            - self.Epoch_ts_day: (1-dimensional array, seconds) Resample Epoch but in days
            
            - "self.Flux_ts": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) each item of this dictionary corresponds to one ion especie and contains a 2-dimensional array with the fluxes (or intensities) averaged over a time step of sampling_factor*original_time_resolution. In these 2-dimensional arrays the columns represent the energy channels and the rows represent the epochs (resampled). 
            
            - "self.Flux_unc_ts": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) uncertainties of "self.Flux_ts". 
            
            - "self.Rate_ts": (dictionary, counts/second) same as "self.Flux_ts" but with rates. 
            
            - "self.Rate_unc_ts": (dictionary, counts/second) same as "self.Flux_unc_ts" but with rates. 
            
            
        """
        # 1) EPOCH:
        # 1.1) we cut the epoch:
        ind1 = select_date_index(self.Epoch, cutoff_date=limits[0]) #index of the first date where you want to cut
        ind2 = select_date_index(self.Epoch, cutoff_date=limits[1]) #index of the second date where you want to cut
        self.ind1_ts, self.ind2_ts = ind1, ind2
        Epoch_cut = self.Epoch[ind1:ind2+1]
        Epoch_sec_cut = self.Epoch_sec[ind1:ind2+1]
        
        # 1.2) We reshample the epoch that we've cut above:
        self.Epoch_ts = resample_time(Epoch_cut, sampling_factor) # (1-dimensional array, dates) Resampled epoch. The format of each item is date as datetime.datetime
        self.Epoch_ts_sec = resample_time(Epoch_sec_cut, sampling_factor) # (1-dimensional array, seconds) Resample Epoch but in seconds
        self.Epoch_ts_day = self.Epoch_ts_sec / 86400 # (1-dimensional array, seconds) Resample Epoch but in days
        
        # 2) FLUXES and RATES
        # 2.0) create empty dictionaries. Each item will correspond to an ion species:
        self.Flux_ts, self.Flux_unc_ts, self.Rate_ts, self.Rate_unc_ts = {},{},{},{} 
        for elem in self.elements: # for the ion species "elem":
            
            # 2.1) Select the rows (which correspond to epochs) chosen by "limits" variable and all energy channels, this is why "[ind1:ind2+1, :]"
            Flux_cut = self.Flux[elem][ind1:ind2+1, :]
            Flux_unc_cut = self.Flux_unc[elem][ind1:ind2+1, :] 
            Rate_cut = self.Rate[elem][ind1:ind2+1, :]
            Rate_unc_cut = self.Rate_unc[elem][ind1:ind2+1, :]
            
            # 2.2) Resample the data
            self.Flux_ts[elem] = resample_flux(Flux_cut, sampling_factor)
            self.Flux_unc_ts[elem] = resample_fluxerr(Flux_unc_cut, sampling_factor)
            self.Rate_ts[elem] = resample_flux(Rate_cut, sampling_factor)
            self.Rate_unc_ts[elem] = resample_fluxerr(Rate_unc_cut, sampling_factor)
            
            # All these above variables have been checked and they are correct. 
    
            
            
    def timeseries_fluxes_rates_sum(self, channels, sampling_factor, limits):
        
        """
        This method sums the intensities of the selected energy channels (whose index are in the input variable "channels"). 
        
        - Input variables:
        
            - "channels": it is an array, list or tuple with the index of the energy channels you want to select. The first energy channels (with lowest energy) has index 0, the second has index 1, the third has index 2,...
            
        - Output variables:
        
            - "self.Flux_ts_sum": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) we sum the fluxes (or intensities) of the selected energy channels in each time step. So each item of this dictionary corresponds to an ion species and is the sum of fluxes (or intensities) of the selected energy channels for each time step (time step of width sampling_factor*original_time_resolution). 
            
            - "self.Flux_unc_ts_sum": (dictionary, [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]) uncertainties of "self.Flux_ts_sum". 
            
            - "self.Rate_ts_sum": (dictionary, counts/second) same as "self.Flux_ts_sum" but with rates. 
            
            - "self.Rate_unc_ts_sum": (dictionary, counts/second) uncertainties of "self.Rate_ts_sum". 
            
        """
        
        self.timeseries_fluxes_rates(sampling_factor, limits) #Initialize this method, thus we have the variables self.Flux_ts, self.Flux_unc_ts, self.Rate_ts and self.Rate_unc_ts that we need. 
        
        self.Flux_ts_sum, self.Flux_unc_ts_sum, self.Rate_ts_sum, self.Rate_unc_ts_sum = {},{},{},{}
        for elem in self.elements: # for the ion species "elem":
            
            self.Flux_ts_sum[elem] = self.Flux_ts[elem][:,channels].sum(axis=1)
            self.Flux_unc_ts_sum[elem] = np.sqrt(((self.Flux_unc_ts[elem][:,channels])**2).sum(axis=1))
            self.Rate_ts_sum[elem] = self.Rate_ts[elem][:,channels].sum(axis=1)
            self.Rate_unc_ts_sum[elem] = np.sqrt(((self.Rate_unc_ts[elem][:,channels])**2).sum(axis=1))
            
    
    def max_impulsive(self, elem, sampling_factor, limits, channels):
        
        """
        - Output variables: self.epoch_max_list, self.epoch_max_low, self.epoch_max_up, self.ind_max, self.flux_max, self.d_flux_max_low, self.flux_max_up, self.epoch_max_days, self.d_epoch_max_days_low, self.d_epoch_max_days_up
        """
        
        ind_max_list, flux_max_list, d_flux_max_low_list = [], [], []
        self.epoch_max, self.epoch_max_low, self.epoch_max_up, epoch_max_days_list, d_epoch_max_days_low_list, d_epoch_max_days_up_list = [], [], [], [], [], []
        j=0
        for ch in channels:
            
            if type(limits[0][0]) is int:
                lim = limits
            elif type(limits[0][0]) is tuple or type(limits[0][0]) is list:
                lim = limits[j]
            self.timeseries_fluxes_rates(sampling_factor=sampling_factor, limits=lim)
            Flux_j = self.Flux_ts[elem][:,ch]
            d_Flux_j = self.Flux_unc_ts[elem][:,ch]
            epoch_days = self.Epoch_ts_sec/86400
            xunc = x_unc(epoch_days, Flux_j, d_Flux_j) 
            ind_max = xunc['ind_max']
            ind_max_list.append(ind_max)
            flux_max_list.append(xunc['y_max'])
            d_flux_max_low_list.append(xunc['dy_max_low']) # the values are the same for up, so we do not create "flux_max_list_low"
            self.epoch_max.append(self.Epoch_ts[ind_max]) # each item of this list is the epoch (or time, in datetime.datetime format) corresponding to the maximum intensity of the corresponding energy channel
            self.epoch_max_low.append(unixdays_to_date(xunc['x1'])) # list of dates corresponding to the lower boundary uncertainties of "self.epoch_max". 
            self.epoch_max_up.append(unixdays_to_date(xunc['x2'])) # list of dates corresponding to the upper boundary uncertainties of "self.epoch_max"
            epoch_max_days_list.append(xunc['x_max'])
            d_epoch_max_days_low_list.append(xunc['x_max'] - xunc['x1'])
            d_epoch_max_days_up_list.append(xunc['x2'] - xunc['x_max'])
            j=j+1

        self.ind_max = np.array(ind_max_list) # 1-dimensional array with the index of the maximum intensity of each energy channel of the ion specie "elem"
        self.flux_max = np.array(flux_max_list) # 1-dimensional array with the values of the maxima intensities of the corresponding energy channels
        self.d_flux_max_low = np.array(d_flux_max_low_list) # 1-dimensional array with the values of upper uncertainties of "self.flux_max"
        self.d_flux_max_up = self.d_flux_max_low # 1-dimensional array with the values of lower uncertainties of "self.flux_max"
        self.epoch_max_days = np.array(epoch_max_days_list) # 1-dimensional array with the same as "self.epoch_max" (see above) but the format of each item is number of day instead of datetime.datetime
        self.d_epoch_max_days_low = np.array(d_epoch_max_days_low_list)  # 1-dimensional array with the lower uncertainties of "self.epoch_max_days" (format: day number)
        self.d_epoch_max_days_up = np.array(d_epoch_max_days_up_list)  # 1-dimensional array with the lower uncertainties of "self.epoch_max_days" (format: day number)
        
    def plot_timeseries_fluxes_rates(self, elements, sampling_factor, limits, channels_list, plot_max_impulsive='no', fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=sub_ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size):
        
        n_elem = len(elements)
        fig, ax = plt.subplots(nrows=n_elem, ncols=1, figsize=fig_size)
        self.timeseries_fluxes_rates(sampling_factor, limits)
        x_ts = self.Epoch_ts
        i_el = 0
        for elem in elements:
            channels = channels_list[i_el]
            self.max_impulsive(elem, sampling_factor, limits, channels)
            for ch in channels:
                y_ts = self.Flux_ts[elem][:, ch]
                dy_ts = self.Flux_unc_ts[elem][:, ch]
                enuc_low = int(1000*self.Enuc_low[elem][ch]) / 1000
                enuc_up = int(1000*self.Enuc_up[elem][ch]) / 1000
                ax[i_el].errorbar(x=x_ts, y=y_ts, yerr=[dy_ts, dy_ts], label=f'{enuc_low} - {enuc_up}'r'MeV n$^{-1}$')
                ax[i_el].set_ylabel(f'{elements_properties(elem)["symbol"]}', fontsize=sub_ylabel_size)
                ax[i_el].set_yscale('log')
                ax[i_el].set_xlim(x_ts[0]-dt.timedelta(hours=0.5), x_ts[len(x_ts)-1]+dt.timedelta(hours=1))
            ax[i_el].legend(loc='upper right', fontsize=legend_size)
            
            #if i_el < n_elem-1: # plot xticks only in the lowest subplot
                #ax[i_el].xaxis.set_ticks_position('bottom')
                #ax[i_el].spines.bottom.set_visible(False)
                
            if plot_max_impulsive == 'yes':
                ax[i_el].plot(self.epoch_max, self.flux_max, '^k')
            if plot_max_impulsive == 'unc':
                ax[i_el].errorbar(self.epoch_max, self.flux_max, yerr=[self.d_flux_max_low, self.d_flux_max_up], fmt='.', color='black')
                i_ch = 0
                for ch in channels:
                    xl = self.epoch_max_low[i_ch]
                    xu = self.epoch_max_up[i_ch]
                    ylu = self.flux_max[i_ch]
                    ax[i_el].plot([xl, xu], [ylu, ylu], color='black')
                    i_ch = i_ch+1
            
            i_el = i_el+1
        
        ax[0].set_title(f'SIS sensor, time series of ions, time resolution: {sampling_factor * self.original_time_resolution[0]/60} min', fontsize=title_size)
        fig.supxlabel(xlabel_text, fontsize=xlabel_size)
        fig.supylabel(r'Intensity [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]', fontsize=ylabel_size)
        plt.show(block=False)
        
    def ratio_fluxes(self, elem_pair, sampling_factor, limits, channels_pair):
        
        # first ion species 
        self.timeseries_fluxes_rates_sum(channels_pair[0], sampling_factor, limits)
        y1 = self.Flux_ts_sum[elem_pair[0]]
        dy1 = self.Flux_unc_ts_sum[elem_pair[0]]
        
        # second ion species
        self.timeseries_fluxes_rates_sum(channels_pair[1], sampling_factor, limits)
        y2 = self.Flux_ts_sum[elem_pair[1]]
        dy2 = self.Flux_unc_ts_sum[elem_pair[1]]
        
        # ratio of the intensities: first ion species / second ion species 
        self.y_ratio = y1 / y2
        self.dy_ratio = 1/y2 * np.sqrt( dy1**2 + (self.y_ratio*dy2)**2 )
        
        # x axis
        self.x_ratio = self.Epoch_ts
        
    def plot_ratio_fluxes(self, elem_pair, sampling_factor, limits, channels_pair, fig_size=fig_size, title_size=title_size, xlabel_text='time [yyyy-mm-dd]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=sub_ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, y_lims='no'):
        
        # create the variables of the ratio abundances
        self.ratio_fluxes(elem_pair, sampling_factor, limits, channels_pair)
        
        # create variables for the energy range of the title
        el_1 = elem_pair[0]
        ch_1  = channels_pair[0]
        enuc1_low = self.Enuc_low[el_1][ch_1[0]]
        enuc1_up = self.Enuc_up[el_1][ch_1[len(ch_1)-1]]
        enuc1_low_label, enuc1_up_label = int(1000*enuc1_low)/1000, int(1000*enuc1_up)/1000
        
        # plot ratio abundances
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=fig_size)
        ax.errorbar(x=self.x_ratio, y=self.y_ratio, yerr=[self.dy_ratio, self.dy_ratio], fmt='.')
        ax.set_title(r'SIS sensor, ratio of intensities [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$], 'f'energy range: {enuc1_low_label} - {enuc1_up_label} MeV/nuc', fontsize=title_size)
        ax.set_xlabel(xlabel_text, fontsize=xlabel_size)
        ax.set_ylabel(f'[{elem_pair[0]}/{elem_pair[1]}]', fontsize=ylabel_size)
        ax.set_yscale('log')
        if y_lims != 'no':
            ax.set_ylim(y_lims)
        ax.set_xlim(self.x_ratio[0]-dt.timedelta(hours=0.5), self.x_ratio[len(self.x_ratio)-1])
        plt.show(block=False)
        
        
    def enuc_to_speed(self, elem, channels): 
        """
        Function that calculate the speed corresponding to each energy channel, for one ion specie. It also calculate the uncertainties considering the width of the channels. The uncerstainties are calculated using error propagation: $\Delta v_{ch} = \sqrt{ ( \frac{ \partial v }{ \partial E } \Delta E )^2 }$, where $\Delta E$ is: 
            - For lower uncertaintie: $\Delta E_{low} = E_{quadratic mean} - E_{lower boundary of the channel}$
            - For upper uncertaintie: $\Delta E_{up} = E_{upper boundary of the channel} - E_{quadratic mean}$
        Of course, there are the lower and the upper speed uncertainties: $\Delta v_{low}$, $\Delta v_{up}$
        - Input variables: Sis object and the id of the ion whose speeds you want to calculate.
        - Output variables: the speed corresponding to each energy channel (calculated with the quadratic mean of the energy of the channel) and their uncertainties (uncertainties are Delta_v, not v-Delta_v). Units: [AU/day]
        """
        elem_prop = elements_properties(elem)
        mass = elem_prop['mass'] #[kg]
        n_nuc = elem_prop['number_of_nucleons'] 

        # energies only of the selected channels
        enuc_ch = self.Enuc[elem][channels]
        enuc_ch_low = self.Enuc_low[elem][channels]
        enuc_ch_up = self.Enuc_up[elem][channels]

        # energy per nucleon in MeV to energy per ion in Joules:
        energy_Joules = n_nuc * enuc_ch * mev_to_J # quadratic average
        energy_Joules_low = n_nuc * enuc_ch_low * mev_to_J # lower boundary
        energy_Joules_up = n_nuc * enuc_ch_up * mev_to_J #upper boundary

        # speed of the element "elem" in Astronomical Units per day (AU/day):
        self.speed_ms = np.sqrt(2 * energy_Joules / mass) #m/s
        self.speed = self.speed_ms * ms_to_AUd #AU/day
        self.d_speed_low = (energy_Joules - energy_Joules_low) / (self.speed_ms * mass) * ms_to_AUd /3 #AU/day
        self.d_speed_up = (energy_Joules_up - energy_Joules) / (self.speed_ms * mass) * ms_to_AUd /3 #AU/day
        
        
    def speed_vs_time_preparation(self, elem, sampling_factor, limits, channels, init_values=[2,2,18595]):
        """
        Prepare the data for the fit of one ion species. 
        """
        
        # load some necessary variables from the methods
        self.max_impulsive(elem, sampling_factor, limits, channels) # load variables for x axis (epochs of the maxima and their uncertainties)
        self.enuc_to_speed(elem, channels) # load variables for y axis (speeds)
        
        # x and y axis, and their uncertainties; for the model
        self.x_odr = self.epoch_max_days - init_values[2] #[days]
        self.dx_odr = (self.d_epoch_max_days_low + self.d_epoch_max_days_up) /2 #[days]
        self.y_odr = self.speed #AU/day
        self.dy_odr = (self.d_speed_low + self.d_speed_up) /2 #AU/day
            

    def plot_speed_vs_time(self, elem, sampling_factor, limits, channels, init_values=[2,2,18595], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=box_size, box_position=box_position):
        """
        Plot the fit of one ion specie "elem". 
        """
        # load and create needed variables
        self.speed_vs_time_preparation(elem, sampling_factor, limits, channels, init_values=[2,2,18595])
        Fit = fit_ODR(self.x_odr, self.y_odr, self.dx_odr, self.dy_odr, init_values=init_values)
        
        # create the figure and the axes
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=fig_size)
        
        # plot measurements (speeds and dates of the maxima intensities) and their uncertainties
        ax.errorbar(x=self.epoch_max, y=self.speed, yerr=[self.d_speed_low, self.d_speed_up], fmt='.b', label='SIS data') # plot x (time [date]), y (speed [AU/day]) and uncertainties of y.
        for i in range(len(self.epoch_max)): # uncertainties of x (time [date])
            xl = self.epoch_max_low[i]
            xu = self.epoch_max_up[i]
            ylu = [self.speed[i], self.speed[i]]
            ax.plot([xl,xu], ylu,'b')
        
        # plot curve of the model
        ax.plot(Fit['xfit'], Fit['yfit'], color='red', label='fit')
        
        # box of annotations
        ax.annotate(Fit['annotation_string'], xy=box_position, xycoords='axes fraction' ,backgroundcolor='w', fontsize=box_size)
        
        # title, labels,...
        ax.set_title(f'SIS sensor, fit of {elements_properties(elem)["symbol"]}, time resolution: {sampling_factor * self.original_time_resolution[0]/60} min.', fontsize=title_size)
        ax.set_xlabel(xlabel_text, fontsize=xlabel_size)
        ax.set_ylabel('speed [AU/day]', fontsize=ylabel_size)
        plt.xticks(fontsize=xticks_size)
        plt.yticks(fontsize=yticks_size)
        plt.legend(fontsize=legend_size)
        plt.show(block=False)
        
        

    def speed_vs_time_total(self, elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, object_ept, init_values_EPT = [2,2,18595]):
        """
        Prepare the data for the fit of several ion species in "elements" list (or tuple). 
        """
        
        x_odr_ions, dx_odr_ions, y_odr_ions, dy_odr_ions = [], [], [], [] 
        x_odrplot_ions_list, x_odrplot_ions_low_list, x_odrplot_ions_up_list = [], [], []
        dy_odrplot_ions_list_low, dy_odrplot_ions_list_up = [], [] 
        self.y_odrplot_ions, self.dy_odrplot_ions_low, self.dy_odrplot_ions_up = {}, {}, {}
        self.xfit_ions, self.yfit_ions, self.s_ions, self.ds_ions, self.t0_ions, self.t0_days_ions, self.dt0_days_ions, self.chi2red_ions, self.annotation_string_ions = {}, {}, {}, {}, {}, {}, {}, {}, {}
        self.x_odrplot_ions, self.x_odrplot_ions_low, self.x_odrplot_ions_up = {}, {}, {} 
        self.y_odrplot_ions, self.dy_odrplot_ions_low, self.dy_odrplot_ions_up = {}, {}, {}
        ii = 0
        for elem in elements: 
            # input variables for each ion species, and load variables from the method "self.speed_vs_time". 
            sampling_factor = sampling_factor_list[ii]
            limits = limits_list[ii]
            channels = channels_list[ii]
            init_values = init_values_list[ii] 
            
            # load some necessary variables from the functions
            self.speed_vs_time_preparation(elem, sampling_factor, limits, channels, init_values=init_values)
            Fit = fit_ODR(self.x_odr, self.y_odr, self.dx_odr, self.dy_odr, init_values=init_values)
            
            # variables and parameters from the model
            self.xfit_ions[elem] = Fit['xfit']
            self.yfit_ions[elem] = Fit['yfit'] 
            self.s_ions[elem] = Fit['s']
            self.ds_ions[elem] = Fit['ds']
            self.t0_ions[elem] = Fit['t0']
            self.t0_days_ions[elem] = Fit['t0_days']
            self.dt0_days_ions[elem] = Fit['dt0_days']
            self.chi2red_ions[elem] = Fit['chi2red']
            self.annotation_string_ions[elem] = Fit['annotation_string']
            
            # variables to do total fit (list to concatenate afterwards)
            x_odr_ions.append(self.x_odr)
            dx_odr_ions.append(self.dx_odr)
            y_odr_ions.append(self.speed)
            dy_odr_ions.append(self.dy_odr)
            
            # variables to plot total fit (list)
            x_odrplot_ions_list.append(self.epoch_max)
            x_odrplot_ions_low_list.append(self.epoch_max_low)
            x_odrplot_ions_up_list.append(self.epoch_max_up)
            dy_odrplot_ions_list_low.append(self.d_speed_low)
            dy_odrplot_ions_list_up.append(self.d_speed_up)
            # variables to plot ion by ion (tuple)
            self.x_odrplot_ions[elem] = self.epoch_max
            self.x_odrplot_ions_low[elem] = self.epoch_max_low
            self.x_odrplot_ions_up[elem] = self.epoch_max_up
            self.y_odrplot_ions[elem] = self.speed
            self.dy_odrplot_ions_low[elem] = self.d_speed_low
            self.dy_odrplot_ions_up[elem] = self.d_speed_up
            
            ii = ii+1
            
        ## EPT
        if object_ept != 'no':
            """
            If you add EPT sensor, you should do this previously:
                1) from EPT_class import EPT
                2) ept = EPT(direction, first_date, last_date, first_date_quiet, last_date_quiet, show_loaded_files)
                3) ept.speed_vs_time_preparation(sampling_factor_EPT, limits_EPT, channels_EPT, init_values_EPT=[2,2,18595])
                4) The input variable "object_ept" of this method "speed_vs_time_total" should be "object_ept = ept" instead of "object_ept = 'no'". Then, the variable "object_ept" is an object (the object of the class "EPT"). 
            """
            # variables to do total fit (list to concatenate afterwards)
            x_odr_ions.append(object_ept.x_odr)
            dx_odr_ions.append(object_ept.dx_odr)
            y_odr_ions.append(object_ept.speed)
            dy_odr_ions.append(object_ept.dy_odr)
            
            # variables to plot total fit (list)
            x_odrplot_ions_list.append(object_ept.epoch_max)
            x_odrplot_ions_low_list.append(object_ept.epoch_max_low)
            x_odrplot_ions_up_list.append(object_ept.epoch_max_up)
            dy_odrplot_ions_list_low.append(object_ept.d_speed_low)
            dy_odrplot_ions_list_up.append(object_ept.d_speed_up)
            # variables to plot ion by ion (tuple)
            self.x_odrplot_ions['H_EPT'] = object_ept.epoch_max
            self.x_odrplot_ions_low['H_EPT'] = object_ept.epoch_max_low
            self.x_odrplot_ions_up['H_EPT'] = object_ept.epoch_max_up
            self.y_odrplot_ions['H_EPT'] = object_ept.speed
            self.dy_odrplot_ions_low['H_EPT'] = object_ept.d_speed_low
            self.dy_odrplot_ions_up['H_EPT'] = object_ept.d_speed_up
            
            
        ## TOTAL FIT:
        
        # x and y axis, and their uncertainties; for the model; of the total fit
        self.x_odr_total = np.concatenate(x_odr_ions)
        self.dx_odr_total = np.concatenate(dx_odr_ions)
        self.y_odr_total = np.concatenate(y_odr_ions)
        self.dy_odr_total = np.concatenate(dy_odr_ions)
        
        # Do the total fit
        Fit = fit_ODR(self.x_odr_total, self.y_odr_total, self.dx_odr_total, self.dy_odr_total, init_values=init_values_total)
        
        # variables and parameters from the model of the total fit
        self.xfit_total = Fit['xfit']
        self.yfit_total = Fit['yfit'] 
        self.s_total = Fit['s']
        self.ds_total = Fit['ds']
        self.t0_total = Fit['t0']
        self.t0_days_total = Fit['t0_days']
        self.dt0_days_total = Fit['dt0_days']
        self.chi2red_total = Fit['chi2red']
        self.annotation_string_total = Fit['annotation_string']
        
        # variables to total plot
        self.x_odrplot_total = np.concatenate(x_odrplot_ions_list)
        self.x_odrplot_total_low = np.concatenate(x_odrplot_ions_low_list)
        self.x_odrplot_total_up = np.concatenate(x_odrplot_ions_up_list)
        self.y_odrplot_total = np.concatenate(y_odr_ions)
        self.dy_odrplot_total_low = np.concatenate(dy_odrplot_ions_list_low)
        self.dy_odrplot_total_up = np.concatenate(dy_odrplot_ions_list_up)
        
        
        
    def plot_speed_vs_time_total(self, object_ept, elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, Plot=['ions', 'total'], init_values_EPT = [2,2,18595], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=box_size, box_position=box_position):
        
        
        self.speed_vs_time_total(elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, object_ept=object_ept)
        
        if object_ept != 'no':
            elements.append('H_EPT')
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=fig_size)
        
        color_list = ['blue', 'green', 'red', 'magenta', 'orange', 'purple', 'yellow']
        iii = 0
        for elem in elements:
            # plot x (times of maximum intensities [date]), y (speeds) and uncertainties of y, for each ion species "elem" 
            if elem != 'H_EPT':
                text_legend = f'{elements_properties(elem)["symbol"]} data'
            else: text_legend = 'H data (EPT)'
            ax.errorbar(x=self.x_odrplot_ions[elem], y=self.y_odrplot_ions[elem], yerr=[self.dy_odrplot_ions_low[elem], self.dy_odrplot_ions_up[elem]], fmt='.', color=color_list[iii], label=text_legend)

            # plot uncertainties of x (time [date]) for each ion species "elem":
            k=0
            for k in range(len(self.x_odrplot_ions_low[elem])): 
                xl = self.x_odrplot_ions_low[elem][k]
                xu = self.x_odrplot_ions_up[elem][k]
                ylu = self.y_odrplot_ions[elem][k]
                ax.plot([xl,xu], [ylu, ylu], color=color_list[iii])
                k=k+1
            
            if 'ions' in Plot: # plot curve of the fit (or model) for each ion species "elem":
                if elem != 'H_EPT':
                    ax.plot(self.xfit_ions[elem], self.yfit_ions[elem], '--', color=color_list[iii], label=f'{elements_properties(elem)["symbol"]} fit') 
                else: 
                    Fit = fit_ODR(object_ept.x_odr, object_ept.y_odr, object_ept.dx_odr, object_ept.dy_odr, init_values=init_values_EPT)
                    ax.plot(Fit['xfit'], Fit['yfit'], '--', color=color_list[iii], label='H fit (EPT)') 
                
            iii = iii+1
        
                
                
        if 'total' in Plot: 
            
            # plot curve of total fit
            ax.plot(self.xfit_total, self.yfit_total, color='black', label='total fit')
            
            # box of annotations (of total fit)
            ax.annotate(self.annotation_string_total, xy=box_position, xycoords='axes fraction' ,backgroundcolor='w', fontsize=box_size)
        
        # title, labels,...
        ax.set_title(f'Velocity dispersion fit', fontsize=title_size)
        ax.set_xlabel(xlabel_text, fontsize=ylabel_size)
        ax.set_ylabel('speed [AU/day]', fontsize=xlabel_size)
        plt.xticks(fontsize=xticks_size)
        plt.yticks(fontsize=yticks_size)
        plt.legend()
        plt.show(block=False)
        
        
    def spectrum(self, limits):
        """
        This method calculates the averaged intensities and rates (and their uncertainties) during a period of time determined by the input variable "limits" (described below) for all the energy channels. It do this for all ion species determined by the input variable "elements" of the "__init__" method (above). But it do this for the period of time determined by "limits" (as we say above), but also substract the quiet time and calculate the spectrum of the quiet time. You should represent these variables against the energy (variable "self.Enuc[elem]").
        
        - Input variables:
        
            - "limits": it is a list or a tuple with 2 items. Each item is a tuple or a list that represent a date in format (yyyy,mm,dd,hh) or (yyyy,mm,dd,hh,min). The second date shoud be larger than the first one because these dates are the limits where you want to restict the time window. Example: "limits=((2022,2,12), (2022,2,6))"
            
        - Output variables:
            
            - "self.Flux_sp": average of fluxes for the all energy channels and during the time window determined by the input variable "limits". 
            - "self.Flux_unc_sp": uncertainties of "self.Flux_sp". 
            - "self.Rate_sp": same as "self.Flux_sp", but with rates.
            - "self.Rate_unc_sp": uncertainties of "self.Rate_sp". 
            
            - "self.Flux_sp_noqt", "self.Flux_unc_sp_noqt", "self.Rate_sp_noqt", "self.Rate_unc_sp_noqt": same as the last 4 variables above but with the quiet time substracted. 
            
            - "self.Flux_sp_quiet", "self.Flux_unc_sp_quiet", "self.Rate_sp_quiet", "self.Rate_unc_sp_quiet": same but only the quiet time, so it is the spectrum of the quiet time. Thus, for example the flux variable: "self.Flux_sp - self.Flux_sp_quiet = self.Flux_sp_noqt". 
        
        """
        
        ind1 = select_date_index(self.Epoch, cutoff_date=limits[0]) #index of the first date where you want to cut
        ind2 = select_date_index(self.Epoch, cutoff_date=limits[1]) #index of the second date where you want to cut
        self.ind1_sp, self.ind2_sp = ind1, ind2
        
        self.Flux_sp_quiet, self.Flux_unc_sp_quiet, self.Rate_sp_quiet, self.Rate_unc_sp_quiet = {},{},{},{} #quiet time
        self.Flux_sp, self.Flux_unc_sp, self.Rate_sp, self.Rate_unc_sp = {},{},{},{} # not substracted quiet time
        self.Flux_sp_noqt, self.Flux_unc_sp_noqt, self.Rate_sp_noqt, self.Rate_unc_sp_noqt = {},{},{},{} # substracted quiet time
        for elem in self.elements:
            
            # Sum over all the channels in the quiet time (Nils said that in spectrum I should sum the counts (rates * time step) or intensities (fluxes) during the selected interval of time, but in time series I should do the average of counts (rates * time step) or intensities (fluxes) for each time step (time resolution). 
            self.Flux_sp_quiet[elem] = self.Flux_quiet[elem].sum(axis=0) 
            self.Flux_unc_sp_quiet[elem] = np.sqrt(((self.Flux_unc_quiet[elem])**2).sum(axis=0))
            self.Rate_sp_quiet[elem] = self.Rate_quiet[elem].sum(axis=0)
            self.Rate_unc_sp_quiet[elem] = np.sqrt(((self.Rate_unc_quiet[elem])**2).sum(axis=0))
            
            # Sum over all the channels in the data without substract the quiet time
            self.Flux_sp[elem] = self.Flux[elem][ind1:ind2+1, :].sum(axis=0) 
            self.Flux_unc_sp[elem] = np.sqrt(((self.Flux_unc[elem][ind1:ind2+1, :])**2).sum(axis=0))
            self.Rate_sp[elem] = self.Rate[elem][ind1:ind2+1, :].sum(axis=0)
            self.Rate_unc_sp[elem] = np.sqrt(((self.Rate_unc[elem][ind1:ind2+1, :])**2).sum(axis=0))

            #Substract quiet time
            self.Flux_sp_noqt[elem] = self.Flux_sp[elem] - self.Flux_sp_quiet[elem]
            self.Flux_unc_sp_noqt[elem] = np.sqrt(self.Flux_unc_sp[elem]**2 + self.Flux_unc_sp_quiet[elem]**2)
            self.Rate_sp_noqt[elem] = self.Rate_sp[elem] - self.Rate_sp_quiet[elem]
            self.Rate_unc_sp_noqt[elem] = np.sqrt(self.Rate_unc_sp[elem]**2 + self.Rate_unc_sp_quiet[elem]**2)
            
            # All these above variables have been checked and they are correct. 
    
    def plot_spectrum(self, elem, limits, pure_nobg_bg=['no','yes','yes']):
        self.spectrum(limits)
        
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=fig_size)
        
        # prepare information for legend
        ## first and last dates of the event (year/month/day) 
        label_event1, label_event2 = f'{limits[0][0]}/{limits[0][1]}/{limits[0][2]}', f'{limits[1][0]}/{limits[1][1]}/{limits[1][2]}'
        ## first minute of the event
        if len(limits[0]) == 4: mm1 = '00' #if limits[0] has not minutes (therefore it's format is (yyyy,mm,dd,hh))
        elif len(limits[0]) == 5: #if limits[0] has minutes (therefore it's format is (yyyy,mm,dd,hh,mm))
            if limits[0][4] > 9: mm1 = f'{limits[0][4]}'
            else: mm1 = f'0{limits[0][4]}' # if the number of minutes is between 1 and 9, add a leading zero (add a zero to the left)
        ## last minute of the event
        if len(limits[1]) == 4: mm2 = '00'
        elif len(limits[1]) == 5: 
            if limits[1][4] > 9: mm2 = f'{limits[1][4]}'
            else: mm2 = f'0{limits[1][4]}'
        ## complete dates and hours of the legend for the event. Format: yyyy/mm/dd hh:mm
        label_event =f'{label_event1} {limits[0][3]}:{mm1} - {label_event2} {limits[1][3]}:{mm2}'
        ## dates of the quiet time (for the background). TODO: hours and minutes have not been added, in fact quiet time in the "spectrum" function is not cut by "limits", but it is directly from the __init__ method using the data from complete days (the tuples with the dates of quiet time are "self.first_date_qt" and "self.last_date_qt")
        label_background = f'{self.first_date_qt[0]}/{self.first_date_qt[1]}/{self.first_date_qt[2]} - {self.last_date_qt[0]}/{self.last_date_qt[1]}/{self.last_date_qt[2]}'
        
        # add legend
        if pure_nobg_bg[0] == 'yes':
            ax.errorbar(x=self.Enuc[elem], y=self.Flux_sp[elem], yerr=[self.Flux_unc_sp[elem],self.Flux_unc_sp[elem]], fmt='^', color='red', label=f'event ({label_event})')
        if pure_nobg_bg[1] == 'yes':
            ax.errorbar(x=self.Enuc[elem], y=self.Flux_sp_noqt[elem], yerr=[self.Flux_unc_sp_noqt[elem],self.Flux_unc_sp_noqt[elem]], fmt='^', color='blue', label=f'event-background ({label_event})')
        if pure_nobg_bg[2] == 'yes':
            ax.errorbar(x=self.Enuc[elem], y=self.Flux_sp_quiet[elem], yerr=[self.Flux_unc_sp_quiet[elem],self.Flux_unc_sp_quiet[elem]], fmt='^', color='green', label=f'background ({label_background})')
        
        # title, labels, x scale, y scale,...
        ax.set_title(f'SIS sensor, {elements_properties(elem)["symbol"]}, spectrum of intensities')
        ax.set_xlabel(r'Energy per nucleon [MeV n$^{-1}$]', fontsize=xlabel_size)
        ax.set_ylabel(r'Intensity [$s^{-1} cm^{-2} sr^{-1} MeV^{-1}$]', fontsize=ylabel_size)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.legend(fontsize=legend_size)
        plt.show(block=False)
    
    
# Falta:
    # Comprobar que funcionan bien todas las funciones de speed vs time
    # que en la última función, la de total fit plot, te devuelva una tabla o algo (o un print) con los resultados (s, t0, ds,...) de cada ion species.
    # cambiar EPT y añadirlo aquí
    # añadir las funciones que faltan (como las de fit spectrum, time series spectra slope,... que están en el Jupyter Notebook SLOPES_SPECTRA__a)
    # comentar bien y poner descripción a todas las funciones
