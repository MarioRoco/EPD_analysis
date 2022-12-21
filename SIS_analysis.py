from SIS_class import SIS
from SIS_auxfuncs import *

# To disable the warnings
import warnings
warnings.filterwarnings("ignore")

quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
quietday_end=(2020,11,12)#inclusive, for ACR background estimation
eventday_start, eventday_end = (2022,3,5), (2022,3,10)

limits=[(2022,3,5,14,7), (2022,3,10,22)]

sis = SIS(direction='sun', elements=['H','He3','He4','O','Fe'], first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

#sis.plot_timeseries_fluxes_rates(elements=['H','He4'], sampling_factor=2*60, limits=limits, channels_list=[[0,1,2,3,4,5,6,7,8,9,10],[1,2,3,4,5,6,7,8,9,10]], plot_max_impulsive='yes', fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=14, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=7)

# EPT
from EPT_class import EPT

init_values_EPT = [2,2,18595]

ept = EPT(direction='sun', first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

ept.speed_vs_time_preparation(sampling_factor=60*20, limits=[(2022,3,5,18), (2022,3,6,12)], channels=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,24,27,29,30], init_values=init_values_EPT)

# SIS ion species


ch_H = [0,1,2,3,4,5]
sf_H = 60
lim_H = [(2022,3,5,9), (2022,3,6,20)]
initv_H = [2, 2, 18595]

ch_He4 = [1,2, 3, 4, 5]
sf_He4 = 60
lim_He4 = [(2022,3,5,9), (2022,3,6,20)]
initv_He4 = [2, 2, 18595]

ch_He3 = [0,1,2,3,4,5,6,7]
sf_He3 = 60
lim_He3 = [(2022,3,5,9), (2022,3,6,20)]
initv_He3 = [2, 2, 18595]

ch_Fe = [0,1,2,3,4,5,6,7]
sf_Fe = 60
lim_Fe = [(2022,3,5,9), (2022,3,6,20)]
initv_Fe = [2, 2, 18595]

elements = ['H', 'He4', 'He3', 'Fe']
sampling_factor_list = [sf_H, sf_He4, sf_He3, sf_Fe]
limits_list = [lim_H, lim_He4, lim_He3, lim_Fe]
channels_list = [ch_H, ch_He4, ch_He3, ch_Fe]
init_values_list = [initv_H, initv_He4, initv_He3, initv_Fe]
init_values_total = [2, 2, 18595]

sis.plot_speed_vs_time_total(object_ept=ept, elements=elements, sampling_factor_list=sampling_factor_list, limits_list=limits_list, channels_list=channels_list, init_values_list=init_values_list, init_values_total=init_values_total, Plot=['ions', 'total'], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=box_size, box_position=box_position)

print(sis.s_ions)
print(sis.ds_ions)
print(sis.t0_ions)
print(sis.chi2red_ions)












#################################################################################### check time series
check_timeseries = 'no'

if check_timeseries == 'yes':
	quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
	quietday_end=(2020,11,12)#inclusive, for ACR background estimation
	eventday_start, eventday_end = (2022,3,5), (2022,3,10)
	sis = SIS(direction='sun', elements=['H','He3','He4','O','Fe'], first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

	sf=98
	sis.timeseries_fluxes_rates(sampling_factor=sf, limits=[(2022,3,5,14,7), (2022,3,10,22)])

	elem='Fe'
	n=34
	ch=4
	print(sum(sis.Flux[elem][sis.ind1_ts+n*sf:sis.ind1_ts+n*sf+sf, ch])/sf) 
	print(sis.Flux_ts[elem][n,ch]) 

	print(np.sqrt(sum((sis.Flux_unc[elem][sis.ind1_ts+n*sf:sis.ind1_ts+n*sf+sf, ch])**2))/sf) 
	print(sis.Flux_unc_ts[elem][n,ch]) 


######################################################################## EXAMPLES:
use_timeseries = 'no'
use_timeseries_ratios = 'no'
use_spectrum = 'no'
use_speed_vs_time = 'no'
use_speed_vs_time_someions = 'no'

if use_timeseries=='yes' or use_timeseries_ratios=='yes' or use_spectrum=='yes' or use_speed_vs_time=='yes' or use_speed_vs_time_someions=='yes':
	
	# variables we need later:
	quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
	quietday_end=(2020,11,12)#inclusive, for ACR background estimation
	eventday_start, eventday_end = (2022,3,5), (2022,3,10)
	direction = 'sun'
	ion_species = ['H','He3','He4','O','Fe', 'C', 'Mg']
	limits=[(2022,3,5,9), (2022,3,6,20,46)]
	
	ch_H = [1,2,3,4,5,6]
	ch_He4 = [2, 3, 4, 5, 6]
	ch_He3 = [1,2,3,4,5,6,7,8]
	ch_Fe = [1,2,3,4,5,6,7,8]

	# load SIS class
	sis = SIS(direction=direction, elements=ion_species, first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')
	
	# Time series 
if use_timeseries == 'yes':
	sis.plot_timeseries_fluxes_rates(elements=['H', 'He4', 'He3', 'Fe'], sampling_factor=60, limits=limits, channels_list=[ch_H, ch_He4, ch_He3, ch_Fe], plot_max_impulsive='yes', fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=14, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=7)

	# time series of ratios
if use_timeseries_ratios == 'yes':
	sis.plot_ratio_fluxes(elem_pair=['H', 'He4'], sampling_factor=120, limits=limits, channels_pair=[ch_H, ch_He4], fig_size=(15,2), title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=sub_ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, y_lims='no')
	
	# spectrum
if use_spectrum == 'yes':
	sis.plot_spectrum(elem='H', limits=limits, pure_nobg_bg=['yes','yes','yes'])


	# fit speed vs time for one ion species
if use_speed_vs_time == 'yes':
	sis.plot_speed_vs_time(elem='H', sampling_factor=60, limits=limits, channels=[0,1,2,3,4,5], init_values=[2,2,18595], fig_size=fig_size, title_size=title_size, xlabel_text='time [dd hh:mm]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=14, box_position=(0.3,0.8))

if use_speed_vs_time_someions == 'yes':
	# fit speed vs time for several ion species from SIS
	ch_H = [1,2,3,4,5,6]
	sf_H = 60
	lim_H = [(2022,3,5,9), (2022,3,6,20)]
	initv_H = [2, 2, 18595]

	ch_He4 = [2, 3, 4, 5, 6]
	sf_He4 = 60
	lim_He4 = [(2022,3,5,9), (2022,3,6,20)]
	initv_He4 = [2, 2, 18595]

	ch_He3 = [1,2,3,4,5,6,7,8]
	sf_He3 = 60
	lim_He3 = [(2022,3,5,9), (2022,3,6,20)]
	initv_He3 = [2, 2, 18595]

	ch_Fe = [1,2,3,4,5,6,7,8]
	sf_Fe = 60
	lim_Fe = [(2022,3,5,9), (2022,3,6,20)]
	initv_Fe = [2, 2, 18595]


	elements = ['H', 'He4', 'He3', 'Fe']
	sampling_factor_list = [sf_H, sf_He4, sf_He3, sf_Fe]
	limits_list = [lim_H, lim_He4, lim_He3, lim_Fe]
	channels_list = [ch_H, ch_He4, ch_He3, ch_Fe]
	init_values_list = [initv_H, initv_He4, initv_He3, initv_Fe]
	init_values_total = [2, 2, 18595]
	
	sis.plot_speed_vs_time_total(elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, Plot=['ions', 'total'], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=15, box_position=(0.3, 0.78))
	sis.plot_speed_vs_time_total(elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, Plot=['ions'], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=15, box_position=(0.3, 0.78))
	sis.plot_speed_vs_time_total(elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, Plot=['total'], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=15, box_position=(0.3, 0.78))
	sis.plot_speed_vs_time_total(elements, sampling_factor_list, limits_list, channels_list, init_values_list, init_values_total, Plot=[], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=15, box_position=(0.3, 0.78))

	
