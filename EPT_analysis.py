from EPT_class import EPT
from SIS_auxfuncs import *

# To disable the warnings
import warnings
warnings.filterwarnings("ignore")


quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
quietday_end=(2020,11,12)#inclusive, for ACR background estimation
eventday_start, eventday_end = (2022,3,5), (2022,3,10)


#ch_Hept = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,27,29,30]
sf_Hept = 60*20
lim_Hept = [(2022,3,5,18), (2022,3,6,12)]
initv_Hept = [2, 2, 18595]
ch_Hept = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,24,27,29,30]


ept = EPT(direction='sun', first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

ept.plot_timeseries_fluxes_rates(sampling_factor=sf_Hept, limits=lim_Hept, channels=ch_Hept, plot_max_impulsive='yes', fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=sub_ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=4) 


ept.plot_speed_vs_time(sampling_factor=sf_Hept, limits=lim_Hept, channels=ch_Hept, init_values=[2,2,18595], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=box_size, box_position=box_position)


















#################################################################################### check time series
check_timeseries = 'no'

if check_timeseries == 'yes':
	quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
	quietday_end=(2020,11,12)#inclusive, for ACR background estimation
	eventday_start, eventday_end = (2022,3,5), (2022,3,10)
	ept = EPT(direction='sun', first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

	sf=60*20
	ept.timeseries_fluxes_rates(sampling_factor=sf, limits=[(2022,3,5,14,7), (2022,3,10,22)])

	n=34
	ch=4
	print(sum(ept.Flux[ept.ind1_ts+n*sf:ept.ind1_ts+n*sf+sf, ch])/sf) 
	print(ept.Flux_ts[n,ch]) 

	print(np.sqrt(sum((ept.Flux_unc[ept.ind1_ts+n*sf:ept.ind1_ts+n*sf+sf, ch])**2))/sf) 
	print(ept.Flux_unc_ts[n,ch]) 



########################################################################## Examples
use_timeseries = 'no'
use_spectrum = 'no'
use_speed_vs_time = 'no'

if use_timeseries=='yes' or use_spectrum=='yes' or use_speed_vs_time=='yes':
	quietday_start=(2020,11,7)#inclusive, for ACR background estimation 
	quietday_end=(2020,11,12)#inclusive, for ACR background estimation
	eventday_start, eventday_end = (2022,3,5), (2022,3,10)
	direction = 'sun'
	limits = [(2022,3,5,9), (2022,3,6,20,46)]
	channels = [0,4,9,12,18,20,60]

	ept = EPT(direction=direction, first_date=eventday_start, last_date=eventday_end, first_date_quiet=quietday_start, last_date_quiet=quietday_end, show_loaded_files='no')

if use_timeseries=='yes':
	ept.plot_timeseries_fluxes_rates(sampling_factor=3600*1, limits=limits, channels=channels, plot_max_impulsive='yes', fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, sub_ylabel_size=sub_ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size)

if use_spectrum == 'yes':
	ept.plot_spectrum(limits=limits, pure_nobg_bg=['no','yes','yes'])

if use_speed_vs_time == 'yes':
	ept.plot_speed_vs_time(sampling_factor=3600, limits=limits, channels=channels, init_values=[2,2,18595], fig_size=fig_size, title_size=title_size, xlabel_text='time [mm-dd hh]', xlabel_size=xlabel_size, ylabel_size=ylabel_size, xticks_size=xticks_size, yticks_size=yticks_size, legend_size=legend_size, box_size=box_size, box_position=box_position)

