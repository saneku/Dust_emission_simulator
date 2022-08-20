from utils import *
import netCDF4 as nc
import os
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_out_file="/wrfout_d01_2016-06-24_00:00:00_gocart"
print wrf_dir+wrf_out_file
nc_fid = nc.MFDataset(wrf_dir+wrf_out_file)
times =nc_fid.variables['Times'][:]

prev_emissions=0

print "processing " +wrf_dir+wrf_out_file
for time_idx in range(1,len(times),1):
	print ''.join(times[time_idx])
	date_time_obj = datetime.strptime(''.join(times[time_idx]), '%Y-%m-%d_%H:%M:%S')

	emissions=nc_fid.variables['EDUST1'][time_idx,0]+nc_fid.variables['EDUST2'][time_idx,0]+nc_fid.variables['EDUST3'][time_idx,0]+nc_fid.variables['EDUST4'][time_idx,0]+nc_fid.variables['EDUST5'][time_idx,0]

	#because emissions are summed here (cumulative, kg/m2)
	#make them flux
	flux=(emissions-prev_emissions)/3600.0
	total_emission_flux=np.sum(surface*flux)
	prev_emissions=emissions

	fig = plt.figure(figsize=(12,12))
	ash_map = Basemap(**basemap_params)
	x, y = ash_map(xlon,xlat)

	decorateMap(ash_map)	
	#plot_cities(ash_map)

	plt.title(date_time_obj.strftime("%d %B, %H:%M %p")+"\n Total emission flux: "+"{:0.1f}".format(total_emission_flux)+" ($kg\ sec^{-1}$)")

	ai_norm = colors.BoundaryNorm(np.logspace(-12.0, -7.0, num=11), ncview_colormap_short.N, clip=False)
	cs=ash_map.pcolormesh(x,y,flux,cmap=ncview_colormap_short, norm=ai_norm)
	
	cbar = fig.colorbar(cs,orientation='horizontal',extend='max')
	cbar.set_label('Original Dust emissions, '+units)#,fontsize=CB_LABEL_TEXT_SIZE)

	plt.savefig("ornlg_"+str(time_idx)+".png",bbox_inches="tight")
	####################################################################################

nc_fid.close()
