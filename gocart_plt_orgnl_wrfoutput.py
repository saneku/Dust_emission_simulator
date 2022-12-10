from utils import *
import netCDF4 as nc
import os
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_dir="./"
wrf_out_file="wrfout_d01_2016-06-24_00:00:00_gocart"
print wrf_dir+wrf_out_file
nc_fid = nc.MFDataset(wrf_dir+wrf_out_file)
times =nc_fid.variables['Times'][:]

print "processing " +wrf_dir+wrf_out_file
for time_idx in range(1,len(times),1):	
	date_time_obj = datetime.strptime(''.join(times[time_idx]), '%Y-%m-%d_%H:%M:%S')

	flux= nc_fid.variables['EDUST1'][time_idx,0]+\
		  nc_fid.variables['EDUST2'][time_idx,0]+\
		  nc_fid.variables['EDUST3'][time_idx,0]+\
		  nc_fid.variables['EDUST4'][time_idx,0]+\
		  nc_fid.variables['EDUST5'][time_idx,0]

	#because emissions are accumulated here (cumulative, kg/m2)
	#make them average (during one hour) flux
	#flux=(emissions-prev_emissions)/3600.0 # (kg/m2/sec) it will be a little different from what calculated in gocart_python.py
	
	#flux (kg/m2/sec)
	total_emission_flux=np.sum(surface*flux) #(kg/sec)

	fig = plt.figure(figsize=(12,12))
	ash_map = Basemap(**basemap_params)
	x, y = ash_map(xlon,xlat)
	decorateMap(ash_map)	

	plt.title(date_time_obj.strftime("%d %B, %H:%M %p")+"\n Integrated Instant. emission flux: "+"{:0.1f}".format(total_emission_flux)+" ($kg\ sec^{-1}$)")

	cs=ash_map.pcolormesh(x,y,flux,cmap=ncview_colormap_short, norm=ai_norm)
	cbar = fig.colorbar(cs,orientation='horizontal',extend='max',format='%.0e')
	cbar.set_label('Orignl. Instant. GOCART Dust emissions, '+units)

	plt.savefig("gocart_orgl_inst_flux_"+str(time_idx)+".png",bbox_inches="tight")

	print (''.join(times[time_idx]),total_emission_flux)
	####################################################################################

nc_fid.close()
