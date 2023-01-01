from utils import *
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_out_file = "gocart.nc"
print (wrf_dir+wrf_out_file)
nc_fid = nc.MFDataset(wrf_dir+wrf_out_file)
times =nc_fid.variables['Times'][:]
#WRF-Chem's flux (kg/m2/sec)
flux= nc_fid.variables['EDUST1'][:,0]+\
	  nc_fid.variables['EDUST2'][:,0]+\
	  nc_fid.variables['EDUST3'][:,0]+\
	  nc_fid.variables['EDUST4'][:,0]+\
	  nc_fid.variables['EDUST5'][:,0]
nc_fid.close()


k=len(times)
fig, axes = plt.subplots(1, k,figsize=(k*6,6))

print ("processing " +wrf_dir+wrf_out_file)
for time_idx in np.arange(0, k):
	#because emissions are accumulated here (cumulative, kg/m2)
	#make them average (during one hour) flux
	#flux=(emissions-prev_emissions)/3600.0 # (kg/m2/sec) it will be a little different from what calculated in gocart_python.py
	
	total_emission_flux=np.sum(surface*flux[time_idx]) #(kg/sec)

	m = Basemap(**basemap_params,ax=axes[time_idx])
	x, y = m(xlon,xlat)
	decorateMap(m)	

	date_time_obj = datetime.strptime(str(b"".join(times[time_idx])), "b'%Y-%m-%d_%H:%M:%S'")
	axes[time_idx].set_title(date_time_obj.strftime("%d %B, %H%M UTC")+"\n Instant dust flux: "+"{:0.1f}".format(total_emission_flux)+" ($kg\ sec^{-1}$)")

	cs = m.pcolormesh(x, y, flux[time_idx], cmap=colmap, norm=ai_norm)


cbar = fig.colorbar(cs,ax=axes.ravel().tolist(),orientation="horizontal",extend='max',format='%.0e')
cbar.set_label('Orignl. GOCART Dust emissions, '+units)
plt.savefig("gocart_orgl_inst_flux.png",bbox_inches="tight",dpi=300)