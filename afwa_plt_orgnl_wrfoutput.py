from utils import *
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_out_file="afwa.nc"
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

	#flux (kg/m2/sec)
	total_emission_flux=np.sum(surface*flux[time_idx]) #(kg/sec)

	m = Basemap(**basemap_params,ax=axes[time_idx])
	x, y = m(xlon,xlat)
	decorateMap(m)	

	date_time_obj = datetime.strptime(str(b"".join(times[time_idx])), "b'%Y-%m-%d_%H:%M:%S'")
	axes[time_idx].set_title(date_time_obj.strftime("%d %B, %H%M UTC")+"\n Instant dust flux: "+"{:0.1f}".format(total_emission_flux)+" ($kg\ sec^{-1}$)")

	cs = m.pcolormesh(x, y, flux[time_idx], cmap=colmap, norm=ai_norm)
	
cbar = fig.colorbar(cs,ax=axes.ravel().tolist(),orientation="horizontal",extend='max',format='%.0e')
cbar.set_label('Orignl. AFWA Dust emissions, '+units)
plt.savefig("afwa_ornlg_inst_flux.png",bbox_inches="tight",dpi=300)



'''
#this works only for AFWA wrfoutput, when wrf_out_file="/wrfout_d01_2016-06-24_00:00:00_awfa"
#plot from WRFOUTPUT	
k=9
fig, axs = plt.subplots(1, k,figsize=(k*5,5))
for l in np.arange(0,k):
	u_ts_orgnl=nc_fid.variables['TOT_DUST'][time_idx,l,:]
	cs=axs[l].pcolormesh(u_ts_orgnl,cmap=cm.get_cmap('rainbow', 21),norm=uts_norm)
	axs[l].set_title('WRFOUT u_ts_'+str(l))
fig.colorbar(cs, ax=axs.ravel().tolist(),orientation='horizontal')
plt.savefig("orgnl_uts_"+str(time_idx)+".png")


#plot u_tres_XXX from WRFOUTPUT
k=9
fig, axs = plt.subplots(1, k,figsize=(k*5,5))
index=0
for l in np.arange(9,18):
	u_tres_orgnl=nc_fid.variables['TOT_DUST'][time_idx,l,:]
	cs=axs[index].pcolormesh(u_tres_orgnl,cmap=cm.get_cmap('rainbow', 21),norm=uts_norm)
	axs[index].set_title('WRFOUT u_tres_'+str(index))
	index=index+1
fig.colorbar(cs, ax=axs.ravel().tolist(),orientation='horizontal')
plt.savefig("orgnl_utres_"+str(time_idx)+".png")
os.system('convert '+"orgnl_utres_"+str(time_idx)+".png "+"orgnl_uts_"+str(time_idx)+".png"+" -append "+"orgnl_u_"+str(time_idx)+".png; rm "+"orgnl_utres_"+str(time_idx)+".png "+"orgnl_uts_"+str(time_idx)+".png")
'''