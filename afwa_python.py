#https://ldas.gsfc.nasa.gov/gldas/soils

from utils import *
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import os

from afwa_source_dust import afwa_source_dust
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_dir="./"
wrf_out_file="/wrfout_d01_2016-06-24_00:00:00_afwa"
print (wrf_dir+wrf_out_file)
nc_fid = nc.MFDataset(wrf_dir+wrf_out_file)
times =nc_fid.variables['Times'][:]
xland=nc_fid.variables['XLAND'][:]

ustar=nc_fid.variables['UST'][:]
smois = nc_fid.variables["SMOIS"][:, 0, :]  # soil moisture of first level
isltyp = nc_fid.variables["ISLTYP"][:]
snowh = nc_fid.variables["SNOWH"][:]
znt = nc_fid.variables["ZNT"][:]

airden = 1.0 / nc_fid.variables["ALT"][:, 0, :]
erodtot = nc_fid.variables["EROD"][:] # Ginoux DSR, there is a DRI also

massfrac=np.zeros(shape=(3,ny,nx))  #MASSFRAC  Fraction of mass in each of 3 soil classes    (clay, silt, sand)
massfrac[0][:]=nc_fid.variables['CLAYFRAC'][0,:]
massfrac[1][:]=1-(nc_fid.variables['CLAYFRAC'][0,:]+nc_fid.variables['SANDFRAC'][0,:])
massfrac[2][:]=nc_fid.variables['SANDFRAC'][0,:]
nc_fid.close()

#print drylimit.shape, volsm.shape, erodtot.shape, ustar.shape, ilwi.shape, gravsm.shape, airden.shape, massfrac.shape  #(15, 44, 44)

flux=np.zeros(shape=(ny,nx))

print "processing " +wrf_dir+wrf_out_file
for time_idx in range(1,len(times),1):

	flux,u_ts,u_tres=afwa_source_dust(nx,ny,ustar[time_idx], massfrac,erodtot[time_idx], isltyp[time_idx], smois[time_idx], airden[time_idx],xland[time_idx],znt[time_idx],snowh[time_idx])
	#flux (kg/m2/sec)
	total_emission_flux=np.sum(surface*flux)

	fig = plt.figure(figsize=(12,12))
	ash_map = Basemap(**basemap_params)
	x, y = ash_map(xlon,xlat)
	decorateMap(ash_map)	

	date_time_obj = datetime.strptime(''.join(times[time_idx]), '%Y-%m-%d_%H:%M:%S')
	plt.title(date_time_obj.strftime("%d %B, %H:%M %p")+"\n Total emission flux: "+"{:0.1f}".format(total_emission_flux)+" ($kg\ sec^{-1}$)")

	cs=ash_map.pcolormesh(x,y,flux,cmap=ncview_colormap_short, norm=ai_norm)
	cbar = fig.colorbar(cs,orientation='horizontal',extend='max',format='%.0e')
	cbar.set_label('Instant. GOCART-AFWA Dust emissions, '+units)

	plt.savefig("afwa_inst_flux_"+str(time_idx)+".png",bbox_inches="tight")

	print("".join(times[time_idx]),total_emission_flux)
	####################################################################################
	
	'''
	############################
	DIAGNOSTIC
	#plot u_ts_XXX's calculated by this script
	k=u_ts.shape[0]
	fig, axs = plt.subplots(1, k,figsize=(k*5,5))
	for l in np.arange(0,k):
		cs=axs[l].pcolormesh(u_ts[l],cmap=cm.get_cmap('rainbow', 21),norm=uts_norm)
		axs[l].set_title('u_ts_'+str(l))
	fig.colorbar(cs, ax=axs.ravel().tolist(),orientation='horizontal')
	plt.savefig("uts_"+str(time_idx)+".png")
	
	#plot u_tres_XXX's calculated by this script
	k=u_tres.shape[0]
	fig, axs = plt.subplots(1, k,figsize=(k*5,5))
	for l in np.arange(0,k):
		cs=axs[l].pcolormesh(u_tres[l],cmap=cm.get_cmap('rainbow', 21),norm=uts_norm)
		axs[l].set_title('u_tres_'+str(l))
	fig.colorbar(cs, ax=axs.ravel().tolist(),orientation='horizontal')
	plt.savefig("u_tres_"+str(time_idx)+".png")
	os.system('convert '+"u_tres_"+str(time_idx)+".png "+"uts_"+str(time_idx)+".png"+" -append "+"u_"+str(time_idx)+".png; rm "+"uts_"+str(time_idx)+".png "+"u_tres_"+str(time_idx)+".png")
	'''
	############################


'''
print ustar.shape, massfrac.shape, erodtot.shape, ilwi.shape, gravsm.shape, volsm.shape, airden.shape, drylimit.shape
#exit()
x=np.arange(0,ny)
y=np.arange(0,ny)
fig, axs = plt.subplots(2, 5,figsize=(25,11))

cs=axs[0, 0].pcolormesh(ustar[time_idx],vmin=0, vmax=1.3)
axs[0, 0].set_title('ustar[time_idx]')
fig.colorbar(cs,ax=axs[0, 0],orientation='horizontal')

cs=axs[0, 1].pcolormesh(massfrac[0])
axs[0, 1].set_title('massfrac[0] clay')
fig.colorbar(cs,ax=axs[0, 1],orientation='horizontal')

cs=axs[0, 2].pcolormesh(massfrac[1])
axs[0, 2].set_title('massfrac[1] silt')
fig.colorbar(cs,ax=axs[0, 2],orientation='horizontal')

cs=axs[0, 3].pcolormesh(massfrac[2])
axs[0, 3].set_title('massfrac[2] sand')
fig.colorbar(cs,ax=axs[0, 3],orientation='horizontal')

cs=axs[0, 4].pcolormesh(erodtot[time_idx])
axs[0, 4].set_title('erodtot[time_idx]')
fig.colorbar(cs,ax=axs[0, 4],orientation='horizontal')
	
cs=axs[1, 0].pcolormesh(ilwi[time_idx])
axs[1, 0].set_title('ilwi[time_idx]')
fig.colorbar(cs,ax=axs[1, 0],orientation='horizontal')

cs=axs[1, 1].pcolormesh(gravsm[time_idx])
axs[1, 1].set_title('gravsm[time_idx]')
fig.colorbar(cs,ax=axs[1, 1],orientation='horizontal')

cs=axs[1, 2].pcolormesh(volsm[time_idx])
axs[1, 2].set_title('volsm[time_idx]')
fig.colorbar(cs,ax=axs[1, 2],orientation='horizontal')

cs=axs[1, 3].pcolormesh(airden[time_idx],vmin=0, vmax=1.2)
axs[1, 3].set_title('airden[time_idx]')
fig.colorbar(cs,ax=axs[1, 3],orientation='horizontal')

cs=axs[1, 4].pcolormesh(drylimit[time_idx])
axs[1, 4].set_title('drylimit[time_idx]')
fig.colorbar(cs,ax=axs[1, 4],orientation='horizontal')

plt.savefig(str(time_idx)+"_all.png")

os.system('convert '+str(time_idx)+"_all.png"+" "+str(time_idx)+".png"+" +append ./"+str(time_idx)+".png; rm "+str(time_idx)+"_all.png")
'''