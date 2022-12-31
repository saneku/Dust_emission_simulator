# https://ldas.gsfc.nasa.gov/gldas/soils

from utils import *
import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import os

from gocart_source_dust import gocart_source_dust
from mpl_toolkits.basemap import Basemap
from datetime import datetime

wrf_dir="./data/"
wrf_out_file = "gocart.nc"
print(wrf_dir + wrf_out_file)
nc_fid = nc.MFDataset(wrf_dir + wrf_out_file)
times = nc_fid.variables["Times"][:]
xland = nc_fid.variables["XLAND"][0, :]

airden = 1.0 / nc_fid.variables["ALT"][:, 0, :]
u10 = nc_fid.variables["U10"][:]
v10 = nc_fid.variables["V10"][:]
w10m = np.sqrt(u10 * u10 + v10 * v10)
isltyp = nc_fid.variables["ISLTYP"][:]
smois = nc_fid.variables["SMOIS"][:, 0, :]  # soil moisture of first level
erod = nc_fid.variables["EROD"][:]
nc_fid.close()

flux = np.zeros(shape=(ny, nx))



print("processing " + wrf_dir + wrf_out_file)
for time_idx in range(1, len(times), 1):
	flux, u_ts, u_tres = gocart_source_dust(nx, ny, w10m[time_idx], isltyp[time_idx], smois[time_idx], erod[time_idx], airden[time_idx], xland)
	#flux (kg/m2/sec)
	total_emission_flux = np.sum(surface*flux) #(kg/sec)

	fig = plt.figure(figsize=(8, 8))
	ash_map = Basemap(**basemap_params)
	x, y = ash_map(xlon, xlat)
	decorateMap(ash_map)

	date_time_obj = datetime.strptime(str(b"".join(times[time_idx])), "b'%Y-%m-%d_%H:%M:%S'")
	plt.title(date_time_obj.strftime("%d %B, %H:%M %p") + "\n Integrated Instant. emission flux: " + "{:0.1f}".format(total_emission_flux) + " ($kg\ sec^{-1}$)")

	cs = ash_map.pcolormesh(x, y, flux, cmap=colmap, norm=ai_norm)
	cbar = fig.colorbar(cs, orientation="horizontal",extend='max',format='%.0e')
	cbar.set_label("Instant. GOCART Dust emissions, " + units)

	plt.savefig("gocart_inst_flux_" + str(time_idx) + ".png", bbox_inches="tight")

	print(date_time_obj.strftime("%d.%m.%Y, %H:%M:%S"),total_emission_flux)

	"""
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
	"""


"""
print ustar.shape, massfrac.shape, erodtot.shape, ilwi.shape, gravsm.shape, volsm.shape, airden.shape, drylimit.shape
#exit()
x=np.arange(0,ny)
y=np.arange(0,ny)
fig, axs = plt.subplots(2, 5,figsize=(25,11))

cs=axs[0, 0].pcolormesh(w10m[time_idx],vmin=0, vmax=1.3)
axs[0, 0].set_title('w10m[time_idx]')
fig.colorbar(cs,ax=axs[0, 0],orientation='horizontal')

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
"""