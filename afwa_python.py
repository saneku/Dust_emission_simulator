#https://ldas.gsfc.nasa.gov/gldas/soils

import netCDF4 as nc
import numpy as np
from matplotlib import pyplot as plt
import os
from netCDF4 import Dataset
from afwa_source_dust import * 
from mpl_toolkits.basemap import Basemap
from datetime import datetime
import matplotlib.colors as colors

wrf_dir='/scratch/ukhova/WRF/run_airquality_oneMonth_100km'
wrf_input_file='wrfinput_d01'
wrf_bdy_file='wrfbdy_d01'
wrf_out_file='/wrfout_d01_2016-06-24_00:00:00'

wrfinput=Dataset(wrf_dir+"/"+wrf_input_file,'r')
xlon=wrfinput.variables['XLONG'][0,:]
xlat=wrfinput.variables['XLAT'][0,:]

dy=wrfinput.getncattr('DY')
dx=wrfinput.getncattr('DX')

MAPFAC_MX=wrfinput.variables['MAPFAC_MX'][0,:]
MAPFAC_MY=wrfinput.variables['MAPFAC_MY'][0,:]

surface=(dx/MAPFAC_MX)*(dy/MAPFAC_MY)       #surface in m2
wrfinput.close()

wrfbddy = Dataset(wrf_dir+"/"+wrf_bdy_file,'r')
nx=len(wrfbddy.dimensions['west_east'])
ny=len(wrfbddy.dimensions['south_north'])
wrfbddy.close()

###############
import ncview_colormap
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import cm

c_map_short = ncview_colormap.default._resample(9)
#print c_map_short.N
c_map = plt.cm.get_cmap(ncview_colormap.default) # select the desired cmap
colorlist=list()
colorlist.append("#ffffff")
for c in np.linspace(0,1,c_map_short.N):
    rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
    clr=colors.rgb2hex(rgba) #convert to hex
    colorlist.append(str(clr)) # create a list of these colors

ncview_colormap_short = LinearSegmentedColormap.from_list('cmap_name', colorlist, N=10)
###############

from matplotlib import cm
ai_norm = colors.BoundaryNorm(np.logspace(-12.0, -7.0, num=11), ncview_colormap_short.N, clip=True)
uts_norm = colors.BoundaryNorm(np.logspace(-1.0, 1.1, num=21), cm.get_cmap('rainbow', 21).N, clip=True)

#####################
massfrac=np.zeros(shape=(3,ny,nx))  #MASSFRAC  Fraction of mass in each of 3 soil classes    (clay, silt, sand)
emissions=np.zeros(shape=(ny,nx))

print wrf_dir+wrf_out_file
nc_fid = nc.MFDataset(wrf_dir+wrf_out_file)
times =nc_fid.variables['Times'][:]
xland=nc_fid.variables['XLAND'][0,:]

drylimit=nc_fid.variables['TOT_DUST'][:,19,:]
volsm=nc_fid.variables['TOT_DUST'][:,20,:]
erodtot=nc_fid.variables['TOT_DUST'][:,21,:]
ustar=nc_fid.variables['TOT_DUST'][:,22,:]
ilwi=nc_fid.variables['TOT_DUST'][:,23,:]
gravsm=nc_fid.variables['TOT_DUST'][:,24,:]
airden=nc_fid.variables['TOT_DUST'][:,25,:]

massfrac[0][:]=nc_fid.variables['CLAYFRAC'][0,:]
massfrac[1][:]=1-(nc_fid.variables['CLAYFRAC'][0,:]+nc_fid.variables['SANDFRAC'][0,:])
massfrac[2][:]=nc_fid.variables['SANDFRAC'][0,:]

#print drylimit.shape, volsm.shape, erodtot.shape, ustar.shape, ilwi.shape, gravsm.shape, airden.shape, massfrac.shape  #(15, 44, 44)
#print massfrac.shape 		#(3, 44, 44)
#exit()

print "processing " +wrf_dir+wrf_out_file
for time_idx in range(1,len(times),1):
	print ''.join(times[time_idx])
	#time_labels.append(''.join(times[time_idx]))

	fig = plt.figure(figsize=(12,12))
	ash_map = Basemap(projection='cyl', llcrnrlat=8.0, urcrnrlat=45.5,llcrnrlon=24.0, urcrnrlon=59.5, resolution='l',area_thresh=10000.)#,ax=axes[k])
	x, y = ash_map(xlon,xlat)

	decorateMap(ash_map)	
	#plot_cities(ash_map)

	date_time_obj = datetime.strptime(''.join(times[time_idx]), '%Y-%m-%d_%H:%M:%S')
	
	emissions,u_ts,u_tres=source_dust(nx,ny,ustar[time_idx], massfrac,erodtot[time_idx], ilwi[time_idx], gravsm[time_idx], volsm[time_idx], airden[time_idx], drylimit[time_idx],xland)
	
	total_emission_flux=np.sum(surface*emissions)
	plt.title(date_time_obj.strftime("%d %B, %H:%M %p")+"\n Total emission flux: "+"{:0.1f}".format(total_emission_flux)+" kg/sec")

	cs=ash_map.pcolormesh(x,y,emissions,cmap=ncview_colormap_short, norm=ai_norm)
	#cs=ash_map.contourf(x,y,emissions,cmap=ncview_colormap_short, norm=ai_norm)
	#OR cs=ash_map.pcolormesh(x,y,emissions,norm=colors.LogNorm(1.0e-7, vmax=0.01))

	cbar = fig.colorbar(cs,orientation='horizontal')
	cbar.set_label('Dust emissions, (kg m^-2 s^-1)')#,fontsize=CB_LABEL_TEXT_SIZE)

	plt.savefig(str(time_idx)+".png")	

	####################################################################################
	####################################################################################
	
	#plot u_ts_XXX's calculated by this script
	k=u_ts.shape[0]
	fig, axs = plt.subplots(1, k,figsize=(k*5,5))
	for l in np.arange(0,k):
		cs=axs[l].pcolormesh(u_ts[l],cmap=cm.get_cmap('rainbow', 21),norm=uts_norm)
		axs[l].set_title('u_ts_'+str(l))
	fig.colorbar(cs, ax=axs.ravel().tolist(),orientation='horizontal')
	plt.savefig("uts_"+str(time_idx)+".png")
	############################
	############################
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
	####################################################################################
	####################################################################################
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

	####################################################################################
	emissions=nc_fid.variables['EDUST1'][time_idx,0]+nc_fid.variables['EDUST2'][time_idx,0]+nc_fid.variables['EDUST3'][time_idx,0]+nc_fid.variables['EDUST4'][time_idx,0]+nc_fid.variables['EDUST5'][time_idx,0]
	total_emission_flux=np.sum(surface*emissions)

	fig = plt.figure(figsize=(12,12))
	ash_map = Basemap(projection='cyl', llcrnrlat=8.0, urcrnrlat=45.5,llcrnrlon=24.0, urcrnrlon=59.5, resolution='l',area_thresh=10000.)#,ax=axes[k])
	x, y = ash_map(xlon,xlat)

	decorateMap(ash_map)	
	#plot_cities(ash_map)

	plt.title("\n Total emission flux: "+"{:0.1f}".format(total_emission_flux)+" kg/sec")

	ai_norm = colors.BoundaryNorm(np.logspace(-12.0, -7.0, num=11), ncview_colormap_short.N, clip=True)
	cs=ash_map.pcolormesh(x,y,emissions,cmap=ncview_colormap_short, norm=ai_norm)
	#cs=ash_map.contourf(x,y,emissions,cmap=ncview_colormap_short, norm=ai_norm)
	
	#OR cs=ash_map.pcolormesh(x,y,emissions,norm=colors.LogNorm(1.0e-7, vmax=0.01))

	cbar = fig.colorbar(cs,orientation='horizontal')
	cbar.set_label('Dust emissions, (kg m^-2 s^-1)')#,fontsize=CB_LABEL_TEXT_SIZE)

	plt.savefig("ornlg_"+str(time_idx)+".png")
	####################################################################################
	'''

nc_fid.close()





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