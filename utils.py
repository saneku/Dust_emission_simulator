import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm

units='($kg\ m^{-2}\ s^{-1}$)'

proj='lcc'
cen_lat = 27.96201
cen_lon = 42.022
true_lat1 = 27.962
true_lat2 = 27.962

dy=100000.0 #meters
dx=100000.0
nx=44 #points
ny=44

nyd=46
nxd=46

basemap_params = dict(width=dx*nxd,height=dy*nyd,resolution='l',area_thresh=100.,projection=proj, lat_1=true_lat1,lat_2=true_lat2,lat_0=cen_lat,lon_0=cen_lon)

wrf_input_file='grid.nc'
wrf_dir='./data/'

c_map = plt.cm.get_cmap('Oranges')
colorlist=list()
colorlist.append("#ffffff")
for c in np.linspace(0,1,plt.cm.Oranges.N):
    rgba=c_map(c) #select the rgba value of the cmap at point c which is a number between 0 to 1
    clr=colors.rgb2hex(rgba) #convert to hex
    colorlist.append(str(clr)) # create a list of these colors

colmap = colors.LinearSegmentedColormap.from_list('cmap_name', colorlist, N=10)
colmap.set_over(color='k')
ai_norm = colors.BoundaryNorm(np.logspace(-12.0, -7.0, num=11), colmap.N, clip=False)


wrfinput=Dataset(wrf_dir+"/"+wrf_input_file,'r')
xlon=wrfinput.variables['XLONG'][0,:]
xlat=wrfinput.variables['XLAT'][0,:]


MAPFAC_MX=wrfinput.variables['MAPFAC_MX'][0,:]
MAPFAC_MY=wrfinput.variables['MAPFAC_MY'][0,:]

surface=(dx/MAPFAC_MX)*(dy/MAPFAC_MY)       #surface in m2
wrfinput.close()

def decorateMap(m):
    m.drawcoastlines(linewidth=0.8)
    m.drawcountries(linewidth=0.2)
    m.drawstates(linewidth=0.2)

    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0], linewidth=0.3,fontsize=10)
    # draw meridians
    meridians = np.arange(0.,180.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1], linewidth=0.3,fontsize=10)

    #domain boundary
    wrf_lons=np.concatenate((xlon[:,0],xlon[ny-1,:],xlon[:,nx-1][::-1],xlon[0,:][::-1]), axis=0)
    wrf_lats=np.concatenate((xlat[:,0],xlat[ny-1,:],xlat[:,nx-1][::-1],xlat[0,:][::-1]), axis=0)
    xx, yy = m(wrf_lons, wrf_lats)
    m.plot(xx, yy, marker=None,color='brown')