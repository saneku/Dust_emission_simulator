import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature

units='($kg\ m^{-2}\ s^{-1}$)'

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

# Lambert Conformal projection matching the original Basemap ('lcc') setup.
map_proj = ccrs.LambertConformal(central_longitude=cen_lon, central_latitude=cen_lat,
                                  standard_parallels=(true_lat1, true_lat2))

# Map extent (in projected meters), centered on (cen_lon, cen_lat), matching
# the width/height Basemap used to draw.
_x0, _y0 = map_proj.transform_point(cen_lon, cen_lat, ccrs.PlateCarree())
map_extent = (_x0 - dx*nxd/2.0, _x0 + dx*nxd/2.0, _y0 - dy*nyd/2.0, _y0 + dy*nyd/2.0)

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

def decorateMap(ax):
    ax.set_extent(map_extent, crs=map_proj)

    ax.coastlines(resolution='50m', linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.2)
    ax.add_feature(cfeature.STATES, linewidth=0.2)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    #domain boundary
    wrf_lons=np.concatenate((xlon[:,0],xlon[ny-1,:],xlon[:,nx-1][::-1],xlon[0,:][::-1]), axis=0)
    wrf_lats=np.concatenate((xlat[:,0],xlat[ny-1,:],xlat[:,nx-1][::-1],xlat[0,:][::-1]), axis=0)
    ax.plot(wrf_lons, wrf_lats, marker=None, color='brown', transform=ccrs.PlateCarree())