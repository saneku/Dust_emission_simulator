import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import matplotlib.colors as colors

units='($kg\ m^{-2}\ s^{-1}$)'
basemap_params = dict(projection='cyl', llcrnrlat=8.0, urcrnrlat=45.5,llcrnrlon=24.0, urcrnrlon=59.5, resolution='l',area_thresh=10000.)

wrf_input_file='grid.nc'
wrf_dir='./data/'

from matplotlib import cm
colmap =cm.get_cmap('gist_ncar_r')
ai_norm = colors.BoundaryNorm(np.logspace(-12.0, -7.0, num=11), colmap.N, clip=False)


wrfinput=Dataset(wrf_dir+"/"+wrf_input_file,'r')
xlon=wrfinput.variables['XLONG'][0,:]
xlat=wrfinput.variables['XLAT'][0,:]

dy=100000.0 #meters
dx=100000.0
nx=44 #points
ny=44

MAPFAC_MX=wrfinput.variables['MAPFAC_MX'][0,:]
MAPFAC_MY=wrfinput.variables['MAPFAC_MY'][0,:]

surface=(dx/MAPFAC_MX)*(dy/MAPFAC_MY)       #surface in m2
wrfinput.close()

MAP_TICKS_TEXT_SIZE=10

'''
def plot_cities(ash_map):
    points=['Jeddah ',' Riyadh','Dammam '," Sana'a"," Muskat"," Mecca"," Abu Dhabi"," Doha"," Ankara"," Damascus","Beirut "," Amman","Jerusalem "," Tehran"," Baghdad"," Kuwait City"," Cairo"," Manama"]
    align={'Jeddah ':"right",' Riyadh':'left','Dammam ':"right",'Yanbu ':"right",'Rabigh ':"right"," Sana'a":'left'," Muskat":"left"," Mecca":"left"," Abu Dhabi":"left"," Medina":"left"," Doha":"left"," Ankara":"left"," Damascus":"left","Beirut ":"right"," Amman":"left","Jerusalem ":"right"," Tehran":"left"," Baghdad":"left"," Kuwait City":"left"," Cairo":"left"," Manama":"left"}
    point_lons=[39.222167,46.891028,49.983083,44.168,58.338,39.831,54.382,51.542,32.863,36.293,35.496,35.923,35.209,51.388,44.358,47.977,31.238,50.582]
    point_lats=[21.416222,24.55266667,26.246639,15.42,23.582,21.427,24.480,25.286,39.938,33.511,33.892,31.972,31.763,35.697,33.330,29.378,30.047,26.221]
    MARKER_SIZE=4
    TEXT_SIZE=12

    x, y = ash_map(point_lons,point_lats)
    ash_map.plot(x, y,'o',markersize=MARKER_SIZE,markeredgewidth=1,markeredgecolor='k',markerfacecolor='none',zorder=11)

    k=0
    for s in points:
        x,y=ash_map(point_lons[k],point_lats[k])
        plt.text(x,y,s,fontsize=TEXT_SIZE,style='italic',ha=align.get(s),va='center',color='k',zorder=11)
        #axes.text(x,y,s,fontsize=TEXT_SIZE,style='italic',ha=align.get(s),va='center',color='k',zorder=11)
        k=k+1
'''

def decorateMap(ash_map):
    ash_map.drawcoastlines(linewidth=0.8)
    ash_map.drawcountries(linewidth=0.2)
    ash_map.drawstates(linewidth=0.2)


    parallels = np.arange(0.,90,10.)
    ash_map.drawparallels(parallels,labels=[1,0,0,0], linewidth=0.3,fontsize=MAP_TICKS_TEXT_SIZE)
    # draw meridians
    meridians = np.arange(0.,180.,10.)
    ash_map.drawmeridians(meridians,labels=[0,0,0,1], linewidth=0.3,fontsize=MAP_TICKS_TEXT_SIZE)
