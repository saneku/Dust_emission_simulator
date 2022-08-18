import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def source_dust(nx,ny,ustar, massfrac,erod, ilwi, gravsm, volsm, airden, drylimit,xland):
  reff_salt=np.array([0.71e-6,1.37e-6,2.63e-6,5.00e-6,9.50e-6,18.1e-6,34.5e-6,65.5e-6,125.0e-6])
  den_salt=np.array([2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650.])
  spoint=np.array([0,1,1,1,1,1,2,2,2])#  ! 0 Clay, 1 Silt, 2 Sand
  frac_salt=np.array([1.0,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333])
  
  #distr_dust=np.array([1.074e-1,1.012e-1,2.078e-1,4.817e-1,1.019e-1]) #sum(distr_dust)=1

  #dt=3600 #seconds
  g0=9.8
  g = g0*1.0E2                          # (cm s^-2)  
  smx=9

  alpha=0.5
  gamma=1.0
  smois_opt=0

  #smtune=1.0
  #volsm=volsm*smtune

  # Global tuning constant, alpha.  
  # Sandblasting mass efficiency, beta.
  # Beta maxes out for clay fractions above 0.2 => betamax.

  betamax=5.25E-4

  # Experimental optional exponential tuning constant for erodibility.
  # 0 < gamma < 1 -> more relative impact by low erodibility regions.
  # 1 < gamma -> accentuates high erodibility regions.  Recommend this
  # be set to 1 (default) unless looking for a way to increase spread
  # within an ensemble framework.

  # Constant of proportionality from Marticorena et al, 1997 (unitless)
  # Arguably more ~consistent~ fudge than alpha, which has many walnuts
  # sprinkled throughout the literature. 

  cmb=1.0    #! Marticorena et al,1997
  #! REAL, PARAMETER :: cmb=2.61   ! White,1979

  # Calculate saltation surface area distribution from sand, silt, and clay
  # mass fractions and saltation bin fraction. This will later become a 
  # modifier to the total saltation flux.  The reasoning here is that the 
  # size and availability of saltators affects saltation efficiency. Based
  # on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

  dsurface=np.zeros(shape=(smx,ny,nx))
  for n in np.arange(0,smx):
    dmass=massfrac[spoint[n]]*frac_salt[n]
    dsurface[n,:]=0.75*dmass/(den_salt[n]*reff_salt[n])  
    
  # The following equation yields relative surface area fraction.  It will only
  # work if you are representing the "full range" of all three soil classes.
  # For this reason alone, we have incorporated particle sizes that encompass
  # the clay class, to account for its relative area over the basal
  # surface, even though these smaller bins would be unlikely to play any large
  # role in the actual saltation process.

  stotal=np.sum(dsurface,axis=0)

  ds_rel=np.zeros(shape=(smx,ny,nx))
  for n in np.arange(0,smx):
    ds_rel[n,:]=dsurface[n,:]/stotal

  # Calculate total dust emission due to saltation of sand sized particles.
  # Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
  # u_ts0 for moisture to get threshold friction velocity (u_ts). Then
  # calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally, 
  # calculate total dust emission (tot_emit), taking into account erodibility. 

  emit=np.zeros(shape=(ny,nx))
  u_ts0=np.zeros(shape=(smx,ny,nx))
  u_ts=np.zeros(shape=(smx,ny,nx))

  for n in np.arange(0,smx):            # Loop over saltation bins
      den = den_salt[n]*1.0e-3          # (g cm^-3)
      diam = 2.0*reff_salt[n]*1.0e2     # (cm)
      rhoa = airden*1.0e-3              # (g cm^-3)

      #print n,den,diam
      #print 1000*np.min(rhoa),1000*np.mean(rhoa),1000*np.max(rhoa)

      # Threshold friction velocity as a function of the dust density and
      # diameter from Bagnold (1941) (m s^-1).

      u_ts0[n] = 0.0013*(np.sqrt(den*g*diam/rhoa) * np.sqrt(1.0+0.006/(den*g*diam**2.5))) / (np.sqrt(1.928*(1331.0*diam**1.56+0.38)**0.092-1.0)) 

      # Friction velocity threshold correction function based on physical
      # properties related to moisture tension. Soil moisture greater than
      # dry limit serves to increase threshold friction velocity (making
      # it more difficult to loft dust). When soil moisture has not reached
      # dry limit, treat as dry (no correction to threshold friction velocity). GAC

      # Calculate threshold friction velocity. If volumetric (gravimetric)
      # water content is less than the drylimit, then the threshold friction
      # velocity (u_ts) will be equal to the dry threshold friction velocity
      # (u_ts0). EDH


      #print u_ts0.shape
      #print volsm.shape
      #print drylimit.shape
      #exit()
      salt=np.zeros(shape=(ny,nx))
      for j in np.arange(0,ny):
        for i in np.arange(0,nx):
          if(xland[j,i]<1.5):
            if (smois_opt==1):
              if (100.0*volsm[j,i] > drylimit[j,i]):
                u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*np.sqrt(1.0+1.21*(100.0*volsm[j,i]-drylimit[j,i])**0.68))
              else:
                u_ts[n,j,i] = u_ts0[n,j,i]
            else:
              if (gravsm[j,i] > drylimit[j,i]):
                u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*np.sqrt(1.0+1.21*(gravsm[j,i]-drylimit[j,i])**0.68))
              else:
                u_ts[n,j,i] = u_ts0[n,j,i]

            # Saltation flux (kg m^-1 s^-1) from MB95
            # ds_rel is the relative surface area distribution

            if ((ustar[j,i] > u_ts[n,j,i]) and (erod[j,i] > 0.0) and (ilwi[j,i]==1.0)):
              salt[j,i] = cmb*ds_rel[n,j,i]*(airden[j,i]/g0)*(ustar[j,i]**3)*(1.0 + u_ts[n,j,i]/ustar[j,i])*(1. - (u_ts[n,j,i]**2)/(ustar[j,i]**2))
            else:
              salt[j,i] = 0.0

            # Calculate total vertical mass flux (note beta has units of m^-1)
            # Beta acts to tone down dust in areas with so few dust-sized particles
            # that the lofting efficiency decreases.  Otherwise, super sandy zones
            # would be huge dust producers, which is generally not the case.
            # Equation derived from wind-tunnel experiments (see MB95).

            beta=10**(13.6*massfrac[0,j,i]-6.0)  # (unitless)
            if (beta > betamax):
              beta=betamax

            emit[j,i]=emit[j,i]+salt[j,i]*(erod[j,i]**gamma)*alpha*beta    # (kg m^-2 sec-1)
          else:
            u_ts0[n,j,i]=0


  return emit,u_ts,u_ts0   #np.sum(emit,axis=0)  # (kg m^-2 s^-1)
  

  '''
   # Loop over saltation bins
  for n in np.arange(0,smx):           
    dsrc = emit*distr_dust(n)*dt1  ! (kg m^-2) per dt1
  IF (dsrc < 0.0) dsrc = 0.0
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


def decorateMap(ash_map):
  ash_map.drawcountries(linewidth=0.25)
  ash_map.drawcoastlines(linewidth=0.25)
  ash_map.drawmapboundary(linewidth=0.25)
  #ash_map.drawrivers()

  # draw parallels.
  parallels = np.arange(0.,90,10.)
  ash_map.drawparallels(parallels,labels=[1,0,0,0],fontsize=12,linewidth=0.2)
  # draw meridians
  meridians = np.arange(0.,180.,10.)
  ash_map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12,linewidth=0.2)
