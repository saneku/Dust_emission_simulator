import numpy as np

def afwa_source_dust(nx, ny, ustar, massfrac, erod,isltyp, smois, airden,xland,znt,snowh,**tuning_params):
  reff_salt=np.array([0.71e-6,1.37e-6,2.63e-6,5.00e-6,9.50e-6,18.1e-6,34.5e-6,65.5e-6,125.0e-6])
  den_salt=np.array([2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650.])
  spoint=np.array([0,1,1,1,1,1,2,2,2])#  ! 0 Clay, 1 Silt, 2 Sand
  frac_salt=np.array([1.0,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333])
  porosity=([0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.000, 0.200, 0.421, 0.468, 0.200,0.339])
  
  #distr_dust=np.array([1.074e-1,1.012e-1,2.078e-1,4.817e-1,1.019e-1]) #sum(distr_dust)=1

  #dt=3600 #seconds
  g0=9.8
  g = g0*1.0E2                          # (cm s^-2)  
  smx=9
  
  alpha=tuning_params.get('alpha')
  gamma=tuning_params.get('gamma')
  smtune=tuning_params.get('smtune')
  ustune=tuning_params.get('ustune')

  erod=np.sum(erod,axis=0)
  volsm=np.zeros(shape=(ny,nx))
  gravsm=np.zeros(shape=(ny,nx))
  drylimit=np.zeros(shape=(ny,nx))
  

  # Choose an LSM option and drylimit option.
  # Drylimit calculations based on look-up table in
  # Clapp and Hornberger (1978) for RUC and PXLSM and
  # Cosby et al. (1984) for Noah LSM.

  '''
  dust_smois=0  
  if dust_smois == 1:
   sfc_select: SELECT CASE(config_flags%sf_surface_physics)
      CASE ( RUCLSMSCHEME, PXLSMSCHEME )
         drylimit(1,1) =0.035*(13.52*clay(i,j)+3.53)**2.68
         smois_opt = 1
      CASE ( LSMSCHEME )
         drylimit(1,1) =0.0756*(15.127*clay(i,j)+3.09)**2.3211
         smois_opt = 1
      CASE DEFAULT

         ! Don't currently support volumetric soil moisture
         ! for this scheme, use drylimit based on gravimetric 

         drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
   END SELECT sfc_select
  else:
    # use drylimit based on gravimetric soil moisture
  '''

  smois_opt = 2         #0 for gravimetric SM, 1 for volumetric SM, 2 for GOCART Simple

  #use drylimit based on gravimetric soil moisture
  if smois_opt==0:
    drylimit=14.0*massfrac[0]*massfrac[0]+17.0*massfrac[0]  #! 0 Clay, 1 Silt, 2 Sand should have gravimetric SM units
  # no substantial difference between Noah LSM scheme and RUCLSMSCHEME, PXLSMSCHEME schemes
  #use drylimit based on volumetric soil moisture
  if smois_opt==1:
    #Noah LSM scheme
    drylimit = 0.0756*(15.127*massfrac[0]+3.09)**2.3211   #! 0 Clay, 1 Silt, 2 Sand should have gravimetric SM units
    #or RUCLSMSCHEME, PXLSMSCHEME   
    #drylimit = 0.035 * (13.52 * massfrac[0] + 3.53)**2.68


  # Friction velocity tuning constant (Note: recommend 0.7 for PXLSM,
  # else use 1.0.  This was created due to make the scheme compatible
  # with the much stronger friction velocities coming out of PXLSM).
  ustar=ustar * ustune

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
  ilwi=np.zeros(shape=(ny,nx))
  fecan=np.zeros(shape=(ny,nx))

  for n in np.arange(0,smx):            # Loop over saltation bins
      den = den_salt[n]*1.0e-3          # (g cm^-3)
      diam = 2.0*reff_salt[n]*1.0e2     # (cm)
      rhoa = airden*1.0e-3              # (g cm^-3)

      #print n,den,diam
      #print 1000*np.min(rhoa),1000*np.mean(rhoa),1000*np.max(rhoa)

      # Threshold friction velocity as a function of the dust density and
      # diameter from ###Bagnold (1941) (m s^-1).
      #from Marticorena and Bergametti (1995) eq.4.

      u_ts0[n] = 0.0013*(np.sqrt(den*g*diam/rhoa) * np.sqrt(1.0+0.006/(den*g*diam**2.5))) / (np.sqrt(1.928*(1331.0*diam**1.56+0.38)**0.092-1.0)) 
      #cs=plt.pcolormesh(u_ts0[n]); plt.colorbar(cs); plt.title(f"{n}"); plt.savefig(f"./{n}.png")
      salt=np.zeros(shape=(ny,nx))
      for j in np.arange(0,ny):
        for i in np.arange(0,nx):

          if(xland[j,i]<1.5):
            ilwi[j,i]=1.0

          if(znt[j,i] > 0.2):
            ilwi[j,i]=0.0

          if(isltyp[j,i]==14 or isltyp[j,i]==15 or isltyp[j,i]==17):
            ilwi[j,i]=0.0

          if(snowh[j,i] > 0.01):
            ilwi[j,i]=0.0

          # Calculate volumetric and gravimetric soil moisture.            
          volsm[j, i]=max(smois[j, i]*smtune,0.0)
          gravsm[j,i]=100.0*volsm[j, i]/((1.0-porosity[isltyp[j, i]])*(2.65*(1-massfrac[0,j,i])+2.50*massfrac[0,j,i]))


          # Friction velocity threshold correction function based on physical
          # properties related to moisture tension. Soil moisture greater than
          # dry limit serves to increase threshold friction velocity (making
          # it more difficult to loft dust). When soil moisture has not reached
          # dry limit, treat as dry (no correction to threshold friction velocity). GAC

          # Calculate threshold friction velocity. If volumetric (gravimetric)
          # water content is less than the drylimit, then the threshold friction
          # velocity (u_ts) will be equal to the dry threshold friction velocity
          # (u_ts0). EDH

          #Use Gravimetric SM
          if (smois_opt==0):
            if (gravsm[j,i] > drylimit[j,i]):
              u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*np.sqrt(1.0+1.21*(gravsm[j,i]-drylimit[j,i])**0.68))
              fecan[j,i] = np.sqrt(1.0+1.21*(gravsm[j,i]-drylimit[j,i])**0.68)
            else:
              u_ts[n,j,i] = u_ts0[n,j,i]
              fecan[j,i] = 1

          #Use Volumetric SM
          if (smois_opt==1):
            if (100.0*volsm[j,i] > drylimit[j,i]):
              u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*np.sqrt(1.0+1.21*(100.0*volsm[j,i]-drylimit[j,i])**0.68))
              fecan[j,i] = np.sqrt(1.0+1.21*(100.0*volsm[j,i]-drylimit[j,i])**0.68)
            else:
              u_ts[n,j,i] = u_ts0[n,j,i]
              fecan[j,i] = 1        

          #Simple GOCART variant
          if (smois_opt==2):
            #  volumetric soil moisture over porosity
            gwet = volsm[j, i] / porosity[isltyp[j, i]]

            # Case of surface dry enough to erode
            if gwet < 0.5:
              u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*(1.2+0.2*np.log10(max(0.001, gwet))))
              fecan[j,i] = 1.2 + 0.2 * np.log10(max(0.001, gwet))
            else:
              # Case of wet surface, no erosion
              u_ts[n, j, i] = 100
              fecan[j,i] = 10


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

  return emit,u_ts,u_ts0