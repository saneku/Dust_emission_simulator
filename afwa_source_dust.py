import numpy as np

def afwa_source_dust(nx, ny, ustar, massfrac, erod,isltyp, smois, airden,xland,znt,snowh,**tuning_params):
  reff_salt=np.array([0.71e-6,1.37e-6,2.63e-6,5.00e-6,9.50e-6,18.1e-6,34.5e-6,65.5e-6,125.0e-6])
  den_salt=np.array([2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650.])
  spoint=np.array([0,1,1,1,1,1,2,2,2])#  ! 0 Clay, 1 Silt, 2 Sand
  frac_salt=np.array([1.0,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333])
  porosity=np.array([0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.000, 0.200, 0.421, 0.468, 0.200,0.339])

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

  # isltyp is a 1-based USDA soil category; porosity is a 0-based lookup table.
  isltyp_idx=isltyp.astype(int)-1
  poros=porosity[isltyp_idx]

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

  drylimit=np.zeros(shape=(ny,nx))
  #use drylimit based on gravimetric soil moisture
  if smois_opt==0:
    drylimit=14.0*massfrac[0]*massfrac[0]+17.0*massfrac[0]  #! 0 Clay, 1 Silt, 2 Sand should have gravimetric SM units
  # no substantial difference between Noah LSM scheme and RUCLSMSCHEME, PXLSMSCHEME schemes
  #use drylimit based on volumetric soil moisture
  elif smois_opt==1:
    #Noah LSM scheme
    drylimit = 0.0756*(15.127*massfrac[0]+3.09)**2.3211   #! 0 Clay, 1 Silt, 2 Sand should have gravimetric SM units
    #or RUCLSMSCHEME, PXLSMSCHEME
    #drylimit = 0.035 * (13.52 * massfrac[0] + 3.53)**2.68
  elif smois_opt!=2:
    raise ValueError(f"Unsupported smois_opt: {smois_opt}")


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

  dmass=massfrac[spoint]*frac_salt[:,None,None]
  dsurface=0.75*dmass/(den_salt[:,None,None]*reff_salt[:,None,None])

  # The following equation yields relative surface area fraction.  It will only
  # work if you are representing the "full range" of all three soil classes.
  # For this reason alone, we have incorporated particle sizes that encompass
  # the clay class, to account for its relative area over the basal
  # surface, even though these smaller bins would be unlikely to play any large
  # role in the actual saltation process.

  stotal=np.sum(dsurface,axis=0)
  ds_rel=dsurface/stotal[None,:,:]

  # Calculate total dust emission due to saltation of sand sized particles.
  # Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
  # u_ts0 for moisture to get threshold friction velocity (u_ts). Then
  # calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally,
  # calculate total dust emission (tot_emit), taking into account erodibility.

  emit=np.zeros(shape=(ny,nx))
  u_ts0=np.zeros(shape=(smx,ny,nx))
  u_ts=np.zeros(shape=(smx,ny,nx))

  # Grid-cell masks and soil-moisture quantities below are independent of the
  # saltation bin (n), so compute them once instead of inside the bin loop.

  ilwi=np.where(xland<1.5, 1.0, 0.0)
  ilwi=np.where(znt>0.2, 0.0, ilwi)
  ilwi=np.where((isltyp==14)|(isltyp==15)|(isltyp==17), 0.0, ilwi)
  ilwi=np.where(snowh>0.01, 0.0, ilwi)

  # Calculate volumetric and gravimetric soil moisture.
  volsm=np.maximum(smois*smtune, 0.0)
  gravsm=100.0*volsm/((1.0-poros)*(2.65*(1-massfrac[0])+2.50*massfrac[0]))

  # Friction velocity threshold correction function based on physical
  # properties related to moisture tension. Soil moisture greater than
  # dry limit serves to increase threshold friction velocity (making
  # it more difficult to loft dust). When soil moisture has not reached
  # dry limit, treat as dry (no correction to threshold friction velocity). GAC

  # Calculate threshold friction velocity. If volumetric (gravimetric)
  # water content is less than the drylimit, then the threshold friction
  # velocity (u_ts) will be equal to the dry threshold friction velocity
  # (u_ts0). EDH

  if smois_opt==0:
    #Use Gravimetric SM
    diff=gravsm-drylimit
    factor=np.where(diff>0.0, np.sqrt(1.0+1.21*np.maximum(diff,0.0)**0.68), 1.0)
    wet_mask=None
  elif smois_opt==1:
    #Use Volumetric SM
    diff=100.0*volsm-drylimit
    factor=np.where(diff>0.0, np.sqrt(1.0+1.21*np.maximum(diff,0.0)**0.68), 1.0)
    wet_mask=None
  else:
    #Simple GOCART variant: volumetric soil moisture over porosity
    gwet=volsm/poros
    wet_mask=gwet>=0.5   # Case of wet surface, no erosion
    factor=1.2+0.2*np.log10(np.maximum(0.001, gwet))

  # Beta acts to tone down dust in areas with so few dust-sized particles
  # that the lofting efficiency decreases.  Otherwise, super sandy zones
  # would be huge dust producers, which is generally not the case.
  # Equation derived from wind-tunnel experiments (see MB95).

  beta=10**(13.6*massfrac[0]-6.0)  # (unitless)
  beta=np.minimum(beta, betamax)

  for n in np.arange(0,smx):            # Loop over saltation bins
      den = den_salt[n]*1.0e-3          # (g cm^-3)
      diam = 2.0*reff_salt[n]*1.0e2     # (cm)
      rhoa = airden*1.0e-3              # (g cm^-3)

      # Threshold friction velocity as a function of the dust density and
      # diameter from ###Bagnold (1941) (m s^-1).
      #from Marticorena and Bergametti (1995) eq.4.

      u_ts0[n] = 0.0013*(np.sqrt(den*g*diam/rhoa) * np.sqrt(1.0+0.006/(den*g*diam**2.5))) / (np.sqrt(1.928*(1331.0*diam**1.56+0.38)**0.092-1.0))

      #Use Gravimetric/Volumetric SM
      if wet_mask is None:
        u_ts[n] = np.maximum(0.0, u_ts0[n]*factor)
      else:
        #Simple GOCART variant
        u_ts_dry = np.maximum(0.0, u_ts0[n]*factor)
        # Case of wet surface, no erosion
        u_ts[n] = np.where(wet_mask, 100.0, u_ts_dry)

      # Saltation flux (kg m^-1 s^-1) from MB95
      # ds_rel is the relative surface area distribution

      erodible = (ustar > u_ts[n]) & (erod > 0.0) & (ilwi==1.0)
      with np.errstate(divide='ignore', invalid='ignore'):
        salt = np.where(erodible,
                        cmb*ds_rel[n]*(airden/g0)*(ustar**3)*(1.0 + u_ts[n]/ustar)*(1. - (u_ts[n]**2)/(ustar**2)),
                        0.0)

      # Calculate total vertical mass flux (note beta has units of m^-1)

      emit = emit + salt*(erod**gamma)*alpha*beta    # (kg m^-2 sec-1)

  return emit,u_ts,u_ts0
