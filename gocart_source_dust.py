import numpy as np


def gocart_source_dust(nx, ny, w10m, isltyp, smois, erod, airden, xland,**tuning_params):
    den_dust = np.array([2500.0, 2650.0, 2650.0, 2650.0, 2650.0])
    reff_dust = np.array([0.73e-6, 1.4e-6, 2.4e-6, 4.5e-6, 8.0e-6])
    ipoint = np.array([2, 1, 1, 1, 1])
    frac_s = np.array([0.15,0.1,0.25,0.4,0.1])
    porosity=np.array([0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.000, 0.200, 0.421, 0.468, 0.200,0.339])
    g = 9.8 * 1.0e2  # (cm s^-2)
    nmx = 5
    ch_dust = 0.8e-9

    C_factor=tuning_params['C_factor']

    # Calculate total dust emission of sand sized particles.
    # Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    # u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    # calculate vertical flux where w10m has exceeded u_ts.  Finally,
    # calculate total dust emission (dsrc).

    emit = np.zeros(shape=(ny, nx))
    u_ts0 = np.zeros(shape=(nmx, ny, nx))
    u_ts = np.zeros(shape=(nmx, ny, nx))

    land_mask = xland < 1.5

    # isltyp is a 1-based USDA soil category; porosity is a 0-based lookup table.
    poros = porosity[isltyp.astype(int)-1]
    gwet = smois / poros                  # volumetric soil moisture over porosity
    wet_mask = gwet >= 0.5                 # Case of wet surface, no erosion

    rhoa = airden*1.0e-3                   # (g cm^-3)

    # Threshold velocity as a function of the dust density and the diameter from Bagnold (1941)
    for n in np.arange(0, nmx):  # Loop over dust bins
        den = den_dust[n] * 1.0e-3  # (g cm^-3)
        diam = 2.0 * reff_dust[n] * 1.0e2  # (cm)

        # Pointer to the 3 classes considered in the source data files
        m = ipoint[n]

        # Threshold friction velocity (u_ts0) as a function of the dust density and
        # diameter from Marticorena and Bergametti (1995) eq.4.
        u_ts0[n] = 0.0013*(np.sqrt(den*g*diam/rhoa) * np.sqrt(1.0+0.006/(den*g*diam**2.5))) / (np.sqrt(1.928*(1331.0*diam**1.56+0.38)**0.092-1.0))

        u_ts_dry = np.maximum(0.0, u_ts0[n]*(1.2+0.2*np.log10(np.maximum(0.001, gwet))))
        u_ts[n] = np.where(land_mask, np.where(wet_mask, 100.0, u_ts_dry), 0.0)
        u_ts0[n] = np.where(land_mask, u_ts0[n], 0.0)

        srce = frac_s[n] * erod[m]  # (kg s^2 m^-5)
        # (kg s^2 m^-5)*(m^3 s^-3) = (kg/m2/sec)
        with np.errstate(invalid='ignore'):
            dsrc = C_factor*ch_dust*srce*w10m**2*(w10m - u_ts[n])
        dsrc = np.maximum(0.0, dsrc)  # (kg m^-2 sec-1)

        emit = emit + np.where(land_mask, dsrc, 0.0)

    return emit, u_ts, u_ts0
