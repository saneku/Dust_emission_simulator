import numpy as np
import matplotlib


def gocart_source_dust(nx, ny, w10m, isltyp, smois, erod, airden, xland):
    den_dust = np.array([2500.0, 2650.0, 2650.0, 2650.0, 2650.0])
    reff_dust = np.array([0.73e-6, 1.4e-6, 2.4e-6, 4.5e-6, 8.0e-6])
    ipoint = np.array([2, 1, 1, 1, 1])
    frac_s = np.array([0.1, 0.25, 0.25, 0.25, 0.25])
    porosity=([0.339, 0.421, 0.434, 0.476, 0.476, 0.439, 0.404, 0.464, 0.465, 0.406, 0.468, 0.468, 0.439, 1.000, 0.200, 0.421, 0.468, 0.200,0.339])

    g = 9.8 * 1.0e2  # (cm s^-2)
    nmx = 5
    ch_dust = 0.8e-9

    # Calculate total dust emission of sand sized particles.
    # Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
    # u_ts0 for moisture to get threshold friction velocity (u_ts). Then
    # calculate vertical flux where w10m has exceeded u_ts.  Finally,
    # calculate total dust emission (dsrc).

    emit = np.zeros(shape=(ny, nx))
    u_ts0 = np.zeros(shape=(nmx, ny, nx))
    u_ts = np.zeros(shape=(nmx, ny, nx))
    dsrc = np.zeros(shape=(nmx, ny, nx))

    # Threshold velocity as a function of the dust density and the diameter from Bagnold (1941)
    for n in np.arange(0, nmx):  # Loop over saltation bins
        den = den_dust[n] * 1.0e-3  # (g cm^-3)
        diam = 2.0 * reff_dust[n] * 1.0e2  # (cm)
        rhoa = airden * 1.0e-3  # (g cm^-3)

        # Pointer to the 3 classes considered in the source data files
        m = ipoint[n]

        # Threshold friction velocity (u_ts0) as a function of the dust density and
        # diameter from Marticorena and Bergametti (1995) eq.4.
        u_ts0[n] = 0.0013*(np.sqrt(den*g*diam/rhoa) * np.sqrt(1.0+0.006/(den*g*diam**2.5))) / (np.sqrt(1.928*(1331.0*diam**1.56+0.38)**0.092-1.0)) 
        
        for j in np.arange(0, ny):
            for i in np.arange(0, nx):
                if xland[j, i] < 1.5:
                    #  volumetric soil moisture over porosity
                    gwet = smois[j, i] / porosity[isltyp[j, i]]

                    # Case of surface dry enough to erode
                    if gwet < 0.5:
                        u_ts[n,j,i] = max(0.0,u_ts0[n,j,i]*(1.2+0.2*np.log10(max(0.001, gwet))))
                    else:
                        # Case of wet surface, no erosion
                        u_ts[n, j, i] = 100

                    srce = frac_s[n] * erod[m, j, i]  # (kg s^2 m^-5)
                    # if (ilwi[j,i]==1.0):
                    # (kg s^2 m^-5)*(m^3 s^-3) = (kg/m2/sec)
                    dsrc = ch_dust*srce*w10m[j,i]**2*(w10m[j,i] - u_ts[n,j,i])

                    dsrc = max(0.0, dsrc)  # (kg m^-2 sec-1)

                    emit[j, i] = emit[j, i] + dsrc
                else:
                    u_ts0[n, j, i] = 0

    return emit, u_ts, u_ts0
