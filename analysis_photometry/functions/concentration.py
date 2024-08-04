# Imports
import math

import numpy as np
from photutils.aperture import EllipticalAperture, aperture_photometry
from scipy.optimize import minimize_scalar


def Median_Set(data, mask):
    nx = np.shape(data)[0]
    ny = np.shape(data)[1]

    for i in range(0, nx):
        for j in range(0, ny):
            if (np.isnan(data[i, j])) & (mask[i, j] == 1):
                try:
                    med_val = np.nanmedian(
                        [
                            data[i - 1, j - 1],
                            data[i, j - 1],
                            data[i + 1, j - 1],
                            data[i - 1, j],
                            data[i + 1, j],
                            data[i - 1, j + 1],
                            data[i, j + 1],
                            data[i + 1, j + 1],
                        ]
                    )
                    data[i, j] = med_val

                except IndexError:
                    0

    return data


def elli_calc(a, theta, cent, e_2, T_flux, data, lim):
    # Using major axis and eccentricity squared calculate the minor axis
    b = np.sqrt((a**2) * (1 - e_2))

    # Fit an ellipse
    elli = EllipticalAperture(cent, a, b, theta)

    # Get data from fit
    phot_table = aperture_photometry(data, elli)
    app_flux = phot_table["aperture_sum"][0]

    ratio = round(app_flux / T_flux, 4)

    # the minimise will find the largest sized ellipse that fits the limiting criteria
    if ratio < lim:
        ratio = np.inf

    else:
        ratio = ratio

    return ratio


def a_error(a, e_2, e_2_err, theta, theta_err, T_flux, T_flux_var):
    b = np.sqrt((a**2) * (1 - e_2))
    b_err = b * (e_2_err / e_2)

    a_err = a * np.sqrt((b_err / b) ** 2 + (theta_err / theta) ** 2 + T_flux_var / (T_flux**2))

    return a_err


def Concentration(df_galaxy_list, centre_asym, gal_image, var_image, back_noise_sum, mask, ARC_SCALE, i):
    # Set nans go to zero
    data = gal_image.copy()
    data = Median_Set(data, mask)
    data[np.isnan(data)] = 0

    cent = centre_asym  # Centre of asymmetry

    # Get position angle and convert to radians with 90 degree python correction
    theta = df_galaxy_list["pa_galfit"][i]
    theta = math.radians(theta + 90)  # 90 degree correction due to python

    theta_err = df_galaxy_list["pa_err_galfit"][i]  # Error in position angle
    theta_err = math.radians(theta_err)  # Error in position angle converted to degrees

    # Get q value so that ellipticity can be obtained and then in turn converted to eccentricity squared
    q_val = df_galaxy_list["q_galfit"][i]
    q_val_err = df_galaxy_list["q_err_galfit"][i]

    eta = 1 - q_val  # Ellipticity
    eta_err = q_val_err

    e_2 = eta * (2 - eta)  # Eccentricity squared
    e_2_err = e_2 * np.sqrt(eta_err**2 + ((eta**2) * np.sqrt(2 * ((eta_err / eta) ** 2))) ** 2)

    # Effective radius of galaxy which can be used as bounds for Minimize Scalar
    re = df_galaxy_list["re_galfit"][i]
    re_pix = re / ARC_SCALE  # Convert to pixels

    # Background noise information
    avg_noise = np.average(back_noise_sum)
    var_noise = np.var(back_noise_sum)

    # Total flux in image with error but not background corrected
    T_flux_unc = np.nansum(data)  # Uncorrected
    T_flux_unc_var = np.nansum(var_image)  # Uncorrected

    # Total flux in image with error and background corrected
    T_flux = T_flux_unc - avg_noise
    T_flux_var = T_flux_unc_var + var_noise

    # Radius of ellipse containing 20 percent of total flux
    res_20 = minimize_scalar(
        elli_calc, bounds=(0, re_pix * 2), method="bounded", args=(theta, cent, e_2, T_flux, data, 0.20)
    )
    a_20 = res_20.x

    # Radius of ellipse containing 80 percent of total flux
    res_80 = minimize_scalar(
        elli_calc, bounds=(a_20, re_pix * 4), method="bounded", args=(theta, cent, e_2, T_flux, data, 0.80)
    )
    a_80 = res_80.x

    # For error analyis
    a_20_err = a_error(a_20, e_2, e_2_err, theta, theta_err, T_flux, T_flux_var)
    a_80_err = a_error(a_80, e_2, e_2_err, theta, theta_err, T_flux, T_flux_var)

    # Concentration with errors
    C = 5 * np.log10(a_80 / a_20)
    C_err = C * np.sqrt((a_80_err / a_80) ** 2 + (a_20_err / a_20) ** 2)

    return C, C_err
