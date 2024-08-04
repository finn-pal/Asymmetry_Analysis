# Imports
import math
import os
import time

import numpy as np
import pandas as pd
from astropy.io import fits
from functions.get_kinematics import get_kinematics
from functions.get_sfr import get_sfr
from photutils.aperture import EllipticalAperture

start_time = time.time()

# directories
FOLDER_DIR_DATA = "data/"  # data folder
CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder
FIELD_DIR = FOLDER_DIR_DATA + "fields/"  # galaxy fields
GALAXY_DIR = FOLDER_DIR_DATA + "galaxy_images/"  # galaxy images
ARC_SCALE = 0.2  # Arcsec / Pixel scale

# Galaxy_Data_Import
GALAXY_DETAIL_DATA_FILE = "galaxy_list_phot.csv"
df_galaxy_list_det = pd.read_csv(CSV_DIR + GALAXY_DETAIL_DATA_FILE)
del df_galaxy_list_det["Unnamed: 0"]

arr = [images for images in os.listdir(GALAXY_DIR) if images.endswith(".fits")]
arr_sorted = sorted(arr)

df_parent = get_kinematics(df_galaxy_list_det)
df_parent = get_sfr(df_parent)


def SNR_calc(image, im_var, cent, a, e_2, theta):
    b = np.sqrt((a**2) * (1 - e_2))
    elli = EllipticalAperture(cent, a, b, theta)
    n = elli.area

    elli_mask = elli.to_mask()

    ima_cut = elli_mask.get_values(image)
    var_cut = elli_mask.get_values(im_var)

    sum_arr = []

    for i in range(0, len(ima_cut)):
        hold = ima_cut[i] / np.sqrt(var_cut[i] + ima_cut[i])

        sum_arr.append(hold)

    SNR = np.nansum(sum_arr) / n

    return SNR


SNR_arr = []

for i in range(0, len(df_parent)):
    # Opens the correct fits file that matches the idx from the dataframe
    hdu1 = fits.open(GALAXY_DIR + arr_sorted[i])

    im_g = hdu1[2].data  # G_mod band image
    st_g = hdu1[3].data  # G_mod band stat
    im_i = hdu1[4].data  # i band image
    st_i = hdu1[5].data  # G_mod band stat
    im_r = hdu1[6].data  # i band image
    st_r = hdu1[7].data  # G_mod band stat
    mask = hdu1[10].data  # undilated_mask

    image = im_g + im_i + im_r  # sum of g mod, i and r band images
    im_var = st_g + st_i + st_r  # sum of g mod, i and r band stats
    im_err = np.sqrt(im_var)

    SEG_ID = df_parent.segID[i]  # Gets SEG ID for galaxy

    # Transforms the mask such that the mask is equal to 1, background is set to 0 and everything else
    # (other galaxies) are nans.
    mask = mask.astype(float)
    mask[(mask != SEG_ID) & (mask != 0)] = np.nan
    mask[mask == SEG_ID] = 1

    theta = df_parent.pa_galfit[i]
    theta = math.radians(theta + 90)  # 90 degree correction due to python

    re = df_parent.re_galfit[i]
    re = re / ARC_SCALE  # Convert to pixels

    q_val = df_parent.q_galfit[i]
    eta = 1 - q_val  # Ellipticity
    e_2 = eta * (2 - eta)  # Eccentricity squared

    cent = [df_parent.centre_asym_x[i], df_parent.centre_asym_y[i]]

    print(arr_sorted[i])
    print()
    SNR = SNR_calc(image, im_var, cent, re, e_2, theta)
    SNR_arr.append(SNR)

df_parent["SNR_Phot"] = SNR_arr

##### Phot_Flags start #########################################################

nan_frac_limit = 0.5
SNR_cut_off = 3
C_cut = 5
Mag_cut = 23


# Stellar kinematic flag formula selection
def Nan_Frac_Flag(row):
    if (
        (row["nan_frac"] < nan_frac_limit)
        & (row["SNR_Phot"] > SNR_cut_off)
        & (~np.isnan(row["P_asym"]))
        & (row["C"] <= C_cut)
        & (row["mag_it"] < Mag_cut)
    ):
        return 1

    else:
        return 0


df_parent = df_parent.assign(Phot_Flag=df_parent.apply(Nan_Frac_Flag, axis=1))
df_parent[df_parent.Phot_Flag.isnull()] = 0

##### Phot_Flags end ###########################################################

##### Kin_Flags start ##########################################################

# Only select galaxies with gas and stellar resolutions above this fraction
stel_kin_frac_lim = 0.80
gas_kin_frac_lim = 0.80
SNR_kin_limit = 5

########## ADD the following to both when get data from Ryan


# Stellar kinematic flag formula selection
def Kinematic_Flag_s(row):
    if (
        (row["stelkinfracRe"] >= stel_kin_frac_lim)
        & (row["SNR_s"] > SNR_kin_limit)
        & (~np.isnan(row["v_asym_s"]))
    ):
        return 1

    else:
        return 0


# Gas kinematic flag formula selection
def Kinematic_Flag_g(row):
    if (
        (row["gaskinfracRe"] >= gas_kin_frac_lim)
        & (row["SNR_g"] > SNR_kin_limit)
        & (~np.isnan(row["v_asym_g"]))
    ):
        return 1

    else:
        return 0


# Apply to magpi parent csv

df_parent = df_parent.assign(Kin_Flag_s=df_parent.apply(Kinematic_Flag_s, axis=1))
df_parent = df_parent.assign(Kin_Flag_g=df_parent.apply(Kinematic_Flag_g, axis=1))

# Redshift limit

df_parent = df_parent[(df_parent["z.x"] > 0.25) & (df_parent["z.x"] < 0.35)]

##### Kin_Flags end ############################################################

# Upadate and print new Galaxy List information to csv
filepath = CSV_DIR + "galaxy_parent.csv"
df_parent.to_csv(filepath)

print("--- %s seconds ---" % (time.time() - start_time))
