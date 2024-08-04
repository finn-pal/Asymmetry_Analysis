# Imports
import os
import time

import numpy as np
import pandas as pd
import skimage
from astropy.io import fits
from functions.concentration import Concentration
from functions.phot_asymmetry import asym_calc_uncorrected, noise_correction
from photutils.centroids import centroid_com
from scipy import optimize as opt

start_time = time.time()

# directories
FOLDER_DIR_DATA = "data/"  # data folder
CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder
FIELD_DIR = FOLDER_DIR_DATA + "fields/"  # galaxy fields
GALAXY_DIR = FOLDER_DIR_DATA + "galaxy_images/"  # galaxy images

# files
MAGPI_PARENT_FILE = "magpi_parent.csv"  # MAGPI parent file

# set constants
RE_LIM = 0.7
ARC_SCALE = 0.2  # Arcsec / Pixel scale
NOISE_IT = 100  # Number of iterations for bacground noise

WALK_LIM_ARCSEC = 0.6  # Limit walk away from the centroid when searching for true centre (arcsec)
WALK_LIM_PIX = WALK_LIM_ARCSEC / ARC_SCALE  # Centre walk limit, expressed in pixels

# Galaxy_Data_Import
GALAXY_DATA_FILE = "galaxy_list.csv"
df_galaxy_list = pd.read_csv(CSV_DIR + GALAXY_DATA_FILE)
del df_galaxy_list["Unnamed: 0"]

arr = [images for images in os.listdir(GALAXY_DIR) if images.endswith(".fits")]
arr_sorted = sorted(arr)


def parent_formula(
    df_galaxy_list: pd.DataFrame,
    GALAXY_DIR: str,
    arr_sorted: list[str],
    i: int,
    WALK_LIM_PIX: float,
    NOISE_IT: int,
):
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

    seg_id = df_galaxy_list["segID"].iloc[i]  # Gets SEG ID for galaxy

    # Transforms the mask such that the mask is equal to 1, background is set to 0 and everything else
    # (other galaxies) are nans.
    mask = mask.astype(float)
    mask[(mask != seg_id) & (mask != 0)] = np.nan
    mask[mask == seg_id] = 1

    # Creates a copy of the mask called segmap and then a masked imaged of the galaxy
    segmap = mask.copy()
    gal_image = image * segmap

    # Creates a masked error array expressed in variance
    var_image = im_var * segmap

    centroid = centroid_com(gal_image)  # Centroid of galaxy
    xc, yc = centroid  # Centroid coordinates

    # Function to determine the minimum asymmetry of the galaxy for optimisation
    # The function is embedded in the code rather than called upon as the opt.fmin function fails when taking
    # inputs other than that which is beign varied

    def asym_calc_opt(centre):
        x, y = centre

        if (
            np.sqrt(abs(xc - x) ** 2 + abs(yc - y) ** 2) <= WALK_LIM_PIX
        ):  # Limits the walk of the centre of the galaxy to be less than walk_limit_pixel from the centroid
            # to prevent the optimiser from using NaNs to minimse the asymmetry
            gal_image_180 = skimage.transform.rotate(
                gal_image, 180.0, center=centre
            )  # rotate galaxy image about stated centre point
            rotated_residual = np.abs(
                gal_image - gal_image_180
            )  # Residuals from subtrating the rotated image from the original

            # Creates a mask so that when normalising against the original galaxy image both gal_image and
            # rotated residuals are summed over the same number of pixels. I.e makes sure there are the same
            # number of NaNs in gal_image and rotated_residuals
            gal_image_div_mask = rotated_residual.copy()
            gal_image_div_mask[~np.isnan(gal_image_div_mask)] = 1
            gal_image_divide = gal_image_div_mask * gal_image

            P_asym_uncorrected = np.nansum(rotated_residual) / np.nansum(gal_image_divide)

        else:
            P_asym_uncorrected = np.inf

        return P_asym_uncorrected

    # the optimised position of the segmap which acts as the point of rotation minimising the asymmetry index
    centre_asym = opt.fmin(asym_calc_opt, centroid)
    xa, ya = centre_asym

    # Finds how far the centre of asymmetry is from the centroid of the galaxy
    walk_dist_pix = np.sqrt(abs(xc - xa) ** 2 + abs(yc - ya) ** 2)

    # Get the background asymmetry value
    noise_asym, err_noise_asym, back_noise_sum, back_noise_arr_avg = noise_correction(
        NOISE_IT, mask, image, centre_asym, gal_image
    )

    # Calculate the P_asym (with errors) before background correction
    P_asym_uncorrected, err_rel_uncorrected = asym_calc_uncorrected(centre_asym, gal_image, var_image)

    # Background Corrected P_asym
    P_asym = P_asym_uncorrected - noise_asym
    P_asym_err = P_asym * np.sqrt(
        (err_rel_uncorrected / P_asym_uncorrected) ** 2 + (err_noise_asym / noise_asym) ** 2
    )

    # Calculate Concentration
    C, C_err = Concentration(
        df_galaxy_list, centre_asym, gal_image, var_image, back_noise_sum, mask, ARC_SCALE, i
    )

    # Count number of nans that exist in resdiual image (i.e. after image has been rotated and subtracted)
    # This indicative of how much of the image is covered in nans, if too high then data could be distorted

    gal_image_180 = skimage.transform.rotate(
        gal_image, 180.0, center=centre_asym
    )  # rotate galaxy image about stated centre point
    rotated_residual = np.abs(
        gal_image - gal_image_180
    )  # Residuals from subtrating the rotated image from the original

    image_size = rotated_residual.shape[0] * rotated_residual.shape[1]  # Image size
    nan_count = np.count_nonzero(np.isnan(rotated_residual))  # Count of nans in image

    nan_frac = nan_count / image_size  # Fraction of image that is nans
    nan_frac = np.round(nan_frac, 2)  # Fraction to 2 decimal places

    return (
        P_asym,
        P_asym_err,
        C,
        C_err,
        centre_asym,
        centroid,
        walk_dist_pix,
        nan_frac,
    )


P_asym_arr = []
P_asym_err_arr = []
C_arr = []
C_err_arr = []
centre_asym_x_arr = []
centre_asym_y_arr = []
centroid_x_arr = []
centroid_y_arr = []
walk_dist_pix_arr = []
nan_frac_arr = []

for i in range(0, len(arr_sorted)):
    print(arr_sorted[i])

    # Extract Data
    val_hold = parent_formula(df_galaxy_list, GALAXY_DIR, arr_sorted, i, WALK_LIM_PIX, NOISE_IT)

    # Save Data to Array
    P_asym_arr.append(val_hold[0])
    P_asym_err_arr.append(val_hold[1])
    C_arr.append(val_hold[2])
    C_err_arr.append(val_hold[3])
    centre_asym_x_arr.append(val_hold[4][0])
    centre_asym_y_arr.append(val_hold[4][1])
    centroid_x_arr.append(val_hold[5][0])
    centroid_y_arr.append(val_hold[5][0])
    walk_dist_pix_arr.append(val_hold[6])
    nan_frac_arr.append(val_hold[7])

df_galaxy_list_det = df_galaxy_list.copy()

# Place data in df_Galaxy_List
df_galaxy_list_det["P_asym"] = P_asym_arr
df_galaxy_list_det["P_asym_err"] = P_asym_err_arr
df_galaxy_list_det["C"] = C_arr
df_galaxy_list_det["C_err"] = C_err_arr
df_galaxy_list_det["centre_asym_x"] = centre_asym_x_arr
df_galaxy_list_det["centre_asym_y"] = centre_asym_y_arr
df_galaxy_list_det["centroid_x"] = centroid_x_arr
df_galaxy_list_det["centroid_y"] = centroid_y_arr
df_galaxy_list_det["walk_dist_pix"] = walk_dist_pix_arr
df_galaxy_list_det["nan_frac"] = nan_frac_arr

# Upadate and print new Galaxy List information to csv
filepath = CSV_DIR + "galaxy_list_phot.csv"
df_galaxy_list_det.to_csv(filepath)

print("--- %s seconds ---" % (time.time() - start_time))
