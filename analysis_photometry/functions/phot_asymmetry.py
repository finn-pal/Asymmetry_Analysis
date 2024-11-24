import random

import numpy as np
import skimage


def asym_calc_uncorrected(centre_asym, gal_image, var_image):
    # rotate galaxy image about stated centre point
    gal_image_180 = skimage.transform.rotate(gal_image, 180.0, center=centre_asym)

    # Residuals from subtrating the rotated image from the original
    rotated_residual = np.abs(gal_image - gal_image_180)

    # Creates a mask so that when normalising against the original galaxy image both gal_image and rotated
    # residuals are summed over the same number of pixels. I.e makes sure there are the same number of NaNs in
    # gal_image and rotated_residuals
    gal_image_div_mask = rotated_residual.copy()
    gal_image_div_mask[~np.isnan(gal_image_div_mask)] = 1
    gal_image_divide = gal_image_div_mask * gal_image

    # P_asym value that has NOT been background corrected
    P_asym_uncorrected = np.nansum(rotated_residual) / np.nansum(gal_image_divide)

    # Error in this value
    var_image_180 = skimage.transform.rotate(var_image, 180.0, center=centre_asym)  # Rotate error image
    var_residuals = var_image + var_image_180  # Add variances

    # Clean error array removing unwanted values
    var_residuals = var_residuals[(~np.isnan(var_residuals)) & (var_residuals != 0)]

    # Sum residual variances (numerator value)
    var_residuals_sum = np.nansum(var_residuals)

    # Clean variances errors of denomiator
    var_image_divide = gal_image_div_mask * var_image
    var_image_divide = var_image_divide[(~np.isnan(var_image_divide)) & (var_image_divide != 0)]

    # Sum variances errors of denomiator
    var_image_divide_sum = np.nansum(var_image_divide)

    # Add relative errors in quadrature

    err_rel_uncorrected = np.sqrt(
        (var_residuals_sum / (np.nansum(rotated_residual) ** 2))
        + (var_image_divide_sum / (np.nansum(gal_image_divide) ** 2))
    )

    return P_asym_uncorrected, err_rel_uncorrected


def noise_correction(noise_it, mask, image, centre_asym, gal_image):
    # Create an array of just back ground noise from the image. Then fill the mask of the galaxy with this
    # back ground noise and rotate it about the same asym_centre point and get back ground asym. Monte carlo
    # with different iterations of noise.
    noise_mask = mask.copy()
    noise_mask[noise_mask == 1] = np.nan
    noise_mask[noise_mask == 0] = 1

    noise = image * noise_mask  # Get just noise values from original image
    noise = noise[~np.isnan(noise)]  # Remove nans and turn into a 1D array

    gal_mask = mask.copy()
    nx_mask, ny_mask = gal_mask.shape

    noise_arr = []
    back_noise_sum = []
    back_noise_arr = []

    # loop that creates an asymmetry value for the noise that can later be subtracted from the measured value

    for _ in range(noise_it):
        np.random.shuffle(noise)  # Randomly shuffles nosie in selection 'bucket' (i.e. array 'noise')
        noise_data = gal_mask.copy()  # Array to be filled and replaced

        # For loops that sift through each pixel of array and assign a random noise value to the segmap of the
        # galaxy (i.e only places noise where galaxy is)
        for i in range(0, nx_mask):
            for j in range(0, ny_mask):
                if gal_mask[i, j] == 1:
                    # sample w/ replacement (i.e noise values can be chosen more than once)
                    k = random.randint(0, len(noise) - 1)
                    noise_data[i, j] = noise[k]  # Assignment noise value to sample

                else:
                    # Skip pixel locations that are not within the segmap of the galaxy

                    pass

        noise_image = noise_data.copy()

        back_noise_sum.append(np.nansum(noise_image))
        back_noise_arr.append(noise_image)

        # Performs asymmetry on the noise
        # rotate noise image about centre of asymmetry
        noise_image_180 = skimage.transform.rotate(noise_image, 180.0, center=centre_asym)
        noise_residuals = np.abs(noise_image - noise_image_180)

        noise_div_mask = noise_residuals.copy()
        noise_div_mask[~np.isnan(noise_div_mask)] = 1
        noise_divide = noise_div_mask * gal_image

        noise_asym = np.nansum(noise_residuals) / np.nansum(noise_divide)

        noise_arr.append(noise_asym)  # array of nosie asymmetry values

    # Noise asymmetry stats
    noise_asym_val = np.average(noise_arr)  # Average noise asymmetry
    err_noise_asym = np.std(noise_arr)  # Variance of noise asymmetry

    back_noise_arr_avg = np.average(back_noise_arr, 0)

    return noise_asym_val, err_noise_asym, back_noise_sum, back_noise_arr_avg
