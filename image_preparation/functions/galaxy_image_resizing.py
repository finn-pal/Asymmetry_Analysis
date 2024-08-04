import math

import pandas as pd
from astropy.io import fits


def image_size_lim(r: float) -> float:
    """
    Summary:
        Formula for the minimum resizing based on analysis completed from StatMorph investigation.

    Args:
        r (float): Radius of galaxy (units: pixels).

    Returns:
        float: Minimum sizing of image for accurate galaxy analysis (units: pixels).
    """
    n = 0.0002 * r**3 - 0.0684 * r**2 + 8.2745 * r + 3.3452
    return n


def get_galaxy_image(
    gal_idx: int,
    df_field_galaxy_info: pd.DataFrame,
    hdu1,
    psf,
    hdu2,
    image_g,
    stat_g,
    hdu3,
    image_i,
    stat_i,
    hdu4,
    image_r,
    stat_r,
    hdu5,
    image_w,
    hdu6,
    image_mask_dilated,
    hdu7,
    image_mask_undilated,
    GALAXY_DIR,
):
    mag_id = df_field_galaxy_info["MAGPIID"].iloc[gal_idx]  # Gets MAGPI ID
    seg_id = df_field_galaxy_info["segID_x"].iloc[gal_idx]  # Gets SED ID

    # everything that is not the seg_id of the galaxy becomes a zero, this is for sizing of the dilated mask
    mask_sizing = (image_mask_dilated == seg_id).astype(int)

    x_len = mask_sizing.shape[0]  # x size of image
    y_len = mask_sizing.shape[1]  # y size of image

    # create empty place holders
    x_min = []
    x_max = []
    y_min = []
    y_max = []

    # x_min position of mask
    for idx in range(0, x_len):
        x_min.append(max(mask_sizing[:, idx]))

        if max(x_min) != 0:
            break

    # x_max position of mask
    for idx in range(0, x_len):
        x_max.append(max(mask_sizing[:, -1 - idx]))

        if max(x_max) != 0:
            break

    # y_min position of mask
    for idx in range(0, y_len):
        y_min.append(max(mask_sizing[idx, :]))

        if max(y_min) != 0:
            break

    # y_max position of mask
    for idx in range(0, y_len):
        y_max.append(max(mask_sizing[-1 - idx, :]))

        if max(y_max) != 0:
            break

    # x and y Size of mask
    seg_x_len = x_len - len(x_max) - len(x_min)
    seg_y_len = y_len - len(y_max) - len(y_min)

    Sizing = max(seg_x_len, seg_y_len)  # maximum dimension of dilated mask

    # for the relevant galaxy obtain the x and y dimensions of the image
    xc = df_field_galaxy_info["xcubecen"].iloc[gal_idx]
    x1 = df_field_galaxy_info["x1"].iloc[gal_idx]
    x2 = df_field_galaxy_info["x2"].iloc[gal_idx]

    yc = df_field_galaxy_info["ycubecen"].iloc[gal_idx]
    y1 = df_field_galaxy_info["y1"].iloc[gal_idx]
    y2 = df_field_galaxy_info["y2"].iloc[gal_idx]

    d_x = x2 - x1  # x sizing of current image
    d_y = y2 - y1  # y sizing of current image

    # r50 of galaxy in pixels
    r50_pixel = df_field_galaxy_info.loc[
        df_field_galaxy_info.index[df_field_galaxy_info["MAGPIID"] == mag_id][0], "r50_pixel"
    ]

    min_sizing = image_size_lim(r50_pixel)

    # takes the largest number between the two sizing methods and converts to an integer
    delta_x = math.ceil(max(d_x, Sizing * 2, min_sizing) / 2) * 2

    # the length from the centre to the x-dim edge of the image
    delta_x = int(delta_x / 2)

    # Check the increase in image size does not go out of the bounds of the reduced image
    if xc + delta_x > image_g.shape[0]:
        delta_x = image_g.shape[0] - xc

    elif xc - delta_x < 0:
        delta_x = xc

    # takes the largest number between the two sizing methods and converts to an integer
    delta_y = math.ceil(max(d_y, Sizing * 2, min_sizing) / 2) * 2

    # the length from the centre to the x-dim edge of the image
    delta_y = int(delta_y / 2)

    # check the increase in image size does not go out of the bounds of the reduced image
    if yc + delta_y > image_g.shape[1]:
        delta_y = image_g.shape[1] - yc

    elif yc - delta_y < 0:
        delta_y = yc

    # new sizing limits on image
    df_field_galaxy_info.at[gal_idx, "x1_d"] = df_field_galaxy_info["xcubecen"][gal_idx] - delta_x
    df_field_galaxy_info.at[gal_idx, "x2_d"] = df_field_galaxy_info["xcubecen"][gal_idx] + delta_x
    df_field_galaxy_info.at[gal_idx, "y1_d"] = df_field_galaxy_info["ycubecen"][gal_idx] - delta_y
    df_field_galaxy_info.at[gal_idx, "y2_d"] = df_field_galaxy_info["ycubecen"][gal_idx] + delta_y

    # conversion to integers for processing
    x1_d = int(df_field_galaxy_info.at[gal_idx, "x1_d"])
    x2_d = int(df_field_galaxy_info.at[gal_idx, "x2_d"])
    y1_d = int(df_field_galaxy_info.at[gal_idx, "y1_d"])
    y2_d = int(df_field_galaxy_info.at[gal_idx, "y2_d"])

    # image resizing
    image_g_resize = image_g[y1_d:y2_d, x1_d:x2_d]
    image_i_resize = image_i[y1_d:y2_d, x1_d:x2_d]
    image_r_resize = image_r[y1_d:y2_d, x1_d:x2_d]
    image_w_resize = image_w[y1_d:y2_d, x1_d:x2_d]
    image_mask_dilated_resize = image_mask_dilated[y1_d:y2_d, x1_d:x2_d]
    image_mask_undilated_resize = image_mask_undilated[y1_d:y2_d, x1_d:x2_d]
    stat_g_resize = stat_g[y1_d:y2_d, x1_d:x2_d]
    stat_i_resize = stat_i[y1_d:y2_d, x1_d:x2_d]
    stat_r_resize = stat_r[y1_d:y2_d, x1_d:x2_d]

    # Conversion of data to fits data type
    res_psf = fits.ImageHDU(psf, name="psf", header=hdu1[4].header)
    res_im_g = fits.ImageHDU(image_g_resize, name="g_mod_band", header=hdu2[1].header)
    res_st_g = fits.ImageHDU(stat_g_resize, name="g_mod_band_stat", header=hdu2[2].header)
    res_im_i = fits.ImageHDU(image_i_resize, name="i_band", header=hdu3[1].header)
    res_st_i = fits.ImageHDU(stat_i_resize, name="i_0band_stat", header=hdu3[2].header)
    res_im_r = fits.ImageHDU(image_r_resize, name="r_band", header=hdu4[1].header)
    res_st_r = fits.ImageHDU(stat_r_resize, name="r_band_stat", header=hdu4[2].header)
    res_im_w = fits.ImageHDU(image_w_resize, name="white", header=hdu5[1].header)
    res_md = fits.ImageHDU(image_mask_dilated_resize, name="dilated_mask", header=hdu6[0].header)
    res_mu = fits.ImageHDU(image_mask_undilated_resize, name="undilated_mask", header=hdu7[0].header)

    hdu1[0].header["MAGPIID"] = mag_id

    newhdul = fits.HDUList(
        [
            hdu1[0],
            res_psf,
            res_im_g,
            res_st_g,
            res_im_i,
            res_st_i,
            res_im_r,
            res_st_r,
            res_im_w,
            res_md,
            res_mu,
        ]
    )

    # Export to new fits file

    newhdul.writeto(GALAXY_DIR + "/MAGPI" + str(mag_id) + "_resized" + ".fits", overwrite=True)


def resize_from_field(
    df_galaxy_list: pd.DataFrame, field_id: int, FIELD_DIR: str, GALAXY_DIR: str, ARC_SCALE: float = 0.2
):
    # Relevant frequency band images to import
    IMAGE_BAND_G = "gmod_SDSS"  # Gmod
    IMAGE_BAND_I = "i_SDSS"  # iBand
    IMAGE_BAND_R = "r_SDSS"  # rBand
    IMAGE_BAND_W = "CollapsedImage"  # WhiteBand
    IMAGE_MASK_DIL = "manual_segmap"  # undilated manual segmap
    IMAGE_MASK_UNDIL = "manual_segmap_undilated"  # dilated manual segmap

    # open relevant fie;d
    field_loc = FIELD_DIR + "MAGPI" + str(field_id) + "/MAGPI" + str(field_id)

    # psf
    hdu1 = fits.open(field_loc + ".fits")
    psf = hdu1[4].data

    # g_mod_band
    hdu2 = fits.open(field_loc + "_" + IMAGE_BAND_G + ".fits")
    image_g = hdu2[1].data
    stat_g = hdu2[2].data

    # i_band
    hdu3 = fits.open(field_loc + "_" + IMAGE_BAND_I + ".fits")
    image_i = hdu3[1].data
    stat_i = hdu3[2].data

    # r_band
    hdu4 = fits.open(field_loc + "_" + IMAGE_BAND_R + ".fits")
    image_r = hdu4[1].data
    stat_r = hdu4[2].data

    # white image
    hdu5 = fits.open(field_loc + "_" + IMAGE_BAND_W + ".fits")
    image_w = hdu5[1].data

    # dilated Mask
    hdu6 = fits.open(field_loc + "_" + IMAGE_MASK_DIL + ".fits")
    image_mask_dilated = hdu6[0].data

    # undilated Mask
    hdu7 = fits.open(field_loc + "_" + IMAGE_MASK_UNDIL + ".fits")
    image_mask_undilated = hdu7[0].data

    # create relevant data frames
    df_edges = pd.read_csv(field_loc + "_segedges.csv")  # import data on current sizing of galaxies
    df_stats = pd.read_csv(field_loc + "_segstats.csv")  # import data on current positon of galaxies

    # combine the above two data frames into one
    df_field_galaxy_info = df_edges.merge(df_stats, left_on="MAGPIID", right_on="MAGPIID")
    df_field_galaxy_info["r50_pixel"] = df_field_galaxy_info["R50"] / ARC_SCALE

    # remove galaxies from df_combine that are not galaxies of interest to be investigated
    cond = ~df_field_galaxy_info["MAGPIID"].isin(df_galaxy_list["MAGPIID"])
    df_field_galaxy_info.drop(df_field_galaxy_info[cond].index, inplace=True)
    df_field_galaxy_info = df_field_galaxy_info.reset_index(drop=True)

    df_field_galaxy_info["x1_d"] = pd.Series(dtype="int")
    df_field_galaxy_info["x2_d"] = pd.Series(dtype="int")
    df_field_galaxy_info["y1_d"] = pd.Series(dtype="int")
    df_field_galaxy_info["y2_d"] = pd.Series(dtype="int")

    for gal_idx in range(0, len(df_field_galaxy_info)):
        # creates resized gaalxy images from a given field
        get_galaxy_image(
            gal_idx,
            df_field_galaxy_info,
            hdu1,
            psf,
            hdu2,
            image_g,
            stat_g,
            hdu3,
            image_i,
            stat_i,
            hdu4,
            image_r,
            stat_r,
            hdu5,
            image_w,
            hdu6,
            image_mask_dilated,
            hdu7,
            image_mask_undilated,
            GALAXY_DIR,
        )


def get_fields_list(df_galaxy_list: pd.DataFrame) -> pd.DataFrame:
    """
    Summary:
        Create a dataframe containing all the relevant galaxy fields to analyse.

    Args:
        df_galaxy_list (pd.DataFrame): Trimmed list of galaxies to analyse.
        ARC_SCALE (float, optional): MUSE pixel to arsecond scale. Defaults to 0.2.

    Returns:
        pd.DataFrame: Dataframe containing all the relevant galaxy fields to analyse.
    """

    # Creates list of reduced cubes of interest from refined list of galaxies
    field_list = df_galaxy_list.iloc[:, 0]
    field_list = field_list.astype(str).str[:4]
    field_list = field_list.drop_duplicates()
    df_field_list = pd.DataFrame(data=field_list)
    df_field_list = df_field_list.reset_index(drop=True)
    df_field_list.columns = ["field_id"]

    return df_field_list


def run_galaxy_extraction(
    df_galaxy_list: pd.DataFrame, FIELD_DIR: str, GALAXY_DIR: str, ARC_SCALE: float = 0.2
):
    df_field_list = get_fields_list(df_galaxy_list)

    for field in df_field_list["field_id"]:
        print("Images extracted from", field)
        resize_from_field(df_galaxy_list, field, FIELD_DIR, GALAXY_DIR, ARC_SCALE)
