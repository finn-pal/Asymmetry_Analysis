import pandas as pd


def sizing_selection(DF_PARENT: pd.DataFrame, RE_LIM: float = 0.7, MAG_LIM: float = 20.0) -> pd.DataFrame:
    """
    Summary:
        Create a trimmed list of galaxies to analyse. These should not have a radius in arseconds less than
        RE_LIM and should have a radius in the i-band that is larger than the corresponding i-band PSF.

    Args:
        DF_PARENT (pd.DataFrame): MAGPI parent data file containing all galaxies across all fields
        RE_LIM (float, optional): Minimum size of galaxies to be analysed (units: arseconds). Defaults to 0.7.

    Returns:
        pd.DataFrame: Trimmed dataframe.
    """
    df_galaxy_list = DF_PARENT.copy()  # create a copy to be filtered

    # radius limit selection
    df_galaxy_list = df_galaxy_list[df_galaxy_list["re_galfit"] > RE_LIM]

    # only select galaxies with an i band r_50 greater than the FWHM i band seeing.
    df_galaxy_list = df_galaxy_list[df_galaxy_list["R50_it"] > df_galaxy_list["fwhm_i"]]

    # exlude failed galfit fits
    df_galaxy_list = df_galaxy_list[df_galaxy_list["n_galfit"] > 0]

    # remove galaxies with i-band magnitudes greater than 20
    df_galaxy_list = df_galaxy_list[df_galaxy_list["mag_it"] < MAG_LIM]

    # rest dataframe index
    df_galaxy_list = df_galaxy_list.reset_index(drop=True)

    return df_galaxy_list
