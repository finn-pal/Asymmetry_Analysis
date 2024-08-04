import numpy as np
import pandas as pd
from astropy.io import ascii


def get_sfr(df_gal):
    # directories
    FOLDER_DIR_DATA = "data/"  # data folder
    CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder

    # Galaxy_Data_Import
    SFR_20_DATA_FILE = "SFR_Galaxy_List_20.csv"
    SFR_80_DATA_FILE = "SFR_Galaxy_List_80.csv"
    SFR_50_DATA_FILE = "MAGPI_kinemetry_sample_s05_BPT.csv"

    SFR_TREV = "MAGPI_Balmer_SFRs_onedspec.tbl"

    df_sfr_20 = pd.read_csv(CSV_DIR + SFR_20_DATA_FILE)
    df_sfr_80 = pd.read_csv(CSV_DIR + SFR_80_DATA_FILE)
    df_sfr_50 = pd.read_csv(CSV_DIR + SFR_50_DATA_FILE)

    hold_trev = ascii.read(CSV_DIR + SFR_TREV)
    df_sfr_trev = hold_trev.to_pandas()

    filepath = "dataframe_print/" + "trev.csv"
    df_sfr_trev.to_csv(filepath)

    # Holding list to fill with information to then be added to dataframe
    SFR_trev_arr = []
    SFR_err_trev_arr = []
    SFR_flag_trev_arr = []

    # Holding list to fill with information to then be added to dataframe
    SFR_50_arr = []
    SFR_50_err_arr = []
    SFR_50_Ha_arr = []
    SFR_50_Ha_err_arr = []

    SFR_20_arr = []
    SFR_20_err_arr = []
    SFR_80_arr = []
    SFR_80_err_arr = []

    SFR_20_Ha_arr = []
    SFR_20_Ha_err_arr = []
    SFR_80_Ha_arr = []
    SFR_80_Ha_err_arr = []

    StarForm_ID_20_arr = []
    StarForm_ID_80_arr = []
    StarForm_ID_50_arr = []

    for i in range(0, len(df_gal)):
        try:
            ID = df_gal["MAGPIID"][i]  # Get MAGPI Id for ith element in dataframe

            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx_20 = df_sfr_20.index[df_sfr_20["MAGPIID"] == ID].tolist()[0]
            SFR_20_val = df_sfr_20["SFR"][idx_20]  # SFR_20
            SFR_err_20_val = df_sfr_20["SFR_err"][idx_20]  # SFR_20 error
            Ha_20 = df_sfr_20["Ha"][idx_20]
            Ha_err_20 = df_sfr_20["Ha_err"][idx_20]
            StarForm_ID_20 = df_sfr_20["type(sf=1, sy=2, ln=3) SII"][idx_20]

            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx_80 = df_sfr_80.index[df_sfr_80["MAGPIID"] == ID].tolist()[0]
            SFR_80_val = df_sfr_80["SFR"][idx_80]  # SFR_80
            SFR_err_80_val = df_sfr_80["SFR_err"][idx_80]  # SFR_80 error
            Ha_80 = df_sfr_80["Ha"][idx_80]
            Ha_err_80 = df_sfr_80["Ha_err"][idx_80]
            StarForm_ID_80 = df_sfr_80["type(sf=1, sy=2, ln=3) SII"][idx_80]

            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx_50 = df_sfr_50.index[df_sfr_50["MAGPIID"] == ID].tolist()[0]
            SFR_50_val = df_sfr_50["SFR"][idx_50]  # SFR_80
            SFR_err_50_val = df_sfr_50["SFR_err"][idx_50]  # SFR_80 error
            Ha_50 = df_sfr_50["Ha"][idx_50]
            Ha_err_50 = df_sfr_50["Ha_err"][idx_50]
            StarForm_ID_50 = df_sfr_50["type(sf=1, sy=2, ln=3) SII"][idx_50]

            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx_t = df_sfr_trev.index[df_sfr_trev["MAGPI_ID"] == ID].tolist()[0]
            SFR_t_val = df_sfr_trev["SFR"][idx_t]  # SFR_80
            SFR_err_t_val = df_sfr_trev["SFR_err"][idx_t]  # SFR_80 error50]
            StarForm_ID_t = df_sfr_trev["SFR_flag"][idx_t]

            SFR_20_arr.append(SFR_20_val)  # Add value to holding list
            SFR_20_err_arr.append(SFR_err_20_val)  # Add value to holding list
            SFR_20_Ha_arr.append(Ha_20)  # Add value to holding list
            SFR_20_Ha_err_arr.append(Ha_err_20)  # Add value to holding list

            SFR_80_arr.append(SFR_80_val)  # Add value to holding list
            SFR_80_err_arr.append(SFR_err_80_val)  # Add value to holding list
            SFR_80_Ha_arr.append(Ha_80)  # Add value to holding list
            SFR_80_Ha_err_arr.append(Ha_err_80)  # Add value to holding list

            SFR_50_arr.append(SFR_50_val)  # Add value to holding list
            SFR_50_err_arr.append(SFR_err_50_val)  # Add value to holding list
            SFR_50_Ha_arr.append(Ha_50)  # Add value to holding list
            SFR_50_Ha_err_arr.append(Ha_err_50)  # Add value to holding list

            StarForm_ID_20_arr.append(StarForm_ID_20)
            StarForm_ID_80_arr.append(StarForm_ID_80)
            StarForm_ID_50_arr.append(StarForm_ID_50)

            SFR_trev_arr.append(SFR_t_val)  # Add value to holding list
            SFR_err_trev_arr.append(SFR_err_t_val)  # Add value to holding list
            SFR_flag_trev_arr.append(StarForm_ID_t)

        except:  # noqa: E722
            SFR_20_arr.append(np.nan)
            SFR_20_err_arr.append(np.nan)
            SFR_20_Ha_arr.append(np.nan)
            SFR_20_Ha_err_arr.append(np.nan)

            SFR_80_arr.append(np.nan)
            SFR_80_err_arr.append(np.nan)
            SFR_80_Ha_arr.append(np.nan)
            SFR_80_Ha_err_arr.append(np.nan)

            SFR_50_arr.append(np.nan)
            SFR_50_err_arr.append(np.nan)
            SFR_50_Ha_arr.append(np.nan)
            SFR_50_Ha_err_arr.append(np.nan)

            StarForm_ID_20_arr.append(np.nan)
            StarForm_ID_80_arr.append(np.nan)
            StarForm_ID_50_arr.append(np.nan)

            SFR_trev_arr.append(np.nan)  # Add value to holding list
            SFR_err_trev_arr.append(np.nan)  # Add value to holding list
            SFR_flag_trev_arr.append(np.nan)

    df_gal["SFR_20"] = SFR_20_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["SFR_20_err"] = SFR_20_err_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_20"] = SFR_20_Ha_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_20_err"] = SFR_20_Ha_err_arr  # Add this array to the df_Galaxy_List dataframe

    df_gal["SFR_80"] = SFR_80_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["SFR_80_err"] = SFR_80_err_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_80"] = SFR_80_Ha_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_80_err"] = SFR_80_Ha_err_arr  # Add this array to the df_Galaxy_List dataframe

    df_gal["SFR_50"] = SFR_50_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["SFR_50_err"] = SFR_50_err_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_50"] = SFR_50_Ha_arr  # Add this array to the df_Galaxy_List dataframe
    df_gal["Ha_50_err"] = SFR_50_Ha_err_arr  # Add this array to the df_Galaxy_List dataframe

    df_gal["StarForm_ID_20"] = StarForm_ID_20_arr
    df_gal["StarForm_ID_80"] = StarForm_ID_80_arr
    df_gal["StarForm_ID_50"] = StarForm_ID_50_arr

    df_gal["SFR_trev"] = SFR_trev_arr
    df_gal["SFR_err_trev"] = SFR_err_trev_arr
    df_gal["SFR_flag_trev"] = SFR_flag_trev_arr

    return df_gal
