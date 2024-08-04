import numpy as np
import pandas as pd


def get_kinematics(df_gal):
    # directories
    FOLDER_DIR_DATA = "data/"  # data folder
    CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder

    # Galaxy_Data_Import
    GALAXY_KINEMATICS_DATA_FILE = "MAGPI_kinemetry_sample_s05.csv"
    df_kins = pd.read_csv(CSV_DIR + GALAXY_KINEMATICS_DATA_FILE)

    # Holding list to fill with information to then be added to dataframe

    v_asym_g_arr = []
    v_asym_g_err_arr = []
    V_rot_g_arr = []
    SNR_g_arr = []
    Sigma_g_arr = []

    v_asym_s_arr = []
    v_asym_s_err_arr = []
    V_rot_s_arr = []
    SNR_s_arr = []
    Sigma_s_arr = []

    for i in range(0, len(df_gal)):
        try:
            ID = df_gal["MAGPIID"][i]  # Get MAGPI Id for ith element in dataframe
            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx = df_kins.index[df_kins["MAGPIID"] == ID].tolist()

            # Find corresponding value
            v_asym_g = df_kins["v_asym_g"][idx].iloc[0]
            v_asym_g_err = df_kins["v_asym_g_err"][idx].iloc[0]
            V_rot_g = df_kins["V_rot_g"][idx].iloc[0]
            SNR_g = df_kins["SNR_g"][idx].iloc[0]
            Sigma_g = df_kins["Sigma_g"][idx].iloc[0]

            v_asym_s = df_kins["v_asym_s"][idx].iloc[0]
            v_asym_s_err = df_kins["v_asym_s_err"][idx].iloc[0]
            V_rot_s = df_kins["V_rot_s"][idx].iloc[0]
            SNR_s = df_kins["SNR_s"][idx].iloc[0]
            Sigma_s = df_kins["Sigma_s"][idx].iloc[0]

            # Add value to holding array
            v_asym_g_arr.append(v_asym_g)
            v_asym_g_err_arr.append(v_asym_g_err)
            V_rot_g_arr.append(V_rot_g)
            SNR_g_arr.append(SNR_g)
            Sigma_g_arr.append(Sigma_g)

            v_asym_s_arr.append(v_asym_s)
            v_asym_s_err_arr.append(v_asym_s_err)
            V_rot_s_arr.append(V_rot_s)
            SNR_s_arr.append(SNR_s)
            Sigma_s_arr.append(Sigma_s)

        except:  # noqa: E722
            # If no correspondign value found then set to nan
            v_asym_g_arr.append(np.nan)
            v_asym_g_err_arr.append(np.nan)
            V_rot_g_arr.append(np.nan)
            SNR_g_arr.append(np.nan)
            Sigma_g_arr.append(np.nan)

            v_asym_s_arr.append(np.nan)
            v_asym_s_err_arr.append(np.nan)
            V_rot_s_arr.append(np.nan)
            SNR_s_arr.append(np.nan)
            Sigma_s_arr.append(np.nan)

    # Add arrays to the df_Galaxy_List dataframe
    df_gal["v_asym_g"] = v_asym_g_arr
    df_gal["v_asym_g_err"] = v_asym_g_err_arr
    df_gal["V_rot_g"] = V_rot_g_arr
    df_gal["SNR_g"] = SNR_g_arr
    df_gal["Sigma_g"] = Sigma_g_arr

    df_gal["v_asym_s"] = v_asym_s_arr
    df_gal["v_asym_s_err"] = v_asym_s_err_arr
    df_gal["V_rot_s"] = V_rot_s_arr
    df_gal["SNR_s"] = SNR_s_arr
    df_gal["Sigma_s"] = Sigma_s_arr

    return df_gal
