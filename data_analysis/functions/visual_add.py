import numpy as np
import pandas as pd


def visual_add(df_gal):
    # directories
    FOLDER_DIR_DATA = "data/"  # data folder
    CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder

    VISUAL_INSPECTION_FILE = "visual_class.csv"
    df_vis = pd.read_csv(CSV_DIR + VISUAL_INSPECTION_FILE)

    vis_class = []
    neighbour = []
    edge_case = []
    obs_int = []

    for i in range(0, len(df_gal)):
        try:
            ID = df_gal["MAGPIID"][i]  # Get MAGPI Id for ith element in dataframe

            # Get dataframe index of this ID in the df_vis_stats dataframe
            idx = df_vis.index[df_vis["MAGPIID"] == ID].tolist()[0]
            vis = df_vis["vis_class"][idx]
            nei = df_vis["neighbours"][idx]
            edg = df_vis["edge"][idx]
            obs = df_vis["obvs_int_neigh"][idx]

            vis_class.append(vis)
            neighbour.append(nei)
            edge_case.append(edg)
            obs_int.append(obs)

        except:  # noqa: E722
            vis_class.append(np.nan)
            neighbour.append(np.nan)
            edge_case.append(np.nan)
            obs_int.append(np.nan)

    df_gal["vis_class"] = vis_class
    df_gal["neighbour"] = neighbour
    df_gal["edge_case"] = edge_case
    df_gal["obs_int"] = obs_int

    return df_gal
