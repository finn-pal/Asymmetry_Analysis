import time

import pandas as pd
from functions.asym_contour_plot import (
    find_s05_dist,
    mass_hist,
    mass_sersic,
    plot_asym_cont,
    plot_mass_dep_int,
    plot_s05_new,
    sfr_conc_new,
    sfr_new,
)
from functions.asym_p_cut import p_asym_cut
from functions.asym_v_cut import v_asym_cut
from functions.conc_p_asym import conc_p_asym
from functions.interaction_stage import int_stage, int_stage_new

# from functions.marker_dict import md
from functions.marker_dict_new import md
from functions.plot_mass_dep import plot_mass_dep
from functions.plot_sfr import plot_sfr_outliers
from functions.sample_sizing import sample_sizing_new
from functions.sersic_cut import sersic_divide
from functions.visual_add import visual_add

start_time = time.time()


# directories
FOLDER_DIR_DATA = "data/"  # data folder
CSV_DIR = FOLDER_DIR_DATA + "csv_data_files/"  # csv data folder

# Galaxy_Data_Import
GALAXY_DATA_FILE = "galaxy_parent.csv"
df_parent = pd.read_csv(CSV_DIR + GALAXY_DATA_FILE)
del df_parent["Unnamed: 0"]

VISUAL_INSPECTION_FILE = "visual_class.csv"
df_vis = pd.read_csv(CSV_DIR + VISUAL_INSPECTION_FILE)

df_child = df_parent[(df_parent.Kin_Flag_s == 1) | (df_parent.Kin_Flag_g == 1) | (df_parent.Phot_Flag == 1)]
df_child = df_child.reset_index(drop=True)

df_child = visual_add(df_child)
df_child = df_child[df_child.edge_case == 0]
df_child = df_child.reset_index(drop=True)

sample_sizing_new(df_child)

df_child = sersic_divide(df_child)
df_child = p_asym_cut(df_child, sig_cut=3)  # sigma cut = 3
df_child = v_asym_cut(df_child, cut_lim=0.841)  # 1 sigma upper limit

int_stage(df_child)
df_child = int_stage_new(df_child)
plot_mass_dep(df_child)

plot_asym_cont(df_child, md, plt_error=False)
df_child = plot_s05_new(df_child, md)
find_s05_dist(df_child)
sfr_new(df_child, md)
sfr_conc_new(df_child, md)
plot_sfr_outliers(df_child)
mass_sersic(df_child, md)
# vis_check(df_child)
mass_hist(df_child, md)
plot_mass_dep_int(df_child, md, plt_error=False)

conc_p_asym(df_child)

# plt.show()

# filepath = "dataframe_print/" + "df_data_child.csv"
# filepath = FOLD
FILE_NAME = "galaxy_child.csv"
df_child.to_csv(CSV_DIR + FILE_NAME)

print("--- %s seconds ---" % (time.time() - start_time))
