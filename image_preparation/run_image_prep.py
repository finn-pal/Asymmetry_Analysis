# Imports
import time

import pandas as pd
from functions.galaxy_image_resizing import run_galaxy_extraction
from functions.sizing_selection import sizing_selection

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
MAG_LIM = 20

DF_PARENT = pd.read_csv(CSV_DIR + MAGPI_PARENT_FILE)  # create a pandas dataframe from csv

# create a copy of the magpi parent file which will be filtered to only contain galaxies of interest
df_galaxy_list = sizing_selection(DF_PARENT, RE_LIM)

run_galaxy_extraction(df_galaxy_list, FIELD_DIR, GALAXY_DIR, ARC_SCALE)

# save dataframe containing galaxy list
df_galaxy_list.to_csv(CSV_DIR + "galaxy_list.csv")

print("--- %s seconds ---" % (time.time() - start_time))
