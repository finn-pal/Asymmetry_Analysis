import matplotlib.pyplot as plt
import numpy as np


def int_stage(df_gal):
    phot_asym = (df_gal.P_asym_Flag).tolist()
    vg_asym = (df_gal.vg_asym_Flag).tolist()
    vs_asym = (df_gal.vs_asym_Flag).tolist()

    df_f_phot = (df_gal.Phot_Flag).tolist()
    df_f_king = (df_gal.Kin_Flag_g).tolist()
    df_f_kins = (df_gal.Kin_Flag_s).tolist()

    int_arr = []

    for i in range(0, len(df_gal)):
        if df_f_phot[i] == 0:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(5)  # potentially disturbed
            else:
                int_arr.append(4)  # potentially undisturbed

        elif phot_asym[i] == 0:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(3)  # post coalescent phase
            elif ((vg_asym[i] == 0) & (df_f_king[i] == 1)) | ((vs_asym[i] == 0) & (df_f_kins[i] == 1)):
                int_arr.append(0)  # non-interacting
            else:
                int_arr.append(4)  # potentially undisturbed

        elif phot_asym[i] == 1:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(1)  # merger phase
            elif ((vg_asym[i] == 0) & (df_f_king[i] == 1)) | ((vs_asym[i] == 0) & (df_f_kins[i] == 1)):
                int_arr.append(2)  # pair phase
            else:
                int_arr.append(5)  # potentially disturbed

    df_gal["int_stage"] = int_arr

    for i in range(0, 6):
        print(str(i) + ": " + str(len(df_gal[df_gal["int_stage"] == i])))

    return df_gal


def int_stage_new(df_gal):
    phot_asym = (df_gal.P_asym_Flag).tolist()
    vg_asym = (df_gal.vg_asym_Flag).tolist()
    vs_asym = (df_gal.vs_asym_Flag).tolist()

    df_f_phot = (df_gal.Phot_Flag).tolist()
    df_f_king = (df_gal.Kin_Flag_g).tolist()
    df_f_kins = (df_gal.Kin_Flag_s).tolist()

    int_arr = []

    for i in range(0, len(df_gal)):
        if df_f_phot[i] == 0:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(5)  # potentially disturbed
            else:
                int_arr.append(4)  # potentially undisturbed

        elif phot_asym[i] == 0:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(3)  # post coalescent phase
            elif ((vg_asym[i] == 0) & (df_f_king[i] == 1)) | ((vs_asym[i] == 0) & (df_f_kins[i] == 1)):
                int_arr.append(0)  # non-interacting
            else:
                int_arr.append(4)  # potentially undisturbed

        elif phot_asym[i] == 1:
            if ((vg_asym[i] == 1) & (df_f_king[i] == 1)) | ((vs_asym[i] == 1) & (df_f_kins[i] == 1)):
                int_arr.append(1)  # interacting
            elif ((vg_asym[i] == 0) & (df_f_king[i] == 1)) | ((vs_asym[i] == 0) & (df_f_kins[i] == 1)):
                int_arr.append(1)  # interacting
            else:
                int_arr.append(5)  # potentially disturbed

    df_gal["int_stage"] = int_arr

    return df_gal


def mag_int_hist(df_gal, md):
    plt.figure(figsize=(12, 6))
    bins = np.linspace(17, 20, 14)

    for i in range(0, 6):
        data = np.array(df_gal[df_gal.int_stage == i].mag_it)
        plt.hist(
            data,
            color=md["mk_c"][i],
            fill=False,
            hatch=md["mk_hs"][i],
            edgecolor=md["mk_ec"][i],
            bins=bins,
            label=md["mk_n"][i],
            alpha=md["mk_al"][i],
        )

    plt.xlabel("i-Band Magnitude")
    plt.ylabel("Bin Count")

    plt.legend(loc="upper left", fontsize=11)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "interaction_hist.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################
