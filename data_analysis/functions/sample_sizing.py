import matplotlib.pyplot as plt
from matplotlib_venn import venn3_circles, venn3_unweighted


def sample_sizing(df_gal, md):
    tot_count = len(df_gal)

    pho_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 0)])
    gas_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 0)])
    ste_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 1)])

    pho_g_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 0)])
    pho_s_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 1)])
    g_s_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 1)])

    pho_g_s_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 1)])

    print("###########################")
    print()
    print("Number of Galaxies in Sample: " + str(tot_count))
    print()
    print("Photometry Only: " + str(pho_co))
    print("Gas Kinematics Only: " + str(gas_co))
    print("Stellar Kinematics Only: " + str(ste_co))
    print()
    print("Photometry and Gas Kinematics: " + str(pho_g_co))
    print("Photometry and Stellar Kinematics: " + str(pho_s_co))
    print("Gas and Stellar Kinematics: " + str(g_s_co))
    print()
    print("Photometry and Gas and Stellar Kinematics " + str(pho_g_s_co))
    print()
    print("###########################")
    print()

    for i in range(0, 6):
        print(str(md["mk_n"][i]) + ": " + str(len(df_gal[df_gal.int_stage == i])))

    print()
    print("###########################")

    labels = ("Gas Kinematics", "Stellar Kinematics", "Photometry")

    plt.figure()

    subsets = (gas_co, ste_co, g_s_co, pho_co, pho_g_co, pho_s_co, pho_g_s_co)

    v = venn3_unweighted(subsets, set_labels=labels, set_colors=("w", "w", "w"))
    venn3_circles((1, 1, 1, 1, 1, 1, 1), linestyle="-", linewidth=0.5, color="k")

    lbl = v.get_label_by_id("C")
    x, y = lbl.get_position()
    lbl.set_position((x, y - 0.05))

    for text in v.set_labels:
        text.set_fontsize(11)

    for text in v.subset_labels:
        text.set_fontsize(11)

    plt.gca().invert_yaxis()

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "venn_sample.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


def sample_sizing_new(df_gal):
    tot_count = len(df_gal)

    pho_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 0)])
    gas_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 0)])
    ste_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 1)])

    pho_g_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 0)])
    pho_s_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 0) & (df_gal.Kin_Flag_s == 1)])
    g_s_co = len(df_gal[(df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 1)])

    pho_g_s_co = len(df_gal[(df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 1) & (df_gal.Kin_Flag_s == 1)])

    print("###########################")
    print()
    print("Number of Galaxies in Sample: " + str(tot_count))
    print()
    print("Photometry Only: " + str(pho_co))
    print("Gas Kinematics Only: " + str(gas_co))
    print("Stellar Kinematics Only: " + str(ste_co))
    print()
    print("Photometry and Gas Kinematics: " + str(pho_g_co))
    print("Photometry and Stellar Kinematics: " + str(pho_s_co))
    print("Gas and Stellar Kinematics: " + str(g_s_co))
    print()
    print("Photometry and Gas and Stellar Kinematics " + str(pho_g_s_co))
    print()
    print("###########################")

    labels = ("Gas Kinematics", "Stellar Kinematics", "Photometry")

    plt.figure()

    subsets = (gas_co, ste_co, g_s_co, pho_co, pho_g_co, pho_s_co, pho_g_s_co)

    v = venn3_unweighted(subsets, set_labels=labels, set_colors=("w", "w", "w"))
    venn3_circles((1, 1, 1, 1, 1, 1, 1), linestyle="-", linewidth=0.5, color="k")

    lbl = v.get_label_by_id("C")
    x, y = lbl.get_position()
    lbl.set_position((x, y - 0.05))

    for text in v.set_labels:
        text.set_fontsize(11)

    for text in v.subset_labels:
        text.set_fontsize(11)

    plt.gca().invert_yaxis()

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "venn_sample.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################
