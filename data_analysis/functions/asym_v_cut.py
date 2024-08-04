import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from matplotlib.ticker import MaxNLocator
from scipy import stats


def v_asym_cut(df_gal, cut_lim=0.8):
    label_lst = ["Gas", "Stellar"]
    code_lst = ["v_asym_g", "v_asym_s"]
    cut_lst = []

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    fig.subplots_adjust(wspace=0.2)

    for i in range(0, 2):
        data = np.array(df_gal[(~(df_gal[code_lst[i]]).isnull())][code_lst[i]])
        data = sorted(data)

        shape_loc, loc_param, scale_param = stats.lognorm.fit(data)

        max_data = np.round(np.max(data) + 0.02, 2)

        sample_lin = np.linspace(0, max_data, 400)
        pdf = sp.stats.lognorm.pdf(sample_lin, shape_loc, loc=loc_param, scale=scale_param)

        sample_dist = stats.lognorm.pdf(data, shape_loc, loc=loc_param, scale=scale_param)
        cut = stats.lognorm.ppf(cut_lim, shape_loc, loc_param, scale_param)

        bins = np.linspace(0, 0.17, 18)

        axs[i].hist(data, density=True, bins=bins, color="silver", edgecolor="black", ls="--", linewidth=0.5)
        axs[i].plot(sample_lin, pdf, "g-", label="Log Nomral Fit to the Distribution")
        axs[i].fill_between(
            sample_lin,
            pdf,
            where=(sample_lin <= cut),
            edgecolor="green",
            color="green",
            alpha=0.15,
            label="Probability interval of: " + str(np.round(cut_lim * 100, 2)) + "%",
        )
        axs[i].plot(
            [cut, cut],
            [-1, 50],
            c="r",
            ls="dashdot",
            marker=None,
            label="Asymmetry Threshold at: " + str(np.round(cut, 3)),
            lw=1.2,
        )

        axs[i].set_xlabel("Kinematic Asymmetry")
        axs[i].set_ylim([0, np.max(sample_dist) + 2])
        axs[i].yaxis.set_major_locator(MaxNLocator(integer=True))
        axs[i].tick_params("y", labelleft=True)
        axs[i].set_yticks(np.linspace(0, 40, 9))
        axs[i].set_title(label_lst[i] + " Kinematics", fontsize=12)

        axs[i].set_xticks(np.linspace(0, 0.175, 8))
        axs[i].set_xlim([0, 0.175])

        axs[i].set_ylabel("Bin Count")

        axs[i].legend(loc="upper right", fontsize=10, framealpha=1)

        axs[i].tick_params(labelsize=12)

        cut_lst.append(cut)

    vg_asym_Flag = [
        np.nan if v_f == 0 else 1 if cut_lst[0] < v_a else 0
        for v_a, v_f in zip(df_gal["v_asym_g"], df_gal["Kin_Flag_g"])
    ]

    vs_asym_Flag = [
        np.nan if v_f == 0 else 1 if cut_lst[1] < v_a else 0
        for v_a, v_f in zip(df_gal["v_asym_s"], df_gal["Kin_Flag_s"])
    ]

    df_gal["vg_asym_Flag"] = vg_asym_Flag
    df_gal["vs_asym_Flag"] = vs_asym_Flag

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "v_asym_threshold.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    return df_gal
