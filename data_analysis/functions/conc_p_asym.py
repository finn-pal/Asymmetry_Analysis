import matplotlib.pyplot as plt


def conc_p_asym(df_gal):
    fs = 20

    _, ax = plt.subplots(figsize=(6, 6))
    plt.rcParams.update({"font.size": fs - 5})

    # x = df_gal[df_gal.Phot_Flag == 1].C
    # y = df_gal[df_gal.Phot_Flag == 1].P_asym

    x_base = df_gal[~(df_gal.n_galfit).isnull()].n_galfit
    y_base = df_gal[~(df_gal.n_galfit).isnull()].P_asym

    x_fail = df_gal[~(df_gal.n_fit).isnull()].n_fit
    y_fail = df_gal[~(df_gal.n_fit).isnull()].P_asym

    ax.scatter(x_base, y_base, marker=".", c="blue")
    ax.scatter(x_fail, y_fail, marker="d", c="green")

    ax.set_xlim([9.5, 0])
    # ax.set_xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5])

    ax.set_ylim([0.7, 0])
    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])

    ax.xaxis.label.set_size(fs)
    ax.yaxis.label.set_size(fs)
    ax.tick_params(axis="both", which="both", labelsize=fs - 3, length=5, width=1.5)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "conc_p_asym.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################
