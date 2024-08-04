import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def fit_func(x, a, b, c):
    return a * np.exp(-b * x) + c


def p_asym_cut(df_gal, sig_cut=3):
    ############# Failed fit to photometric asymmetry

    fs = 18

    p_asym_dict = {
        "n_02": df_gal[
            (df_gal.n_galfit > 0)
            & (df_gal.n_galfit <= 2)
            & (~(df_gal.n_galfit).isnull())
            & (df_gal.Phot_Flag == 1)
        ].P_asym,
        "n_24": df_gal[
            (df_gal.n_galfit > 2)
            & (df_gal.n_galfit <= 4)
            & (~(df_gal.n_galfit).isnull())
            & (df_gal.Phot_Flag == 1)
        ].P_asym,
        "n_46": df_gal[
            (df_gal.n_galfit > 4)
            & (df_gal.n_galfit <= 6)
            & (~(df_gal.n_galfit).isnull())
            & (df_gal.Phot_Flag == 1)
        ].P_asym,
        "n_68": df_gal[
            (df_gal.n_galfit > 6)
            & (df_gal.n_galfit <= 8)
            & (~(df_gal.n_galfit).isnull())
            & (df_gal.Phot_Flag == 1)
        ].P_asym,
    }

    n_base = df_gal[(~(df_gal.n_galfit).isnull()) & (df_gal.Phot_Flag == 1)].n_galfit
    p_base = df_gal[(~(df_gal.n_galfit).isnull()) & (df_gal.Phot_Flag == 1)].P_asym

    n_fail = df_gal[(~(df_gal.n_fit).isnull()) & (df_gal.Phot_Flag == 1)].n_fit
    p_fail = df_gal[(~(df_gal.n_fit).isnull()) & (df_gal.Phot_Flag == 1)].P_asym

    bins = np.linspace(0, 8, 5)

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))
    fig.subplots_adjust(hspace=0)
    plt.rcParams.update({"font.size": fs - 5})

    box_yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    box_ylabels = ["", 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    # hist_yticks = np.linspace(0, 22, 12)
    hist_yticks = [0, 4, 8, 12, 16, 20]

    axs[0].boxplot(
        p_asym_dict.values(),
        patch_artist=True,
        widths=0.5,
        medianprops=dict(color="blue", linewidth=2),
        boxprops=dict(facecolor="silver"),
        flierprops={"markerfacecolor": "r", "markeredgecolor": "r", "markersize": 5},
    )

    axs[0].set_xticks([])
    axs[0].set_yticks(box_yticks)
    axs[0].set_yticklabels(box_ylabels)
    axs[0].set_ylim([0, 0.7])

    axs[1].hist(n_base, bins, edgecolor="black", ls="--", color="w")
    axs[1].set_xticks([0, 2, 4, 6, 8])
    axs[1].set_xlim([0, 8])
    axs[1].set_ylim([0, 20])
    axs[1].set_yticks(hist_yticks)
    axs[1].invert_yaxis()

    axs[1].set_xlabel("S$\\acute{e}$rsic Index (n)")

    axs[0].set_ylabel("Photometric Asymmetry")
    axs[1].set_ylabel("Bin Count")

    for ax in axs:
        ax.xaxis.label.set_size(fs)
        ax.yaxis.label.set_size(fs)
        ax.tick_params(axis="both", which="both", labelsize=fs - 3, length=5, width=1.5)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "outlier_identification.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    upper_lst = []
    for arr in p_asym_dict.values():
        # finding the 1st quartile
        q1 = np.quantile(arr, 0.25)

        # finding the 3rd quartile
        q3 = np.quantile(arr, 0.75)

        # finding the iqr region
        iqr = q3 - q1

        # finding upper and lower whiskers
        upper_bound = q3 + (1.5 * iqr)
        upper_lst.append(upper_bound)

    data_n = []
    data_p = []
    out_n = []
    out_p = []
    for n, p in zip(n_base, p_base):
        if n <= 2:
            if p < upper_lst[0]:
                data_n.append(n)
                data_p.append(p)
            else:
                out_n.append(n)
                out_p.append(p)

        elif n <= 4:
            if p < upper_lst[1]:
                data_n.append(n)
                data_p.append(p)
            else:
                out_n.append(n)
                out_p.append(p)

        elif n <= 6:
            if p < upper_lst[2]:
                data_n.append(n)
                data_p.append(p)
            else:
                out_n.append(n)
                out_p.append(p)

        elif n <= 8:
            if p < upper_lst[3]:
                data_n.append(n)
                data_p.append(p)
            else:
                out_n.append(n)
                out_p.append(p)

    data_p_sort = [p for _, p in sorted(zip(data_n, data_p))]
    data_n_sort = sorted(data_n)

    popt, pcov = curve_fit(fit_func, data_n_sort, data_p_sort)
    lw = 1.2
    err_cov = np.sqrt(np.diag(pcov))

    x_new = np.linspace(0, 10, 200)
    y_new = [fit_func(xi, *popt) for xi in x_new]

    plt.figure(figsize=(8, 4.5))
    plt.scatter(data_n_sort, data_p_sort, s=10, c="b", label="Not Obviously Asymmetric")
    plt.scatter(out_n, out_p, s=10, c="r", label="Obvious Outliers")

    plt.scatter(
        n_fail,
        p_fail,
        s=22,
        marker="d",
        c="w",
        edgecolors="g",
        lw=1.2,
        label="Estimated S$\\acute{e}$rsic Indexes",
    )

    plt.plot(x_new, y_new, "b", lw=lw, linestyle="dashdot", zorder=0)
    plt.plot(
        x_new,
        y_new + sig_cut * err_cov[2],
        c="r",
        linestyle="dashdot",
        lw=lw,
        zorder=0,
    )

    plt.xticks(np.linspace(0, 10, 11))
    plt.yticks(np.linspace(0, 0.7, 8))
    plt.xlim([0, 10])
    plt.ylim([0, 0.7])

    plt.xlabel("S$\\acute{e}$rsic Index (n)")
    plt.ylabel("Photometric Asymmetry")
    plt.legend(fontsize=10, markerscale=1.3)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "p_asym_threshold.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    P_asym_Flag = [
        1 if p > (fit_func(n, *popt) + sig_cut * err_cov[2]) else 0
        for p, n in zip(df_gal["P_asym"], df_gal["n_use"])
    ]

    df_gal["P_asym_Flag"] = P_asym_Flag

    print("\n######################\n")
    print("Photometric Asymmetry Threshold")
    print(
        "p_{asym} = "
        + "("
        + str(np.round(popt[0], 2))
        + r" \pm "
        + str(np.round(err_cov[0], 2))
        + ") "
        + "e^{"
        + "-("
        + str(np.round(popt[1], 2))
        + r" \pm "
        + str(np.round(err_cov[1], 2))
        + " ) n} + "
        + "("
        + str(np.round(popt[2], 2))
        + r" \pm "
        + str(np.round(err_cov[2], 2))
        + ") "
    )
    print("\n######################\n")
    # a * np.exp(-b * x) + c

    return df_gal
