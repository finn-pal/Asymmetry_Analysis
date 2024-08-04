import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull

# Asym Param #############################################


def plot_asym_cont(df_gal, md, plt_error=True):
    ew = 0.6
    lw = 0.3
    z = 0
    # ms = 30
    fs = 8
    x_offset = 0.007
    y_offset = 0.003

    fig, axs = plt.subplots(nrows=2, ncols=1)

    for i in range(0, 6):
        #################################### Gas ####################################
        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 1)]

        x = np.array(df_plot.P_asym)
        x_err = np.array(df_plot.P_asym_err)
        y = np.array(df_plot.v_asym_g)
        y_err = np.array(df_plot.v_asym_g_err)

        axs[0].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=0.1, s=10)
        if plt_error:
            axs[0].errorbar(x, y, xerr=x_err, yerr=y_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

        x_hull = np.array([xi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])
        y_hull = np.array([yi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])

        length = len(x_hull)

        if length >= 3:
            points = np.concatenate([x_hull, y_hull]).reshape((2, length)).T
            hull = ConvexHull(points)

            axs[0].add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_g == 0)]
        x = np.array(df_plot.P_asym)
        x_err = np.array(df_plot.P_asym_err)
        y = np.array([y_offset] * len(x))
        axs[0].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=lw, s=12, marker="v")
        if plt_error:
            axs[0].errorbar(x, y, xerr=x_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_g == 1)]
        y = np.array(df_plot.v_asym_g)
        y_err = np.array(df_plot.v_asym_g_err)
        x = np.array([x_offset] * len(y))
        axs[0].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=lw, s=12, marker="<")
        if plt_error:
            axs[0].errorbar(x, y, yerr=y_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

        ################################### Stars ###################################
        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_s == 1)]

        x = np.array(df_plot.P_asym)
        x_err = np.array(df_plot.P_asym_err)
        y = np.array(df_plot.v_asym_s)
        y_err = np.array(df_plot.v_asym_s_err)

        axs[1].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=0.1, s=10)
        if plt_error:
            axs[1].errorbar(x, y, xerr=x_err, yerr=y_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

        x_hull = np.array([xi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])
        y_hull = np.array([yi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])

        length = len(x_hull)

        if length >= 3:
            points = np.concatenate([x_hull, y_hull]).reshape((2, length)).T
            hull = ConvexHull(points)

            axs[1].add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 1) & (df_gal.Kin_Flag_s == 0)]
        x = np.array(df_plot.P_asym)
        x_err = np.array(df_plot.P_asym_err)
        y = np.array([y_offset] * len(x))
        axs[1].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=lw, s=12, marker="v")
        if plt_error:
            axs[1].errorbar(x, y, xerr=x_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

        df_plot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 0) & (df_gal.Kin_Flag_s == 1)]
        y = np.array(df_plot.v_asym_s)
        y_err = np.array(df_plot.v_asym_s_err)
        x = np.array([x_offset] * len(y))
        axs[1].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=lw, s=12, marker="<")
        if plt_error:
            axs[1].errorbar(x, y, yerr=y_err, color="grey", fmt="none", elinewidth=ew, zorder=z)

    ######################## Formatting ##############################
    fig.subplots_adjust(hspace=0, wspace=0)

    axs[0].set_ylabel("Gas Kinematic Asymmetry", fontsize=fs)
    axs[0].set_xticklabels([])
    axs[0].set_yticks(np.linspace(0.02, 0.14, 7))

    axs[1].set_xlabel("Photometric Asymmetry", fontsize=fs)
    axs[1].set_ylabel("Stellar Kinematic Asymmetry", fontsize=fs)
    axs[1].set_yticks(np.linspace(0.0, 0.14, 8))

    for i in range(0, 2):
        axs[i].tick_params(right=True, top=True, direction="inout")
        axs[i].set(xlim=[0, 0.7], ylim=[0, 0.14])
        axs[i].set_xticks(np.linspace(0.0, 0.70, 8))

        axs[i].tick_params(axis="both", labelsize=fs)

    ######################## Legend ##############################

    # create legend

    handles, _ = plt.gca().get_legend_handles_labels()

    for i in range(0, 6):
        if i != 2:
            # create manual symbols for legend
            point = Line2D(
                [0],
                [0],
                label=md["mk_n"][i],
                marker=".",
                markersize=7,
                markeredgecolor=md["mk_ec"][i],
                markerfacecolor=md["mk_c"][i],
                linestyle="",
                mew=0.1,
            )

            handles.append(point)

    axs[1].legend(handles=handles, fontsize=fs - 1, loc="upper right", frameon=False)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "asym_cont.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


# S05 Linear Fit #########################################


def fit_func(x, m, b):
    return m * x + b


def lin_fit(Mass, S05):
    popt, pcov = curve_fit(fit_func, np.log10(Mass), np.log10(S05))
    m, b = popt
    err_cov_log = np.sqrt(np.diag(pcov))
    # m, b = np.polyfit(np.log10(Mass), np.log10(S05), 1)

    b_err_log = err_cov_log[1]
    b_err = 10 ** (b_err_log)

    Mass_fit = np.arange(3 * 10**9, 2 * 10**12, 10**10)

    S05_log_fit = m * np.log10(Mass_fit) + b
    S05_fit = 10 ** (S05_log_fit)

    print(
        r"\log_{10}(S_{0.5}) = "
        + "("
        + str(np.round(m, 2))
        + r" \pm "
        + str(np.round(err_cov_log[0], 2))
        + ")"
        + r"\log_{10}(M_{*} / M_{\odot})"
        + " + ("
        + str(np.round(b, 2))
        + r" \pm "
        + str(np.round(err_cov_log[1], 2))
        + ")"
    )

    return Mass_fit, S05_fit, b_err


# S05 ####################################################


def plot_s05_new(df_gal, md):
    z = 0
    fs = 10
    ew = 0.2

    #### QUESTIONNNNNNNS WHAT ABOUT NANS

    # calculate s05 and error
    df_gal["S05_g"] = np.sqrt(0.5 * (df_gal.V_rot_g) ** 2 + (df_gal.Sigma_g) ** 2)
    df_gal["S05_s"] = np.sqrt(0.5 * (df_gal.V_rot_s) ** 2 + (df_gal.Sigma_s) ** 2)

    Error_arr_s = []
    Error_arr_g = []

    Default_s = np.nanmean(df_gal["v_asym_s_err"] / df_gal["v_asym_s"])
    Default_g = np.nanmean(df_gal["v_asym_g_err"] / df_gal["v_asym_g"])

    for i in range(0, len(df_gal)):
        if (~np.isnan(df_gal["v_asym_s"][i])) & (~np.isnan(df_gal["v_asym_s_err"][i])):
            err_rel_s = df_gal["v_asym_s_err"][i] / df_gal["v_asym_s"][i]

        else:
            err_rel_s = Default_s

        if (~np.isnan(df_gal["v_asym_g"][i])) & (~np.isnan(df_gal["v_asym_g_err"][i])):
            err_rel_g = df_gal["v_asym_g_err"][i] / df_gal["v_asym_g"][i]

        else:
            err_rel_g = Default_g

        Error_arr_s.append(err_rel_s)
        Error_arr_g.append(err_rel_g)

    # plot s05 - mass relation

    df_gal["S05_g_err"] = Error_arr_g * df_gal["S05_g"]
    df_gal["S05_s_err"] = Error_arr_s * df_gal["S05_s"]

    fig, axs = plt.subplots(2, 1, figsize=(8, 8), sharey=True, sharex=True)

    for i in range(0, 6):
        # kinematics of the gas
        df_gal_g = df_gal[(df_gal.Kin_Flag_g == 1) & (df_gal.int_stage == i)]

        xg = (df_gal_g.StellarMass_median).tolist()
        xg_el = (df_gal_g.StellarMass_median - df_gal_g.StellarMass_16).tolist()
        xg_eu = (df_gal_g.StellarMass_84 - df_gal_g.StellarMass_median).tolist()

        yg = (df_gal_g.S05_g).tolist()
        yg_e = (df_gal_g.S05_g_err).tolist()

        # kinematics of the stars
        df_gal_s = df_gal[(df_gal.Kin_Flag_s == 1) & (df_gal.int_stage == i)]

        xs = (df_gal_s.StellarMass_median).tolist()
        xs_el = (df_gal_s.StellarMass_median - df_gal_s.StellarMass_16).tolist()
        xs_eu = (df_gal_s.StellarMass_84 - df_gal_s.StellarMass_median).tolist()

        ys = (df_gal_s.S05_s).tolist()
        ys_e = (df_gal_s.S05_s_err).tolist()

        axs[0].errorbar(xg, yg, xerr=[xg_el, xg_eu], yerr=yg_e, color="grey", fmt="none", lw=ew, zorder=z)
        axs[0].scatter(xg, yg, c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=10)

        axs[1].errorbar(xs, ys, xerr=[xs_el, xs_eu], yerr=ys_e, color="grey", fmt="none", lw=ew, zorder=z)
        axs[1].scatter(xs, ys, c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=10)

    print("\n########################")
    print("S05 Fit - Gas")

    xg = (df_gal[(df_gal.Kin_Flag_g == 1) & (df_gal.int_stage == 0)].StellarMass_median).tolist()
    yg = (df_gal[(df_gal.Kin_Flag_g == 1) & (df_gal.int_stage == 0)].S05_g).tolist()
    xfit_g, yfit_g, b_err_g = lin_fit(xg, yg)
    axs[0].plot(xfit_g, yfit_g, c=md["mk_c"][0], linestyle="-.")
    axs[0].fill_between(
        xfit_g, yfit_g - b_err_g, yfit_g + b_err_g, fc="tab:blue", alpha=0.1, ec=None, zorder=z
    )

    print()
    print("S05 Fit - Stars")

    xs = (df_gal[(df_gal.Kin_Flag_s == 1) & (df_gal.int_stage == 0)].StellarMass_median).tolist()
    ys = (df_gal[(df_gal.Kin_Flag_s == 1) & (df_gal.int_stage == 0)].S05_s).tolist()
    xfit_s, yfit_s, b_err_s = lin_fit(xs, ys)
    axs[1].plot(xfit_s, yfit_s, c=md["mk_c"][0], linestyle="-.")
    axs[1].fill_between(
        xfit_s, yfit_s - b_err_s, yfit_s + b_err_s, fc="tab:blue", alpha=0.1, ec=None, zorder=z
    )

    print("########################")

    # formatting
    for i in range(0, 2):
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")

        axs[i].tick_params(which="both", right=True, top=True, direction="inout")

    axs[0].set_ylabel("S05$_{Gas}$ (km/s)", fontsize=fs)
    axs[1].set_ylabel("S05$_{Stars}$ (km/s)", fontsize=fs)
    axs[1].set_xlabel(r"M$_{*}$ / M$_{\odot}$", fontsize=fs)

    # Remove space between axes
    fig.subplots_adjust(hspace=0, wspace=0)

    handles, _ = plt.gca().get_legend_handles_labels()

    for i in range(0, 6):
        if i != 2:
            # create manual symbols for legend
            point = Line2D(
                [0],
                [0],
                label=md["mk_n"][i],
                marker=".",
                markersize=7,
                markeredgecolor=md["mk_ec"][i],
                markerfacecolor=md["mk_c"][i],
                linestyle="",
                mew=0.1,
            )

            handles.append(point)

    line = Line2D(
        [0],
        [0],
        label="S05 Undisturbed",
        marker=None,
        color="tab:blue",
        linestyle="-.",
        lw=1,
    )

    handles.append(line)

    axs[0].legend(handles=handles, loc="upper center", fontsize=9, ncol=3)

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "s05.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    return df_gal


def find_s05_dist(df_gal):
    # Gas #######################################################
    df_king = df_gal[(df_gal.Kin_Flag_g == 1)]

    xg0 = np.array((df_king[df_king.int_stage == 0].StellarMass_median).tolist())
    yg0 = np.array((df_king[df_king.int_stage == 0].S05_g).tolist())

    popt, _ = curve_fit(fit_func, np.log10(xg0), np.log10(yg0))
    m, b = popt

    print("\n########################")
    print("S05 Dist - Gas")
    for i in range(0, 6):
        xg = np.array((df_king[df_king.int_stage == i].StellarMass_median))
        yg = np.array((df_king[df_king.int_stage == i].S05_g))

        yg_fit_log = np.log10(xg) * m + b
        yg_fit = 10 ** (yg_fit_log)
        abs_dist_g = np.abs(yg - yg_fit)

        if len(yg_fit) > 3:
            g_m = np.mean(abs_dist_g)
            g_e = np.std(abs_dist_g) / np.sqrt(len(abs_dist_g))
            print("int stage: " + str(i) + " - " + str(np.round(g_m, 2)) + " +/- " + str(np.round(g_e, 2)))

    # Stars #######################################################
    df_kins = df_gal[(df_gal.Kin_Flag_s == 1)]

    xs0 = np.array((df_kins[df_kins.int_stage == 0].StellarMass_median).tolist())
    ys0 = np.array((df_kins[df_kins.int_stage == 0].S05_s).tolist())

    popt, _ = curve_fit(fit_func, np.log10(xs0), np.log10(ys0))
    m, b = popt

    print()
    print("S05 Dist - Stars")
    for i in range(0, 6):
        xs = np.array((df_kins[df_kins.int_stage == i].StellarMass_median))
        ys = np.array((df_kins[df_kins.int_stage == i].S05_s))

        ys_fit_log = np.log10(xs) * m + b
        ys_fit = 10 ** (ys_fit_log)
        abs_dist_s = np.abs(ys - ys_fit)

        if len(ys_fit) > 3:
            s_m = np.mean(abs_dist_s)
            s_e = np.std(abs_dist_s) / np.sqrt(len(abs_dist_s))
            print("int stage: " + str(i) + " - " + str(np.round(s_m, 2)) + " +/- " + str(np.round(s_e, 2)))

    print("########################")


# SFR New ################################################


def sfr_new(df_gal, md):
    # lw = 1
    z = 0
    # fs = 9

    df_sf = df_gal[(df_gal.StarForm_ID_50 == 1) & (df_gal.SFR_50 >= 10 ** (-3))]

    plt.figure()

    for i in range(0, 6):
        df_plot = df_sf[(df_sf.int_stage == i)]

        x = df_plot.StellarMass_median
        y = df_plot.SFR_50

        plt.scatter(x, y, c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=10)

    x_full = df_sf[df_sf.SFR_50 >= 10 ** (-3)].StellarMass_median

    x_lin = np.linspace(np.min(x_full), np.max(x_full), 100)
    y_lin = [x_l ** (0.712) * (10 ** (-7.293)) for x_l in x_lin]
    y_lin_max = [x_l ** (0.712 + 0.04) * (10 ** (-7.293 + 0.06)) for x_l in x_lin]
    y_lin_min = [x_l ** (0.712 - 0.04) * (10 ** (-7.293 - 0.06)) for x_l in x_lin]

    plt.plot(x_lin, y_lin, c="magenta", ls="dashdot", zorder=z)
    plt.fill_between(x_lin, y_lin_min, y_lin_max, fc="magenta", alpha=0.1, ec=None, zorder=z)

    plt.xscale("log")
    plt.yscale("log")

    plt.xlim([2 * 10**9, 10**12])

    plt.xlabel(r"M$_{*}$ / M$_{\odot}$")
    plt.ylabel("SFR")

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "sfr.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


# SFR Conc ###############################################


def sfr_conc_new(df_gal, md):
    lw = 0.5
    # z = 0
    # fs = 9

    df_sf = df_gal[(df_gal.StarForm_ID_50 == 1) & (df_gal.SFR_50 >= 10 ** (-3))]

    plt.figure()

    for i in range(0, 6):
        df_plot = df_sf[df_sf.int_stage == i]

        ha_20 = df_plot.Ha_20
        ha_50 = df_plot.Ha_50

        y = ha_20 / ha_50
        x = df_plot.StellarMass_median
        # x = df_plot.n_use

        if len(df_plot) > 3:
            y_mean = np.mean(y)
            y_std_err = np.std(y) / np.sqrt(len(y))
            y1 = y_mean + y_std_err
            y2 = y_mean - y_std_err
            # plt.axhspan(y1, y2, color=md["mk_c"][i], alpha=0.1, lw=0)
            plt.plot([2 * 10**9, 10**12], [y_mean, y_mean], c=md["mk_c"][i], ls="--", lw=lw)

        plt.scatter(x, y, c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=10)
        plt.xscale("log")
        plt.xlim([2 * 10**9, 10**12])

        plt.xlabel(r"M$_{*}$ / M$_{\odot}$")
        plt.ylabel("H\u03b1$_{20}$ / H\u03b1$_{50}$")

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "sfr_conc.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


# Mass v Sersic ##########################################


def mass_sersic(df_gal, md):
    lw = 1
    # z = 0
    # fs = 9

    # plt.figure()
    fig, ax = plt.subplots()

    for i in range(0, 6):
        df_plot = df_gal[df_gal.int_stage == i]

        y = df_plot.n_use
        x = df_plot.StellarMass_median

        ax.scatter(x, y, c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=10)

        x_hull = np.array([xi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])
        y_hull = np.array([yi for xi, yi in zip(x, y) if (~np.isnan(xi)) & (~np.isnan(yi))])

        length = len(x_hull)

        if length >= 3:
            points = np.concatenate([x_hull, y_hull]).reshape((2, length)).T
            hull = ConvexHull(points)

            ax.add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

    ax.set_xscale("log")
    # plt.xlim([2 * 10**9, 10**12])

    ax.set_xlabel(r"M$_{*}$ / M$_{\odot}$")
    ax.set_ylabel("Sersic")

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "mass_sersic.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


# Mass Histogram #########################################


def mass_hist(df_gal, md):
    _, ax = plt.subplots()

    print()
    for i in range(0, 6):
        print(md["mk_n"][i] + ": " + str(len(df_gal[df_gal.int_stage == i])))
    print()

    mass = df_gal.StellarMass_median
    logbins = np.logspace(np.log10(np.min(mass)) - 1, np.log10(np.max(mass)) + 1, 10)

    for i in range(0, 4):
        if i != 2:
            mass_i = df_gal[df_gal.int_stage == i].StellarMass_median

            ax.hist(
                mass_i,
                bins=logbins,
                color=md["mk_c"][i],
                fill=False,
                hatch=md["mk_hs"][i],
                edgecolor=md["mk_ec"][i],
                label=md["mk_n"][i],
                alpha=md["mk_al"][i],
                zorder=i,
            )

    ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    ax.set_xscale("log")
    ax.set_xlim([10**9, 10**13])
    ax.set_xlabel(r"M$_{*}$ / M$_{\odot}$")
    ax.set_ylabel("Bin Count")
    ax.legend()

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "mass_int_hist.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


# Mass Dependence ##########################################


def plot_mass_dep_int(df_gal, md, plt_error=True):
    z = 0
    ew = 0.5

    fig, axs = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    fig.subplots_adjust(hspace=0, wspace=0.2)

    # Photometry ################################################
    for i in range(0, 6):
        df_phot = df_gal[(df_gal.int_stage == i) & (df_gal.Phot_Flag == 1)]

        y = df_phot.P_asym
        x = (df_phot.StellarMass_median).tolist()
        axs[0].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=0.1, s=10)

        if plt_error:
            y_err = df_phot.P_asym_err
            x_el = (df_phot.StellarMass_median - df_phot.StellarMass_16).tolist()
            x_eu = (df_phot.StellarMass_84 - df_phot.StellarMass_median).tolist()
            axs[0].errorbar(x, y, xerr=[x_el, x_eu], yerr=y_err, color="grey", fmt="none", lw=ew, zorder=z)

        length = len(df_phot)

        if (length >= 3) & (i < 4):
            points = np.concatenate([x, y]).reshape((2, length)).T
            hull = ConvexHull(points)

            axs[0].add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

    # Gas Kinematics ################################################
    for i in range(0, 6):
        df_king = df_gal[(df_gal.int_stage == i) & (df_gal.Kin_Flag_g == 1)]

        y = df_king.v_asym_g
        x = (df_king.StellarMass_median).tolist()
        axs[1].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=0.1, s=10)

        if plt_error:
            y_err = df_king.P_asym_err
            x_el = (df_king.StellarMass_median - df_king.StellarMass_16).tolist()
            x_eu = (df_king.StellarMass_84 - df_king.StellarMass_median).tolist()
            axs[1].errorbar(x, y, xerr=[x_el, x_eu], yerr=y_err, color="grey", fmt="none", lw=ew, zorder=z)

        length = len(df_king)

        if (length >= 3) & (i < 4):
            points = np.concatenate([x, y]).reshape((2, length)).T
            hull = ConvexHull(points)

            axs[1].add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

    # Stellar Kinematics ##############################################
    for i in range(0, 6):
        df_kins = df_gal[(df_gal.int_stage == i) & (df_gal.Kin_Flag_s == 1)]

        y = df_kins.v_asym_s
        x = (df_kins.StellarMass_median).tolist()
        axs[2].scatter(x, y, c=md["mk_c"][i], ec=md["mk_ec"][i], lw=0.1, s=10)

        if plt_error:
            y_err = df_kins.P_asym_err
            x_el = (df_kins.StellarMass_median - df_kins.StellarMass_16).tolist()
            x_eu = (df_kins.StellarMass_84 - df_kins.StellarMass_median).tolist()
            axs[2].errorbar(x, y, xerr=[x_el, x_eu], yerr=y_err, color="grey", fmt="none", lw=ew, zorder=z)

        length = len(df_kins)

        if (length >= 3) & (i < 4):
            points = np.concatenate([x, y]).reshape((2, length)).T
            hull = ConvexHull(points)

            axs[2].add_collection(
                PolyCollection(
                    [points[hull.vertices, :]],
                    edgecolors=md["mk_c"][i],
                    facecolors=md["mk_c"][i],
                    alpha=0.2,
                    linewidths=2,
                    zorder=-1,
                )
            )

    # Formatting #########################################

    ylabel = ["Photometric Asymmetry", "Kinematic Asymmetry (Gas)", "Kinematic Asymmetry (Stars)"]

    axs[2].set_xlabel(r"M$_{*}$ / M$_{\odot}$")

    for i in range(0, 3):
        axs[i].set(xscale="log", yscale="log", ylabel=ylabel[i])
        axs[i].tick_params(which="both", right=True, top=True, direction="inout", length=3)

    for i in range(1, 3):
        axs[i].set_ylim([8 * 10**-3, 2 * 10**-1])

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "mass_dep_int.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################
