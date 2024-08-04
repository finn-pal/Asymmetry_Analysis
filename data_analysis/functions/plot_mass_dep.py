import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy import stats


def plot_mass_dep(df_gal, plt_error=False):
    z = 0
    ew = 0.5
    fs = 18

    custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "custom", ["darkblue", "blue", "magenta", "red", "crimson", "darkred"]
    )
    cmap = custom_cmap

    # photometry
    df_phot = df_gal[(df_gal.Phot_Flag == 1) & (~(df_gal.n_galfit).isnull())]
    yp = df_phot.P_asym
    yp_err = df_phot.P_asym_err
    cp = df_phot.n_galfit

    xp = (df_phot.StellarMass_median).tolist()
    xp_el = (df_phot.StellarMass_median - df_phot.StellarMass_16).tolist()
    xp_eu = (df_phot.StellarMass_84 - df_phot.StellarMass_median).tolist()

    # gas kinematics
    df_king = df_gal[(df_gal.Kin_Flag_g == 1) & (~(df_gal.n_galfit).isnull())]
    yg = df_king.v_asym_g
    yg_err = df_king.v_asym_g_err
    cg = df_king.n_galfit

    xg = (df_king.StellarMass_median).tolist()
    xg_el = (df_king.StellarMass_median - df_king.StellarMass_16).tolist()
    xg_eu = (df_king.StellarMass_84 - df_king.StellarMass_median).tolist()

    # stellar kinematics
    df_kins = df_gal[(df_gal.Kin_Flag_s == 1) & (~(df_gal.n_galfit).isnull())]
    ys = df_kins.v_asym_s
    ys_err = df_kins.v_asym_s_err
    cs = df_kins.n_galfit

    xs = (df_kins.StellarMass_median).tolist()
    xs_el = (df_kins.StellarMass_median - df_kins.StellarMass_16).tolist()
    xs_eu = (df_kins.StellarMass_84 - df_kins.StellarMass_median).tolist()

    fig, axs = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    fig.subplots_adjust(hspace=0, wspace=0.2)
    plt.rcParams.update({"font.size": fs - 5})

    c_plot = axs[0].scatter(xp, yp, c=cp, cmap=cmap, s=10, vmin=0, vmax=8)
    axs[1].scatter(xg, yg, c=cg, cmap=cmap, s=10, vmin=0, vmax=8)
    axs[2].scatter(xs, ys, c=cs, cmap=cmap, s=10, vmin=0, vmax=8)

    if plt_error:
        axs[0].errorbar(xp, yp, xerr=[xp_el, xp_eu], yerr=yp_err, color="grey", fmt="none", lw=ew, zorder=z)
        axs[1].errorbar(xg, yg, xerr=[xg_el, xg_eu], yerr=yg_err, color="grey", fmt="none", lw=ew, zorder=z)
        axs[2].errorbar(xs, ys, xerr=[xs_el, xs_eu], yerr=ys_err, color="grey", fmt="none", lw=ew, zorder=z)

    ################ Failed Sersic Fits

    # photometry
    df_phot = df_gal[(df_gal.Phot_Flag == 1) & (~(df_gal.n_fit).isnull())]
    yp = df_phot.P_asym
    yp_err = df_phot.P_asym_err
    cp = df_phot.n_use

    xp = (df_phot.StellarMass_median).tolist()
    xp_el = (df_phot.StellarMass_median - df_phot.StellarMass_16).tolist()
    xp_eu = (df_phot.StellarMass_84 - df_phot.StellarMass_median).tolist()

    # gas kinematics
    df_king = df_gal[(df_gal.Kin_Flag_g == 1) & (~(df_gal.n_fit).isnull())]
    yg = df_king.v_asym_g
    yg_err = df_king.v_asym_g_err
    cg = df_king.n_use

    xg = (df_king.StellarMass_median).tolist()
    xg_el = (df_king.StellarMass_median - df_king.StellarMass_16).tolist()
    xg_eu = (df_king.StellarMass_84 - df_king.StellarMass_median).tolist()

    # stellar kinematics
    df_kins = df_gal[(df_gal.Kin_Flag_s == 1) & (~(df_gal.n_fit).isnull())]
    ys = df_kins.v_asym_s
    ys_err = df_kins.v_asym_s_err
    cs = df_kins.n_use

    xs = (df_kins.StellarMass_median).tolist()
    xs_el = (df_kins.StellarMass_median - df_kins.StellarMass_16).tolist()
    xs_eu = (df_kins.StellarMass_84 - df_kins.StellarMass_median).tolist()

    axs[0].scatter(xp, yp, c="w", edgecolor="k", s=12)
    axs[1].scatter(xg, yg, c="w", edgecolor="k", s=12)
    axs[2].scatter(xs, ys, c="w", edgecolor="k", s=12)

    if plt_error:
        axs[0].errorbar(xp, yp, xerr=[xp_el, xp_eu], yerr=yp_err, color="grey", fmt="none", lw=ew, zorder=z)
        axs[1].errorbar(xg, yg, xerr=[xg_el, xg_eu], yerr=yg_err, color="grey", fmt="none", lw=ew, zorder=z)
        axs[2].errorbar(xs, ys, xerr=[xs_el, xs_eu], yerr=ys_err, color="grey", fmt="none", lw=ew, zorder=z)

    ######### PEARSON RANK STATS ###########
    df_phot = df_gal[(df_gal.Phot_Flag == 1)]
    yp = df_phot.P_asym
    xp = df_phot.StellarMass_median

    df_king = df_gal[(df_gal.Kin_Flag_g == 1)]
    yg = df_king.v_asym_g
    xg = df_king.StellarMass_median

    df_kins = df_gal[(df_gal.Kin_Flag_s == 1)]
    ys = df_kins.v_asym_s
    xs = df_kins.StellarMass_median

    res_p = stats.pearsonr(np.log(xp), np.log(yp))
    res_g = stats.pearsonr(np.log(xg), np.log(yg))
    res_s = stats.pearsonr(np.log(xs), np.log(ys))

    res = [res_p, res_g, res_s]

    for i in range(0, 3):
        axs[i].text(
            0.85,
            0.8,
            "r = " + str(round(res[i][0], 2)) + "\n" + r"$\rho$ = " + str(round(res[i][1], 2)),
            transform=axs[i].transAxes,
            fontsize=12,
            ha="left",
            color="k",
        )

    # Formatting
    cb_ax = fig.add_axes([0.91, 0.12, 0.015, 0.754])
    fig.colorbar(c_plot, orientation="vertical", cax=cb_ax, label="S$\\acute{e}$rsic Index")

    ylabel = ["Photometric\nAsymmetry", "Kinematic\nAsymmetry (Gas)", "Kinematic\nAsymmetry (Stars)"]

    axs[2].set_xlabel(r"M$_{*}$ / M$_{\odot}$", fontsize=fs - 3)

    for i in range(0, 3):
        axs[i].set(xscale="log", yscale="log")
        axs[i].set_ylabel(ylabel[i], fontsize=fs - 5)
        axs[i].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        axs[i].tick_params(axis="both", which="both", labelsize=fs - 5)

    for i in range(1, 3):
        axs[i].set_ylim([0.008, 0.2])
        axs[i].set_yticks([0.01, 0.1])
        axs[i].tick_params(which="both", top=True, direction="inout", length=3)

    axs[0].tick_params(which="both", bottom=True, direction="inout", length=3)

    axs[0].set_ylim([0.008, 0.8])

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "mass_dep.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################
