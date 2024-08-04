import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
from matplotlib.lines import Line2D


def sfr_conc(df_gal, md, sersic_cut=2.5):
    lw = 1
    z = 0
    fs = 9

    n_cut = [0 if sersic <= sersic_cut else 1 for sersic in df_gal.n_use]
    df_gal["n_cut"] = n_cut

    df_sf = df_gal[df_gal.StarForm_ID_80 == 1]

    fig, axs = plt.subplots(3, 2, figsize=(16, 8), sharex=True)

    df_phot = df_sf[(df_sf.Phot_Flag == 1) & (df_sf.Ha_20 != 0)]
    df_gas = df_sf[(df_sf.Kin_Flag_g == 1) & (df_sf.Ha_20 != 0)]
    df_ste = df_sf[(df_sf.Kin_Flag_s == 1) & (df_sf.Ha_20 != 0)]

    # photometry
    for i in range(0, 2):
        df_p = df_phot[df_phot.n_cut == i]
        a0 = axs[0, i]

        for j in range(0, 6):
            df_pj = df_p[df_p.int_stage == j]

            ha_20 = df_pj.Ha_20
            ha_50 = df_pj.Ha_50

            x = ha_20 / ha_50
            y = df_pj.P_asym

            a0.scatter(x, y, marker=md["mk"][j], c=md["mk_c"][j], edgecolor=md["mk_ec"][j], s=md["mk_s"][j])

            if (j == 0) or (j == 1) or (j == 3):
                x_mean = np.mean(x)
                x_std_err = np.std(x) / np.sqrt(len(x))
                xmin = x_mean + x_std_err
                xmax = x_mean - x_std_err
                a0.plot([x_mean, x_mean], [0, 1], c=md["mk_c"][j], ls="--", lw=lw, zorder=z)
                a0.axvspan(xmin, xmax, color=md["mk_c"][j], alpha=0.1, lw=0)

    # gas
    for i in range(0, 2):
        df_g = df_gas[df_gas.n_cut == i]
        a1 = axs[1, i]

        for j in range(0, 6):
            df_gj = df_g[df_g.int_stage == j]

            ha_20 = df_gj.Ha_20
            ha_50 = df_gj.Ha_50

            x = ha_20 / ha_50
            y = df_gj.v_asym_g

            a1.scatter(x, y, marker=md["mk"][j], c=md["mk_c"][j], edgecolor=md["mk_ec"][j], s=md["mk_s"][j])

            if (j == 0) or (j == 1) or (j == 3):
                x_mean = np.mean(x)
                x_std_err = np.std(x) / np.sqrt(len(x))
                xmin = x_mean + x_std_err
                xmax = x_mean - x_std_err
                a1.plot([x_mean, x_mean], [0, 1], c=md["mk_c"][j], ls="--", lw=lw, zorder=z)
                a1.axvspan(xmin, xmax, color=md["mk_c"][j], alpha=0.1, lw=0)

    # stars
    for i in range(0, 2):
        df_s = df_ste[df_ste.n_cut == i]
        a2 = axs[2, i]

        for j in range(0, 6):
            df_sj = df_s[df_s.int_stage == j]

            ha_20 = df_sj.Ha_20
            ha_50 = df_sj.Ha_50

            x = ha_20 / ha_50
            y = df_sj.v_asym_s

            a2.scatter(x, y, marker=md["mk"][j], c=md["mk_c"][j], edgecolor=md["mk_ec"][j], s=md["mk_s"][j])

            if (j == 0) or (j == 1) or (j == 3):
                x_mean = np.mean(x)
                x_std_err = np.std(x) / np.sqrt(len(x))
                xmin = x_mean + x_std_err
                xmax = x_mean - x_std_err
                a2.plot([x_mean, x_mean], [0, 1], c=md["mk_c"][j], ls="--", lw=lw, zorder=z)
                a2.axvspan(xmin, xmax, color=md["mk_c"][j], alpha=0.1, lw=0)

    # formatting
    fig.subplots_adjust(hspace=0, wspace=0)

    ylabels = ["Photometric Asymmetry", "Kinematic Asymmetry \n (Gas)", "Kinematic Asymmetry \n (Stars)"]

    for _, ax in enumerate(axs.flat):
        ax.tick_params(which="both", right=True, top=True, direction="inout", length=3)
        ax.set_xlim([0, 0.8])
        # ax.set(xscale="log", yscale="log")

    for i in range(0, 2):
        axs[0, i].set_ylim([0, 0.55])
        axs[1, i].set_ylim([0, 0.12])
        axs[2, i].set_ylim([0, 0.12])

    for i in range(0, 3):
        axs[i, 1].set_yticklabels([])

    axs[2, 0].set_xlabel("H\u03b1$_{20}$ / H\u03b1$_{50}$")
    axs[2, 1].set_xlabel("H\u03b1$_{20}$ / H\u03b1$_{50}$")

    for i in range(0, 3):
        axs[i, 0].set_ylabel(ylabels[i])

    p_y_ticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    p_y_tick_labels = ["", "0.10", "0.20", "0.30", "0.40", "0.50", "0.60"]
    vg_y_ticks = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
    vg_y_tick_labels = ["", 0.02, 0.04, 0.06, 0.08, "0.10", 0.12]
    vs_y_ticks = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
    vs_y_tick_labels = [0, 0.02, 0.04, 0.06, 0.08, "0.10", 0.12]

    for i in range(0, 2):
        axs[0, i].set_yticks(p_y_ticks)
        axs[0, 0].set_yticklabels(p_y_tick_labels)
        axs[1, i].set_yticks(vg_y_ticks)
        axs[1, 0].set_yticklabels(vg_y_tick_labels)
        axs[2, i].set(yticks=vs_y_ticks)
        axs[2, 0].set_yticklabels(vs_y_tick_labels)

    axs[2, 1].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    axs[2, 1].set_xticklabels(["0.10", "0.20", "0.30", "0.40", "0.50", "0.60", "0.70", "0.80"])

    ####### Text Positioning #########
    n_low = r"n $\leq$ " + str(sersic_cut)
    n_high = str(sersic_cut) + " < n"
    n_label = [n_low, n_high]

    n_ID_colour_samples = ["Low S$\\acute{e}$rsic Sample", "High S$\\acute{e}$rsic Sample"]

    left, width = 0.32, 0.5
    bottom, height = 0.38, 0.5
    right = left + width
    top = bottom + height

    for i in range(0, 2):
        for j in range(0, 3):
            axs[j, i].text(
                right + 0.05,
                top,
                n_label[i] + "\n" + n_ID_colour_samples[i],
                transform=axs[j, i].transAxes,
                fontsize=fs - 1,
                ha="center",
                va="center",
            )

    # create legend
    handles, _ = plt.gca().get_legend_handles_labels()

    for i in range(0, 6):
        # create manual symbols for legend
        point = Line2D(
            [0],
            [0],
            label=md["mk_n"][i],
            marker=md["mk"][i],
            markersize=md["mk_ls"][i],
            markeredgecolor=md["mk_ec"][i],
            markerfacecolor=md["mk_c"][i],
            linestyle="",
        )

        # add manual symbols to auto legend
        handles.append(point)

    axs[0, 0].legend(
        loc="upper center",
        bbox_to_anchor=(1, 1.2),
        ncol=6,
        fancybox=True,
        shadow=True,
        handles=handles,
        fontsize=fs,
    )

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "sfr_conc.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################


def sfr_conc_sersic(df_gal, md):
    lw = 1
    # z = 0
    # fs = 9

    df_sf = df_gal[df_gal.StarForm_ID_50 == 1]

    plt.figure()

    for i in range(0, 6):
        df_plot = df_sf[df_sf.int_stage == i]

        ha_20 = df_plot.Ha_20
        ha_50 = df_plot.Ha_50

        y = ha_20 / ha_50
        x = df_plot.StellarMass_median
        # x = df_plot.n_use

        if (i < 2) or (i == 3):
            y_mean = np.mean(y)
            y_std_err = np.std(y) / np.sqrt(len(y))
            y1 = y_mean + y_std_err
            y2 = y_mean - y_std_err

            plt.axhspan(y1, y2, color=md["mk_c"][i], alpha=0.1, lw=0)
            plt.plot([2 * 10**9, 10**12], [y_mean, y_mean], c=md["mk_c"][i], ls="--", lw=lw)

        plt.scatter(x, y, marker=md["mk"][i], c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=md["mk_s"][i])
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


def sfr(df_gal, md, sersic_cut=2.5):
    # lw = 1
    z = 0
    # fs = 9

    # n_cut = [0 if sersic <= sersic_cut else 1 for sersic in df_gal.n_use]
    # df_gal["n_cut"] = n_cut

    df_sf = df_gal[df_gal.StarForm_ID_50 == 1]

    plt.figure()

    for i in range(0, 6):
        df_plot = df_sf[(df_sf.int_stage == i) & (df_sf.SFR_50 >= 10 ** (-3))]

        x = df_plot.StellarMass_median
        y = df_plot.SFR_50

        plt.scatter(x, y, marker=md["mk"][i], c=md["mk_c"][i], edgecolor=md["mk_ec"][i], s=md["mk_s"][i])

    x_full = df_sf[df_sf.SFR_50 >= 10 ** (-3)].StellarMass_median

    x_lin = np.linspace(np.min(x_full), np.max(x_full), 100)
    y_lin = [x_l ** (0.712) * (10 ** (-7.293)) for x_l in x_lin]
    y_lin_max = [x_l ** (0.712 + 0.04) * (10 ** (-7.293 + 0.06)) for x_l in x_lin]
    y_lin_min = [x_l ** (0.712 - 0.04) * (10 ** (-7.293 - 0.06)) for x_l in x_lin]

    plt.plot(x_lin, y_lin, c="magenta", ls="dashdot", zorder=z)
    plt.fill_between(x_lin, y_lin_min, y_lin_max, fc="magenta", alpha=0.1, ec=None, zorder=z)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"M$_{*}$ / M$_{\odot}$")
    plt.ylabel("SFR")

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "sfr.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    # plt.figure()


def plot_sfr_outliers(df_gal):
    df_sf = df_gal[df_gal.StarForm_ID_50 == 1]

    for _, row in df_sf.iterrows():
        if (row.int_stage == 1) & (row.StellarMass_median < 2 * 10**10):
            # directories
            FOLDER_DIR_DATA = "data/"  # data folder
            GALAXY_DIR = FOLDER_DIR_DATA + "galaxy_images/"  # galaxy images

            mag_id = row.MAGPIID
            file_name = "MAGPI" + str(mag_id) + "_resized.fits"

            hdu1 = fits.open(
                GALAXY_DIR + file_name
            )  # Opens the correct fits file that matches the idx from the dataframe

            im_g = hdu1[2].data  # G_mod band image
            im_i = hdu1[4].data  # i band image
            im_r = hdu1[6].data  # i band image

            im = make_lupton_rgb(im_i, im_r, im_g, Q=10, stretch=0.5)

            plt.figure(figsize=(8, 8))

            plt.imshow(im, origin="lower")
            plt.axis("off")

            #########################################
            # save figure
            PRINT_DIR = "plot_print/"
            FIG_NAME = "sfr_outlier_int_" + str(mag_id) + ".pdf"
            plt.savefig(PRINT_DIR + FIG_NAME)
            #########################################

    for _, row in df_sf.iterrows():
        if (row.SFR_50 < 1 * 10 ** (-1)) & (row.SFR_50 >= 10 ** (-3)):
            # directories
            FOLDER_DIR_DATA = "data/"  # data folder
            GALAXY_DIR = FOLDER_DIR_DATA + "galaxy_images/"  # galaxy images

            mag_id = row.MAGPIID
            file_name = "MAGPI" + str(mag_id) + "_resized.fits"

            hdu1 = fits.open(
                GALAXY_DIR + file_name
            )  # Opens the correct fits file that matches the idx from the dataframe

            im_g = hdu1[2].data  # G_mod band image
            im_i = hdu1[4].data  # i band image
            im_r = hdu1[6].data  # i band image

            im = make_lupton_rgb(im_i, im_r, im_g, Q=10, stretch=0.5)

            plt.figure(figsize=(8, 8))

            plt.imshow(im, origin="lower")
            plt.axis("off")

            #########################################
            # save figure
            PRINT_DIR = "plot_print/"
            FIG_NAME = "sfr_outlier_low_" + str(mag_id) + ".pdf"
            plt.savefig(PRINT_DIR + FIG_NAME)
            #########################################


#############################################


# log10(SFR) = 0.712±0.04 × log10(M∗) - 7.293±0.06

# plt.fill_between(x, y-error, y+error)

# d = np.log10(y) - np.log10( x ** (0.712) * ( 10 ** (-7.293) ))
