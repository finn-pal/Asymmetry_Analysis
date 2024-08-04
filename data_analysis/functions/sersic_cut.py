import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import curve_fit


def fit_func(x, m, b):
    return m * x + b


def sersic_divide(df_gal):
    fs = 27
    z = 0

    # if bad sersic fit then set to nan, these will later be defined by this relationship
    df_gal.loc[(df_gal.n_galfit == 8) | (df_gal.n_galfit == 0.5), "n_galfit"] = np.nan
    df_gal.loc[(df_gal.n_err_galfit > 0.3), "n_galfit"] = np.nan

    x = np.array(df_gal[~(df_gal.n_galfit).isnull()].n_galfit)
    x_err = np.array(df_gal[~(df_gal.n_galfit).isnull()].n_err_galfit)

    y = np.array(df_gal[~(df_gal.n_galfit).isnull()].C)
    y_err = np.array(df_gal[~(df_gal.n_galfit).isnull()].C_err)

    x_log = np.log10(x)
    y_log = np.log10(y)

    x_log_p = np.delete(x_log, np.where(x_log < 0))
    y_log_p = np.delete(y_log, np.where(x_log < 0))

    # m, b = np.polyfit(x_log_p, y_log_p, 1)

    popt, pcov = curve_fit(fit_func, x_log_p, y_log_p)

    m, b = popt
    err_cov_log = np.sqrt(np.diag(pcov))
    b_err_log = err_cov_log[1]
    b_err = 10 ** (b_err_log)

    xfit = np.linspace(0.5, 9, 200)
    xfitlog = np.log10(xfit)
    yfitlog = np.array([fit_func(xi, m, b) for xi in xfitlog])

    # xfit = np.arange(np.min(x), np.max(x))
    # yfitlog = m * np.log10(xfit) + b
    yfit = 10 ** (yfitlog)

    _, ax = plt.subplots(figsize=(14, 8))
    plt.rcParams.update({"font.size": fs - 5})

    zord_err = 0
    ax.errorbar(x, y, xerr=x_err, yerr=y_err, color="grey", fmt="none", elinewidth=1.5, zorder=zord_err)

    ax.scatter(x, y, s=40, c="blue", label="MAGPI Galaxies")
    ax.plot(xfit, yfit, "k-.", lw=2.5, label="Fit to MAGPI Galaxies")

    ax.fill_between(xfit, yfit - b_err, yfit + b_err, fc="grey", alpha=0.1, ec=None, zorder=z)

    ax.set(xscale="log", yscale="log", xlim=(0.6, 9))
    ax.set_xlabel("S$\\acute{e}$rsic Index (n)")
    ax.set_ylabel("Concentration (C)")

    ax.xaxis.label.set_size(fs)
    ax.yaxis.label.set_size(fs)
    ax.tick_params(axis="both", which="both", labelsize=fs - 3, length=5, width=1.5)

    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(mticker.ScalarFormatter())

    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    # ax.legend(loc="upper left")

    #########################################
    # save figure
    PRINT_DIR = "plot_print/"
    FIG_NAME = "sersic_con.pdf"
    plt.savefig(PRINT_DIR + FIG_NAME)
    #########################################

    ############## Calculate new sersic for poor fits

    C_val = df_gal.C
    n_val = 10 ** ((np.log10(C_val) - b) / m)

    n_val_lst = [fit if np.isnan(galfit) else np.nan for fit, galfit in zip(n_val, df_gal.n_galfit)]

    df_gal["n_fit"] = n_val_lst

    n_use = []
    for i in range(0, len(df_gal)):
        if np.isnan(df_gal["n_galfit"][i]):
            n_use.append(df_gal["n_fit"][i])
        else:
            n_use.append(df_gal["n_galfit"][i])

    df_gal["n_use"] = n_use

    ############### Print relationship

    print("\n########################")
    print("Sersic - Concentration")

    print(
        r"\log_{10}(C) = "
        + "("
        + str(np.round(m, 2))
        + r" \pm "
        + str(np.round(err_cov_log[0], 2))
        + ")"
        + r"\log_{10}(n)"
        + " + ("
        + str(np.round(b, 2))
        + r" \pm "
        + str(np.round(err_cov_log[1], 2))
        + ")"
    )

    return df_gal
