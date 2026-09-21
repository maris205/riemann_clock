import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
from matplotlib.ticker import FormatStrFormatter

DATA_DIR = Path("data")
EXP_DIR = DATA_DIR / "exp"
THEORY_DIR = DATA_DIR / "theory"

exp_03 = np.load(EXP_DIR / "tscan_beta03_exp.npz")
exp_05 = np.load(EXP_DIR / "tscan_beta05_exp.npz")
exp_beta = np.load(EXP_DIR / "beta_scan_exp.npz")

theo_03 = np.load(THEORY_DIR / "tscan_beta03_theory.npz")
theo_05 = np.load(THEORY_DIR / "tscan_beta05_theory.npz")
theo_beta = np.load(THEORY_DIR / "beta_scan_theory.npz")

M, N = exp_05["x_mean"].shape

LINE_WIDTH_THEO = 1.40
LINE_WIDTH_tick = 1.0
MARKER_SIZE_A = 3.8
MARKER_SIZE_ZERO = 3.2
CAP_SIZE = 2.6

plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath, bm} \boldmath"
mpl.rcParams["mathtext.fontset"] = "custom"
mpl.rcParams["mathtext.rm"] = "Arial"
mpl.rcParams["mathtext.it"] = "Computer Modern:italic"
mpl.rcParams["mathtext.bf"] = "Computer Modern:bold"

plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"

FONT_NAME = "Arial"
FONT_SIZE_AXES_LABEL = 10.4
FONT_SIZE_TICK_LABEL = 8.8
FONT_SIZE_TITLE = 10.4
FONT_SIZE_LEGEND = 9.8
LETTER_FONT_SIZE = 12.0
LINE_WIDTH_ERR = 0.95
LINE_WIDTH_EXP_POINT = 0.95

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "mathtext.fontset": "stix",
    "axes.linewidth": LINE_WIDTH_tick,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 3.2,
    "ytick.major.size": 3.2,
    "xtick.major.width": 0.85,
    "ytick.major.width": 0.85,
})


def hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i+2], 16) for i in (0, 2, 4)]) / 255.0


COLOR_THEO_BLUE = hex2rgb("#0072B2")
COLOR_EXP_ORANGE = hex2rgb("#D55E00")
COLOR_SHADE_CYAN = hex2rgb("#DFF3FF")
COLOR_LOOP_BLUE = hex2rgb("#005EA8")

col_t = {
    "A_th": COLOR_THEO_BLUE,
    "A_exp": COLOR_EXP_ORANGE,
    "A_shade": COLOR_SHADE_CYAN,
    "ErrBar_A": COLOR_EXP_ORANGE,
    "X_exp": COLOR_EXP_ORANGE,
    "XY_th": COLOR_LOOP_BLUE,
}

rows_zeros = exp_05["rows_zeros"]

# ---------- beta = 0.5 t-scan ----------
abs_mean_exp_05 = exp_05["abs_mean"]
std_abs_05 = exp_05["abs_std"]
abs_m_05 = abs_mean_exp_05.T.reshape(-1)
lAbs_05 = std_abs_05.T.reshape(-1)
uAbs_05 = lAbs_05.copy()

x_exp_05 = exp_05["t_flat"]
x_th_05 = theo_05["t_flat"]

X_t_05 = theo_05["curve_x_full"].T.reshape(-1)
Y_t_05 = theo_05["curve_y_full"].T.reshape(-1)
Abs_th_05 = np.sqrt(np.real(X_t_05)**2 + np.real(Y_t_05)**2)

X_mean_s_05 = exp_05["x_mean"]
Y_mean_s_05 = exp_05["y_mean"]
std_X_05 = exp_05["x_std"]
std_Y_05 = exp_05["y_std"]

X_theo_full_xy_05 = theo_05["xy_x_full"]
Y_theo_full_xy_05 = theo_05["xy_y_full"]

# ---------- beta = 0.3 t-scan ----------
abs_mean_exp_03 = exp_03["abs_mean"]
std_abs_03 = exp_03["abs_std"]
abs_m_03 = abs_mean_exp_03.T.reshape(-1)
lAbs_03 = std_abs_03.T.reshape(-1)
uAbs_03 = lAbs_03.copy()

x_exp_03 = exp_03["t_flat"]
x_th_03 = theo_03["t_flat"]

X_t_03 = theo_03["curve_x_full"].T.reshape(-1)
Y_t_03 = theo_03["curve_y_full"].T.reshape(-1)
Abs_th_03 = np.sqrt(np.real(X_t_03)**2 + np.real(Y_t_03)**2)

X_mean_s_03 = exp_03["x_mean"]
Y_mean_s_03 = exp_03["y_mean"]
std_X_03 = exp_03["x_std"]
std_Y_03 = exp_03["y_std"]

X_theo_full_xy_03 = theo_03["xy_x_full"]
Y_theo_full_xy_03 = theo_03["xy_y_full"]

# ---------- selected rows for d/e ----------
X_data_03_exp = X_mean_s_03[rows_zeros, :]
X_data_03_th = np.real(X_theo_full_xy_03[rows_zeros, 200:600+1])
Y_data_03_exp = Y_mean_s_03[rows_zeros, :]
Y_data_03_th = np.real(Y_theo_full_xy_03[rows_zeros, 200:600+1])

X_data_05_exp = X_mean_s_05[rows_zeros, :]
X_data_05_th = np.real(X_theo_full_xy_05[rows_zeros, 200:600+1])
Y_data_05_exp = Y_mean_s_05[rows_zeros, :]
Y_data_05_th = np.real(Y_theo_full_xy_05[rows_zeros, 200:600+1])

std_X_03_sub = std_X_03[rows_zeros, :]
std_Y_03_sub = std_Y_03[rows_zeros, :]
std_X_05_sub = std_X_05[rows_zeros, :]
std_Y_05_sub = std_Y_05[rows_zeros, :]

# ---------- beta scan ----------
X_theo_beta = theo_beta["curve_x"]
Y_theo_beta = theo_beta["curve_y"]
abs_theo_beta = np.sqrt(X_theo_beta**2 + Y_theo_beta**2)

BETA_VALUES = exp_beta["beta"]
X_mean_exp_beta = exp_beta["x_mean"]
Y_mean_exp_beta = exp_beta["y_mean"]
abs_mean_exp_beta = exp_beta["abs_mean"]
err_X_beta = exp_beta["x_std"]
err_Y_beta = exp_beta["y_std"]
err_abs_beta = exp_beta["abs_std"]

fig = plt.figure(figsize=(7, 3.2))
fig.patch.set_facecolor("white")

gs = fig.add_gridspec(
    2, 3,
    width_ratios=[1.5, 1.5, 0.82],
    height_ratios=[1.45, 1.02],
    left=0.060, right=0.995,
    bottom=0.108, top=0.988,
    wspace=0.18, hspace=0.25
)

ax_A = fig.add_subplot(gs[0, 0])
ax_B = fig.add_subplot(gs[0, 1], sharey=ax_A)
ax_C = fig.add_subplot(gs[0, 2])

gs_D = gs[1, 0].subgridspec(1, 3, wspace=0.16)
gs_E = gs[1, 1].subgridspec(1, 3, wspace=0.16)

ax_D_zeros = [fig.add_subplot(gs_D[0, k]) for k in range(3)]
ax_E_zeros = [fig.add_subplot(gs_E[0, k], sharey=ax_D_zeros[0]) for k in range(3)]
ax_F = fig.add_subplot(gs[1, 2])


def add_panel_label(fig, ax, label, dx=0.028, dy=0.020):
    pos = ax.get_position()
    fig.text(
        pos.x0 - dx, pos.y1 + dy, label,
        fontsize=LETTER_FONT_SIZE,
        fontweight="bold",
        ha="left",
        va="top"
    )


for ax in [ax_A, ax_B, ax_C, ax_F] + ax_D_zeros + ax_E_zeros:
    ax.set_facecolor("white")
    ax.tick_params(labelsize=FONT_SIZE_TICK_LABEL)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)

for ax in ax_D_zeros + ax_E_zeros + [ax_F]:
    ax.set_aspect("auto")

# ---------- panel a ----------
ax_A.fill_between(
    x_th_03, Abs_th_03, 0.0,
    facecolor=col_t["A_shade"], alpha=0.9, linewidth=0
)
h_th_A, = ax_A.plot(
    x_th_03, Abs_th_03,
    color=col_t["A_th"], linewidth=LINE_WIDTH_THEO
)
h_exp_A = ax_A.errorbar(
    x_exp_03, abs_m_03, yerr=[lAbs_03, uAbs_03],
    fmt="o",
    ecolor=col_t["ErrBar_A"],
    elinewidth=LINE_WIDTH_ERR,
    capsize=CAP_SIZE,
    capthick=LINE_WIDTH_ERR,
    markerfacecolor="none",
    markeredgecolor=col_t["A_exp"],
    markeredgewidth=LINE_WIDTH_EXP_POINT,
    markersize=MARKER_SIZE_A,
)

xmin_A = min(x_th_03.min(), x_th_05.min()) - 0.4
xmax_A = max(x_th_03.max(), x_th_05.max()) + 0.4

ax_A.set_xlim(xmin_A, xmax_A)
ax_A.set_ylim(0, 0.63)
ax_A.set_xticks([10, 15, 20, 25, 30, 35])
ax_A.set_yticks([0, 0.25, 0.5])
ax_A.set_ylabel(r"$|\mathcal{L}|$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.8)
ax_A.set_xlabel(r"$t$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.5)
ax_A.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_A.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_A.text(
    0.50, 0.92, r"$\beta=0.3$",
    transform=ax_A.transAxes,
    fontsize=FONT_SIZE_TITLE,
    fontweight="bold",
    ha="center",
    va="bottom"
)
ax_A.legend(
    [h_th_A, h_exp_A],
    [r"$\mathcal{L}_{\mathrm{theo}}$", r"$\mathcal{L}_{\mathrm{exp}}$"],
    loc="upper left",
    bbox_to_anchor=(0.02, 0.92),
    frameon=False,
    ncol=2,
    fontsize=FONT_SIZE_LEGEND,
    columnspacing=0.25,
    handlelength=1.6,
    handletextpad=0.5,
    borderaxespad=0.0
)

# ---------- panel b ----------
ax_B.fill_between(
    x_th_05, Abs_th_05, 0.0,
    facecolor=col_t["A_shade"], alpha=0.9, linewidth=0
)
ax_B.plot(
    x_th_05, Abs_th_05,
    color=col_t["A_th"], linewidth=LINE_WIDTH_THEO
)
ax_B.errorbar(
    x_exp_05, abs_m_05, yerr=[lAbs_05, uAbs_05],
    fmt="o",
    ecolor=col_t["ErrBar_A"],
    elinewidth=LINE_WIDTH_ERR,
    capsize=CAP_SIZE,
    capthick=LINE_WIDTH_ERR,
    markerfacecolor="none",
    markeredgecolor=col_t["A_exp"],
    markeredgewidth=LINE_WIDTH_EXP_POINT,
    markersize=MARKER_SIZE_A,
)
ax_B.set_xlim(xmin_A, xmax_A)
ax_B.set_ylim(0, 0.63)
ax_B.set_xticks([10, 15, 20, 25, 30, 35])
ax_B.set_yticks([0, 0.25, 0.5])
ax_B.tick_params(axis="y", labelleft=False)
ax_B.set_xlabel(r"$t$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.5)
ax_B.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_B.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_B.text(
    0.50, 0.92, r"$\beta=0.5$",
    transform=ax_B.transAxes,
    fontsize=FONT_SIZE_TITLE,
    fontweight="bold",
    ha="center",
    va="bottom"
)

# ---------- panel c ----------
idx_beta = slice(1, 18)
ax_C.fill_between(
    BETA_VALUES[idx_beta],
    abs_theo_beta[idx_beta],
    0.0,
    facecolor=col_t["A_shade"],
    alpha=0.9,
    linewidth=0,
)
ax_C.plot(
    BETA_VALUES[idx_beta],
    abs_theo_beta[idx_beta],
    "-",
    color=col_t["A_th"],
    linewidth=LINE_WIDTH_THEO,
)
ax_C.errorbar(
    BETA_VALUES[idx_beta],
    abs_mean_exp_beta[idx_beta],
    yerr=err_abs_beta[idx_beta],
    fmt="o",
    ecolor=col_t["ErrBar_A"],
    elinewidth=LINE_WIDTH_ERR,
    capsize=CAP_SIZE,
    capthick=LINE_WIDTH_ERR,
    markerfacecolor="none",
    markeredgecolor=col_t["X_exp"],
    markeredgewidth=LINE_WIDTH_EXP_POINT,
    markersize=MARKER_SIZE_A,
)
ax_C.axvline(
    0.5,
    ymin=0,
    ymax=0.82,
    color=(0.5, 0.5, 0.5),
    linestyle=":",
    linewidth=0.8,
)
ax_C.set_xlim(0, 1)
ax_C.set_ylim(0, 0.18)
ax_C.set_xticks([0.1, 0.5, 0.9])
ax_C.set_yticks([0, 0.08, 0.16])
ax_C.set_xlabel(r"$\beta$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.5)
ax_C.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_C.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_C.text(
    0.52, 0.92, r"$t=14.134725$",
    transform=ax_C.transAxes,
    fontsize=FONT_SIZE_TITLE,
    fontweight="bold",
    ha="center",
    va="bottom"
)

# ---------- panels d / e ----------
t_labels_pos = np.array([0.43, 0.71, 0.41])
t_labels_text = [
    r"$t \in [12,16]$",
    r"$t \in [23,27]$",
    r"$t \in [31,35]$"
]
F_label_y = np.array([0.05, 0.83, 0.11])

xlims_d = [
    (-0.12, 0.46),
    (-0.22, 0.54),
    (-0.12, 0.52),
]
xlims_e = [
    (-0.10, 0.46),
    (-0.22, 0.54),
    (-0.10, 0.55),
]

for k, axd in enumerate(ax_D_zeros):
    axd.plot(
        X_data_03_th[k, :], Y_data_03_th[k, :],
        "-", linewidth=LINE_WIDTH_THEO, color=col_t["XY_th"]
    )
    axd.errorbar(
        X_data_03_exp[k, :], Y_data_03_exp[k, :],
        yerr=std_Y_03_sub[k, :], xerr=std_X_03_sub[k, :],
        fmt="o",
        ecolor=col_t["ErrBar_A"],
        elinewidth=LINE_WIDTH_ERR,
        capsize=CAP_SIZE,
        capthick=LINE_WIDTH_ERR,
        markerfacecolor="none",
        markeredgecolor=col_t["X_exp"],
        markeredgewidth=LINE_WIDTH_ERR,
        markersize=MARKER_SIZE_ZERO,
    )
    axd.axvline(0, linestyle=":", color="k", linewidth=0.55)
    axd.axhline(0, linestyle=":", color="k", linewidth=0.55)
    axd.set_xlim(*xlims_d[k])
    axd.set_ylim(-0.42, 0.36)
    axd.set_xticks([0.0, 0.2, 0.4])
    axd.set_yticks([-0.25, 0, 0.25])
    axd.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    axd.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    if k == 0:
        axd.set_ylabel(r"$\langle\sigma_y\rangle$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.3)
    else:
        axd.tick_params(axis="y", labelleft=False)
    if k == 1:
        axd.set_xlabel(r"$\langle\sigma_x\rangle$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.4)
    axd.text(
        t_labels_pos[k], F_label_y[k], t_labels_text[k],
        transform=axd.transAxes,
        fontsize=FONT_SIZE_TICK_LABEL,
        ha="center",
        va="bottom",
        color="black"
    )

for k, axe in enumerate(ax_E_zeros):
    axe.plot(
        X_data_05_th[k, :], Y_data_05_th[k, :],
        "-", linewidth=LINE_WIDTH_THEO, color=col_t["XY_th"]
    )
    axe.errorbar(
        X_data_05_exp[k, :], Y_data_05_exp[k, :],
        yerr=std_Y_05_sub[k, :], xerr=std_X_05_sub[k, :],
        fmt="o",
        ecolor=col_t["ErrBar_A"],
        elinewidth=LINE_WIDTH_ERR,
        capsize=CAP_SIZE,
        capthick=LINE_WIDTH_ERR,
        markerfacecolor="none",
        markeredgecolor=col_t["X_exp"],
        markeredgewidth=LINE_WIDTH_ERR,
        markersize=MARKER_SIZE_ZERO,
    )
    axe.axvline(0, linestyle=":", color="k", linewidth=0.55)
    axe.axhline(0, linestyle=":", color="k", linewidth=0.55)
    axe.set_xlim(*xlims_e[k])
    axe.set_ylim(-0.42, 0.32)
    axe.set_xticks([0.0, 0.2, 0.4])
    axe.set_yticks([-0.25, 0, 0.25])
    axe.tick_params(axis="y", labelleft=False)
    axe.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    axe.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    if k == 1:
        axe.set_xlabel(r"$\langle\sigma_x\rangle$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.4)
    axe.text(
        t_labels_pos[k], F_label_y[k], t_labels_text[k],
        transform=axe.transAxes,
        fontsize=FONT_SIZE_TICK_LABEL,
        ha="center",
        va="bottom",
        color="black"
    )

# ---------- panel f ----------
idx_beta = slice(6, 13)
ax_F.plot(
    X_theo_beta[idx_beta], Y_theo_beta[idx_beta],
    "-", linewidth=LINE_WIDTH_THEO, color=col_t["XY_th"]
)
ax_F.errorbar(
    X_mean_exp_beta[idx_beta], Y_mean_exp_beta[idx_beta],
    xerr=err_X_beta[idx_beta], yerr=err_Y_beta[idx_beta],
    fmt="o",
    ecolor=col_t["ErrBar_A"],
    elinewidth=LINE_WIDTH_ERR,
    capsize=CAP_SIZE,
    capthick=LINE_WIDTH_ERR,
    markerfacecolor="none",
    markeredgecolor=col_t["X_exp"],
    markeredgewidth=LINE_WIDTH_ERR,
    markersize=MARKER_SIZE_ZERO,
)
ax_F.axvline(0, linestyle=":", color="k", linewidth=0.55)
ax_F.axhline(0, linestyle=":", color="k", linewidth=0.55)

x_min = min(X_theo_beta[idx_beta].min(), X_mean_exp_beta[idx_beta].min())
x_max = max(X_theo_beta[idx_beta].max(), X_mean_exp_beta[idx_beta].max())

ax_F.set_xlim(x_min - 0.012, x_max + 0.012)
ax_F.set_ylim(-0.01, 0.045)
ax_F.set_xticks([0, 0.05])
ax_F.set_yticks([0, 0.03])
ax_F.set_xlabel(r"$\langle\sigma_x\rangle$", fontsize=FONT_SIZE_AXES_LABEL, labelpad=1.4)
ax_F.xaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_F.yaxis.set_major_formatter(FormatStrFormatter("%g"))
ax_F.text(
    0.72, 0.86, r"$\beta \in [0.35,0.65]$",
    transform=ax_F.transAxes,
    fontsize=FONT_SIZE_TICK_LABEL,
    ha="center",
    va="bottom",
    color="black"
)

# ---------- aligned panel letters ----------
add_panel_label(fig, ax_A, "a")
add_panel_label(fig, ax_B, "b")
add_panel_label(fig, ax_C, "c")
add_panel_label(fig, ax_D_zeros[0], "d", dy=0.04)
add_panel_label(fig, ax_E_zeros[0], "e", dy=0.04)
add_panel_label(fig, ax_F, "f", dy=0.04)

plt.savefig("Fig4.pdf", bbox_inches="tight", pad_inches=0.003)