# -*- coding: utf-8 -*-
"""
3D simulation of hydraulic classification (elutriation) in an upflow column.

Merges the three previous studies into a single model:

  * stokes.py  -> Stokes law (laminar regime, Re < 1)
  * newton.py  -> Newton law (turbulent regime, Cd ~ 0.44)
  * b.py       -> intermediate regime with Cd(Re) and a shape factor

The integrated dynamics (drag, shape factor and boundary conditions) follow
b.py; the Stokes and Newton laws appear only as reference curves in the
comparison panel.

Output: a GIF with the 3D column, the VT(d) panel and the live particle
balance.
"""

import math
import os
import random

import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
from PIL import Image

# =====================================================================
# 1. PHYSICAL CONSTANTS
# =====================================================================

G = 9.81            # gravity                          [m/s^2]
RHO_F = 1000.0      # water density                    [kg/m^3]
MU = 1.0e-3         # water dynamic viscosity          [Pa.s]

RHO_QUARTZ = 2650.0     # kg/m^3
RHO_HEMATITE = 5200.0   # kg/m^3

CD_NEWTON = 0.44        # Cd of the Newton regime (used only in the curve)

# =====================================================================
# 2. COLUMN GEOMETRY AND THE 6 LPM FEED
# =====================================================================

SIDE = 0.10                      # square column section              [m]
HEIGHT = 0.15                    # useful column height               [m]
AREA = SIDE * SIDE               # cross-sectional area               [m^2]

Q_LPM = 6.0                      # imposed flow rate                  [L/min]
Q = Q_LPM * 1.0e-3 / 60.0        # volumetric flow rate               [m^3/s]
V_FLUID = Q / AREA               # superficial upflow velocity        [m/s]

# =====================================================================
# 3. NUMERICAL PARAMETERS (short run: 4 s of process)
# =====================================================================

DT = 0.002           # integration time step      [s]
T_TOTAL = 4.0        # simulated time             [s]
N_STEPS = int(round(T_TOTAL / DT))
STEPS_PER_FRAME = 16
FPS = 25

N_PARTICLES = 240
D_MIN, D_MAX = 20e-6, 250e-6     # particle diameter range            [m]

SEED = 7

# =====================================================================
# 4. PALETTE (validated against the dark surface)
# =====================================================================

BACKGROUND = "#131313"
PANEL = "#1a1a19"
INK_1 = "#ffffff"
INK_2 = "#c3c2b7"
INK_3 = "#8a8a80"

C_QUARTZ = "#c98500"
C_HEMATITE = "#d55181"

C_STOKES = "#3987e5"
C_NEWTON = "#d95926"
C_INTERMEDIATE = "#199e70"

C_FLOW = "#3987e5"


# =====================================================================
# 5. TERMINAL VELOCITY LAWS
# =====================================================================

def reynolds(v, d):
    """Particle Reynolds number (b.py)."""
    return RHO_F * abs(v) * d / MU


def drag_cd(Re):
    """Cd(Re) correlation, valid from the laminar to the turbulent regime (b.py)."""
    Re = max(Re, 1e-9)
    return (0.63 + 4.8 / math.sqrt(Re)) ** 2


def vt_intermediate(rho_p, d, shape_factor=1.0, n_iter=200, tol=1e-12):
    """
    Terminal velocity in the intermediate regime, by fixed-point iteration:
    VT -> Re -> Cd(Re) -> VT ...

    shape_factor (psi): 1.0 = nearly spherical, 0.4 = flat/irregular.
    A less spherical particle has a larger effective drag -> Cd_eff = Cd / psi.
    """
    if rho_p <= RHO_F:
        return 0.0

    v = 0.01
    for _ in range(n_iter):
        cd_eff = drag_cd(reynolds(v, d)) / shape_factor
        v_new = math.sqrt(4.0 * (rho_p - RHO_F) * G * d / (3.0 * RHO_F * cd_eff))
        if abs(v_new - v) < tol:
            v = v_new
            break
        v = v_new
    return v


def vt_stokes(rho_p, d):
    """Stokes law (stokes.py), written in terms of diameter."""
    return (rho_p - RHO_F) * G * d * d / (18.0 * MU)


def vt_newton(rho_p, d):
    """Newton law (newton.py), with Cd = 0.44 and radius r = d/2."""
    r = 0.5 * d
    return math.sqrt((8.0 / (3.0 * CD_NEWTON)) * G * ((rho_p - RHO_F) / RHO_F) * r)


def cut_diameter(rho_p, shape_factor=1.0, v_target=V_FLUID):
    """Diameter at which VT_intermediate equals the fluid velocity (theoretical d50)."""
    lo, hi = 1e-7, 5e-3
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if vt_intermediate(rho_p, mid, shape_factor) < v_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# =====================================================================
# 6. PARTICLE POPULATION
# =====================================================================

def make_particles(n=N_PARTICLES, seed=SEED):
    rng = random.Random(seed)
    p = {}

    d = np.array([rng.uniform(D_MIN, D_MAX) for _ in range(n)])
    rho = np.array([rng.choice([RHO_QUARTZ, RHO_HEMATITE]) for _ in range(n)])
    psi = np.array([rng.uniform(0.4, 1.0) for _ in range(n)])

    p["d"] = d
    p["rho"] = rho
    p["psi"] = psi
    p["vt"] = np.array([vt_intermediate(rho[i], d[i], psi[i]) for i in range(n)])
    p["re"] = np.array([reynolds(p["vt"][i], d[i]) for i in range(n)])

    p["x"] = np.array([rng.uniform(-0.45, 0.45) * SIDE for _ in range(n)])
    p["z"] = np.array([rng.uniform(-0.45, 0.45) * SIDE for _ in range(n)])
    p["y"] = np.array([rng.uniform(0.05, 0.95) * HEIGHT for _ in range(n)])

    p["hematite"] = rho > 3000.0
    p["flat"] = psi < 0.65
    p["rng"] = rng
    return p


def step(p):
    """
    One integration step.

    Vertical velocity: vy = v_fluid - VT  (fluid rises, particle settles).
    Lateral dispersion: proportional to (1 - psi) -- an irregular particle
    wanders more, as in b.py.

    Boundary conditions (b.py):
      * side walls in x and z : periodic (the particle re-enters on the
                                opposite face);
      * bottom (y < 0)        : the particle has settled -> reinjected at the top;
      * top    (y > HEIGHT)   : the particle was carried over -> reinjected at
                                the bottom.
    """
    rng = p["rng"]
    n = len(p["d"])

    vy = V_FLUID - p["vt"]
    spread = (1.0 - p["psi"]) * 0.4 * V_FLUID
    vx = np.array([rng.uniform(-1.0, 1.0) for _ in range(n)]) * spread
    vz = np.array([rng.uniform(-1.0, 1.0) for _ in range(n)]) * spread

    p["x"] += vx * DT
    p["y"] += vy * DT
    p["z"] += vz * DT

    half = 0.5 * SIDE
    p["x"] = np.where(p["x"] > half, p["x"] - SIDE, p["x"])
    p["x"] = np.where(p["x"] < -half, p["x"] + SIDE, p["x"])
    p["z"] = np.where(p["z"] > half, p["z"] - SIDE, p["z"])
    p["z"] = np.where(p["z"] < -half, p["z"] + SIDE, p["z"])

    settled = p["y"] < 0.0
    carried = p["y"] > HEIGHT

    for idx in np.flatnonzero(settled):
        p["y"][idx] = HEIGHT
        p["x"][idx] = rng.uniform(-0.45, 0.45) * SIDE
        p["z"][idx] = rng.uniform(-0.45, 0.45) * SIDE

    for idx in np.flatnonzero(carried):
        p["y"][idx] = 0.0
        p["x"][idx] = rng.uniform(-0.45, 0.45) * SIDE
        p["z"][idx] = rng.uniform(-0.45, 0.45) * SIDE

    return int(settled.sum()), int(carried.sum())


# =====================================================================
# 7. FIGURE
# =====================================================================

def _style_3d_axes(ax):
    ax.set_facecolor(BACKGROUND)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.07, 0.07, 0.07, 1.0))
        axis._axinfo["grid"]["color"] = (1, 1, 1, 0.06)
        axis.line.set_color(INK_3)
    ax.tick_params(colors=INK_2, labelsize=7.5)


def _column_box(ax):
    """Column edges, in mm."""
    a = 0.5 * SIDE * 1000.0
    h = HEIGHT * 1000.0
    base = [(-a, -a), (a, -a), (a, a), (-a, a), (-a, -a)]

    for y in (0.0, h):
        xs = [pt[0] for pt in base]
        zs = [pt[1] for pt in base]
        ax.plot(xs, zs, [y] * 5, color=C_FLOW, lw=1.0, alpha=0.45)

    for (x, z) in base[:4]:
        ax.plot([x, x], [z, z], [0.0, h], color=C_FLOW, lw=0.8, alpha=0.25)


def _flow_arrows(ax):
    """Upflow feed of 6 LPM at the bottom of the column."""
    a = 0.5 * SIDE * 1000.0
    grid = np.linspace(-a * 0.7, a * 0.7, 3)
    for x in grid:
        for z in grid:
            ax.quiver(x, z, 2.0, 0, 0, 22.0, color=C_FLOW,
                      alpha=0.55, linewidth=0.9, arrow_length_ratio=0.35)


def _vt_panel(ax, d50_q, d50_h):
    """Comparison panel: VT(d) for the three laws (quartz) + Cd(Re) for hematite."""
    ax.set_facecolor(PANEL)

    d = np.linspace(5e-6, 2e-3, 400)
    y_st = np.array([vt_stokes(RHO_QUARTZ, di) for di in d]) * 1000.0
    y_nw = np.array([vt_newton(RHO_QUARTZ, di) for di in d]) * 1000.0
    y_in = np.array([vt_intermediate(RHO_QUARTZ, di) for di in d]) * 1000.0
    y_in_h = np.array([vt_intermediate(RHO_HEMATITE, di) for di in d]) * 1000.0
    d_um = d * 1e6

    ax.loglog(d_um, y_st, color=C_STOKES, lw=2.0, label="Stokes — quartz")
    ax.loglog(d_um, y_nw, color=C_NEWTON, lw=2.0, label="Newton — quartz")
    ax.loglog(d_um, y_in, color=C_INTERMEDIATE, lw=2.0,
              label="Cd(Re) — quartz 2650")
    ax.loglog(d_um, y_in_h, color=C_INTERMEDIATE, lw=2.0, ls=(0, (5, 2)),
              label="Cd(Re) — hematite 5200")

    u = V_FLUID * 1000.0
    ax.axhline(u, color=INK_2, lw=1.2, ls=(0, (4, 3)), alpha=0.9)
    ax.text(d_um[-1], u * 1.15, f"u = {u:.0f} mm/s  ({Q_LPM:.0f} LPM)",
            color=INK_2, fontsize=7.5, ha="right", va="bottom")

    ax.axvspan(D_MIN * 1e6, D_MAX * 1e6, color=INK_2, alpha=0.07)
    ax.text(math.sqrt(D_MIN * D_MAX) * 1e6, 4.5e2,
            "simulated range", color=INK_3, fontsize=7, ha="center", va="top")

    for d50 in (d50_q, d50_h):
        ax.plot([d50 * 1e6], [u], marker="o", ms=6, color=C_INTERMEDIATE,
                mec=PANEL, mew=1.5, zorder=5)
    ax.text(d50_q * 1e6 * 1.2, u * 0.45, f"d50 quartz\n{d50_q*1e6:.0f} µm",
            color=INK_2, fontsize=7, ha="left", va="top")
    ax.text(d50_h * 1e6 * 0.8, u * 0.45, f"d50 hematite\n{d50_h*1e6:.0f} µm",
            color=INK_2, fontsize=7, ha="right", va="top")

    ax.set_xlabel("particle diameter  [µm]", color=INK_2, fontsize=8)
    ax.set_ylabel("terminal velocity  [mm/s]", color=INK_2, fontsize=8)
    ax.set_title("Terminal velocity (spherical particle, ψ = 1)",
                 color=INK_1, fontsize=9, loc="left", pad=6)
    ax.tick_params(colors=INK_2, labelsize=7.5, which="both")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_3)
        ax.spines[side].set_linewidth(0.8)
    ax.grid(True, which="major", color=INK_1, alpha=0.08, lw=0.7)
    ax.set_ylim(1e-2, 1e3)

    leg = ax.legend(loc="upper left", frameon=False, fontsize=7.5,
                    labelcolor=INK_2, handlelength=1.6, borderaxespad=0.2)
    for text in leg.get_texts():
        text.set_color(INK_2)


def build_figure(p):
    fig = plt.figure(figsize=(11.2, 6.3), dpi=86, facecolor=BACKGROUND)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0],
                          height_ratios=[1.0, 0.58],
                          left=0.02, right=0.965, top=0.82, bottom=0.06,
                          wspace=0.18, hspace=0.55)

    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_vt = fig.add_subplot(gs[0, 1])
    ax_txt = fig.add_subplot(gs[1, 1])
    ax_txt.axis("off")

    _style_3d_axes(ax3d)
    a = 0.5 * SIDE * 1000.0
    ax3d.set_xlim(-a, a)
    ax3d.set_ylim(-a, a)
    ax3d.set_zlim(0, HEIGHT * 1000.0)
    ax3d.set_box_aspect((1, 1, 1.7))
    ax3d.set_xlabel("x [mm]", color=INK_2, fontsize=8, labelpad=-2)
    ax3d.set_ylabel("z [mm]", color=INK_2, fontsize=8, labelpad=-2)
    ax3d.set_zlabel("height y [mm]", color=INK_2, fontsize=8, labelpad=-2)

    _column_box(ax3d)
    _flow_arrows(ax3d)

    size = 6.0 + 46.0 * (p["d"] - D_MIN) / (D_MAX - D_MIN)
    groups = {}
    for name, mask, colour, marker in (
        ("quartz, spherical", (~p["hematite"]) & (~p["flat"]), C_QUARTZ, "o"),
        ("quartz, flat", (~p["hematite"]) & p["flat"], C_QUARTZ, "D"),
        ("hematite, spherical", p["hematite"] & (~p["flat"]), C_HEMATITE, "o"),
        ("hematite, flat", p["hematite"] & p["flat"], C_HEMATITE, "D"),
    ):
        sc = ax3d.scatter(p["x"][mask] * 1000.0, p["z"][mask] * 1000.0,
                          p["y"][mask] * 1000.0,
                          s=size[mask], c=colour, marker=marker,
                          edgecolors=BACKGROUND, linewidths=0.4,
                          depthshade=True, label=name)
        groups[name] = (sc, mask)

    leg = ax3d.legend(loc="upper left", bbox_to_anchor=(-0.02, 0.98),
                      frameon=False, fontsize=7.5, scatterpoints=1,
                      handletextpad=0.4, labelspacing=0.35)
    for text in leg.get_texts():
        text.set_color(INK_2)

    d50_q = cut_diameter(RHO_QUARTZ)
    d50_h = cut_diameter(RHO_HEMATITE)
    _vt_panel(ax_vt, d50_q, d50_h)

    fig.text(0.02, 0.955,
             f"Elutriation in an upflow column — {Q_LPM:.0f} LPM · section "
             f"{SIDE*1000:.0f} × {SIDE*1000:.0f} mm · H = {HEIGHT*1000:.0f} mm · "
             f"water at 20 °C",
             color=INK_1, fontsize=12.5, ha="left", va="center")
    fig.text(0.02, 0.918,
             "Drag Cd(Re) = (0.63 + 4.8/√Re)² · shape correction Cd_eff = Cd/ψ",
             color=INK_3, fontsize=8.5, ha="left", va="center")
    fig.text(0.02, 0.888,
             "Boundaries: periodic side walls · reinjection at bottom and top",
             color=INK_3, fontsize=8.5, ha="left", va="center")

    return fig, ax3d, ax_txt, groups, size, (d50_q, d50_h)


# =====================================================================
# 8. ANIMATION AND GIF EXPORT
# =====================================================================

def _frame_to_image(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return Image.fromarray(buf[:, :, :3].copy())


def run(gif_name="elutriation_3d.gif", gif_colours=96):
    p = make_particles()
    fig, ax3d, ax_txt, groups, size, (d50_q, d50_h) = build_figure(p)

    fines_q = float(np.mean(p["vt"][~p["hematite"]] < V_FLUID)) * 100.0
    fines_h = float(np.mean(p["vt"][p["hematite"]] < V_FLUID)) * 100.0

    box = ax_txt.text(
        0.0, 1.06, "", color=INK_2, fontsize=9, family="monospace",
        ha="left", va="top", transform=ax_txt.transAxes, linespacing=1.6)

    frames = []
    n_settled = n_carried = 0

    for k in range(N_STEPS + 1):
        if k > 0:
            s, c = step(p)
            n_settled += s
            n_carried += c

        if k % STEPS_PER_FRAME:
            continue

        t = k * DT
        for name, (sc, mask) in groups.items():
            sc._offsets3d = (p["x"][mask] * 1000.0,
                             p["z"][mask] * 1000.0,
                             p["y"][mask] * 1000.0)

        ax3d.view_init(elev=16.0, azim=-62.0 + 0.28 * (k / STEPS_PER_FRAME))

        box.set_text(
            f"t = {t:5.2f} s   of {T_TOTAL:.1f} s   (dt = {DT*1000:.0f} ms)\n"
            f"u  = {V_FLUID*1000:5.2f} mm/s   Q = {Q_LPM:.1f} L/min\n"
            f"\n"
            f"carried over the top ..... {n_carried:4d}\n"
            f"settled at the bottom .... {n_settled:4d}\n"
            f"\n"
            f"VT < u  quartz   {fines_q:5.1f} %\n"
            f"VT < u  hematite {fines_h:5.1f} %\n"
            f"population Re: {p['re'].min():.2g} – {p['re'].max():.2g}"
        )

        frames.append(_frame_to_image(fig))

    plt.close(fig)

    base = frames[0].quantize(colors=gif_colours, method=Image.MEDIANCUT)
    paletted = [f.quantize(palette=base, dither=Image.FLOYDSTEINBERG)
                for f in frames]

    target = os.path.join(os.path.dirname(os.path.abspath(__file__)), gif_name)
    paletted[0].save(
        target, save_all=True, append_images=paletted[1:],
        duration=int(1000 / FPS), loop=0, optimize=True)

    return target, len(frames), n_settled, n_carried, d50_q, d50_h


def console_summary(d50_q, d50_h):
    print(f"Column .............. {SIDE*1000:.0f} x {SIDE*1000:.0f} x "
          f"{HEIGHT*1000:.0f} mm")
    print(f"Flow rate ........... {Q_LPM:.1f} L/min = {Q:.2e} m3/s")
    print(f"Upflow velocity ..... {V_FLUID*1000:.2f} mm/s")
    print(f"d50 quartz (psi=1) .. {d50_q*1e6:.1f} um")
    print(f"d50 hematite (psi=1)  {d50_h*1e6:.1f} um")
    print()
    print(" d [um]   rho      VT Stokes   VT Newton   VT Cd(Re)   Re")
    for d in (20e-6, 50e-6, 100e-6, 150e-6, 250e-6):
        for rho in (RHO_QUARTZ, RHO_HEMATITE):
            vi = vt_intermediate(rho, d)
            print(f"{d*1e6:7.0f}{rho:8.0f}{vt_stokes(rho,d)*1000:12.3f}"
                  f"{vt_newton(rho,d)*1000:12.1f}{vi*1000:12.3f}"
                  f"{reynolds(vi,d):9.3f}")


if __name__ == "__main__":
    target, n_frames, n_settled, n_carried, d50_q, d50_h = run()
    console_summary(d50_q, d50_h)
    print()
    print(f"GIF: {target}  ({n_frames} frames, {FPS} fps)")
    print(f"Balance over {T_TOTAL:.1f} s: {n_carried} carried over / "
          f"{n_settled} settled")
