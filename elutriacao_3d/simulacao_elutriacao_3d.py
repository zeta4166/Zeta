# -*- coding: utf-8 -*-
"""
Simulacao 3D de elutriacao / classificacao hidraulica em coluna ascendente.

Reune os tres estudos anteriores em um unico modelo:

  * stokes.py  -> lei de Stokes (regime laminar, Re < 1)
  * newton.py  -> lei de Newton (regime turbulento, Cd ~ 0.44)
  * b.py       -> regime intermediario com Cd(Re) e fator de forma

O modelo dinamico (arrasto, fator de forma e condicoes de contorno) segue o
b.py; as leis de Stokes e Newton entram apenas como curvas de referencia no
painel comparativo.

Saida: um GIF com a coluna em 3D, o painel VT(d) e os contadores de balanco.
"""

import math
import os
import random

import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registra a projecao 3d)
from PIL import Image

# =====================================================================
# 1. CONSTANTES FISICAS
# =====================================================================

G = 9.81            # gravidade                         [m/s^2]
RHO_F = 1000.0      # massa especifica da agua          [kg/m^3]
MU = 1.0e-3         # viscosidade dinamica da agua      [Pa.s]

RHO_QUARTZO = 2650.0    # kg/m^3
RHO_HEMATITA = 5200.0   # kg/m^3

CD_NEWTON = 0.44        # Cd do regime de Newton (usado so na curva teorica)

# =====================================================================
# 2. GEOMETRIA DA COLUNA E VAZAO DE 6 LPM
# =====================================================================

LADO = 0.10                      # secao quadrada da coluna          [m]
ALTURA = 0.15                    # altura util da coluna             [m]
AREA = LADO * LADO               # area da secao transversal         [m^2]

Q_LPM = 6.0                      # vazao imposta                     [L/min]
Q = Q_LPM * 1.0e-3 / 60.0        # vazao volumetrica                 [m^3/s]
V_FLUIDO = Q / AREA              # velocidade superficial ascendente [m/s]

# =====================================================================
# 3. PARAMETROS NUMERICOS (tempo curto: 4 s de processo)
# =====================================================================

DT = 0.002           # passo de integracao        [s]
T_TOTAL = 4.0        # tempo simulado             [s]
N_PASSOS = int(round(T_TOTAL / DT))
PASSOS_POR_QUADRO = 16
FPS = 25

N_PARTICULAS = 240
D_MIN, D_MAX = 20e-6, 250e-6     # faixa de diametros                [m]

SEED = 7

# =====================================================================
# 4. PALETA (validada para fundo escuro)
# =====================================================================

COR_FUNDO = "#131313"
COR_PAINEL = "#1a1a19"
INK_1 = "#ffffff"
INK_2 = "#c3c2b7"
INK_3 = "#8a8a80"

COR_QUARTZO = "#c98500"
COR_HEMATITA = "#d55181"

COR_STOKES = "#3987e5"
COR_NEWTON = "#d95926"
COR_INTERM = "#199e70"

COR_FLUXO = "#3987e5"


# =====================================================================
# 5. LEIS DE VELOCIDADE TERMINAL
# =====================================================================

def reynolds(v, d):
    """Reynolds da particula (b.py)."""
    return RHO_F * abs(v) * d / MU


def drag_cd(Re):
    """Correlacao Cd(Re) valida do regime laminar ao turbulento (b.py)."""
    Re = max(Re, 1e-9)
    return (0.63 + 4.8 / math.sqrt(Re)) ** 2


def vt_intermediario(rho_p, d, fator_forma=1.0, n_iter=200, tol=1e-12):
    """
    VT no regime intermediario por ponto fixo: Re -> Cd(Re) -> VT -> Re ...

    fator_forma (psi): 1.0 = quase esferica, 0.4 = achatada/irregular.
    Particula menos esferica tem arrasto efetivo maior -> Cd_eff = Cd / psi.
    """
    if rho_p <= RHO_F:
        return 0.0

    v = 0.01
    for _ in range(n_iter):
        Cd_eff = drag_cd(reynolds(v, d)) / fator_forma
        v_novo = math.sqrt(4.0 * (rho_p - RHO_F) * G * d / (3.0 * RHO_F * Cd_eff))
        if abs(v_novo - v) < tol:
            v = v_novo
            break
        v = v_novo
    return v


def vt_stokes(rho_p, d):
    """Lei de Stokes (stokes.py), escrita em diametro."""
    return (rho_p - RHO_F) * G * d * d / (18.0 * MU)


def vt_newton(rho_p, d):
    """Lei de Newton (newton.py), com Cd = 0.44 e raio r = d/2."""
    r = 0.5 * d
    return math.sqrt((8.0 / (3.0 * CD_NEWTON)) * G * ((rho_p - RHO_F) / RHO_F) * r)


def diametro_de_corte(rho_p, fator_forma=1.0, v_alvo=V_FLUIDO):
    """Diametro em que VT_intermediario = velocidade do fluido (d50 teorico)."""
    lo, hi = 1e-7, 5e-3
    for _ in range(80):
        meio = 0.5 * (lo + hi)
        if vt_intermediario(rho_p, meio, fator_forma) < v_alvo:
            lo = meio
        else:
            hi = meio
    return 0.5 * (lo + hi)


# =====================================================================
# 6. POPULACAO DE PARTICULAS
# =====================================================================

def criar_particulas(n=N_PARTICULAS, seed=SEED):
    rng = random.Random(seed)
    p = {}

    d = np.array([rng.uniform(D_MIN, D_MAX) for _ in range(n)])
    rho = np.array([rng.choice([RHO_QUARTZO, RHO_HEMATITA]) for _ in range(n)])
    psi = np.array([rng.uniform(0.4, 1.0) for _ in range(n)])

    p["d"] = d
    p["rho"] = rho
    p["psi"] = psi
    p["vt"] = np.array([vt_intermediario(rho[i], d[i], psi[i]) for i in range(n)])
    p["re"] = np.array([reynolds(p["vt"][i], d[i]) for i in range(n)])

    p["x"] = np.array([rng.uniform(-0.45, 0.45) * LADO for _ in range(n)])
    p["z"] = np.array([rng.uniform(-0.45, 0.45) * LADO for _ in range(n)])
    p["y"] = np.array([rng.uniform(0.05, 0.95) * ALTURA for _ in range(n)])

    p["hematita"] = rho > 3000.0
    p["achatada"] = psi < 0.65
    p["rng"] = rng
    return p


def passo(p):
    """
    Um passo de integracao.

    Velocidade vertical: vy = v_fluido - VT  (fluido sobe, particula afunda).
    Dispersao lateral: proporcional a (1 - psi) -- particula irregular desvia
    mais, como no b.py.

    Condicoes de contorno (b.py):
      * paredes laterais x e z: periodicas (a particula reentra pelo lado oposto);
      * base  (y < 0)      : particula sedimentou -> reinjetada no topo;
      * topo  (y > ALTURA) : particula foi arrastada -> reinjetada na base.
    """
    rng = p["rng"]
    n = len(p["d"])

    vy = V_FLUIDO - p["vt"]
    desvio = (1.0 - p["psi"]) * 0.4 * V_FLUIDO
    vx = np.array([rng.uniform(-1.0, 1.0) for _ in range(n)]) * desvio
    vz = np.array([rng.uniform(-1.0, 1.0) for _ in range(n)]) * desvio

    p["x"] += vx * DT
    p["y"] += vy * DT
    p["z"] += vz * DT

    meia = 0.5 * LADO
    p["x"] = np.where(p["x"] > meia, p["x"] - LADO, p["x"])
    p["x"] = np.where(p["x"] < -meia, p["x"] + LADO, p["x"])
    p["z"] = np.where(p["z"] > meia, p["z"] - LADO, p["z"])
    p["z"] = np.where(p["z"] < -meia, p["z"] + LADO, p["z"])

    sedimentou = p["y"] < 0.0
    arrastou = p["y"] > ALTURA

    for idx in np.flatnonzero(sedimentou):
        p["y"][idx] = ALTURA
        p["x"][idx] = rng.uniform(-0.45, 0.45) * LADO
        p["z"][idx] = rng.uniform(-0.45, 0.45) * LADO

    for idx in np.flatnonzero(arrastou):
        p["y"][idx] = 0.0
        p["x"][idx] = rng.uniform(-0.45, 0.45) * LADO
        p["z"][idx] = rng.uniform(-0.45, 0.45) * LADO

    return int(sedimentou.sum()), int(arrastou.sum())


# =====================================================================
# 7. FIGURA
# =====================================================================

def _estilo_eixo_3d(ax):
    ax.set_facecolor(COR_FUNDO)
    for eixo in (ax.xaxis, ax.yaxis, ax.zaxis):
        eixo.set_pane_color((0.07, 0.07, 0.07, 1.0))
        eixo._axinfo["grid"]["color"] = (1, 1, 1, 0.06)
        eixo.line.set_color(INK_3)
    ax.tick_params(colors=INK_2, labelsize=7.5)


def _caixa_coluna(ax):
    """Arestas da coluna, em mm."""
    a = 0.5 * LADO * 1000.0
    h = ALTURA * 1000.0
    base = [(-a, -a), (a, -a), (a, a), (-a, a), (-a, -a)]

    for y in (0.0, h):
        xs = [pt[0] for pt in base]
        zs = [pt[1] for pt in base]
        ax.plot(xs, zs, [y] * 5, color=COR_FLUXO, lw=1.0, alpha=0.45)

    for (x, z) in base[:4]:
        ax.plot([x, x], [z, z], [0.0, h], color=COR_FLUXO, lw=0.8, alpha=0.25)


def _setas_fluxo(ax):
    """Fluxo ascendente de 6 LPM na base da coluna."""
    a = 0.5 * LADO * 1000.0
    passo_grade = np.linspace(-a * 0.7, a * 0.7, 3)
    for x in passo_grade:
        for z in passo_grade:
            ax.quiver(x, z, 2.0, 0, 0, 22.0, color=COR_FLUXO,
                      alpha=0.55, linewidth=0.9, arrow_length_ratio=0.35)


def _painel_vt(ax, d50_q, d50_h):
    """Painel comparativo: VT(d) pelas tres leis (quartzo) + Cd(Re) da hematita."""
    ax.set_facecolor(COR_PAINEL)

    d = np.linspace(5e-6, 2e-3, 400)
    y_st = np.array([vt_stokes(RHO_QUARTZO, di) for di in d]) * 1000.0
    y_nw = np.array([vt_newton(RHO_QUARTZO, di) for di in d]) * 1000.0
    y_in = np.array([vt_intermediario(RHO_QUARTZO, di) for di in d]) * 1000.0
    y_in_h = np.array([vt_intermediario(RHO_HEMATITA, di) for di in d]) * 1000.0
    d_um = d * 1e6

    ax.loglog(d_um, y_st, color=COR_STOKES, lw=2.0, label="Stokes — quartzo")
    ax.loglog(d_um, y_nw, color=COR_NEWTON, lw=2.0, label="Newton — quartzo")
    ax.loglog(d_um, y_in, color=COR_INTERM, lw=2.0,
              label="Cd(Re) — quartzo 2650")
    ax.loglog(d_um, y_in_h, color=COR_INTERM, lw=2.0, ls=(0, (5, 2)),
              label="Cd(Re) — hematita 5200")

    u = V_FLUIDO * 1000.0
    ax.axhline(u, color=INK_2, lw=1.2, ls=(0, (4, 3)), alpha=0.9)
    ax.text(d_um[-1], u * 1.15, f"u = {u:.0f} mm/s  ({Q_LPM:.0f} LPM)",
            color=INK_2, fontsize=7.5, ha="right", va="bottom")

    ax.axvspan(D_MIN * 1e6, D_MAX * 1e6, color=INK_2, alpha=0.07)
    ax.text(math.sqrt(D_MIN * D_MAX) * 1e6, 4.5e2,
            "faixa simulada", color=INK_3, fontsize=7, ha="center", va="top")

    for d50 in (d50_q, d50_h):
        ax.plot([d50 * 1e6], [u], marker="o", ms=6, color=COR_INTERM,
                mec=COR_PAINEL, mew=1.5, zorder=5)
    ax.text(d50_q * 1e6 * 1.2, u * 0.45, f"d50 quartzo\n{d50_q*1e6:.0f} µm",
            color=INK_2, fontsize=7, ha="left", va="top")
    ax.text(d50_h * 1e6 * 0.8, u * 0.45, f"d50 hematita\n{d50_h*1e6:.0f} µm",
            color=INK_2, fontsize=7, ha="right", va="top")

    ax.set_xlabel("diâmetro da partícula  [µm]", color=INK_2, fontsize=8)
    ax.set_ylabel("velocidade terminal  [mm/s]", color=INK_2, fontsize=8)
    ax.set_title("Velocidade terminal (partícula esférica, ψ = 1)",
                 color=INK_1, fontsize=9, loc="left", pad=6)
    ax.tick_params(colors=INK_2, labelsize=7.5, which="both")
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(INK_3)
        ax.spines[lado].set_linewidth(0.8)
    ax.grid(True, which="major", color=INK_1, alpha=0.08, lw=0.7)
    ax.set_ylim(1e-2, 1e3)

    leg = ax.legend(loc="upper left", frameon=False, fontsize=7.5,
                    labelcolor=INK_2, handlelength=1.6, borderaxespad=0.2)
    for texto in leg.get_texts():
        texto.set_color(INK_2)


def construir_figura(p):
    fig = plt.figure(figsize=(11.2, 6.3), dpi=86, facecolor=COR_FUNDO)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1.0],
                          height_ratios=[1.0, 0.58],
                          left=0.02, right=0.965, top=0.82, bottom=0.06,
                          wspace=0.18, hspace=0.55)

    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_vt = fig.add_subplot(gs[0, 1])
    ax_txt = fig.add_subplot(gs[1, 1])
    ax_txt.axis("off")

    _estilo_eixo_3d(ax3d)
    a = 0.5 * LADO * 1000.0
    ax3d.set_xlim(-a, a)
    ax3d.set_ylim(-a, a)
    ax3d.set_zlim(0, ALTURA * 1000.0)
    ax3d.set_box_aspect((1, 1, 1.7))
    ax3d.set_xlabel("x [mm]", color=INK_2, fontsize=8, labelpad=-2)
    ax3d.set_ylabel("z [mm]", color=INK_2, fontsize=8, labelpad=-2)
    ax3d.set_zlabel("altura y [mm]", color=INK_2, fontsize=8, labelpad=-2)

    _caixa_coluna(ax3d)
    _setas_fluxo(ax3d)

    tamanho = 6.0 + 46.0 * (p["d"] - D_MIN) / (D_MAX - D_MIN)
    grupos = {}
    for nome, mask, cor, marcador in (
        ("quartzo esférico", (~p["hematita"]) & (~p["achatada"]), COR_QUARTZO, "o"),
        ("quartzo achatado", (~p["hematita"]) & p["achatada"], COR_QUARTZO, "D"),
        ("hematita esférica", p["hematita"] & (~p["achatada"]), COR_HEMATITA, "o"),
        ("hematita achatada", p["hematita"] & p["achatada"], COR_HEMATITA, "D"),
    ):
        sc = ax3d.scatter(p["x"][mask] * 1000.0, p["z"][mask] * 1000.0,
                          p["y"][mask] * 1000.0,
                          s=tamanho[mask], c=cor, marker=marcador,
                          edgecolors=COR_FUNDO, linewidths=0.4,
                          depthshade=True, label=nome)
        grupos[nome] = (sc, mask)

    leg = ax3d.legend(loc="upper left", bbox_to_anchor=(-0.02, 0.98),
                      frameon=False, fontsize=7.5, scatterpoints=1,
                      handletextpad=0.4, labelspacing=0.35)
    for texto in leg.get_texts():
        texto.set_color(INK_2)

    d50_q = diametro_de_corte(RHO_QUARTZO)
    d50_h = diametro_de_corte(RHO_HEMATITA)
    _painel_vt(ax_vt, d50_q, d50_h)

    fig.text(0.02, 0.955,
             f"Elutriação em coluna ascendente — {Q_LPM:.0f} LPM · seção "
             f"{LADO*1000:.0f} × {LADO*1000:.0f} mm · H = {ALTURA*1000:.0f} mm · "
             f"água a 20 °C",
             color=INK_1, fontsize=12.5, ha="left", va="center")
    fig.text(0.02, 0.918,
             "Arrasto Cd(Re) = (0,63 + 4,8/√Re)² · correção de forma Cd_ef = Cd/ψ",
             color=INK_3, fontsize=8.5, ha="left", va="center")
    fig.text(0.02, 0.888,
             "Contorno: paredes laterais periódicas · base e topo com reinjeção",
             color=INK_3, fontsize=8.5, ha="left", va="center")

    return fig, ax3d, ax_txt, grupos, tamanho, (d50_q, d50_h)


# =====================================================================
# 8. ANIMACAO E EXPORTACAO DO GIF
# =====================================================================

def _quadro_para_imagem(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return Image.fromarray(buf[:, :, :3].copy())


def simular(saida_gif="elutriacao_3d.gif", cores_gif=96):
    p = criar_particulas()
    fig, ax3d, ax_txt, grupos, tamanho, (d50_q, d50_h) = construir_figura(p)

    elutr_q = float(np.mean(p["vt"][~p["hematita"]] < V_FLUIDO)) * 100.0
    elutr_h = float(np.mean(p["vt"][p["hematita"]] < V_FLUIDO)) * 100.0

    caixa_txt = ax_txt.text(
        0.0, 1.06, "", color=INK_2, fontsize=9, family="monospace",
        ha="left", va="top", transform=ax_txt.transAxes, linespacing=1.6)

    quadros = []
    n_sed = n_arr = 0

    for k in range(N_PASSOS + 1):
        if k > 0:
            s, a = passo(p)
            n_sed += s
            n_arr += a

        if k % PASSOS_POR_QUADRO:
            continue

        t = k * DT
        for nome, (sc, mask) in grupos.items():
            sc._offsets3d = (p["x"][mask] * 1000.0,
                             p["z"][mask] * 1000.0,
                             p["y"][mask] * 1000.0)

        ax3d.view_init(elev=16.0, azim=-62.0 + 0.28 * (k / PASSOS_POR_QUADRO))

        caixa_txt.set_text(
            f"t = {t:5.2f} s   de {T_TOTAL:.1f} s   (dt = {DT*1000:.0f} ms)\n"
            f"u  = {V_FLUIDO*1000:5.2f} mm/s   Q = {Q_LPM:.1f} L/min\n"
            f"\n"
            f"arrastadas pelo topo .... {n_arr:4d}\n"
            f"sedimentadas na base .... {n_sed:4d}\n"
            f"\n"
            f"VT < u  quartzo  {elutr_q:5.1f} %\n"
            f"VT < u  hematita {elutr_h:5.1f} %\n"
            f"Re da população: {p['re'].min():.2g} – {p['re'].max():.2g}"
        )

        quadros.append(_quadro_para_imagem(fig))

    plt.close(fig)

    base = quadros[0].quantize(colors=cores_gif, method=Image.MEDIANCUT)
    paletizados = [q.quantize(palette=base, dither=Image.FLOYDSTEINBERG)
                   for q in quadros]

    destino = os.path.join(os.path.dirname(os.path.abspath(__file__)), saida_gif)
    paletizados[0].save(
        destino, save_all=True, append_images=paletizados[1:],
        duration=int(1000 / FPS), loop=0, optimize=True)

    return destino, len(quadros), n_sed, n_arr, d50_q, d50_h


def resumo_console(d50_q, d50_h):
    print(f"Coluna .............. {LADO*1000:.0f} x {LADO*1000:.0f} x "
          f"{ALTURA*1000:.0f} mm")
    print(f"Vazao ............... {Q_LPM:.1f} L/min = {Q:.2e} m3/s")
    print(f"Velocidade ascend. .. {V_FLUIDO*1000:.2f} mm/s")
    print(f"d50 quartzo (psi=1) . {d50_q*1e6:.1f} um")
    print(f"d50 hematita (psi=1)  {d50_h*1e6:.1f} um")
    print()
    print(" d [um]   rho      VT Stokes   VT Newton   VT Cd(Re)   Re")
    for d in (20e-6, 50e-6, 100e-6, 500e-6):
        for rho in (RHO_QUARTZO, RHO_HEMATITA):
            vi = vt_intermediario(rho, d)
            print(f"{d*1e6:7.0f}{rho:8.0f}{vt_stokes(rho,d)*1000:12.3f}"
                  f"{vt_newton(rho,d)*1000:12.1f}{vi*1000:12.3f}"
                  f"{reynolds(vi,d):9.3f}")


if __name__ == "__main__":
    destino, n_quadros, n_sed, n_arr, d50_q, d50_h = simular()
    resumo_console(d50_q, d50_h)
    print()
    print(f"GIF: {destino}  ({n_quadros} quadros, {FPS} fps)")
    print(f"Balanco em {T_TOTAL:.1f} s: {n_arr} arrastadas / {n_sed} sedimentadas")
