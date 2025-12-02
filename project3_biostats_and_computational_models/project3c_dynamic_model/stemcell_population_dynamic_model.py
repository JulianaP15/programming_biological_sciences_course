"""
Script: stemcell_population_dynamic_model.py
Author: Juliana Patiño Gallego
Description:
    Simulation of a 3-compartment stem cell system using ordinary differential equations (ODEs).
    The model represents hierarchical tissue organization with:
        - S(t): Stem cells
        - P(t): Progenitor cells
        - D(t): Differentiated cells

    The script solves the system under two biological scenarios:
        1. Near-homeostasis (balanced self-renewal, differentiation, and loss)
        2. Dysregulated dynamics (increased self-renewal and proliferation, reduced differentiation)

    It generates time-course plots for both scenarios and saves them as PNG files.

Model:
    The system is governed by:

        dS/dt = r_s * S * (1 - S/K) - α*S - δ_s*S
        dP/dt = α*S + r_p*P - β*P - δ_p*P
        dD/dt = β*P - δ_d*D

    where:
        r_s      = stem-cell self-renewal rate
        K        = carrying capacity of stem-cell niche
        α        = differentiation rate from S → P
        δ_s      = stem-cell loss / death rate
        r_p      = proliferation rate of progenitors
        β        = differentiation rate from P → D
        δ_p      = progenitor death rate
        δ_d      = differentiated-cell turnover rate

Features:
    - Flexible parameter dictionaries for biological scenarios
    - Automatic figure saving at high resolution (300 DPI)
    - Linear-scale plots of the full dynamics and early-time zoom
    - Clean seaborn/Matplotlib workflow
    - Modular code with helper functions for saving and parameter handling

Output:
    Figures saved to ./figures :
        - SPD-hom.png       (Homeostatic scenario)
        - SPD-dys-100.png   (Dysregulated scenario, full time)
        - SPD-dys-10.png    (Early dysregulated dynamics, first 10 units)

Notes / assumptions:
    * Time units are arbitrary (dimensionless)
    * Initial hierarchy assumes S < P < D
    * dysregulated parameters mimic early tumor-like behavior (expansion of S and P)
    * Uses scipy.integrate.odeint for numerical integration

"""


# ==== Library imports ====
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import odeint
import os

# ==== Configuration ====
# Output folder for saved figures
output_dir = "figures"

# ==== Functions ====
# Save figure
def save_figure(fig, name, output_dir):
    """
    Save a Matplotlib figure as .png

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure object to be saved.
    name : str
        Base name of the output file (without extension).
    output_dir : str, optional
        Path to the output directory where figures will be stored.

    Notes
    -----
    * Uses tight bounding boxes to avoid excessive white margins.
    * Saves 300 DPI PNG.
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    # Define output path
    path_png = os.path.join(output_dir, f"{name}.png")
    # Save
    fig.savefig(path_png, format="png", bbox_inches="tight", dpi=300)

# Params dictionary to tuple
def params_to_tuple(pdict):
    return (pdict["r_s"], pdict["K"], pdict["alpha"], pdict["delta_s"],
            pdict["r_p"], pdict["beta"], pdict["delta_p"],
            pdict["delta_d"])

# ODE system
def stem_cell_system(y, t, params):
    """
    Simple 3-compartment model:
    S(t): stem cells
    P(t): progenitor cells
    D(t): differentiated cells
    """
    S, P, D = y
    (r_s, K, alpha, delta_s,
     r_p, beta, delta_p,
     delta_d) = params

    # dS/dt: self-renewal (logistic) - differentiation - death
    dSdt = r_s * S * (1 - S / K) - alpha * S - delta_s * S

    # dP/dt: input from S + proliferation - differentiation - death
    dPdt = alpha * S + r_p * P - beta * P - delta_p * P

    # dD/dt: input from P - loss of differentiated cells
    dDdt = beta * P - delta_d * D

    return np.array([dSdt, dPdt, dDdt])

# ==== Define two scenarios: "homeostasis" vs "dysregulated" ====
# Scenario 1: near-homeostatic tissue
params_homeo = {
    "r_s": 0.5,       # stem-cell self-renewal rate (moderate)
    "K": 1000.0,      # carrying capacity of the stem-cell niche
    "alpha": 0.3,     # differentiation rate S -> P (high enough to prevent overgrowth)
    "delta_s": 0.1,   # stem-cell death / loss rate

    "r_p": 0.2,       # progenitor proliferation rate (low; prevents exponential expansion)
    "beta": 0.25,     # differentiation P -> D (slightly higher than r_p)
    "delta_p": 0.05,  # progenitor death rate (helps stabilize P levels)

    "delta_d": 0.2    # death / turnover rate of differentiated cells
}

# Scenario 2: dysregulated stem cell compartment (higher self-renewal, lower differentiation)
params_dys = {
    "r_s": 0.8,       # stem-cell self-renewal rate (much higher → overgrowth)
    "K": 2000.0,      # expanded niche capacity (e.g., microenvironment permissive)

    "alpha": 0.1,     # differentiation S -> P (reduced; typical in dysregulated tissues)
    "delta_s": 0.02,  # lower stem-cell death → accumulation of S

    "r_p": 0.5,       # high proliferation rate of progenitors (drives exponential expansion)
    "beta": 0.1,      # low differentiation P -> D (cells stay proliferative)
    "delta_p": 0.02,  # low progenitor death → amplifies overgrowth

    "delta_d": 0.1    # normal turnover of differentiated cells
}


# ==== Initial conditions and time span ====
# Initial cell counts
# These values reflect a biologically plausible hierarchy: S < P < D
S0 = 100.0
P0 = 200.0
D0 = 500.0
y0 = [S0, P0, D0]  # initial state vector passed to the ODE system

Nt = 500   # number of time points for the simulation
tmax = 100.0   # total duration of the simulation (arbitrary units)
t = np.linspace(0.0, tmax, Nt)  # create an evenly spaced time vector where the system will be evaluated


# ==== Solve ODEs for both scenarios ====
res_homeo = odeint(stem_cell_system, y0, t,
                   args=(params_to_tuple(params_homeo),))

res_dys = odeint(stem_cell_system, y0, t,
                 args=(params_to_tuple(params_dys),))

# Extract each compartment from the resulting solution
S_home, P_home, D_home = res_homeo.T
S_dys,  P_dys,  D_dys  = res_dys.T


# ==== Plots ====
# ---- Homeostatic scenario ----
plt.figure(figsize=(8,5))
sns.lineplot(x=t, y=S_home, label="Stem cells (S)", linewidth=2)
sns.lineplot(x=t, y=P_home, label="Progenitors (P)", linewidth=2)
sns.lineplot(x=t, y=D_home, label="Differentiated (D)", linewidth=2)

plt.title("Stem cell system – Near homeostasis")
plt.xlabel("Time")
plt.ylabel("Cell counts")
plt.grid(alpha=0.3)
plt.tight_layout()

fig = plt.gcf()  # get the current active figure
save_figure(fig, "SPD-hom", output_dir)

plt.show()

# ---- Dysregulated scenario ----
plt.figure(figsize=(8,5))
sns.lineplot(x=t, y=S_dys, label="Stem cells (S)", linewidth=2)
sns.lineplot(x=t, y=P_dys, label="Progenitors (P)", linewidth=2)
sns.lineplot(x=t, y=D_dys, label="Differentiated (D)", linewidth=2)

plt.title("Stem cell system – Dysregulated")
plt.xlabel("Time")
plt.ylabel("Cell counts")
plt.grid(alpha=0.3)
plt.tight_layout()

fig = plt.gcf()  # get the current active figure
save_figure(fig, "SPD-dys-100", output_dir)

plt.show()

# Mask to select only time points up to t = 20
mask = t <= 10.0

# Slice time and solutions
t_10   = t[mask]
S_10   = S_dys[mask]
P_10   = P_dys[mask]
D_10   = D_dys[mask]

plt.figure(figsize=(8, 5))

# Plot only the early-time dynamics in linear scale
sns.lineplot(x=t_10, y=S_10, label="Stem cells (S)", linewidth=2)
sns.lineplot(x=t_10, y=P_10, label="Progenitors (P)", linewidth=2)
sns.lineplot(x=t_10, y=D_10, label="Differentiated (D)", linewidth=2)

plt.title("Stem cell system – Dysregulated (first 10 time units)")
plt.xlabel("Time (arbitrary units)")
plt.ylabel("Cell counts")
plt.grid(alpha=0.3)
plt.tight_layout()

fig = plt.gcf()  # get the current active figure
save_figure(fig, "SPD-dys-10", output_dir)

plt.show()