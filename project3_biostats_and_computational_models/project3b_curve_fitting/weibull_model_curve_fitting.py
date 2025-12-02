"""
Script: weibull_model_curve_fitting.py
Author: Juliana Patiño Gallego

Description:
    This script simulates a "true" three-parameter Weibull survival curve
    and evaluates how accurately its parameters can be recovered under
    different levels of Gaussian noise.

    For each noise level tested, the script:
        1. Generates noisy survival observations based on a known
           ground-truth Weibull model.
        2. Fits the model using three strategies:
               - With informed initial guesses (stable baseline)
               - Without initial guesses but WITH parameter bounds
               - Without initial guesses AND WITHOUT bounds (unstable baseline)
        3. Computes predicted curves, residuals, and R² values.
        4. Produces comparative survival plots and residual plots.
        5. Stores all estimated parameters in a summary table
           exported as a CSV file.

Model:
    Weibull survival function with time delay:
        S(t) = exp( - ((max(t - delay, 0) / lambda) ** k) )

    Parameters:
        - lambda (scale): controls curve width
        - k (shape): determines hazard trend
        - delay: latency before survival starts decreasing

Outputs:
    * High-resolution PNG figures for:
         - Model fits vs ground truth (one per noise level)
         - Residuals for all fitting strategies
    * A 'param_summary.csv' file summarizing:
         - True parameters
         - Estimated parameters for each method
         - R² for each method
         - Convergence status

Notes / assumptions:
    * curve_fit is sensitive to starting conditions; this script
      demonstrates how initialization and parameter bounds affect
      convergence, stability, and accuracy.
    * Noise is simulated using Gaussian distributions with increasing
      standard deviations to test robustness.
    * Survival values are clipped to the interval [0,1] to maintain
      biological interpretability.
    * The analysis is intended for instructional purposes in the context
      of survival modeling and parameter recovery under noise.
"""

# ==== Library imports ====
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
import os


# ==== Functions ====
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


# ==== Configuration ====
# Output folder for saved figures
output_dir = "survival_curves"


# ==== Global aesthetics ====
# Base palette
color_true = "#000000"  # black charcoal
color_init = "gray"  # gray
color_noinit = "#E69F00"  # orange
color_nobounds = "red"  # red
color_observed = "#999999"  # neutral gray

# Style: light grid, no top/right spines
custom_params = {'axes.spines.right': False, 'axes.spines.top': False}
sns.set_theme(style='ticks', rc=custom_params)


# ==== Define Weibull survival model (3 parameters) ====
def weibull_survival(t, lambda_scale, k_shape, t_delay):
    """
    Three-parameter Weibull survival model with a time delay.

    S(t) = exp( - ((max(t - delay, 0) / lambda)^k) )
    """
    t_effective = np.maximum(t - t_delay, 0)
    return np.exp(- (t_effective / lambda_scale)**k_shape)


# ==== Generate synthetic "true" survival curve ====
# High-resolution time axis
time_months = np.linspace(0, 36, 200)  # 200 time points between 0 and 36 months

# True (underlying) Weibull parameters
true_lambda = 15.0   # scale (months)
true_k      = 1.5    # shape (>1 = increasing hazard)
true_delay  = 2.0    # delay (months)

S_true = weibull_survival(time_months, true_lambda, true_k, true_delay)


# ==== Define noise levels to test ====
noise_levels = [0.05, 0.15, 0.25, 0.50]  # standard deviation of Gaussian noise
results = []  # to store summary per noise level


# ==== Fit model with different noise level
for sigma in noise_levels:
    # ---- Generate noisy observations ----
    noise = np.random.normal(loc=0.0, scale=sigma, size=time_months.size)
    S_observed = np.clip(S_true + noise, 0, 1)  # keep in [0,1]

    # Common bounds for both fits
    bounds = ([1e-3, 0.1, -5.0],   # lower bounds [lambda, k, delay]
              [100.0, 5.0, 20.0])  # upper bounds

    # ---- Fit WITH informed initial guesses ----
    # Try to give the optimizer a reasonable starting point
    lambda_guess = 10.0    # around the middle of follow-up
    k_guess      = 1.0     # close to exponential shape
    delay_guess  = 0.0     # no delay initially

    initial_guess = [lambda_guess, k_guess, delay_guess]

    params_with_init, cov_with_init = curve_fit(weibull_survival, time_months, S_observed,
                                                p0=initial_guess,
                                                bounds=bounds
                                                )

    lambda_est_init, k_est_init, delay_est_init = params_with_init

    # Predicted survival and R² with initial guesses
    S_pred_init = weibull_survival(time_months, lambda_est_init, k_est_init, delay_est_init)
    residuals_init = S_observed - S_pred_init
    ss_res_init = np.sum(residuals_init**2)
    ss_tot = np.sum((S_observed - S_observed.mean())**2)
    R2_init = 1 - ss_res_init / ss_tot

    # ---- Fit WITHOUT providing initial guesses ----
    # Here curve_fit will choose default initial parameters.
    # - Without bounds, it uses p0 = 1 for each parameter.
    # - With bounds, it uses the middle of the bounds as starting point.
    # We keep the same bounds for a fair comparison.
    try:
        params_no_init, cov_no_init = curve_fit(weibull_survival, time_months, S_observed,
                                                bounds=bounds
                                                )

        lambda_est_ni, k_est_ni, delay_est_ni = params_no_init

        S_pred_ni = weibull_survival(time_months, lambda_est_ni, k_est_ni, delay_est_ni)
        residuals_ni = S_observed - S_pred_ni
        ss_res_ni = np.sum(residuals_ni**2)
        R2_ni = 1 - ss_res_ni / ss_tot

        success_no_init = True

    except RuntimeError:
        # Fitting failed to converge
        lambda_est_ni = np.nan
        k_est_ni = np.nan
        delay_est_ni = np.nan
        R2_ni = np.nan
        success_no_init = False

    # ---- Fit WITHOUT bounds AND WITHOUT initial guesses ----
    # This is the "raw" fit: no parameter restrictions and no starting values.
    # curve_fit will start from p0 = 1 for all parameters by default.
    try:
        # No p0, no bounds → fully unconstrained optimization
        params_no_bounds, cov_no_bounds = curve_fit(weibull_survival, time_months, S_observed)

        lambda_est_nb, k_est_nb, delay_est_nb = params_no_bounds

        S_pred_nb = weibull_survival(time_months, lambda_est_nb, k_est_nb, delay_est_nb)
        residuals_nb = S_observed - S_pred_nb
        ss_res_nb = np.sum(residuals_nb ** 2)
        R2_nb = 1 - ss_res_nb / ss_tot

        success_no_bounds = True

    except (RuntimeError, ValueError, OverflowError):
        # Fitting fails more often here
        lambda_est_nb = np.nan
        k_est_nb = np.nan
        delay_est_nb = np.nan
        R2_nb = np.nan
        success_no_bounds = False

    # ---- Store summary for this noise level ----
    results.append({
        "noise_sigma": sigma,

        "lambda_true": true_lambda,
        "k_true": true_k,
        "delay_true": true_delay,

        # With init
        "lambda_init": lambda_est_init,
        "k_init": k_est_init,
        "delay_init": delay_est_init,
        "R2_init": R2_init,

        # Without init, with bounds
        "lambda_no_init": lambda_est_ni,
        "k_no_init": k_est_ni,
        "delay_no_init": delay_est_ni,
        "R2_no_init": R2_ni,
        "success_no_init": success_no_init,

        # NEW: Without init AND without bounds
        "lambda_nobounds": lambda_est_nb,
        "k_nobounds": k_est_nb,
        "delay_nobounds": delay_est_nb,
        "R2_nobounds": R2_nb,
        "success_nobounds": success_no_bounds
    })

    # ---- Survival curves for this noise level ----
    plt.figure(figsize=(7, 5))

    # True curve
    sns.lineplot(x=time_months, y=S_true,
                 label="True survival", linewidth=2, color=color_true)

    # Observed noisy data
    sns.scatterplot(x=time_months, y=S_observed,
                    s=12, alpha=0.5, color=color_observed)

    # Fit with initial guesses
    sns.lineplot(x=time_months, y=S_pred_init,
                 linewidth=2.5, label="Fit (with init)", color=color_init)

    # Fit without init (bounded)
    if success_no_init:
        sns.lineplot(x=time_months, y=S_pred_ni,
                     linestyle="--", label="Fit (bounded, no init)", linewidth=2.5, color=color_noinit)

    # Fit without bounds
    if success_no_bounds:
        sns.lineplot(x=time_months, y=S_pred_nb,
                     linestyle=":", label="Fit (no bounds, no init)", linewidth=3, color=color_nobounds)

    # Build R² text block (only for successful fits)
    r2_lines = [f"R² (with init) = {R2_init:.3f}"]
    if success_no_init:
        r2_lines.append(f"R² (bounded, no init) = {R2_ni:.3f}")
    if success_no_bounds:
        r2_lines.append(f"R² (no bounds) = {R2_nb:.3f}")

    r2_text = "\n".join(r2_lines)

    ax = plt.gca()
    ax.text(
        0.02, 0.02,
        r2_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3",
                  facecolor="white",
                  edgecolor="grey",
                  alpha=0.85)
    )

    plt.title(f"Weibull Fit Comparison – Noise σ = {sigma}")
    plt.xlabel("Time (months)")
    plt.ylabel("Survival Probability S(t)")
    plt.ylim(-0.05, 1.05)
    plt.legend()
    plt.tight_layout()

    # Save figure
    fig = plt.gcf()  # get the current active figure
    save_figure(fig, f"sigma={sigma} curve", output_dir)

    plt.show()


    # ---- Residual plots for all methods for this noise level ----
    plt.figure(figsize=(8, 5))

    # Residuals for with init
    sns.scatterplot(x=time_months, y=residuals_init,
                    label="Residuals (with init)", s=15, alpha=0.6, color=color_init)

    # Residuals for no init (bounded)
    if success_no_init:
        sns.scatterplot(x=time_months, y=residuals_ni,
                        label="Residuals (no init, bounded)", s=15, alpha=0.6, color=color_noinit)

    # Residuals for no bounds
    if success_no_bounds:
        sns.scatterplot(x=time_months, y=residuals_nb,
                        label="Residuals (no bounds, no init)", s=15, alpha=0.6, color=color_nobounds)

    # Add reference line at zero residuals
    plt.axhline(0, color="#333333", linewidth=1, linestyle="--")

    plt.title(f"Residuals Comparison – Noise σ = {sigma}")
    plt.xlabel("Time (months)")
    plt.ylabel("Residual (Observed – Predicted)")
    plt.legend(loc='lower right', framealpha=0.9)
    plt.ylim(-1, 1)
    plt.tight_layout()

    # Save figure
    fig = plt.gcf()  # get the current active figure
    save_figure(fig, f"sigma={sigma} residuals", output_dir)

    plt.show()


# ==== Summary table as a DataFrame ====
results_df = pd.DataFrame(results)

# Parameter summary
param_summary = results_df[[
    "noise_sigma",

    # True (for reference); With initial guesses; Without init, bounded; Without bounds
    # Lambda
    "lambda_true", "lambda_init", "lambda_no_init", "lambda_nobounds",
    # k
    "k_true", "k_init", "k_no_init", "k_nobounds",
    # Delay
    "delay_true", "delay_init", "delay_no_init", "delay_nobounds",
    # R²
    "R2_init", "R2_no_init", "R2_nobounds"

]]

print("\n===== COMPACT PARAMETER SUMMARY =====\n")
print(param_summary)
param_summary.to_csv("param_summary.csv", index=False)


