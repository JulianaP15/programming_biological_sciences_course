# Project 3b – Weibull Survival Model Curve Fitting

## Description
This project implements a **three-parameter Weibull survival model** to explore how well its parameters can be recovered under different levels of noise.  
The script simulates a ground-truth survival curve, adds controlled Gaussian noise, and compares three fitting strategies:

- **Fit with informed initial guesses**
- **Fit without initial guesses but with parameter bounds**
- **Fit without initial guesses and without bounds**

The goal is to evaluate **model stability**, **sensitivity to noise**, and the **importance of initialization and bounds** in nonlinear optimization.

---

## Features

- **Simulates a ground-truth Weibull survival curve** with time delay.
- **Adds Gaussian noise** at four noise levels to test robustness.
- **Compares three fitting strategies**:
  - With informed initial parameters
  - Without initial parameters (bounded)
  - Without initial parameters (unbounded)
- **Computes and visualizes**:
  - Fitted survival curves
  - Residual distributions
  - R² metrics
- **Automatically saves figures** (`.png`) for each noise level.
- **Exports a parameter summary** (`param_summary.csv`) containing:
  - True parameters
  - Estimated parameters
  - Convergence flags
  - R² for each method

---

## Requirements

- **Python 3.8+**
- Python libraries:
  - `numpy`
  - `pandas`
  - `scipy`
  - `matplotlib`
  - `seaborn`

Install all dependencies with:

```bash
pip install numpy pandas scipy matplotlib seaborn
```

---

## Project Files

- `weibull_survival_analysis.py`  
  Main script that performs simulation, model fitting, plotting, and summary export.

- `survival_curves/`  
  Output directory (auto-created) containing:
  - Survival comparison plots  
  - Residual plots  

- `param_summary.csv`  
  Summary table with estimated parameters and R² values for each noise level.

---

## Model Description

The survival function used is a **three-parameter Weibull model with delay**, defined as:

$$
S(t) = \exp\left[-\left(\frac{\max(t - \text{delay},\, 0)}{\lambda}\right)^{k}\right]
$$

Where:

- `lambda` → scale (controls decay rate)  
- `k` → shape (k > 1 means increasing hazard)  
- `delay` → latency before decay begins  

This model captures scenarios with **latency**, **accelerating hazard**, and **smooth survival dynamics**.

---

## Usage

Run the script directly from a terminal:

```bash
python weibull_survival_analysis.py
```

The script will:

1. Simulate a true survival curve  
2. Add noisy observations  
3. Fit the model using three strategies  
4. Generate plots for each noise level  
5. Export results into `param_summary.csv`

---

## Output Files

### Figures
Stored inside the `survival_curves/` folder:

- `sigma=X curve.png`  
  Comparison of:
  - Ground-truth curve  
  - Noisy observations  
  - Fitted curves from all strategies

- `sigma=X residuals.png`  
  Residual scatterplot comparing all fitting methods.

### Summary Table
- `param_summary.csv`  
  Contains:
  - Noise level  
  - True λ, k, and delay  
  - Estimated parameters (all methods)  
  - R² values  
  - Convergence flags

---

## Interpretation Highlights

- The **fit with informed initial guesses** is consistently the most stable and accurate.
- The **bounded fit without initial guesses** performs reasonably well but deteriorates under high noise.
- The **unbounded fit** often fails to converge, produces unrealistic parameter values, or shows poor curve behavior.
- Residual plots reveal systematic deviations and help diagnose model instability.
- Increasing noise reduces R² and amplifies parameter uncertainty.

---

## License

This project is distributed under the **MIT License**.

---

## Author

Juliana Patiño Gallego  
jpatinoga@unal.edu.co
