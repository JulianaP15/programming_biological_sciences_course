# Project 3c – Stem Cell Population Dynamic Model (S–P–D System)

## Description
This project implements and simulates a **three-compartment ordinary differential equation (ODE) model** representing a hierarchical tissue composed of:

- **Stem cells (S)**
- **Progenitor cells (P)**
- **Differentiated cells (D)**

The model evaluates how tissues behave under two contrasting biological scenarios:

1. **Near-homeostasis** (balanced renewal and differentiation)
2. **Dysregulated dynamics** (increased self-renewal, reduced differentiation, high proliferation)

The script solves the ODE system, visualizes the time-course behavior, and saves high-quality figures for interpretation.

---
## Model

The system is described by the following ODEs:

$$
\frac{dS}{dt} = r_s S\left(1 - \frac{S}{K}\right) - \alpha S - \delta_s S
$$

$$
\frac{dP}{dt} = \alpha S + r_p P - \beta P - \delta_p P
$$

$$
\frac{dD}{dt} = \beta P - \delta_d D
$$

Where:

- \( r_s \): stem-cell self-renewal rate  
- \( K \): stem-cell niche carrying capacity  
- \( \alpha \): differentiation S → P  
- \( r_p \): progenitor proliferation rate  
- \( \beta \): differentiation P → D  
- \( \delta_s, \delta_p, \delta_d \): death/turnover rates

---
## Features

- **Two biological scenarios**
  - Near-homeostasis  
  - Dysregulated stem/progenitor expansion  

- **High-quality visualizations**
  - Full dynamics for both scenarios  
  - Early-time zoom (first 10 time units)  
  - Seaborn-based smooth line plots  
  - Automatic PNG saving at 300 DPI  

- **Modular code**
  - Clean ODE system definition  
  - Parameter dictionaries for easy tuning  
  - Helper for figure saving  
  - Parameter-to-tuple conversion for `odeint`

---
## Requirements

- **Python 3.8+**

Install required libraries:
```bash
pip install numpy matplotlib seaborn scipy
```

---
## Project Files

- `stem_cell_dynamics_SPD.py` — Main simulation script  
- `figures/` — Automatically generated folder for output plots  

---
## Usage

Run the simulation with:

```bash
python stem_cell_dynamics_SPD.py
```

The script will:

1. Define the ODE system  
2. Simulate two biological scenarios  
3. Generate full and early-time plots  
4. Save output figures under `figures/`

---
## Output Files

The script generates the following figures:

| File | Description |
|------|-------------|
| `SPD-hom.png`        | Homeostatic scenario (full dynamics) |
| `SPD-dys-100.png`    | Dysregulated scenario (full dynamics) |
| `SPD-dys-10.png`     | Dysregulated scenario (first 10 time units) |

All images are saved in **300 DPI** with tight bounding boxes for clean formatting.

---
## Notes

- Initial conditions reflect a typical hierarchy: **S < P < D**.  
- Time units are arbitrary but internally consistent.  
- Logistic growth in the stem-cell compartment prevents unbounded expansion.  
- Dysregulated parameters mimic hyperproliferation or early tumor-like dynamics.  
- Numerical integration is performed using `scipy.integrate.odeint`.  

---
## License

Distributed under the **MIT License**.

---
## Author

Juliana Patiño Gallego  
jpatinoga@unal.edu.co
