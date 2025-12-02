# Project 3 – Biostatistics and Computational Models

This project consists of three components, each focused on a different type of quantitative analysis applied to biological data.

# 🗂 Project Structure

```txt
project3a_statistical_analysis_replication/      → KM survival analysis (GSE14520)
project3b_curve_fitting/                         → Weibull survival model under noise
project3c_dynamic_model/                         → Stem cell ODE system (S–P–D model)
README.md                                        → General overview
```

---
# 1. Kaplan–Meier Survival Analysis
This module includes:

- Loading and preprocessing GEO data (GSE14520)  
- Mapping probes → gene symbols and building the 13-gene TEX signature  
- Computing a **risk score** using Cox coefficients  
- Kaplan–Meier curves for High vs Low risk groups  
- Log-rank test and number-at-risk table

---
# 2. Weibull Parametric Survival Model

This module includes:

- Generation of a “true” survival curve using a 3-parameter Weibull model  
- Simulation of Gaussian noise (σ = 0.05–0.50)  
- Comparison of three fitting strategies:
  - With initial guesses  
  - Without initial guesses (bounded)  
  - Without bounds and no initial guesses  
- Survival and residual plots for each noise level  
- Summary table of fitted parameters and \(R^2\)

---
# 3. Stem Cell Dynamics (ODE S–P–D Model)

This module simulates a three-compartment ordinary differential equation model:

- **S**: stem cells  
- **P**: progenitor cells  
- **D**: differentiated cells  

Two biological scenarios are included:

1. **Homeostasis** — controlled growth and stable differentiation  
2. **Dysregulation** — excessive proliferation and reduced differentiation  

The script generates full time-course plots and an early-time zoom.

---
## License

Distributed under the **MIT License**.

---
## Author

Juliana Patiño Gallego  
jpatinoga@unal.edu.co
