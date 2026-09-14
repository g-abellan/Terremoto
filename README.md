# 2026 Northern Venezuela Doublet Earthquake Sequence: Data & Code

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2609.XXXXX-b31b1b.svg)](https://arxiv.org)

This repository contains the consolidated seismic catalog, network topology algorithms, anisotropic field propagator routines, forecast evaluation benchmarks, and the research paper:

> **Segmented Rupture Dynamics, Causal Network Topology, and Anisotropic Stress Transfer of the 2026 Northern Venezuela Doublet Earthquake Sequence**  
> *G. Abellán* (Escuela de Física, Facultad de Ciencias, Universidad Central de Venezuela)

---

## 📁 Repository Contents

| File | Description |
| :--- | :--- |
| `main_v03.pdf` | Research paper (6-page preprint, RevTeX 4-1 format). |
| `consolidated_catalog_2026.csv` | Relocated 60-day seismic catalog with doublet hypocenters ($M_w \ge 2.5$). |
| `generate_consolidated_catalog.py` | Pipeline for fetching/assembling regional broadband and teleseismic data. |
| `network_topology.py` | Computes Baiesi--Paczuski causal tree, exact Forman--Ricci curvature ($F$), and weighted curvature ($F_W$). |
| `anisotropic_field_theory.py` | Calibrates two-source anisotropic space-time propagator (MLE + bootstrap CIs). |
| `forecast_evaluation.py` | Prospective 30-day and 60-day forecast testing against CSEP Poisson $N$-test. |
| `fig1_network_topology.png` | Figure 1: Network curvature evolution and degree distribution. |
| `fig2_anisotropic_field.png` | Figure 2: Anisotropic cumulative susceptibility field along $\text{N}80^\circ\text{E}$ strike. |
| `fig3_forecast_validation.png` | Figure 3: Daily aftershock rate and 60-day forecast validation. |
| `mle_parameters_table.csv` | Maximum-likelihood point-process parameters with 95% bootstrap confidence intervals. |
| `network_metrics_summary.csv` | Event-by-event table of causal parents, metric distances, and Forman curvatures. |
| `requirements.txt` | Python dependencies. |
| `.gitignore` | Standard exclusion list for caches and environments. |

---

## ⚙️ Quick Start & Reproduction

### 1. Installation
Install dependencies via `requirements.txt`:
```bash
pip install -r requirements.txt
```

### 2. Network Topology Analysis (Generates Figure 1)
```bash
python network_topology.py
```
*Outputs:* `fig1_network_topology.png`, `network_metrics_summary.csv`.

### 3. Calibrate Anisotropic Propagator (Generates Figure 2 & Table I)
```bash
python anisotropic_field_theory.py
```
*Outputs:* `fig2_anisotropic_field.png`, `mle_parameters_table.csv`.

### 4. Prospective Forecast Validation (Generates Figure 3)
```bash
python forecast_evaluation.py
```
*Outputs:* `fig3_forecast_validation.png` and CSEP $N$-test statistics.

---

## 🔬 Key Doublet Sequence Parameters
- **Shock A (Western Segment):** Lat $10.371^\circ\text{N}$, Lon $-68.556^\circ\text{W}$, depth $10.0\text{ km}$, $M_w 7.2$, Origin: 22:04:31 UTC (Yaracuy corridor).
- **Shock B (Eastern Segment):** Lat $10.596^\circ\text{N}$, Lon $-67.221^\circ\text{W}$, depth $10.0\text{ km}$, $M_w 7.5$, Origin: 22:05:04 UTC (San Sebastián offshore fault).
- **Spatial Separation:** $\Delta L \approx 145.7\text{ km}$ along $\text{N}80^\circ\text{E}$ strike.
- **Temporal Lag:** $\Delta t = 33\text{ s} \implies v_{\text{trig}} \approx 4.41\text{ km/s}$ ($S$/Rayleigh wave dynamic rupture triggering).

---

## 📜 Citation & License
Open access under the MIT License. If using this code or dataset, please cite the corresponding paper.
