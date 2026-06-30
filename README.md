# A Hybrid Information-Geometric and Stochastic Field Theory Framework for the 2026 Caracas Doublet Earthquake Sequence

This repository contains the replication code, dataset, and LaTeX manuscript draft for the article *"A Hybrid Information-Geometric and Stochastic Field Theory Framework for the 2026 Caracas Doublet Earthquake Sequence"*. 

The framework bridges discrete information geometry (via causal network topology) and continuous stochastic field theory (via maximum-likelihood field propagators) to model, analyze, and forecast the relaxation sequence of the devastating doublet earthquake of June 24, 2026, in northern Venezuela.

---

## Directory Structure

*   `Sources/`: Operational Python scripts that execute the calculations, model fits, and plotting.
*   `Drafts/`: LaTeX source files, citations (`references.bib`), and style files used to compile the manuscript.
*   `Reports/`: Outputs of the pipeline, including parameter tables, statistical reports, and publication-quality figures.
*   `requirements.txt`: List of Python library dependencies required to run the scripts.

---

## Python Scripts (in `Sources/`)

1.  `run_pipeline.py`: The master execution script. It runs the entire database building, simulation, estimation, correlation, and forecasting pipeline in sequence.
2.  `fetch_data.py`: Connects to the USGS ComCat API to download the seed earthquake catalog ($M_w \ge 4.3$) within the regional bounding box for the plate boundary corridor (June 24 to June 29, 2026).
3.  `etas_simulation.py`: Implements a branching Epidemic-Type Aftershock Sequence (ETAS) simulation to generate a realistic local micro-seismicity catalog down to $M_c = 2.0$ using the USGS events as seeds.
4.  `geometry_analysis.py`: Constructs the causal network using the Baiesi-Paczuski metric distance, extracts the directed triggering tree, and computes the Forman-Ricci Curvature (FRC) of nodes as a function of time (**generates Figure 1**).
5.  `field_theory_analysis.py`: Computes the continuous susceptibility field, estimates the stochastic field propagator parameters using Maximum Likelihood Estimation (MLE), and plots the 2D cumulative hazard field map (**generates Figure 2**).
6.  `coupling_analysis.py`: Calculates Pearson and Spearman correlations between network FRC and the continuous field rate (log-rate proxy) as well as the exact analytical spatial Laplacian $\nabla^2 \Phi(\vec{x})$, plotting both relationships in a publication-quality two-panel format (**generates Figure 3**).
7.  `run_forecast.py`: Computes the forward prospective 30-day forecast of seismicity rates (regional total and Caracas metropolitan sub-basin) across different magnitude thresholds.
8.  `summarize_data.py`: A helper script that parses the seed catalog and prints a summary of events, locations, and released seismic energy.

---

## Installation & Setup

All scripts are written in standard Python 3. The pipeline requires the standard scientific stack (`numpy`, `pandas`, `scipy`, `matplotlib`, `networkx`).

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/<username>/Terremoto.git
    cd Terremoto
    ```

2.  **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

---

## Running the Pipeline

To run all analysis steps, fit the propagator parameters, and regenerate all figures:
```bash
python3 Sources/run_pipeline.py
```
Upon successful execution, the following outputs will be updated in the `Reports/` directory:
*   `simulated_aftershocks_geometry_evolution.png` (Figure 1: curvature relaxation curve).
*   `simulated_aftershocks_susceptibility_field.png` (Figure 2: 2D spatial susceptibility map).
*   `simulated_aftershocks_coupling.png` (Figure 3: two-panel continuous-discrete coupling plot).
*   `simulated_aftershocks_field_parameters.csv` (Table I: MLE parameters fit).
*   `simulated_aftershocks_correlations.csv` (Statistical summary of coupling correlations).

