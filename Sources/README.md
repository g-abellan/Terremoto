# Replication Code: 2026 Caracas Doublet Earthquake Seismology Paper

This folder contains the complete Python pipeline to replicate the results, parameter estimations, and figures presented in the article *"A Hybrid Information-Geometric and Stochastic Field Theory Framework for the 2026 Caracas Doublet Earthquake Sequence"*.

## Pipeline Directory Structure

*   `fetch_data.py`: Downloads the seed catalog ($M_w \ge 4.3$) for the active plate boundary corridor from the USGS ComCat API between June 24 and June 29, 2026.
*   `etas_simulation.py`: Performs a branching Epidemic-Type Aftershock Sequence (ETAS) simulation from the seeds to generate the micro-seismicity catalog ($M_c = 2.0$).
*   `geometry_analysis.py`: Constructs the information-geometric causal network using the Baiesi-Paczuski metric, extracts the directed spanning tree (stress geodesics), computes the Forman-Ricci curvature (FRC) for all events, and saves the curvature relaxation metrics (**Figure 1**).
*   `field_theory_analysis.py`: Fits the stochastic field propagator parameters using Maximum Likelihood Estimation (MLE) and generates the 2D spatial susceptibility field map with coastline and city overlays (**Figure 2**).
*   `coupling_analysis.py`: Computes the Pearson and Spearman correlations between network curvature and field rate, performs the linear regression fit, and generates the bin-averaged scaling coupling plot (**Figure 3**).
*   `run_forecast.py`: Computes the forward prospective 30-day forecast (total events and Caracas metropolitan area) as a function of the magnitude threshold $M_{\text{th}}$.
*   `run_pipeline.py`: Master script that runs all steps in sequence.

## Dependencies

The code requires standard Python science libraries:
*   `numpy`
*   `pandas`
*   `matplotlib`
*   `scipy`

Install dependencies using `pip`:
```bash
pip install -r requirements.txt
```

## Running the Replication Pipeline

To run the entire pipeline and regenerate all datasets and figures in one click:
```bash
python Sources/run_pipeline.py
```

All generated figures and fitted parameters will be saved directly into the `Reports/` directory.

## Compiling the Manuscript

The LaTeX sources are in the `Drafts/` directory. To compile the PDF:
```bash
pdflatex Drafts/main.tex
bibtex Drafts/main
pdflatex Drafts/main.tex
pdflatex Drafts/main.tex
```
This compiles the final manuscript as `Drafts/main.pdf`.
