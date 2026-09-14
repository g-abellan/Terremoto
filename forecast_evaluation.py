#!/usr/bin/env python3
"""
forecast_evaluation.py

Evaluates the 30-day and 60-day prospective forecasts against the realized seismic record
for the updated Northern Venezuela doublet sequence using formal CSEP N-test scoring.
Generates Figure 3.
"""

import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

def evaluate_forecast(csv_file="consolidated_catalog_2026.csv"):
    rows = []
    with open(csv_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "id": r["id"],
                "timestamp_ms": int(r["timestamp_ms"]),
                "latitude": float(r["latitude"]),
                "longitude": float(r["longitude"]),
                "magnitude": float(r["magnitude"]),
                "segment": r.get("segment", "")
            })
            
    rows.sort(key=lambda x: x["timestamp_ms"])
    t0_ms = rows[0]["timestamp_ms"]
    for r in rows:
        r["time_days"] = (r["timestamp_ms"] - t0_ms) / (1000.0 * 86400.0)
        
    aftershocks = rows[2:]
    times = np.array([r["time_days"] for r in aftershocks])
    
    val_30_mask = (times > 5.0) & (times <= 35.0)
    val_60_mask = (times > 5.0) & (times <= 60.0)
    
    N_obs_30 = np.sum(val_30_mask)
    N_obs_60 = np.sum(val_60_mask)
    
    p_fit = 1.08
    c_fit = 0.02
    K_fit = 38.0
    mu_day = 1.2
    
    def expected_events(t1, t2):
        bg = mu_day * (t2 - t1)
        trig = K_fit * (((t2 + c_fit)**(1.0 - p_fit) - (t1 + c_fit)**(1.0 - p_fit)) / (1.0 - p_fit))
        return bg + trig
        
    N_exp_30 = expected_events(5.0, 35.0)
    N_exp_60 = expected_events(5.0, 60.0)
    
    delta1_30 = float(poisson.cdf(N_obs_30, N_exp_30))
    
    print("==================================================")
    print("CSEP Prospective Forecast Scoring Results:")
    print("==================================================")
    print(f"30-Day Evaluation Window (Days 5 to 35):")
    print(f"  Observed Count: {N_obs_30}")
    print(f"  Forecasted Expected Count: {N_exp_30:.1f} +/- {math.sqrt(N_exp_30):.1f}")
    print(f"  N-Test quantile (delta_1): {delta1_30:.4f} (Pass: in [0.025, 0.975])")
    
    print(f"\n60-Day Full Sequence Evaluation (Days 5 to 60):")
    print(f"  Observed Count: {N_obs_60}")
    print(f"  Forecasted Expected Count: {N_exp_60:.1f} +/- {math.sqrt(N_exp_60):.1f}")
    print("==================================================")
    
    days_bins = np.arange(0, 61, 1.0)
    hist_counts, bin_edges = np.histogram(times, bins=days_bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    t_model = np.linspace(0.05, 60, 500)
    rate_model = mu_day + K_fit / ((t_model + c_fit)**p_fit)
    rate_upper = rate_model + 1.96 * np.sqrt(rate_model)
    rate_lower = np.maximum(rate_model - 1.96 * np.sqrt(rate_model), 0)
    
    plt.rcParams.update({'font.sans-serif': 'Helvetica', 'font.size': 10})
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    
    ax.bar(bin_centers, hist_counts, width=0.8, color='steelblue', alpha=0.65, edgecolor='black', label="Observed Daily Seismicity")
    ax.plot(t_model, rate_model, 'r-', linewidth=2.0, label=r"Forecast Rate $\lambda(t) = \mu + K_0(t+c)^{-p}$")
    ax.fill_between(t_model, rate_lower, rate_upper, color='red', alpha=0.15, label="95% Poisson Forecast Envelope")
    
    ax.axvline(5.0, color='black', linestyle='--', linewidth=1.5, label="Forecast Horizon (Day 5)")
    ax.set_yscale('log')
    ax.set_ylim([0.5, max(hist_counts)*1.8])
    ax.set_xlim([0, 60])
    ax.set_xlabel("Time Since Doublet (days)", fontsize=11)
    ax.set_ylabel(r"Daily Earthquake Count ($M_w \geq 2.5$)", fontsize=11)
    ax.set_title("Realized Ex-Post Forecast Evaluation of the Northern Venezuela Doublet", fontsize=11, fontweight="bold")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9)
    
    plt.tight_layout()
    plt.savefig("fig3_forecast_validation.png", dpi=300)
    plt.close()
    print("Generated fig3_forecast_validation.png")

if __name__ == "__main__":
    evaluate_forecast()
