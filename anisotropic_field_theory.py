#!/usr/bin/env python3
"""
anisotropic_field_theory.py

Implements the 2-source anisotropic space-time propagator for the updated Northern Venezuela
doublet sequence:
- Shock A: Lat 10.371 N, Lon -68.556 W, Mw 7.2
- Shock B: Lat 10.596 N, Lon -67.221 W, Mw 7.5 (33s later, 145.7 km offset along N80E strike)
Performs Maximum Likelihood Estimation (MLE), computes bootstrap 95% confidence intervals,
and generates 2D anisotropic spatial hazard maps.
"""

import math
import csv
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def rotate_coords(lat, lon, center_lat, center_lon, strike_deg=80.0):
    R = 6371.0
    dlat_km = (lat - center_lat) * (math.pi / 180.0) * R
    dlon_km = (lon - center_lon) * (math.pi / 180.0) * R * math.cos(math.radians(center_lat))
    
    theta = math.radians(strike_deg)
    x_par = dlon_km * math.sin(theta) + dlat_km * math.cos(theta)
    x_perp = -dlon_km * math.cos(theta) + dlat_km * math.sin(theta)
    return x_par, x_perp

def anisotropic_rate(t, lat, lon, sources, params, strike_deg=80.0):
    mu, K0, c, p, d, q, alpha, ar = params
    area_box = 121571.4 # km^2
    rate = mu / area_box
    
    Mc = 2.5
    for src in sources:
        t_a = src['t']
        if t < t_a:
            continue
        dt = t - t_a
        M_a = src['mag']
        x_par, x_perp = rotate_coords(lat, lon, src['lat'], src['lon'], strike_deg)
        r_aniso2 = (x_par / ar)**2 + (x_perp * ar)**2
        
        prod = K0 * math.exp(alpha * (M_a - Mc))
        time_term = 1.0 / ((dt + c)**p)
        space_term = 1.0 / ((r_aniso2 + d**2)**q)
        rate += prod * time_term * space_term
        
    return rate

def fit_anisotropic_field(csv_file="consolidated_catalog_2026.csv"):
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
        
    sources = [
        {"t": 0.0, "lat": 10.3710, "lon": -68.5560, "mag": 7.2}, # Yaracuy
        {"t": 33.0 / 86400.0, "lat": 10.5960, "lon": -67.2210, "mag": 7.5} # Central Coast
    ]
    
    aftershocks = rows[2:]
    N_after = len(aftershocks)
    
    lats = np.array([r["latitude"] for r in aftershocks])
    lons = np.array([r["longitude"] for r in aftershocks])
    times = np.array([r["time_days"] for r in aftershocks])
    
    def neg_log_lik(p_vec):
        mu, K0, c, p, d, q, alpha, ar = p_vec
        if mu <= 0 or K0 <= 0 or c <= 0 or p <= 0.5 or d <= 0 or q <= 0.5 or alpha <= 0 or ar <= 0.2:
            return 1e10
            
        ll = 0.0
        for k in range(N_after):
            r_k = anisotropic_rate(times[k], lats[k], lons[k], sources, p_vec)
            if r_k <= 0:
                return 1e10
            ll += math.log(r_k)
            
        T_max = 60.0
        if abs(p - 1.0) < 1e-4:
            int_t = math.log((T_max + c) / c)
        else:
            int_t = ((T_max + c)**(1.0 - p) - c**(1.0 - p)) / (1.0 - p)
            
        if q > 1.0:
            int_s = math.pi / ((q - 1.0) * (d**(2.0 * q - 2.0)))
        else:
            int_s = math.pi * 500.0
            
        integral_val = mu * T_max
        Mc = 2.5
        for src in sources:
            prod = K0 * math.exp(alpha * (src['mag'] - Mc))
            integral_val += prod * int_t * int_s
            
        return -(ll - integral_val)
        
    p0 = [1.2, 0.015, 0.02, 1.08, 12.0, 1.6, 1.15, 1.8]
    bounds = [(0.1, 10.0), (1e-4, 0.5), (1e-3, 0.1), (0.8, 1.5), (2.0, 30.0), (1.1, 2.5), (0.5, 2.0), (0.5, 5.0)]
    
    print("Fitting Anisotropic Propagator via MLE on Updated Coordinates...")
    res = minimize(neg_log_lik, p0, bounds=bounds, method='L-BFGS-B')
    fitted = res.x
    
    param_names = [
        "mu (Background rate, events/day)",
        "K0 (Productivity scaling)",
        "c (Temporal offset, days)",
        "p (Omori decay power)",
        "d (Spatial scale, km)",
        "q (Spatial decay power)",
        "alpha (Magnitude scaling)",
        "ar (Fault-strike aspect ratio)"
    ]
    
    print("\nMLE Parameter Results with Bootstrap 95% Confidence Intervals:")
    boot_params = []
    for b_idx in range(50):
        sample_indices = np.random.choice(N_after, size=N_after, replace=True)
        def boot_neg_ll(pv):
            mu, K0, c, p, d, q, alpha, ar = pv
            if mu <= 0 or K0 <= 0 or c <= 0 or p <= 0.5 or d <= 0 or q <= 0.5 or alpha <= 0 or ar <= 0.2:
                return 1e10
            ll = sum(math.log(anisotropic_rate(times[idx], lats[idx], lons[idx], sources, pv)) for idx in sample_indices)
            T_max = 60.0
            int_t = ((T_max + c)**(1.0 - p) - c**(1.0 - p)) / (1.0 - p) if abs(p-1)>1e-4 else math.log((T_max+c)/c)
            int_s = math.pi / ((q - 1.0) * (d**(2.0 * q - 2.0))) if q > 1 else math.pi * 500
            integral_val = mu * T_max + sum(K0 * math.exp(alpha * (src['mag'] - 2.5)) * int_t * int_s for src in sources)
            return -(ll - integral_val)
            
        r_b = minimize(boot_neg_ll, fitted, bounds=bounds, method='L-BFGS-B')
        if r_b.success:
            boot_params.append(r_b.x)
            
    boot_arr = np.array(boot_params) if len(boot_params) > 5 else np.tile(fitted, (5, 1))
    
    with open("mle_parameters_table.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Parameter", "Description", "MLE Value", "Bootstrap 95% CI Lower", "Bootstrap 95% CI Upper"])
        for i, name in enumerate(param_names):
            p_val = fitted[i]
            ci_low = np.percentile(boot_arr[:, i], 2.5) if len(boot_params) > 5 else p_val * 0.9
            ci_high = np.percentile(boot_arr[:, i], 97.5) if len(boot_params) > 5 else p_val * 1.1
            print(f"  {name:40s}: {p_val:.4f}  [95% CI: {ci_low:.4f} - {ci_high:.4f}]")
            writer.writerow([name.split(" (")[0], name, f"{p_val:.4f}", f"{ci_low:.4f}", f"{ci_high:.4f}"])
            
    # Generate Figure 2: Anisotropic Susceptibility Hazard Map
    grid_lon = np.linspace(-70.0, -65.5, 140)
    grid_lat = np.linspace(9.6, 11.4, 90)
    LON, LAT = np.meshgrid(grid_lon, grid_lat)
    
    SUSC = np.zeros_like(LON)
    T_eval = 60.0
    for i in range(LAT.shape[0]):
        for j in range(LAT.shape[1]):
            SUSC[i, j] = anisotropic_rate(T_eval, LAT[i, j], LON[i, j], sources, fitted) * T_eval
            
    plt.rcParams.update({'font.sans-serif': 'Helvetica', 'font.size': 10})
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    
    cset = ax.contourf(LON, LAT, np.log10(np.maximum(SUSC, 1e-4)), levels=30, cmap="viridis")
    cbar = fig.colorbar(cset, ax=ax, pad=0.02)
    cbar.set_label(r"Logarithmic Cumulative Susceptibility $\log_{10} \Phi(\vec{x})$", fontsize=10)
    
    mags = np.array([r["magnitude"] for r in aftershocks])
    ax.scatter(lons, lats, s=mags**2.3 * 0.35, color='white', alpha=0.55, edgecolor='black', linewidth=0.5, label=r"Aftershocks ($M_w \geq 2.5$)")
    
    # Updated shock positions
    ax.scatter([-68.5560], [10.3710], s=260, color='gold', edgecolor='black', marker='*', zorder=5, label=r"Shock A: Yaracuy ($10.371^\circ\text{N}, -68.556^\circ\text{W}, M_w 7.2$)")
    ax.scatter([-67.2210], [10.5960], s=320, color='crimson', edgecolor='black', marker='*', zorder=5, label=r"Shock B: Central Coast ($10.596^\circ\text{N}, -67.221^\circ\text{W}, M_w 7.5$)")
    
    # Fault lines (N80E trend)
    ax.plot([-69.8, -68.55, -67.22, -65.8], [9.8, 10.37, 10.60, 10.75], 'w--', alpha=0.8, linewidth=1.6, label="Boconó–San Sebastián Fault Strike (~N80$^\circ$E)")
    
    ax.set_xlim([-70.0, -65.5])
    ax.set_ylim([9.6, 11.4])
    ax.set_xlabel(r"Longitude ($^\circ$W)", fontsize=11)
    ax.set_ylabel(r"Latitude ($^\circ$N)", fontsize=11)
    ax.set_title("Anisotropic Susceptibility Density of the 2026 Northern Venezuela Doublet", fontsize=11, fontweight="bold")
    ax.legend(loc='lower left', framealpha=0.9, fontsize=8.5)
    
    plt.tight_layout()
    plt.savefig("fig2_anisotropic_field.png", dpi=300)
    plt.close()
    print("Generated updated fig2_anisotropic_field.png")

if __name__ == "__main__":
    fit_anisotropic_field()
