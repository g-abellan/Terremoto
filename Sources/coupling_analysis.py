import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import pearsonr, spearmanr
import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Haversine distance in km
def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0 # Earth radius in km
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2)**2 + 
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2)**2)
    c_val = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c_val

# Time integral helper
def time_integral(t_j, T_end, c, p):
    if T_end <= t_j:
        return 0.0
    if abs(p - 1.0) < 1e-5:
        return math.log((T_end - t_j + c) / c)
    else:
        return ((T_end - t_j + c)**(1.0 - p) - c**(1.0 - p)) / (1.0 - p)

def analyze_coupling():
    # Paths
    annotated_network_path = os.path.join(BASE_DIR, "../Reports/simulated_aftershocks_annotated_network.csv")
    field_params_path = os.path.join(BASE_DIR, "../Reports/simulated_aftershocks_field_parameters.csv")
    
    if not os.path.exists(annotated_network_path) or not os.path.exists(field_params_path):
        print("Required files do not exist.")
        return
        
    df = pd.read_csv(annotated_network_path)
    df_params = pd.read_csv(field_params_path)
    
    # Map parameters
    params = dict(zip(df_params['parameter'], df_params['value']))
    mu = params['mu']
    K0 = params['K0']
    c = params['c']
    p = params['p']
    d = params['d']
    q = params['q']
    alpha_val = params['alpha']
    Mc = 2.0
    
    # Find doublet mainshocks
    mainshocks = df[df['magnitude'] >= 7.0].sort_values(by='time_utc').reset_index(drop=True)
    if len(mainshocks) < 2:
        mainshocks = df.iloc[:2].copy()
        
    t0_ms = mainshocks.loc[0, 'timestamp_ms']
    df['time_days'] = (df['timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
    
    t1 = 0.0
    t2 = (mainshocks.loc[1, 'timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
    
    lat1, lon1 = mainshocks.loc[0, 'latitude'], mainshocks.loc[0, 'longitude']
    lat2, lon2 = mainshocks.loc[1, 'latitude'], mainshocks.loc[1, 'longitude']
    
    M1, M2 = mainshocks.loc[0, 'magnitude'], mainshocks.loc[1, 'magnitude']
    
    # Calculate area
    lat_min, lat_max = 9.5, 11.5
    lon_min, lon_max = -70.0, -65.0
    width_km = haversine_distance((lat_min + lat_max)/2.0, lon_min, (lat_min + lat_max)/2.0, lon_max)
    height_km = haversine_distance(lat_min, (lon_min + lon_max)/2.0, lat_max, (lon_min + lon_max)/2.0)
    area_km2 = width_km * height_km
    
    T_end = df['time_days'].max()
    
    # Compute continuous rate, susceptibility, and exact Laplacian
    rates = []
    susceptibilities = []
    laplacians = []
    
    # Pre-calculate integrals
    I1 = time_integral(t1, T_end, c, p)
    I2 = time_integral(t2, T_end, c, p)
    
    K1_coeff = math.exp(alpha_val * (M1 - Mc)) * K0 * I1
    K2_coeff = math.exp(alpha_val * (M2 - Mc)) * K0 * I2 if T_end > t2 else 0.0
    
    for idx, row in df.iterrows():
        t = row['time_days']
        lat = row['latitude']
        lon = row['longitude']
        
        # 1. Continuous field rate
        rate = mu / area_km2
        if t > t1:
            r1 = haversine_distance(lat1, lon1, lat, lon)
            rate += K0 * math.exp(alpha_val * (M1 - Mc)) * (1.0 / ((t - t1 + c)**p)) * (1.0 / ((r1**2 + d**2)**q))
        if t > t2:
            r2 = haversine_distance(lat2, lon2, lat, lon)
            rate += K0 * math.exp(alpha_val * (M2 - Mc)) * (1.0 / ((t - t2 + c)**p)) * (1.0 / ((r2**2 + d**2)**q))
        rates.append(rate)
        
        # 2. Cumulative susceptibility field
        phi = (mu / area_km2) * T_end
        r1 = haversine_distance(lat1, lon1, lat, lon)
        phi += math.exp(alpha_val * (M1 - Mc)) * K0 * I1 * (1.0 / ((r1**2 + d**2)**q))
        if T_end > t2:
            r2 = haversine_distance(lat2, lon2, lat, lon)
            phi += math.exp(alpha_val * (M2 - Mc)) * K0 * I2 * (1.0 / ((r2**2 + d**2)**q))
        susceptibilities.append(phi)
        
        # 3. Exact analytical spatial Laplacian
        lap1 = 4.0 * q * (q * r1**2 - d**2) / ((r1**2 + d**2)**(q + 2.0))
        lap2 = 4.0 * q * (q * r2**2 - d**2) / ((r2**2 + d**2)**(q + 2.0)) if T_end > t2 else 0.0
        laplacians.append(K1_coeff * lap1 + K2_coeff * lap2)
        
    df['field_rate'] = rates
    df['susceptibility'] = susceptibilities
    df['laplacian'] = laplacians
    
    # Exclude mainshocks to focus on aftershocks
    df_after = df[~df['id'].isin(mainshocks['id'].values)].copy()
    
    # Calculate correlations (Rate vs FRC)
    p_r, p_p = pearsonr(df_after['forman_curvature'], df_after['field_rate'])
    s_r, s_p = spearmanr(df_after['forman_curvature'], df_after['field_rate'])
    
    # Calculate correlations (Laplacian vs FRC)
    lap_p_r, lap_p_p = pearsonr(df_after['forman_curvature'], df_after['laplacian'])
    lap_s_r, lap_s_p = spearmanr(df_after['forman_curvature'], df_after['laplacian'])
    
    print(f"Coupling Correlations calculated:")
    print(f"  Rate vs FRC: Pearson r = {p_r:.4f} (p = {p_p:.2e}), Spearman rho = {s_r:.4f} (p = {s_p:.2e})")
    print(f"  Laplacian vs FRC: Pearson r = {lap_p_r:.4f} (p = {lap_p_p:.2e}), Spearman rho = {lap_s_r:.4f} (p = {lap_s_p:.2e})")
    
    # Generate publication-quality plot (2 panels side-by-side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Use serif font to match LaTeX style
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['text.usetex'] = False
    
    # ------------------ PANEL A: FRC vs Log-Rate ------------------
    x_a = np.log10(df_after['field_rate'])
    y_a = df_after['forman_curvature']
    
    scatter_a = ax1.scatter(x_a, y_a, c=df_after['magnitude'], cmap='plasma', alpha=0.15, edgecolors='none', s=df_after['magnitude']**2 * 4)
    
    # Compute bin averages for Panel A
    num_bins = 12
    bin_edges = np.linspace(min(x_a), max(x_a), num_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    binned_means_a = []
    binned_sems_a = []
    for i in range(num_bins):
        mask = (x_a >= bin_edges[i]) & (x_a < bin_edges[i+1])
        if i == num_bins - 1:
            mask = mask | (x_a == bin_edges[i+1])
        bin_y = y_a[mask]
        if len(bin_y) > 0:
            binned_means_a.append(np.mean(bin_y))
            binned_sems_a.append(np.std(bin_y) / np.sqrt(len(bin_y)) if len(bin_y) > 1 else 0.0)
        else:
            binned_means_a.append(np.nan)
            binned_sems_a.append(0.0)
            
    valid_a = ~np.isnan(binned_means_a)
    ax1.errorbar(bin_centers[valid_a], np.array(binned_means_a)[valid_a], yerr=np.array(binned_sems_a)[valid_a], 
                 fmt='o', color='black', ecolor='red', elinewidth=2, capsize=4, markersize=8, label='Bin-Averaged Curvature $\\langle F \\rangle$')
    
    m_a, b_a = np.polyfit(x_a, y_a, 1)
    x_fit_a = np.linspace(min(x_a), max(x_a), 100)
    ax1.plot(x_fit_a, m_a * x_fit_a + b_a, 'r--', linewidth=2, label=f'Linear fit (Slope: {m_a:.2f})')
    
    ax1.set_xlabel('$\\log_{10}$ Continuous Field Rate $\\lambda(t, \\vec{x})$ ($km^{-2} \\cdot day^{-1}$)', fontsize=12)
    ax1.set_ylabel('Forman-Ricci Curvature $F(u)$', fontsize=12)
    ax1.set_title('(a) FRC vs Logarithmic Field Rate (Proxy)', fontsize=13, fontweight='bold', pad=10)
    ax1.grid(True, linestyle=':', alpha=0.5)
    
    stats_a = f"Spearman $\\rho = {s_r:.4f}$ ($p = {s_p:.2e}$)\nPearson $r = {p_r:.4f}$ ($p = {p_p:.2e}$)"
    ax1.text(0.05, 0.05, stats_a, transform=ax1.transAxes, fontsize=10,
             verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    ax1.legend(loc='upper right')
    
    # ------------------ PANEL B: FRC vs Laplacian ------------------
    x_b = df_after['laplacian']
    y_b = df_after['forman_curvature']
    
    scatter_b = ax2.scatter(x_b, y_b, c=df_after['magnitude'], cmap='plasma', alpha=0.15, edgecolors='none', s=df_after['magnitude']**2 * 4)
    
    # Compute quantile bin averages for Panel B (to deal with highly skewed Laplacian values)
    num_bins_b = 8
    df_after['lap_bin'] = pd.qcut(df_after['laplacian'], q=num_bins_b, labels=False, duplicates='drop')
    binned_b = df_after.groupby('lap_bin').agg(
        mean_lap=('laplacian', 'mean'),
        mean_frc=('forman_curvature', 'mean'),
        std_frc=('forman_curvature', 'std'),
        count=('id', 'count')
    ).reset_index()
    
    binned_b['sem_frc'] = binned_b['std_frc'] / np.sqrt(binned_b['count'])
    
    ax2.errorbar(binned_b['mean_lap'], binned_b['mean_frc'], yerr=binned_b['sem_frc'], 
                 fmt='s', color='black', ecolor='blue', elinewidth=2, capsize=4, markersize=8, label='Quantile-Averaged Curvature $\\langle F \\rangle$')
    
    m_b, b_b = np.polyfit(x_b, y_b, 1)
    x_fit_b = np.linspace(min(x_b), max(x_b), 100)
    ax2.plot(x_fit_b, m_b * x_fit_b + b_b, 'b--', linewidth=2, label=f'Linear fit (Slope: {m_b:.2e})')
    
    ax2.set_xlabel('Analytical Susceptibility Laplacian $\\nabla^2 \\Phi(\\vec{x})$ ($km^{-2}$)', fontsize=12)
    ax2.set_ylabel('Forman-Ricci Curvature $F(u)$', fontsize=12)
    ax2.set_title('(b) FRC vs Spatial Susceptibility Laplacian (Direct)', fontsize=13, fontweight='bold', pad=10)
    ax2.grid(True, linestyle=':', alpha=0.5)
    
    stats_b = f"Spearman $\\rho = {lap_s_r:.4f}$ ($p = {lap_s_p:.2f}$)\nPearson $r = {lap_p_r:.4f}$ ($p = {lap_p_p:.2e}$)"
    ax2.text(0.05, 0.05, stats_b, transform=ax2.transAxes, fontsize=10,
             verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))
    ax2.legend(loc='upper right')
    
    # Combined colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(scatter_b, cax=cbar_ax)
    cbar.set_label('Aftershock Magnitude ($M_w$)', fontsize=12)
    
    fig.suptitle('Continuous-Discrete Coupling: Discrete Network Curvature vs. Continuous Field Representation', fontsize=15, fontweight='bold', y=0.98)
    
    plot_path = os.path.join(BASE_DIR, "../Reports/simulated_aftershocks_coupling.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Coupling plot saved to: {plot_path}")
    
    # Save correlation statistics
    stats_df = pd.DataFrame({
        'correlation_pair': ['curvature_vs_field_rate', 'curvature_vs_laplacian'],
        'pearson_r': [p_r, lap_p_r],
        'pearson_p': [p_p, lap_p_p],
        'spearman_rho': [s_r, lap_s_r],
        'spearman_p': [s_p, lap_s_p]
    })
    stats_df.to_csv(os.path.join(BASE_DIR, "../Reports/simulated_aftershocks_correlations.csv"), index=False)
    print("Updated correlations CSV.")

if __name__ == "__main__":
    analyze_coupling()
