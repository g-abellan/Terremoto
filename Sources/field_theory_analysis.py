import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os
import argparse
import math
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Haversine distance in km
def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0 # Earth radius in km
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2)**2 + 
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2)**2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

def time_integral(t_j, T_end, c, p):
    if T_end <= t_j:
        return 0.0
    if abs(p - 1.0) < 1e-5:
        return math.log((T_end - t_j + c) / c)
    else:
        return ((T_end - t_j + c)**(1.0 - p) - c**(1.0 - p)) / (1.0 - p)

def analyze_field_theory(csv_path, Mc=2.0):
    if not os.path.exists(csv_path):
        print(f"Error: File {csv_path} does not exist.")
        return
        
    df = pd.read_csv(csv_path)
    df = df.sort_values(by='time_utc').reset_index(drop=True)
    
    # Let's extract the two main doublet events
    # They are the two largest events in the catalog (M7.2 and M7.5 on June 24)
    # Let's find them
    mainshocks = df[df['magnitude'] >= 7.0].sort_values(by='time_utc').reset_index(drop=True)
    if len(mainshocks) < 2:
        print("Warning: Could not find both mainshocks (M >= 7.0) in the CSV. Using the first two events as seeds.")
        mainshocks = df.iloc[:2].copy()
        
    print("Detected doublet mainshocks:")
    for idx, row in mainshocks.iterrows():
        place_str = f" ({row['place']})" if 'place' in row else ""
        print(f"  Source {idx+1}: M {row['magnitude']} at {row['time_utc']}{place_str}")
        
    # Get times in days since the first mainshock
    t0_ms = mainshocks.loc[0, 'timestamp_ms']
    df['time_days'] = (df['timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
    
    t1 = 0.0 # first mainshock time is 0
    t2 = (mainshocks.loc[1, 'timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
    
    lat1, lon1 = mainshocks.loc[0, 'latitude'], mainshocks.loc[0, 'longitude']
    lat2, lon2 = mainshocks.loc[1, 'latitude'], mainshocks.loc[1, 'longitude']
    
    M1, M2 = mainshocks.loc[0, 'magnitude'], mainshocks.loc[1, 'magnitude']
    
    # We will exclude the mainshocks from the aftershock list when fitting the propagator
    aftershocks = df[~df['id'].isin(mainshocks['id'].values)].copy().reset_index(drop=True)
    n_after = len(aftershocks)
    print(f"Number of aftershocks for fitting: {n_after}")
    
    # Time window for fitting: from t1 to the end of the catalog
    T_start = 0.0
    T_end = df['time_days'].max()
    
    # Bounding box spatial area calculation (approximate in km^2)
    lat_min, lat_max = 9.5, 11.5
    lon_min, lon_max = -70.0, -65.0
    width_km = haversine_distance((lat_min + lat_max)/2.0, lon_min, (lat_min + lat_max)/2.0, lon_max)
    height_km = haversine_distance(lat_min, (lon_min + lon_max)/2.0, lat_max, (lon_min + lon_max)/2.0)
    area_km2 = width_km * height_km
    print(f"Seismic zone area: {area_km2:.1f} km^2")
    
    # Define the log-likelihood function for MLE
    # Parameters to fit: theta = [log10(mu), log10(K0), log10(c), p, log10(d), q, alpha]
    # We use log-transformed parameters to enforce positivity and improve convergence
    def neg_log_likelihood(params):
        mu = 10**params[0]
        K0 = 10**params[1]
        c = 10**params[2]
        p = params[3]
        d = 10**params[4]
        q = params[5]
        alpha_val = params[6]
        
        # Enforce physical bounds
        if p <= 1.0 or q <= 1.0 or p > 2.5 or q > 3.0:
            return 1e10
            
        # 1. Calculate rate lambda(t, x) for each aftershock
        total_log_lambda = 0.0
        for idx, row in aftershocks.iterrows():
            t = row['time_days']
            lat = row['latitude']
            lon = row['longitude']
            
            # Rate contribution from background
            rate = mu / area_km2
            
            # Contribution from Source 1
            if t > t1:
                r1 = haversine_distance(lat1, lon1, lat, lon)
                temp_factor1 = K0 * math.exp(alpha_val * (M1 - Mc))
                g_t1 = 1.0 / ((t - t1 + c)**p)
                g_s1 = 1.0 / ((r1**2 + d**2)**q)
                rate += temp_factor1 * g_t1 * g_s1
                
            # Contribution from Source 2
            if t > t2:
                r2 = haversine_distance(lat2, lon2, lat, lon)
                temp_factor2 = K0 * math.exp(alpha_val * (M2 - Mc))
                g_t2 = 1.0 / ((t - t2 + c)**p)
                g_s2 = 1.0 / ((r2**2 + d**2)**q)
                rate += temp_factor2 * g_t2 * g_s2
                
            if rate <= 1e-15:
                rate = 1e-15
            total_log_lambda += math.log(rate)
            
        # Time integral
        int_t1 = time_integral(t1, T_end, c, p)
        int_t2 = time_integral(t2, T_end, c, p)
        
        # Spatial integral over R^2: Int_R^2 (r^2 + d^2)^(-q) dx dy = pi * d^(2 - 2q) / (q - 1)
        int_s = math.pi * (d**(2.0 - 2.0 * q)) / (q - 1.0)
        
        integral_trigger1 = math.exp(alpha_val * (M1 - Mc)) * K0 * int_t1 * int_s
        integral_trigger2 = math.exp(alpha_val * (M2 - Mc)) * K0 * int_t2 * int_s if T_end > t2 else 0.0
        
        integral_bg = mu * (T_end - T_start)
        
        total_integral = integral_bg + integral_trigger1 + integral_trigger2
        
        return total_integral - total_log_lambda

    # Initial parameter guess
    # mu_init = 10 events/day, K0 = 0.01, c = 0.01 days, p = 1.1, d = 5 km, q = 1.5, alpha = 1.4
    init_params = [
        math.log10(10.0),   # log10(mu)
        math.log10(0.01),   # log10(K0)
        math.log10(0.01),   # log10(c)
        1.15,               # p
        math.log10(5.0),    # log10(d)
        1.5,                # q
        1.4                 # alpha
    ]
    
    print("Fitting field-theoretic parameters using MLE...")
    res = minimize(neg_log_likelihood, init_params, method='Nelder-Mead', options={'maxiter': 1000})
    
    if res.success:
        fit_params = res.x
        mu = 10**fit_params[0]
        K0 = 10**fit_params[1]
        c = 10**fit_params[2]
        p = fit_params[3]
        d = 10**fit_params[4]
        q = fit_params[5]
        alpha_val = fit_params[6]
        
        print("\n--- MLE Parameters Fit Results ---")
        print(f"Background rate (mu): {mu:.4f} events/day")
        print(f"Productivity constant (K0): {K0:.6f}")
        print(f"Omori time offset (c): {c:.4f} days ({c*24*60:.1f} minutes)")
        print(f"Omori decay exponent (p): {p:.4f}")
        print(f"Spatial scale (d): {d:.4f} km")
        print(f"Spatial decay exponent (q): {q:.4f}")
        print(f"Magnitude scaling (alpha): {alpha_val:.4f}")
    else:
        print("Error: MLE did not converge. Using initial guess for plotting.")
        mu, K0, c, p, d, q, alpha_val = 10.0, 0.01, 0.01, 1.15, 5.0, 1.5, 1.4
        
    # Save the parameters to a report file
    basename = os.path.basename(csv_path).replace('.csv', '')
    param_path = os.path.join(BASE_DIR, f"../Reports/{basename}_field_parameters.csv")
    pd.DataFrame({
        'parameter': ['mu', 'K0', 'c', 'p', 'd', 'q', 'alpha'],
        'value': [mu, K0, c, p, d, q, alpha_val],
        'description': [
            'Background rate (events/day)',
            'Productivity scaling factor',
            'Omori time offset (days)',
            'Omori time decay power',
            'Spatial scale (km)',
            'Spatial decay power',
            'Magnitude scaling parameter'
        ]
    }).to_csv(param_path, index=False)
    print(f"Parameters saved to: {param_path}")
    
    # 3. Generate 2D Spatial Susceptibility Map (Contour Plot)
    # Grid of latitudes and longitudes
    lat_grid = np.linspace(lat_min, lat_max, 100)
    lon_grid = np.linspace(lon_min, lon_max, 100)
    LON, LAT = np.meshgrid(lon_grid, lat_grid)
    
    # We will compute the accumulated field strength (integrated over time from T_start to T_end)
    # Field strength at (x, y) = Int_T_start^T_end lambda(t, x, y) dt
    # = mu/area_km2 * (T_end - T_start) + Sum_{j=1}^2 e^(alpha*(M_j-Mc)) * K0 * Int_t_j^T_end (t - t_j + c)^(-p) dt * (r_j^2 + d^2)^(-q)
    int_t1 = time_integral(t1, T_end, c, p)
    int_t2 = time_integral(t2, T_end, c, p)
    
    FIELD = np.zeros_like(LON)
    
    for i in range(LON.shape[0]):
        for j in range(LON.shape[1]):
            lat = LAT[i, j]
            lon = LON[i, j]
            
            # Background contribution (integrated over time)
            val = (mu / area_km2) * (T_end - T_start)
            
            # Source 1 contribution
            r1 = haversine_distance(lat1, lon1, lat, lon)
            val += math.exp(alpha_val * (M1 - Mc)) * K0 * int_t1 * (1.0 / ((r1**2 + d**2)**q))
            
            # Source 2 contribution
            if T_end > t2:
                r2 = haversine_distance(lat2, lon2, lat, lon)
                val += math.exp(alpha_val * (M2 - Mc)) * K0 * int_t2 * (1.0 / ((r2**2 + d**2)**q))
                
            FIELD[i, j] = val
            
    # Plotting
    plt.figure(figsize=(10, 8))
    # Plot field as filled contours (log scale for high dynamic range)
    contour = plt.contourf(LON, LAT, np.log10(FIELD), levels=20, cmap='inferno')
    cbar = plt.colorbar(contour)
    cbar.set_label('Log10 Cumulative Seismicity Density ($km^{-2}$)', fontsize=12)
    
    # Plot simplified coastline
    coast_lons = [-70.0, -69.2, -68.4, -68.3, -68.2, -68.0, -67.77, -67.6, -66.93, -66.07, -65.9, -65.4, -65.0]
    coast_lats = [11.3,  11.1,  10.9,  10.7,  10.5,  10.48, 10.5,   10.5,  10.6,   10.6,   10.3,  10.15, 10.10]
    plt.plot(coast_lons, coast_lats, color='lightgray', linewidth=2.0, linestyle='-', label='Coastline')
    
    # Plot aftershocks
    plt.scatter(aftershocks['longitude'], aftershocks['latitude'], 
                s=aftershocks['magnitude']**2 * 2, color='white', alpha=0.4, edgecolors='black', label='Aftershocks')
                
    # Plot mainshocks
    plt.scatter([lon1], [lat1], s=300, marker='*', color='cyan', edgecolors='black', linewidth=1.5, label='Mainshock 1 (M7.2)')
    plt.scatter([lon2], [lat2], s=400, marker='*', color='yellow', edgecolors='black', linewidth=1.5, label='Mainshock 2 (M7.5)')
    
    # Plot cities for reference
    cities = {
        'San Felipe': (10.48, -68.74),
        'Valencia': (10.16, -68.00),
        'Maracay': (10.25, -67.60),
        'Caracas': (10.48, -66.90),
        'La Guaira': (10.60, -66.93)
    }
    for city, coords in cities.items():
        plt.plot(coords[1], coords[0], 's', color='white', markersize=5, markeredgecolor='black')
        plt.text(coords[1] + 0.05, coords[0] + 0.03, city, color='white', fontsize=9, fontweight='bold')
        
    plt.xlabel('Longitude (°W)', fontsize=12)
    plt.ylabel('Latitude (°N)', fontsize=12)
    plt.title('2D Seismicity Susceptibility Field $\Phi(x, y)$', fontsize=14)
    plt.xlim(lon_min, lon_max)
    plt.ylim(lat_min, lat_max)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='upper right')
    
    plot_path = os.path.join(BASE_DIR, f"../Reports/{basename}_susceptibility_field.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Susceptibility field plot saved to: {plot_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze seismicity as a stochastic field.')
    parser.add_argument('--input', type=str, default=os.path.join(BASE_DIR, 'simulated_aftershocks.csv'),
                        help='Path to the seismic data CSV file')
    parser.add_argument('--mc', type=float, default=2.0, help='Cutoff magnitude Mc')
    
    args = parser.parse_args()
    
    # Make sure Reports directory exists
    os.makedirs(os.path.join(BASE_DIR, "../Reports"), exist_ok=True)
    
    analyze_field_theory(args.input, args.mc)
