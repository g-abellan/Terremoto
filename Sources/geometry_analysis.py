import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
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

def analyze_geometry(csv_path, b_val=1.0, df_val=2.0, Mc=2.0, window_size=100, step_size=10):
    if not os.path.exists(csv_path):
        print(f"Error: File {csv_path} does not exist.")
        return
        
    df = pd.read_csv(csv_path)
    df = df.sort_values(by='time_utc').reset_index(drop=True)
    n_events = len(df)
    print(f"Loaded {n_events} events for geometric analysis...")
    
    # Calculate time in days since the first event
    t0_ms = df['timestamp_ms'].min()
    df['time_days'] = (df['timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
    
    # 1. Construct Baiesi-Paczuski Causal Network
    # For each event j > 0, find parent i < j that minimizes eta_ij
    # eta_ij = (t_j - t_i) * 10^(-b*(M_i - Mc)) * (L_ij)^df
    # We set a small epsilon spatial distance to avoid 0 km issues
    eps_dist = 0.1 # 100 meters
    
    parents = []
    parent_distances = []
    
    print("Constructing Baiesi-Paczuski causal network...")
    # Seed event (first event) has no parent
    parents.append(None)
    parent_distances.append(np.nan)
    
    for j in range(1, n_events):
        t_j = df.loc[j, 'time_days']
        lat_j = df.loc[j, 'latitude']
        lon_j = df.loc[j, 'longitude']
        
        min_eta = float('inf')
        best_parent = None
        
        # Look at all previous events i < j
        for i in range(j):
            t_i = df.loc[i, 'time_days']
            lat_i = df.loc[i, 'latitude']
            lon_i = df.loc[i, 'longitude']
            M_i = df.loc[i, 'magnitude']
            
            dt = t_j - t_i
            # If events occur at the exact same millisecond, force dt to be small but non-zero
            if dt <= 0:
                dt = 1e-6 
                
            L_ij = haversine_distance(lat_i, lon_i, lat_j, lon_j)
            if L_ij < eps_dist:
                L_ij = eps_dist
                
            eta = dt * (10**(-b_val * (M_i - Mc))) * (L_ij**df_val)
            
            if eta < min_eta:
                min_eta = eta
                best_parent = i
                
        parents.append(best_parent)
        parent_distances.append(min_eta)
        
    df['parent_idx'] = parents
    df['parent_eta'] = parent_distances
    
    # Create the complete NetworkX graph
    G_full = nx.Graph()
    G_full.add_nodes_from(range(n_events))
    for j in range(1, n_events):
        if parents[j] is not None:
            G_full.add_edge(parents[j], j)
            
    # 2. Sliding Window Geometric Analysis
    # We track how the curvature of the network changes over time.
    # In each window, we extract the induced subgraph of events in that window.
    # We calculate the Forman-Ricci Curvature (FRC) of the subgraph.
    
    # Pre-calculate FRC on the full graph G_full
    degrees_full = dict(G_full.degree())
    node_curvatures_full = {}
    for u in G_full.nodes():
        deg_u = degrees_full[u]
        if deg_u == 0:
            node_curvatures_full[u] = 0.0
        else:
            curv_sum = 0.0
            for v in G_full.neighbors(u):
                deg_v = degrees_full[v]
                curv_sum += (4.0 - deg_u - deg_v)
            node_curvatures_full[u] = curv_sum

    times = []
    avg_curvatures = []
    std_curvatures = []
    negative_curv_frac = []
    
    for start in range(0, n_events - window_size + 1, step_size):
        end = start + window_size
        nodes_in_window = list(range(start, end))
        
        # Get pre-calculated curvatures for nodes in this window (excluding mainshock seeds 0 and 1 to prevent scaling distortion)
        curv_vals = [node_curvatures_full[u] for u in nodes_in_window if u >= 2]
        
        avg_c = np.mean(curv_vals)
        std_c = np.std(curv_vals)
        neg_frac = np.mean([1.0 if c < 0 else 0.0 for c in curv_vals])
        
        # Midpoint time of window
        t_mid = df.loc[start:end-1, 'time_days'].mean()
        
        times.append(t_mid)
        avg_curvatures.append(avg_c)
        std_curvatures.append(std_c)
        negative_curv_frac.append(neg_frac)
        
    # Save geometry report
    basename = os.path.basename(csv_path).replace('.csv', '')
    report_path = os.path.join(BASE_DIR, f"../Reports/{basename}_geometry_report.csv")
    report_df = pd.DataFrame({
        'time_days': times,
        'avg_forman_curvature': avg_curvatures,
        'std_forman_curvature': std_curvatures,
        'negative_curvature_fraction': negative_curv_frac
    })
    report_df.to_csv(report_path, index=False)
    print(f"Geometry analysis report saved to: {report_path}")
    
    # 3. Plotting Curvature Evolution
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Average curvature plot
    ax1.plot(times, avg_curvatures, 'o-', color='#d62728', linewidth=2, label='Mean Forman-Ricci Curvature $\langle F \\rangle$')
    ax1.fill_between(times, 
                     np.array(avg_curvatures) - np.array(std_curvatures), 
                     np.array(avg_curvatures) + np.array(std_curvatures), 
                     color='#d62728', alpha=0.15, label='$\pm 1$ Std Dev')
    ax1.set_ylabel('Forman-Ricci Curvature', fontsize=12)
    ax1.set_title(f'Geometric Curvature Evolution of Seismic Network: {basename.replace("_", " ").title()}', fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='upper right')
    
    # Negative curvature fraction plot
    ax2.plot(times, negative_curv_frac, 'D-', color='#9467bd', linewidth=2, label='Fraction of Negatively Curved Nodes')
    ax2.set_xlabel('Time (days since doublet mainshock)', fontsize=12)
    ax2.set_ylabel('Fraction $F(u) < 0$', fontsize=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='lower right')
    
    plt.tight_layout()
    plot_path = os.path.join(BASE_DIR, f"../Reports/{basename}_geometry_evolution.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Geometry evolution plot saved to: {plot_path}")
    
    # 4. Save the network structure to visualize spatial curvature
    # We will compute the node curvature on the full graph and save it in a csv for mapping
    degrees_full = dict(G_full.degree())
    node_curvatures_full = {}
    for u in G_full.nodes():
        deg_u = degrees_full[u]
        if deg_u == 0:
            node_curvatures_full[u] = 0.0
        else:
            curv_sum = 0.0
            for v in G_full.neighbors(u):
                deg_v = degrees_full[v]
                curv_sum += (4.0 - deg_u - deg_v)
            node_curvatures_full[u] = curv_sum
            
    df['forman_curvature'] = [node_curvatures_full[i] for i in range(n_events)]
    df['degree'] = [degrees_full[i] for i in range(n_events)]
    
    annotated_csv_path = os.path.join(BASE_DIR, f"../Reports/{basename}_annotated_network.csv")
    df.to_csv(annotated_csv_path, index=False)
    print(f"Annotated catalog with curvature saved to: {annotated_csv_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze seismic network geometry.')
    parser.add_argument('--input', type=str, default=os.path.join(BASE_DIR, 'simulated_aftershocks.csv'),
                        help='Path to the seismic data CSV file')
    parser.add_argument('--b', type=float, default=1.0, help='Gutenberg-Richter b-value')
    parser.add_argument('--df', type=float, default=2.0, help='Fractal dimension')
    parser.add_argument('--mc', type=float, default=2.0, help='Cutoff magnitude Mc')
    parser.add_argument('--window', type=int, default=100, help='Sliding window size')
    parser.add_argument('--step', type=int, default=10, help='Sliding step size')
    
    args = parser.parse_args()
    
    # Make sure Reports directory exists
    os.makedirs(os.path.join(BASE_DIR, "../Reports"), exist_ok=True)
    
    analyze_geometry(args.input, args.b, args.df, args.mc, args.window, args.step)
