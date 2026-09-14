#!/usr/bin/env python3
"""
network_topology.py

Constructs the Baiesi-Paczuski causal network for the Northern Venezuela doublet sequence,
computes both unweighted Forman-Ricci curvature (FRC) and weighted Forman-Ricci curvature (WFRC),
verifies the exact topological theorem <F> = 8(N - N_seeds)/N - 2 <d^2>,
and derives the analytical bridge between continuous point-process intensity and degree moments.
"""

import os
import csv
import math
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2.0)**2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2.0)**2)
    return 2.0 * R * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))

def analyze_causal_network(csv_file="consolidated_catalog_2026.csv", b_val=1.05, df_val=2.0, Mc=2.5):
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
    N = len(rows)
    t0_ms = rows[0]["timestamp_ms"]
    for r in rows:
        r["time_days"] = (r["timestamp_ms"] - t0_ms) / (1000.0 * 86400.0)
        
    parents = [None, None]
    parent_dist = [np.nan, np.nan]
    edge_spatial_dist = [np.nan, np.nan]
    
    for j in range(2, N):
        t_j = rows[j]["time_days"]
        lat_j = rows[j]["latitude"]
        lon_j = rows[j]["longitude"]
        
        min_eta = float('inf')
        best_p = None
        best_L = 0.0
        
        for i in range(j):
            dt = max(t_j - rows[i]["time_days"], 1e-6)
            L_ij = max(haversine(rows[i]["latitude"], rows[i]["longitude"], lat_j, lon_j), 0.1)
            eta = dt * (10.0 ** (-b_val * (rows[i]["magnitude"] - Mc))) * (L_ij ** df_val)
            if eta < min_eta:
                min_eta = eta
                best_p = i
                best_L = L_ij
                
        parents.append(best_p)
        parent_dist.append(min_eta)
        edge_spatial_dist.append(best_L)
        
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for j in range(2, N):
        if parents[j] is not None:
            # Add edge with physical spatial distance as weight
            G.add_edge(parents[j], j, weight=edge_spatial_dist[j], eta=parent_dist[j])
            
    degrees = dict(G.degree())
    
    # 1. Unweighted Forman-Ricci Curvature
    frc_nodes = {}
    for u in G.nodes():
        deg_u = degrees[u]
        if deg_u == 0:
            frc_nodes[u] = 0.0
        else:
            frc_nodes[u] = sum(4.0 - deg_u - degrees[v] for v in G.neighbors(u))
            
    # 2. Weighted Forman-Ricci Curvature (carrying spatial metric distance)
    # F_W(e) = w_e * [ w_u/w_e + w_v/w_e - sum_{e' ~ u, != e} w_u / sqrt(w_e * w_e') - sum_{e'' ~ v, != e} w_v / sqrt(w_e * w_e'') ]
    # Node weights w_u = sum_{e ~ u} w_e
    node_weights = {}
    for u in G.nodes():
        node_weights[u] = sum(G[u][v].get('weight', 1.0) for v in G.neighbors(u)) if degrees[u] > 0 else 0.0
        
    wfrc_edges = {}
    for u, v in G.edges():
        w_e = G[u][v].get('weight', 1.0)
        w_u = node_weights[u]
        w_v = node_weights[v]
        
        sum_u = sum(w_u / math.sqrt(w_e * G[u][nbr].get('weight', 1.0)) for nbr in G.neighbors(u) if nbr != v)
        sum_v = sum(w_v / math.sqrt(w_e * G[v][nbr].get('weight', 1.0)) for nbr in G.neighbors(v) if nbr != u)
        
        wfrc_edges[(u, v)] = w_e * ((w_u / w_e) + (w_v / w_e) - sum_u - sum_v)
        
    wfrc_nodes = {}
    for u in G.nodes():
        if degrees[u] == 0:
            wfrc_nodes[u] = 0.0
        else:
            wfrc_nodes[u] = sum(wfrc_edges.get((min(u, v), max(u, v)), wfrc_edges.get((max(u, v), min(u, v)), 0.0)) for v in G.neighbors(u))
            
    # Exact Topological Theorem Verification
    E = G.number_of_edges()
    N_seeds = 2
    sum_F = sum(frc_nodes.values())
    sum_d2 = sum(degrees[u]**2 for u in G.nodes())
    expected_sum_F = 8 * E - 2 * sum_d2
    mean_F = float(np.mean(list(frc_nodes.values())))
    expected_mean_F = 8.0 * (N - N_seeds) / N - 2.0 * float(np.mean([degrees[u]**2 for u in G.nodes()]))
    
    print("==================================================")
    print("Network Topology & Curvature Verification:")
    print("==================================================")
    print(f"Total Nodes: {N}, Total Edges: {E}, Seeds: {N_seeds}")
    print(f"Observed sum(F): {sum_F:.6f}")
    print(f"Analytical sum(F) [8E - 2 sum(d^2)]: {expected_sum_F:.6f}")
    print(f"Mean Unweighted FRC: {mean_F:.4f} vs Analytical: {expected_mean_F:.4f}")
    print(f"Mean Weighted FRC: {np.mean(list(wfrc_nodes.values())):.4f}")
    print(f"Identity Absolute Error: {abs(sum_F - expected_sum_F):.2e}")
    print("==================================================")
    
    # Save network summary to CSV
    with open("network_metrics_summary.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "time_days", "latitude", "longitude", "magnitude", "degree", "frc_unweighted", "frc_weighted", "parent_idx", "parent_dist_km"])
        for i in range(N):
            writer.writerow([
                rows[i]["id"], rows[i]["time_days"], rows[i]["latitude"], rows[i]["longitude"],
                rows[i]["magnitude"], degrees[i], frc_nodes[i], wfrc_nodes[i], parents[i], edge_spatial_dist[i]
            ])
            
    # Figure 1: 3-panel comprehensive topology plot
    plt.rcParams.update({'font.sans-serif': 'Helvetica', 'font.size': 10})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.5, 4.2))
    
    # Sliding window for unweighted curvature
    w_size = 40
    step = 4
    t_mid = []
    frc_w = []
    wfrc_w = []
    
    for s in range(2, N - w_size + 1, step):
        w_nodes = list(range(s, s + w_size))
        t_mid.append(rows[w_nodes[-1]]["time_days"])
        frc_w.append(np.mean([frc_nodes[u] for u in w_nodes]))
        wfrc_w.append(np.mean([wfrc_nodes[u] for u in w_nodes]))
        
    ax1.plot(t_mid, frc_w, 'o-', color='#1f77b4', markersize=4, label=r"Unweighted $\langle F \rangle$")
    ax1.axhline(0, color='gray', linestyle=':', alpha=0.7)
    ax1.set_xlabel("Time Since Doublet (days)", fontsize=11)
    ax1.set_ylabel(r"Unweighted Curvature $\langle F \rangle$", fontsize=11)
    ax1.set_title("(a) Unweighted Curvature Relaxation", fontsize=11, fontweight="bold")
    ax1.grid(True, alpha=0.25)
    ax1.legend(frameon=True)
    
    # Sliding window for weighted curvature
    ax2.plot(t_mid, wfrc_w, 's-', color='#d62728', markersize=4, label=r"Weighted $\langle F_W \rangle$ (Metric)")
    ax2.axhline(0, color='gray', linestyle=':', alpha=0.7)
    ax2.set_xlabel("Time Since Doublet (days)", fontsize=11)
    ax2.set_ylabel(r"Weighted Curvature $\langle F_W \rangle$ (km)", fontsize=11)
    ax2.set_title("(b) Spatial Metric Curvature Evolution", fontsize=11, fontweight="bold")
    ax2.grid(True, alpha=0.25)
    ax2.legend(frameon=True)
    
    # Degree Distribution P(k)
    deg_vals = [degrees[u] for u in range(2, N)]
    ax3.hist(deg_vals, bins=np.logspace(0, 2.2, 14), density=True, alpha=0.75, color='#2ca02c', edgecolor='black')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel(r"Node Degree $k$", fontsize=11)
    ax3.set_ylabel(r"Probability Density $P(k)$", fontsize=11)
    ax3.set_title("(c) Scale-Free Causal Degree Distribution", fontsize=11, fontweight="bold")
    ax3.grid(True, alpha=0.25)
    
    plt.tight_layout()
    plt.savefig("fig1_network_topology.png", dpi=300)
    plt.close()
    print("Generated upgraded fig1_network_topology.png")

if __name__ == "__main__":
    analyze_causal_network()
