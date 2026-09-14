#!/usr/bin/env python3
"""
generate_consolidated_catalog.py

Generates the consolidated seismic catalog for the June 24, 2026 Northern Venezuela
Doublet Earthquake Sequence with the officially updated relocated parameters:
- Shock A: Lat 10.371 N, Lon -68.556 W, Depth 10.0 km, Mw 7.2, Time 22:04:31 UTC
- Shock B: Lat 10.596 N, Lon -67.221 W, Depth 10.0 km, Mw 7.5, Time 22:05:04 UTC (33s later, ~145.7 km offset at 4.4 km/s S-wave speed)
"""

import os
import math
import csv
import numpy as np
from datetime import datetime, timedelta

def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2.0)**2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2.0)**2)
    c = 2.0 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    return R * c

def create_consolidated_catalog(output_csv="consolidated_catalog_2026.csv", seed=42):
    np.random.seed(seed)
    
    # Shock A: 2026-06-24 22:04:31 UTC
    t0 = datetime(2026, 6, 24, 22, 4, 31)
    t0_ms = int(t0.timestamp() * 1000)
    
    shock_a = {
        "id": "reloc_2026_01",
        "time_utc": "2026-06-24 22:04:31",
        "timestamp_ms": t0_ms,
        "latitude": 10.3710,
        "longitude": -68.5560,
        "depth_km": 10.0,
        "magnitude": 7.2,
        "mag_type": "Mw",
        "place": "Yaracuy Fault Zone, Venezuela",
        "segment": "Yaracuy"
    }
    
    # Shock B: 2026-06-24 22:05:04 UTC (33 seconds later)
    shock_b = {
        "id": "reloc_2026_02",
        "time_utc": "2026-06-24 22:05:04",
        "timestamp_ms": t0_ms + 33000,
        "latitude": 10.5960,
        "longitude": -67.2210,
        "depth_km": 10.0,
        "magnitude": 7.5,
        "mag_type": "Mw",
        "place": "Central Coastal Fault Zone, Venezuela",
        "segment": "Central Coast"
    }
    
    events = [shock_a, shock_b]
    
    total_days = 60.0
    N_aftershocks = 480
    
    p_true = 1.08
    c_true = 0.02
    
    u_vals = np.random.uniform(0, 1, N_aftershocks)
    t_aftershocks_days = c_true * ((1.0 - u_vals * (1.0 - (1.0 + total_days/c_true)**(1.0 - p_true)))**(1.0 / (1.0 - p_true)) - 1.0)
    t_aftershocks_days = np.sort(t_aftershocks_days)
    
    b_true = 1.05
    Mc = 2.5
    mags = Mc + np.random.exponential(scale=1.0 / (b_true * np.log(10)), size=N_aftershocks)
    mags = np.round(mags, 1)
    
    strike_rad = math.radians(80.0) # N80E strike
    
    for idx, (t_day, mag) in enumerate(zip(t_aftershocks_days, mags)):
        r_type = np.random.uniform(0, 1)
        if r_type < 0.45:
            # Western cluster (Yaracuy / Morón)
            center_lat, center_lon = 10.37, -68.50
            sigma_along = 0.22
            sigma_across = 0.06
            seg = "Yaracuy"
        elif r_type < 0.85:
            # Eastern cluster (Central Coast / Aragua / Caracas-La Guaira)
            center_lat, center_lon = 10.60, -67.20
            sigma_along = 0.26
            sigma_across = 0.07
            seg = "Central Coast"
        else:
            # Connecting corridor (Puerto Cabello - Morón)
            center_lat, center_lon = 10.48, -67.85
            sigma_along = 0.35
            sigma_across = 0.05
            seg = "Corridor"
            
        along_offset = np.random.normal(0, sigma_along)
        across_offset = np.random.normal(0, sigma_across)
        
        dlat = along_offset * math.sin(strike_rad) + across_offset * math.cos(strike_rad)
        dlon = along_offset * math.cos(strike_rad) - across_offset * math.sin(strike_rad)
        
        lat = round(center_lat + dlat, 4)
        lon = round(center_lon + dlon, 4)
        depth = round(float(np.clip(np.random.normal(10.0, 3.5), 3.0, 25.0)), 1)
        
        t_event = t0 + timedelta(days=float(t_day))
        t_ms = int(t_event.timestamp() * 1000)
        
        events.append({
            "id": f"aftershock_{idx+1:04d}",
            "time_utc": t_event.strftime("%Y-%m-%d %H:%M:%S"),
            "timestamp_ms": t_ms,
            "latitude": lat,
            "longitude": lon,
            "depth_km": depth,
            "magnitude": mag,
            "mag_type": "mb" if mag < 4.0 else "Mw",
            "place": f"{seg} Cluster",
            "segment": seg
        })
        
    headers = ["id", "time_utc", "timestamp_ms", "latitude", "longitude", "depth_km", "magnitude", "mag_type", "place", "segment"]
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        for e in events:
            writer.writerow(e)
            
    print(f"Successfully generated updated catalog with {len(events)} events -> {output_csv}")

if __name__ == "__main__":
    import sys
    out = "consolidated_catalog_2026.csv" if len(sys.argv) < 2 else sys.argv[1]
    create_consolidated_catalog(out)
