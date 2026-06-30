import pandas as pd
import numpy as np
import math

import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(BASE_DIR, "caracas_seismic_data.csv")
output_file = os.path.join(BASE_DIR, "simulated_aftershocks.csv")

# ETAS Model Parameters (typical for tectonic regions)
Mc = 2.0        # Cutoff magnitude for simulation
b_val = 1.0     # Gutenberg-Richter b-value
alpha = 1.4     # Magnitude efficiency parameter (efficiency of triggering aftershocks)
K0 = 0.08       # Productivity constant (calibrated to get ~500 events in branching cascade)
c_val = 0.01    # Omori time offset (in days, ~15 mins)
p_val = 1.15    # Omori decay exponent (must be > 1)
d_val = 0.08    # Spatial scaling parameter (in degrees)
q_val = 1.5     # Spatial decay exponent (must be > 1)

np.random.seed(42)  # For reproducibility

try:
    df_seeds = pd.read_csv(csv_file)
    print(f"Loaded {len(df_seeds)} seed events from USGS.")
    
    # We will represent events as a list of dictionaries
    # Keys: id, time_days, magnitude, latitude, longitude, parent_id, generation
    events = []
    
    # Convert seed times to days since the first event
    t0_ms = df_seeds['timestamp_ms'].min()
    
    for idx, row in df_seeds.iterrows():
        t_days = (row['timestamp_ms'] - t0_ms) / (1000.0 * 60 * 60 * 24)
        events.append({
            'id': row['id'],
            'time_days': t_days,
            'timestamp_ms': row['timestamp_ms'],
            'latitude': row['latitude'],
            'longitude': row['longitude'],
            'magnitude': row['magnitude'],
            'parent_id': 'NONE',
            'generation': 0
        })
        
    # Branching process simulation
    current_gen = 0
    max_generations = 10
    
    # We loop through events and generate offspring
    # To do a branching process, we can keep a queue of events to process
    queue = list(events)
    
    simulated_events = []
    
    event_counter = 1
    
    while queue:
        parent = queue.pop(0)
        simulated_events.append(parent)
        
        gen = parent['generation']
        if gen >= max_generations:
            continue
            
        # Calculate mean number of offspring: mu = K0 * exp(alpha * (M - Mc))
        M = parent['magnitude']
        mu = K0 * math.exp(alpha * (M - Mc))
        
        # Sample actual number of offspring from a Poisson distribution
        num_offspring = np.random.poisson(mu)
        
        for _ in range(num_offspring):
            # 1. Magnitude: Gutenberg-Richter law (exponential distribution)
            # P(M) = b * ln(10) * 10^(-b*(M-Mc))
            u_mag = np.random.uniform(0, 1)
            mag = Mc - math.log10(1.0 - u_mag) / b_val
            
            # Round magnitude to 1 decimal place
            mag = round(mag, 1)
            
            # 2. Time: Omori-Utsu law (power-law distribution)
            # f(dt) = (p-1) * c^(p-1) * (dt + c)^(-p)
            u_time = np.random.uniform(0, 1)
            # Inverse CDF of Omori law:
            dt = c_val * ((1.0 - u_time)**(1.0 / (1.0 - p_val)) - 1.0)
            
            # Limit dt to 10 days max to prevent extremely long delays
            dt = min(dt, 10.0)
            t_days = parent['time_days'] + dt
            
            # 3. Location: Isotropic spatial decay f(r) ~ r * (r^2 + d^2)^(-q)
            # We sample r using inverse CDF of f(r)
            # CDF: F(r) = 1 - (d^2 / (r^2 + d^2))^(q-1) for q > 1
            u_loc = np.random.uniform(0, 1)
            r = d_val * math.sqrt((1.0 - u_loc)**(1.0 / (1.0 - q_val)) - 1.0)
            
            # Limit r to 1.5 degrees max to keep aftershocks in the regional fault zone
            r = min(r, 1.5)
            
            theta = np.random.uniform(0, 2 * math.pi)
            dlon = r * math.cos(theta) / math.cos(math.radians(parent['latitude']))
            dlat = r * math.sin(theta)
            
            lat = parent['latitude'] + dlat
            lon = parent['longitude'] + dlon
            
            child_id = f"sim_{event_counter}"
            event_counter += 1
            
            child = {
                'id': child_id,
                'time_days': t_days,
                'timestamp_ms': t0_ms + int(t_days * 24 * 60 * 60 * 1000),
                'latitude': lat,
                'longitude': lon,
                'magnitude': mag,
                'parent_id': parent['id'],
                'generation': gen + 1
            }
            
            # Add to queue for next generation triggering, and append to database
            # We only allow triggering if magnitude is large enough (e.g. M >= Mc)
            if mag >= Mc:
                queue.append(child)
                
    # Create DataFrame and sort by time
    df_sim = pd.DataFrame(simulated_events)
    df_sim = df_sim.sort_values(by='time_days')
    
    # Re-calculate datetime strings
    df_sim['time_utc'] = df_sim['timestamp_ms'].apply(
        lambda x: pd.to_datetime(x, unit='ms').strftime('%Y-%m-%d %H:%M:%S')
    )
    
    # Save to file
    df_sim.to_csv(output_file, index=False)
    
    print("\n--- Simulation Complete ---")
    print(f"Generated a total of {len(df_sim)} events in the catalog.")
    print(f"Number of simulated aftershocks: {len(df_sim) - len(df_seeds)}")
    print(f"Cutoff magnitude: M {Mc}")
    print(f"Max generation reached: {df_sim['generation'].max()}")
    print(f"Simulated catalog saved to: {output_file}")
    
except Exception as e:
    print(f"Error during ETAS simulation: {e}")
