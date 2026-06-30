import os
import csv
import math

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
csv_file = os.path.join(BASE_DIR, "caracas_seismic_data.csv")

try:
    events = []
    with open(csv_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            row['magnitude'] = float(row['magnitude'])
            row['depth_km'] = float(row['depth_km']) if row['depth_km'] else 0.0
            row['energy_joules'] = 10**(4.8 + 1.5 * row['magnitude'])
            events.append(row)
            
    if not events:
        print("No events found in the CSV.")
    else:
        print("--- Seismic Catalog Summary ---")
        print(f"Total events recorded: {len(events)}")
        
        # Sort by time
        events_sorted = sorted(events, key=lambda x: x['time_utc'])
        
        print(f"Date range: {events_sorted[0]['time_utc']} to {events_sorted[-1]['time_utc']}")
        
        magnitudes = [e['magnitude'] for e in events]
        print(f"Magnitude range: {min(magnitudes)} to {max(magnitudes)}")
        
        total_energy = sum(e['energy_joules'] for e in events)
        print(f"Total seismic energy released: {total_energy:.2e} Joules")
        
        print("\nEvents list:")
        for e in events_sorted:
            energy_pct = (e['energy_joules'] / total_energy) * 100
            print(f"- UTC: {e['time_utc']} | M {e['magnitude']} ({e['mag_type']}) | Depth: {e['depth_km']}km | Energy: {e['energy_joules']:.2e} J ({energy_pct:.2f}%) | Location: {e['place']}")
            
except Exception as e:
    print(f"Error reading or parsing CSV: {e}")
