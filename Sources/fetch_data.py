import urllib.request
import json
import csv
from datetime import datetime

# Define parameters for the query
# Bounding box for northwestern/central Venezuela (San Felipe to Caracas)
min_lat = 9.5
max_lat = 11.5
min_lon = -70.0
max_lon = -65.0

start_time = "2026-06-24T00:00:00"
# Current date is June 29, 2026
end_time = "2026-06-30T00:00:00"

url = (
    f"https://earthquake.usgs.gov/fdsnws/event/1/query?"
    f"format=geojson&"
    f"starttime={start_time}&"
    f"endtime={end_time}&"
    f"minlatitude={min_lat}&"
    f"maxlatitude={max_lat}&"
    f"minlongitude={min_lon}&"
    f"maxlongitude={max_lon}"
)

print(f"Fetching data from: {url}")

try:
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read().decode('utf-8'))
        
    features = data.get('features', [])
    print(f"Successfully retrieved {len(features)} events.")
    import os
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    csv_file = os.path.join(BASE_DIR, "caracas_seismic_data.csv")
    
    with open(csv_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # Header row
        writer.writerow([
            "id", "time_utc", "timestamp_ms", "latitude", "longitude", 
            "depth_km", "magnitude", "mag_type", "place"
        ])
        
        for feat in features:
            props = feat.get('properties', {})
            geom = feat.get('geometry', {})
            coords = geom.get('coordinates', [0, 0, 0]) # lon, lat, depth
            
            event_id = feat.get('id')
            time_ms = props.get('time') # UTC timestamp in ms
            time_str = datetime.utcfromtimestamp(time_ms / 1000.0).strftime('%Y-%m-%d %H:%M:%S') if time_ms else ""
            
            magnitude = props.get('mag')
            mag_type = props.get('magType')
            place = props.get('place')
            
            lon = coords[0]
            lat = coords[1]
            depth = coords[2] if len(coords) > 2 else None
            
            writer.writerow([
                event_id, time_str, time_ms, lat, lon, depth, magnitude, mag_type, place
            ])
            
    print(f"Data saved to {csv_file}")
    
except Exception as e:
    print(f"Error fetching or parsing data: {e}")
