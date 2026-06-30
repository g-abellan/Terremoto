import os
import pandas as pd
import numpy as np
import math

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
df_params = pd.read_csv(os.path.join(BASE_DIR, '../Reports/simulated_aftershocks_field_parameters.csv'))
params = dict(zip(df_params['parameter'], df_params['value']))

mu = params['mu']
K0 = params['K0']
c = params['c']
p = params['p']
d = params['d']
q = params['q']
alpha_val = params['alpha']
Mc = 2.0

t1 = 0.0
t2 = 39.0 / (24 * 3600)
M1, M2 = 7.2, 7.5
lat1, lon1 = 10.436, -68.5277
lat2, lon2 = 10.4351, -68.4716

T1 = 5.0  # June 29
T2 = 35.0 # July 29

def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2)**2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

# Bounding box for entire catalog
cat_lat_min, cat_lat_max = 9.5, 11.5
cat_lon_min, cat_lon_max = -70.0, -65.0
cat_width = haversine_distance((cat_lat_min + cat_lat_max)/2.0, cat_lon_min, (cat_lat_min + cat_lat_max)/2.0, cat_lon_max)
cat_height = haversine_distance(cat_lat_min, (cat_lon_min + cat_lon_max)/2.0, cat_lat_max, (cat_lon_min + cat_lon_max)/2.0)
cat_area = cat_width * cat_height

# Caracas sub-basin
lat_min, lat_max = 10.45, 10.55
lon_min, lon_max = -67.0, -66.8
sub_area = haversine_distance((lat_min+lat_max)/2.0, lon_min, (lat_min+lat_max)/2.0, lon_max) * haversine_distance(lat_min, (lon_min+lon_max)/2.0, lat_max, (lon_min+lon_max)/2.0)

# Time integrals
if abs(p - 1.0) < 1e-5:
    int_t1 = math.log((T2 - t1 + c) / (T1 - t1 + c))
    int_t2 = math.log((T2 - t2 + c) / (T1 - t2 + c))
else:
    int_t1 = ((T2 - t1 + c)**(1.0 - p) - (T1 - t1 + c)**(1.0 - p)) / (1.0 - p)
    int_t2 = ((T2 - t2 + c)**(1.0 - p) - (T1 - t2 + c)**(1.0 - p)) / (1.0 - p)

# 1. Regional Triggered (over R^2)
int_s = math.pi * (d**(2.0 - 2.0 * q)) / (q - 1.0)
N_trig1_reg = math.exp(alpha_val * (M1 - Mc)) * K0 * int_t1 * int_s
N_trig2_reg = math.exp(alpha_val * (M2 - Mc)) * K0 * int_t2 * int_s
N_bg_reg = mu * (T2 - T1)
N_total_reg = N_bg_reg + N_trig1_reg + N_trig2_reg

# 2. Caracas Triggered (spatial grid integration)
lat_points = np.linspace(lat_min, lat_max, 100)
lon_points = np.linspace(lon_min, lon_max, 100)
dlat_km = haversine_distance(lat_min, (lon_min+lon_max)/2.0, lat_max, (lon_min+lon_max)/2.0) / 100.0
dlon_km = haversine_distance((lat_min+lat_max)/2.0, lon_min, (lat_min+lat_max)/2.0, lon_max) / 100.0
dA = dlat_km * dlon_km

int_s1_car = 0.0
int_s2_car = 0.0
for lat in lat_points:
    for lon in lon_points:
        r1 = haversine_distance(lat1, lon1, lat, lon)
        r2 = haversine_distance(lat2, lon2, lat, lon)
        int_s1_car += (1.0 / ((r1**2 + d**2)**q)) * dA
        int_s2_car += (1.0 / ((r2**2 + d**2)**q)) * dA

N_trig1_car = math.exp(alpha_val * (M1 - Mc)) * K0 * int_t1 * int_s1_car
N_trig2_car = math.exp(alpha_val * (M2 - Mc)) * K0 * int_t2 * int_s2_car
N_bg_car = (mu / cat_area) * sub_area * (T2 - T1)
N_total_car = N_bg_car + N_trig1_car + N_trig2_car

print(f"Regional total area: {cat_area:.1f} km^2")
print(f"Expected regional total events: {N_total_reg:.2f}")
print(f"Expected regional background: {N_bg_reg:.2f}")
print(f"Expected regional triggered: {N_trig1_reg + N_trig2_reg:.2f}")
print(f"Expected Caracas total events: {N_total_car:.4f}")
print(f"Expected Caracas background: {N_bg_car:.4f}")
print(f"Expected Caracas triggered: {N_trig1_car + N_trig2_car:.4f}")
