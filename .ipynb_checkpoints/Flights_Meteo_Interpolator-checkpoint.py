# Latest Run 01/17/2024
import Nio
import xarray as xr
import math
import glob
import pandas as pd
import numpy as np
from metpy.calc import dewpoint_from_relative_humidity
from metpy.calc import dewpoint_from_specific_humidity
from metpy.units import units
import metpy.calc
from metpy.units import units
from scipy.interpolate import interp1d, RegularGridInterpolator
from datetime import datetime, timezone, timedelta
from metpy.calc import height_to_geopotential
from tqdm import trange, tqdm
# from joblib import Parallel, delayed

import os
import shutil
import sys
from ecmwfapi import *

###-------------------------------------------------------------------------------------------------------###

def utc_to_unix(date_str):
    # decode byte object into string
    date_str = date_str.decode()
    # Parse the date string into a datetime object
    date_format = '%m/%d/%Y (%H:%M)'
    date_utc = datetime.strptime(date_str, date_format).replace(tzinfo=timezone.utc)
    # Get the timezone offset in seconds
    timezone_offset = date_utc.utcoffset().total_seconds()
    # Convert the datetime object to Unix time (in seconds)
    unix_time = int((date_utc - timedelta(seconds=timezone_offset)).timestamp())
    return unix_time

def normalized_distance(x, x_0, x_1):
    return (x-x_0)/(x_1-x_0)


def MIND_BLOWING_INTERPOLATION(flight_time, flight_latitude, flight_longitude, flight_geopotential, desired_index):
    
    # flight_time, flight_latitude, flight_longitude, flight_geopotential must be arrays
    
    # !!!!!!!!!!   VERY IMPORTANT NOTES:
    
    # lat, lon, time, plvl must exist in global environment
    
    # lat must be in descending order [90:-90]
    # plvl levels must be in descending order
    # lon must be in ascending order [0:360]
    # time must be in ascending order in UNIX format
    
    
    if len(flight_time) == len(flight_latitude) == len(flight_longitude) == len(flight_geopotential):
        pass
    else:
        raise ValueError('Given array are of different lengths')
    
    
    index = desired_index
    i = flight_time[index]
    # j = unkown
    k = flight_latitude[index]
    w = flight_longitude[index]
    t = flight_geopotential[index]
    
    # for i - Time axis interpolation
    
    id_i_0 = np.searchsorted(time, i)-1
    id_i_1 = np.searchsorted(time, i)
    nd_i = normalized_distance(i, time[id_i_0], time[id_i_1])
    data_reduced_i = z_data[id_i_0,:,:] + (z_data[id_i_1,:,:] - z_data[id_i_0,:,:])*nd_i
    
    
    # for k - latitude axis interpolation
    
     # BE VEEERY CAREFUL
        
    lat_0 = np.sort(lat)
    id_k_0 = abs(np.searchsorted(lat_0, k)-1 - (len(lat)-1))
    id_k_1 = abs(np.searchsorted(lat_0, k) - (len(lat)-1))
    nd_k = normalized_distance(k, lat[id_k_0], lat[id_k_1])
    data_reduced_ik = data_reduced_i[:,id_k_0,:] + (data_reduced_i[:,id_k_1,:] - data_reduced_i[:,id_k_0,:])*nd_k
    
    
    # for w - longitude axis interpolation
    
    id_w_0 = np.searchsorted(lon, w)-1
    id_w_1 = np.searchsorted(lon, w)
    nd_w = normalized_distance(w, lon[id_w_0], lon[id_w_1])
    data_reduced_ikw = data_reduced_ik[:,id_w_0] + (data_reduced_ik[:,id_w_1] - data_reduced_ik[:,id_w_0])*nd_w
    
    
    # for t - reversed interpolation of pressure level for a given geopotential
    
     # BE VEEERY CAREFUL
        
    data_reduced_ikw_0 = np.sort(data_reduced_ikw)
    id_t_0 = abs(np.searchsorted(data_reduced_ikw_0, t)-1 - (len(data_reduced_ikw)-1))
    id_t_1 = abs(np.searchsorted(data_reduced_ikw_0, t) - (len(data_reduced_ikw)-1))
    
    nd_t = normalized_distance(t, data_reduced_ikw[id_t_0], data_reduced_ikw[id_t_1])
    result = plvl[id_t_0] + (plvl[id_t_1] - plvl[id_t_0])*nd_t
    
    return result

reverse_map_function = lambda lon: (lon + 360) if (lon < 0) else lon


def extend_dim(arr):
    new_element = np.copy(arr[:,:,:,0])
    new_element = np.expand_dims(new_element, axis=-1)
    extended_array = np.concatenate((arr, new_element), axis=-1)
    
    return extended_array

###-------------------------------------------------------------------------------------------------------###
# Getting Flight data paths
if len(sys.argv) != 2:
    print("Usage: python script_name.py <directory_path>")
    sys.exit(1)

# Get the directory path from the command-line argument
directory_path = sys.argv[1]

save_path = f'interpolated_flights_{directory_path}/'

# Get a list of file names in the directory
file_names = [filename for filename in os.listdir(directory_path) if filename.endswith('.pickle')]

# Making Reuired Directories
if not os.path.exists(save_path):
    os.makedirs(save_path) 
    
if not os.path.exists('meteo/'):
    os.makedirs('meteo/') 

# Your list of flight dates
for file_name in tqdm(file_names):
    file_name_exts = file_name.split('.')[0]
    print(f'Starting Operation for {file_name_exts} \n')
    print('       __|__')
    print('--o--o--(_)--o--o--\n')
    flight_date = pd.read_pickle(f'{directory_path}/{file_name}')['datetime'][0].strftime('%Y-%m-%d')
    flight_date_obj = datetime.strptime(flight_date, '%Y-%m-%d')
    next_day_date_obj = flight_date_obj + timedelta(days=1)
    next_day_date_str = next_day_date_obj.strftime('%Y-%m-%d')
    
    flight_date = str(flight_date_obj).split(" ")[0]
    next_day_date = str(next_day_date_obj).split(" ")[0]


    if not os.path.exists(f"meteo/{flight_date}.grib"):
        # Parsing ECMWF DATA 
        # yyyy/mm/dd   
        print(f'Getting ECMWF DATA {flight_date}...\n')
            
        my_param = "129.128/157.128"
        all_param = "129.128/130.128/131/132/133.128/157.128"
        server = ECMWFService("mars")
    
        server.execute(
            {
            "class": "od",
            "date": f"{flight_date}",
            "expver": "1",
            "levelist": "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/800/850/900/925/950/1000",
            "levtype": "pl",
            "param": all_param,
            "stream": "oper",
            "time": "00:00:00/06:00:00/12:00:00/18:00:00",
            "type": "an",
            "grid": "0.1/0.1",
            },
            f"meteo/{flight_date}.grib")
    
        server.execute(
            {
            "class": "od",
            "date": f"{next_day_date}", 
            "expver": "1",
            "levelist": "1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/800/850/900/925/950/1000",
            "levtype": "pl",
            "param": all_param,
            "stream": "oper",
            "time": "00:00:00",
            "type": "an",
            "grid": "0.1/0.1",
            },
            f"meteo/{next_day_date}_.grib")
        
        print('ECMWF DATA downloaded!\n')

###-------------------------------------------------------------------------------------------------------###
    print("Loading data to memory...\n")
    #Import weather data grid
    path = f"meteo/{flight_date}.grib"
    path_2 = f"meteo/{next_day_date}_.grib"
    ds = xr.open_dataset(path, engine="pynio")
    ds_2 = xr.open_dataset(path_2, engine="pynio")

    #Input flight data
    flight = pd.read_pickle(f"{directory_path}/{file_name}") #First 3 points
    flight = flight[(flight['geoaltitude'] <= 15000) & (flight['geoaltitude'] >= 9000)]
    unique_dates = flight['datetime'].dt.date.unique()
    flight = flight[flight['datetime'].dt.date==unique_dates[0]] # Temp solution 
    flight.dropna(subset=['time','lat','lon','geoaltitude'],inplace=True)
    ftime = flight["time"].values
    flat = flight["lat"].values
    flon = flight["lon"].values
    flon = np.vectorize(reverse_map_function)(flon)
    falt = flight["geoaltitude"].values

    height = falt * units.m
    geopot = height_to_geopotential(height).magnitude

    #Convert data to array
    time_data =ds.isel()["initial_time0"].to_numpy()
    z_data = ds.isel()["Z_GDS0_ISBL"].to_numpy()
    z_data_4 = ds_2.isel()["Z_GDS0_ISBL"].to_numpy() #Next day
    z_data = np.concatenate((z_data,z_data_4.reshape(1,25,1801,3600)),axis=0)

    del z_data_4

    t_data = ds.isel()["T_GDS0_ISBL"].to_numpy()
    t_data_4 = ds_2.isel()["T_GDS0_ISBL"].to_numpy()
    t_data = np.concatenate((t_data,t_data_4.reshape(1,25,1801,3600)),axis=0) 

    del t_data_4

    u_data = ds.isel()["U_GDS0_ISBL"].to_numpy()
    u_data_4 = ds_2.isel()["U_GDS0_ISBL"].to_numpy()
    u_data = np.concatenate((u_data,u_data_4.reshape(1,25,1801,3600)),axis=0)

    del u_data_4

    v_data = ds.isel()["V_GDS0_ISBL"].to_numpy()
    v_data_4 = ds_2.isel()["V_GDS0_ISBL"].to_numpy()
    v_data = np.concatenate((v_data,v_data_4.reshape(1,25,1801,3600)),axis=0)

    del v_data_4


    q_data = ds.isel()["Q_GDS0_ISBL"].to_numpy()
    q_data_4 = ds_2.isel()["Q_GDS0_ISBL"].to_numpy()
    q_data = np.concatenate((q_data,q_data_4.reshape(1,25,1801,3600)),axis=0)

    del q_data_4


    r_data = ds.isel()["R_GDS0_ISBL"].to_numpy()
    r_data_4 = ds_2.isel()["R_GDS0_ISBL"].to_numpy()
    r_data = np.concatenate((r_data,r_data_4.reshape(1,25,1801,3600)),axis=0)

    del r_data_4
      
    z_data = extend_dim(z_data)
    t_data = extend_dim(t_data)
    u_data = extend_dim(u_data)
    v_data = extend_dim(v_data)
    q_data = extend_dim(q_data)
    r_data = extend_dim(r_data)

    #Import time data
    byte_literal = time_data[0]
    string_value = byte_literal.decode('ascii')
    datetime_format = "%m/%d/%Y (%H:%M)"
    datetime_value = datetime.strptime(string_value, datetime_format)
    one_day = timedelta(days=1)
    new_datetime = datetime_value + one_day
    new_string_value = new_datetime.strftime(datetime_format)
    new_byte_literal = new_string_value.encode('ascii')
    time_data = np.append(time_data,new_byte_literal)

    #Allocate lon lat plvl time
    # map_function = lambda lon: (lon - 360) if (lon > 180) else lon
    time = []
    for i in range(len(time_data)) : 
        time.append(utc_to_unix(time_data[i]))

    # lon =  np.vectorize(map_function)(lon_)
    lon= np.round(np.arange(0,360.1,0.1,dtype='float32'),1)
    lat = np.round(np.arange(90,-90.1,-0.1,dtype='float32'),1)
    plvl = np.array([   1,    2,    3,    5,    7,   10,   20,   30,   50,   70,  100,  150,
            200,  250,  300,  400,  500,  600,  700,  800,  850,  900,  925,  950,
           1000], dtype='int32')

    print("Data Loading Done !")

###-------------------------------------------------------------------------------------------------------###
    print('Interpolating the data...')
    
    # interp_z = RegularGridInterpolator((time, plvl, lat, lon),z_data)
    interp_t = RegularGridInterpolator((time, plvl, lat, lon),t_data)
    interp_u = RegularGridInterpolator((time, plvl, lat, lon),u_data)
    interp_v = RegularGridInterpolator((time, plvl, lat, lon),v_data)
    interp_q = RegularGridInterpolator((time, plvl, lat, lon),q_data)
    interp_r = RegularGridInterpolator((time, plvl, lat, lon),r_data)
    
    met_at_flight = []
    
    for i in trange(len(flight)):
            plvl_point = MIND_BLOWING_INTERPOLATION(ftime, flat, flon, geopot, i)
    
            s =[plvl_point,
            interp_t((ftime[i],plvl_point,flat[i],flon[i])),
            interp_u((ftime[i],plvl_point,flat[i],flon[i])),
            interp_v((ftime[i],plvl_point,flat[i],flon[i])),
            interp_q((ftime[i],plvl_point,flat[i],flon[i])),
            interp_r((ftime[i],plvl_point,flat[i],flon[i]))]
    
            met_at_flight.append(s)
        
    df_temp = pd.DataFrame(np.array(met_at_flight),columns=['p', 'temp', 'u', 'v', 'spfh', 'rh'])
    result_df = pd.concat([flight.reset_index(drop=True), df_temp.reset_index(drop=True)], axis=1)
    
    result_df.to_csv(fr'{save_path}/{file_name_exts}_meteo.csv')
    print('Interpolation Done !')
    # Deleting Meteo Data
    # os.remove(f"meteo/{flight_date}.grib")
    # os.remove(f"meteo/{next_day_date}.grib")

folder_path = "meteo"
# Use shutil.rmtree to remove the folder and its contents
try:
    shutil.rmtree(folder_path)
    print(f"Deleted folder: {folder_path}")
except Exception as e:
    print(f"Error deleting folder {folder_path}: {e}")