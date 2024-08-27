import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

# Function to read HDF5 file and extract rw_obs data
def read_hdf5_file(filename, station_name=None, radar_tilt=None):
    with h5py.File(filename, 'r') as f:
        obs = f['ObsValue']['equivalentReflectivityFactor'][:]
        bak = f['hofx0']['equivalentReflectivityFactor'][:]
        ana = f['hofx3']['equivalentReflectivityFactor'][:]
        lat = f['MetaData']['latitude'][:]  # Assuming latitudes are stored here
        lon = f['MetaData']['longitude'][:]  # Assuming longitudes are stored here
        azi = f['MetaData']['radarAzimuth'][:]
        til = f['MetaData']['radarTilt'][:]
        nam = f['MetaData']['stationIdentification'][:].astype(str) 
    
    n = len(obs)
    # For specific station name and tilt
    if station_name is not None:
        station_mask = (nam == station_name)
    else:
        station_mask = np.ones(n, dtype=bool)
        
    if radar_tilt is not None:
        tilt_mask = (np.abs(til - radar_tilt) < 0.15)
    else:
        tilt_mask = np.ones(n, dtype=bool)

    combined_mask = station_mask & tilt_mask

    # Apply the combined mask to all arrays
    obs = obs[combined_mask]
    bak = bak[combined_mask]
    ana = ana[combined_mask]
    lat = lat[combined_mask]
    lon = lon[combined_mask]

    return obs, bak, ana, lat, lon

# Create a custom colormap
def create_custom_colormap():
    colors = ['white', 'blue', 'green', 'orange', 'red', 'purple']
    cmap_name = 'custom_cmap'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

def get_tick_interval(data_range, num_ticks=5):
    """Determine a suitable tick interval based on the data range."""
    min_val, max_val = data_range
    range_val = max_val - min_val
    
    # Define possible tick intervals
    intervals = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    
    # Choose an interval based on the range
    for interval in intervals:
        if range_val / interval < num_ticks:
            return interval
    return intervals[-1]

def set_flexible_ticks(ax, lat, lon, num_ticks=5):
    """Set flexible ticks based on data range."""
    # Latitude
    min_val = np.min(lat)
    max_val = np.max(lat)
    interval = get_tick_interval((min_val, max_val), num_ticks)
    ticks = np.arange(np.floor(min_val / interval) * interval,
                      np.ceil(max_val / interval) * interval + interval,
                      interval)
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{abs(int(y))}°N' if y >= 0 else f'{abs(int(y))}°S'))

    # Longitude
    min_val = np.min(lon)
    max_val = np.max(lon)
    interval = get_tick_interval((min_val, max_val), num_ticks)
    ticks = np.arange(np.floor(min_val / interval) * interval,
                      np.ceil(max_val / interval) * interval + interval,
                      interval)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{abs(int(y))}°E' if y >= 0 else f'{abs(int(y))}°W'))


# Plotting the dbz data on a map
def plot_dbz_on_map(obs, bak, ana, lat, lon, output_filename):
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Define the color map and normalization
    cmap = create_custom_colormap()
    levels = np.arange(0, 75, 5)  # Intervals of 5 from 0 to 70
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # Define features for each subplot
    for ax, data, title in zip(axes, [obs, bak, ana], ['OBS', 'BAK', 'ANA']):
        # sort the values of data
        sort_indices = np.argsort(data)
        data_sort = data[sort_indices]
        lat_sort = lat[sort_indices]
        lon_sort = lon[sort_indices]

        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, edgecolor='black')
        ax.add_feature(cfeature.LAKES, alpha=0.5)
        ax.add_feature(cfeature.RIVERS) 
        ax.add_feature(cfeature.STATES, linestyle=':')
        sc = ax.scatter(lon_sort, lat_sort, c=data_sort, cmap=cmap, norm=norm, s=10, alpha=0.7, transform=ccrs.PlateCarree())
        ax.set_title(title)

        set_flexible_ticks(ax, lat, lon)

    # Shared colorbar
    cbar = fig.colorbar(sc, ax=axes, orientation='horizontal', fraction=0.05, pad=0.1)
    cbar.set_label('Radar reflectivity (dBZ)')
    cbar.set_ticks(levels)
    cbar.set_ticklabels(levels)

    plt.savefig(output_filename, bbox_inches='tight')

# Main script
filename = 'obsout_da_radar_dbz.h5'
station_name = "KRLX" # Set both station_name and radar_tilt to None for composite reflectivity graphic.
radar_tilt = 0.5
if station_name is not None:
   if radar_tilt is not None:
       output_filename = f'radar_dbz_{station_name}_tilt{radar_tilt}.png'
   else:
       output_filename = f'radar_dbz_{station_name}_alltilts.png'
else:
   if radar_tilt is not None:
       output_filename = f'all_radar_dbz_tilt{radar_tilt}.png'
   else:
       output_filename = 'composite_reflectivity.png'

obs, bak, ana, lat, lon = read_hdf5_file(filename, station_name, radar_tilt)

plot_dbz_on_map(obs, bak, ana, lat, lon, output_filename)
