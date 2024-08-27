import h5py
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

# Function to read HDF5 file and extract rw_obs data
def read_hdf5_file(filename, station_name=None, radar_tilt=None):
    with h5py.File(filename, 'r') as f:
        obs = f['ObsValue']['radialVelocity'][:]
        bak = f['hofx0']['radialVelocity'][:]
        ana = f['hofx3']['radialVelocity'][:]
        lat = f['MetaData']['latitude'][:]  # Assuming latitudes are stored here
        lon = f['MetaData']['longitude'][:]  # Assuming longitudes are stored here
        azi = f['MetaData']['radarAzimuth'][:]
        til = f['MetaData']['radarTilt'][:]
        nam = f['MetaData']['stationIdentification'][:].astype(str) 
    
        sin_til = f['MetaData']['sinTilt'][:]
        sinazi_costil =  f['MetaData']['sinAzimuthCosTilt'][:]
        sinazi_sintil =  f['MetaData']['sinAzimuthSinTilt'][:]


    nobs = len(sin_til)
    for iobs in range()

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
    colors = ['navy', 'blue', 'green', 'white', 'orange', 'red', 'purple']
    cmap_name = 'custom_cmap'
    return LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

# Plotting the dbz data on a map
def plot_dbz_on_map(obs, bak, ana, lat, lon, output_filename):
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Define the color map and normalization
    cmap = create_custom_colormap()
    levels = np.arange(-40, 45, 5)  # Intervals of 5 from 0 to 70
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # Define features for each subplot
    for ax, data, title in zip(axes, [obs, bak, ana], ['OBS', 'BAK', 'ANA']):
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, edgecolor='black')
        ax.add_feature(cfeature.LAKES, alpha=0.5)
        ax.add_feature(cfeature.RIVERS) 
        ax.add_feature(cfeature.STATES, linestyle=':')
        sc = ax.scatter(lon, lat, c=data, cmap=cmap, norm=norm, s=10, alpha=0.7, transform=ccrs.PlateCarree())
        ax.set_title(title)

        # Set latitude and longitude ticks
        ax.set_xticks(np.arange(min(lon), max(lon), step=1), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(min(lat), max(lat), step=1), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}'))

    # Shared colorbar
    cbar = fig.colorbar(sc, ax=axes, orientation='horizontal', fraction=0.05, pad=0.1)
    cbar.set_label('Composite radial velocity (m/s)')
    cbar.set_ticks(levels)
    cbar.set_ticklabels(levels)

    plt.suptitle('Composite Reflectivity Observations', fontsize=16)
    plt.savefig(output_filename, bbox_inches='tight')

# Main script
filename = 'obsout_da_radar_rw.h5'
output_filename = 'obs_rw.png'  # Specify the output file name
station_name = "KRLX"
radar_tilt = 0.5

obs, bak, ana, lat, lon = read_hdf5_file(filename, station_name, radar_tilt)

plot_dbz_on_map(obs, bak, ana, lat, lon, output_filename)
