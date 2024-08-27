import h5py
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_obs_loc(obsout_dir, diag_obs, diag_var, invariant, if_regional=False, if_radiance=False):

   # Open the HDF5 file
   hdf5_file = f"{obsout_dir}/obsout_da_{diag_obs}.h5"
   print(f"processing {hdf5_file}")
   with h5py.File(hdf5_file, 'r') as f:
      # Extract QC values
      if not if_radiance:
         qc_values = f['EffectiveQC0'][diag_var][:]
      else:
         qc_values = f['EffectiveQC0']['brightnessTemperature'][:, diag_var-1]
      # Extract latitude and longitude from the MetaData group
      latitude = f['MetaData']['latitude'][:]
      longitude = f['MetaData']['longitude'][:]
    
      # Filter the data where QC values are equal to 0
      qc_zero_mask = (qc_values == 0)
      latitude_qc_zero = latitude[qc_zero_mask]
      longitude_qc_zero = longitude[qc_zero_mask]
      qc_nonzero_mask = (qc_values > 0)
      latitude_qc_nonzero = latitude[qc_nonzero_mask]
      longitude_qc_nonzero = longitude[qc_nonzero_mask]

   # Open MPAS file for terrain data
   ds = xr.open_dataset(invariant)

   # Convert longitude and latitude from radians to degrees for terrain data
   lonData = np.degrees(ds.lonCell)
   latData = np.degrees(ds.latCell)
   lonData = ((lonData + 180) % 360) - 180  # Adjust longitude to [-180, 180]

   # Prepare triangulation for terrain
   triang = tri.Triangulation(lonData, latData)

   # Extract terrain data
   terrain = ds["ter"]

   # Determine the map extent based on the observation data
   if if_regional:
      lat_min, lat_max = latData.min()-2.0, latData.max()+2.0
      lon_min, lon_max = lonData.min()-2.0, lonData.max()+2.0

   # Set up the plot with a map projection
   plt.figure(figsize=(12, 8))
   ax = plt.axes(projection=ccrs.PlateCarree())
   if if_regional:
      ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
   else:
      ax.set_global()
   ax.coastlines()
   ax.add_feature(cfeature.BORDERS, linestyle=':')
   ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
   ax.add_feature(cfeature.OCEAN)

   # Plot the terrain data using the triangulated mesh as the background
   contour = plt.tricontourf(triang, terrain, cmap='terrain', transform=ccrs.PlateCarree())
   # Add a colorbar that shows the actual terrain values
   cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05, shrink=0.7)
   cbar.set_label('Terrain Elevation (m)')

   # Overlay the scatter plot of observations on the terrain map
   plt.scatter(longitude_qc_nonzero, latitude_qc_nonzero, color='blue', marker='o', s=2, transform=ccrs.PlateCarree())
   plt.scatter(longitude_qc_zero, latitude_qc_zero, color='red', marker='o', s=2, transform=ccrs.PlateCarree())
   plt.title(f'{diag_obs} {diag_var} locations (blue: QC > 0, red: QC = 0)')

   # Save the figure with the map and scatter plot
   if not if_radiance:
      output_file = f'obsloc_{diag_obs}_{diag_var}_withTerrain.png'
   else:
      output_file = f'obsloc_{diag_obs}_channel{diag_var}_withTerrain.png'
   plt.savefig(output_file, dpi=300, bbox_inches='tight')


if __name__ == "__main__":

   sensors = ['satwind','sondes','aircraft','sfc','gnssrobndropp1d']
   #sensors = ['amsua_n15','amsua_n18','amsua_n19','amsua_metop-b','amsua_metop-c']
   obsout_dir = '/glade/u/home/taosun/scratch/pandac/taosun_3denvar-60-60-iter_amsua_metop-c_OR3km-SINGV/CyclingDA/2024061700/dbOut'
   invariant = '/glade/u/home/taosun/scratch/pandac/taosun_3denvar-60-60-iter_OR3km-SINGV/CyclingDA/2024061700/invariant.522172.nc'
   if_regional = True

   for diag_obs in sensors:

     if_radiance = False
     # Determine how specific diag_var to plot
     if diag_obs == "aircraft" or diag_obs == "sondes":
        variables = ['windNorthward', 'windEastward', 'airTemperature', 'specificHumidity']
     elif diag_obs == "satwind" or diag_obs == "satwnd":
        variables = ['windNorthward', 'windEastward']
     elif diag_obs == "sfc":
        variables = ['stationPressure']
     elif diag_obs == "gnssrobndropp1d":
        variables = ['bendingAngle']
     elif diag_obs == "amsua_metop-a":
        variables = [5, 6, 9]
        if_radiance = True
     elif diag_obs == "amsua_metop-b":
        variables = [8, 9]
        if_radiance = True
     elif diag_obs == "amsua_metop-c":
        variables = [5, 6, 9]
        if_radiance = True
     elif diag_obs == "amsua_n15":
        variables = [5, 7, 8, 9]
        if_radiance = True
     elif diag_obs == "amsua_n18":
        variables = [5, 6, 7, 8, 9]
        if_radiance = True
     elif diag_obs == "amsua_n19":
        variables = [5, 6, 7, 9]
        if_radiance = True
     elif diag_obs.startswith("mhs"):
        variables = [3, 4, 5]
        if_radiance = True

     for diag_var in variables:
        # plot conventional obs locations 
        plot_obs_loc(obsout_dir, diag_obs, diag_var, invariant, if_regional, if_radiance)
