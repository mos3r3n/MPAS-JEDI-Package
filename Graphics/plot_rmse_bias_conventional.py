import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def plot_data_level(ax, data, levels, label, diag_var, var_units, linecolor='black', linestyle='-', ):
    ax.plot(data, levels, label=f'{label}', color=linecolor, linestyle=linestyle, linewidth=2.5)

    #ax.set_ylim(levels[0], levels[-1])
    #ax.set_yticks(levels)
    if levels[1] > levels[0]:
      ax.set_ylabel('Height (km)', fontsize=18)
    else:
      ax.set_ylabel('Pressure (hPa)', fontsize=18)
      ax.invert_yaxis()

    ax.set_xlabel(f'Unit ({var_units})', fontsize=18)
    ax.set_title(f'{diag_var}', fontsize=20)
    ax.tick_params(axis='both', labelsize=18)
    ax.grid(True)

def plot_data_time(ax, times, data, label, diag_var, var_units, linecolor='black', linestyle='-', ):

    # Modify x-axis labels to include date and month
    x_labels = []
    x_ticks = []

    for itime in times:
        if (itime.day == 1 and itime.hour == 0) or itime == times[0]:
            #x_labels.append(f'{itime.month}-{itime.day}')
            x_labels.append(f'{itime.day}\n{itime.strftime("%b")}')
            x_ticks.append(itime)
            #count_day = 1
        elif itime.hour == 0:
            if (itime.day - times[0].day) % 2 ==0:
               x_labels.append(itime.day)
               x_ticks.append(itime)
            #count_day += 1
        elif itime == times[-1]:
            x_labels.append('')
            x_ticks.append(itime)

    # Set x-axis properties with a step size of 7 to skip every 7th label
    ticks_ini = times[0] - timedelta(hours=3)
    ticks_end = times[-1] + timedelta(hours=3)
    ax.set_xlim(ticks_ini, ticks_end)
    ax.set_xticks(x_ticks)  # Set ticks at indices of first day of each month
    ax.set_xticklabels(x_labels, ha='left')  # Use the same step size for labels

    ax.plot(times, data, label=f'{label}', color=linecolor, linestyle=linestyle, linewidth=2.5)
    ax.set_xlabel('Assimilation cycles', fontsize=22)
    ax.set_ylabel(f'Unit ({var_units})', fontsize=22)
    ax.set_title(f'{diag_var}', fontsize=24)
    ax.tick_params(axis='both', labelsize=22)
    ax.grid(True)
  
def process_data(exps_name, exps_label, diag_obs, ini_time, end_time, interval, outerloop, plot_rmse=True, plot_bias=True, plot_oma=False, apply_mask=False, mask_domain=None):

  # Determine how specific diag_var to plot
  if diag_obs == "aircraft" or diag_obs == "sondes":
     variables = ['windNorthward', 'windEastward', 'airTemperature', 'specificHumidity']
     units = ['$\mathrm{m \cdot s^{-1}}$', '$\mathrm{m \cdot s^{-1}}$', '$\mathrm{K}$', '$\mathrm{g \cdot kg^{-1}}$']
     labels = ['(a) U', '(b) V', '(c) T', '(d) Qv']
  elif diag_obs == "satwind" or diag_obs == "satwnd":
     variables = ['windNorthward', 'windEastward']
     units = ['$\mathrm{m \cdot s^{-1}}$', '$\mathrm{m \cdot s^{-1}}$']
     labels = ['(a) U', '(b) V']
  elif diag_obs == "sfc":
     variables = ['stationPressure']
     units = ['Pa']
     labels = ['PS']
  elif diag_obs == "gnssrobndropp1d":
     variables = ['bendingAngle']
     units = ['radian']
     labels = ['(a) Bending angle']

  time_format = "%Y%m%d%H"
  exp_dir = f"/glade/u/home/taosun/scratch/pandac"
  colors = ['blue', 'red', 'darkgreen', 'purple']
  if apply_mask:
    lat_min = mask_domain[0]
    lat_max = mask_domain[1]
    lon_min = mask_domain[2]
    lon_max = mask_domain[3]

  # System function to get directories
  times = []
  dirs = []
  run_time = ini_time
  while run_time <= end_time:
     date_str = run_time.strftime(time_format)
     dirs.append(date_str)
     times.append(run_time)
     run_time += interval

  cycle_num = len(times)
  exp_num = len(exps_name)
  var_num = len(variables)

  # Plotting
  IT = np.arange(1, cycle_num + 1)

  if var_num == 1 or var_num == 2:
    num_rows = 1
    num_cols = var_num
  elif var_num == 3 or var_num == 4:
    num_rows = 2
    num_cols = 2
  elif var_num == 5 or var_num == 6:
    num_rows = 2
    num_cols = 3
  elif var_num == 7 or var_num == 8 or var_num == 9:
    num_rows = 3
    num_cols = 3
  else:
    num_rows = int(np.ceil(np.sqrt(var_num)))
    num_cols = int(np.ceil(var_num / num_rows))

  # Vertical layers
  if diag_obs != "satwnd":
    vert_num = 30
    interval = 1.0
    KLEV = np.arange(1.0, vert_num*1.0+1.0, interval)
  else:
    KLEV = [1000, 925, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50]
    levels_half = [1050, 975, 850, 750, 650, 550, 450, 350, 275, 225, 175, 125, 75, 10]
    vert_num = len(KLEV)

  # Initialize arrays for each experiment  
  # For vertical profils
  if plot_rmse:
     rmse_omb_vert = np.zeros((exp_num, cycle_num, vert_num, var_num))
     if plot_oma:
        rmse_oma_vert = np.zeros((exp_num, cycle_num, vert_num, var_num))
  if plot_bias:
     mean_omb_vert = np.zeros((exp_num, cycle_num, vert_num, var_num))
     if plot_oma:
        mean_oma_vert = np.zeros((exp_num, cycle_num, vert_num, var_num))
  # For time series
  if plot_rmse:
     rmse_omb = np.zeros((exp_num, cycle_num, var_num))
     if plot_oma:
        rmse_oma = np.zeros((exp_num, cycle_num, var_num))
  if plot_bias:
     mean_omb = np.zeros((exp_num, cycle_num, var_num))
     if plot_oma:
        mean_oma = np.zeros((exp_num, cycle_num, var_num))
       
   # Loop over cycles
  for iexp in range(exp_num):

      for icycle in range(cycle_num):

          filename = f"{exp_dir}/{exps_name[iexp]}/CyclingDA/{dirs[icycle]}/dbOut/obsout_da_{diag_obs}.h5"
          print(filename)

          with h5py.File(filename, 'r') as f1:

              for ivar, diag_var in enumerate(variables):

                  # Meta data
                  if diag_obs != "satwnd":
                    hgt_omb = f1['MetaData']['height'][:] / 1000.
                    hgt_omb = np.where(hgt_omb > 0., hgt_omb, np.nan)
                  else:
                    hgt_omb = f1['MetaData']['pressure'][:] / 100.
                    hgt_omb = np.where(hgt_omb > 0., hgt_omb, np.nan)

                  lat_omb = f1['MetaData']['latitude'][:]
                  lon_omb = f1['MetaData']['longitude'][:]
                  lon_omb = np.where(lon_omb > 180., lon_omb-360., lon_omb)

                  # Assuming QC information is available in datasets named 'EffectiveQC0'
                  omb_qc = f1['EffectiveQC0'][diag_var][:]
                  if apply_mask:
                     omb_qc = np.where((lat_omb >= lat_min) & (lat_omb <= lat_max), omb_qc, np.nan)
                     omb_qc = np.where((lon_omb >= lon_min) & (lon_omb <= lon_max), omb_qc, np.nan)
                  if plot_oma:
                     oma_qc = f1[f'EffectiveQC{outerloop}'][diag_var][:]
                     if apply_mask:
                       oma_qc = np.where((lat_omb >= lat_min) & (lat_omb <= lat_max), oma_qc, np.nan)
                       oma_qc = np.where((lon_omb >= lon_min) & (lon_omb <= lon_max), oma_qc, np.nan)

                  ##### RMSE 
                  omb_var = f1['ombg'][diag_var][:]
                  if plot_oma:
                     oma_var = f1['oman'][diag_var][:]

                  # Filtering based on QC
                  omb_var = np.where(omb_qc == 0, omb_var, np.nan)
                  if plot_oma:
                     oma_var = np.where(oma_qc == 0, oma_var, np.nan)

                  # handing bending angeles
                  if (diag_var == "bendingAngle"):
                      omb_var = np.where(np.abs(omb_var) < 100., omb_var, np.nan)
                      if plot_oma:
                          oma_var = np.where(np.abs(oma_var) < 100., oma_var, np.nan)

                  if (diag_var == "specificHumidity"):
                      omb_var *= 1000.
                      if plot_oma:
                          oma_var *=1000.
                        
                  # Vertical profiles
                  for ilev in range(0, vert_num):

                     if diag_obs == "satwnd":
                        bot_lev = levels_half[ilev+1]
                        top_lev = levels_half[ilev]
                     else:
                        bot_lev = KLEV[ilev] - 0.5
                        top_lev = KLEV[ilev] + 0.5

                     omb_var1 = np.where((hgt_omb > bot_lev) & (hgt_omb <= top_lev), omb_var, np.nan)
                     if plot_oma:
                        oma_var1 = np.where((hgt_omb > bot_lev) & (hgt_omb <= top_lev), oma_var, np.nan)

                     if plot_rmse:
                        rmse_omb_vert[iexp, icycle, ilev, ivar] = np.sqrt(np.nanmean(omb_var1**2.0))
                        if plot_oma:
                           rmse_oma_vert[iexp, icycle, ilev, ivar] = np.sqrt(np.nanmean(oma_var1**2.0))

                     if plot_bias:
                        mean_omb_vert[iexp, icycle, ilev, ivar] = np.nanmean(omb_var1)
                        if plot_oma:
                           mean_oma_vert[iexp, icycle, ilev, ivar] = np.nanmean(oma_var1)

                  # Total level
                  if plot_rmse:
                     rmse_omb[iexp, icycle, ivar] = np.sqrt(np.nanmean(omb_var**2.0))
                     if plot_oma:
                        rmse_oma[iexp, icycle, ivar] = np.sqrt(np.nanmean(oma_var**2.0))
                  if plot_bias:
                     mean_omb[iexp, icycle, ivar] = np.nanmean(omb_var)
                     if plot_oma:
                        mean_oma[iexp, icycle, ivar] = np.nanmean(oma_var)

                  # Cleanup
                  del omb_var, omb_qc, omb_var1
                  if plot_oma:
                     del oma_var, oma_qc, oma_var1

# Vertical profiles
  if apply_mask:
    plot_name = f"global_verify_{diag_obs}_vert_height.pdf"
  else:
    plot_name = f"regional_verify_{diag_obs}_vert_height.pdf"

  fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 6, num_rows * 6))

  if plot_rmse:
     rmse_omb_level = np.nanmean(rmse_omb_vert, axis=1)
     if plot_oma:
        rmse_oma_level = np.nanmean(rmse_oma_vert, axis=1)
  if plot_bias:
     mean_omb_level = np.nanmean(mean_omb_vert, axis=1)
     if plot_oma:
        mean_oma_level = np.nanmean(mean_oma_vert, axis=1)

  for iexp in range(exp_num):

      exp_label = exps_label[iexp]

      for ivar, diag_var in enumerate(variables):

          # Determine the localtion of panels.
          row_position = ivar // num_cols
          col_position = ivar % num_cols

          if var_num == 1:
             current_ax = axes
          else:
             if num_rows > 1:
                current_ax = axes[row_position, col_position]
             else:
                current_ax = axes[col_position]

          if plot_rmse:
              if plot_oma:
                  plot_data_level(current_ax, rmse_omb_level[iexp, :, ivar], KLEV, f'{exp_label}: RMSE - OMB', labels[ivar], units[ivar], colors[iexp], '-')
                  plot_data_level(current_ax, rmse_oma_level[iexp, :, ivar], KLEV, f'{exp_label}: RMSE - OMA', labels[ivar], units[ivar], colors[iexp], '--')
              else:
                  plot_data_level(current_ax, rmse_omb_level[iexp, :, ivar], KLEV, f'RMSE: {exp_label}', labels[ivar], units[ivar], colors[iexp], '-')
          if plot_bias:
              if plot_oma:
                  plot_data_level(current_ax, mean_omb_level[iexp, :, ivar], KLEV, f'{exp_label}: MEAN - OMB', labels[ivar], units[ivar], colors[iexp], '-')
                  plot_data_level(current_ax, mean_oma_level[iexp, :, ivar], KLEV, f'{exp_label}: MEAN - OMA', labels[ivar], units[ivar], colors[iexp], '--')
              else:
                  plot_data_level(current_ax, mean_omb_level[iexp, :, ivar], KLEV, f'BIAS: {exp_label}', labels[ivar], units[ivar], colors[iexp], '-.')

          if ivar == var_num - 1:
             # Add legend to the last panel
             current_ax.legend(fontsize=15)
            
  for irow in range(num_rows):
     for icol in range(num_cols):
         if (irow * num_cols + icol) > (var_num - 1):
             axes[irow, icol].axis('off')

  plt.subplots_adjust(hspace=0.30, wspace=0.30)

  plt.subplots_adjust(bottom=0.15, top=0.85, left=0.15, right=0.85)

  # Save the plot
  plt.savefig(plot_name,dpi=600)

  ###########################################################################################  
  # Time series
  if apply_mask:
    plot_name = f"regional_{diag_obs}_time_series.pdf"
  else:
    plot_name = f"global_{diag_obs}_time_series.pdf"
  fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 13, num_rows * 9))

  # End of processing data  
  for iexp in range(exp_num):

      exp_label = exps_label[iexp]

      for ivar, diag_var in enumerate(variables):

          # Determine the localtion of panels.
          row_position = ivar // num_cols
          col_position = ivar % num_cols

          if var_num == 1:
             current_ax = axes
          else:
             if num_rows > 1:
                current_ax = axes[row_position, col_position]
             else:
                current_ax = axes[col_position]

          if plot_rmse:
              if plot_oma:
                  plot_data_time(current_ax, times, rmse_omb[iexp, :, ivar], f'{exp_label}: RMSE - OMB', labels[ivar], units[ivar], colors[iexp], '-')
                  plot_data_time(current_ax, times, rmse_oma[iexp, :, ivar], f'{exp_label}: RMSE - OMA', labels[ivar], units[ivar], colors[iexp], '--')
              else:
                  plot_data_time(current_ax, times, rmse_omb[iexp, :, ivar], f'RMSE: {exp_label}', labels[ivar], units[ivar], colors[iexp], '-')
          if plot_bias:
              if plot_oma:
                  plot_data_time(current_ax, times, mean_omb[iexp, :, ivar], f'{exp_label}: MEAN - OMB', labels[ivar], units[ivar], colors[iexp], '-')
                  plot_data_time(current_ax, times, mean_oma[iexp, :, ivar], f'{exp_label}: MEAN - OMA', labels[ivar], units[ivar], colors[iexp], '--')
              else:
                  plot_data_time(current_ax, times, mean_omb[iexp, :, ivar], f'BIAS: {exp_label}', labels[ivar], units[ivar], colors[iexp], '-.')

          if ivar == var_num - 1:
             # Add legend to the last panel
             current_ax.legend(fontsize=15)

  plt.subplots_adjust(hspace=0.30, wspace=0.25)

  plt.subplots_adjust(bottom=0.2, top=0.8, left=0.2, right=0.8)

  for irow in range(num_rows):
     for icol in range(num_cols):
         if (irow * num_cols + icol) > (var_num - 1):
             axes[irow, icol].axis('off')

  # Save the plot
  plt.savefig(plot_name,dpi=600)

if __name__ == "__main__":
    # Get command-line arguments
    #exps_name = ['taosun_3dhybrid_O30kmIE60km_USGS_MPAS-modelv8.2.1','taosun_3dhybrid_O60-3kmIE60km_TEST_Modelv8.2.1_samePP']
    #exps_label = ['30km','60-3km']
    exps_name = ['taosun_3dhybrid_O30kmIE60km_AllSky_MPAS-modelv8.2.1','taosun_3dhybrid_O60-3kmIE60km_AllSky_Modelv8.2.1_samePP']
    exps_label = ['AllSky_30km','AllSky_60-3km']
    sensors = ['satwind','satwind','sondes','aircraft','sfc','gnssrobndropp1d']
    plot_rmse = True
    plot_bias = False
    plot_oma = False
    apply_mask = False
    mask_domain = [25., 50., -125., -65.]
    outerloop = 2

    ini_time = datetime(2018, 4, 17, 0)
    end_time = datetime(2018, 5, 10, 0)

    for diag_obs in sensors:
        if diag_obs == 'sondes':
           interval = timedelta(hours=12)
        else:
           interval = timedelta(hours=6)
        process_data(exps_name, exps_label, diag_obs, ini_time, end_time, interval, outerloop, plot_rmse, plot_bias, plot_oma, apply_mask, mask_domain)
                                                                                                                                                          
