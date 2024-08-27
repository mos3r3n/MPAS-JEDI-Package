import os
import sys
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def plot_data_level(ax, xvalue, yvalue, label, diag_var, var_units, linecolor='black', linestyle='-', ):
    
    ax.plot(xvalue, yvalue, label=f'{label}', color=linecolor, linestyle=linestyle, linewidth=2.5)
    ax.set_ylabel('Model level', fontsize=18)
    ax.set_xlabel(f'RMS ({var_units})', fontsize=18)
    ax.set_title(f'{diag_var}', fontsize=18)
    ax.tick_params(axis='both', labelsize=18)
    max_value = np.float32(np.ceil(np.max(np.abs(xvalue))*2.0)/2.0)
    ax.set_xlim(0., max_value)
    ax.grid(True)

def plot_data_time(ax, times, data, label, diag_var, var_units, linecolor='black', linestyle='-', ):

    # Modify x-axis labels to include date and month
    x_labels = []
    x_ticks = []
    
    for itime in times:
        if (itime.day == 1 and itime.hour == 0) or itime == times[0]:
            x_labels.append(f'{itime.day}\n{itime.strftime("%b")}')
            x_ticks.append(itime)
            count_day = 1
        elif itime.hour == 0:
            x_labels.append(itime.day)
            x_ticks.append(itime)
            count_day += 1
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
    ax.set_ylabel(f'RMS ({var_units})', fontsize=22)
    ax.set_title(f'{diag_var}', fontsize=24)
    ax.tick_params(axis='both', labelsize=22)
    ax.grid(True)

def process_data(exps_name, exps_label, num_lev, grid_mesh, ensemble, variables, labels, units, ini_time, end_time, interval):

  time_format = "%Y%m%d%H"
  exp_dir = f"/glade/u/home/taosun/scratch/pandac"
  colors = ['black', 'blue', 'red']

  # System function to get directories
  times = []
  dirs = []
  run_time = ini_time
  while run_time <= end_time:
     date_str = run_time.strftime(time_format)
     dirs.append(date_str)
     times.append(run_time)
     run_time += interval

  num_exp = len(exps_name)
  num_var = len(variables)
  num_cycle = len(times)
  num_num = len(exps_name)
  
  if num_var == 1 or num_var == 2:
    num_rows = 1
    num_cols = num_var
  elif num_var == 3 or num_var == 4:
    num_rows = 2
    num_cols = 2
  elif num_var == 5 or num_var == 6:
    num_rows = 2
    num_cols = 3
  elif num_var == 7 or num_var == 8 or num_var == 9:
    num_rows = 3
    num_cols = 3
  else:
    num_rows = int(np.ceil(np.sqrt(num_var)))
    num_cols = int(np.ceil(num_var / num_rows))

  # Initialize arrays for each experiment
  rmse_bak = np.zeros((num_exp, num_cycle, num_lev, num_var))
  
  # Loop over cycles
  for iexp in range(num_exp):

      if grid_mesh[iexp] == '30km':
          ratio = 1
          cells = 655362
      elif grid_mesh[iexp] == '60-3km':
          ratio = 20
          cells = 835586

      for icycle in range(num_cycle):

          cycle_time = dirs[icycle]
          year  = cycle_time[:4]
          month = cycle_time[4:6]
          day   = cycle_time[6:8]
          hour  = cycle_time[8:]
          
          date_str = f"{year}-{month}-{day}_{hour}"

          # RMSE
          if ensemble[iexp]:
             filename1 = f"{exp_dir}/{exps_name[iexp]}/CyclingDA/{cycle_time}/bg/mean/mpasout.{date_str}.00.00.nc"
          else:
             filename1 = f"{exp_dir}/{exps_name[iexp]}/CyclingDA/{cycle_time}/bg/bg.{date_str}.00.00.nc"  

          filename2 = f"{exp_dir}/{exps_name[iexp]}/ExternalAnalyses/{grid_mesh[iexp]}/{cycle_time}/x{ratio}.{cells}.init.{date_str}.00.00.nc"
          print(filename1)

          with nc4.Dataset(filename1, 'r') as f1, nc4.Dataset(filename2, 'r') as f2:

             for ivar, diag_var in enumerate(variables):

                if diag_var == "U":
                    var_bak = f1['uReconstructZonal'][0,:,:]
                    var_gfs = f2['uReconstructZonal'][0,:,:]

                elif diag_var == "V": 
                    var_bak = f1['uReconstructMeridional'][0,:,:]
                    var_gfs = f2['uReconstructMeridional'][0,:,:]
                          
                elif diag_var == "Tv":
                    var_bak = f1['theta'][0,:,:]
                    var_gfs = f2['theta'][0,:,:]

                elif diag_var == "T":
                    var_p = f1['pressure_base'][0,:,:] + f1['pressure_p'][0,:,:] # pressure
                    var_t = f1['theta'][0,:,:]
                    var_bak = var_t * (var_p / 100000.)**0.286
                    var_p = f2['pressure_base'][0,:,:] + f2['pressure_p'][0,:,:] # pressure
                    var_t = f2['theta'][0,:,:]
                    var_gfs = var_t * (var_p / 100000.)**0.286

                elif diag_var == "Qv":
                    var_bak = f1['qv'][0,:,:] * 1000.  # kg/kg to g/kg
                    var_gfs = f2['qv'][0,:,:] * 1000.
          
                rmse_bak[iexp, icycle, :, ivar] = np.sqrt(np.nanmean((var_bak - var_gfs)**2.0, axis=0)) 

  # End of processing data 

  # Vertical profiles
  plot_name = f"verify_gfs_model_vert.pdf"
  fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 6, num_rows * 6))
    
  KLEV = np.arange(1, num_lev+1, 1)  
  rmse_bak_level = np.nanmean(rmse_bak, axis=1)
   
  for iexp in range(num_exp):

      exp_label = exps_label[iexp]

      for ivar, diag_var in enumerate(variables):

          # Determine the localtion of panels.
          row_position = ivar // num_cols
          col_position = ivar % num_cols

          if num_var == 1:
             current_ax = axes
          else:
             if num_rows > 1:
                current_ax = axes[row_position, col_position]
             else:
                current_ax = axes[col_position]

          # Plot RMSE
          plot_data_level(current_ax, rmse_bak_level[iexp, :, ivar], KLEV, f'RMSE: {exp_label}', labels[ivar], units[ivar], colors[iexp], '--')
      
          if ivar == num_var - 1:
             # Add legend to the last panel
             current_ax.legend(fontsize=15)

  plt.subplots_adjust(hspace=0.25, wspace=0.2) 

  for irow in range(num_rows):
     for icol in range(num_cols):
         if (irow * num_cols + icol) > (num_var - 1):
             axes[irow, icol].axis('off')

  # Save the plot
  plt.savefig(plot_name,dpi=600)
  
  ###########################################################################################  
  # Time series
  plot_name = f"verify_gfs_model_time.pdf"
  fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 12, num_rows * 9))

  rmse_bak_time = np.nanmean(rmse_bak[:,:,:,:], axis=2)

  # End of processing data  
  for iexp in range(num_exp):

      exp_label = exps_label[iexp]

      for ivar, diag_var in enumerate(variables):

          # Determine the localtion of panels.
          row_position = ivar // num_cols
          col_position = ivar % num_cols

          if num_var == 1:
             current_ax = axes
          else:
             if num_rows > 1:
                current_ax = axes[row_position, col_position]
             else:
                current_ax = axes[col_position]

          # Plot RMSE
          plot_data_time(current_ax, times, rmse_bak_time[iexp, :, ivar], f'RMSE: {exp_label}', labels[ivar], units[ivar], colors[iexp], '--')
      
          if ivar == num_var - 1:
             # Add legend to the last panel
             current_ax.legend(fontsize=15)

  plt.subplots_adjust(hspace=0.24, wspace=0.15) 

  for irow in range(num_rows):
     for icol in range(num_cols):
         if (irow * num_cols + icol) > (num_var - 1):
             axes[irow, icol].axis('off')

  # Save the plot
  plt.savefig(plot_name,dpi=600)

if __name__ == "__main__":
    # Get command-line arguments
    exps_name = ['taosun_hybrid_3denvar_O30kmIE60km_WarmStart_allsky_new','taosun_hybrid_3denvar_O30kmIE60km_WarmStart_allsky_noBCupdate']
    exps_label = ['HYB','HYB-noBCupdate', 'HYB-LOC']
    grid_mesh = ['30km','30km','30km']
    ensemble = [False, False, False]
    num_lev = 55
    variables = ['U','V','T','Qv']
    units = ['$\mathrm{m \cdot s^{-1}}$', '$\mathrm{m \cdot s^{-1}}$', '$\mathrm{K}$', '$\mathrm{g \cdot kg^{-1}}$']

    labels = ['(a) U','(b) V','(c) T','(d) Qv']

    ini_time = datetime(2023, 7, 7, 0)
    end_time = datetime(2023, 7, 20, 0)
    interval = timedelta(hours=6)

    process_data(exps_name, exps_label, num_lev, grid_mesh, ensemble, variables, labels, units, ini_time, end_time, interval)
