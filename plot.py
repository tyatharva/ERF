import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from util.visualizer import safe_fss


def plot_precipitation(ax, ds, precip_var, title, elevation_data, letter, intensity, cmap='turbo'):
    """
    ax: Matplotlib axis
    ds: xarray Dataset
    precip_var: Precipitation variable name
    title: Title for the subplot
    elevation_data: Elevation data to overlay
    cmap: Colormap for precipitation
    show_grid: Whether to show lat/lon grid
    letter: Letter to place in the top left corner of the plot
    """
    # Extract data
    precip = ds[precip_var].values
    actual = ds['Target'].values
    lon = ds['lon'].values
    lat = ds['lat'].values

    fss = "{:.3f}".format(0)
    if intensity == 'light': fss = "{:.3f}".format(safe_fss(precip, actual, 0.25, 5))
    elif intensity =='moderate': fss = "{:.3f}".format(safe_fss(precip, actual, 1, 5))
    elif intensity == 'heavy': fss = "{:.3f}".format(safe_fss(precip, actual, 2, 5))

    # Plot data
    c = ax.contourf(lon, lat, precip, levels=np.linspace(0.25, 6.25, 25), cmap=cmap, extend='max', alpha=0.9, zorder=3)
    ax.contourf(lon, lat, elevation_data, levels=np.linspace(-750, 2000, 23), cmap='terrain', alpha=0.5, zorder=1)
    ax.add_feature(cfeature.STATES, edgecolor='black', linewidth=0.7, zorder=2)
    ax.add_feature(cfeature.LAKES, alpha=0.7, zorder=2)
    ax.text(0.05, 0.95, letter, transform=ax.transAxes, fontsize=16, fontweight='bold', color='black', va='top', ha='left')
    ax.text(0.85, 0.95, fss, transform=ax.transAxes, fontsize=16, fontweight='bold', color='black', va='top', ha='left')

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, xpadding=-2, ypadding=-1.5, linewidth=0.5, color='gray', alpha=0.7)
    gl.top_labels = True
    gl.left_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.xformatter = LongitudeFormatter(number_format='.1f')
    gl.yformatter = LatitudeFormatter(number_format='.1f')
    gl.xlabel_style = {'size': 12}
    gl.ylabel_style = {'size': 12}

    return c


# Load data
l1 = 'results/0.001_wl1_25_18.0_9_e/exp_latest/20170129_09e.nc'
l2 = 'results/0.001_wl1_25_18.0_9_e/exp_latest/20250103_15e.nc'
l3 = 'results/0.001_wl1_25_18.0_9_e/exp_latest/20221117_06e.nc'

col1 = xr.open_dataset('results/0.001_l1_25_18.0_9_e/exp_latest/20170129_09e.nc')
col2 = xr.open_dataset('results/0.001_l1_25_18.0_9_e/exp_latest/20250103_15e.nc')
col3 = xr.open_dataset('results/0.001_l1_25_18.0_9_e/exp_latest/20221117_06e.nc')


main1 = xr.open_dataset(l1)
main2 = xr.open_dataset(l2)
main3 = xr.open_dataset(l3)


# Load elevation
elev_ds = xr.open_dataset('../dem/dem_e.nc')
elev_subset = elev_ds['elev'].values


# Set up figure and axes
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(25, 20), subplot_kw={'projection': ccrs.PlateCarree()})
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.rc('pdf', fonttype=42)


# Create labels
# Columns are timestamps
cols = [
    datetime.strptime(os.path.splitext(os.path.basename(l1))[0][:-1], '%Y%m%d_%H').strftime('%H%M UTC\n%d %b. %Y') + ' (Light)',
    datetime.strptime(os.path.splitext(os.path.basename(l2))[0][:-1], '%Y%m%d_%H').strftime('%H%M UTC\n%d %b. %Y') + ' (Moderate)',
    datetime.strptime(os.path.splitext(os.path.basename(l3))[0][:-1], '%Y%m%d_%H').strftime('%H%M UTC\n%d %b. %Y') + ' (Heavy)',
]

# Rows are data sources
rows = ['HRRR', 'LESNet-A', 'LESNet-B', 'Target']

# Annotate top columns with dates
for ax, col in zip(axes[0], cols):
    ax.annotate(col, xy=(0.5, 1.05), xycoords='axes fraction',
                ha='center', va='bottom', fontsize=16)

# Annotate left rows with source labels
for ax, row in zip(axes[:, 0], rows):
    ax.annotate(row, xy=(-0.01, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                xycoords='axes fraction', textcoords='offset points',
                size=16, ha='right', va='center', rotation=90)



# Plot data
# Column 0: l1
plot_precipitation(axes[0, 0], main1, 'HRRR', '', elev_subset, 'a', 'light')
plot_precipitation(axes[1, 0], main1, 'Model', '', elev_subset, 'd', 'light')
plot_precipitation(axes[2, 0], col1, 'Model', '', elev_subset, 'g', 'light')
plot_precipitation(axes[3, 0], main1, 'Target', '', elev_subset, 'j', 'light')

# Column 1: l2
plot_precipitation(axes[0, 1], main2, 'HRRR', '', elev_subset, 'b', 'moderate')
plot_precipitation(axes[1, 1], main2, 'Model', '', elev_subset, 'e', 'moderate')
plot_precipitation(axes[2, 1], col2, 'Model', '', elev_subset, 'h', 'moderate')
plot_precipitation(axes[3, 1], main2, 'Target', '', elev_subset, 'k', 'moderate')

# Column 2: l3
plot_precipitation(axes[0, 2], main3, 'HRRR', '', elev_subset, 'c', 'heavy')
plot_precipitation(axes[1, 2], main3, 'Model', '', elev_subset, 'f', 'heavy')
plot_precipitation(axes[2, 2], col3, 'Model', '', elev_subset, 'i', 'heavy')
c1 = plot_precipitation(axes[3, 2], main3, 'Target', '', elev_subset, 'l', 'heavy')



# Create colorbar
cbar = fig.colorbar(c1, ax=axes, orientation='horizontal', fraction=0.05, pad=0.025)
cbar.set_label('Next-hour precipitation (mm)', fontsize=16) # (mm h$^{-1}$)
cbar.ax.tick_params(labelsize=14)
# cbar.set_ticks(np.arange(0.5, 5.0, 0.5))


# Save and display
# plt.savefig('e.pdf', dpi=600, bbox_inches='tight')
plt.show()
