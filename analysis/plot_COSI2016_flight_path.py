import cartopy.crs as ccrs
import pandas as pd
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import matplotlib.cm as cm


### Import aspect data

# Path to csv
file = '/Volumes/JB/COSI/COSI2016_csvfiles/aspect_dso.csv'

# Load data. One row per second. Time in UTC, latitude and longitude in degrees
df = pd.read_csv(file, 
					header=0,
        			names=['time', 'latitude', 'longitude', 'altitude'], 
        			usecols=[0, 1, 2, 3])

# Remove rows which read 0 aspect info
# "2016-06-13 08:26:39" to "2016-06-13 11:43:26" have 0 values, among other times
df = df[df.latitude != 0.00]


### Use cartopy Orthographic projection as base image
projection = ccrs.Orthographic(
				central_longitude=-110.0, 
				central_latitude=-60.0, 
				globe=None)
geo = ccrs.Geodetic()

ax = plt.axes(projection=projection)
ax.stock_img()
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Plot location ~every hour (approximate because you removed all zero values)
cmap = cm.plasma

# Plot thick yellow line along the path (max color of cmap)
ax.plot(df.longitude[::3600], df.latitude[::3600],
         color=cmap(1.0),
         alpha=1,
         linewidth=3, 
         linestyle='-',
         transform=geo,
         zorder=0)

# Overlay colored points to indicate balloon altitude
plt.scatter(df.longitude[::3600], df.latitude[::3600],
				c=df.altitude[::3600]/33000.,
				cmap=cmap,
				marker='.',
				transform=geo,
				zorder=1)

# Plot thin black line over the whole path
ax.plot(df.longitude[::3600], df.latitude[::3600],
         color='k', 
         alpha=0.5,
         linewidth=1, 
         linestyle='-',
         transform=geo,
         zorder=2)

cbar = plt.colorbar()
cbar.set_label('Balloon altitude [33 km]')

plt.show()

