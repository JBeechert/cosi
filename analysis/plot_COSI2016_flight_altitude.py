import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Path to csv
file = '/Volumes/JB/COSI/COSI2016_csvfiles/COSI2016_FlightAltitude.csv'


# Load data. One row per second. Time in UTC, altitude in meters
df = pd.read_csv(file, skiprows=1, header=None, names=['time', 'altitude'])


# Remove rows which read 0 altitude
# "2016-06-13 08:26:39" to "2016-06-13 11:43:26" have 0 altitude, among others
df = df[df.altitude != 0.00]


# Make the figure
fig, ax = plt.subplots(1, 1) 

# Plot altitude [km] ~every hour (approximate because you removed non-zero values)
ax.plot(df['time'][::3600], df['altitude'][::3600]/1000, 'k-')
ax.set_xlabel('Time [UTC]', fontsize=14)
ax.set_ylabel('Altitude [km]', fontsize=14)

# Time labels: ~one label per week (approximate because you removed non-zero values)
times_per_week = np.array(df['time'][::604800])
x_labels = []
for x in times_per_week:
	x = x.split(' ')[0]
	x_labels.append(x)
ax.set_xticks(times_per_week)
ax.set_xticklabels(x_labels)
ax.tick_params(axis='both', size=14)

# Plot nominal 33 km flight altitude
plt.axhline(y=33, c='k', linestyle='--', label='33 km')

# Legend
plt.legend(loc='best')

plt.show()
