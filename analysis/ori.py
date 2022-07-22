# Written by T. Siegert to create a custom orientation file with proper instrument coordiantes.
# Example: write a file which moves the instrument across the Galactic Plane at specified time and 
#  lon/lat intervals.

import numpy as np
from numpy import sin,cos,tan,arctan2,arctan,arccos,arcsin,pi
from scipy.spatial.transform import Rotation as R


def sph2car(l, b):
	# make cartesian coordinates from spherical
	return(np.array([cos(np.deg2rad(b))*cos(np.deg2rad(l)), cos(np.deg2rad(b))*sin(np.deg2rad(l)), sin(np.deg2rad(b))]))

def car2sph(vec):
	# make spherical from cartesian

	# ensure arcsin is well-defined
	if vec[2] > 1.0:
		vec[2] = 1.0
	elif vec[2] < -1.0:
		vec[2] = -1.0

	return(np.rad2deg(np.array([arctan2(vec[1], vec[0]), arcsin(vec[2])])))

def rotation(l, b):
    # initial coordinates:
    # x pointing toward galactic centre: l_x = 0, b_x = 0 (right hand middle finger toward screen)
    # z pointing toward the galactic north pole, l_z = 0, b_z = 90 (right hand thumb up)
    # y pointing toward the left, l_y = 90, b_y = 0 (right hand index to the left)
	l_x, b_x = 0., 0.
	l_z, b_z = 0, 90.
	l_y, b_y = 90., 0.

    # making unit vectors out of that with spherical coordinates
	z = sph2car(l_z, b_z)
	x = sph2car(l_x, b_x)
	y = sph2car(l_y, b_y)

    # now x is pointing toward the centre...

    # rotation about y (index) gives latitude change
    # rotation about z (middle) gives longitude change
    # need -b intead of b because of left-handed Galactic coords? weird but whatever
	r = R.from_euler('yz', [-b, l], degrees=True)

	x_new, y_new, z_new = r.apply(x), r.apply(y), r.apply(z)

	l_z_new, b_z_new = car2sph(z_new)
	l_x_new, b_x_new = car2sph(x_new)
	l_y_new, b_y_new = car2sph(y_new)

	return(np.array([l_x_new, b_x_new, l_y_new, b_y_new, l_z_new, b_z_new]))

 

# Swap x and z in the above for the ori file since
# Type OrientationsGalactic
# OG 0 90 0 0 0
# gives you one pointing at 0 0 with some orientation of the instrument.
# It’s (time) gb_x gl_x gb_z gl_z, where
# z is the optical axis (zenith direction) and
# x is one of the other axes of COSI. From the right-hand rule, y follows automatically which is why you don't need y.
# So, 90 0 0 0 points with x to the Galactic north pole, z to the Galactic center, and y “to the right”, i.e. at b=0, l=90.
# In rotation(l, b) function, the inital x and z coordinates are swapped from the ori definition.

lines = []

min_time = 0 # s
sim_time = 10000
Delta_T = 100
times = np.arange(min_time, sim_time+Delta_T, Delta_T)
num_steps = len(times)

min_lon = -100 # deg
max_lon = 100 

# zenith latitude = 0, zenith longitude moves along the Galactic Plane from min to max lon zenith
lons = np.linspace(min_lon, max_lon, int(np.floor(num_steps/2)))

for i in range(0, int(np.floor(num_steps/2))):
	time = times[i]

	l = lons[i]
	b = 0

	#print('time, l, b: ', time, l, b, '\n')
	l_x, b_x, l_y, b_y, l_z, b_z = rotation(l, b)

	# swap x and z
	gb_x = b_z
	gl_x = l_z
	gb_z = b_x
	gl_z = l_x

	line = f'OG {time} {gb_x} {gl_x} {gb_z} {gl_z}'
	lines.append(line)


# zenith latitude = 15, zenith longitude moves along the Galactic Plane from max to min
#  lon zenith 
lons_2 = np.linspace(max_lon, min_lon, int(np.ceil(num_steps/2)))

for i in range(int(np.floor(num_steps/2)), num_steps):
	time = times[i]

	l = lons_2[i-int(np.floor(num_steps/2))]
	b = 15

	#print('time, l, b: ', time, l, b, '\n')
	l_x, b_x, l_y, b_y, l_z, b_z = rotation(l, b)

	# swap x and z
	gb_x = b_z
	gl_x = l_z
	gb_z = b_x
	gl_z = l_x

	line = f'OG {time} {gb_x} {gl_x} {gb_z} {gl_z}'
	lines.append(line)


with open('GalCenter.ori', 'w') as f:
	f.write('Type OrientationsGalactic \n')
	f.writelines('\n'.join(lines))
	f.write('\nEN')
f.close()

