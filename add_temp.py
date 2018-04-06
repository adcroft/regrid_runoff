import re
import numpy as np
import netCDF4
import sys
import pdb
from scipy.interpolate import interp1d
import pyroms

outfile = sys.argv[1]

# Read the river temperatures
f = open('KenaiRiverTemps.dat', 'r')
# Eat first two lines
f.readline()
f.readline()

# These are for the ROMS sources file
ttime = []
temp = []
salt = []

#pdb.set_trace()

for line in f:
    nul, a, b, c = re.split('\s+', line)
    ttime.append(float(a))
    temp.append(float(b))
    salt.append(0.0)

grd = pyroms.grid.get_ROMS_grid('ARCTIC4')
lat_rho = grd.hgrid.lat_rho

time2 = [0, 46, 76, 107, 137, 167, 198, 228, 259, 289, 320, 365]
t2 = [0.0, 0.0, 0.0, 0.0, 2.0, 10.0, 15.0, 13.0, 5.0, 2.0, 0.0, 0.0]
f2 = interp1d(time2, t2, kind='cubic')
temp2 = f2(ttime)
for t in range(len(temp2)):
    temp2[t] = max(temp2[t], 0.0)

# create file with all the objects
out = netCDF4.Dataset(outfile, 'a', format='NETCDF3_64BIT')
xi = out.variables['river_Xposition'][:]
eta = out.variables['river_Eposition'][:]
Nr = xi.shape[0]
river_tr = np.zeros((len(ttime), Nr))

for k in range(Nr):
    lat = lat_rho[eta[k],xi[k]]
    if lat > 65.0:
        river_tr[:,k] = temp2
    else:
        river_tr[:,k] = temp


out.createDimension('river_tracer_time', len(ttime))

times = out.createVariable('river_tracer_time', 'f8', ('river_tracer_time'))
times.units = 'day'
times.cycle_length = 365.25
times.long_name = 'river tracer time'

temp = out.createVariable('river_temp', 'f8', ('river_tracer_time', 'river'))
temp.long_name = 'river runoff potential temperature'
temp.units = 'Celsius'
temp.time = 'river_tracer_time'

salt = out.createVariable('river_salt', 'f8', ('river_tracer_time', 'river'))
salt.long_name = 'river runoff salinity'
salt.time = 'river_tracer_time'

out.variables['river_tracer_time'][:] = ttime
out.variables['river_temp'][:] = river_tr
river_tr = np.zeros((len(ttime), Nr))
out.variables['river_salt'][:] = river_tr

out.close()
