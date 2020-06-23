from subprocess import call
import xarray
import sys


year = sys.argv[1]
# Open output from regrid_runoff
# Fill any missing/nan values with zero.
base = xarray.open_dataset('Arctic5_runoff_'+year+'.nc').fillna(0)

# Choose where to save the results from this script
output_file = 'Arctic5_runoff_'+year+'_a.nc'

# Save file, making sure that time is the unlimited dimension.
compress = True
if compress:
   fileformat, zlib, complevel, area_dtype = 'NETCDF4', True, 1, 'f4'
else:
   fileformat, zlib, complevel, area_dtype = 'NETCDF3_64BIT_OFFSET', None, None, 'd'

comp = dict(zlib=zlib, complevel=complevel)
encoding = {var: comp for var in base.data_vars}
base.to_netcdf(
    output_file,
    unlimited_dims=['time'],
    format=fileformat,
    encoding=encoding
)

# Delete the _FillValue attribute.
# MOM will crash if this attribute is present.
# I can't find how to get xarray to not write it.
# Need to have NCO in your path or module load nco before running this script.
# Get rid of all of them, or just that of friver:
call(f'ncatted -h -O -a _FillValue,,d,, {output_file}', shell=True)
#call(f'ncatted -h -O -a _FillValue,friver,d,, {output_file}', shell=True)
