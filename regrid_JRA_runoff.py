#!/usr/bin/env python

import argparse
import glob
import netCDF4
import numpy
from numba import jit
import os.path
import scipy.sparse

# Process command-line arguments
parser = argparse.ArgumentParser(description='''
    Re-grids runoff data from the JRA55-do dataset
    to an ocean model-grid of comparable or coarser resolution.''',
    epilog='Written by A.Adcroft, 2017.')
parser.add_argument('runoff_files', type=str, nargs='+',
    help='''Runoff files to process. Arguments will be wildcard expanded.''')
parser.add_argument('-sg','--super_grid', required=True,
    help='The ocean model supergrid file.')
parser.add_argument('-om','--ocean_mask', required=True,
    help='The ocean model mask file.')
parser.add_argument('-ra','--river_cell_area', required=True,
    help='The area of (river) cells for the gridded runoff data.')
parser.add_argument('-cd','--coastal_diagonals', action='store_true',
    help='Wet cells that touch a land cell diagonally are considered coastal cells.')
parser.add_argument('-wid','--write_intermediate_data', action='store_true',
    help='Write intermediate data to files.')
parser.add_argument('-q','--quiet', action='store_true',
    help='Turn off progress information.')
parser.add_argument('-d','--debug', action='store_true',
    help='Turn on debugging information.')
args = parser.parse_args()

if args.debug: print('args=',args)

# Read super-grid
if not args.quiet: print('Reading ocean grid data')
# (olon,olat) will be cell centers
olon = netCDF4.Dataset(args.super_grid).variables['x'][1::2,1::2]
olat = netCDF4.Dataset(args.super_grid).variables['y'][1::2,1::2]
oarea = netCDF4.Dataset(args.super_grid).variables['area'][:] # Returns 4 sub-cells
oarea = ( oarea[::2,::2] + oarea[1::2,1::2] ) + ( oarea[::2,1::2] + oarea[1::2,::2] ) # Area of model grid

# Read ocean mask
if not args.quiet: print('Reading ocean mask data')
vars = netCDF4.Dataset(args.ocean_mask).variables
omsk = None
for v in vars:
  if vars[v].shape==oarea.shape:
    if omsk is None:
      omsk = vars[v][:]
      break
    else:
      raise Exception('More than one variable in %s matches the expected shape'%(args.ocean_mask))
del vars, v
if omsk is None:    
  raise Exception('No variables found in %s with the expected shape'%(args.ocean_mask))

# Find coastal cells on model grid
def define_coastal(omsk):
  cmsk = 0 * omsk # All land should be 0
  cmsk[ (omsk>0) & (numpy.roll(omsk,1,axis=1)==0) ] = 1 # Land to the west
  cmsk[ (omsk>0) & (numpy.roll(omsk,-1,axis=1)==0) ] = 1 # Land to the east
  cmsk[ (omsk>0) & (numpy.roll(omsk,1,axis=0)==0) ] = 1 # Land to the south
  nom = numpy.roll(omsk,-1,axis=0) # Shift southward
  nom[-1,:] = omsk[-1,::-1] # Tri-polar fold
  cmsk[ (omsk>0) & (nom==0) ] = 1 # Land to the north
  if args.coastal_diagonals:
    cmsk[ (omsk>0) & (numpy.roll(numpy.roll(omsk,1,axis=0),1,axis=1)==0) ] = 1 # Land to the south-west
    cmsk[ (omsk>0) & (numpy.roll(numpy.roll(omsk,1,axis=0),-1,axis=1)==0) ] = 1 # Land to the south-east
    cmsk[ (omsk>0) & (numpy.roll(nom,1,axis=1)==0) ] = 1 # Land to the north-west
    cmsk[ (omsk>0) & (numpy.roll(nom,-1,axis=1)==0) ] = 1 # Land to the north-east
  return cmsk
if not args.quiet: print('Calculating coastal cells on model grid')
cmsk = define_coastal(omsk)

if args.write_intermediate_data:
  print('Write cmsk here')

# Ice-9 the nearest coastal cell id to all model cells
if not args.quiet: print('Calculating nearest coastal cells for all model cells')
def find_nearest_costal_cell(cmsk):
  oid = numpy.arange(cmsk.size,dtype=numpy.int32).reshape(cmsk.shape) # Ocean mode cell id
  ocid = oid * cmsk - (1 - cmsk)# Will have nearest oid of coastal cells (-1 for unassigned)
  ocidm = 1 * cmsk # Mask for assigned field
  while (ocid<0).sum()>0:
    # Look east
    difm = numpy.roll( ocidm, -1, axis=1) - ocidm
    ocid[ difm>0 ] = numpy.roll( ocid, -1, axis=1)[ difm>0 ]
    ocidm[ ocid>=0 ] = 1 # Flag all that have been assigned
    # Look west
    difm = numpy.roll( ocidm, 1, axis=1) - ocidm
    ocid[ difm>0 ] = numpy.roll( ocid, 1, axis=1)[ difm>0 ]
    ocidm[ ocid>=0 ] = 1 # Flag all that have been assigned
    # Look south
    difm = numpy.roll( ocidm, 1, axis=0) - ocidm
    difm[0,:] = 0 # Non-periodic across south
    ocid[ difm>0 ] = numpy.roll( ocid, 1, axis=0)[ difm>0 ]
    ocidm[ ocid>=0 ] = 1 # Flag all that have been assigned
    # Look north
    difm = numpy.roll( ocidm, -1, axis=0) - ocidm
    difm[-1,:] = 0 # THIS DOES NOT DO THE TRI-POLAR FOLD PROPERLY YET ***************
    ocid[ difm>0 ] = numpy.roll( ocid, -1, axis=0)[ difm>0 ]
    ocidm[ ocid>=0 ] = 1 # Flag all that have been assigned
  if args.debug: print('oid, ocid shapes:', oid.shape, ocid.shape)
  return oid,ocid
oid,ocid = find_nearest_costal_cell(cmsk)

if args.write_intermediate_data:
  print('Write ocid here')

# Read river grid
if not args.quiet: print('Reading river grid')
rlon = netCDF4.Dataset(args.river_cell_area).variables['longitude'][:]
rlat = netCDF4.Dataset(args.river_cell_area).variables['latitude'][:]
rarea = netCDF4.Dataset(args.river_cell_area).variables['areacello'][0]
if args.debug: print('rlon,rlat,rarea shapes:', rlon.shape, rlat.shape, rarea.shape)

# Connect river cells to model cells 
if not args.quiet: print('Calculating river-cell model-cell mapping (guess phase)')
roid = numpy.zeros((rlat.size,rlon.size),dtype=numpy.int32) - 1 # Will hold id of ocean cell overlying river grid
ir = (4*numpy.mod( olon, 360)).astype(int) # i-index in river-space of ocean cell
jr = (4*( olat+90.)).astype(int) # j-index in river-space of ocean cell
if args.debug: print('roid, ir, jr shapes:', roid.shape, ir.shape, jr.shape, jr.flatten().shape)
roid[ jr.flatten(), ir.flatten() ] = oid.flatten() # Associate river cell to ocean cell id
del ir, jr
# (note in the above there could be multiple ocean cells per river cell for which the result is ill-defined)
# The following fills holes that should really be filled with a nearest neighbor search...
while (roid<0).sum():
    roid[ roid<0 ] = numpy.roll(roid, 1, axis=1)[roid<0] # Assign from west
    roid[ roid<0 ] = numpy.roll(roid, -1, axis=1)[roid<0] # Assign from east
    tmp = numpy.roll(roid, 1, axis=0); tmp[0,:] = -1
    roid[ roid<0 ] = tmp[roid<0] # Assign from south
    tmp = numpy.roll(roid, -1, axis=0); tmp[-1,:] = -1
    roid[ roid<0 ] = tmp[roid<0] # Assign from north
del tmp

# Correction phase
@jit
def iter_ij(rx, ry, nio, njo, olonp, olatp, io, jo, it=None):
  if it is None: it=nio
  oi, oj = 0+io, 0+jo
  while it>0:
    it -= 1
    # d2 is the cost function we need to minimize
    # oi,oj are indices in model-grid space which in the padded arrays correspond to oi+1,oj+1
    dx = numpy.mod((rx - olonp[oj:oj+3,oi:oi+3]) + 540., 360.) - 180.
    d2 = dx**2 + (ry - olatp[oj:oj+3,oi:oi+3])**2
    # Move down gradient
    if d2[1,1]>d2[1,0]: oi -= 1
    elif d2[1,1]>d2[1,2]: oi += 1
    elif d2[1,1]>d2[0,1] and oj>1: oj -= 1
    elif d2[1,1]>d2[2,1]: oj += 1
    else: return oi,oj
    # Periodic in i
    if oi<=0: oi += nio
    elif oi>=nio: oi -= nio
    # Non-periodic in south, folded in north
    if oj<0: oj=0
    elif oj>=njo:
        oj,oi = njo-1,nio+1-oi
  raise Exception('Did not find minimum')
def correct_connections(rlon, rlat, rarea, olon, olat, omsk, oarea, roid):
  njr, nir = rarea.shape
  njo, nio = omsk.shape
  olatp = numpy.zeros((njo+2,nio+2))
  olatp[1:-1,1:-1] = olat
  olatp[0,:] = olatp[1,:] # Bottom edge
  olatp[-1,:] = olatp[-1,::-1] # Tri-polar fold
  olatp[:,0] = olatp[:,-2]; olatp[:,-1] = olatp[:,1] # Periodic in i
  olonp = numpy.zeros((njo+2,nio+2))
  olonp[1:-1,1:-1] = olon
  olonp[0,:] = olonp[1,:] # Bottom edge
  olonp[-1,:] = olonp[-1,::-1] # Tri-polar fold
  olonp[:,0] = olonp[:,-2] - 360. ; olonp[:,-1] = olonp[:,1] + 360. # Periodic in i
  if args.debug: print('olonp,olatp shapes:', olonp.shape, olatp.shape)
  ominlat = olat.min() # Used to avoid wasted interactions in Antarctica
  for rj in range(600,njr):
    ry = max( rlat[rj], ominlat )
    for ri in range(nir):
      rx = numpy.mod( rlon[ri] + 300., 360) - 300.
      oi = numpy.mod( roid[rj,ri], nir)
      oj = int( ( roid[rj,ri] - oi ) / nir )
      oi, oj = iter_ij( rx, ry, nio, njo, olonp, olatp, oi, oj )
      #if args.debug: print(oi,oj)
      roid[rj,ri] = oj * nio + oi
  return roid
if not args.quiet: print('Calculating river-cell model-cell mapping (correction phase)')
roid = correct_connections(rlon, rlat, rarea, olon, olat, omsk, oarea, roid)

# Construct sparse matrix
if not args.quiet: print('Constructing sparse matrix')
roid = roid.flatten()
rocid = ocid.flatten()[roid]
A = scipy.sparse.lil_matrix( (omsk.size, rarea.size), dtype=numpy.double )
A[rocid[numpy.arange(rarea.size,dtype=int)],numpy.arange(rarea.size,dtype=int)] = \
    rarea.flatten() #/ oarea.flatten()[roid[numpy.arange(rarea.size,dtype=int)]]
A = scipy.sparse.csr_matrix( A )

# Read each file and regrid data
def copy_meta(nc_in, nc_out, osize):
  for a in nc_in.ncattrs():
    v = nc_in.getncattr(a)
    if a == 'title': v = v + ', regridded on to OM4 0.25-degree grid'
    nc_out.setncattr(a, v)
  for d in nc_in.dimensions:
    dim = nc_in.dimensions[d]
    nm,sz = d,dim.size
    if dim.isunlimited(): sz = None
    if d == 'latitude': nm,sz = 'jh',osize[0]
    if d == 'longitude': nm,sz = 'ih',osize[1]
    nc_out.createDimension(nm, sz)
  for name in nc_in.variables:
    var = nc_in.variables[name]
    dims = list(var.dimensions)
    for n in range(len(dims)):
      if dims[n] == 'latitude': dims[n] = 'jh'
      if dims[n] == 'longitude': dims[n] = 'ih'
    nm = name
    if name == 'latitude': nm = 'jh'
    if name == 'longitude': nm = 'ih'
    n = nc_out.createVariable(nm, var.datatype, dimensions=(dims))
    for a in var.ncattrs():
      v = var.getncattr(a)
      if v == 'degrees_east': v = 'cell'
      if v == 'degrees_north': v = 'cell'
      if not a == 'point_spacing': n.setncattr(a, v)

def regrid(A, oarea, in_file, var_name='friver'):
  if args.debug: print('in_file =', in_file)
  out_file = os.path.splitext( os.path.basename(in_file) )[0] + '.%ix%i.nc'%(oarea.shape)
  nc_in = netCDF4.Dataset(in_file)
  nc_out = netCDF4.Dataset(out_file, 'w', format='NETCDF3_CLASSIC', clobber=True)
  copy_meta(nc_in, nc_out, (1080,1440))
  in_var = nc_in.variables[var_name]
  time = nc_out.variables['time']
  time_bnds = nc_out.variables['time_bnds']
  out_var = nc_out.variables[var_name]
  nc_out.variables['jh'][:] = numpy.arange(0.5, 0.5+omsk.shape[0] )
  nc_out.variables['ih'][:] = numpy.arange(0.5, 0.5+omsk.shape[1] )
  for n in range(in_var.shape[0]):
    rflux = in_var[n].astype(numpy.double)
    Rtot = (rarea * rflux ).sum()
    oflux = ( A * rflux.flatten() ).reshape( omsk.shape ) * omsk / oarea
    Otot = (oarea * oflux ).sum()
    if abs(Otot - Rtot) > 3e-15 * Rtot:
      print('record',n+1,'river total',Rtot,'kg s-1','ocean total',Otot,'kg s-1','error',(Otot - Rtot)/Rtot)
    out_var[n,:,:] = oflux[:,:]
    time[n] = nc_in.variables['time'][n]
    time_bnds[n,:] = nc_in.variables['time_bnds'][n,:]
  nc_in.close()
  nc_out.close()
for arg in args.runoff_files:
  for a_file in glob.glob(arg):
    if not args.quiet: print('Regridding', a_file)
    regrid(A, oarea, a_file)
