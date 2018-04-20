#!/bin/env python
# Add bio tracers to a rivers file.
# This one is for COBALT.

import re
import numpy as np
import netCDF4
import sys
import pdb

outfile = sys.argv[1]

# Read the river temperatures
f = open('NGOA_River.dat', 'r')
# Eat first two lines
f.readline()
f.readline()

# These are for the ROMS sources file
ttime = []
alk = []
dic = []
doc = []
sdon = []
ldon = []
din = []
do = []
po4 = []
dsi = []

# Converting from mmol/m^3 to mol/kg
scale = 1.0e-6
#pdb.set_trace()

for line in f:
    a, b, c, d, e, f, g, h, i, j, k = re.split('\s+', line)
    ttime.append(float(a))
    alk.append(scale * float(b))
    dic.append(scale * float(c))
    doc.append(scale * float(d))
    sdon.append(scale * float(e))
    ldon.append(scale * float(f))
    din.append(scale * float(g))
    do.append(scale * float(h))
    po4.append(scale * float(i))
    dsi.append(scale * float(j))

# create file with all the objects
out = netCDF4.Dataset(outfile, 'w', format='NETCDF3_64BIT')
jra = netCDF4.Dataset('/center1/AKWATERS/kshedstrom/Arctic2/Files/JRA_Arctic_rivers_2011.nc', 'r')
rivers = jra.variables['river'][:]

new = True
#new = False
if new:
    out.createDimension('river_bgc_time', len(ttime))
    out.createDimension('river', len(rivers))

    times = out.createVariable('river_bgc_time', 'f8', ('river_bgc_time'))
    times.units = 'day'
    times.cycle_length = 365.25
    times.long_name = 'river bgc time'

    alky = out.createVariable('river_alk', 'f8', ('river_bgc_time', 'river'))
    alky.long_name = 'river runoff alkalinity'
    alky.units = 'mol/kg'
    alky.time = 'river_bgc_time'

    cadet_arag = out.createVariable('river_cadet_arag', 'f8', ('river_bgc_time', 'river'))
    cadet_arag.long_name = 'river runoff detrital CaCO3'
    cadet_arag.units = 'mol/kg'
    cadet_arag.time = 'river_bgc_time'

    cadet_calc = out.createVariable('river_cadet_calc', 'f8', ('river_bgc_time', 'river'))
    cadet_calc.long_name = 'river runoff detrital CaCO3'
    cadet_calc.units = 'mol/kg'
    cadet_calc.time = 'river_bgc_time'

    rdic = out.createVariable('river_dic', 'f8', ('river_bgc_time', 'river'))
    rdic.long_name = 'river runoff dissolved inorganic carbon'
    rdic.units = 'mol/kg'
    rdic.time = 'river_bgc_time'

    fed = out.createVariable('river_fed', 'f8', ('river_bgc_time', 'river'))
    fed.long_name = 'river runoff dissolved iron'
    fed.units = 'mol/kg'
    fed.time = 'river_bgc_time'

    fedet = out.createVariable('river_fedet', 'f8', ('river_bgc_time', 'river'))
    fedet.long_name = 'river runoff detrital iron'
    fedet.units = 'mol/kg'
    fedet.time = 'river_bgc_time'

    fedi = out.createVariable('river_fedi', 'f8', ('river_bgc_time', 'river'))
    fedi.long_name = 'river runoff diazotroph iron'
    fedi.units = 'mol/kg'
    fedi.time = 'river_bgc_time'

    felg = out.createVariable('river_felg', 'f8', ('river_bgc_time', 'river'))
    felg.long_name = 'river runoff large phytoplankton iron'
    felg.units = 'mol/kg'
    felg.time = 'river_bgc_time'

    femd = out.createVariable('river_femd', 'f8', ('river_bgc_time', 'river'))
    femd.long_name = 'river runoff medium phytoplankton iron'
    femd.units = 'mol/kg'
    femd.time = 'river_bgc_time'

    fesm = out.createVariable('river_fesm', 'f8', ('river_bgc_time', 'river'))
    fesm.long_name = 'river runoff small phytoplankton iron'
    fesm.units = 'mol/kg'
    fesm.time = 'river_bgc_time'

    ldon = out.createVariable('river_ldon', 'f8', ('river_bgc_time', 'river'))
    ldon.long_name = 'river runoff labile DON'
    ldon.units = 'mol/kg'
    ldon.time = 'river_bgc_time'

    ldop = out.createVariable('river_ldop', 'f8', ('river_bgc_time', 'river'))
    ldop.long_name = 'river runoff labile DOP'
    ldop.units = 'mol/kg'
    ldop.time = 'river_bgc_time'

    lith = out.createVariable('river_lith', 'f8', ('river_bgc_time', 'river'))
    lith.long_name = 'river runoff lithogenic aluminosilicate'
    lith.units = 'mol/kg'
    lith.time = 'river_bgc_time'

    lithdet = out.createVariable('river_lithdet', 'f8', ('river_bgc_time', 'river'))
    lithdet.long_name = 'river runoff lithdet'
    lithdet.units = 'mol/kg'
    lithdet.time = 'river_bgc_time'

    nbact = out.createVariable('river_nbact', 'f8', ('river_bgc_time', 'river'))
    nbact.long_name = 'river runoff bacterial N'
    nbact.units = 'mol/kg'
    nbact.time = 'river_bgc_time'

    ndet = out.createVariable('river_ndet', 'f8', ('river_bgc_time', 'river'))
    ndet.long_name = 'river runoff detrital N'
    ndet.units = 'mol/kg/m3'
    ndet.time = 'river_bgc_time'

    ndi = out.createVariable('river_ndi', 'f8', ('river_bgc_time', 'river'))
    ndi.long_name = 'river runoff diazotroph N'
    ndi.units = 'mol/kg'
    ndi.time = 'river_bgc_time'

    nlg = out.createVariable('river_nlg', 'f8', ('river_bgc_time', 'river'))
    nlg.long_name = 'river runoff large phytoplankton N'
    nlg.units = 'mol/kg'
    nlg.time = 'river_bgc_time'

    nmd = out.createVariable('river_nmd', 'f8', ('river_bgc_time', 'river'))
    nmd.long_name = 'river runoff medium phytoplankton N'
    nmd.units = 'mol/kg'
    nmd.time = 'river_bgc_time'

    nsm = out.createVariable('river_nsm', 'f8', ('river_bgc_time', 'river'))
    nsm.long_name = 'river runoff small phytoplankton N'
    nsm.units = 'mol/kg'
    nsm.time = 'river_bgc_time'

    nh4 = out.createVariable('river_nh4', 'f8', ('river_bgc_time', 'river'))
    nh4.long_name = 'river runoff ammonia'
    nh4.units = 'mol/kg'
    nh4.time = 'river_bgc_time'

    no3 = out.createVariable('river_no3', 'f8', ('river_bgc_time', 'river'))
    no3.long_name = 'river runoff nitrate'
    no3.units = 'mol/kg'
    no3.time = 'river_bgc_time'

    o2 = out.createVariable('river_o2', 'f8', ('river_bgc_time', 'river'))
    o2.long_name = 'river runoff oxygen'
    o2.units = 'mol/kg'
    o2.time = 'river_bgc_time'

    pdet = out.createVariable('river_pdet', 'f8', ('river_bgc_time', 'river'))
    pdet.long_name = 'river runoff detrital phosphorus'
    pdet.units = 'mol/kg'
    pdet.time = 'river_bgc_time'

    rpo4 = out.createVariable('river_po4', 'f8', ('river_bgc_time', 'river'))
    rpo4.long_name = 'river runoff phosphate'
    rpo4.units = 'mol/kg'
    rpo4.time = 'river_bgc_time'

    srdon = out.createVariable('river_srdon', 'f8', ('river_bgc_time', 'river'))
    srdon.long_name = 'river runoff semi-refractory DON'
    srdon.units = 'mol/kg'
    srdon.time = 'river_bgc_time'

    srdop = out.createVariable('river_srdop', 'f8', ('river_bgc_time', 'river'))
    srdop.long_name = 'river runoff semi-refractory DOP'
    srdop.units = 'mol/kg'
    srdop.time = 'river_bgc_time'

    sldon = out.createVariable('river_sldon', 'f8', ('river_bgc_time', 'river'))
    sldon.long_name = 'river runoff semilabile DON'
    sldon.units = 'mol/kg'
    sldon.time = 'river_bgc_time'

    sldop = out.createVariable('river_sldop', 'f8', ('river_bgc_time', 'river'))
    sldop.long_name = 'river runoff semilabile DOP'
    sldop.units = 'mol/kg'
    sldop.time = 'river_bgc_time'

    sidet = out.createVariable('river_sidet', 'f8', ('river_bgc_time', 'river'))
    sidet.long_name = 'river runoff detrital silicon'
    sidet.units = 'mol/kg'
    sidet.time = 'river_bgc_time'

    silg = out.createVariable('river_silg', 'f8', ('river_bgc_time', 'river'))
    silg.long_name = 'river runoff large phytoplankton silicon'
    silg.units = 'mol/kg'
    silg.time = 'river_bgc_time'

    simd = out.createVariable('river_simd', 'f8', ('river_bgc_time', 'river'))
    simd.long_name = 'river runoff medium phytoplankton silicon'
    simd.units = 'mol/kg'
    simd.time = 'river_bgc_time'

    sio4 = out.createVariable('river_sio4', 'f8', ('river_bgc_time', 'river'))
    sio4.long_name = 'river runoff silicate'
    sio4.units = 'mol/kg'
    sio4.time = 'river_bgc_time'

    nsmz = out.createVariable('river_nsmz', 'f8', ('river_bgc_time', 'river'))
    nsmz.long_name = 'river runoff small zooplankton N'
    nsmz.units = 'mol/kg'
    nsmz.time = 'river_bgc_time'

    nmdz = out.createVariable('river_nmdz', 'f8', ('river_bgc_time', 'river'))
    nmdz.long_name = 'river runoff medium zooplankton N'
    nmdz.units = 'mol/kg'
    nmdz.time = 'river_bgc_time'

    nlgz = out.createVariable('river_nlgz', 'f8', ('river_bgc_time', 'river'))
    nlgz.long_name = 'river runoff large zooplankton N'
    nlgz.units = 'mol/kg'
    nlgz.time = 'river_bgc_time'

nrivers = len(rivers)
zoo = np.zeros((len(alk), nrivers))
foo = np.ones((len(alk), nrivers))
twee = len(alk)

print('alk', alk)
print('dic', dic)
print('no3', din)
print('po4', po4)
print('dsi', dsi)

alk_a = np.broadcast_to(np.array(alk).reshape((twee,1)), foo.shape)
dic_a = np.broadcast_to(np.array(dic).reshape((twee,1)), foo.shape)
din_a = np.broadcast_to(np.array(din).reshape((twee,1)), foo.shape)
po4_a = np.broadcast_to(np.array(po4).reshape((twee,1)), foo.shape)
dsi_a = np.broadcast_to(np.array(dsi).reshape((twee,1)), foo.shape)
do_a = np.broadcast_to(np.array(do).reshape((twee,1)), foo.shape)

out.variables['river_bgc_time'][:] = ttime
out.variables['river_alk'][:] = alk_a
out.variables['river_cadet_arag'][:] = zoo
out.variables['river_cadet_calc'][:] = zoo
out.variables['river_dic'][:] = dic_a
out.variables['river_fed'][:] = 3.0e-8 * foo
out.variables['river_fedet'][:] = zoo
out.variables['river_fedi'][:] = zoo
out.variables['river_felg'][:] = zoo
out.variables['river_femd'][:] = zoo
out.variables['river_fesm'][:] = zoo
out.variables['river_ldon'][:] = foo * scale * 10.0
out.variables['river_ldop'][:] = foo * scale * 0.25
out.variables['river_lith'][:] = zoo
out.variables['river_lithdet'][:] = zoo
out.variables['river_nbact'][:] = zoo
out.variables['river_ndet'][:] = zoo
out.variables['river_ndi'][:] = zoo
out.variables['river_nlg'][:] = zoo
out.variables['river_nmd'][:] = zoo
out.variables['river_nsm'][:] = zoo
out.variables['river_nh4'][:] = 0.02 * din_a
out.variables['river_no3'][:] = din_a
out.variables['river_o2'][:] = do_a
out.variables['river_pdet'][:] = zoo
out.variables['river_po4'][:] = po4_a
out.variables['river_srdon'][:] = foo * scale * 5.0
out.variables['river_srdop'][:] = foo * scale * 0.125
out.variables['river_sldon'][:] = foo * scale * 5.0
out.variables['river_sldop'][:] = foo * scale * 0.125
out.variables['river_sidet'][:] = zoo
out.variables['river_silg'][:] = zoo
out.variables['river_simd'][:] = zoo
out.variables['river_sio4'][:] = dsi_a
out.variables['river_nsmz'][:] = zoo
out.variables['river_nmdz'][:] = zoo
out.variables['river_nlgz'][:] = zoo

out.close()
