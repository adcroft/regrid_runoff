TESTDATA = ocean_hgrid.nc ocean_static.nc runoff.daitren.clim.10FEB2011.nc
REGRID_RUNOFF = ../regrid_runoff.py
TARGETS = test.nc test-f.nc test-fz.nc

all: $(TARGETS)
testdata: $(TESTDATA)
ocean_hgrid.nc ocean_static.nc:
	wget ftp://ftp.gfdl.noaa.gov/perm/Alistair.Adcroft/MOM6-testing/OM4_05/$@
runoff.daitren.clim.10FEB2011.nc:
	wget http://data1.gfdl.noaa.gov/~nnz/mom4/COREv2/data_IAF/CORRECTED/calendar_years/$@
help:
	$(REGRID_RUNOFF) --help
test.nc: $(TESTDATA)
	$(REGRID_RUNOFF) ocean_hgrid.nc ocean_static.nc -m wet runoff.daitren.clim.10FEB2011.nc $@ -p
test-f.nc: $(TESTDATA)
	$(REGRID_RUNOFF) ocean_hgrid.nc ocean_static.nc -m wet runoff.daitren.clim.10FEB2011.nc $@ -p -f
test-fz.nc: $(TESTDATA)
	$(REGRID_RUNOFF) ocean_hgrid.nc ocean_static.nc -m wet runoff.daitren.clim.10FEB2011.nc $@ -p -f -z
check:
	md5sum test*.nc
	md5sum -c check.md5
