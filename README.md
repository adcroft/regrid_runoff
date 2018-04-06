# regrid_runoff

Conservatively re-grids runoff data on a regular latitude-longitude grid to the nearest coastal wet-point
of a curvilinear ocean grid.

The MOM6 ocean model does not require river input to be smoothed out - something that is sometimes done
in other models for stability.
Runoff data is often provided on a uniform/regular latitude-longitude grid but has some runoff in two lines of
cells on some coasts, presumably because an actual river discharge point spanned cell boundaries.
The regridding algorithm here regrids runoff data to the ocean model grid and at the same time assures
that the discharge is "sharpened" in the direction normal to coastlines so that the discharge is only
ever one cell away from and land-masked cell.

## Usage

```
./regrid_runoff.py -h
./regrid_runoff.py ocean_hgrid.nc ocean_mask.nc -m mask river_runoff.nc new_file.nc
```

## Algorithm

An [example notebook](https://github.com/adcroft/regrid_runoff/blob/master/Regrid%20runoff%20data.ipynb)
(used for development) illustrates the procedure:
1. Define and find the "coastal cells" on the target ocean model grid
   - I currently use the definition "Any wet cell with a land cell to the north, south, east, or west".
   - Optionally this can include the wet cells diagonally-adjacent to land but this tends to broaden
   the discharge regions at coastal head lands.
2. Use an "[ice-9](https://en.wikipedia.org/wiki/Ice-nine)" algorithm to label each cell on the model
grid with the nearest "coastal cell".
   - The ice-9 algorithm implemented here is imperfect for expediency in that it as directional biases
   built-in.
3. Find the models cells among which to divide a river cell flux when the river cells are coarse
   relative to the ocean cells.
4. Find the model cell into which to send all flux in a river-grid cell when the river cells are fine
   relative to the ocean cells.
   There is no area-weight re-mapping.
   - We can take advantage of the regular nature of the river-grid, rather than generating a general
   intersection of grids
   1. Using the latitude-longitude of each model cell-center to enter a model-cell id in the river cells.
   This will be incomplete but cover a good fraction of the river cells.
   2. Use an "ice-9" algorithm to fill in the model-cell id for river cells not touched by in step 3.i.
   We consider this a first guess.
   3. For each river-grid cell, use a down-gradient search algorithm to correct the model-cell id.
5. Construct a sparse matrix, **A**, of river-areas from the indirect-indexing information calculated
   in step 3.
   - The matrix will have a column for each river-grid cell, and a row for each ocean model-grid cell.
   - The matrix might have multiple entries on a row but only one entry in a column.
6. The flux, o, on the model grid is "o = **A**.r / areacello" where "areacello" is the model cell area,
   and r is the flux on the river grid.

## Usage for ROMS

The --regional switch attempts to turn off rivers outside of your regional domain, but it instead
moves them all to the NE corner (-1,-1) which routes it to a coastal point near there. I then
set those values to zero in the make\_rivers\_clim script, not a completely satisfactory fix.
```
./regrid_runoff.py -h
./regrid_runoff.py ARCTIC4 --regional -f runoff_all.1980.nc Arctic4_runoff_1980.nc

```

Then some weird old scripts I have to get it into a rivers file:

ROMS SIDE PREP

1. First, the maskedge program needs to be run on the ROMS grid file.
It produces a list of the i,j pairs along the edge of the land mask.
Save it in a file and edit out the parts of the land mask you don't
want included (for NWGOA, i<209). Call this file maskedge.out.

    python maskedge.py roms_grid.nc > maskedge.out
    STDERR:  New body! 1, 449
    STDERR:  There are 1 bodies of water.

Note: if you get 0 bodies of water, check the \_FillValue attribute
on your masks. Some old grids have them set to 1.0, which is no
longer desired. Fix with:

```
ncatted -O -a _FillValue,mask_rho,d,, -a _FillValue,mask_u,d,, -a _FillValue,mask_v,d,, -a _FillValue,mask_psi,d,, CGOA_grid_4.nc
```

2. Run the add\_rivers program to create the (empty) ROMS file of i,j
river locations:

```
python add_rivers.py JRA_Arctic_1980.nc
```

3. Getting the remapped runoff into the ROMS rivers file:

```
python make_river_clim.py Arctic4_runoff_1980.nc JRA_Arctic_1980.nc
```

The runoff is in kg/m2/sec and ROMS wants m^3/sec. Also, the rivers outside
the domain all end up in the far NE coastal cell that maps from (-1,-1).
Someone should probably fix this in regrid_runoff.py.

4. Just because we can:

```
python squeeze_rivers.py JRA_Arctic_1980.nc squeeze.nc
mv squeeze.nc JRA_Arctic_1980.nc
```

5. Then there's the tracers... Seth gave me one year of daily values to
be reused year after year. We need a second, cyclic river\_time.
Plus we need to hack ROMS to accept only one value per day - the new
ONE\_TRACER\_SOURCE flag (or one value per horizontal location).

```
python add_temp.py JRA_Arctic_1980.nc
```

6. Finally, setting Vshape:

```
python set_vshape.py JRA_Arctic_1980.nc
```

7. If using dyes:

```
python add_dye.py
```

Or just use the run\_regrid script...
