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
