# regrid_JRA_runoff

Conservatively re-grids JRA55-do runoff data to the nearest coastal wet-point of a MOM6 ocean grid.

The MOM6 ocean model does not require river input to be smoothed out - something that is sometimes done
in other models for stability.
The JRA55-do runoff data is on a 1/4-degree regular lat-lon grid but has some runoff in two lines of
cells on some coasts, presumably because an actual river discharge point spanned cell boundaries.
The regridding algorithm here regrids runoff data to the ocean model grid and at the same time assures
that the discharge is "sharpened" in the direction normal to coastlines so that the discharge is only
ever one cell away from and land-masked cell.

## Usage

```
./regrid_JRA_runoff.py -h
./regrid_JRA_runoff.py -sg ocean_hgrid.nc -om ocean_mask.nc -ra runoff_cell_area.15Dec2016.nc runoff_all.\*.15Dec2016.padded.nc 
```

## Algorithm

An [example notebook](https://gist.github.com/adcroft/8894494254777d7da0863b6a3cf69722) (used for
development) illustrates the procedure:
1. Define and find the "coastal cells" on the target ocean model grid
   - I currently use the definition "Any wet cell with a land cell to the north, south, east, or west".
   - Optionally this can include the wet cells diagonally-adjacent to land but this tends to broaden
   the discharge regions at coastal head lands.
2. Use an "[ice-9](https://en.wikipedia.org/wiki/Ice-nine)" algorithm to label each cell on the model
grid with the nearest "coastal cell".
   - The ice-9 algorithm implemented here is imperfect for expediency in that it as directional biases
   built-in.
3. Find the model cell into which to send any flux in a river-grid cell.
   - Note: this currently is all-or-nothing.
   There is no partial re-mapping, or splitting of a river cell into multiple ocean cells.
   This is a reasonable approximation so long as the model cells are generally broader than actual
   river mouths, and broader/comparable to the river-grid cells.
   - We can take advantag of the regular nature of the river-grid, rather than generating a general
   intersection of grids
   1. Using the latitude-longitude of each model cell-center to enter a model-cell id in the river cells.
   This will be incomplete but cover a good fraction of the river cells.
   2. Use an "ice-9" algorithm to fill in the model-cell id for river cells not touched by in step 3.i.
   We consider this a first guess.
   3. For each river-grid cell, use a down-gradient search algorithm to correct the model-cell id.
4. Construct a sparse matrix, **B**, of river-areas from the indirect-indexing information calculated
   in step 3.
   - The matrix will have a column for each river-grid cell, and a row for each ocean model-grid cell.
   - The matrix might have multiple entries on a row but only one entry in a column.
5. The flux, o, on the model grid is "o = **B**.r / areacello" where "areacello" is the model cell area,
   and r is the flux on the river grid.
