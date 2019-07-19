DESI Tile Picker
================

A prototype interactive tile picker for [DESI](https://desi.lbl.gov)
commissioning.

## Usage ##

From the command line:

    python tilepicker.py -i tiles.fits -o tiles.html

Or from within a Jupyter notebook:

    import tilepicker
    from astropy.io import fits
    tiles = fits.getdata('tiles.fits')
    tilepicker.plot_visibility(tiles, airmass=1.5)

## Input / Output ##

Input tiles table must have columns TILEID, RA, DEC, PROGRAM.
This repository includes an example `tiles.fits` file, used by
the DESI Commissioning Instrument observing program (v7).

Output is an interactive tile selector like [tiles-example.html](tiles-example.html).

## To do ##

* Moon location
* Ecliptic
* Configurable selection of which columns are included in table

<hr/>
**Stephen Bailey**<br/>
Lawrence Berkeley National Lab<br/>
Summer 2019
