#!/usr/bin/env python

"""
Prototype code for a visual DESI Tile Picker
"""

import sys, os
import numpy as np

import math
import datetime
import time

def xyzrot(axis, v, theta):
    '''Rotate `v [x,y,z]` counter-clockwise around `axis` by `theta` in radians'''
    
    #- From https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    #- Uses Euler-Rodrigues formula for constructing 3D rotation matrix
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / np.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    return np.dot(rotation_matrix(axis, theta), v)

def radec2xyz(ra, dec):
    '''ra,dec in degrees -> x,y,z on unit sphere'''
    theta = np.radians(90-dec)
    phi = np.radians(ra)
    x = np.sin(theta)*np.cos(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(theta)
    return x, y, z

def xyz2radec(x,y,z):
    '''x,y,z on unit sphere -> ra,dec in degrees'''
    phi = np.arctan2(y, x)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(r, z)
    ra = (np.degrees(phi) + 360) % 360
    dec = 90 - np.degrees(theta)
    return ra, dec

def visibility(lst, lat, airmass, n=500):
    '''lst, lat in degrees; returns (ra,dec) arrays'''
    dec = 90.0 - np.ones(n) * np.degrees(np.arccos(1/airmass))
    ra = np.linspace(0, 360, n)
    
    x, y, z = radec2xyz(ra, dec)
    axis = [-np.sin(np.radians(lst)), np.cos(np.radians(lst)), 0]
    xx, yy, zz = xyzrot(axis, [x,y,z], np.radians(90-lat))
    ra, dec = xyz2radec(xx, yy, zz)
    return ra, dec

def kpnotime2lst(year, month, day, hour, minute, second=0):
    '''Returns LST at KPNO in degrees given KPNO local date and time'''
    #- Local time -> Universal time
    t = datetime.datetime(year, month, day, hour, minute, second) + datetime.timedelta(hours=7)

    #- time relative to noon UT on 2000-01-01
    dt = t - datetime.datetime(2000, 1, 1, 12, 0, 0)
    days = dt.days + dt.seconds / (24*3600)
    
    #- Approximation
    mayall_longitude_degrees = -(111 + 35/60. + 59.6/3600)
    LST_hours = ((18.697374558 + 24.06570982441908 * days) + mayall_longitude_degrees/15) % 24
    LST_degrees = LST_hours * 15
    
    return LST_degrees


def plot_visibility(tiles, airmass=1.5, width=800, outfile=None, title=None):
    '''
    Plot tile visibility using bokeh, optionally outputting to html file
                    
    Args:
        tiles: np.array of tiles with columns TILEID, RA, DEC, PROGRAM

    Options:
        airmass (float): airmass visibility window to highlight
        width (float): plot width in pixels
        outfile (str): output filename
        title (str): plot title

    If using within a Jupyter notebook, keep outfile=None and it will
    display within the notebook instead.
    '''

    import bokeh.plotting as bk
    from bokeh.models.widgets import tables as bktables
    from bokeh.models import CustomJS, Slider, DateSlider, HoverTool

    if title is None:
        title = 'DESI Tile Picker'

    lat = 31.964  #- KPNO Latitude
    now = time.localtime()
    # lst = kpnotime2lst(2019, 1, 1, 16, 0)
    lst = kpnotime2lst(now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min)

    ra0_1, dec_1 = visibility(0, lat, airmass, n=350)
    ra0_2, dec_2 = visibility(0, lat, 2.0, n=350)
    
    ra0 = np.concatenate([ra0_1, ra0_2])
    dec = np.concatenate([dec_1, dec_2])

    # ra0, dec = visibility(0, lat, airmass)
    ra = (ra0 + lst) % 360

    tiledata = dict(
        ra=tiles['RA'],
        dec=tiles['DEC'],
        tileid=tiles['TILEID'],
        program=tiles['PROGRAM'],
        selected=np.ones(len(tiles), dtype=bool),
    )
    for colname in ['STAR_DENSITY', 'EBV_MED']:
        if colname in tiles.dtype.names:
            tiledata[colname] = tiles[colname]
    
    source = bk.ColumnDataSource(data=dict(ra0=ra0, ra=ra, dec=dec))
    tile_source = bk.ColumnDataSource(data=tiledata)
    
    colformat = bktables.NumberFormatter(format='0,0.00')
    columns = [
        bktables.TableColumn(field='tileid', title='TILEID', width=60),
        bktables.TableColumn(field='ra', title='RA', formatter=colformat),
        bktables.TableColumn(field='dec', title='DEC', formatter=colformat),
    ]
    
    for colname in ['STAR_DENSITY', 'EBV_MED']:
        if colname in tiledata:
            columns.append(bktables.TableColumn(field=colname, title=colname, formatter=colformat))

    
    columns.append(bktables.TableColumn(field='selected', title='Selected'))

    tiletable = bktables.DataTable(columns=columns, source=tile_source, width=width)

    tile_source.selected.js_on_change('indices', CustomJS(args=dict(s1=tile_source), code="""
        var inds = cb_obj.indices;
        var d1 = s1.data;

        for (var i=0; i<d1['selected'].length; i++) {
            d1['selected'][i] = false;
        }

        for (var i = 0; i < inds.length; i++) {
            d1['selected'][inds[i]] = true;
        }
        s1.change.emit();
    """)
    )
    
    if outfile is not None:
        if os.path.exists(outfile):
            os.remove(outfile)
        bk.output_file(outfile, title=title)
    else:
        bk.output_notebook()
    
    fig = bk.figure(width=width, height=width//2,
                    tools=['pan', 'box_zoom', 'box_select', 'save', 'undo', 'reset'], active_drag='box_select')

    #- TODO: adjust size depending upon number of circles drawn
    circlesize = 2.0

    tile_circles = fig.circle('ra', 'dec', source=tile_source, radius=circlesize, alpha=0.8)
    fig.circle('ra', 'dec', source=source, size=0.1, color='black')
    fig.xaxis.axis_label = 'RA [degrees]'
    fig.yaxis.axis_label = 'Declination [degrees]'
    fig.title.text = '{}; airmass < {:.1f}, 2 @ LST={:.1f} deg'.format(
        title, airmass, lst)

    #- KPNO Date and local time sliders
    start = datetime.datetime(2019, 1, 1, 0, 0, 0)
    today = datetime.datetime(now.tm_year, now.tm_mon, now.tm_mday, 0, 0, 0)
    date_slider = DateSlider(start=start, end = start+datetime.timedelta(days=365), value=today,
        format = "%B %d", title='Date of sunset', width=width)

    start = datetime.datetime(2019, 1, 1, 16, 0, 0)
    timenow = datetime.datetime(2019, 1, 1, now.tm_hour, now.tm_min, now.tm_sec)
    localtime_slider = DateSlider(start = start, end = start+datetime.timedelta(hours=16), value=timenow,
        format = "%H:%M", title='KPNO local time', width=width)

    callback = CustomJS(
        args=dict(source=source, date_slider=date_slider, localtime_slider=localtime_slider,
                  airmass=airmass, fig=fig),
        code="""
            // First set times as if they were UTC
            var t = new Date(localtime_slider.value);
            var d = new Date(date_slider.value);
            if (t.getUTCHours() < 12) {
                d.setTime(date_slider.value + 24*3600*1000);
            } else {
                d.setTime(date_slider.value);
            }

            d.setUTCHours(t.getUTCHours());
            d.setUTCMinutes(t.getUTCMinutes());
            d.setUTCSeconds(0);

            // Correct to KPNO local time
            // d object still thinks in UTC, which is 7 hours ahead of KPNO
            d.setTime(d.getTime() + 7*3600*1000);

            // noon UT on 2000-01-01
            var reftime = new Date();
            reftime.setUTCFullYear(2000);
            reftime.setUTCMonth(0);   // Months are 0-11 (!)
            reftime.setUTCDate(1);    // Days are 1-31 (!)
            reftime.setUTCHours(12);
            reftime.setUTCMinutes(0);
            reftime.setUTCSeconds(0);

            // time difference in days (starting from milliseconds)
            var dt = (d.getTime() - reftime.getTime()) / (24*3600*1000);

            // Convert to LST
            var mayall_longitude_degrees = -(111 + 35/60. + 59.6/3600);
            var LST_hours = ((18.697374558 + 24.06570982441908 * dt) + mayall_longitude_degrees/15) % 24;
            var LST_degrees = LST_hours * 15;
            
            var data = source.data;
            var ra = data['ra'];
            var ra0 = data['ra0'];
            for (var i = 0; i < ra.length; i++) {
                ra[i] = (ra0[i] + LST_degrees) % 360;
            }
            // DESI tile picker; airmass < {:.1f} @ LST={:.1f} deg
            fig.title.text = "DESI tile picker; airmass < " + airmass + " @ LST="+LST_degrees.toFixed(2);
            source.change.emit();

            // console.log(d, reftime, LST_degrees);
        """)

    date_slider.js_on_change('value', callback)
    localtime_slider.js_on_change('value', callback)

    fig.x_range.start = 360
    fig.x_range.end = 0
    
    hovertips = [
        ("Tile ID", "@tileid"),
        ("RA,dec", "@ra, @dec"),
        ("Program", "@program"),
    ]
    fig.add_tools(HoverTool(tooltips=hovertips, callback=None, renderers=[tile_circles,]))
    text_opts = dict(text_align='center', text_baseline='middle', text_alpha=0.3, text_font_size=dict(value='10pt'))
    fig.text(180, 84, ['North',], **text_opts)
    fig.text(180, -20, ['South',], **text_opts)
    fig.text(352, 50, ['East',], angle=np.pi/2, **text_opts)
    fig.text(7, 50, ['West',], angle=np.pi/2, **text_opts)

    if outfile is None:
        bk.show(bk.Column(fig, date_slider, localtime_slider, tiletable))
    else:
        bk.save(bk.Column(fig, date_slider, localtime_slider, tiletable))

    # bk.reset_output()

def main():
    import argparse
    
    parser = argparse.ArgumentParser(usage = "{prog} [options]")
    parser.add_argument("-i", "--input", type=str, required=True, help="input tiles file (FITS)")
    parser.add_argument("-o", "--output", type=str, required=True, help="output html file")
    parser.add_argument("-t", "--title", type=str,  help="title to include on plot")
    # parser.add_argument("-v", "--verbose", action="store_true", help="some flag")
    
    args = parser.parse_args()
    
    from astropy.io import fits
    citiles = fits.getdata(args.input)
    plot_visibility(citiles, airmass=1.5, outfile=args.output, title=args.title)

if __name__ == '__main__':
    main()











