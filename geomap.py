# Functions for plotting geographic data

# DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.ticker as mticker
from datetime import datetime, timedelta
from metpy.units import units
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.colors
from matplotlib.offsetbox import AnchoredText
from shapely import wkt



# FUNCTIONS:
#---------------------------------------------------------------------
def fix_cartopy_vectors(u, v, uvlats):
    
    """Function to output vector components for plotting in cartopy. 
    
    Reads in vectors and associated latitudes, return fixed vectors. 
    
    Cartopy doesn't know meters per degree increase toward the pole 
    in zonal direction as cosine(lat) when reprojecting for vectors 
    given as m/s (where  a meter covers different amount of a degree 
    depending on latitude), we need to rescale the u (east-west) speeds. 
    otherwise at high latitudes, u will be drawn much larger than 
    they should which will give incorrect angles

INPUT: 
- u: (N x M) array of eastward velocity component (m/s)
- v: (N x M) array of northward velocity component (m/s)
- uvlats: (N x M) array latitudes associated with u,v vectors

OUTPUT:
- u_fixed: (N x M) array of u with correct angle and magnitude for plotting
- v_fixed: (N x M) array of v with correct angle and magnitude for plotting

Latest recorded update:
12-17-2024
    """
    
    # FIX ANGLE
    # fix u scale to be correct relative to v scale
    #----------------------------------------------
    # for u (m/s), m/degree decreases as cos(lat)
    u_fixed = u/np.cos(uvlats/180*np.pi) 
    v_fixed = v  # v does not change with latitude

    # FIX MAGNITUDE
    # scale fixed u,v to have correct magnitude 
    #-----------------------
    # original magnitude 
    orig_mag = ma.sqrt(u**2+v**2)
    # new magnitude
    fixed_mag = ma.sqrt(u_fixed**2+v_fixed**2)
    u_fixed = u_fixed*(orig_mag/fixed_mag)
    v_fixed = v_fixed*(orig_mag/fixed_mag)

    return u_fixed, v_fixed


def add_land(ax, scale = '50m', color='gray', alpha=1, fill_dateline_gap = True, zorder=2):
    
    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature land feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = land color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- alpha = land opacity (default: 1)
- zorder: drawing order of land layer (default: 2)
- fill_dateline_gap: specify whether to fill gap in cartopy land feature along 
   dateline that crosses Russia and Wrangel Island (default: True)

Latest recorded update:
12-17-2024

    """
    
        
    # grab land from cfeat.NaturalEarthFeature
    #-----------------------------------------
    ax.add_feature(cfeat.NaturalEarthFeature(category='physical', name='land', 
                                             scale=scale, facecolor=color),
                                             alpha = alpha, zorder = zorder)

    # if specified, fill dateline gap in land feature with shapely polygons
    if fill_dateline_gap == True:
        # generate polygon to fill line across Wrangel Island and line across Russia
        WKT_fill_Wrangel = 'POLYGON ((-180.1 71.51,-180.1 71.01,-179.9 71.01,-179.9 71.51,-180.1 71.51))'
        poly1 = wkt.loads(WKT_fill_Wrangel)
        ax.add_geometries([poly1], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)
        WKT_fill_Russia = 'POLYGON ((-180.1 65.1,-180.1 68.96,-179.9 68.96,-179.9 65.1,-180.1 65.1))'
        poly2 = wkt.loads(WKT_fill_Russia)
        ax.add_geometries([poly2], crs=ccrs.PlateCarree(), 
              facecolor=color, edgecolor=color, alpha = alpha, zorder=zorder)



def add_coast(ax, scale = '50m', color='gray', linewidth = 1, alpha=1, zorder=3):

    """Add land feature to cartopy figure
    
INPUT:
- ax: cartopy figure axis
- scale = NaturalEarthFeature coast feature scale (e.g. '10m', '50m', '110m')
        (default: '50m')
- color = coastline color (e.g. 'k' or [0.9,0.6,0.5]) (default: 'gray')
- linewidth = coastline linewidth (default: 1)
- alpha = coastline opacity (default: 1)
- zorder: drawing order of coast layer (default: 3)

Latest recorded update:
12-17-2024

    """

    # coastline
    #----------
    ax.coastlines(scale, color=color, linewidth=linewidth, alpha = alpha, zorder = zorder)



def add_grid(ax, lats = None, lons = None, linewidth = 1, color = 'gray', alpha=0.5, zorder = 4): 
    
    """Add specified gridlines to cartopy figure.
    
INPUT:
- ax: cartopy figure axis
- lats: None or array of latitudes to plot lines (default: None)
- lons: None or array of latitudes to plot lines (default: None)
- linewdith: grid line linewidths (default: 1)
- color: grid line color (default: 'gray')
- alpha: line transparency (default: 0.5)
- zorder: drawing order of gridlines layer (default: 4)


Latest recorded update:
12-17-2024

    """
        
    # give gridline specifications
    #-----------------------------
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=linewidth, color=color, alpha=alpha, zorder = zorder)

    # add the longitude gridlines
    #----------------------------
    if lons is None:
        gl.xlocator = mticker.FixedLocator([])
    else:
        # shift all longitudes from [0,360] to [180,-180]
        lons = np.concatenate((lons[(lons>180)]-360,lons[(lons<=180)]))
        gl.xlocator = mticker.FixedLocator(lons)

        
    # add the latitude gridlines
    #----------------------------
    if lats is None:
        gl.ylocator = mticker.FixedLocator([])
    else:
        gl.ylocator = mticker.FixedLocator(lats)


def add_date(fig, ax, dt_obj, date_format = '%b %d, %Y (%H:%M UTC)', method = 'anchor',
             boxstyle="round,pad=0.,rounding_size=0.2", facecolor = 'black', edgecolor = 'black',
             zorder = 10,
             
             anchor_loc = 4, anchor_prop = {'size': 20, 'color':'white'},
             

             x = 0.02, y= 0.05, textcolor = 'white',fontsize=15): 
    
    """Add date label to cartopy plot.
    
INPUT:
- fig: cartopy figure
- ax: cartopy figure axis
- dt_obj: datetime object of date for plotted data 
            OR
          string with text to show (date format already provided (e.g. 'Dec 20, 2018 (6:00 UTC)')
          
IF dt_obj IS DATETIME OBJECT:
- date_format: str, format to display date (default: '%b %d, %Y (%H:%M UTC)')
    - example 1: '%b %d, %Y (%H:%M UTC)' could give 'Dec 20, 2018 (6:00 UTC)'
    - example 2: '%m-%d-%Y' could give '12-20-2018'
    
- method: method to place the date label (either 'anchor' for AnchoredText or 'manual' to place manually).
        (default: 'anchor')
- boxstyle: anchor box shape style (default: "round,pad=0.,rounding_size=0.2")
- facecolor: color of bounding box (default: 'black')
- edgecolor: color of bounding box edge (default: 'black')
- zorder: drawing order of date layer (default: 10)

IF METHOD = 'anchor':
- anchor_loc: anchor text location (default: 4)
- anchor_prop: anchor properties dictionary (default: {'size': 20, 'color':'white'})

IF METHOD = 'manual':
- x: x-location of figure extent to place date
- y: y-location of figure extent to place date
- textcolor: color oftext (default: 'white')
- fontsize: fontsize of text (defult: 15)

OUTPUT:
- input plot with added date label

Latest recorded update:
12-17-2024
    """

    
    assert method in ['anchor', 'manual'], f">>> method should be 'manual' or 'anchor', given: '{method}'"
    
    assert str(type(dt_obj)) in ["<class 'datetime.datetime'>", "<class 'str'>"], f">>> dt_obj should be datetime object or string, given: {str(type(dt_obj))}"
    
    
    # if given as datetime object, convert to specified date format
    if str(type(dt_obj)) == "<class 'datetime.datetime'>":
        date_text = dt_obj.strftime(date_format)
    
    # else, set date directly to given string object
    else:
        date_text = dt_obj
    
    
    # add text
    #---------
    if str(method) == 'anchor':
        at = AnchoredText(date_text, loc=anchor_loc, prop=anchor_prop)
        at.patch.set_boxstyle(boxstyle)
        at.patch.set_facecolor(facecolor)
        at.patch.set_edgecolor(edgecolor)
        at.zorder = zorder
        ax.add_artist(at)
    
    elif str(method) == 'manual':
        ax.text(x, y, date_text, 
                bbox=dict(boxstyle = boxstyle, facecolor=facecolor, edgecolor = edgecolor), 
                transform=ax.transAxes, fontsize=fontsize, 
                c=textcolor, verticalalignment='top', zorder = zorder);





def scalebar(ax, corner = (0.1, 0.1), stepsize = 50,  numsteps = 4, unit = 'km', 
             unit_label = None,
             colors=['k','w'], textsize=9, lw = 3, 
             unitpad = 0.025, labelpad = 0.025):

    """Add scalebar to to cartopy plot.
    
INPUT:
- ax: cartopy figure axis
- corner = (x,y) of scalebar's lower left corner, in axes coordinates [0,1] (default: (0.1,0.1))
- stepsize: distance between scalebar ticks (default: 50)
- numsteps: number of scalebar ticks (default: 4)
- unit: unit of scalebar (default: 'km')
- unit_label: unit label to display (if None, will use unit) (default: None)
- colors: list of 2 alternating colors for scalebar ticks (default: ['k','w'])
- lw: linewidth of scalebar (default: 3)
- textsize: size of scalebar text (default: 9)
- unitpad: padding between unit label and scalebar (default: 0.025)
- labelpad: padding between scalebar ticks and labels (default: 0.025)

Latest recorded update:
01-30-2025
    """

    # convert step size to m
    step = (stepsize * units(unit)).to('m').magnitude

    # convert lower left corner of scalebar
    # from figure display coordinates to projected coordinates
    # follow: https://stackoverflow.com/questions/56662941/cartopy-convert-point-from-axes-coordinates-to-lat-lon-coordinates
    xi, yi = corner
    # convert from Axes coordinates to display coordinates
    (xd, yd) = ax.transAxes.transform((xi, yi))
    # convert from display coordinates to data coordinates
    (x0, y0) = ax.transData.inverted().transform((xd, yd))

    # top and bottom label locations

    # number labels
    # convert label pad to axes coordinates
    # convert from Axes coordinates to display coordinates
    (xd, yd) = ax.transAxes.transform((xi, yi+labelpad))
    # convert from display coordinates to data coordinates
    (xxx, labely) = ax.transData.inverted().transform((xd, yd))

    # unit labels
    unitx = x0 + ((numsteps/2-1)*step)
    # convert label pad to axes coordinates
    # convert from Axes coordinates to display coordinates
    (xd, yd) = ax.transAxes.transform((xi, yi-unitpad))
    # convert from display coordinates to data coordinates
    (xxx, unity) = ax.transData.inverted().transform((xd, yd))

    # print(step_m)
    # figure_projection = ax.projection
    
    # plot scalebar
    x1 = x0 - step - step/1000
    if numsteps%2==0:
        x2 = x0 + ((numsteps-1)*step) + step/1000  
    else:
        x2 = x0 + ((numsteps-1)*step) 
    ax.plot([x1, x2], [y0, y0], lw=lw+1, c=colors[0], zorder=100)
    
    for ii in range(numsteps+1):
        
        x1 = x0 + (step*(ii-1))
        x2 = x0 + step*ii
        
        if ii < numsteps:
            ax.plot([x1, x2], [y0, y0], lw=lw, c=colors[ii%2], zorder=100)
        
        
        # label ticks
        dist = f'{step*ii/1000:.0f}'
        ax.text(x1, labely, dist, va='bottom', ha = 'center', size=textsize, zorder=100)
    
    
    if unit_label is not None:
        LABEL = unit_label
    else:
        LABEL = unit
    
    ax.text(unitx, unity, 
            LABEL, va='top', ha = 'center', size=textsize, zorder=100)
    



def northarrow(ax, loc = (0.1, 0.1), textsize=9,):

    """Add north arrow to to cartopy plot.
    
INPUT:
- ax: cartopy figure axis
- loc = (x,y) of arrow in axes coordinates [0,1] (default: (0.1,0.1))
- textsize: size of arrow label text (default: 9)

Latest recorded update:
01-30-2025
    """


    # convert desired arrow coordinates to projection
    # from figure display coordinates to projected coordinates
    # follow: https://stackoverflow.com/questions/56662941/cartopy-convert-point-from-axes-coordinates-to-lat-lon-coordinates
    xi, yi = loc
    # convert from Axes coordinates to display coordinates
    (xd, yd) = ax.transAxes.transform((xi, yi))
    # convert from display coordinates to data coordinates
    (x0, y0) = ax.transData.inverted().transform((xd, yd))

    # convert from data to cartesian coordinates
    proj_cart = ccrs.PlateCarree()
    (lon, lat) = proj_cart.transform_point(*(x0, y0), src_crs=ax.projection)

    # arrow_lon = -153.5
    # arrow_lat = 69.75
    
    ax.text(np.array([lon]), np.array([lat]), 'N', weight='bold', va = 'center', ha='center',
            size = textsize,    
             transform=ccrs.PlateCarree(), zorder=300)
    
    ax.quiver(np.array([lon]), np.array([lat+0.05]), 
                np.array([0]), np.array([1]), scale=28,
                width = 0.0005, 
                headaxislength = 650, headlength = 1000, headwidth=1000,
                transform=ccrs.PlateCarree(), zorder=300)
    