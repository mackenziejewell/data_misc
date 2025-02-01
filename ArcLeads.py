# Opening locally-stored Willmes et al Al (2023) sea ice lead data
# Willmes, Sascha; Heinemann, Günther; Reiser, Fabian (2023): ArcLeads: Daily sea-ice lead maps for the Arctic, 2002-2021, NOV-APR [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.955561

# DEPENDENCIES:
import xarray as xr
import numpy as np
import numpy.ma as ma
# import cartopy
# import cartopy.crs as ccrs
# import cartopy.feature as cfeat
# import matplotlib.ticker as mticker
from datetime import datetime, timedelta
# from metpy.units import units
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# import matplotlib.colors
# from matplotlib.offsetbox import AnchoredText
# from shapely import wkt



# FUNCTIONS:
#---------------------------------------------------------------------
def open_local_file(year, 
                        crop_i = [1400, 2200], 
                        crop_j = [600, 1200], 
                        main_dir = '/Volumes/Seagate_Jewell/KenzieStuff/ArcLeads/'):
    
    """Read in local lead file from Willmes et al. (2023) dataset.
        Dates are adjusted since they are off by one year in second portion of timeseries. 
        Default is to crop spatially, which signifcantly saves time

INPUT: 
- year: (int) year2 of file to open (files are stored Nov. year1 -- Apr. year2)
- crop_i: (2 x 1) array of indices to crop i dimension (axis 0)
- crop_j: (2 x 1) array of indices to crop j dimension (axis 1)
- main_dir: (str) directory where data files are stored

OUTPUT:
- ds2: (xarray Dataset) dataset with cropped data and correct dates

Latest recorded update:
01-31-2025
    """

    # crop_i is range to crop i index
    # crop_j is range to crop j index
    
    # JFMA year corresponds to second year of file
    # Nov. year1 -- Apr. year2
    
    # open data file
    #0=clouds, 1=land, 2=sea ice, 3=artefacts, 4=leads, 5=open water.
    ds = xr.open_dataset(main_dir+f'{year-1}{str(year)[-2:]}_ArcLeads.nc')
    ds.close()
    
    # interpret file dates
    # note that seocond set of dates (those after Jan 1) list the wrong (previous) year
    dates1 = [datetime.strptime(time, '%Y%m%d') for time in ds.time.values[:61].astype(str)]
    dates2 = [datetime(int(time[:4])+1, int(time[4:6]), int(time[6:])) for time in ds.time.values[61:].astype(str)]
    correct_dates = np.append(dates1, dates2)

    # crop
    ai, bi = crop_i[0], crop_i[1]
    aj, bj = crop_j[0], crop_j[1]
    
    # grab coordinates
    lats = ds.lat[ai:bi, aj:bj]
    lons = ds.lon[ai:bi, aj:bj]
    
    # lead map
    lead_map = ds.leadmap[:, ai:bi, aj:bj]

    # create new dataset (cropped and with correct dates)
    variables = {'lon': (('nrow', 'ncol'), lons.values, {'units': 'degreeE'}),
                'lat': (('nrow', 'ncol'), lats.values, {'units': 'degreeN'}),
                'leadmap': (('time', 'nrow', 'ncol'), lead_map.values, 
                            {'units': ds.leadmap.units, 'CRS': ds.leadmap.CRS}),
                }
    
    coordinates = {'nrow': np.arange(ds.dims['nrows'])[ai:bi],
                    'ncol': np.arange(ds.dims['ncols'])[aj:bj],
                    'time': correct_dates 
                    }

    # Create the Dataset
    ds2 = xr.Dataset(data_vars = variables,
                     coords = coordinates,
                     attrs=ds.attrs)

    return ds2
