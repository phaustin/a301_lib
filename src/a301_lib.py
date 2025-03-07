"""
utility code for ATSC 301
"""

from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from pathlib import Path
import pyproj
from pyproj import Transformer, CRS
import xarray as xr

g=9.8  #m/s^2 don't worry about g(z) for this class
Rd=287.  #kg/m^3



def rowcol2latlon(tiffile,row,col):
    """
    return the latitude and longitude of a row and column in a geotif

    Parameters
    ----------

    tiffile: Path object
        path to a clipped tiffile
    row: float
       row of the pixel
    col: float
       column of the pixel

    Returns
    -------

    (lon, lat):  (float, float)
       longitude (deg east) and latitude (deg north) on
       a WGS84 datum
    
    """
    has_file = tiffile.exists()
    if not has_file:
        raise IOError(f"can't find {filename}, something went wrong above") 
    the_band = rioxarray.open_rasterio(tif_filename,masked=True)
    epsg_code = the_band.rio.crs.to_epsg()
    p_utm = CRS.from_epsg(epsg_code)
    p_latlon = CRS.from_epsg(4326)
    affine_transform = the_band.rio.transform()
    x, y = affine_transform*(row, col)
    crs_transform = Transformer.from_crs(p_utm, p_latlon)
    lat, lon = transform.transform(x,y) 
    return (lon, lat)


def make_pal(ax = None,vmin = None, vmax = None, palette = "viridis"):
    """
    return a dictionary containing a palette and a Normalize object

    Parameters
    ----------

    vmin: minimum colorbar value (optional, defaults to minimum data value)
    vmax: maximum colorbar value (optional, defaults to maximum data value)
    palette: string (optional, defaults to "viridis")

    Returns
    -------

    out_dict: dict
      dictionary with key:values  cmap:colormap, norm:Normalize
      

    """
    the_norm = Normalize(vmin=vmin, vmax=vmax, clip=False)
    pal = plt.get_cmap(palette)
    pal.set_bad("0.75")  # 75% grey for out-of-map cells
    pal.set_over("w")  # color cells > vmax white
    pal.set_under("k")  # color cells < vmin black
    out_dict=dict(cmap=pal,norm=the_norm)
    return out_dict

def calcScaleHeight(df):
    """
    Calculate the pressure scale height H_p
    
    Parameters
    ----------

    df: pd.DataFrame

    columns:
    
        T: vector (float)
          temperature (K)

        p: vector (float) of len(T)
          pressure (pa)

        z: vector (float) of len(T
          height (m)

    Returns
    -------
    
    Hbar: vector (float) of len(T)
      pressure scale height (m)
    
    """
    z=df['z'].values
    Temp=df['temp'].values
    dz=np.diff(z)
    TLayer=(Temp[1:] + Temp[0:-1])/2.
    oneOverH=g/(Rd*TLayer)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar

def calcDensHeight(df):
    """
    Calculate the density scale height H_rho
    
    Parameters
    ----------

    df: pd.DataFrame

    df columns:
    
        T: vector (float)
          temperature (K)

        p: vector (float) of len(T)
          pressure (pa)

        z: vector (float) of len(T
          height (m)
      
    Returns
    -------
    
    Hbar: vector (float) of len(T)
      density scale height (m)
    """
    z=df['z'].values
    Temp=df['temp'].values
    dz=np.diff(z)
    TLayer=(Temp[1:] + Temp[0:-1])/2.
    dTdz=np.diff(Temp)/np.diff(z)
    oneOverH=g/(Rd*TLayer) + (1/TLayer*dTdz)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar

c, h, k = 299792458.0, 6.62607004e-34, 1.38064852e-23
c1 = 2.0 * h * c ** 2.0
c2 = h * c / k

def calc_radiance(wavel, Temp):
    """
    Calculate the blackbody radiance
    
    Parameters
    ----------

      wavel: float or array
           wavelength (meters)

      Temp: float
           temperature (K)

    Returns
    -------

    Llambda:  float or arr
           monochromatic radiance (W/m^2/m/sr)
    """
    Llambda_val = c1 / (wavel**5. * (np.exp(c2 / (wavel * Temp)) - 1))
    return Llambda_val

def radiance_invert(wavel, L):
    """
    Calculate the brightness temperature
    
    Parameters
    ----------
      wavel: float
           wavelength (meters)
      L: float or array
           radiance (W/m^2/m/sr)
    
    Returns
    -------
    Tbright:  float or arr
           brightness temperature (K)
    """
    c, h, k = 299792458.0, 6.62607004e-34, 1.38064852e-23
    c1 = 2.0 * h * c ** 2.0
    c2 = h * c / k
    Tbright = c2 / (wavel * np.log(c1 / (wavel ** 5.0 * L) + 1.0))
    return Tbright


from pydantic import BaseModel

class Row_col_offsets(BaseModel):
    """
    row/col offsets needed to define
    rows and column ranges in a region
    """
    west: int
    south: int
    east: int
    north: int

class Clip_point(BaseModel):
    """
    a location within the clipped region 
    in geodetic coordinates
    """
    lon: float
    lat: float

class Bounding_box(BaseModel):
    """
    bounding box for rioxarray.clip_box in map coordinates
    """
    xmin: float
    ymin: float
    xmax: float
    ymax: float

class Cartopy_extent(BaseModel):
    """
    cartopy plotting extent in map coords
    """
    xmin: float
    xmax: float
    ymin: float
    ymax: float


def find_bounds(
    scene_ds: xr.DataArray,
    offsets: Row_col_offsets,
    image_point: Clip_point,
    scene_crs: pyproj.CRS | None = None ) -> Bounding_box: 
    """
    Given a point in x,y map coordinates and
    a set of row, column offsets, return a
    bounding box that can be used by rio.clip_box to clip
    the image

    Parameters
    ----------
    scene_ds: rioxarray DataArray
    offsets: west, south, east, north
    image_point: lon, lat for point in image
    scene_crs: Optional -- point in scene, if missing use scene_ds.rio.crs

    Returns
    -------

    bounds: west, south, east, north bounds for rio.clip_box
    """ 
    if scene_crs is None:
        #
        # turn the rasterio CRS into a pyproj.CRS
        #
        scene_crs = pyproj.CRS(scene_ds.rio.crs)
    #
    # step 1: find the image point in utm coords
    # need to transform the latlon point to x,y utm zone 10
    #
    full_affine = scene_ds.rio.transform()
    p_latlon = CRS.from_epsg(4326)
    transform = Transformer.from_crs(p_latlon, scene_crs)
    image_x, image_y = transform.transform(image_point.lat, image_point.lon)
    #
    # step 2: find the row and column of the image point
    #
    image_col, image_row = ~full_affine * (image_x, image_y)
    image_col, image_row = int(image_col), int(image_row)
    #
    # step 3 find the top, bottom, left and right rows and columns
    # given the specified offsets
    #
    col_left, col_right = (image_col + offsets.west, image_col + 
                           offsets.east)
    row_bot, row_top = image_row + offsets.south, image_row + offsets.north
    row_max, col_max = scene_ds.shape
    #
    # step 4: sanity check to make sure we haven't gone outside the image
    #
    if (col_left < 0 or col_right > (col_max -1) or
         row_top < 0 or row_bot > (row_max -1)):
        raise ValueError((f"bounds error with "
                          f"{(col_left, col_right, row_bot, row_top)=}"
                          f"\n{scene_ds.shape=}"))
    #
    # step 5: find the x,y corners of the bounding box
    #
    ul_x, ul_y = full_affine * (col_left, row_top)
    lr_x, lr_y = full_affine * (col_right, row_bot)
    #
    # step 6 create a new instance of the Bounding_box object and return it
    #
    bounds = Bounding_box(xmin=ul_x,ymin=lr_y,xmax=lr_x,ymax=ul_y)
    return bounds
    
