
from astropy import units as u
# from astropy.coordinates.distances import Distance
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import numpy as np
# from difflib import SequenceMatcher

from modules.IO import readINI, readData, write2File
from modules import makePlot, crossMatch


def main():
    """
    """
    max_sep, defFlag, plotFlag, dpi, mode = readINI()

    # Read databases.
    print("\nReading databases:")
    allDatabases = readData()

    # Cross-match all databases.
    if len(allDatabases.keys()) > 1:
        print("Perform cross-match (max_sep={})".format(max_sep))
        allData = crossMatch.match(allDatabases, max_sep)
    else:
        raise ValueError("At least two databases must be present in 'input/")


    import pdb; pdb.set_trace()  # breakpoint faa7716a //

    # Define Galactocentric frame
    gc_frame = frameGalactocentric(defFlag)
    # Add Cartesian data.
    allData = dist2plane(allData, gc_frame)

    # Write output file.
    write2File(allData)

    if plotFlag:
        print("Plot mode: {}".format(mode))
        makePlot.plot(dpi, mode, allData['crossMdata'], gc_frame)


def frameGalactocentric(defFlag=True):
    """
    Transform to Galactocentric coordinate
    http://docs.astropy.org/en/stable/api/
         astropy.coordinates.Galactocentric.html

    defFlag: use astropy's default values.

    """
    if defFlag:
        # Default Galactic Center is 8.3 kpc (Gillessen et al. 2009)
        return coord.Galactocentric()
    else:
        # Sun's distance to galactic center from Camargo et al (2013)
        # (taken from Bica et al. 2006)
        return coord.Galactocentric(galcen_distance=7.2 * u.kpc)


# def similar(a, b):
#     return 1. - SequenceMatcher(None, a, b).ratio()


def dist2plane(allData, gc_frame):
    """
    Convert equatorial coordinates, and obtain the
    Cartesian coordinates with 'z_pc' the vertical distance.
    """
    for data in allData.values():

        eq = SkyCoord(ra=data['ra'], dec=data['dec'], frame='icrs')
        lb = eq.transform_to('galactic')
        data['lon'] = lb.l.wrap_at(180 * u.deg).radian * u.radian
        data['lat'] = lb.b.radian * u.radian

        try:
            data['dist_pc'] = data['dist_pc'].filled(np.nan)
        except AttributeError:
            pass

        # Galactic coordinates.
        coords = SkyCoord(
            l=data['lon'], b=data['lat'], distance=data['dist_pc'] * u.pc,
            frame='galactic')
        # Galactocentric coordinates.
        c_glct = coords.transform_to(gc_frame)
        data['x_pc'], data['y_pc'], data['z_pc'] = c_glct.x, c_glct.y, c_glct.z

    return allData


if __name__ == '__main__':
    main()
