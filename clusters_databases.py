
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
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
        # raise ValueError("At least two databases must be present in 'input/")
        print("Perform cross-match (max_sep={})".format(max_sep))
        crossMdata = crossMatch.match(allDatabases, max_sep)
    else:
        # TODO
        crossMdata = allDatabases

    # Define Galactocentric frame
    gc_frame = frameGalactocentric(defFlag)
    # Add Cartesian data.
    crossMdata = dist2plane(allDatabases, crossMdata, gc_frame)

    dtBs_names = list(allDatabases.keys())

    # Write output file.
    write2File(crossMdata, dtBs_names)

    if plotFlag:
        print("Plot mode: {}".format(mode))
        makePlot.plot(dpi, mode, crossMdata, gc_frame)


def frameGalactocentric(defFlag=True):
    """
    Transform to Galactocentric coordinate5
    http://docs.astropy.org/en/stable/api/
         astropy.coordinates.Galactocentric.html

    defFlag: use astropy's default values.

    """
    if defFlag:
        # Default Galactic Center is 8.3 kpc (Gillessen et al. 2009)
        return coord.Galactocentric()
    else:
        # Sun's distance to galactic center from Camargo et al (2013)
        # (taken from Bica et al. 2006)a
        return coord.Galactocentric(galcen_distance=7.2 * u.kpc)


def dist2plane(allDatabases, crossMdata, gc_frame):
    """
   Obtain the Cartesian coordinates
    """
    # Galactic coordinates.
    eq = SkyCoord(ra=crossMdata['ra'], dec=crossMdata['dec'], frame='icrs')
    lb = eq.transform_to('galactic')
    crossMdata['lon'] = lb.l.wrap_at(180 * u.deg).radian * u.radian
    crossMdata['lat'] = lb.b.radian * u.radian
    coords = SkyCoord(
        l=crossMdata['lon'], b=crossMdata['lat'],
        distance=crossMdata['dist_pc'], frame='galactic')

    # Galactocentric coordinates.
    c_glct = coords.transform_to(gc_frame)
    crossMdata['x_pc'], crossMdata['y_pc'], crossMdata['z_pc'] =\
        c_glct.x, c_glct.y, c_glct.z

    return crossMdata


if __name__ == '__main__':
    # Create /output dir if it does not exist.
    if not os.path.exists('output/'):
        os.makedirs('output/')
    main()
