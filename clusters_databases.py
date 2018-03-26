
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np


def main():
    """


    MWSC - Milky Way Star Clusters Catalog
    https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html

    Camargo et al (2010-2013); Table 1, 2
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J%2FMNRAS%2F432%2F3349

    Example plots:

    * http://www.astropy.org/astropy-tutorials/plot-catalog.html
    * http://docs.astropy.org/en/stable/coordinates/skycoord.html
    * http://docs.astropy.org/en/stable/visualization/wcsaxes/overlaying_coordinate_systems.html
    * http://docs.astropy.org/en/stable/visualization/wcsaxes/
    * http://adass.org/adass/proceedings/adass94/greisene.html

    """

    # Read databases.
    openclst = readData()
    import pdb; pdb.set_trace()  # breakpoint 783dfc96 //


    # ra_random = np.random.rand(100)*360.0 * u.degree
    # dec_random = (np.random.rand(100)*180.0-90.0) * u.degree

    # ra_rad = np.random.uniform(-np.pi, np.pi, 100)
    # dec_rad = np.random.uniform(-np.pi / 2., np.pi / 2., 100)

    l_rad = np.random.uniform(-np.pi, np.pi, 100)
    b_rad = np.random.uniform(-np.pi / 2., np.pi / 2., 100)

    plt.figure()
    plt.subplot(111, projection="aitoff")
    plt.grid(True)
    plt.plot(l_rad, b_rad, 'o', markersize=2, alpha=0.3)
    plt.show()


def readData():
    """
    * OPENCLUST - New Optically Visible Open Clusters and Candidates Catalog
      (https://heasarc.gsfc.nasa.gov/W3Browse/all/openclust.html)
      Table with proper format from:
      http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/txt?B/ocl
      Commented the header.
    """

    openclst = ascii.read('input/OPENCLUST.dat')
    openclst['col1'].name = 'name'
    openclst['col2'].name = 'ra_dec'
    openclst['col5'].name = 'dist_pc'
    openclst['col18'].name = 'fe_h'

    return openclst


if __name__ == '__main__':
    main()
