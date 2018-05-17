
from astropy.io import ascii
from astropy import units as u
# from astropy.coordinates.distances import Distance
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import numpy as np

from clusters_databases import frameGalactocentric
from clusters_databases import readData
from clusters_databases import dist2plane


def main(defFlag=True):
    """

    max_sep: match radius in arcseconds.
    plotDBs: generate plots for each database.
    plotCM: generate plots for the cross-matched set.

    """
    # Define Galactocentric frame
    gc_frame = frameGalactocentric(defFlag)

    # Read databases.
    allData = readData('openclust')

    # Add Cartesian data.
    allData = dist2plane(allData, gc_frame)
    openclust = allData['OPENCLUST']

    makePlot(openclust)


def makePlot(openclust):
    z_lim, N_lim = 200, 100

    m_z = openclust['z_pc'] > z_lim
    m_N = openclust[m_z]['N_m'] > N_lim
    N, N1, N2 = len(openclust), len(openclust[m_z]),\
        len(openclust[m_z][m_N])

    plt.subplot(121)
    names = ['all', '>' + str(z_lim) + 'pc', 'N_m>' + str(N_lim)]
    plt.bar(names, [N, N1, N2])

    plt.subplot(122)
    R_GC = np.sqrt(
        openclust[m_z][m_N]['x_pc']**2 + openclust[m_z][m_N]['x_pc']**2)
    plt.xlabel('R_GC (pc)')
    plt.ylabel('z (pc)')
    plt.scatter(
        R_GC, openclust[m_z][m_N]['z_pc'], c=openclust[m_z][m_N]['log_age'])
    cbar = plt.colorbar()
    cbar.set_label('log(age)')

    plt.show()


if __name__ == '__main__':
    main()
