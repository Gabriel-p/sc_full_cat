
import os
import configparser
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def readINI():
    """
    Read .ini config file
    """
    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    cmpars = in_params['Cross-match parameters']
    max_sep, defFlag = cmpars.getfloat('max_sep'), cmpars.getboolean('defFlag')

    outpars = in_params['Output parameters']
    plotFlag, dpi, mode = outpars.getboolean('plotFlag'),\
        outpars.getfloat('dpi'), outpars.get('mode')

    if mode not in ('z_dist', 'd_dist'):
        raise ValueError("Unrecognized 'mode' value: {}".format(mode))

    return max_sep, defFlag, plotFlag, dpi, mode


def readData():
    """
    Read *all* the databases in 'input/ folder. Prepare so that all have the
    required columns with matching names.
    """
    def oc_read():
        # OPENCLUST - New Optically Visible Open Clusters and Candidates
        # Catalog
        # https://heasarc.gsfc.nasa.gov/W3Browse/all/openclust.html
        # Table with proper format from:
        # http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/txt?B/ocl
        #
        # Manual edits: Commented the header, added a column separator between
        # ra and dec.
        openclst = ascii.read('input/OPENCLUST.dat')
        openclst['col1'].name = 'name'
        openclst['col2'].name = 'ra'
        openclst['col3'].name = 'dec'
        # Add units to (ra, dec)
        eq = SkyCoord(
            ra=openclst['ra'], dec=openclst['dec'], unit=(u.hour, u.deg),
            frame='icrs')
        openclst['ra'] = eq.ra
        openclst['dec'] = eq.dec
        openclst['col6'].name = 'dist_pc'
        openclst['dist_pc'] = openclst['dist_pc'].astype(float)
        openclst['col8'].name = 'log_age'
        openclst['col18'].name = 'fe_h'
        openclst['col13'].name = 'N_m'

        return openclst

    def mwsc_read():
        # MWSC - Milky Way Star Clusters Catalog
        # https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html
        #
        # Manual edits: removed initial '|' chars.
        mwsc = ascii.read('input/MWSC.dat', format='fixed_width')
        mwsc['distance'].name = 'dist_pc'
        mwsc['dist_pc'] = mwsc['dist_pc'].astype(float)
        mwsc['metallicity'].name = 'fe_h'
        # Add units to (ra, dec)
        eq = SkyCoord(
            ra=mwsc['ra'], dec=mwsc['dec'], unit=(u.hour, u.deg),
            frame='icrs')
        mwsc['ra'] = eq.ra
        mwsc['dec'] = eq.dec

        # Only use objects classified as open clusters.
        op_msk = mwsc['class'] == 'OPEN STAR CLUSTER'
        mwsc = mwsc[op_msk]
        # m10 = abs(mwsc['z_pc']) > 10000
        # print(mwsc[m10])

        return mwsc

    def camargo_read():
        # Camargo et al (2010-2013); Table 1, 2
        # http://vizier.u-strasbg.fr/viz-bin/
        # VizieR?-source=J%2FMNRAS%2F432%2F3349
        #
        # Manual edits: stitched together both tables.
        camargo = ascii.read('input/Camargo.dat')
        camargo['Cluster'].name = 'name'
        camargo['RAJ2000'].name = 'ra'
        camargo['DEJ2000'].name = 'dec'
        # Add units to (ra, dec)
        eq = SkyCoord(
            ra=camargo['ra'], dec=camargo['dec'], unit=(u.hour, u.deg),
            frame='icrs')
        camargo['ra'] = eq.ra
        camargo['dec'] = eq.dec
        camargo['Dist'].name = 'dist_pc'
        camargo['dist_pc'] = camargo['dist_pc'] * 1000.
        camargo['Age'].name = 'log_age'
        camargo['log_age'] = np.log10(camargo['log_age'] * 1000000.)
        camargo['fe_h'] = np.array([np.nan for _ in camargo])

        return camargo

    def webda_read():
        # WEBDA - http://www.univie.ac.at/webda/
        #
        # Manual edits: downloaded from the parameters form, removed a line
        # 'Back to WEBDA home page'
        webda = ascii.read('input/WEBDA.dat')
        webda['Cluster_name'].name = 'name'
        webda['RA_2000'].name = 'ra'
        webda['Dec_2000'].name = 'dec'
        # Add units to (ra, dec)
        eq = SkyCoord(
            ra=webda['ra'], dec=webda['dec'], unit=(u.hour, u.deg),
            frame='icrs')
        webda['ra'] = eq.ra
        webda['dec'] = eq.dec
        webda['Age'].name = 'log_age'
        webda['Dist'].name = 'dist_pc'
        webda['dist_pc'] = webda['dist_pc'].astype(float)
        webda['Fe/H'].name = 'fe_h'

        return webda

    files = os.listdir('input/')
    files = [_.replace('.dat', '').upper() for _ in files]
    dbs_dict = {
        'OPENCLUST': oc_read, 'MWSC': mwsc_read, 'CAMARGO': camargo_read,
        'WEBDA': webda_read}

    db_dict = {}
    for db in files:
        db_data = dbs_dict[db]()
        print("  {}, N={}".format(db, len(db_data)))
        db_dict[db] = db_data

    return db_dict


def write2File(allData):
    """
    Write cross-matched data to 'crossMdata.dat' file.
    """
    # Create /output dir if it does not exist.
    if not os.path.exists('output/'):
        os.makedirs('output/')

    crossMdata = allData['crossMdata']

    # Equatorial to degrees (from radians)
    eq = SkyCoord(crossMdata['ra'], crossMdata['dec'])
    crossMdata['ra_h'] = eq.ra.to_string(unit=u.hour)
    crossMdata['dec_d'] = eq.dec.to_string(unit=u.degree)
    # Galactic to degrees (from radians)
    gl = SkyCoord(l=crossMdata['lon'], b=crossMdata['lat'], frame='galactic')
    crossMdata['lon_d'] = gl.l.wrap_at(360 * u.deg)
    crossMdata['lat_d'] = gl.b

    # Order by 'ra' (if I attempt to order the table as is, a ValueError
    # is raised) TODO
    bb = Table([[_ for _ in range(len(crossMdata))], crossMdata['ra']])
    bb.sort('ra')
    mask = bb['col0'].data
    crossMdata = crossMdata[mask]

    dtBs_names = list(allData.keys())
    dtBs_names.remove('crossMdata')
    col_order = ['name', 'ra_h', 'dec_d'] + dtBs_names +\
        ['N_d', 'dist_pc', 'lon_d', 'lat_d', 'x_pc', 'y_pc', 'z_pc']
    ascii.write(
        crossMdata[col_order], 'output/crossMdata.dat', format='fixed_width',
        formats={'N_d': '%5.0f',
                 'dist_pc': '%10.2f', 'lon_d': '%10.4f', 'lat_d': '%10.4f',
                 'x_pc': '%10.2f', 'y_pc': '%10.2f', 'z_pc': '%10.2f'},
        overwrite=True)
    print("Cross-matched data written to file.")
