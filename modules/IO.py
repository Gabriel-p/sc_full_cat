
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

    WEBDA - https://webda.physics.muni.cz/name_selection.html
    To download full list go to the link above and submit an empty query.

    OPENCLUST - New Optically Visible Open Clusters and Candidates Catalog
    https://heasarc.gsfc.nasa.gov/W3Browse/all/openclust.html
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=B/ocl

    MWSC - Milky Way Star Clusters Catalog
    https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html
    Manual edits: removed initial '|' chars.

    """

    dbs_dict = {
        # 'BOSSINI': (' ', 'Cluster', 'RA_ICRS', 'DE_ICRS', 'Dist_mod', 'logA',
        #             'Fe/H', False, 'mod'),
        # 'CAMARGO': ('|', 'Cluster', 'RAJ2000', 'DEJ2000', 'Dist', 'Age',
        #             None, True, 'pc', True),
        # 'CG2018': (',', 'Cluster', 'RAJ2000', 'DEJ2000', 'dmode',
        #            None, None),
        'CG20': ('Name', 'RA', 'DEC', 'D_pc', 'Age', None),
        'WEBDA': ('Cluster_name', 'RA_2000', 'Dec_2000', 'Dist', 'Age',
                  'Fe/H', True),
        'DIAS21': ('Cluster', 'RA_ICRS', 'DE_ICRS', 'Dist', 'logage',
                   '[Fe/H]'),
        'OPENCLUST': ('Cluster', 'RAJ2000', 'DEJ2000', 'Dist', 'Age',
                      '[Fe/H]', True),
        'MWSC': ('name', 'ra', 'dec', 'distance', 'log_age',
                 'metallicity', True),
        'BICA19': ('Name', 'RAJ2000', 'DEJ2000', None, None, None, True),
        'LIUPANG': ('ID', '_RA.icrs', '_DE.icrs', 'plx', 'Age',
                    None, False, 'plx')
    }

    files = os.listdir('input/')
    files = [_.replace('.dat', '').upper() for _ in files]

    db_dict = {}
    for db_name in files:
        # print(db_name)
        db = readDB(db_name, *dbs_dict[db_name])
        db_dict[db_name] = db

    return db_dict


def readDB(
    db_name, c_name, c_ra, c_de, c_dist, c_age, c_feh,
        ra_h=False, dist_type='pc', age_years=False, hs=0):
    """
    Read a database and return a properly formatted table
    """
    db = ascii.read(
        'input/' + db_name + '.dat', delimiter=',', header_start=hs)

    db[c_name].name = 'name'
    db['name'] = db['name'].astype(str)
    if db_name == 'BICA19':
        for i, nm in enumerate(db['name']):
            db['name'][i] = nm.split(',')[0]

    db[c_ra].name = 'ra'
    db[c_de].name = 'dec'
    # Add units to (ra, dec)
    if ra_h is True:
        eq = SkyCoord(
            ra=db['ra'], dec=db['dec'], unit=(u.hour, u.deg), frame='icrs')
    else:
        eq = SkyCoord(
            ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    db['ra'] = eq.ra
    db['dec'] = eq.dec

    if c_dist is not None:
        db[c_dist].name = 'dist_pc'
        db['dist_pc'] = db['dist_pc'].astype(float)
        if dist_type == 'mod':
            db['dist_pc'] = 10**(.2 * (db['dist_pc'] + 5))
        elif dist_type == 'plx':
            db['dist_pc'] = 1000 / db['dist_pc']

        try:
            db['dist_pc'].fill_value = np.nan
            db['dist_pc'] = db['dist_pc'].filled()
        except AttributeError:
            # Column is not masked
            pass
    else:
        db['dist_pc'] = np.nan * np.ones(len(db))

    if c_age is not None:
        db[c_age].name = 'log_age'
        if age_years is True:
            db['log_age'] = np.log10(db['log_age'] * 1000000.)

    if c_feh is not None:
        db[c_feh].name = 'fe_h'

    if db_name == 'MWSC':
        # Only use objects classified as open clusters.
        msk = db['class'] == 'OPEN STAR CLUSTER'
        db = db[msk]

    if db_name == 'BICA19':
        # Only use objects classified as open clusters.
        msk = (db['Class1'] == 'OC') #| (db['Class1'] == 'OCC')
        db = db[msk]

    no_dist = np.isnan(db['dist_pc']).sum()
    print("  {}, N_tot={}, N_nodist={}".format(db_name, len(db), no_dist))

    return db


def write2File(crossMdata, dtBs_names):
    """
    Write cross-matched data to 'crossMdata.dat' file.
    """
    # Equatorial to degrees (from radians)
    # eq = SkyCoord(crossMdata['ra'], crossMdata['dec'])
    # crossMdata['ra_h'] = eq.ra.to_string(unit=u.hour)
    # crossMdata['dec_d'] = eq.dec.to_string(unit=u.degree)
    # Galactic to degrees (from radians)
    gl = SkyCoord(l=crossMdata['lon'], b=crossMdata['lat'], frame='galactic')
    crossMdata['lon_d'] = gl.l.wrap_at(360 * u.deg)
    crossMdata['lat_d'] = gl.b

    # Order by 'lon' (if I attempt to order the table as is, a ValueError
    # is raised) TODO
    bb = Table([[_ for _ in range(len(crossMdata))], crossMdata['lon_d']])
    bb.sort('lon_d')
    mask = bb['col0'].data
    crossMdata = crossMdata[mask]

    col_order = ['name', 'lon_d', 'lat_d', 'N_m', 'dist_pc', 'x_pc', 'y_pc',
                 'z_pc']
    ascii.write(
        crossMdata[col_order], 'output/crossMdata.dat', format='fixed_width',
        formats={'N_m': '%5.0f',
                 'dist_pc': '%10.2f', 'lon_d': '%10.4f', 'lat_d': '%10.4f',
                 'x_pc': '%10.2f', 'y_pc': '%10.2f', 'z_pc': '%10.2f'},
        overwrite=True)
    print("Cross-matched data written to file.")
