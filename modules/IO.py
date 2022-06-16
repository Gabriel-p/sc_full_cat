
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

    OPENCLUST - New Optically Visible Open Clusters and Candidates Catalog
    https://heasarc.gsfc.nasa.gov/W3Browse/all/openclust.html
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=B/ocl

    MWSC - Milky Way Star Clusters Catalog
    https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html
    Manual edits: removed initial '|' chars.

    """

    dbs_dict = {
        'BOSSINI': (' ', 'Cluster', 'RA_ICRS', 'DE_ICRS', 'Dist_mod', 'logA',
                    'Fe/H', False, True),
        'CAMARGO': ('|', 'Cluster', 'RAJ2000', 'DEJ2000', 'Dist', 'Age',
                    None, True, False, True),
        'CANTAT_GAUDIN': (',', 'Cluster', 'RA_ICRS', 'DE_ICRS', 'DistPc',
                          'AgeNN', None),
        'CANTAT_GAUDIN_2018': (' ', 'Cluster', 'RAJ2000', 'DEJ2000', 'dmode',
                               None, None),
        'DIAS': (' ', 'Cluster', 'RA_ICRS', 'DE_ICRS', 'Dist', 'logage',
                 '[Fe/H]'),
        # 'LIU_PANG': (';', 'Name', 'GLON', 'GLAT', 'Dist', 'Age',
        #               '[Fe/H]', True),
        'OPENCLUST': (';', 'Cluster', 'RAJ2000', 'DEJ2000', 'Dist', 'Age',
                      '[Fe/H]', True),
        'MWSC': ('|', 'name', 'ra', 'dec', 'distance', 'log_age',
                 'metallicity', True),
        'WEBDA': (',', 'Cluster_name', 'RA_2000', 'Dec_2000', 'Dist', 'Age',
                  'Fe/H', True)
    }

    # def bossini_read():
    #     """
    #     Bossini et al. (2019)
    #     https://ui.adsabs.harvard.edu/abs/2019A%26A...623A.108B

    #     Source (Vizier):
    #     https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A%2bA/623/A108/
    #     """
    #     db = ascii.read('input/BOSSINI.dat')
    #     db['Cluster'].name = 'name'
    #     db['RA_ICRS'].name = 'ra'
    #     db['DE_ICRS'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    #     db['ra'] = eq.ra
    #     db['dec'] = eq.dec
    #     db['logA'].name = 'log_age'
    #     db['Dist_mod'].name = 'dist_pc'
    #     # To parsec
    #     db['dist_pc'] = 10**(.2 * (db['dist_pc'] + 5))
    #     db['dist_pc'] = db['dist_pc'].astype(float)
    #     db['Fe/H'].name = 'fe_h'

    #     return db

    # def oc_read():
    #     openclst = ascii.read('input/OPENCLUST.dat', delimiter=';')
    #     openclst['Cluster'].name = 'name'
    #     openclst['RAJ2000'].name = 'ra'
    #     openclst['DEJ2000'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=openclst['ra'], dec=openclst['dec'], unit=(u.hour, u.deg),
    #         frame='icrs')
    #     openclst['ra'] = eq.ra
    #     openclst['dec'] = eq.dec
    #     openclst['Dist'].name = 'dist_pc'
    #     openclst['dist_pc'] = openclst['dist_pc'].astype(float)
    #     openclst['Age'].name = 'log_age'
    #     openclst['[Fe/H]'].name = 'fe_h'

    #     return openclst

    # def mwsc_read():
    #     mwsc = ascii.read('input/MWSC.dat', format='fixed_width')
    #     mwsc['distance'].name = 'dist_pc'
    #     mwsc['dist_pc'] = mwsc['dist_pc'].astype(float)
    #     mwsc['metallicity'].name = 'fe_h'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=mwsc['ra'], dec=mwsc['dec'], unit=(u.hour, u.deg),
    #         frame='icrs')
    #     mwsc['ra'] = eq.ra
    #     mwsc['dec'] = eq.dec

    #     # Only use objects classified as open clusters.
    #     op_msk = mwsc['class'] == 'OPEN STAR CLUSTER'
    #     mwsc = mwsc[op_msk]
    #     # m10 = abs(mwsc['z_pc']) > 10000
    #     # print(mwsc[m10])

    #     return mwsc

    # def camargo_read():
    #     # Camargo et al (2010-2013); Table 1, 2
    #     # http://vizier.u-strasbg.fr/viz-bin/
    #     # VizieR?-source=J%2FMNRAS%2F432%2F3349
    #     #
    #     # Manual edits: stitched together both tables.
    #     camargo = ascii.read('input/Camargo.dat')
    #     camargo['Cluster'].name = 'name'
    #     camargo['RAJ2000'].name = 'ra'
    #     camargo['DEJ2000'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=camargo['ra'], dec=camargo['dec'], unit=(u.hour, u.deg),
    #         frame='icrs')
    #     camargo['ra'] = eq.ra
    #     camargo['dec'] = eq.dec
    #     camargo['Dist'].name = 'dist_pc'
    #     camargo['dist_pc'] = camargo['dist_pc'] * 1000.
    #     camargo['Age'].name = 'log_age'
    #     camargo['log_age'] = np.log10(camargo['log_age'] * 1000000.)
    #     camargo['fe_h'] = np.array([np.nan for _ in camargo])

    #     return camargo

    # def webda_read():
    #     # WEBDA - http://www.univie.ac.at/webda/
    #     #
    #     # Manual edits: downloaded from the parameters form, removed a line
    #     # 'Back to WEBDA home page'
    #     webda = ascii.read('input/WEBDA.dat')
    #     webda['Cluster_name'].name = 'name'
    #     webda['RA_2000'].name = 'ra'
    #     webda['Dec_2000'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=webda['ra'], dec=webda['dec'], unit=(u.hour, u.deg),
    #         frame='icrs')
    #     webda['ra'] = eq.ra
    #     webda['dec'] = eq.dec
    #     webda['Age'].name = 'log_age'
    #     webda['Dist'].name = 'dist_pc'
    #     webda['dist_pc'] = webda['dist_pc'].astype(float)
    #     webda['Fe/H'].name = 'fe_h'

    #     return webda

    # def cantat2018_read():
    #     """
    #     Cantat-Gaudin et al. (2018)

    #     """
    #     db = ascii.read('input/CANTAT_GAUDIN_2018.dat')
    #     db['Cluster'].name = 'name'
    #     db['RAJ2000'].name = 'ra'
    #     db['DEJ2000'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    #     db['ra'] = eq.ra
    #     db['dec'] = eq.dec
    #     # db['AgeNN'].name = 'log_age'
    #     db['dmode'].name = 'dist_pc'
    #     db['dist_pc'] = db['dist_pc'].astype(float)
    #     # db['Fe/H'].name = 'fe_h'

    #     return db

    # def cantat_read():
    #     """
    #     Cantat-Gaudin et al. (2020)
    #     https://ui.adsabs.harvard.edu/abs/2020A%26A...640A...1C

    #     Source (Vizier):
    #     https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/640/A1
    #     """
    #     db = ascii.read('input/CANTAT_GAUDIN.dat')
    #     db['Cluster'].name = 'name'
    #     db['RA_ICRS'].name = 'ra'
    #     db['DE_ICRS'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    #     db['ra'] = eq.ra
    #     db['dec'] = eq.dec
    #     db['AgeNN'].name = 'log_age'
    #     db['DistPc'].name = 'dist_pc'
    #     db['dist_pc'] = db['dist_pc'].astype(float)
    #     # db['Fe/H'].name = 'fe_h'

    #     return db

    # def dias_read():
    #     """
    #     Dias et al. (2021)
    #     https://ui.adsabs.harvard.edu/abs/2021MNRAS.504..356D

    #     Source (Vizier):
    #     https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/504/356
    #     """
    #     db = ascii.read('input/DIAS.dat')
    #     db['Cluster'].name = 'name'
    #     db['RA_ICRS'].name = 'ra'
    #     db['DE_ICRS'].name = 'dec'
    #     # Add units to (ra, dec)
    #     eq = SkyCoord(
    #         ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    #     db['ra'] = eq.ra
    #     db['dec'] = eq.dec
    #     db['logage'].name = 'log_age'
    #     db['Dist'].name = 'dist_pc'
    #     db['dist_pc'] = db['dist_pc'].astype(float)
    #     db['[Fe/H]'].name = 'fe_h'

    #     return db

    files = os.listdir('input/')
    files = [_.replace('.dat', '').upper() for _ in files]

    db_dict = {}
    for db_name in files:
        # print(db_name)
        db = readDB(db_name, *dbs_dict[db_name])
        db_dict[db_name] = db

    return db_dict


def readDB(
    db_name, deli, c_name, c_ra, c_de, c_dist, c_age, c_feh,
        ra_h=False, dist_mod=False, age_years=False, hs=0):
    """
    Read a database and return a properly formatted table
    """
    db = ascii.read(
        'input/' + db_name + '.dat', delimiter=deli, header_start=hs)
    db[c_name].name = 'name'
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

    db[c_dist].name = 'dist_pc'
    db['dist_pc'] = db['dist_pc'].astype(float)
    if dist_mod is True:
        db['dist_pc'] = 10**(.2 * (db['dist_pc'] + 5))

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

    try:
        no_dist = db['dist_pc'].mask.sum()
    except AttributeError:
        no_dist = 0
    print("  {}, N_tot={}, N_nodist={}".format(db_name, len(db), no_dist))

    return db


def write2File(crossMdata, dtBs_names):
    """
    Write cross-matched data to 'crossMdata.dat' file.
    """
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

    # + dtBs_names +\
    col_order = ['name', 'ra_h', 'dec_d'] +\
        ['N_m', 'dist_pc', 'lon_d', 'lat_d', 'x_pc', 'y_pc', 'z_pc']
    ascii.write(
        crossMdata[col_order], 'output/crossMdata.dat', format='fixed_width',
        formats={'N_m': '%5.0f',
                 'dist_pc': '%10.2f', 'lon_d': '%10.4f', 'lat_d': '%10.4f',
                 'x_pc': '%10.2f', 'y_pc': '%10.2f', 'z_pc': '%10.2f'},
        overwrite=True)
    print("Cross-matched data written to file.")
