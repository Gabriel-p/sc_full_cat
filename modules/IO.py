
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

    # Read info on DBs files
    dbs = dict(in_params.items('DBs column names'))
    dbs_dict = {}
    for key, value in dbs.items():
        temp = []
        for elem in value.split():
            if elem == 'None':
                elem = None
            if elem == 'False':
                elem = False
            if elem == 'True':
                elem = True
            temp.append(elem)
        dbs_dict[key.upper()] = temp

    cmpars = in_params['Cross-match parameters']
    max_sep, defFlag = cmpars.getfloat('max_sep'), cmpars.getboolean('defFlag')

    outpars = in_params['Output parameters']
    plotFlag, dpi, mode = outpars.getboolean('plotFlag'),\
        outpars.getfloat('dpi'), outpars.get('mode')

    if mode not in ('z_dist', 'd_dist'):
        raise ValueError("Unrecognized 'mode' value: {}".format(mode))

    return dbs_dict, max_sep, defFlag, plotFlag, dpi, mode


def readData(dbs_dict):
    """
    Read *all* the databases in 'input/ folder. Prepare so that all have the
    required columns with matching names.
    """

    files = os.listdir('input/')
    files = [_.replace('.dat', '').upper() for _ in files]

    db_dict = {}
    for db_name in files:
        try:
            db = readDB(db_name, *dbs_dict[db_name])
            db_dict[db_name] = db
        except KeyError:
            pass

    return db_dict


def readDB(
    db_name, c_name, c_ra, c_de, c_dist, c_age, c_feh, ra_h, xy_type,
        dist_type, age_years, hs=0):
    """
    Read a database and return a properly formatted table
    """
    db = ascii.read(
        'input/' + db_name + '.dat', delimiter=',', header_start=int(hs))

    db[c_name].name = 'name'
    db['name'] = db['name'].astype(str)

    db[c_ra].name = 'ra'
    db[c_de].name = 'dec'
    # Add units to (ra, dec)
    if xy_type == 'equat':
        if ra_h is True:
            eq = SkyCoord(
                ra=db['ra'], dec=db['dec'], unit=(u.hour, u.deg), frame='icrs')
        else:
            eq = SkyCoord(
                ra=db['ra'], dec=db['dec'], unit=(u.deg, u.deg), frame='icrs')
    elif xy_type == 'galac':
        eq = SkyCoord(
            frame="galactic", l=db['ra'], b=db['dec'], unit=(u.deg, u.deg))

    if xy_type == 'equat':
        db['ra'] = eq.ra
        db['dec'] = eq.dec
    elif xy_type == 'galac':
        db['ra'] = eq.icrs.ra
        db['dec'] = eq.icrs.dec

    if c_dist is not None:
        db[c_dist].name = 'dist_pc'
        db['dist_pc'] = db['dist_pc'].astype(float)
        if dist_type == 'mod':
            db['dist_pc'] = 10**(.2 * (db['dist_pc'] + 5))
        elif dist_type == 'plx':
            db['dist_pc'][db['dist_pc'] <= 0] = np.nan
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

    no_dist = np.isnan(db['dist_pc']).sum()
    print("  {}, N_tot={}, N_nodist={}".format(db_name, len(db), no_dist))

    return db


def write2File(crossMdata):
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

    col_order = ['name', 'DBs', 'lon_d', 'lat_d', 'Nm', 'dist_pc', 'x_pc',
                 'y_pc', 'z_pc', 'ra_all', 'dec_all', 'dist_all']
    ascii.write(
        crossMdata[col_order], 'output/crossMdata.dat', format='fixed_width',
        formats={'Nm': '%5.0f',
                 'dist_pc': '%10.0f', 'lon_d': '%10.4f', 'lat_d': '%10.4f',
                 'x_pc': '%10.2f', 'y_pc': '%10.2f', 'z_pc': '%10.2f'},
        overwrite=True)
    print("\nCross-matched data written to file.")
