
# import warnings
import re
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def match(allDatabases, max_sep):
    """
    """

    # crossMdata = nameMatch(allDatabases)
    r_max_deg = (max_sep / 3600) * u.deg

    # Initial Table
    DB_names = list(allDatabases.keys())
    crossMdata = allDatabases[DB_names[0]]
    print("\nInitial database: {}".format(DB_names[0]))

    new_names = []
    for i, cl in enumerate(crossMdata['name']):
        new_names.append(cl + ' (' + DB_names[0] + ')')
    crossMdata['name'] = new_names

    crossMdata['N_m'] = np.ones(len(crossMdata))

    # Process the remaining DBs
    for DB_name in DB_names[1:]:
        print("  (processing {})".format(DB_name))

        data = allDatabases[DB_name]

        # Initial full list of observed and queried catalogs.
        c1_ids = np.arange(len(crossMdata))
        c2_ids = np.arange(len(data))

        cm_ra, cm_de = crossMdata['ra'].value * u.deg,\
            crossMdata['dec'].value * u.deg
        da_ra, da_de = data['ra'].value * u.deg, data['dec'].value * u.deg

        # Store the indexes in 'c1' and c2'
        idx1_unq, idx2_unq, N_old = [], [], 0
        while True:

            # Define catalogs to be matched.
            c1 = SkyCoord(ra=cm_ra, dec=cm_de)
            c2 = SkyCoord(ra=da_ra, dec=da_de)

            # 'idx2' are indices into 'c2' that are the closest objects to each
            # of the coordinates in 'c1'.
            idx2, d2d, _ = c1.match_to_catalog_sky(c2)

            N_new = (d2d < r_max_deg).sum()
            if N_new <= 0 or N_old == N_new:
                break
            N_old = N_new

            # Sort by smallest distance first
            di = np.argsort(d2d)
            for i in di:
                if d2d[i] < r_max_deg:
                    if idx2[i] not in idx2_unq:
                        # This is an acceptable match
                        idx1_unq.append(c1_ids[i])
                        idx2_unq.append(c2_ids[idx2[i]])

            # Remove the matched elements from the catalogs
            cm_ra[idx1_unq] = 0.
            cm_de[idx1_unq] = 0.
            da_ra[idx2_unq] = 0.
            da_de[idx2_unq] = 0.

        # Indexes for elements with no match found.
        idx1_ncm = [_ for _ in c1_ids if _ not in idx1_unq]
        idx2_ncm = [_ for _ in c2_ids if _ not in idx2_unq]

        crossMdata = updtCrossMData(
            crossMdata, data, DB_name, idx1_unq, idx2_unq, idx1_ncm, idx2_ncm)

    # Store mean distance values
    crossMdata['dist_pc'] = crossMdata['dist_pc'] / crossMdata['N_m']

    print("Databases cross-matched")
    for Nm in np.arange(crossMdata['N_m'].max(), 0, -1):
        if Nm == 1:
            txt = 'No cross-match found'
        else:
            txt = '{} matches found'.format(Nm)
        print("  {}: {}".format(txt, (crossMdata['N_m'] == Nm).sum()))

    print("\nUnique clusters in DBS")
    all_dbs = []
    for i, cl in enumerate(crossMdata['name']):
        if crossMdata['N_m'][i] == 1:
            dbs = re.findall('\(.*?\)',cl)
            dbs = [_.replace('(', '').replace(')', '') for _ in dbs]
            all_dbs += dbs
    all_dbs = np.array(all_dbs)
    for db in DB_names:
        msk = db == all_dbs
        print("{}: {}".format(db, msk.sum()))

    return crossMdata


def nameMatch(allDatabases):
    """
    """

    # Initial Table
    DB_names = list(allDatabases.keys())
    crossMdata = allDatabases[DB_names[0]]

    idx_names = {}
    for DB_name in DB_names[1:]:
        data = allDatabases[DB_name]

        cm_n = [_.replace(' ', '').replace('_', '').lower() for _ in crossMdata['name']]
        da_n = [_.replace(' ', '').replace('_', '').lower() for _ in data['name']]

        idx_cm, idx_da = [], []
        for i, cl in enumerate(cm_n):
            try:
                j = da_n.index(cl)
                idx_cm.append(i)
                idx_da.append(j)
            except ValueError:
                pass

        idx_names[DB_name] = (idx_cm, idx_da)

    breakpoint()

    new_names = []
    for i, cl in enumerate(crossMdata['name']):
        new_names.append(cl + ' (' + DB_names[0] + ')')
    crossMdata['name'] = new_names

    # for DB_name, idxs in idx_names.items():


    return crossMdata


def updtCrossMData(
        crossMdata, data, DB_name, idx1_unq, idx2_unq, idx1_ncm, idx2_ncm):
    """
    """
    name_m, ra_m, dec_m, dist_m, Nm_m = [], [], [], [], []

    # For each crossed-match element
    for j, i2 in enumerate(idx2_unq):

        i1 = idx1_unq[j]

        # Store all names.
        name_m.append(crossMdata['name'][i1] + ', ' + data['name'][i2]
                      + ' (' + DB_name + ')')
        # Store averaged (ra, dec) values.
        ra_m.append(np.mean([crossMdata['ra'][i1], data['ra'][i2]]))
        dec_m.append(np.mean([crossMdata['dec'][i1], data['dec'][i2]]))
        # Sum individual distances
        dist_m.append(np.nansum([
            crossMdata['dist_pc'][i1], data['dist_pc'][i2]]))
        Nm_m.append(crossMdata['N_m'][i1] + 1.)

    name_nm, ra_nm, dec_nm, idx_nm, dist_nm, Nm_nm = [], [], [], [], [], []
    for i1 in idx1_ncm:
        name_nm.append(crossMdata['name'][i1])
        ra_nm.append(crossMdata['ra'][i1])
        dec_nm.append(crossMdata['dec'][i1])
        idx_nm.append('--')
        dist_nm.append(crossMdata['dist_pc'][i1])
        Nm_nm.append(crossMdata['N_m'][i1])

    for i2 in idx2_ncm:
        name_nm.append(data['name'][i2] + ' (' + DB_name + ')')
        ra_nm.append(data['ra'][i2])
        dec_nm.append(data['dec'][i2])
        idx_nm.append(i2)
        dist_nm.append(data['dist_pc'][i2])
        Nm_nm.append(1.)

    # Combine matched and not matched data.
    n_cmb = Column(name_m + name_nm, name='name')
    ra_cmb = Column(ra_m + ra_nm, name='ra', unit=u.degree)
    dec_cmb = Column(dec_m + dec_nm, name='dec', unit=u.degree)
    dist_cmb = Column(dist_m + dist_nm, name='dist_pc', unit=u.pc)
    Nm_cmb = Column(Nm_m + Nm_nm, name='N_m')
    tempData = Table([n_cmb, ra_cmb, dec_cmb, dist_cmb, Nm_cmb])

    # Update final database.
    crossMdata = tempData

    return crossMdata
