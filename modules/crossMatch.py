
# import warnings
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def match(allDatabases, max_sep):
    """
    """
    print("\nPerform cross-match (max_sep={})".format(max_sep))

    # crossMdata = nameMatch(allDatabases)
    r_max_deg = (max_sep / 3600) * u.deg

    # Initial Table with dummy entry
    DB_names = list(allDatabases.keys())
    # Eq coordinates of (l, b) = (0, 90)
    ra_l_0, dec_b_90 = 192.85948121 * u.deg, 27.12825118 * u.deg
    crossMdata = Table(
        [['dummy'], ['noname'], [ra_l_0], [dec_b_90]],
        names=('IDs', 'name', 'ra', 'dec'))

    # Process all DBs
    for DB_name in DB_names:
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

            # If no new matches were found, break out
            N_new = (d2d < r_max_deg).sum()
            if N_new == 0 or N_old == N_new:
                break
            N_old = N_new

            # Sort by smallest distance first
            di = np.argsort(d2d)
            for i in di:
                if d2d[i] < r_max_deg:
                    if idx2[i] not in idx2_unq:
                        # breakpoint()
                        # This is an acceptable match
                        idx1_unq.append(c1_ids[i])
                        idx2_unq.append(c2_ids[idx2[i]])

            # Remove the matched elements from the catalogs
            cm_ra[idx1_unq] = ra_l_0
            cm_de[idx1_unq] = dec_b_90
            da_ra[idx2_unq] = ra_l_0
            da_de[idx2_unq] = dec_b_90

        # Indexes for elements with no match found.
        idx1_ncm = [_ for _ in c1_ids if _ not in idx1_unq]
        idx2_ncm = [_ for _ in c2_ids if _ not in idx2_unq]

        crossMdata = updtCrossMData(
            crossMdata, data, DB_name, idx1_unq, idx2_unq, idx1_ncm, idx2_ncm)

    return crossMdata


def updtCrossMData(
        crossMdata, data, DB_name, idx1_unq, idx2_unq, idx1_ncm, idx2_ncm):
    """
    """
    idx_m, name_m, ra_m, dec_m = [], [], [], []

    # For each crossed-match element
    for j, i2 in enumerate(idx2_unq):
        i1 = idx1_unq[j]
        # Store averaged (ra, dec) values.
        ra_m.append(np.mean([crossMdata['ra'][i1], data['ra'][i2]]))
        dec_m.append(np.mean([crossMdata['dec'][i1], data['dec'][i2]]))
        idx_m.append(crossMdata['IDs'][i1] + '; ' + str(i2) + ' ' + DB_name)
        # name_m.append(
        #     crossMdata['name'][i1].lower().replace('_', '').replace(' ', '')\
        #     + ';' + data['name'][i2].lower().replace('_', '').replace(' ', ''))

    name_nm, ra_nm, dec_nm, idx_nm, dist_nm, Nm_nm = [[] for _ in range(6)]
    for i1 in idx1_ncm:
        ra_nm.append(crossMdata['ra'][i1])
        dec_nm.append(crossMdata['dec'][i1])
        idx_nm.append(crossMdata['IDs'][i1] + '; -- ' + DB_name)
        # name_nm.append(crossMdata['name'][i1].lower().replace('_', '').replace(' ', '') + ';')
    for i2 in idx2_ncm:
        ra_nm.append(data['ra'][i2])
        dec_nm.append(data['dec'][i2])
        idx_nm.append(str(i2) + ' ' + DB_name)
        # name_nm.append(data['name'][i2].lower().replace('_', '').replace(' ', '') + ';')

    # Combine matched and not matched data.
    ra_cmb = Column(ra_m + ra_nm, name='ra', unit=u.degree)
    dec_cmb = Column(dec_m + dec_nm, name='dec', unit=u.degree)
    idx_cmb = Column(idx_m + idx_nm, name='IDs')
    # name_cmb = Column(name_m + name_nm, name='name')

    # Update final database.
    # crossMdata = Table([idx_cmb, name_cmb, ra_cmb, dec_cmb])
    crossMdata = Table([idx_cmb, ra_cmb, dec_cmb])

    return crossMdata


# def nameMatch(allDatabases):
#     """
#     """

#     # Initial Table
#     DB_names = list(allDatabases.keys())
#     crossMdata = allDatabases[DB_names[0]]

#     idx_names = {}
#     for DB_name in DB_names[1:]:
#         data = allDatabases[DB_name]

#         cm_n = [_.replace(' ', '').replace('_', '').lower() for _ in crossMdata['name']]
#         da_n = [_.replace(' ', '').replace('_', '').lower() for _ in data['name']]

#         idx_cm, idx_da = [], []
#         for i, cl in enumerate(cm_n):
#             try:
#                 j = da_n.index(cl)
#                 idx_cm.append(i)
#                 idx_da.append(j)
#             except ValueError:
#                 pass

#         idx_names[DB_name] = (idx_cm, idx_da)

#     breakpoint()

#     new_names = []
#     for i, cl in enumerate(crossMdata['name']):
#         new_names.append(cl + ' (' + DB_names[0] + ')')
#     crossMdata['name'] = new_names

#     # for DB_name, idxs in idx_names.items():

#     return crossMdata
