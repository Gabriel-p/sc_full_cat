
# import warnings
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def match(allDatabases, max_sep):
    """
    Cross-match the databases in 'allDatabases', using the 'max_sep' value
    as the limiting match radius in arcsec.

    At each iteration the (ra, dec) values are averaged and these updated
    values are compared with the next database.
    """

    # # Ignore Warning converting nan to masked element.
    # warnings.filterwarnings("ignore", category=UserWarning)

    # Initial Table with proper format and a single dummy entry.
    crossMdata = Table(
        [['NaN'], [-179.9] * u.deg, [-89.9] * u.deg, [0.]],
        names=('name', 'ra', 'dec', 'dist_pc'))

    procDatabases = []
    for DB_name, data in allDatabases.items():
        print("  (processing {})".format(DB_name))

        # idx2_unq, d2d_unq, idx2_ncm = nameMatch(crossMdata, data)

        # Find matches between this database, and the cross-matched data
        idx2_unq, d2d_unq, idx2_ncm = unqCrossMatch(crossMdata, data)

        n_m, ra_m, dec_m, dist_m = [], [], [], []
        idx1_ncm, idx1_m, idx2_m = [], [], []

        # For each element in crossMdata
        for i1, i2 in enumerate(idx2_unq):

            # No match found for this element.
            if np.isnan(i2):
                idx1_ncm.append(i1)
            else:
                # Found match
                # sim = similar(crossMdata['name'][i1], data['name'][i2])
                if d2d_unq[i1] < max_sep * u.arcsec:  # or sim:
                    # Store all names.
                    n_m.append(
                        # + ' ({})'.format(DB_name[0]))
                        crossMdata['name'][i1] + ', ' + data['name'][i2])

                    # Store averaged (ra, dec) values.
                    ra_m.append(np.mean([
                        crossMdata['ra'][i1], data['ra'][i2]]))
                    dec_m.append(np.mean([
                        crossMdata['dec'][i1], data['dec'][i2]]))
                    # Sum individual distances
                    dist_m.append(np.nansum([
                        crossMdata['dist_pc'][i1], data['dist_pc'][i2]]))

                    # Store indexes
                    idx1_m.append(i1)
                    idx2_m.append(i2)
                else:
                    # Rejected match, 'max_sep' criteria.
                    idx1_ncm.append(i1)
                    idx2_ncm.append(i2)

        n_nm, ra_nm, dec_nm, idx_nm, dist_nm = [], [], [], [], []
        for i1 in idx1_ncm:
            n_nm.append(crossMdata['name'][i1])
            ra_nm.append(crossMdata['ra'][i1])
            dec_nm.append(crossMdata['dec'][i1])
            idx_nm.append('--')
            dist_nm.append(crossMdata['dist_pc'][i1])

        for i2 in idx2_ncm:
            n_nm.append(data['name'][i2])
            ra_nm.append(data['ra'][i2])
            dec_nm.append(data['dec'][i2])
            idx_nm.append(i2)
            dist_nm.append(data['dist_pc'][i2])

        # Combine matched and not matched data.
        n_cmb = Column(n_m + n_nm, name='name')
        ra_cmb = Column(ra_m + ra_nm, name='ra', unit=u.degree)
        dec_cmb = Column(dec_m + dec_nm, name='dec', unit=u.degree)
        dist_cmb = Column(dist_m + dist_nm, name='dist_pc', unit=u.pc)
        tempData = Table([n_cmb, ra_cmb, dec_cmb, dist_cmb])

        # Add indexes pointing to the clusters in each database
        for n in allDatabases.keys():
            if n == DB_name:
                tempData.add_column(Column(
                    idx2_m + idx_nm, name=DB_name))
            elif n in procDatabases:
                tempData.add_column(Column(
                    list(crossMdata[n][idx1_m])
                    + list(crossMdata[n][idx1_ncm])
                    + ['--'] * (len(tempData) - len(idx1_m) - len(idx1_ncm)),
                    name=n))
            else:
                tempData.add_column(Column(
                    ['--'] * len(tempData), name=n))

        # Signal that this database was already processed.
        procDatabases.append(DB_name)
        # Update final database.
        crossMdata = tempData

    # Count the number of matches found for each cluster.
    all_msk = []
    for db_name in allDatabases.keys():
        all_msk.append(crossMdata[db_name] != '--')
    all_msk = np.array(all_msk).T
    crossMdata['N_m'] = all_msk.sum(1)

    # Store mean distance values
    crossMdata['dist_pc'] = crossMdata['dist_pc'] / crossMdata['N_m']

    # Remove initial dummy 'NaN' entry.
    crossMdata.add_index('name')
    nan_idx = crossMdata.loc['NaN'].index
    crossMdata.remove_row(nan_idx)

    print("Databases cross-matched")

    for Nm in np.arange(crossMdata['N_m'].max(), 0, -1):
        if Nm == 1:
            txt = 'No cross-match found'
        else:
            txt = '{} matches found'.format(Nm)
        print("  {}: {}".format(txt, (crossMdata['N_m'] == Nm).sum()))

    return crossMdata


def unqCrossMatch(crossMdata, data):
    """
    """
    # Define catalogs to be matched.
    c1 = SkyCoord(ra=crossMdata['ra'], dec=crossMdata['dec'])
    c2 = SkyCoord(ra=data['ra'], dec=data['dec'])

    # 'idx2' are indices into 'c2' that are the closest objects to each of
    # the coordinates in 'c1'.
    idx2, d2d, _ = c1.match_to_catalog_sky(c2)

    idx2_unq, d2d_unq = [], []
    # Duplicate matches: keep the closest match and replace the other with
    # a 'nan' value.
    for j, i2 in enumerate(idx2):

        # If this match is already stored.
        if i2 in idx2_unq:
            # Find the index for this index.
            i = idx2_unq.index(i2)
            # If this is a better match than the one stored.
            if d2d[j] < d2d_unq[i]:
                # Store 'nan' values in the 'old' position.
                idx2_unq[i] = np.nan
                d2d_unq[i] = np.nan
                # Replace the old index (and distance) with this new one.
                idx2_unq.append(i2)
                d2d_unq.append(d2d[j])
            else:
                # The index is already stored, but this is not a better
                # match. Store 'nan' values.
                idx2_unq.append(np.nan)
                d2d_unq.append(np.nan)
        else:
            # Index not in unique list, store.
            idx2_unq.append(i2)
            d2d_unq.append(d2d[j])

    # Indexes in 'data' with no match found.
    idx2_ncm = [_ for _ in range(len(c2)) if _ not in idx2_unq]

    return idx2_unq, d2d_unq, idx2_ncm
