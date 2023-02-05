
import os
import numpy as np
from astropy import units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from modules.IO import readINI, readData, write2File
from modules import makePlot, crossMatch


def main():
    """
    """
    dbs_dict, max_sep, defFlag, plotFlag, dpi, mode = readINI()

    # Read databases.
    print(f"\nReading {len(dbs_dict)} databases:")
    allDatabases = readData(dbs_dict)
    # Cross-match all databases.
    if len(allDatabases.keys()) <= 1:
        raise ValueError("At least two databases must be present in 'input/")

    N_tot = sum([len(v) for k, v in allDatabases.items()])
    print(f"Total number of clusters: {N_tot}")

    # Generate cross-match for all the databases
    crossMdata = crossMatch.match(allDatabases, max_sep)
    # Use 'len()-1' to exclude dummy entry
    print(f"\nDatabases cross-matched, clusters found: {len(crossMdata) - 1}")

    # Define Galactocentric frame
    gc_frame = frameGalactocentric(defFlag)

    # Generate the output table
    outTable = createTable(allDatabases, crossMdata, gc_frame)

    for Nm in np.arange(outTable['Nm'].max(), 0, -1):
        if Nm == 1:
            txt = 'No cross-match found'
        else:
            txt = '{} matches found'.format(Nm)
        print("  {}: {}".format(txt, (outTable['Nm'] == Nm).sum()))

    print("\nUnique clusters in DBS")
    all_dbs = []
    for i, cl in enumerate(outTable['name']):
        if outTable['Nm'][i] == 1:
            all_dbs.append(outTable['DBs'][i])
    all_dbs = np.array(all_dbs)
    DB_names = list(allDatabases.keys())
    for db in DB_names:
        msk = db == all_dbs
        print("  {}: {}".format(db, msk.sum()))

    # Write output file.
    write2File(outTable)

    if plotFlag:
        print("Plot mode: {}".format(mode))
        makePlot.plot(dpi, mode, outTable, gc_frame)


def createTable(allDatabases, crossMdata, gc_frame):
    """
    Create table to write to file
    """

    all_names, all_DBs, all_ra, all_de, all_dist, all_Nm, all_ra_avrg,\
        all_de_avrg = [[] for _ in range(8)]
    # For each matched cluster
    for i, ids in enumerate(crossMdata['IDs']):
        ids = ids.split(';')

        # For each DB
        r_names, r_DBs, r_ra, r_de, r_dist, r_Nm, r_ra_avrg, r_de_avrg =\
            [[] for _ in range(8)]
        Nm = 0
        for _id in ids:
            # Skip dummy entry
            if _id == 'dummy':
                continue
            idx, DB = _id.strip().split(' ')
            # Add this DB to the row
            if idx != '--':
                clust = allDatabases[DB][int(idx)]
                r_names.append(str(clust['name']) + ';')
                r_DBs.append(DB + ';')
                ra, de = round(clust['ra'], 4), round(clust['dec'], 4)
                r_ra.append(str(ra) + ';')
                r_de.append(str(de) + ';')
                dist = clust['dist_pc']
                if not np.isnan(dist):
                    dist = int(dist)
                r_dist.append(str(dist) + ';')
                Nm += 1

        if _id != 'dummy':
            r_ra_avrg = crossMdata['ra'][i]
            r_de_avrg = crossMdata['dec'][i]

        if not r_names:
            continue
        # Add row (remove trailing ';')
        all_names.append(''.join(r_names)[:-1])
        all_DBs.append(''.join(r_DBs)[:-1])
        all_ra.append(''.join(r_ra)[:-1])
        all_de.append(''.join(r_de)[:-1])
        all_dist.append(''.join(r_dist)[:-1])
        all_Nm.append(Nm)
        all_ra_avrg.append(r_ra_avrg)
        all_de_avrg.append(r_de_avrg)

    name = Column(all_names, name='name')
    DBs = Column(all_DBs, name='DBs')
    ra_all = Column(all_ra, name='ra_all')
    dec_all = Column(all_de, name='dec_all')
    dist = Column(all_dist, name='dist_all')
    Nm = Column(all_Nm, name='Nm')
    ra_avrg = Column(all_ra_avrg, name='ra', unit=u.degree)
    de_avrg = Column(all_de_avrg, name='dec', unit=u.degree)
    outTable = Table([name, DBs, ra_all, dec_all, dist, Nm, ra_avrg, de_avrg])

    # Add Cartesian data.
    outTable = dist2plane(outTable, gc_frame)

    return outTable


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


def dist2plane(outTable, gc_frame):
    """
   Obtain the Cartesian coordinates
    """
    dist_pc = []
    for i, cl_d in enumerate(outTable['dist_all']):
        ds = cl_d.split(';')
        if len(ds) == 1:
            ds = float(ds[0])
        else:
            msk = np.array(ds) == 'nan'
            if msk.sum() == len(ds):
                ds = np.nan
            else:
                ds = np.nanmean([float(_) for _ in ds])
        dist_pc.append(ds)
    outTable['dist_pc'] = Column(dist_pc, name='dist_pc', unit=u.pc)

    # Galactic coordinates.
    eq = SkyCoord(ra=outTable['ra'], dec=outTable['dec'], frame='icrs')
    lb = eq.transform_to('galactic')
    outTable['lon'] = lb.l.wrap_at(180 * u.deg).radian * u.radian
    outTable['lat'] = lb.b.radian * u.radian
    coords = SkyCoord(
        l=outTable['lon'], b=outTable['lat'],
        distance=outTable['dist_pc'], frame='galactic')

    # Galactocentric coordinates.
    c_glct = coords.transform_to(gc_frame)
    outTable['x_pc'], outTable['y_pc'], outTable['z_pc'] =\
        c_glct.x, c_glct.y, c_glct.z

    return outTable


if __name__ == '__main__':
    # Create /output dir if it does not exist.
    if not os.path.exists('output/'):
        os.makedirs('output/')
    main()
