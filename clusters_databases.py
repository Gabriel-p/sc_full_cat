
import warnings
from astropy.io import ascii
from astropy import units as u
# from astropy.coordinates.distances import Distance
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import numpy as np


def main(max_sep=600., plotFlag=True):
    """
    Example plots:

    * http://www.astropy.org/astropy-tutorials/plot-catalog.html
    * http://docs.astropy.org/en/stable/coordinates/skycoord.html
    * http://docs.astropy.org/en/stable/visualization/wcsaxes/
      overlaying_coordinate_systems.html
    * http://docs.astropy.org/en/stable/visualization/wcsaxes/
    * http://adass.org/adass/proceedings/adass94/greisene.html

    """

    # Read databases.
    print("Read all databases.")
    allData = readData()

    # Cross-match all databases.
    print("Perform cross-match (max_sep={}).".format(max_sep))
    crossMdata = crossMatch(allData, max_sep)
    write2File(allData, crossMdata)

    if plotFlag:
        print('Plot each database.')
        # Plot each catalog separately.
        for name, data in allData.items():
            makePlot(name, data)

        print('Plot Cross-matched data.')
        makePlot('Crossmatch', crossMdata)

        # Cross-match with no MWSC database.
        del allData['MWSC']
        crossMdata = crossMatch(allData, max_sep)
        print('Plot CrossM (no MWSC) data.')
        makePlot('CrossM (no MWSC)', crossMdata)


def readData():
    """
    Read all databases.
    """

    # OPENCLUST - New Optically Visible Open Clusters and Candidates Catalog
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
    openclst['col6'].name = 'dist_pc'
    openclst['dist_pc'] = openclst['dist_pc'].astype(float)
    openclst['col8'].name = 'log_age'
    openclst['col18'].name = 'fe_h'

    openclst = dist2plane(openclst)

    # MWSC - Milky Way Star Clusters Catalog
    # https://heasarc.gsfc.nasa.gov/W3Browse/all/mwsc.html
    #
    # Manual edits: removed initial '|' chars.
    mwsc = ascii.read('input/MWSC.dat', format='fixed_width')
    mwsc['distance'].name = 'dist_pc'
    mwsc['dist_pc'] = mwsc['dist_pc'].astype(float)
    mwsc['metallicity'].name = 'fe_h'

    mwsc = dist2plane(mwsc)
    # m10 = abs(mwsc['z_pc']) > 10000
    # print(mwsc[m10])

    # Camargo et al (2010-2013); Table 1, 2
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J%2FMNRAS%2F432%2F3349
    #
    # Manual edits: stitched together both tables.
    camargo = ascii.read('input/Camargo.dat')
    camargo['Cluster'].name = 'name'
    camargo['RAJ2000'].name = 'ra'
    camargo['DEJ2000'].name = 'dec'
    camargo['Dist'].name = 'dist_pc'
    camargo['dist_pc'] = camargo['dist_pc'] * 1000.
    camargo['Age'].name = 'log_age'
    camargo['log_age'] = np.log10(camargo['log_age'] * 1000000.)
    camargo['fe_h'] = np.array([np.nan for _ in camargo])

    camargo = dist2plane(camargo)

    # WEBDA - http://www.univie.ac.at/webda/
    #
    # Manual edits: downloaded from the parameters form, removed a line
    # 'Back to WEBDA home page'
    webda = ascii.read('input/WEBDA.dat')
    webda['Cluster_name'].name = 'name'
    webda['RA_2000'].name = 'ra'
    webda['Dec_2000'].name = 'dec'
    webda['Age'].name = 'log_age'
    webda['Dist'].name = 'dist_pc'
    webda['dist_pc'] = webda['dist_pc'].astype(float)
    webda['Fe/H'].name = 'fe_h'

    webda = dist2plane(webda)

    return {'Camargo': camargo, 'OPENCLUST': openclst, 'MWSC': mwsc,
            'WEBDA': webda}


def dist2plane(data, eq2rad=True):
    """
    Convert equatorial coordinates (if flag is True), and obtain the
    Cartesian coordinates with 'z_pc' the vertical distance.
    """
    if eq2rad:
        eq = SkyCoord(
            ra=data['ra'], dec=data['dec'], unit=(u.hour, u.deg), frame='icrs')
        data['ra'] = eq.ra.radian
        data['dec'] = eq.dec.radian
    else:
        eq = SkyCoord(
            ra=data['ra'] * u.radian, dec=data['dec'] * u.radian, frame='icrs')

    lb = eq.transform_to('galactic')
    data['lon'] = lb.l.wrap_at(180 * u.deg).radian
    data['lat'] = lb.b.radian

    try:
        data['dist_pc'] = data['dist_pc'].filled(np.nan)
    except AttributeError:
        pass

    # Cartesian coordinates.
    coords = SkyCoord(
        l=data['lon'] * u.radian, b=data['lat'] * u.radian,
        distance=data['dist_pc'] * u.pc, frame='galactic')
    data['x_pc'], data['y_pc'], data['z_pc'] = coords.cartesian.x,\
        coords.cartesian.y, coords.cartesian.z

    return data


def crossMatch(allData, max_sep=5.):
    """
    Cross-match all the databases in 'allData', using the 'max_sep' value
    as the limiting match radius in arcsec.
    """

    # Initial Table with proper format and a single dummy entry.
    crossMdata = Table(
        [['NaN'], [-np.pi / 2.], [-np.pi / 2.], [np.nan]],
        names=('name', 'ra', 'dec', 'dist_pc'))

    procDatabases = []
    for data_name, data in allData.items():

        idx2_unq, d2d_unq, idx2_ncm = unqCrossMatch(crossMdata, data)

        n_m, ra_m, dec_m, s_dist_m = [], [], [], []
        idx1_ncm, idx1_m, idx2_m = [], [], []
        # For each element in crossMdata
        for i1, i2 in enumerate(idx2_unq):

            # No match found for this 'crossMdata' element.
            if np.isnan(i2):
                idx1_ncm.append(i1)
            else:
                # Found match
                if d2d_unq[i1] < max_sep * u.arcsec:
                    # Store all names.
                    n_m.append(
                        crossMdata['name'][i1] + ', ' + data['name'][i2])
                    # Store averaged (ra, dec) values.
                    ra_m.append(np.mean([
                        crossMdata['ra'][i1], data['ra'][i2]]))
                    dec_m.append(np.mean([
                        crossMdata['dec'][i1], data['dec'][i2]]))
                    # Sum distances.
                    s_dist_m.append(
                        nansumwrapper([
                            crossMdata['dist_pc'][i1], data['dist_pc'][i2]]))
                    idx1_m.append(i1)
                    idx2_m.append(i2)
                else:
                    # Rejected match, 'max_sep' criteria.
                    idx1_ncm.append(i1)
                    idx2_ncm.append(i2)

        n_nm, ra_nm, dec_nm, dist_nm, idx_nm = [], [], [], [], []
        for i1 in idx1_ncm:
            n_nm.append(crossMdata['name'][i1])
            ra_nm.append(crossMdata['ra'][i1])
            dec_nm.append(crossMdata['dec'][i1])
            dist_nm.append(crossMdata['dist_pc'][i1])
            idx_nm.append('--')

        for i2 in idx2_ncm:
            n_nm.append(data['name'][i2])
            ra_nm.append(data['ra'][i2])
            dec_nm.append(data['dec'][i2])
            dist_nm.append(data['dist_pc'][i2])
            idx_nm.append(i2)

        # Combine matched and not matched data.
        n_cmb = Column(n_m + n_nm, name='name')
        ra_cmb = Column(ra_m + ra_nm, name='ra')
        dec_cmb = Column(dec_m + dec_nm, name='dec')
        d_cmb = Column(s_dist_m + dist_nm, name='dist_pc')
        tempData = Table([n_cmb, ra_cmb, dec_cmb, d_cmb])

        for n in allData.keys():
            if n == data_name:
                tempData.add_column(Column(
                    idx2_m + idx_nm, name=data_name))
            elif n in procDatabases:
                tempData.add_column(Column(
                    list(crossMdata[n][idx1_m]) +
                    list(crossMdata[n][idx1_ncm]) +
                    ['--'] * (len(tempData) - len(idx1_m) - len(idx1_ncm)),
                    name=n))
            else:
                tempData.add_column(Column(
                    ['--'] * len(tempData), name=n))

        # Signal that this database was already processed.
        procDatabases.append(data_name)
        # Update final database.
        crossMdata = tempData

    # Count the number of valid distance values stored for each cross-matched
    # cluster.
    N_d = Column([0.] * len(crossMdata))
    for n in allData.keys():
        m1 = crossMdata[n] != '--'
        m2 = crossMdata[n] == '--'
        N_d[m1] = N_d[m1] + 1.
        N_d[m2] = N_d[m2]
    crossMdata['N_d'] = N_d
    # Replace sum of distances for its mean.
    crossMdata['dist_pc'] = crossMdata['dist_pc'] / N_d

    # Add Cartesian data.
    crossMdata = dist2plane(crossMdata, eq2rad=False)

    # Remove dummy 'NaN' entry.
    crossMdata.add_index('name')
    nan_idx = crossMdata.loc['NaN'].index
    crossMdata.remove_row(nan_idx)

    return crossMdata


def nansumwrapper(a):
    """
    Source: https://stackoverflow.com/a/48405633/1391441
    """
    if np.isnan(a).all():
        return np.nan
    else:
        return np.nansum(a)


def unqCrossMatch(crossMdata, data):
    """
    """
    # Define catalogs to be matched.
    c1 = SkyCoord(
        ra=crossMdata['ra'] * u.radian, dec=crossMdata['dec'] * u.radian)
    c2 = SkyCoord(ra=data['ra'] * u.radian, dec=data['dec'] * u.radian)

    # 'idx' are indices into 'c2' that are the closest objects to each of
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
            # If this is a better match then the one stored.
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


def write2File(allData, crossMdata):
    """
    Write cross-matched data to file.
    """
    # Equatorial to degrees (from radians)
    eq = SkyCoord(
        ra=crossMdata['ra'] * u.radian, dec=crossMdata['dec'] * u.radian)
    crossMdata['ra'] = eq.ra
    crossMdata['dec'] = eq.dec

    dtBs_names = list(allData.keys())
    col_order = ['name', 'ra', 'dec'] + dtBs_names +\
        ['N_d', 'dist_pc', 'lon', 'lat', 'x_pc', 'y_pc', 'z_pc']
    ascii.write(
        crossMdata[col_order], 'output/crossMdata.dat', format='fixed_width',
        formats={'ra': '%12.7f', 'dec': '%12.7f', 'N_d': '%5.0f',
                 'dist_pc': '%10.2f', 'lon': '%10.4f', 'lat': '%10.4f',
                 'x_pc': '%10.2f', 'y_pc': '%10.2f', 'z_pc': '%10.2f'},
        overwrite=True)
    print("Cross-matched data written to file.")


def makePlot(name, data):
    """
    Gridspec idea: http://www.sc.eso.org/~bdias/pycoffee/codes/20160407/
                   gridspec_demo.html
    """

    # Ignore colorbar warning
    warnings.filterwarnings(
        "ignore", category=UserWarning, module="matplotlib")
    # Ignore RuntimeWarning when creating the masks
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # Vertical distance masks.
    mnan = np.isnan(data['z_pc'])
    m600 = abs(data['z_pc']) <= 600
    m1000 = (600 < abs(data['z_pc'])) & (abs(data['z_pc']) <= 1000)
    m2500 = (1000 < abs(data['z_pc'])) & (abs(data['z_pc']) <= 2500)
    minf = abs(data['z_pc']) > 2500

    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(
        6, 6, height_ratios=[1, 1, 1, 0.04, .7, .7],
        width_ratios=[1, 1, 1, 1, 1, 1])
    gs.update(hspace=0.03, wspace=.3)

    ax = plt.subplot(gs[0:3, 0:6], projection="aitoff")
    plt.title('Database: {} (N={})'.format(name, len(data)), y=1.02)
    ax.set_xticklabels([
        r'210$^{\circ}$', r'240$^{\circ}$', r'270$^{\circ}$', r'300$^{\circ}$',
        r'330$^{\circ}$', r'0$^{\circ}$', r'30$^{\circ}$', r'60$^{\circ}$',
        r'90$^{\circ}$', r'120$^{\circ}$', r'150$^{\circ}$'])
    ax.grid(True)

    plt_data = {
        '0': [3., .3, 'red', r'$No\,dist\; (N={})$'],
        '1': [4., .3, 'grey', r'$z\leq 600\,[pc]\; (N={})$'],
        '2': [4., .3, 'orange', r'$600<z\leq1000\,[pc]\;(N={})$'],
        '3': [6., .4, 'blue', r'$1000<z\leq2500\,[pc]\;(N={})$'],
        '4': [8., .5, 'green', r'$z>2500\,[pc]\;(N={})$']
    }
    for i, m in enumerate([mnan, m600, m1000, m2500, minf]):
        # Only plot if there are objects to plot.
        if sum(m) > 0:
            ms, a, c, lab = plt_data[str(i)]
            ax.plot(
                data['lon'][m], data['lat'][m], 'o', markersize=ms, alpha=a,
                color=c, label=lab.format(len(data['lon'][m])))
            # Only plot cluster names for those with the largest 'z' values.
            if i == 4:
                for _, (lon, lat) in enumerate(
                        zip(*[data['lon'][minf], data['lat'][minf]])):
                    ax.annotate(data['name'][minf][_], (lon, lat),
                                xycoords='data')
                    # print(data['name'][minf][_], data['lon'][minf][_],
                    #       data['lat'][minf][_], data['z_pc'][minf][_])
    plt.legend()

    plt.style.use('seaborn-darkgrid')

    data['x_pc'].convert_unit_to('kpc')
    data['y_pc'].convert_unit_to('kpc')
    data['z_pc'].convert_unit_to('kpc')

    # Sun's distance to galactic center from Camargo et al (2013) (taken from
    # Bica et al. 2006)
    R_0 = 7.2 * u.kpc
    # Move Galactic center to the origin of the system.
    x_gc = data['x_pc'] - R_0
    # Sun's coords
    s_xys = SkyCoord(-R_0, 0., 0., unit='kpc', representation_type='cartesian')

    ax = plt.subplot(gs[4:6, 0:2])
    plt.xlabel(r"$x_{GC}\, [kpc]$")
    plt.ylabel(r"$y_{GC}\, [kpc]$")
    vmin, vmax = max(min(data['z_pc']), -2.5), min(max(data['z_pc']), 2.5)
    plt1 = plt.scatter(
        x_gc, data['y_pc'].data, alpha=.5, c=data['z_pc'], cmap='viridis',
        vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.x, s_xys.y, c='yellow', s=50, edgecolor='k')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # Plot spiral arms
    spiral_arms = spiralArms()
    for sp_name, vals in spiral_arms.items():
        if checkInRange(vals, xmin, xmax, ymin, ymax):
            xy_arm = np.array(list(zip(*vals)))
            if sp_name == 'Orion-Cygnus':
                plt.plot(xy_arm[0], xy_arm[1], c='k', label=sp_name)
            elif sp_name == 'Carina-Sagittarius':
                plt.plot(xy_arm[0], xy_arm[1], c='b', ls='--', label=sp_name)
            elif sp_name == 'Crux-Scutum':
                plt.plot(
                    xy_arm[0], xy_arm[1], c='purple', ls='--', label=sp_name)
            else:
                plt.plot(xy_arm[0], xy_arm[1], ls='--', label=sp_name)
    plt.xlim(max(xmin, -20.), min(xmax, 20.))
    plt.ylim(max(ymin, -15.), min(ymax, 15.))
    plt.legend()
    # colorbar
    cbax = plt.subplot(gs[3:4, 0:2])
    cb = Colorbar(
        ax=cbax, mappable=plt1, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < z_{{GC}}\, [kpc] < {:.1f}$".format(
        min(data['z_pc']), max(data['z_pc'])), labelpad=10)

    ax = plt.subplot(gs[4:6, 2:4])
    plt.xlabel(r"$x_{GC}\, [kpc]$")
    plt.ylabel(r"$z_{GC}\, [kpc]$")
    vmin, vmax = max(ymin, -15.), min(ymax, 15.)
    plt2 = plt.scatter(
        x_gc, data['z_pc'], alpha=.5, c=data['y_pc'], cmap='viridis',
        vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.x, s_xys.z, c='yellow', s=50, edgecolor='k')
    plt.xlim(max(xmin, -20.), min(xmax, 20.))
    plt.ylim(max(min(data['z_pc']), -3.), min(max(data['z_pc']), 3.))
    # colorbar
    cbax = plt.subplot(gs[3:4, 2:4])
    cb = Colorbar(
        ax=cbax, mappable=plt2, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < y_{{GC}}\, [kpc] < {:.1f}$".format(
        min(data['y_pc']), max(data['y_pc'])), labelpad=10)

    ax = plt.subplot(gs[4:6, 4:6])
    plt.xlabel(r"$y_{GC}\, [kpc]$")
    plt.ylabel(r"$z_{GC}\, [kpc]$")
    vmin, vmax = max(xmin, -20.), min(xmax, 20.)
    plt3 = plt.scatter(
        data['y_pc'], data['z_pc'], alpha=.5, c=x_gc, cmap='viridis',
        vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.y, s_xys.z, c='yellow', s=50, edgecolor='k')
    plt.xlim(max(ymin, -15.), min(ymax, 15.))
    plt.ylim(max(min(data['z_pc']), -3.), min(max(data['z_pc']), 3.))
    # colorbar
    cbax = plt.subplot(gs[3:4, 4:6])
    cb = Colorbar(
        ax=cbax, mappable=plt3, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < x_{{GC}}\, [kpc] < {:.1f}$".format(
        min(x_gc.data), max(x_gc.data)), labelpad=10)

    fig.tight_layout()
    fig.savefig('output/' + name + '.png', dpi=150, bbox_inches='tight')
    plt.style.use('default')


def checkInRange(xy_data, xmin, xmax, ymin, ymax):
    """
    Check if at least one point is between the given ranges.
    """
    for x, y in xy_data:
        if xmin <= x <= xmax and ymin <= y <= ymax:
            return True
    return False


def spiralArms():
    """
    Obtained from Momany et al. (2006) "Outer structure of the..."
    """
    spiral_arms = {
        'Outer': (
            (7.559999999999999, 2.514124293785308),
            (7.237333333333332, 3.4180790960452008),
            (6.738666666666667, 4.46327683615819),
            (6.1519999999999975, 5.480225988700564),
            (5.535999999999998, 6.384180790960453),
            (4.861333333333334, 7.1186440677966125),
            (4.479999999999997, 7.514124293785311),
            (3.8053333333333335, 7.909604519774014),
            (2.8373333333333335, 8.474576271186443),
            (1.751999999999999, 9.011299435028253),
            (0.6373333333333306, 9.35028248587571),
            (-0.4480000000000022, 9.63276836158192),
            (-1.2693333333333356, 9.830508474576273),
            (-1.8560000000000016, 9.943502824858761),
            (-2.677333333333337, 10),
            (-3.4986666666666686, 9.943502824858761),
            (-4.085333333333336, 9.915254237288138),
            (-5.024000000000003, 9.830508474576273),
            (-5.845333333333336, 9.689265536723166),
            (-6.4906666666666695, 9.519774011299436),
            (-7.253333333333337, 9.322033898305087),
            (-7.986666666666669, 8.926553672316388),
            (-8.485333333333337, 8.559322033898308),
            (-9.160000000000004, 8.050847457627121),
            (-9.864000000000004, 7.570621468926557),
            (-10.333333333333336, 7.203389830508474),
            (-10.89066666666667, 6.525423728813561),
            (-11.389333333333337, 5.9604519774011315),
            (-11.77066666666667, 5.395480225988699),
            (-12.240000000000004, 4.661016949152543),
            (-12.65066666666667, 3.9830508474576263),
            (-13.061333333333337, 3.361581920903955),
            (-13.442666666666671, 2.7401129943502838)),
        'Perseus': (
            (5.2719999999999985, 2.372881355932204),
            (4.9786666666666655, 3.27683615819209),
            (4.773333333333333, 3.6723163841807924),
            (4.3039999999999985, 4.46327683615819),
            (3.776, 5.141242937853107),
            (3.3359999999999985, 5.593220338983052),
            (2.690666666666665, 6.073446327683616),
            (2.074666666666662, 6.468926553672318),
            (1.4879999999999995, 6.751412429378529),
            (0.607999999999997, 7.005649717514125),
            (-0.15466666666666917, 7.146892655367235),
            (-0.8293333333333344, 7.259887005649716),
            (-1.621333333333336, 7.316384180790962),
            (-2.4426666666666694, 7.288135593220339),
            (-3.14666666666667, 7.1186440677966125),
            (-3.909333333333336, 7.005649717514125),
            (-4.6426666666666705, 6.666666666666668),
            (-5.3173333333333375, 6.3559322033898304),
            (-5.962666666666669, 5.903954802259886),
            (-6.57866666666667, 5.395480225988699),
            (-6.989333333333336, 5.028248587570623),
            (-7.63466666666667, 4.350282485875709),
            (-8.104000000000003, 3.7570621468926575),
            (-8.51466666666667, 3.10734463276836),
            (-8.896000000000004, 2.4011299435028235),
            (-9.218666666666671, 1.6384180790960414),
            (-9.424000000000003, 1.073446327683616),
            (-9.65866666666667, 0.1977401129943459),
            (-9.805333333333337, -0.5084745762711869),
            (-9.92266666666667, -0.8192090395480243),
            (-10.010666666666669, -1.2711864406779654),
            (-10.128000000000004, -2.005649717514128),
            (-10.186666666666671, -2.711864406779661),
            (-10.186666666666671, -2.909604519774014),
            (-10.128000000000004, -3.3615819209039586),
            (-9.952000000000004, -4.2372881355932215),
            (-9.83466666666667, -5),
            (-9.688000000000002, -5.310734463276839),
            (-9.48266666666667, -5.734463276836161),
            (-9.189333333333337, -6.29943502824859),
            (-8.86666666666667, -7.005649717514126),
            (-8.456000000000003, -7.598870056497178)),
        'Orion-Cygnus': (
            (-7.341333333333337, 3.2485875706214706),
            (-7.63466666666667, 2.909604519774014),
            (-7.9280000000000035, 2.485875706214692),
            (-8.280000000000003, 1.9209039548022595),
            (-8.456000000000003, 1.4124293785310726),
            (-8.60266666666667, 1.1016949152542352),
            (-8.808000000000003, 0.5649717514124291),
            (-9.013333333333335, -0.197740112994353),
            (-9.13066666666667, -0.7627118644067821),
            (-9.160000000000004, -1.2146892655367267)),
        'Carina-Sagittarius': (
            (2.8373333333333335, 3.6723163841807924),
            (2.5146666666666633, 4.152542372881356),
            (2.338666666666665, 4.350282485875709),
            (1.9280000000000008, 4.774011299435031),
            (1.2239999999999966, 5.16949152542373),
            (0.49066666666666237, 5.451977401129945),
            (-0.12533333333333552, 5.593220338983052),
            (-0.7706666666666688, 5.621468926553675),
            (-1.7680000000000025, 5.53672316384181),
            (-2.5600000000000023, 5.310734463276834),
            (-2.970666666666668, 5.112994350282488),
            (-3.645333333333337, 4.576271186440678),
            (-3.8800000000000026, 4.350282485875709),
            (-4.2613333333333365, 3.9830508474576263),
            (-4.613333333333337, 3.5593220338983045),
            (-4.789333333333337, 3.192090395480225),
            (-5.1413333333333355, 2.627118644067796),
            (-5.493333333333336, 2.033898305084744),
            (-5.845333333333336, 1.4124293785310726),
            (-6.22666666666667, 0.6497175141242906),
            (-6.608000000000002, -0.14124293785311082),
            (-6.842666666666668, -0.7627118644067821),
            (-7.048000000000004, -1.5536723163841835),
            (-7.19466666666667, -2.25988700564972),
            (-7.253333333333337, -3.1073446327683634),
            (-7.136000000000003, -3.6440677966101696),
            (-7.048000000000004, -3.98305084745763),
            (-6.813333333333338, -4.519774011299436),
            (-6.461333333333336, -5.254237288135597),
            (-6.05066666666667, -5.875706214689268),
            (-5.6106666666666705, -6.440677966101699),
            (-5.024000000000003, -6.977401129943505),
            (-4.466666666666669, -7.485875706214692),
            (-3.8800000000000026, -7.909604519774014),
            (-3.352000000000002, -8.163841807909607),
            (-3.0880000000000027, -8.361581920903957)),
        'Crux-Scutum': (
            (1.663999999999998, 3.1355932203389827),
            (1.3119999999999976, 3.4180790960452008),
            (0.6666666666666643, 3.8135593220338997),
            (-0.09600000000000186, 3.9830508474576263),
            (-0.858666666666668, 4.0395480225988685),
            (-1.5626666666666686, 3.926553672316384),
            (-2.2666666666666693, 3.5875706214689274),
            (-2.9413333333333362, 3.192090395480225),
            (-3.5280000000000022, 2.6553672316384187),
            (-3.909333333333336, 2.033898305084744),
            (-4.29066666666667, 1.3276836158192076),
            (-4.584000000000003, 0.6214689265536713),
            (-4.760000000000003, -0.11299435028248794),
            (-4.906666666666668, -0.9322033898305087),
            (-4.965333333333335, -1.6666666666666714),
            (-4.994666666666669, -2.344632768361585),
            (-4.8480000000000025, -3.2203389830508478),
            (-4.672000000000002, -3.8983050847457648),
            (-4.320000000000004, -4.576271186440678),
            (-3.8800000000000026, -5.225988700564974),
            (-3.4106666666666694, -5.734463276836161),
            (-2.765333333333336, -6.158192090395483),
            (-2.032000000000002, -6.610169491525427),
            (-1.3866666666666685, -6.920903954802263),
            (-0.6826666666666696, -7.06214689265537),
            (0.05066666666666464, -7.25988700564972),
            (0.8719999999999999, -7.401129943502829),
            (1.575999999999997, -7.372881355932208),
            (2.2799999999999976, -7.316384180790964)),
        'Norma': (
            (-3.14666666666667, 0.8474576271186436),
            (-3.3226666666666684, 0.1977401129943459),
            (-3.2933333333333366, -0.7627118644067821),
            (-3.2346666666666692, -1.4406779661016955),
            (-2.970666666666668, -2.1186440677966125),
            (-2.5600000000000023, -2.824858757062149),
            (-2.1200000000000028, -3.4463276836158236),
            (-1.6506666666666696, -3.9548022598870105),
            (-0.9760000000000026, -4.378531073446332),
            (-0.2720000000000038, -4.689265536723166),
            (0.4319999999999986, -4.858757062146896),
            (0.9600000000000009, -4.887005649717519),
            (1.370666666666665, -4.858757062146896),
            (2.0453333333333283, -4.717514124293789),
            (2.6613333333333316, -4.519774011299436),
            (3.042666666666662, -4.406779661016952),
            (3.4826666666666632, -4.152542372881356),
            (4.010666666666662, -3.7288135593220346),
            (4.538666666666664, -3.050847457627121))
    }

    return spiral_arms


if __name__ == '__main__':
    main()
