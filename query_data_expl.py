
from astropy.io import ascii
import numpy as np
from scipy.spatial import distance
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.offsetbox as offsetbox
import matplotlib.gridspec as gridspec


def main(babusiaux_filters=False):
    """
    Explore data downloaded via the 'query.py' script.

    Vizier Gaia DR2 column names:
    http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2

    About the negative parallaxes in Gaia DR2 data:
    https://astronomy.stackexchange.com/q/26250/354
    https://astronomy.stackexchange.com/q/26071/354

    Vizier Pan-STARRS1 column names:
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/349&-to=3

    """

    # Define cluster name, center & radius
    # name = 'NGC6791'
    # center, d_max = (290.2208333, 37.7716667), .15
    # name = 'RUP44'
    # center, d_max = (119.7125, -28.5833), .15

    name = 'GAIA1'
    center, d_max = (101.47, -16.75), .15
    # name = 'GAIA1_ps1'
    # center, d_max = (101.47, -16.75), .15
    # name = 'GAIA1_gaia2_ps1'
    # center, d_max = (101.47, -16.75), .15

    print(name)
    data = ascii.read('input_expl/' + name + '.dat')
    N_old = len(data)
    print("Data read, {} sources".format(N_old))
    # for i, _ in enumerate(data.columns):
    #     print(i, _, data[_].dtype)

    if babusiaux_filters:
        data = babusiaux_filt(data)
        print("Filters applied, {:.1f}% of data lost".format(
            100. - (len(data) * 100.) / N_old))

    ra, dec, mag, e_mag, mag_n = data['RA_ICRS'], data['DE_ICRS'],\
        data['Gmag'], data['e_Gmag'], 'G'

    # # Two colors (Pan-STARRS1)
    # col1, e_col1, col2, e_col2, col1_n, col2_n = data['rmag'] - data['zmag'],\
    #     data['e_rmag'] + data['e_zmag'], data['imag'] - data['zmag'],\
    #     data['e_imag'] + data['e_zmag'], 'r-z', 'i-z'
    # Two colors (Gaia DR2)
    col1, e_col1, col2, e_col2, col1_n, col2_n = data['BP-G'],\
        data['e_BPmag'] + data['e_Gmag'], data['G-RP'],\
        data['e_RPmag'] + data['e_Gmag'], 'BP-G', 'G-RP'

    # Parallax, PMs, rad velocity.
    plx, pm_ra, e_pm_ra, pm_dec, e_pm_dec, radv = data['Plx'], data['pmRA'],\
        data['e_pmRA'], data['pmDE'], data['e_pmDE'], data['RV']

    # Max photometric errors.
    e_mmax, e_c1max, e_c2max = 0.05, 0.1, 0.1
    # Plx range [mas]
    plx_min, plx_max = 0., 5.

    makePlot(
        name, center, d_max, ra, dec, mag, e_mag, col1, e_col1, col2, e_col2,
        plx, pm_ra, e_pm_ra, pm_dec, e_pm_dec, e_mmax, e_c1max, e_c2max,
        plx_min, plx_max, radv, mag_n, col1_n, col2_n)


def babusiaux_filt(data):
    """
    Babusiaux et al. (2018) HRD filters.
    """
    m1 = (data['RPlx'] > 0.)  # 10.
    m2 = (data['RFG'] > 50.)
    m3 = (data['RFBP'] > 20.)  # 20.
    m4 = (data['RFRP'] > 20.)  # 20.
    m5 = (data['E_BR_RP_'] > 1. + 0.015 * (data['BPmag'] - data['RPmag']) ** 2)
    # m6 = (data['Gmag'] < 1.e6)
    m6 = (data['E_BR_RP_'] < 1.3 + 0.06 * (data['BPmag'] - data['RPmag']) ** 2)
    m7 = (data['Nper'] > 8)
    m8 = (data['chi2AL'] / (data['NgAL'] - 5.) < 1.44 * np.clip(
          np.exp(-.4 * (data['Gmag'] - 19.5)), a_min=None, a_max=1.))
    mask = m1 & m2 & m3 & m4 & m5 & m6 & m7 & m8
    for i, m in enumerate([m1, m2, m3, m4, m5, m6, m7, m8]):
        print("  m" + str(i + 1) + " removes {} sources".format(
            len(data) - m.data.sum()))

    return data[mask]


def makePlot(name, center, d_max, ra, dec, mag, e_mag, col1, e_col1, col2,
             e_col2, plx, pmRA, e_pmRA, pmDE, e_pmDE, e_mmax, e_c1max,
             e_c2max, plx_min, plx_max, radv, mag_n, col1_n, col2_n):
    """
    """
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(30, 25))
    gs = gridspec.GridSpec(10, 12)

    col_sizes = [[], []]
    col_names = ['RA'] + [mag_n, col1_n, col2_n] + ['Plx', 'pmRA', 'RV']
    for i, col in enumerate([ra, mag, col1, col2, plx, pmRA, radv]):
        # Count valid data for each column
        col_sizes[0].append(col_names[i])
        col_sizes[1].append(col[~col.mask].size)

    ax1 = plt.subplot(gs[0:2, 0:2])
    ax1.bar(col_sizes[0], col_sizes[1])
    fig.autofmt_xdate()

    ax2 = plt.subplot(gs[0:2, 2:4])
    ax2.set_title("N_T={}, d_max={:.2f} [deg]".format(
        ra.size, d_max), fontsize=8)
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    ax2.scatter(ra, dec, s=star_size(mag), c='k')
    # Radius
    circle = plt.Circle(center, d_max, color='red', lw=1.5, fill=False)
    fig.gca().add_artist(circle)
    ax2.invert_xaxis()

    ax4 = plt.subplot(gs[0:2, 4:6])
    ax4.scatter(mag, e_col1, label='e' + col1_n, s=5, lw=0., alpha=0.5)
    ax4.scatter(mag, e_col2, label='e' + col2_n, s=5, lw=0., alpha=0.5)
    ax4.scatter(mag, e_mag, label='e' + mag_n, s=5, lw=0., alpha=0.5)
    ax4.axhline(e_mmax, ls='--', c='g')
    ax4.axhline(e_c1max, ls='--', c='r')
    ax4.axhline(e_c2max, ls='--', c='r')
    plt.xlabel(mag_n)
    plt.ylim(-0.01, .5)
    plt.legend()

    # distance to cluster's center
    d_col = distance.cdist([center], np.array([ra, dec]).T)[0]

    # Mask for photometric diagrams
    msk = (d_col < d_max) & (e_mag < e_mmax) & (e_col1 < e_c1max) &\
        (e_col2 < e_c2max)

    ax5 = plt.subplot(gs[2:4, 0:2])
    ax5.set_title(
        "N(r<{}, e_mag<{}, e_c1<{}, e_c2<{})={}".format(
            d_max, e_mmax, e_c1max, e_c2max, msk.data.sum()), fontsize=8)
    plt.xlabel(col1_n)
    plt.ylabel(mag_n)
    ax5.scatter(col1[msk], mag[msk], s=4, lw=.1, edgecolor='w')
    # no masked elements
    msk4 = (~col1[msk].mask) & (~mag[msk].mask)
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'mag', col1[msk][msk4], mag[msk][msk4])
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)

    ax6 = plt.subplot(gs[2:4, 2:4])
    plt.xlabel(col2_n)
    plt.ylabel(mag_n)
    ax6.scatter(col2[msk], mag[msk], s=4, lw=.1, edgecolor='w')
    # no masked elements
    msk4 = (~col2[msk].mask) & (~mag[msk].mask)
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'mag', col2[msk][msk4], mag[msk][msk4])
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)

    ax9 = plt.subplot(gs[2:4, 4:6])
    plt.xlabel(col1_n)
    plt.ylabel(col2_n)
    ax9.scatter(col1[msk], col2[msk], s=4, lw=.1, edgecolor='w')
    # no masked elements
    msk4 = (~col1[msk].mask) & (~col2[msk].mask)
    x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = diag_limits(
        'col', col1[msk][msk4], col2[msk][msk4])
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)

    ax8 = plt.subplot(gs[4:6, 0:2])
    ax8.set_title("{} < Plx [mas] < {}".format(plx_min, plx_max), fontsize=8)
    plt.xlabel("Plx [mas]")
    msk2 = (d_col < d_max) & (e_mag < e_mmax) &\
        (e_col1 < e_c1max) & (e_col2 < e_c2max) &\
        (plx > plx_min) & (plx < plx_max)
    y, x, _ = ax8.hist(plx[msk2], bins=75)
    p_max_mas = (.5 * (x[y.argmax()] + x[y.argmax() + 1]))
    d_max_pc = 1000. / p_max_mas
    ax8.axvline(p_max_mas, ls='--', c='r')
    plx_lt_zero = 100. * plx[plx < 0.].size / plx.size
    ob = offsetbox.AnchoredText(
        "Plx_max={:.0f} [pc] ({:.3f} [mas])\nPlx<0: {:.1f}%".format(
            d_max_pc, p_max_mas, plx_lt_zero),
        pad=0.2, loc=1, prop=dict(size=9))
    ob.patch.set(alpha=0.85)
    ax8.add_artist(ob)
    plt.xlim(0., 3.)

    ax7 = plt.subplot(gs[4:6, 2:4])
    ax7.set_title("N(r<, e_G<, e_XP<, <Plx<)={}".format(
        msk2.data.sum()), fontsize=8)
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    cmap = cm.viridis_r
    norm = Normalize(vmin=0., vmax=p_max_mas)
    ax7.scatter(ra[msk2], dec[msk2], s=4, c=cmap(norm(plx[msk2])))
    ax7.invert_xaxis()
    im = plt.scatter(ra[msk2], dec[msk2], s=0, c=plx[msk2], cmap=cmap)
    cbar_ax = fig.add_axes([0.313, 0.53, 0.005, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=5)
    plt.clim(0., 2. * p_max_mas)

    ax3 = plt.subplot(gs[4:6, 4:6])
    plt.xlabel("pmRA [mas/yr]")
    plt.ylabel("pmDEC [mas/yr]")
    msk3 = (pmRA > -30) & (pmDE > -30) & (pmRA < 30) &\
        (pmDE < 30) & (d_col < d_max)
    pmRA_f, pmDE_f, epmRA_f, epmDE_f = pmRA[msk3], pmDE[msk3],\
        e_pmRA[msk3], e_pmDE[msk3]
    ax3.set_title("N(r<d_max, |pmX|<30)={}".format(pmRA_f.size), fontsize=8)
    cmap = cm.viridis
    norm = Normalize(vmin=d_col[msk3].min(), vmax=d_max)
    ax3.errorbar(
        pmRA_f, pmDE_f, yerr=epmDE_f, xerr=epmRA_f, fmt='none', elinewidth=.35,
        ecolor=cmap(norm(d_col[msk3])))
    im = plt.scatter(pmRA_f, pmDE_f, s=0, c=d_col[msk3], cmap=cmap)
    cbar_ax = fig.add_axes([0.48, 0.53, 0.005, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=5)

    fig.tight_layout()
    plt.savefig('output/gaia_' + name + '.png', dpi=150, bbox_inches='tight')


def star_size(mag, N=None, min_m=None):
    '''
    Convert magnitudes into intensities and define sizes of stars in
    finding chart.
    '''
    # Scale factor.
    if N is None:
        N = len(mag)
    if min_m is None:
        min_m = min(mag)
    factor = 500. * (1 - 1 / (1 + 150 / N ** 0.85))
    return 0.1 + factor * 10 ** ((np.array(mag) - min_m) / -2.5)


def kde_limits(phot_x, phot_y):
    '''
    Return photometric diagram limits taken from a 2D KDE.
    '''

    xmin, xmax = min(phot_x), max(phot_x)
    ymin, ymax = min(phot_y), max(phot_y)
    # Stack photometric data.
    values = np.vstack([phot_x, phot_y])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)
    # Grid density (number of points).
    gd = 25
    gd_c = complex(0, gd)
    # Define x,y grid.
    x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x.ravel(), y.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)

    # Generate 30 contour lines.
    plt.figure()
    cs = plt.contour(x, y, np.reshape(k_pos, x.shape), 30)
    plt.close()
    # Extract (x,y) points delimiting each line.
    x_v, y_v = np.asarray([]), np.asarray([])
    # Only use the outer curve.
    col = cs.collections[0]
    # If more than one region is defined by this curve (ie: the main sequence
    # region plus a RC region or some other detached region), obtain x,y from
    # all of them.
    for lin in col.get_paths():
        x_v = np.append(x_v, lin.vertices[:, 0])
        y_v = np.append(y_v, lin.vertices[:, 1])

    return x_v, y_v


def diag_limits(yaxis, phot_x, phot_y):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    x_v, y_v = kde_limits(phot_x, phot_y)

    # Define diagram limits.
    x_min_cmd, x_max_cmd = min(x_v) - 1.25, max(x_v) + 1.25
    y_min_cmd = max(y_v) + 1.25
    # If photometric axis y is a magnitude, make sure the brightest star
    # is always plotted.
    if yaxis == 'mag':
        y_max_cmd = min(phot_y) - 1.
    else:
        y_max_cmd = min(y_v) - 1.

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


if __name__ == '__main__':
    main()
