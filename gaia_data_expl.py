
from astropy.io import ascii
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.offsetbox as offsetbox
import matplotlib.gridspec as gridspec


def main():
    """
    Explore Gaia DR2 data downloaded via the 'query.py' script.

    Vizier Gaia DR2 column names:
    http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2

    About the negative parallax:
    https://astronomy.stackexchange.com/q/26250/354
    https://astronomy.stackexchange.com/q/26071/354
    """

    # Define cluster name, center & radius
    # name = 'NGC6791'
    # center, d_max = (290.2208333, 37.7716667), .15
    # name = 'GAIA1'
    # center, d_max = (101.47, -16.75), .15
    name = 'RUP44'
    center, d_max = (119.7125, -28.5833), .15

    print(name)
    data = ascii.read('input/' + name + '.dat')
    N_old = len(data)
    print("Data read, {} sources".format(N_old))

    data = babusiaux_filt(data)
    print("Filters applied, {:.1f}% of data lost".format(
        100. - (len(data) * 100.) / N_old))

    # Max photometric errors for mask
    e_Gmax, e_XPmax = 0.05, 0.1
    # Plx range [mas]
    plx_min, plx_max = 0., 5.

    makePlot(
        name, center, d_max, e_Gmax, e_XPmax, plx_min, plx_max, data)


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


def makePlot(name, center, d_max, e_Gmax, e_XPmax, plx_min, plx_max, data):
    """
    """
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure(figsize=(30, 25))
    gs = gridspec.GridSpec(10, 12)

    col_sizes = [[], []]
    for col in ['RA_ICRS', 'Gmag', 'BP-RP', 'Plx', 'pmRA', 'RV']:
        # Count valid data for each column
        col_sizes[0].append(col)
        col_sizes[1].append(data[col][~data[col].mask].size)

    ax1 = plt.subplot(gs[0:2, 0:2])
    ax1.bar(col_sizes[0], col_sizes[1])
    fig.autofmt_xdate()

    ax2 = plt.subplot(gs[0:2, 2:4])
    ax2.set_title("N_T={}, d_max={:.2f} [deg]".format(
        data['RA_ICRS'].size, d_max), fontsize=8)
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    ax2.scatter(
        data['RA_ICRS'], data['DE_ICRS'], s=star_size(data['Gmag']), c='k')
    # Radius
    circle = plt.Circle(center, d_max, color='red', lw=1.5, fill=False)
    fig.gca().add_artist(circle)
    ax2.invert_xaxis()

    ax4 = plt.subplot(gs[0:2, 4:6])
    ax4.scatter(
        data['Gmag'], data['e_BPmag'], label='eBP', s=5, lw=0., alpha=0.5)
    ax4.scatter(
        data['Gmag'], data['e_RPmag'], label='eRP', s=5, lw=0., alpha=0.5)
    ax4.scatter(
        data['Gmag'], data['e_Gmag'], label='eG', s=5, lw=0., alpha=0.5)
    ax4.axhline(e_Gmax, ls='--', c='g')
    ax4.axhline(e_XPmax, ls='--', c='r')
    plt.xlabel('G')
    plt.ylim(-0.01, .5)
    plt.legend()

    # distance to cluster's center
    d_col = distance.cdist(
        [center], np.array([data['RA_ICRS'], data['DE_ICRS']]).T)[0]

    # Mask for photometric diagrams
    msk = (d_col < d_max) & (data['e_Gmag'] < e_Gmax) &\
        (data['e_BPmag'] < e_XPmax) & (data['e_RPmag'] < e_XPmax)

    ax5 = plt.subplot(gs[2:4, 0:2])
    ax5.set_title(
        "N(r<{}, e_Gmag<{} , e_XPmag<{})={}".format(
            d_max, e_Gmax, e_XPmax, msk.data.sum()), fontsize=8)
    plt.xlabel("(BP-RP)")
    plt.ylabel("G")
    ax5.scatter(
        data['BP-RP'][msk], data['Gmag'][msk], s=4, lw=.1, edgecolor='w')
    # ax5.set_xlim(np.nanmin(data['BP-RP'][msk]), 4.)
    ax5.invert_yaxis()

    ax6 = plt.subplot(gs[2:4, 2:4])
    plt.xlabel("(BP-RP)")
    plt.ylabel("(G-RP)")
    ax6.scatter(
        data['BP-RP'][msk], data['G-RP'][msk], s=4, lw=.1, edgecolor='w')
    ax6.invert_yaxis()

    ax9 = plt.subplot(gs[2:4, 4:6])
    plt.xlabel("(BP-RP)")
    plt.ylabel("(BP-G)")
    ax9.scatter(
        data['BP-RP'][msk], data['BP-G'][msk], s=4, lw=.1, edgecolor='w')
    ax9.invert_yaxis()

    ax8 = plt.subplot(gs[4:6, 0:2])
    ax8.set_title("{} < Plx [mas] < {}".format(plx_min, plx_max), fontsize=8)
    plt.xlabel("Plx [mas]")
    msk2 = (d_col < d_max) & (data['e_Gmag'] < e_Gmax) &\
        (data['e_BPmag'] < e_XPmax) & (data['e_RPmag'] < e_XPmax) &\
        (data['Plx'] > plx_min) & (data['Plx'] < plx_max)
    y, x, _ = ax8.hist(data['Plx'][msk2], bins=75)
    p_max_mas = (.5 * (x[y.argmax()] + x[y.argmax() + 1]))
    d_max_pc = 1000. / p_max_mas
    ax8.axvline(p_max_mas, ls='--', c='r')
    ob = offsetbox.AnchoredText(
        "Plx_max={:.0f} [pc] ({:.3f} [mas])".format(d_max_pc, p_max_mas),
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
    ax7.scatter(
        data['RA_ICRS'][msk2], data['DE_ICRS'][msk2], s=4,
        c=cmap(norm(data['Plx'][msk2])))
    ax7.invert_xaxis()
    im = plt.scatter(
        data['RA_ICRS'][msk2], data['DE_ICRS'][msk2], s=0,
        c=data['Plx'][msk2], cmap=cmap)
    cbar_ax = fig.add_axes([0.313, 0.53, 0.005, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=5)
    plt.clim(0., 2. * p_max_mas)

    ax3 = plt.subplot(gs[4:6, 4:6])
    plt.xlabel("pmRA [mas/yr]")
    plt.ylabel("pmDEC [mas/yr]")
    msk3 = (data['pmRA'] > -30) & (data['pmDE'] > -30) & (data['pmRA'] < 30) &\
        (data['pmDE'] < 30) & (d_col < d_max)
    pmRA, pmDE, epmRA, epmDE = data['pmRA'][msk3], data['pmDE'][msk3],\
        data['e_pmRA'][msk3], data['e_pmDE'][msk3]
    ax3.set_title("N(r<d_max, |pmX|<30)={}".format(pmRA.size), fontsize=8)
    cmap = cm.viridis
    norm = Normalize(vmin=d_col[msk3].min(), vmax=d_max)
    ax3.errorbar(
        pmRA, pmDE, yerr=epmDE, xerr=epmRA, fmt='none', elinewidth=.35,
        ecolor=cmap(norm(d_col[msk3])))
    # pmRA_m, pmRA_s, pmDE_m, pmDE_s = np.nanmean(pmRA), np.nanstd(pmRA),\
    #     np.nanmean(pmDE), np.nanstd(pmDE)
    # plt.xlim(pmRA_m - 2. * abs(pmRA_s), pmRA_m + 2. * abs(pmRA_s))
    # plt.ylim(pmDE_m - 2. * abs(pmDE_s), pmDE_m + 2. * abs(pmDE_s))
    im = plt.scatter(pmRA, pmDE, s=0, c=d_col[msk3], cmap=cmap)
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


if __name__ == '__main__':
    main()
