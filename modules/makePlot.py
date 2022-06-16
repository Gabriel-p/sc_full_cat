
import warnings
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import numpy as np
from . import spiralArms


def plot(dpi, mode, crossMdata, gc_frame):
    """
    Gridspec idea: http://www.sc.eso.org/~bdias/pycoffee/codes/20160407/
                   gridspec_demo.html
    """

    # Ignore colorbar warning
    warnings.filterwarnings(
        "ignore", category=UserWarning, module="matplotlib")
    # Ignore RuntimeWarning when creating the masks
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if mode == 'z_dist':
        # Vertical distance masks.
        mnan = np.isnan(crossMdata['z_pc'])
        m200 = abs(crossMdata['z_pc']) <= 200
        m600 = (200 < abs(crossMdata['z_pc'])) & (
            abs(crossMdata['z_pc']) <= 600)
        m1000 = (600 < abs(crossMdata['z_pc'])) & (
            abs(crossMdata['z_pc']) <= 1000)
        m1500 = (1000 < abs(crossMdata['z_pc'])) & (
            abs(crossMdata['z_pc']) <= 1500)
        m2500 = (1500 < abs(crossMdata['z_pc'])) & (
            abs(crossMdata['z_pc']) <= 2500)
        minf = abs(crossMdata['z_pc']) > 2500
        plt_data = {
            '0': [20, .5, 'x', 'k', r'$No\,dist\; (N={})$'],
            '1': [20., .3, 'o', 'grey', r'$z\leq 200\,[pc]\; (N={})$'],
            '2': [25., .3, 'o', 'red', r'$200<z\leq600\,[pc]\;(N={})$'],
            '3': [25., .3, 'o', 'orange', r'$600<z\leq1000\,[pc]\;(N={})$'],
            '4': [30., .4, 'o', 'green', r'$1000<z\leq1500\,[pc]\;(N={})$'],
            '5': [50., .5, 'o', 'cyan', r'$1500<z\leq2500\,[pc]\;(N={})$'],
            '6': [70., .6, 'o', 'blue', r'$z>2500\,[pc]\;(N={})$']
        }
        mask_list = [mnan, m200, m600, m1000, m1500, m2500, minf]

    elif mode == 'd_dist':
        mnan = np.isnan(crossMdata['z_pc'])
        m1000 = abs(crossMdata['dist_pc']) <= 1000
        m2000 = (1000 < abs(crossMdata['dist_pc'])) & (
            abs(crossMdata['dist_pc']) <= 2000)
        m3000 = (2000 < abs(crossMdata['dist_pc'])) & (
            abs(crossMdata['dist_pc']) <= 3000)
        m4000 = (3000 < abs(crossMdata['dist_pc'])) & (
            abs(crossMdata['dist_pc']) <= 4000)
        m5000 = (4000 < abs(crossMdata['dist_pc'])) & (
            abs(crossMdata['dist_pc']) <= 5000)
        m6000 = (5000 < abs(crossMdata['dist_pc'])) & (
            abs(crossMdata['dist_pc']) <= 6000)
        minf = abs(crossMdata['dist_pc']) > 6000
        plt_data = {
            '0': [20, .5, 'x', 'k', r'$No\,dist\; (N={})$'],
            '1': [25., .4, 'o', 'grey', r'$d\leq 1000\,[pc]\; (N={})$'],
            '2': [25., .4, 'o', 'red', r'$1000<d\leq2000\,[pc]\;(N={})$'],
            '3': [25., .4, 'o', 'orange', r'$2000<d\leq3000\,[pc]\;(N={})$'],
            '4': [25., .4, 'o', 'green', r'$3000<d\leq4000\,[pc]\;(N={})$'],
            '5': [25., .4, 'o', 'cyan', r'$4000<d\leq5000\,[pc]\;(N={})$'],
            '6': [25., .4, 'o', 'blue', r'$5000<d\leq6000\,[pc]\;(N={})$'],
            '7': [25., .4, 'o', 'purple', r'$d>6000\,[pc]\;(N={})$']
        }
        mask_list = [mnan, m1000, m2000, m3000, m4000, m5000, m6000, minf]

    fig = plt.figure(figsize=(25, 25))
    gs = gridspec.GridSpec(
        6, 6, height_ratios=[1, 1, 1, 0.04, .7, .7],
        width_ratios=[1, 1, 1, 1, 1, 1])
    gs.update(hspace=0.03, wspace=.3)

    ax = plt.subplot(gs[0:3, 0:6], projection="aitoff")
    plt.title('Database N={}'.format(len(crossMdata)), y=1.02)
    ax.set_xticklabels([
        r'210$^{\circ}$', r'240$^{\circ}$', r'270$^{\circ}$', r'300$^{\circ}$',
        r'330$^{\circ}$', r'0$^{\circ}$', r'30$^{\circ}$', r'60$^{\circ}$',
        r'90$^{\circ}$', r'120$^{\circ}$', r'150$^{\circ}$'], fontsize=18)
    plt.yticks(fontsize=18)
    ax.grid(True)

    cl_plots1, cl_plots2 = [[], []], [[], []]
    for i, m in enumerate(mask_list):
        # Only plot if there are objects to plot.
        if sum(m) > 0:
            ms, a, mrk, c, lab = plt_data[str(i)]

            pl = ax.scatter(
                crossMdata['lon'][m], crossMdata['lat'][m], marker=mrk, s=ms,
                alpha=a, c=c, label=lab.format(len(crossMdata['lon'][m])))

            if i in [0, 1, 2, 3]:
                cl_plots1[0].append(pl)
                cl_plots1[1].append(lab.format(len(crossMdata['lon'][m])))
            else:
                cl_plots2[0].append(pl)
                cl_plots2[1].append(lab.format(len(crossMdata['lon'][m])))

            # # Name annotate
            # lon_lat_cust = (
            #     (227.33827932, -8.74737685), (132.14541698, -8.73629164),
            #     (147.48401456, -1.47915429), (242.89042522, -6.86623902),
            #     (242.74214783, 4.89128908), (177.89050533, -0.43351877))
            # name_cust = ('GAIA1', 'GAIA2', 'GAIA4', 'GAIA5', 'GAIA6', 'GAIA7')

            # if db_read == 'custom':
            #     for _, (lon, lat) in enumerate(lon_lat_cust):

            #         lb = SkyCoord(lon, lat, unit='degree', frame='galactic')
            #         rad_l = lb.l.wrap_at(180 * u.deg).radian
            #         rad_b = lb.b.radian

            #         ax.scatter(rad_l, rad_b, marker='o', s=50, c='k', zorder=5)
            #         xshift = -0.25 if name_cust[_] == 'GAIA7' else 0.015
            #         yshift = -0.05 if name_cust[_] == 'GAIA4' else 0.015
            #         ax.annotate(
            #             name_cust[_], (rad_l + xshift, rad_b + yshift),
            #             xycoords='data', color='b', fontsize=16)

            if mode == 'z_dist':
                # Only plot names for those with the two largest 'z' values.
                fs = [6, 8, 10]
                if i in [4, 5, 6]:
                    for _, (lon, lat) in enumerate(
                            zip(*[
                                crossMdata['lon'][m], crossMdata['lat'][m]])):
                        ax.annotate(crossMdata['name'][m][_].split(',')[0],
                                    (lon, lat), xycoords='data',
                                    fontsize=fs[i - 4])

    l1 = plt.legend(cl_plots1[0], cl_plots1[1], loc=1, fontsize=12)
    plt.legend(cl_plots2[0], cl_plots2[1], loc=4, fontsize=12)
    plt.gca().add_artist(l1)

    plt.style.use('seaborn-darkgrid')

    crossMdata['x_pc'].convert_unit_to('kpc')
    crossMdata['y_pc'].convert_unit_to('kpc')
    crossMdata['z_pc'].convert_unit_to('kpc')

    # Sun's coords according to the Galactocentric frame.
    x_sun, z_sun = gc_frame.galcen_distance, gc_frame.z_sun
    s_xys = SkyCoord(
        -x_sun, 0., z_sun, unit='kpc', representation_type='cartesian')

    ax = plt.subplot(gs[4:6, 0:2])
    plt.xlabel(r"$x_{GC}\, [kpc]$", fontsize=12)
    plt.ylabel(r"$y_{GC}\, [kpc]$", fontsize=12)
    vmin, vmax = max(min(crossMdata['z_pc']), -2.5),\
        min(max(crossMdata['z_pc']), 2.5)
    plt1 = plt.scatter(
        crossMdata['x_pc'], crossMdata['y_pc'].data, alpha=.5,
        c=crossMdata['z_pc'], cmap='viridis', vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.x, s_xys.y, c='yellow', s=50, edgecolor='k')
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    plt.scatter(0., 0., c='k', marker='x', s=70)
    # Plot spiral arms
    spiral_arms = spiralArms.momany()
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
    plt.legend(fontsize=10)
    # colorbar
    cbax = plt.subplot(gs[3:4, 0:2])
    cb = Colorbar(
        ax=cbax, mappable=plt1, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < z_{{GC}}\, [kpc] < {:.1f}$".format(
        min(crossMdata['z_pc']), max(crossMdata['z_pc'])), labelpad=10)

    ax = plt.subplot(gs[4:6, 2:4])
    plt.xlabel(r"$x_{GC}\, [kpc]$", fontsize=12)
    plt.ylabel(r"$z_{GC}\, [kpc]$", fontsize=12)
    vmin, vmax = max(ymin, -15.), min(ymax, 15.)
    plt2 = plt.scatter(
        crossMdata['x_pc'], crossMdata['z_pc'], alpha=.5, c=crossMdata['y_pc'],
        cmap='viridis', vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.x, s_xys.z, c='yellow', s=50, edgecolor='k')
    plt.scatter(0., 0., c='k', marker='x', s=70)
    plt.xlim(max(xmin, -20.), min(xmax, 20.))
    plt.ylim(max(min(crossMdata['z_pc']), -3.),
             min(max(crossMdata['z_pc']), 3.))
    # colorbar
    cbax = plt.subplot(gs[3:4, 2:4])
    cb = Colorbar(
        ax=cbax, mappable=plt2, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < y_{{GC}}\, [kpc] < {:.1f}$".format(
        min(crossMdata['y_pc']), max(crossMdata['y_pc'])), labelpad=10)

    ax = plt.subplot(gs[4:6, 4:6])
    plt.xlabel(r"$y_{GC}\, [kpc]$", fontsize=12)
    plt.ylabel(r"$z_{GC}\, [kpc]$", fontsize=12)
    vmin, vmax = max(xmin, -20.), min(xmax, 20.)
    plt3 = plt.scatter(
        crossMdata['y_pc'], crossMdata['z_pc'], alpha=.5, c=crossMdata['x_pc'],
        cmap='viridis', vmin=vmin, vmax=vmax)
    plt.scatter(s_xys.y, s_xys.z, c='yellow', s=50, edgecolor='k')
    plt.xlim(max(ymin, -15.), min(ymax, 15.))
    plt.ylim(max(min(crossMdata['z_pc']), -3.),
             min(max(crossMdata['z_pc']), 3.))
    # colorbar
    cbax = plt.subplot(gs[3:4, 4:6])
    cb = Colorbar(
        ax=cbax, mappable=plt3, orientation='horizontal', ticklocation='top')
    cb.set_label(r"${:.1f} < x_{{GC}}\, [kpc] < {:.1f}$".format(
        min(crossMdata['x_pc'].data), max(crossMdata['x_pc'].data)),
        labelpad=10)

    plt.suptitle(r'$[d_{{GC}}={},\;\;z_{{\odot}}={}]$'.format(
        x_sun, z_sun), x=.52, y=.4, fontsize=14)

    fig.tight_layout()
    fig.savefig('output/crossMdata.png', dpi=dpi, bbox_inches='tight')
    plt.style.use('default')


def checkInRange(xy_data, xmin, xmax, ymin, ymax):
    """
    Check if at least one point is between the given ranges.
    """
    for x, y in xy_data:
        if xmin <= x <= xmax and ymin <= y <= ymax:
            return True
    return False
