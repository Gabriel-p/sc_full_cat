
import numpy as np
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import matplotlib.pyplot as plt


def stats(db, zc, ac, ec, dc):
    """
    """
    Nt = len(db)
    N_zaed = []
    for c in (zc, ac, ec, dc):
        if c == '?':
            N_zaed.append(0)
        else:
            try:
                N_zaed.append((~db[c].mask).sum())
            except:
                N_zaed.append(len(db))
    N_zaed = np.array(N_zaed) / Nt

    print("N={}, z={:.2f}, a={:.2f}, e={:.2f}, d={:.2f}".format(Nt, *N_zaed))


def main():
    """
    """

    db = ascii.read('input/CANTAT_GAUDIN_2018.dat')
    print('CANTAT_GAUDIN_2018')
    stats(db, '?', 'AgeNN', '?', 'dmode')

    db = ascii.read('input/CANTAT_GAUDIN.dat')
    print('CANTAT_GAUDIN')
    stats(db, '?', 'AgeNN', 'AVNN', 'DistPc')

    db = ascii.read('input/DIAS.dat')
    print('DIAS')
    stats(db, '[Fe/H]', 'logage', 'Av', 'Dist')

    db = ascii.read('input/MWSC.dat', format='fixed_width')
    print('MWSC')
    stats(db, 'metallicity', 'log_age', 'e_bv', 'distance')

    db = ascii.read('input/OPENCLUST.dat', delimiter=';')
    print('OPENCLUST')
    stats(db, '[Fe/H]', 'Age', 'E(B-V)', 'Dist')

    db = ascii.read('input/WEBDA.dat')
    print('WEBDA')
    stats(db, 'Fe/H', 'Age', 'EB-V', 'Dist')

    db = ascii.read('input/DIAS.dat')
    # Galactic coordinates.
    eq = SkyCoord(ra=db['RA_ICRS'] * u.deg, dec=db['DE_ICRS'] * u.deg, frame='icrs')
    lb = eq.transform_to('galactic')
    db['lon'] = lb.l.wrap_at(180 * u.deg).radian * u.radian
    db['lat'] = lb.b.radian * u.radian
    coords = SkyCoord(
        l=db['lon'], b=db['lat'],
        distance=db['Dist'] * u.pc, frame='galactic')
    # Galactocentric coordinates.
    c_glct = coords.transform_to(coord.Galactocentric())
    R_GC = np.sqrt(c_glct.x**2 + c_glct.y**2 + c_glct.z**2)
    plt.scatter(R_GC, db['[Fe/H]'])
    # Trend taken from Donor et al. (2020), Fig 7
    # https://ui.adsabs.harvard.edu/abs/2020AJ....159..199D/abstract
    plt.plot((8000., 13900), (0., -0.4), ls=':', lw=2, c='k')
    plt.plot((13900, 22000), (-0.4, -0.473), ls=':', lw=2, c='k')
    plt.show()


if __name__ == '__main__':
    main()
