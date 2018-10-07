
import numpy as np
# from astropy.coordinates import Angle
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
# import astropy.units as u
# import astropy.coordinates as coord
# from astroquery.vizier import Vizier


# catalog_list = Vizier.find_catalogs('Pan-STARRS')
# catalogs = Vizier.get_catalogs(catalog_list.keys())
# print(catalogs)

# # Select catalog
# #
# GaiaDR2 = 'I/345/gaia2'
# Pan_STARRS1 = 'II/349/ps1'
# #
# cat = Pan_STARRS1

# # Select cluster
# #
# # name = 'GAIA1'
# # center, box_s = (101.47, -16.75), "1deg"
# name = 'GAIA1_ps1'
# center, box_s = (101.47, -16.75), "1deg"
# # name = 'NGC6791'
# # center, box_s = (290.2208333, 37.7716667), "1deg"
# # name = 'RUP44_ps1'
# # center, box_s = (119.7125, -28.5833), "1deg"

# # Unlimited rows, all columns
# v = Vizier(row_limit=-1, columns=['all'])
# result = v.query_region(coord.SkyCoord(
#     ra=center[0], dec=center[1], unit=(u.deg, u.deg), frame='icrs'),
#     width=box_s, catalog=[cat])

# ascii.write(result[cat], 'input_expl/' + name + ".dat", overwrite=True)

#
# Format data file
#
data = ascii.read('input_expl/GAIA1_gaia2_ps1.dat')

# colors
gr = MaskedColumn(data['gmag'] - data['rmag'], name='g-r')
e_gr = MaskedColumn(data['e_gmag'] + data['e_rmag'], name='e_gr')
rz = MaskedColumn(data['rmag'] - data['zmag'], name='r-z')
e_rz = MaskedColumn(data['e_rmag'] + data['e_zmag'], name='e_rz')
iz = MaskedColumn(data['imag'] - data['zmag'], name='i-z')
e_iz = MaskedColumn(data['e_imag'] + data['e_zmag'], name='e_iz')

data2 = Table([
    data['DR2Name'], data['RA_ICRS'], data['DE_ICRS'], data['Gmag'],
    data['e_Gmag'], gr, e_gr, rz, e_rz, iz, e_iz, data['Plx'],
    data['e_Plx'], data['pmRA'], data['e_pmRA'], data['pmDE'],
    data['e_pmDE'], data['RV'], data['e_RV']])
import pdb; pdb.set_trace()  # breakpoint 2348d463 //

ascii.write(data2, 'input_expl/GAIA1_2.dat')
