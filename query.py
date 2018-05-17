
# from astropy.coordinates import Angle
from astropy.io import ascii
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier


# catalog_list = Vizier.find_catalogs('Pan-STARRS')
# catalogs = Vizier.get_catalogs(catalog_list.keys())
# print(catalogs)

# Select catalog
#
GaiaDR2 = 'I/345/gaia2'
Pan_STARRS1 = 'II/349/ps1'
#
cat = Pan_STARRS1

# Select cluster
#
# name = 'GAIA1'
# center, box_s = (101.47, -16.75), "1deg"
# name = 'NGC6791'
# center, box_s = (290.2208333, 37.7716667), "1deg"
name = 'RUP44_ps1'
center, box_s = (119.7125, -28.5833), ".1deg"

# Unlimited rows, all columns
v = Vizier(row_limit=-1, columns=['all'])
result = v.query_region(coord.SkyCoord(
    ra=center[0], dec=center[1], unit=(u.deg, u.deg), frame='icrs'),
    width=box_s, catalog=[cat])

ascii.write(result[cat], 'input/' + name + ".dat")
