
# Gaia DR2 questions:

2. Is there a completeness curve for the photometry?


# Gaia 1 & Gaia 2

1. Koposov et al. (2017)

 * DR1 + algorithm for detection --> Gaia 1 & Gaia 2
 * Analysis: combine Gaia 'G' + 2MASS 'JHK' + Pan-STARRS1 'rz'

## Gaia 1
 * (alpha=102.47, delta=-16.75) ; (l=227.34, b=-8.75)
 * Was (is) hidden by Sirius
 * [Fe/H]=-0.7 +- 0.2; age=6.3 +- 2 Gyr (~9.8); d=4.6 Kpc (13.3);
   M=22000 Mo (ball-park estimates)
 * Most likely globular (given its age and metallicity) or an old open cluster.
 * Gaia 1 --> z ~ -640 pc

## Gaia 2 
 * [Fe/H]=-0.6 +- 0.2; age=8 +- 2 Gyr (~9.9); d=5.2 Kpc (13.6)
 * Most likely globular (given its age and metallicity) or an old open cluster.
 * Gaia 2 --> z ~ -760 pc


# Gaia 1

2. Simpson et al. (2017)

 * Spectral analysis (HERMES + AAOmega) --> confirm Gaia 1
 * Radial velocity --> membership analysis
 * [Fe/H]=-0.13 +- 0.13; age=3 Gyr (~9.5); d=4.46 Kpc (13.2), E(B-V)~0.66
 * Dynamical mass ~ 13000 Mo
 * R_GC~11.8 Kpc; z=-640 pc
 * Orbital analysis using 'galpy' + 'MWPotential2014'
 * Forward 1 Gyr --> z_max~1.1 Kpc (+0.4, -0.3 Kpc)
 * "In its present orbit, it makes nine plane crossings every gigayear, for a
    total of over 30 in its 3 Gyr lifetime"
 * "surprising to find Gaia 1 in its present orbit at the present day"


3. Mucciarelli et al. (2017)

 * Observed 6 He-clump stars with MIKE/MAGELLAN spectrograph
 * [Fe/H]~0.00 +- 0.01
 * Rejects possible extra-galactic origin, favors a Galactic open cluster
 * "S17 stars do not define an RGB in the theoretical plane, suggesting that
    their parameters are not correct"
 * Atmospheric parameters from Simpson et al. (2017) are **wrong**
 * "Unremarkable standard Galactic open cluster"
 * "the precise value of R_GC can be affected by the value of E(B-V), which
    still remains uncertain for this cluster."
 * "orbital parameters inferred by S17 may suggest that the cluster has been
    accreted by the Milky Way, they are still, within the uncertainties, fully 
    compatible with the majority of known Galactic open clusters"
 * Poor proper motions --> Uncertainty in the orbit --> wait for Gaia DR2


4. Koch et al. (2017)

 * Chemical abundances for 4 red giant members
 * [Fe/H]~-0.62 +- 0.13
 * Dynamical mass: 30000 - 180000 Mo --> "dwarf galaxy"??
 * galpy + MWPotential2015 backwards 10 Gyr --> z_max~1 Kpc for the 4 stars
 * alpha abundances compatible with the thick disk
 * "distance to the cluster is poorly constrained"
 * "firmly establish Gaia 1 as an object associated with the thick disk"
 * "unclear which mechanisms put it in that place"


5. Carraro (2018) (RN AAS)

 * "(Koch et al. 2017) conclude that Gaia 1 is a thick disk cluster because
    of its actual location, almost 1 kpc above the formal (b=0°) Galactic
    plane" <-- Koch et al. do not state this
 * Carraro says: for the *real* (warped) disk: z<0.5 Kpc
 * NOT a thick disk cluster
