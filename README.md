This is a Python program to design an optimized layout of a (Concentrated Solar
Power (CSP) Plant.  Specifically, it designs the layout of the heliostats
(mirrors) in the solar field, which reflect sunlight up the the receiver tower.

Several parameters can be defined at the top of the script, and these are used
to achieve a desired net power output of the plant at a specified time.
(NOTE: Currently the time is set to be *solar noon*, but the day of the year can
be set to any day of the year.)

# ASSUMPTIONS
## NOMINAL CONDITIONS TO BE USED:
    - solar noon on autumnal equnox     --> Day 80, t_s = 12.0
    - highest measured I_b of the year  --> 1089 [W/m2]
    - visibility loss coefficients used (taken from a lecture):
        a0 = 0.679 / 100  [-]
        a1 = 10.46 / 100  [-/km]
        a2 = -1.70 / 100  [-/km^2]
        a3 = 0.0   / 100  [-/km^3]

## LAND AND TOPOGRAPHY
    - topography assumed perfectly flat
    - curvature of Earth is neglected in terms of altitude (z-axis) change
        --> only about 0.2m difference
            40000*(1-sin(arccos(4/40000))) = 0.00019km (0.19m)
    - maintenance roads are not included
    - power transmission cables are not included
    - (only the solar field is designed.  Power block etc. are not included)

## HELIOSTATS
    - mirror area / total area = 0.96
    - mirror clean factor = 0.95
    - shadowing & blocking --> assumed to be 5% (combined)
    - mirrors are assumed flat (no focusing, therefore there is excess spillage)

## RECEIVER
    - Optical height is center of receiver height
    - receiver is cylindrical       --> Diameter = 10 [m] (horizontal)
                                    --> Length   = 15 [m] (vertical)
    - spillage is approximated as purely horizontal, and straight on
        * since receiver is cylinderical, it appears the same from
          all around, if you neglect vertical changes
    - receiver length is assumed long enough to avoid vertical spillage
    - all receiver losses are assumed to be included the receiver efficiency

## POWER
    - parasitic losses (self-consumption) are not considered
        - molten salt pumps
        - heliostat tracking system
        - monitoring system
    - energy losses from molten salt storage system are neglected
