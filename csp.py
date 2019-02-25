#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from operator import attrgetter

# make some references to np functions, for ease of coding
pi  = np.pi
cos = np.cos
sin = np.sin
tan = np.tan
arccos = np.arccos
arcsin = np.arcsin
arctan = np.arctan

# make a global variable for holding the net power
net_power_thermal = 0.0
# global list for all heliostats in solar field
solar_field = []
# global variable to hold number of rows
rows = 0

print("============================================")

#-------------------------------------------------------------------------------
# Design a CSP plant in Chile
#-------------------------------------------------------------------------------
# (d) For a chosen design point, propose a configuration:
#       - tower height
#       - number of heliostats
#       - layout/ heliostat placement
#       -->  (Consider the available land and land shape)
#
# (e) Estimate total land area required
#
# (f) Estimate field costs, assuming component costs can be scaled linearly
#
# (g) Consider now including an 8-hr storage system
#     Repropose the configuration.  List differences between both
#     configurations (with and without storage)
#-------------------------------------------------------------------------------

#===============================================================================
#-------------------------------- ASSUMPTIONS ----------------------------------
"""
NOMINAL CONDITIONS TO BE USED:
    - solar noon on summer solstice     --> Day 255
                                        --> Hour 12:xx  *(correct for meridian)
    - highest measured I_b of the year  --> 1089 [W/m2]
    - highest temp of the year          --> 294.55 [K]
    - highest pressure of the year      --> 76.98 [kPa]
    - visibility (atmosphere)           --> 15 km

GEOMETRY
    - topography assumed perfectly flat
    - curvature of Earth is neglected in terms of altitude (z-axis) change
        --> only about 0.2m difference
            40000*(1-sin(arccos(4/40000))) = 0.00019km (0.19m)
    - maintenance roads are not accounted for

HELIOSTATS
    - mirror area / total area = 0.96
    - mirror clean factor = 0.95
    - shadowing & blocking --> assumed to be 5%
    - mirrors are flat

RECEIVER
    - Optical height is center of receiver height
    - receiver is cylindrical       --> Diameter = 10 [m]
                                    --> Length   = 15 [m] (vertical)
    - spillage is approximated as purely horizontal, and straight on
        * since receiver is cylinderical, it appears the same from
          all around, if you neglect vertical changes
    - receiver length is assumed long enough to avoid vertical spillage

POWER
    - parasitic loss (self-consumption) is not considered
    - heat loss is not considered
"""
#===============================================================================
# GIVEN PARAMETERS

NET_POWER_MW            = 50                 #[MW]
NET_POWER_KW            = NET_POWER_MW * 1000 #[kW]
NET_POWER_W             = NET_POWER_KW * 1000 #[W]

STORAGE_HOURS           = 8 #[hrs]

EFF_PB                  = 0.43
EFF_RECEIVER_THERMAL    = 0.88
EFF_GENERATOR           = 0.98

NET_POWER_TH_W          = NET_POWER_W / ( EFF_PB
                                        * EFF_RECEIVER_THERMAL
                                        * EFF_GENERATOR )
                        
MIRROR_AREA             = 115  #[m2]
MIRROR_AREA_TO_TOTAL    = 0.96 # (assumption)
MIRROR_REFLECTIVITY     = 0.95
MIRROR_CLEAN_FACTOR     = 0.95 # (assumption)
                        
TOWER_HEIGHT_TOTAL      = 250 #[m]
TOWER_HEIGHT_OPTICAL    = 220 #[m]

RECEIVER_DIAMETER       = 10 #[m] # (assumption)
RECEIVER_LENGTH         = 15 #[m] # (assumption)
                        
INNER_RADIUS            = 100 #[m]
                        
COST_HELIOSTAT          = 140    #[$/m2]
COST_RECEIVER           = 100    #[$/kW_thermal]
COST_TOWER              = 170000 #[$/m]
COST_STORAGE            = 25     #[$/kWh_thermal]


#===============================================================================
# Define the nominal conditions
NOM_DAY         = 173#355       #[day]
NOM_T_S         = 12        #[hr] (solar time)
NOM_TEMP        = 294.55    #[K]
NOM_HOUR        = 12.0      #[hr]   TODO: adjust for meridian difference
NOM_Ib          = 1089      #[W/m2]
NOM_Ig          = 1277      #[W/m2]
NOM_Id          = 547.0879  #[W/m2]
NOM_PRESSURE    = 76.98     #[kPa]
NOM_VISIBILITY  = 15000     #[m]

SUN_ANGULAR_DIAMETER  = 0.5334 #[deg]


#===============================================================================
# Define the polygonal area we can build on
#Point      East (x)    North (y)   Altitude (z)
P1      = ( -68.73613 , -22.57351 , 0 )     # most South
P2      = ( -68.69902 , -22.54856 , 0 )     # most East
P3      = ( -68.70427 , -22.53609 , 0 )
P4      = ( -68.7691  , -22.51801 , 0 )     # most North
P5      = ( -68.7717  , -22.55289 , 0 )     # most West
LAT_LON_POINTS  = [P1, P2, P3, P4, P5]

# determine "center" of polygon
p_sum = (0,0,0)
for p in LAT_LON_POINTS:
    p_sum = (p[0] + p_sum[0], p[1] + p_sum[1], p[2] + p_sum[2])
LAT_LON_CENTER = tuple(np.divide(p_sum, len(LAT_LON_POINTS)))
LONGITUDE = LAT_LON_CENTER[0]   # (x)
LATITUDE  = LAT_LON_CENTER[1]   # (y)
ALTITUDE  = LAT_LON_CENTER[2]   # (z)

MERIDIAN  = -3*15 #[deg E]


#-------------------------------------------------------------------------------
# translate lat/long points into (x,y) coordinates [m]
def lat_lon_to_dx_dy(p1, p2):
    (lon1, lat1, alt1) = p1
    (lon2, lat2, alt2) = p2
    dx = 1000 * (lon2-lon1)*40000*cos((lat1+lat2)*pi/360)/360
    dy = 1000 * (lat2-lat1)*40000/360
    dz = 1000 * (alt2-alt1)
    return (dx, dy, dz)

print("Perimeter (x,y) coordinates:")
perimeter_xyz_points = []
for p in LAT_LON_POINTS:
    (dx,dy,dz) = lat_lon_to_dx_dy(LAT_LON_CENTER, p)
    perimeter_xyz_points.append((dx,dy,dz))
    print("    ({: 6.0f}, {: 6.0f} ) [m]".format(dx,dy))
xyz1, xyz2, xyz3, xyz4, xyz5 = perimeter_xyz_points
print()


#-------------------------- Plot Property Perimeter ----------------------------
# create plot (turn on interactive mode)
plt.ion()

# plot perimeter points
pp = list(perimeter_xyz_points)
pp.append(perimeter_xyz_points[0]) # append [0] element, so it draws full perimeter
ppx, ppy, ppz = zip(*pp)
plot = plt.scatter(ppx, ppy, color='k', marker='')

# draw perimeter border
plt.plot(ppx, ppy, 'k-', color='k', linewidth=1)

# plot center point (tower position)
ppx, ppy, ppz = [0], [0], [0]
plt.scatter(ppx, ppy, color='r', marker='o', s=50)
# plt.scatter([1000], [1000], color='g')

# show plot
plt.show(block=False)


#------------------------------ Solar parameters -------------------------------
lat = np.deg2rad(LATITUDE)
lon = np.deg2rad(LONGITUDE)
alt = np.deg2rad(ALTITUDE)

t_s = NOM_T_S           # solar time
n = NOM_DAY             # day of the year
w = (pi/12)*(t_s - 12)  # hour angle

# solar angles
declination = arcsin( 0.39795 * cos(2*pi*(n-173)/365) )
zenith      = arccos( cos(lat)*cos(declination)*cos(w) +
                      sin(lat)*sin(declination) )
elevation = (pi/2) - zenith

sgn     = w/abs(w) if w else 0
azimuth = sgn * abs( arccos(
                        (cos(zenith)*sin(lat) - sin(declination)) /
                        (sin(zenith)*cos(lat))
                ))

# length of day
N = (24/pi)*arccos( -1*tan(lat) * tan(declination) )

# Convert solar unit vector to cartesian coordinates
# TODO: Extend this to take into account times other than solar noon
# v_S = np.array(( sin(zenith)*sin(azimuth),
#                  sin(zenith)*cos(azimuth),
#                  cos(zenith) ))
v_S = np.array((0, sin(zenith), cos(zenith)))


print("Solar parameters:")
print("    t_s:          {: } hour".format(t_s))
print("    w:            {: 04.2f} deg".format(w))
print("    n:            {: } (day of year)".format(n))
print("    N:            {: 03.1f} hours (length of day)".format(N))
print("    declination:  {: 04.2f} deg".format(np.degrees(declination)))
print("    zenith:       {: 04.2f} deg".format(np.degrees(zenith)))
print("    azimuth:      {: 04.2f} deg".format(np.degrees(azimuth)))
print("    elevation:    {: 04.2f} deg".format(np.degrees(elevation)))
print()


#===============================================================================
class Heliostat(object):
    """
    |-------| Tower height                                      
    |       |                                                   
    |#######|                                                   
    |## + ##|- Receiver height - - - - - - - - - - - - - - - - -
    |#######|  (optical height)                            |    
    |       |                       elevation_to_receiver  | 
    |       |                                              |    
    |       |                                              |    
    |       |                      width                   |    
    |       |                -----------------             |    
    |       |                |               |             |     
    |       |                |      area     | depth   - - - - -
    |       |                |               |             |
    |       |                -----------------             | height
    |       |                        |                     |
    |-------|                       / \    _ _ _ _ _ _ _ _ _ _ _
    """
    area         = MIRROR_AREA * MIRROR_AREA_TO_TOTAL
    reflectivity = MIRROR_REFLECTIVITY
    clean_factor = MIRROR_CLEAN_FACTOR
    depth        = np.floor(np.sqrt(area))
    width        = area / depth
    height       = width / 2    # gives some space between ground and mirror
                                # since width is greater than depth
    elevation_to_receiver = np.array((0, 0, TOWER_HEIGHT_OPTICAL - height))

    def __init__(self, position):
        """
        param: position [np.array] - (x, y, z=0) coordinate from tower base
        """
        # vector from tower base
        self.position = position    # (x,y,z) coordinate, with tower at (0,0,0)

        # initiate vector attributes
        self.vT  = np.ones(3)   # vector to tower (receiver)
        self.nH  = np.ones(3)   # normal vector of mirror surface
        self.v_T = np.ones(3)   # UNIT vector to tower (receiver)
        self.n_H = np.ones(3)   # UNIT normal vector of mirror surface
        self.dR  = 0.0          # distance to receiver [m]

        # angles
        self.slope_angle    = 0.0
        self.azimuth_angle  = 0.0
        self.theta          = 0.0

        # loss factors
        self.cosine_eff     = 1.0
        self.f_shadow_block = 0.0
        self.f_att          = 0.0
        self.f_spill        = 0.0

        # total energy contributed to the receiver
        self.total_contribution = 0.0   # (total W_th to receiver)

        # calculate total contribution
        self.calculate_energy_contribution()

    def calculate_n_H(self):
        """calculate the surface normal vector"""
        self.nH = v_S + self.v_T
        self.n_H = self.nH / np.linalg.norm(self.nH)

    def set_angle(self):
        """Calculate angle to receiver, then surface normal vector"""
        # calculate vector to receiver
        self.vT = Heliostat.elevation_to_receiver - self.position
        # dR = magnitude(self.vT)
        self.dR = np.linalg.norm(self.vT)
        # convert to unit vector
        self.v_T = self.vT / self.dR
        # calculate the surface normal vector
        self.calculate_n_H()

    def determine_cosine_effectiveness(self):
        # (same as theta)
        self.cosine_eff = np.dot(self.v_T, self.n_H)
        self.theta = self.cosine_eff

    def determine_shadow_block(self):
        # (assumed to be ~5%, for simplicity)
        self.f_shadow_block = 0.05

    def determine_attenuation(self):
        """
        Use one of two approaches:
            (1) Visibility distance [km]
            (2) Loss coefficients (ex: [%/km])
        """
        # (1)
        #--------------- Visibility Approach ---------------
        # attenuation due to visibility in atmosphere
        # self.f_att = 1.0 - np.exp(-(3*self.dR) / NOM_VISIBILITY)

        # (2)
        #----------- Loss Coefficiencts Approach -----------
        # visibility loss coefficients
        a0 = 0.679 / 100  #[-]
        a1 = 10.46 / 100  #[-/km]
        a2 = -1.70 / 100  #[-/km^2]
        a3 = 0.0   / 100  #[-/km^3]

        # attenuation using visibility coefficients
        dR = self.dR / 1000
        self.f_att = (   a0
                       + a1 * (dR)
                       + a2 * (dR**2)
                       + a3 * (dR**3)
                     )

    def determine_spillage(self):
        # convert solar angular spread to radians
        delta_s = np.deg2rad(SUN_ANGULAR_DIAMETER/2)
        
        # rename variables to math equation
        AH = Heliostat.width
        AR = RECEIVER_DIAMETER

        # calculate A.I (width of reflected light at receiver)
        cos_theta_M = np.dot(self.n_H, self.v_T)
        AI = AH * cos_theta_M + (2*self.dR * sin(delta_s))

        # calculate f_spill
        if AI > AR:
            self.f_spill = 1.0 - (AR/AI)
        else:
            self.f_spill = 0.0

    def calculate_energy_contribution(self, verbose=False):
        if verbose:
            print(".................................................")
            print("calculating heliostat parameters at position:")
            print("  ({},{},{}) [m]".format(*self.position), end="")

        # calculate angle for aiming to receiver
        self.set_angle()
        if verbose:
            print("  --> dR = {:,.1f} [m]".format(self.dR))
            print(".................................................")

        # calculate all efficiency parameters (losses)
        self.determine_cosine_effectiveness()
        self.determine_shadow_block()
        self.determine_attenuation()
        self.determine_spillage()

        #-----------------------------------------------------------------------
        if verbose:
            print("=============================")
            print("NOM_Ib: {}  [W]".format(NOM_Ib))
            print("area:   {:04.1f} [m2]".format(Heliostat.area))
            print("    _________________________")
            print("    reflectivity:    {: 04.3f}".format(Heliostat.reflectivity))
            print("    clean_factor:    {: 04.3f}".format(Heliostat.clean_factor))
            print("    cosine_eff:      {: 04.3f}".format(self.cosine_eff))
            print("    f_shadow_block:  {: 04.3f}".format(self.f_shadow_block))
            print("    f_att:           {: 04.3f}".format(self.f_att))
            print("    f_spill:         {: 04.3f}".format(self.f_spill))
            print("    _________________________")

        # calculate total efficiency of heliostat
        total_efficiency = (
              Heliostat.reflectivity
            * Heliostat.clean_factor
            * self.cosine_eff
            * (1.0 - self.f_shadow_block)
            * (1.0 - self.f_att)
            * (1.0 - self.f_spill)
        )
        if verbose:
            print("    TOTAL EFFICIENCY: {:04.1f} %".format(total_efficiency*100))
            print()
        
        # calculate total solar thermal energy delivered to the receiver [kW_th]
        self.total_contribution = NOM_Ib * Heliostat.area * total_efficiency

        # print result
        if verbose:
            print("--> total_contribution: {:,.0f} [W]".format(
                self.total_contribution))
            print("                        {:.1f}  [kW]".format(
                self.total_contribution/1000))
            print()
        
        return self.total_contribution
            

#===============================================================================
def place_row_of_heliostats(radius, d_theta, theta_0=0):
    """
    Places a single circle of heliostats, evenly distributed around the circle

    :param: radius  [float] Radius to plot heliostats [m]
    :param: d_theta [float] Angle between adjacent heliostats [rad]
    :param: theta   [float] Starting angle [rad]

    return: heliostat_row [list]  List of heliostat objects
            theta_0       [float] Angle between adjacent heliostats [rad]
            d_theta       [float] Starting angle [rad]
    """
    global net_power_thermal
    global solar_field
    global rows

    # number of heliostats
    number = int((2*pi) / d_theta)

    # message to inform
    print()
    print("Placing row: {:>2}, Radius: {:.1f} m, {} heliostats".format(
        rows, radius, number))

    # generate list of coordinates, and instantiate heliostats
    theta = theta_0
    heliostat_row = []
    rows += 1
    row_contribution = 0
    while theta < 2*pi:
        # generate coordinate
        x_coord = radius * cos(theta)
        y_coord = radius * sin(theta)
        z_coord = 0.0
        coordinate = np.array((x_coord, y_coord, z_coord))

        # instantiate heliostat, and add its contribution
        h = Heliostat(coordinate)
        heliostat_row.append(h)
        solar_field.append(h)
        theta += d_theta
        net_power_thermal += h.total_contribution
        row_contribution  += h.total_contribution

    # message to inform contribution of row
    net_thermal_MW = net_power_thermal / 1000000
    row_contribution_MW = row_contribution / 1000000
    percent_of_goal = (net_power_thermal / NET_POWER_TH_W) * 100
    print("--> Contribution of row: {:.1f} [MW]".format(row_contribution_MW))
    print("--> NET POWER THERMAL: {:.1f} [MW]".format(net_thermal_MW), end=", ")
    print(" ({:.2f} %)".format(percent_of_goal))
            
    return heliostat_row, theta_0, d_theta


def place_layers_of_heliostats(r_min=INNER_RADIUS, row_margin=None,
    margin_min=None, margin_max=None, oversize=1.0):
    """
    Places successive rows of heliostats, with a defined row_margin of spacing
    between each row.

    :param: r_min      [float] Inner heliostat radius from tower (to start) [m]
    :param: row_margin [float] Space between rows [m]
    :param: margin_min [float] Min distance between heliostats in same row [m]
    :param: margin_max [float] Max distance between heliostats in same row [m]
    :param: oversize   [float] Factor to oversize field by, before trimming [-]

    return heliostat_rows_list [list] List of rows
    """
    global net_power_thermal

    # stagger mirrors in each row
    theta_0 = 0
    d_theta = 0

    # default values
    if row_margin is None:
        row_margin = Heliostat.depth
    if margin_min is None:
        margin_min = Heliostat.width * 0.5
    if margin_max is None:
        margin_max = Heliostat.width * 3

    # oversize, so that least efficient heliostats can be removed afterwards
    OVERSIZE_POWER_TH_W = NET_POWER_TH_W * oversize
    
    # calculate starting parameters
    r = r_min
    margin = margin_min
    d_s = Heliostat.width + margin
    d_theta = arcsin(d_s/r)
    number = int((2*pi) / d_theta)
    d_theta = (2*pi) / number 

    heliostat_rows_list = []
    while net_power_thermal < OVERSIZE_POWER_TH_W:
        # update the margin, to see if it's getting to high
        margin = r * sin(d_theta)
        if margin > margin_max:
            # reset margin
            margin = margin_min
            d_s = Heliostat.width + margin
            d_theta = arcsin(d_s/r)
            # make d_theta so it is equal all the way around
            number = np.floor((2*pi) / d_theta)
            d_theta = (2*pi) / number 
            # increase radius (a little) for next section, and stagger mirrors
            r += row_margin / 2
            theta_0 = d_theta / 2 if theta_0 == 0 else 0

        # place heliostat row
        heliostat_row, theta_0, d_theta = \
            place_row_of_heliostats(r, d_theta, theta_0=theta_0)
        # stagger mirrors in the next row
        theta_0 = d_theta / 2 if theta_0 == 0 else 0

        # add row, and increase radius for the next row
        heliostat_rows_list.append(heliostat_row)
        r += row_margin

    return heliostat_rows_list


#===============================================================================
# Generate list of heliostats in the solar field.
# Place concentric circles of heliostats, staggered in each successive row to
# minimize blocking.
#-------------------------------------------------------------------------------
OVERSIZE_FACTOR = 1.5
print("Placing all heliostats...")
place_layers_of_heliostats(oversize=OVERSIZE_FACTOR)

print("--------> Done!")
print("          Placed {} heliostats".format(len(solar_field)))
print()


#-------------------------------------------------------------------------------
# trim away the least effective heliostats
if OVERSIZE_FACTOR > 1:
    print(".........................................")
    print("Oversize factor = {}".format(OVERSIZE_FACTOR))
    print("...Trimming least effective heliostats...")
    print("... ", end="")

    # sort list of heliostats (most to least contribution)
    solar_field.sort(key=lambda x: x.total_contribution, reverse=True)

    # remove heliostats from list, and subtract their thermal energy contribution
    while net_power_thermal > NET_POWER_TH_W:
        h = solar_field[-1]
        if net_power_thermal - h.total_contribution > NET_POWER_TH_W:
            solar_field.remove(h)
            net_power_thermal -= h.total_contribution
        else:
            print("--------> Done!")
            print(".........................................")
            break


#-------------------------------------------------------------------------------
# plot final list of heliostats
for h in solar_field:
    ppx, ppy, ppz = zip(list(h.position))
    plt.scatter(ppx, ppy, color='c', marker='s', s=5)


#===============================================================================
#-------------------------------- Final Values ---------------------------------
# calculate final values
net_thermal_KW = net_power_thermal / 1000
net_thermal_MW = net_power_thermal / 1000000
percent_of_goal = (net_power_thermal / NET_POWER_TH_W) * 100
net_power_electrical = ( net_power_thermal
                       * EFF_PB
                       * EFF_RECEIVER_THERMAL
                       * EFF_GENERATOR )
net_electrical_MW = net_power_electrical / 1000000

number_heliostats = len(solar_field)
total_mirror_area = MIRROR_AREA * number_heliostats 


#---------------------------------- Storage ------------------------------------
# TODO: Include 8hr storage
storage_capacity_MWH = ( NET_POWER_MW
                       * EFF_PB
                       * EFF_GENERATOR
                       * STORAGE_HOURS )
storage_capacity_KWH = storage_capacity_MWH * 1000


#----------------------------------- Costs -------------------------------------
# calculate costs
total_heliostat_cost = COST_HELIOSTAT * total_mirror_area
total_tower_cost     = COST_TOWER * TOWER_HEIGHT_TOTAL
total_receiver_cost  = COST_RECEIVER * net_thermal_KW
total_storage_cost   = COST_STORAGE * storage_capacity_KWH

final_cost = ( total_heliostat_cost
             + total_tower_cost
             + total_receiver_cost
             + total_storage_cost )


#--------------------------------- Land Area -----------------------------------
# get major and minor axes (assume ellipse)
coordinates = [h.position for h in solar_field]
min_x = min(map(lambda c: c[0], coordinates))
max_x = max(map(lambda c: c[0], coordinates))
min_y = min(map(lambda c: c[1], coordinates))
max_y = max(map(lambda c: c[1], coordinates))

# calculate area of ellipse
x_radius = (max_x - min_x) / 2
y_radius = (max_x - min_x) / 2
elliptical_area_m = pi * x_radius * y_radius

# convert to square kilometers
land_area_km2 = elliptical_area_m / 1000000


#===============================================================================
#-------------------------------------------------------------------------------
# (final output)
print()
print("====================================================")
print("                      SUMMARY                       ")
print("====================================================")
print("POWER OUTPUT:")
print("    NET POWER Thermal:    {:.1f} [MW_th]".format(net_thermal_MW))
print("    NET POWER Electrical: {:.1f} [MW_e]".format(net_electrical_MW), end=", ")
print(" ({:.2f} %)".format(percent_of_goal))
print()
print("LAND AREA:")
print("    (min_x: {: .2f} m)".format(min_x ))
print("    (max_x: {: .2f} m)".format(max_x ))
print("    (min_y: {: .2f} m)".format(min_y ))
print("    (max_y: {: .2f} m)".format(max_y ))
print("    (--> x_span:   {:.2f} m)".format(x_radius*2))
print("    (--> y_span:   {:.2f} m)".format(y_radius*2))
print("    Total Land Area:       {:,.2f} [km2]".format(land_area_km2))
print()
print("HELIOSTATS:")
print("    Total Number:          {:,d} heliostats".format(number_heliostats))
print("    Total Mirror Area:     {:,.1f} [m2]".format(total_mirror_area))
print()
print("TOWER:")                   
print("    Tower Height           {:,.1f} [m]".format(TOWER_HEIGHT_TOTAL))
print("    Receiver Height (mid)  {:,.1f} [m]".format(TOWER_HEIGHT_OPTICAL))
print()
print("STORAGE:")                 
print("    Total Storage:         {:,.1f} [MWh]".format(storage_capacity_MWH))
print()
print("INVESTMENT COSTS:")
print("    Total Tower Cost:     $ {:,.0f}".format(total_tower_cost))
print("    Total Receiver Cost:  $ {:,.0f}".format(total_receiver_cost))
print("    Total Helistat Cost:  $ {:,.0f}".format(total_heliostat_cost))
print("    Total Storage Cost:   $ {:,.0f}".format(total_storage_cost))
print("----------------------------------------------------")
print("            FINAL COST:   $ {:,.0f}".format(final_cost))
print("----------------------------------------------------")
print()


#-------------------------------------------------------------------------------
# show plot
plt.show(block=False)
# (wait for user to press enter to close plot and exit)
input()
