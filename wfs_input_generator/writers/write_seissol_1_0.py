#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input file writer for Seissol.

:copyright:
    Marek Simon (marek.simon@geophysik.uni-muenchen.de), 2013
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)


Requires Metis to be installed (see class Gambit2Metis)
"""

EARTH_RADIUS = 6371 * 1000
deviation_from_sphere = 0.0  # EARTH_RADIUS/200

import numpy as np
import scipy as sp
import obspy
import os
import timeit
from scipy.spatial import Delaunay
import math
import json

# Define the required configuration items. The key is always the name of the
# configuration item and the value is a tuple. The first item in the tuple is
# the function or type that it will be converted to and the second is the
# documentation.

REQUIRED_CONFIGURATION = {
    'model': (str, "e.g.: model_eucrust_small_new.dat,\
        Set constant, PREM, linear, or filename of\
        Material file. See Seissol Docu 4.2.1."),
    'mesh': (str, "e.g.: eucrust_small_new. Default: Gambit neutral\
        file (*.neu)"),
    'max_time': (int, "maximum simulated time. can also be restricted\
        by max_iteration, max_wallclocktime"),
    'working_directory': (str, "working directory"),
    'number_of_processors': (int, 'Number of Processors. Will Partition the mesh accordingly')
}

DEFAULT_CONFIGURATION = {
    #<----------------------------------->
    #<  S E I S S O L - V E R S I O N    >
    #<----------------------------------->
    'version': (18, int, "Version number. See Seissol Docu 4.1"),
    #<----------------------------------->
    #<  E Q U A T I O N S                >
    #<----------------------------------->
    'dimension': (3, int, "Number of Dimensions. Either 2D or 3D"),
    'advection': (0, int, "Advection. 1 or 0 for On or Off. Only \
    available for 2D. See Seissol Docu 4.2"),
    # Units of Advection velocity??
    'advection_velocity': ((1.0, 1.0, 1.0), tuple, "If Advection \
    is switched on, give Advection velocity values as tuple \
    (v_1, v_2, v_3). See Seissol Docu 4.2"),
    'anisotropy': (0, int, "Anisotropy. 1 or 0 for On or Off. \
    See Seissol Docu 4.2"),
    'anelasticity': (0, int, "Anelasticity. 1 or 0 for On or Off. \
    Can Increase Computational cost dramatically if switched on! \
    See Seissol Docu 4.2"),
    # Will change output!
    'poroelasticity': (0, int, "Poroelasticity. 1 or 0 for On or Off.\
     Can Increase Computational cost dramatically if switched on! \
     See Seissol Docu 4.2"),
    # What does this actually mean? Documentation??
    # only implemented for 2d, not yet for 3d.
    'adjoint': (0, int, "Adjoint Calculation. 1 or 0 for On or Off."),
    # What is the unit, meaning of material reference values?? kg/m^3. rho
    # mü=Pascal lambda=Pascal
    'material_reference_values': ((3600, 9.0e10, 1.11e11), tuple,
    "Material reference values, in case model=0 or ?? as tupel \
    (rho, mü, lambda). See Seissol Docu 4.2.1"),
    'randomfield': (0, int, "Randomfield. 1 or 0 for On or Off. \
    See Seissol Docu 4.2.2."),
    #<----------------------------------->
    #<  INITIAL CONDITION                >
    #<----------------------------------->
    # Documentation missing!!!
    'initial_values': ('Gauss_Puls_Rad', str, "Documentation missing"),
    'ic_parameter': (1, int, "Documentation missing"),
    'center_of_gaussspulse': ((0.0, 0.0, 0.0), tuple, "Center of \
    Gauss-pulse as tuple (X,Y,Z). Documentation missing"),
    'amplitudes': ((0.0, 0.0, 0.0), tuple, "Amplitudes as tuple \
    (A_1, A_2, A_3). Documentation missing"),
    'halfwidths': ((5.0e3,  5.0e3,  5.0e3), tuple, "Halfwidths as \
    tuple (h_x, h_y, h_z). Documentation missing"),
    #<----------------------------------->
    #<  BOUNDARIES                       >
    #<----------------------------------->
    # Documentation
    'number_free_surfaces': (1, int, "101 in Gambit Neutral File"),
    'number_non_conforming_boundaries': (0, int, "102 in Gambit Neutral File"),
    'number_inflow_boundaries': (0, int, "104 in Gambit Neutral File"),
    'number_open_boundaries': (1, int, "105 in Gambit Neutral File"),
    'number_periodic_boundaries': (0, int, "106 in Gambit Neutral File"),
    'number_dynamic_rupture_fault_boundaries': (0, int, "103 in Gamit Neutral File"),
    #<----------------------------------->
    #<  S O U R C E    T E R M S         >
    #<----------------------------------->
    #
    'sourcetype': (50, int, "Defines source type. See Seissol \
    Docu 4.5."),
    'source_file': ('source.dat', str, "Name of the Sourcefile, \
    if sourcetype=50. Do not change if you want to use automatically\
     generated source!"),
    'source_time_function': ('gauss', str, "Choose the shape of the source\
     time function. Either 'gauss' or 'dirac'"),
    #<----------------------------------->
    #<  S P O N G E    L A Y E R         >
    #<----------------------------------->
    'sponge': (0, int, "Not implemented in Inputfile generator at \
     the moment. Will not generate Input-files for anything other \
     then 0. See Seissol Docu 4.6."),
    #<----------------------------------->
    #<  M E S H                          >
    #<----------------------------------->
    'meshgenerator': ("Gambit3D-Tetra", str,  "For the moment please\
     only use Gambit *.neu files. For Documentation see \
     SeiSsol_tutorial.ppt"),
    #<----------------------------------->
    #<  D I S C R E T I S A T I O N      >
    #<----------------------------------->
    'fine_output': (0, int, "For Documentation see Seissol Docu 4.8"),
    'restartfile': (0, int, "Simulation is restated from previous \
    simulation. For Documentation see Seissol Docu 4.1"),
    'DGMethod': (1, int, "1 = GLobal Time Stepping, 3 = Local Time \
    Stepping. 4 =  See Seissol Docu 4.8"),
    'CK': (0, int, "0 = Standard ADER Cauchy Kowaleski, 1 = Space- \
    Time Discontinuous ADER time integration. See Seissol Docu 4.8"),
    'fluxmethod': (0, int, "0 = godunov, 1 = rusanov"),
    'DGCycle': (1, int, "# 1: one iteration = one cycle; 2: one \
    iteration = update of all elements"),
    # Documentation missing. Old structure?
    'basisfunction_degree': (0, int, ""),
    'reconstructed_basisfunction_degree': (0, int, ""),
    'stencil_security_factor': (0, int, ""),
    'reconstruction_type': (0, int, "0=linear reconstruction or ELSE\
     (nonlinear central reconstruction)"),
    'exponent_r': (0, int, ""),
    'coefficient_epsilon': (0, int, ""),
    'linear_weight': (0, int, ""), 'limiter_security_factor':
    (0, int, ""),
    'minspace_order': (3, int, "Min order basis functions. For 3D \
    from 1-7.\
    For 2D from 1-10. See Seissol Docu 4.8"),
    'maxspace_order': (3, int, "Max order basis functions. For 3D \
    from 1-7.\
    For 2D from 1-10. See Seissol Docu 4.8"),
    'zone_order_flag': (0, int, "Sort order basis \
    function by cell size. See Seissol Docu 4.8"),
    'pAdaptivity_file_name': ('pAdaptivity_file_name', str,
    "If  minspace_order<maxspace_order, required filename! See \
    Seissol Docu 4.8"),
    # Documentation missing. Old structure?
    'material_basis_function_order': (1, int, ""),
    'courant_number': (0.5, float, "Courant number. See Seissol Docu\
    4.8"),
    'min_time_step': (10000, int, "Leave large for automatic pick. \
    See Seissol Docu 4.8"),
    #<----------------------------------->
    #<  O U T P U T                      >
    #<----------------------------------->
    'rotational_output': (0, int, "Flag for rotational output. 0=no,\
    1=rotational, 2=moment tensor contribution. See Seissol Docu \
    4.9"),
    'rotation_components': ((1, 1, 1), tuple, "Tupel (r_1,r_2,r_3) \
    for rotational rate components. See Seissol Docu 4.9"),
    'variable_output': ((0, 0, 0, 0, 0, 0, 1, 1, 1), tuple, "Tupel \
    of variable output. (o_1, ... o_9) with o_n either 0 or 1"),
    'material_parameters_output': ((1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                                   tuple, "Tupel \
    of material_parameter output. (o_1, ... o_9) with o_n either 0 \
    or 1. See Seissol Docu 4.9"),
    'output_character': (0, int, "0 = geometry based; 1= calculation\
    based; See Seissol Docu 4.9.1"),
    'output_format': (1, int, "Currently only 1 = Tecplot is \
    implemented"),
    'timestep_output': (500, int, "Print volume output every N timesteps. \
    See Seissol Docu 4.9.2"),
    'time_output': (500.0, float, "Print volume output every T seconds. "),
    'output_index': (1, int, "Switch: 1 = timestep_output, 2 = \
    time_output, 3 = timesteps+time"),
    'sampling': (0.1, float, "Pickpoint sampling-rate"),
    #<----------------------------------->
    #<  A B O R T - C R I T E R I A      >
    #<----------------------------------->
    'max_time': (2000, int, "Set maximum time in seconds that you \
    want to simulate"),
    'max_iteration': (100000000, int, "Set maximum number of \
    iterations that you would want. Maxiumum = 9digits"),
    'max_wallclocktime': (1e20, str, "Set maximimum wallclock time"),
    # Ask what the delay means!
    'delay': (0, int, "Delay"),
    'check_mesh': (True, bool, "Check if stations are in the mesh or not. \
    If they are not, move them incrementally towards the coordinate system \
    origin until they are in the mesh. If they are much lower then the \
    lowest point of the free surface boundary, but still not in the mesh, \
    they will be considered outside of the mesh and dropped")}

    # The default configuration item. Contains everything that can sensibly
    # be set to some default value. The syntax is very similar to the
    # REQUIRED_CONFIGURATION except that the tuple now has three items, the
    # first one being the actual default value.

def _check_configuration(config):
    
    """
    first draft of configuration checker. can still be extended. \
    checks if options are compatible.
    """
    
    if config.advection is 1:
        config.advection_velocity = config.advection_velocity + '\n'
    else:
        config.advection_velocity = ''
    if config.dimension is 2:
        config.anisotropy = 0
    if config.dimension is 3:
        config.anelasticity = 0
    # self.poroelasticity=poroelasticity #0=no poroelasticity, 1=Biot's theory
    # to describe wave propagation through porous materials filled with a
    # fluid. For the high-frequency range, option 1 should be chosen which
    # uses the standard ADERDG
    # For the low-frequency case, option 2 uses the ADER-DG(ST)
    # and option 3 the ADER-DG(FS) method.
    if config.model in 'constant':
        config.material = 0
    elif config.model in 'PREM':
        config.material = 2  # 0=constant, 1=from file, 2=PREM, 3=linear
    elif config.model in 'linear':
        config.material = 3
    else:
        config.material = 1
        if os.path.exists(os.path.join(config.working_directory,
                                       config.model)):
            print 'Use Material file as Model: ', config.model
            config.materialfile_name = config.model
        else:
            raise IOError, ('Cannot find Material file: %s' % config.model)
    if config.anisotropy is 1 and config.material is not 1:
        raise IOError, 'To include anisotropy choose material file'
    if config.poroelasticity is not 0 and config.material is not 1:
        raise IOError, 'To include poroelasticity choose material file'
    try:
        config.material_reference_values = "%-10.3e %-10.3e %-10.3e"\
            % config.material_reference_values
    except ValueError:
        print 'Cannot parse material reference values'
        print "default material reference values!!"
        print "mv_1=3600."
        print "mv_2=9.0e10"
        print "mv_3=1.116e11"
        mv_1 = 3600.
        mv_2 = 9.0e10
        mv_3 = 1.116e11
        config.material_reference_values = "%-10.3e %-10.3e %-10.3e" \
            % (mv_1, mv_2, mv_3)
    if config.material is 1:
        if len(config.model) < 32:
            config.materialfile_name = "%-37s \n" % config.model
        else:
            raise Error, "material file name too long!! <32char!!"
    else:
        config.materialfile_name = ''
    if config.source_time_function not in ['gauss', 'dirac']:
        raise Exception, 'source_time_function must be either gauss or dirac'
    if config.poroelasticity is 3:
        config.CK = 0
    elif config.poroelasticity is 2:
        config.CK = 1
    if config.DGMethod is 3:
        config.DGcycle = "%-37s \n" % config.DGcycle
    else:
        config.DGcycle = ''
    if config.DGMethod is 4:
        config.basisfunction_degree = "%-37s \n" % \
            config.basisfunction_degree
        config.reconstructed_basisfunction_degree = "%-37s \n" % \
            config.reconstructed_basisfunction_degree
        config.stencil_security_factor = "%-37s \n" % \
            config.stencil_security_factor
        config.reconstruction_type = "%-37s \n" % \
            config.reconstruction_type  # 0 (linear reconstruction) or ELSE (nonlinear central reconstruction)
        config.exponent_r = "%-37s \n" % config.exponent_r
        config.coefficient_epsilon = "%-37s \n" % config.coefficient_epsilon
        config.linear_weight = "%-37s \n" % config.linear_weight
        config.limiter_security_factor = "%-37s \n" % \
            config.limiter_security_factor
    else:
        config.basisfunction_degree = ''
        config.reconstructed_basisfunction_degree = ''
        config.stencil_security_factor = ''
        config.reconstruction_type = ''
        config.exponent_r = ''
        config.coefficient_epsilon = ''
        config.linear_weight = ''
        config.limiter_security_factor = ''
    if config.minspace_order is not config.maxspace_order:
        config.zone_order_flag = '%s \n' % config.zone_order_flag
        if str(1) in config.zone_order_flag:
            try:
                config.pAdaptivity_file_name = "%-37s \n" % \
                    config.pAdaptivity_file_name
            except:
                IOError, 'Zone order flag = 1, but \
                pAdaptivity_filename: %s does not work' % \
                    config.pAdaptivity_file_name
        else:
            config.pAdaptivity_file_name = ''
    else:
        config.zone_order_flag = ''
        config.pAdaptivity_file_name = ''
    if config.rotational_output is 1:
        config.rotation_components = ("%s %s %s \n") % \
            config.rotation_components
        # rotation_components=(1, 1, 1) #AXES X, Y, Z
    else:
        config.rotation_components = ''
    config.variable_output = ("%s %s %s %s %s %s %s %s %s") % \
        config.variable_output
    # variable_output=(0, 0, 0, 0, 0, 0, 1, 1, 1) # up to 9
    config.material_parameters_output = ("%s %s %s %s %s %s %s %s %s %s %s") \
        % config.material_parameters_output

    # material_parameters_output=(1, 1, 1) # up to 11
    # PARAMS['number_of_pickpoints']=len(self.station_coordinatesXYZ.split('\n'))-1
    # Abort Criteria


def _prepare_station_coordinates(config, stations):
    
    """
    Is called from within _write_input()
    converts lat/lon/depth to XYZ with origin at the center of the earth
    """

    station_coordinatesXYZ = str()
    coordinate_format = ('   %1.7e   %1.7e   %1.7e   \n')
    for station in stations:
        lat = station['latitude'] * np.pi / 180
        lon = station['longitude'] * np.pi / 180
        dep = - station['elevation_in_m'] * 1000
        r = 6371000 - dep
        x = r * sp.cos(lat) * sp.cos(lon)
        y = r * sp.cos(lat) * sp.sin(lon)
        z = r * sp.sin(lat)
        station_coordinatesXYZ += coordinate_format % (x, y, z)

    if config.check_mesh == True:
        # Check if station coordinates in Mesh
        csm = CheckStationsMesh(filename_mesh=os.path.join(config.
                                                       working_directory, config.mesh + '.neu'), station_coordinates=
                            station_coordinatesXYZ)
        points_in_mesh = csm.move_points_to_mesh()
        station_coordinatesXYZ = str()
        used_stations = []
        # add only those stations that could successfully be placed within the mesh
        for point in points_in_mesh:
            station_coordinatesXYZ += coordinate_format % (point[0], point[1],
                                                       point[2])
            used_stations.append(_xyz_to_station(point, stations))
    else:
        used_stations=stations

    config.used_stations = _write_used_stations(used_stations)

    return station_coordinatesXYZ


def _xyz_to_station(xyz, stations):
    """
    checks if point xyz can be matched by latitude and longitude to a station in stations

    :type xyz: class: 'numpy.ndarray'
    :param: X, Y, Z coordinate in a 3X1 dim ndarray
    :type stations: class: 'list(dict())'
    :param: stations as defined as a list of dictionaries with the shape of
    [{"id": "BW.FURT",
    "latitude": 48.162899,
    "longitude": 11.2752,
    "elevation_in_m": 565.0,
    "local_depth_in_m": 0.0},..]

    :rtype station: class: 'dict'
    :param returns dictionary with station information of matching station

    to be used within _prepare station coordinates function
    """

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    r = np.sqrt(x * x + y * y + z * z)

    lat = np.arcsin(z / r) / np.pi * 180
    lon = np.arctan(y / x) / np.pi * 180

    for station in stations:
        if math.ceil(float(station['latitude']) * 100) / 100 == \
            math.ceil(lat * 100) / 100 and \
            math.ceil(float(station['longitude']) * 100) / 100 ==\
                math.ceil(lon * 100) / 100:
            return station


def _write_used_stations(used_stations):
    """
    write out the stations used in json format
    """

    try:
        used_stations_json = json.dumps(
            used_stations, sort_keys=True, indent=4,
            separators=(',', ': '))

        return used_stations_json
    except:
        raise Exception, 'Cannot dump used stations into JSON'


def _check_mesh(config):
    """
    Is called from within _write_input()
    checks only the number of open surfaces, and assumes 1 free surface
    works only for Gambit.neu files
    """

    if 'Gambit' in config.meshgenerator:
        f = open(os.path.join(config.working_directory, config.mesh + '.neu'))
    else:
        print 'cant open meshfile'
        return 1, 0
    for i in range(6):
        junk = f.readline()
    tmp = f.readline()
    f.close()
    tmp = [int(s) for s in tmp.split() if s.isdigit()]
    nbsets = tmp[3]
    nbfreesurface = 1
    nbopensurface = nbsets - nbfreesurface
    config.number_free_surfaces = nbfreesurface
    config.number_open_surface = nbopensurface


def _partition_mesh(config):
    """
    Partitions the mesh to the number of cores
    """

    if 'Gambit' in config.meshgenerator:
        Gambit2Metis(filename=os.path.join(config.working_directory,
                                           config.mesh + '.neu'), number_proc=config.number_of_processors)


def _write_source_dat(config, event):
    """
    writes a gaussian as a source time function. requires at the moment \
    to have a half-duration (e.g. from global cmt). The shape of the \
    gauss function for now is somewhat arbitrary. feel free to change \
    according to your demands!
    Rotates the moment tensor to XYZ
    Be carefull, for very shallow events, event can be out of mesh, if \
    the mesh is too coarse on the surface.
    """

    # event
    lat = event['latitude'] * np.pi / 180
    lon = event['longitude'] * np.pi / 180
    dep = event['depth_in_km'] * 1000
    r = EARTH_RADIUS - abs(dep)
    x = r * sp.cos(lat) * sp.cos(lon)
    y = r * sp.cos(lat) * sp.sin(lon)
    z = r * sp.sin(lat)

    if 'gauss' in config.source_time_function:
        # source
        # STF = gauss )
        hd = 10.0
        # hd=float(event['halfduration'])
        # ???
        sd = 6.0 * hd
        ssampling = config.sampling / 100
        src = np.arange(0, sd + ssampling, ssampling)

        def gauss(x):
            x = 1.0 / hd / \
                np.sqrt(2 * np.pi) * np.exp(
                    -1 / 2 * (x - sd / 2) / hd * (x - sd / 2) / hd)
            return x
        vecgauss = np.vectorize(gauss, otypes=[float])
        src = vecgauss(src)
        src_nsamp = len(src)
        src_dt = ssampling

    else:
        src = [0.0, 1.0, 0.0]
        src_nsamp = len(src)
        src_dt = ssampling

    # Coordinate transformation matrix R_ij:
    R_11 = sp.cos(lat) * sp.cos(lon)
    R_12 = sp.cos(lat) * sp.sin(lon)
    R_13 = sp.sin(lat)
    R_21 = sp.sin(lat) * sp.cos(lon)
    R_22 = sp.sin(lat) * sp.sin(lon)
    R_23 = -sp.cos(lat)
    R_31 = -sp.sin(lon)
    R_32 = sp.cos(lon)
    R_33 = 0.0
    R = np.matrix([[R_11, R_12, R_13],
                   [R_21, R_22, R_23],
                   [R_31, R_32, R_33]])
    M = np.matrix([[event['m_rr'], event['m_rt'], event['m_rp']],
                   [event['m_rt'], event['m_tt'], event['m_tp']],
                   [event['m_rp'], event['m_tp'], event['m_pp']]])
    # Coordinate Transformation of Momenttensor
    M_rot = R.getI() * M * R
    # SAVE IN SEISSOL FORMAT
    source_dat = str()
    source_dat += 'Seismic Moment Tensor\n'
    for row in M_rot:
        source_dat += '%2.9e %2.9e %2.9e \n' % (row.tolist()[0][0],
                                                row.tolist()[0][1], row.tolist()[0][2])
    source_dat += 'Number of subfaults\n'
    source_dat += '1\n'
    source_dat += '               x                y                z\
    strike        dip             rake            area       Onset time\n'
    source_dat += '%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\
     \n' % (x,     y,     z,   0.0,   0.0,   0.0,   1.0,   0.0)
    source_dat += 'source time function\n'
    source_dat += '%16.9e %16.9e \n' % (src_dt, src_nsamp)
    source_dat += 'samples \n'
    for value in src:
        source_dat += '%16.9e\n' % value
    return source_dat


def _write_input(config, stations):
    """
    writes input.par file
    """

    # check mesh and stations
    _check_mesh(config)
    station_coordinatesXYZ = _prepare_station_coordinates(config, stations)

    input_par = ''.join([
                        "<-----------------------------------> \n",
                        "<  S E I S S O L - V E R S I O N    > \n",
                        "<-----------------------------------> \n",
                        "%-37s! seissol-version \n" % config.version,
                        "\n",
                        "<-----------------------------------> \n",
                        "<  E Q U A T I O N S                > \n",
                        "<-----------------------------------> \n",
                        "%-37s! number of dimensions ( 2 or 3 ) \n" % config.dimension,
                        "%-37s! Advection \n" % config.advection,
                        str(config.advection_velocity),
                        "%-37s! Anisotropy (only in 3D) \n" % config.anisotropy,
                        "%-37s! Anelasticity (only in 2D) \n" % config.anelasticity,
                        "%-37s! Poroelasticity \n" % config.poroelasticity,
                        "%-37s! adjoint \n" % config.adjoint,
                        "%-37s! Material (0=constant;1=from file;2=PREM;3=linear) \n" %
                        config.material,
                        "%-37s! Reference values \n" % config.material_reference_values,
                        config.materialfile_name,
                        "%-37s! random field \n" % config.randomfield,
                        " \n",
                        "<-----------------------------------> \n",
                        "<  INITIAL CONDITION                > \n",
                        "<-----------------------------------> \n",
                        "%-37s! problem-dependent initial values \n" % config.initial_values,
                        "%-37s! IC Parameter \n" % config.ic_parameter,
                        "%-37s! Center of Gaussspulse (x,y,z) \n" % str("%1.1f  %1.1f  %1.1f" % config.center_of_gaussspulse),
                        "%-37s! Amplitudes \n" % str("%1.1f  %1.1f  %1.1f" % config.amplitudes),
                        "%-37s! Halfwidths (hx>0,hy>0,hz>0) \n" % str("%1.3e  %1.3e  %1.3e" % config.halfwidths),
                        " \n",

                        "<-----------------------------------> \n",
                        "<  BOUNDARIES                       > \n",
                        "<-----------------------------------> \n",

                        "%-37s! number of free surfaces \n" % config.number_free_surfaces,
                        "%-37s! number of non_conforming boundaries \n" %
                        config.number_non_conforming_boundaries,
                        "%-37s! number of dynamic rupture fault boundaries \n" %
                        config.number_dynamic_rupture_fault_boundaries,
                        "%-37s! number of inflow boundaries \n" %
                        config.number_inflow_boundaries,
                        "%-37s! number of open (outflow) boundaries \n" %
                        config.number_open_boundaries,
                        "%-37s! periodic boundaries \n" % config.number_periodic_boundaries,
                        " \n",
                        "<-----------------------------------> \n",
                        "<  S O U R C E    T E R M S         > \n",
                        "<-----------------------------------> \n",
                        "%-37s! Rupture plane source \n" % config.sourcetype,
                        "%-37s! rupture plane input file \n" % config.source_file,
                        " \n",
                        "<-----------------------------------> \n",
                        "<  S P O N G E    L A Y E R         > \n",
                        "<-----------------------------------> \n",
                        "%-37s! \n" % config.sponge,
                        " \n",
                        "<-----------------------------------> \n",
                        "<  M E S H                          > \n",
                        "<-----------------------------------> \n",
                        "%-37s! name of mesh file \n" % config.mesh,
                        "%-37s! name of meshgenerator (format) \n" % config.meshgenerator,
                        "0                                    ! Periodic boundaries \n",
                        "0.        0.       0.                ! Mesh displacement   \
    (2D or 3D) \n",
                        "1.        0.       0.                ! Mesh scaling matrix \
    (2D or 3D) \n",
                        "0.        1.       0.                ! Mesh scaling matrix \
    (2D or 3D) \n",
                        "0.        0.       1.                ! Mesh scaling matrix \
    (2D or 3D) \n",
                        " \n",
                        "<-----------------------------------> \n",
                        "<  D I S C R E T I S A T I O N      > \n",
                        "<-----------------------------------> \n",
                        "%-37s! fine output \n" % config.fine_output,
                        "%-37s! Restart file \n" % config.restartfile,
                        "%-37s! DGmethod \n" % config.DGMethod,
                        "%-37s! (0=CK; 1=STDG) \n" % config.CK,
                        "%-37s! flux 0-godunov 1-rusanov \n" % config.fluxmethod,
                        config.DGcycle,
                        config.basisfunction_degree,
                        config.reconstructed_basisfunction_degree,
                        config.stencil_security_factor,
                        config.reconstruction_type,
                        config.exponent_r,
                        config.coefficient_epsilon,
                        config.linear_weight,
                        config.limiter_security_factor,
                        "%-37s! Min Space accuracy \n" % config.minspace_order,
                        "%-37s! Max Space accuracy \n" % config.maxspace_order,
                        config.zone_order_flag,
                        config.pAdaptivity_file_name,
                        "%-37s! Material basis function order \n" %
                        config.material_basis_function_order,
                        "%-37s! courant number (<=1.0) \n" % config.courant_number,
                        "%-37s! Manualy chosen minimum time \n" % config.min_time_step,
                        ' \n',
                        "<-----------------------------------> \n",
                        "<  O U T P U T                      > \n",
                        "<-----------------------------------> \n",
                        'out', '\n',
                        "%-37s! rotational output ( 0=no / 1=yes ) \n" %
                        config.rotational_output,
                        config.rotation_components,
                        "%-37s! variables \n" % config.variable_output,
                        "%-37s! material values (rho / mu / lambda) \n" %
                        config.material_parameters_output,
                        "%-37s! character of output (0=geo.based, 1=calc. based) \n" %
                        config.output_character,
                        "%-37s! format of output (0=IDL, 1=TECPLOT, 2=IBM DX, 4=GiD)) \n" %
                        config.output_format,
                        "%-37s! index of printed info at timesteps \n" %
                        config.timestep_output,
                        "%-37s! index of printed info at time \n" % config.time_output,
                        "%-37s! Criterion for index of printed info: 1=timesteps,2=time,\
    3=timesteps+time \n" % config.output_index,
                        "%-37s! Pickpoint Sampling \n" % config.sampling,
                        "%-37s! # of coordinates to pick \n" %
                       (len(station_coordinatesXYZ.split('\n')) - 1),
                        station_coordinatesXYZ,
                        str(0) + ' \n',
                        ' \n',
                        "<-----------------------------------> \n",
                        "<  A B O R T - C R I T E R I A      > \n",
                        "<-----------------------------------> \n",
                        "%-37s! maximum computing time \n" % config.max_time,
                        "%-37s! number of maximum iteration \n" % config.max_iteration,
                        "%-37s! Max. wallclock time \n" % config.max_wallclocktime,
                        "%-37s! Delay \n" % config.delay,
                        "\n",
                        "<-----------------------------------> \n",
                        "<  A N A L Y S I S  O F   D A T A   > \n",
                        "<-----------------------------------> \n",
                        "0                                    ! Do not analyse data \n",
                        "\n",
                        "<-----------------------------------> \n",
                        "<  D E B U G G I N G     M O D U S  > \n",
                        "<-----------------------------------> \n",
                        "0                                     ! No debugging \n",
                        "\n"])

    return input_par


def write(config, events, stations):

    if len(events) > 1:
        raise TypeError, 'Can only simulate one event at a time!'
    event = events[0]
    output_files = {}
    _check_configuration(config)
    _partition_mesh(config)

    output_files['input.par'] = _write_input(config, stations)
    output_files['used_stations.json'] = config.used_stations
    output_files['source.dat'] = _write_source_dat(config, event)

    for key in output_files.iterkeys():
        output_files[key] += "\n"

    return output_files

# Mesh partitioning


class Gambit2Metis(object):

    def __init__(self, filename='vercehpc.neu', number_proc=40):
        if filename[-4:] not in '.neu':
            print 'filename not in *.neu format!!'
        filename = filename[:-4]
        proc = str(number_proc)
        # rewritten to python by marek simon
        print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print '     %%                                                             %%'
        print '     %%                        Gambit2Metis                         %%'
        print '     %%                                                             %%'
        print '     %%          (c) Martin Kaeser, Uni Trento, 22.01.2008          %%'
        print '     %%                                                             %%'
        print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print '     %%                                                             %%'
        print '     %% Gambit2Metis converts a GAMBIT-meshfile of                  %%'
        print '     %% tetrahedrons or hexahedrons      %%'
        print '     %% stored as "filename.neu"                                    %%'
        print '     %% into a METIS-file and calls METIS for mesh partitioning.    %%'
        print '     %% METIS mesh partitioner:                                     %%'
        print '     %% http://www-users.cs.umn.edu/~karypis/metis/metis/index.html %%'
        print '     %%                                                             %%'
        print '     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        print '\n \n'

        t = timeit.time.time()
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%        Read Gambit Data
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = open(filename + '.neu')
        print '\n'
        print '--------------------------------------------------------',\
            '---------------------------'
        print ' Reading data from: %s' % filename + '.neu'
        for i in range(6):
            junk = f.readline()
        tmp = f.readline()
        tmp = [int(s) for s in tmp.split() if s.isdigit()]
        NUMNP = tmp[0]  # %Number of Vertices
        NELEM = tmp[1]  # %Number of Elements
        NGRPS = tmp[2]  # %Number of Element Groups
        NBSETS = tmp[3]
            #                     %Number of Boundary Condition Sets
        NDFCD = tmp[4]  # %Number of Dimensions (2 or 3)
        try:
            if 3 == NDFCD:
                f.readline()
                f.readline()
                vertices = np.fromfile(file=f, count=NUMNP * 4, sep=" ")
                vertices = vertices.reshape(NUMNP, 4)
                vertices = vertices.transpose()
                f.readline()
                f.readline()
                position = f.tell()
                first_line = f.readline()
                first_line = [int(s)
                              for s in first_line.split() if s.isdigit()]
                if first_line[1] is 6:
                    n = 7
                if first_line[1] is 4:
                    n = 11
                f.seek(position, 0)
                elements = np.fromfile(file=f, count=NELEM * n, sep=" ")\
                    .reshape(NELEM, n).transpose()
                if elements[1][0] is 6:
                    print 'Read all tetrahedrons successfully!'
                if elements[1][0] is 4:
                    print 'Read all hexahedrons successfully!'

            elif 2 == NDFCD:
                f.readline()
                f.readline()
                vertices = np.fromfile(file=f, count=NUMNP * 3, sep=" ")
                vertices = vertices.reshape(NUMNP, 3)
                vertices = vertices.transpose()
                f.readline()
                f.readline()
                position = f.tell()
                first_line = f.readline()
                first_line = [int(s)
                              for s in first_line.split() if s.isdigit()]
                if first_line[1] is 3:
                    n = 6
                if first_line[1] is 2:
                    n = 7
                f.seek(position, 0)
                elements = np.fromfile(file=f, count=NELEM * n, sep=" ").\
                    reshape(NELEM, n).transpose()
                if elements[1][0] is 3:
                    print 'Read all triangles successfully!'
                if elements[1][0] is 2:
                    print 'Read all quadrilaterals successfully!'
        except:
            print '##ERROR: The GAMBIT input file shows a problem\n\n ##'
        #   return
        f.close()

        print '-----------------------------------------------------------------------------------'

        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%        Writing Metis Data
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        print 'Write METIS-file: %s' % filename + '.met'
        f = open(filename + '.met', 'w')
        if 3 == NDFCD:
            if int(elements[1][0]) is 6:
                f.write('%12i %5i \n' % (NELEM, 2))
                # np.savetxt(f,  elements[3:n].transpose())
                for [a, b, c, d] in elements[3:n].transpose():
                    f.write('%12i %12i %12i %12i \n' % (a, b, c, d))
            if int(elements[1][0]) is 4:
                f.write('%12i %5i \n' % (NELEM, 3))
                # np.savetxt(f,  elements[3:n].transpose())    <------- is
                # somehow slower ??
                for [a, b, c, d, e, ff, g, h] in elements[3:n].transpose():
                    f.write('%12i %12i %12i %12i %12i %12i %12i %12i\n' %
                           (a, b, c, d, e, ff, g, h))
        elif 2 == NDFCD:
            if int(elements[1][0]) is 3:
                f.write('%12i %5i \n' % (NELEM, 1))
                for [a, b, c] in elements[3:n].transpose():
                    f.write('%12i %12i %12i\n' % (a, b, c))
            if int(elements[1][0]) is 2:
                f.write('%12i %5i \n' % (NELEM, 4))
                for [a, b, c, d] in elements[3:n].transpose():
                    f.write('%12i %12i %12i %12i\n' % (a, b, c, d))
        f.close()
        t = t - timeit.time.time()
        print '-----------------------------------------------------------------------------------'
        print ' Conversion finished successfully!  (%s CPU sec)\n' % t
        print '-----------------------------------------------------------------------------------'
        print '\n'

        print '-----------------------------------------------------------------------------------'
        print 'METIS partition starts!\n'
        print '-----------------------------------------------------------------------------------'
        print '\n'
        os.system('mesh2dual ' + filename + '.met')
        os.system('pmetis ' + filename + '.met.dgraph ' + proc)

        # clean up of files
        os.system('rm ' + filename + '.met')
        os.system('rm ' + filename + '.met.dgraph')
        os.system('mv ' + filename + '.met.dgraph.part.' +
                  proc + '  ' + filename + '.met.epart.' + proc)

        print '-----------------------------------------------------------------------------------'
        print 'METIS partition done!'
        print '-----------------------------------------------------------------------------------'
        print '\n'
        print '-----------------------------------------------------------------------------------'
        print 'Wrote final METIS-file: %s' % filename + '.met.epart.' + proc
        print '-----------------------------------------------------------------------------------'
        print '\n'

# Check Stations


class CheckStationsMesh(object):

    def __init__(
        self, filename_mesh='/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/eucrust_small_new.neu',
            station_coordinates=str()):
        """
        Checks if stations given by station_coordinates are inside the mesh.
        If not, the function "move_points_to_mesh" will move the points
        into the mesh and return a numpy ndarray of the coordinates of
        all the stations that were successfully placed into the mesh

        :type filename: class:'~str'
        :param: Fullpath, or relative path of the mesh
        :type station_coordinates: class: '~str'
        :param: String containing station coordinats as in Input.par
        e.g.:
          '3.9463538e+06   1.9802029e+06   4.5486094e+06 \n
           1.9542831e+06   9.9590736e+05   5.9475798e+06 \n
           4.1565145e+06   1.4021579e+06   4.5762539e+06
           .
           .
           .'

        """

        self.filename_mesh = filename_mesh
        try:
            self.station_coordinates = self.__read_station_coordinates(
                station_coordinates)
        except:
            print 'cannot load station coordinates'
            raise Exception, ('cannot load station coordinates %s' %
                              station_coordinates)
        try:
            self.vertices, self.elements, self.boundary_elements = \
                self.__read_mesh(filename_mesh=self.filename_mesh)
        except:
            print 'cannot load mesh'
            raise Exception, ('cannot load mesh %s' % self.filename_mesh)
        try:
            self.surface_volume = self.__construct_surface_volume(
                self.vertices, self.elements, self.boundary_elements,
                surface_volume_thickness=0.97)
        except:
            print 'cannot load construct surface volume'
            raise Exception('cannot construct surface volume')

    def __read_station_coordinates(self, station_coordinates):
        """
        Converts a String of station coordinates as given in Seissol-Input.par
        Input-file into a numpy array.

        :type station_coordinates: class:'str()'
        :param: Converts
        :rtype points: class:'numpy.ndarray()'
        :param: 3xN-dim array containing X,Y,Z coordinates of N stations
        """

        station_coordinates = station_coordinates.split('\n')
        points = []
        for station in station_coordinates:
            if station in '':
                break
            points.append([float(coord) for coord in station.split()])
        points = np.asanyarray(points)
        return points

    def __read_mesh(self, filename_mesh='/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/eucrust_small_new.neu'):
        """
        From mesh-file read vertices, elements and boundary elements

        :type: filename_mesh: 'str'
        :param: filename of the mesh. Full path or relative path
        :rtype: vertices: '~numpy.ndarray'
        :param: 4xN-dim array (N=Number of vertices). 1st-dim=Index;
        2nd-dim= x-coordinate; 3rd-dim = y-coordinate;
        4th-dim = z-coordinate
        :rtype: elements: '~numpy.dnarray'
        :param: 7xN-dim array (N=Number of elements). 1st-dim=Index;
        2nd-dim= irrelevant; 3rd-dim = irrelevant;
        4th-7th-dim = Indices of vertices constructing the element
        :rtype: boundary_elements: 'dict(boundary_condition: numpy.array)'
        :param: Dictionary Refering to the Indices of Elements adhering \
        to boundary conditions. Keys are Boundary conditions type
        identifiers
        Identifier of boundary element type:
        101 = free surface 	    boundary
        102 = non-conforming 	boundary
        103 = dynamic rupture	boundary
        104 = inflow 		    boundary
        105 = absorbing 	 	boundary
        106 = periodic 		    boundary
        """

        filename = filename_mesh
        if filename[-4:] not in '.neu':
            print 'filename not in *.neu format!!'
        filename = filename[:-4]
        t = timeit.time.time()
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%        Read Gambit Data
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = open(filename + '.neu')
        print '\n'
        print '--------------------------------------------------------',\
            '---------------------------'
        print ' Reading data from: %s' % filename + '.neu'
        for i in range(6):
            junk = f.readline()
        tmp = f.readline()
        tmp = [int(s) for s in tmp.split() if s.isdigit()]
        NUMNP = tmp[0]  # %Number of Vertices
        NELEM = tmp[1]  # %Number of Elements
        NGRPS = tmp[2]  # %Number of Element Groups
        NBSETS = tmp[3]  # %Number of Boundary Condition Sets
        NDFCD = tmp[4]  # %Number of Dimensions (2 or 3)
        if 3 == NDFCD:

            # Vertices
            f.readline()
            f.readline()
            vertices = np.fromfile(file=f, count=NUMNP * 4, sep=" ")
            vertices = vertices.reshape(NUMNP, 4)
            vertices = vertices.transpose()

            # Elements
            f.readline()
            f.readline()
            position = f.tell()
            first_line = f.readline()
            first_line = [int(s) for s in first_line.split() if s.isdigit()]
            n = 7
            f.seek(position, 0)
            elements = np.fromfile(file=f, count=NELEM * n, sep=" ",
                                   dtype=int).reshape(NELEM, n).transpose()

            # loop over element groups. irrelevant for this routine, but i\
            # left it in in case someone wants to work with it
            for group in range(NGRPS):
                f.readline()
                f.readline()
                first_line = f.readline()
                first_line = [int(s) for s in first_line.split() if
                              s.isdigit()]
                nelem = first_line[1]
                f.readline()
                f.readline()
                tmp = np.fromfile(file=f, count=nelem, sep=" ")
                position = f.tell()

            # Boundary Elements
            f.seek(position, 0)
            boundary_elements = {}
            for boundaries in range(NBSETS):
                f.readline()
                f.readline()
                first_line = f.readline()
                first_line = [int(s) for s in first_line.split() if
                              s.isdigit()]
                nelem = first_line[2]
                boundary_elements[str(
                    first_line[0])] = np.fromfile(
                        file=f, count=nelem * 3, sep=" ", dtype=int).reshape(3,
                                                                             nelem)[0]
                position = f.tell()
        else:
            print 'only 3d tetrahedral meshes'
        self.vertices = vertices
        self.elements = elements
        self.boundary_elements = boundary_elements
        return vertices, elements, boundary_elements

    def __construct_surface_volume(self, vertices, elements,
                                   boundary_elements, surface_volume_thickness=0.05):
        """
        From Free surface boundary construct convex hull of
        surface volume

        :type vertices: class:'~numpy.ndarray'
        :param: 4xN-dim array (N=Number of vertices). 1st-dim=Index;
        2nd-dim= x-coordinate; 3rd-dim = y-coordinate;
        4th-dim = z-coordinate
        :type elements: class:'~numpy.dnarray'
        :param: 7xN-dim array (N=Number of elements). 1st-dim=Index;
        2nd-dim= irrelevant; 3rd-dim = irrelevant;
        4th-7th-dim = Indices of vertices constructing the element
        :type boundary_elements: class:'dict(boundary_condition: numpy.array)'
        :param: Dictionary Refering to the Indices of Elements adhering \
        to boundary conditions.
        :type surface_volume_thickness: class:'~float'
        :param: Given in percent. Defines the volume for which to
        calculate the convex hull below the surface.
        :rtype surface_volume: class:'~numpy.ndarray'
        :param: Collection of point coordinates of surface element-vertices
        and the same points shifted towards coordinate system origin by
        some percent defined in surface_volume_thickness
        """

        # select elements for boundary condition '101' = free surface
        boundary_vertices = elements[3:].transpose().take(
            boundary_elements['101']-1, axis=0)
        # extract coordinates of vertices
        boundary_points_coordinates = vertices[1:].transpose().take(
            boundary_vertices.flatten()-1, axis=0)

        # create a surface volume by adding the same points, just moved
        # slightly towards the coordinate system origin
        lower_bound_boundary_points = boundary_points_coordinates *\
            (1.0 - surface_volume_thickness)
        surface_volume = np.concatenate((boundary_points_coordinates,
                                         lower_bound_boundary_points), axis=0)
        return surface_volume

    def __in_hull(self, points, hull):
        """
        hull is Delaunay of scipy.spatial

        :type points: class:'~numpy.ndarray'
        :params: 3xN-dim numpy array containing the X,Y,Z coordinates of the
        N stations.
        :type hull: class:'~scipy.spatial.Delaunay'
        :params: convex hull around mesh surface
        :rtype: 1d-Boolean array stating if Point is inside the convex
        hull or not
        """

        return hull.find_simplex(points) >= 0

    def move_points_to_mesh(self, steplength=1.01):
        """
        Checks if the stations are located within the mesh.

        The stations outside the mesh will be moved incrementally towards
        the coordinate system origin, until they are in the mesh again.
        The initial position of each station is moved away from the origin
        by a factor of "steplength"
        Then moving back towards the origin, is done in 10ths of steplength

        If some coordinates could not be moved into the mesh, they will
        be ignored.

        :type steplength: class:'~float'
        :param: sets the inital offset of the station coordinates away
        from the origin
        :rtype points_in_mesh: class:'~numpy.ndarray'
        :param: 3xN-dim numpy array containing all coordinates that could
         successfully be moved into the mesh.
        """
        
      
        
        try:
            self.hull = Delaunay(self.surface_volume)
        except:
            raise Exception, "Cant create convex hull for mesh"

        # set inital values:
        finalstation_coordinates = np.copy(self.station_coordinates *
                                           steplength)
        inhull = self.__in_hull(finalstation_coordinates, self.hull)
        min_radius_finalstation_coordinates = 1.0
        min_radius_surface_volume = 0.0

        # move stations into the convex hull of the surface volume. stop,
        # in case their distance to the origin is smaller then the the
        # bottom layer of the surface volume
        while min_radius_finalstation_coordinates > \
            min_radius_surface_volume and \
                sum(inhull) != len(self.station_coordinates):
            inhull = self.__in_hull(finalstation_coordinates, self.hull)
            inv_in3d = np.column_stack((inhull, inhull, inhull))
            finalstation_coordinates = np.where(inv_in3d,
                                                finalstation_coordinates, finalstation_coordinates *
                                               (1 - steplength / 10))
            min_radius_finalstation_coordinates = min([np.sqrt(x * x + y * y + z * z)
                                                       for x, y, z
                                                       in finalstation_coordinates])
            min_radius_surface_volume = min([np.sqrt(x * x + y * y + z * z) for
                                             x, y, z
                                             in self.surface_volume])

        # only keep the points that could be moved into the mesh
        points_in_mesh = np.compress(inv_in3d.flatten(),
                                     finalstation_coordinates.
                                     flatten(
                                     )).reshape(
                                         len(np.compress(inv_in3d.flatten(),
                                                         finalstation_coordinates.flatten())) / 3, 3)

        return points_in_mesh
        
if __name__ == '__main__':
    from wfs_input_generator import InputFileGenerator
    from obspy.core import UTCDateTime
    gen = InputFileGenerator()
    seissol_example_path = '/home/msimon/git/wfs_input_generator_msimon00/wfs_input_generator/tests/data/seissol_example/'
    data_dir = '/home/msimon/git/wfs_input_generator_msimon00/wfs_input_generator/tests/data/'
    import glob
    gen.add_stations(glob.glob(data_dir+'*dataless*'))
    event = {"latitude": 48.9,"longitude": -2.3,"depth_in_km": 13.0,"origin_time":\
    UTCDateTime(2012, 4, 12, 7, 15, 48, 500000),"m_rr": -2.11e+18,"m_tt": \
    -4.22e+19,"m_pp": 4.43e+19,"m_rt": -9.35e+18,"m_rp": -8.38e+18,"m_tp": -6.44e+18}
    gen.add_events(event)
    gen.config.mesh = 'vercehpc'
    gen.config.model = 'PREM'
    gen.config.working_directory = seissol_example_path
    gen.config.max_time = 1000.0
    gen.config.number_of_processors = 16
    gen.write(format = 'seissol_1_0', output_dir = seissol_example_path)
