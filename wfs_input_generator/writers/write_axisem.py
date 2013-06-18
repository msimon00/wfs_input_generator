#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input file writer for PyAxi (Python wrapper to execute Axisem by Kasra Hosseini) or regular Axisem inputfiles

:copyright:
    Marek Simon (marek.simon@geophysik.uni-muenchen.de), 2013
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""


# Define the required configuration items. The key is always the name of the
# configuration item and the value is a tuple. The first item in the tuple is
# the function or type that it will be converted to and the second is the
# documentation.


REQUIRED_CONFIGURATION = {
    'seismogram_length': (float, "Duration of the simulated waveforms in seconds"),
    'background_model': (str, "prem : Isotropic continental PREM model; \
                        prem_solid: like 'prem', replace fluid outer \
                        core with vs=vp/sqrt(3);\
                        prem_onecrust: like 'prem' but extend lower crust\
                        to surface;\
                        prem_light: PREM without crust;\
                        prem_solid_light: like 'prem_light', but in fluid\
                        outer core vs=vp/sqrt(3);\
                        iasp91: Isotropic IASP91 model with PREM density;\
                        twolayer_solid: discontinuity at r0/2: vel./dens.\
                         increase by 20%;\
                        homomodel: homogeneous solid;"),
    'dominant_period': (float, "DOMINANT period [s]"),
    'number_of_processors': (int, "Number of processors to be used")}


DEFAULT_CONFIGURATION = {
    "address": ('..', str, "Path to where AXISEM is located"),
    "mesh_name": ('DEFAULT_MESHER', str, "Name for the mesh folder"),
    "solver_name": ('DEFAULT_SOLVER' , str, "Name for the solver folder"),
    "vtk_output": ('.false.', str, "Set to .true. to generate VTK output"),
    "time_step": (0.0, float, ""),
    "time_scheme": ('newmark2', str, "pick time scheme from: newmark2,\
    symplec4,ML_SO4m5,ML_SO6m7,KL_O8m17,SS_35o10"),
    "STF": ('heavis', str, "pick source time function from 'dirac_0',\
    'quheavi', 'gauss_0', 'gauss_1' (1st deriv),'gauss_2' (2nd), \
    'heavis'"),
    'receiver_components': ('enz', str, 'Determine which receiver \
                    components to rotate to. ENZ, SPH, CYL, XYZ, SRC'),
    'postprocessingHalfDuration': (5.0, float, 'Half duration for \
    convolution. 0.0 if no convolution'),
    'postprocessingSTF': ('gauss0', str, "Pick source time function to\
    convolve with in post processing from 'dirac_0', 'quheavi', 'gauss_0',\
    'gauss_1' (1st deriv),'gauss_2' (2nd),'heavis'"),
    'dv': ('disp', str, 'Select if you want to use displacement or \
    velocity seismograms'),
    'make_flag': ('gfortran', str, 'add required flags to \
    makemake.pl'),
    'mpi_compiler': ('mpif90 -O3 -fbacktrace', str, 'mpi compiler of \
    your local machine')}


def _write_PyAxi_config(config, event):
    
    """
    Returns a configuration file for PyAxi as String
    """
    
    PyAxiCfg = "# AXISEM input file for python interface (PyAxi.py)\n \
    # This file should be located in the same folder as PyAxi.py\n \
    \n \
    [general]\n \
    address = %-31s ; address where AXISEM is located\n \
    mesh_name = %-29s ; Name for the mesh folder\n \
    solver_name = %-27s ; Name for the solver folder\n \
    verbose = N                              ; for debugging purposes\n \
    " % (config.address, config.mesh_name, config.solver_name)
    # Irrelevant options
    PyAxiCfg += "\n \
    # IMPORTANT:\n \
    # 'new_mesh' flag controls all the required steps that you want for your\
     simulation: \n \
    # (except for post_processing and test)\n \
    \n \
    # Y: means that you want to run the AXISEM for a new mesh! (so MESHER \
    and SOLVER)\n \
    # N: to re-run the solver with different parameters but the same mesh \
    (so SOLVER) ---> obviously you need an old mesh!\n \
    # M: if you change the flag to M, you have a full control on all steps \
    (please refer to the bottom of this section!)\n \
    \n \
    new_mesh = Y                             \n \
    post_processing = Y                      ; perform the post processing \
    step\n \
    \n \
    # ATTENTION: if you change the 'new_mesh' flag to 'M' [Manual] \n \
    # you could change these flags in a way that you want\n \
    mesher = N                               ; to run the mesher\n \
    solver = N                               ; to run the solver\n \
    mesher_makefile = Y                      ; Create Mesher Makefile? (just\
    once!)\n \
    mesher_make = Y                          ; Run make clean; make for \
    mesher\n \
    mesher_move = Y                          ; Run movemesh.csh\n \
    solver_makefile = Y                      ; Create Solver Makefile?\n \
    solver_cp = Y                            ; defines the mesh as a header\n \
    solver_make = Y                          ; Run make clean; make for \
    mesher\n \
    \n \
    "
    # Define Options for compiling
    PyAxiCfg += "[mpi_netCDF]\n \
    make_flag = %-29s; add required flags to \
    makemake.pl\n \
    mpi_compiler = %-26s; mpi compiler of your local \
    machine\n \
    \n \
    netCDF = N                            \n \
    # include folder of where the libnetcdff.a is located (type locate \
    libnetcdff.a)\n \
    netCDF_LIBS = -lm -L $(HOME)/local/lib -lnetcdff \
    -Wl,-rpath,$(HOME)/local/lib\n \
    netCDF_INCLUDE = -I $(HOME)/local/include -I /usr/include\n \
    \n \
    " % (config.make_flag, config.mpi_compiler)
    
    # Meshing Options
    PyAxiCfg += "[mesher]\n \
    model = %-33s ; Background model: prem,\
    prem_solid etc\n \
    period = %-33s ; DOMINANT period [s]\n \
    no_proc = %-33s ; Number of processors to be used\n \
    vtk_output = %-33s ; write vtk output (may cause \
    memory problem for high frequencies (~2s))\n \
    \n \
    [solver]\n \
    no_simu = 1                              ; number of simulations. 1: \
    single Mij/f_i; 2: forces; 4: moment tensor\n \
    seis_length = %-33s ; seismogram length [s]\n \
    time_step = %-33s ; time step [s]. Put to 0.0 to \
    use mesher's suggestion (mesh_params.h)\n \
    time_scheme = %-33s ; time scheme: newmark2,symplec4,\
    ML_SO4m5,ML_SO6m7,KL_O8m17,SS_35o10\n \
    source_type = cmtsolut                   ; source file type: 'sourceparams',\
    'cmtsolut'\n \
    receiver_type = stations                 ; receiver file type: 'colatlon',\
    'stations','database'\n \
    save_XDMF = .false.                      ; save XDMF files (high resolution\
    2D wavefields), more options in inparam_xdmf\n \
    force_aniso = .false.                    ; force anisotropic model \
    handling\n \
    viscoelastic_attenuation = .false.       ; include viscoelastic \
    attenuation\n \
    \n \
    " % ("'"+config.background_model+"'", config.dominant_period, \
    config.number_of_processors, "'"+config.vtk_output+"'", \
    config.seismogram_length,config.time_step, config.time_scheme)
    #There are two options on how to implement a Source, but we only use the\
    #CMT solution option.
    PyAxiCfg += "#--------------------------------------------------------\
    ---------------\n \
    # sourceparams parameters: (not valid for CMTSOLUTION)\n \
    sourceparams_type = 'monopole'      ; excitation type: 'monopole', \
    'dipole', 'quadpole'\n \
    sourceparams_MDQ = 'explosion'\n \
    # 'explosion','mxx_p_myy','mzz','vertforce' (MONOPOLE)\n \
    # 'mxz', 'myz', 'xforce', 'yforce'          (DIPOLE)\n \
    # 'mxy', 'mxx_m_myy'                        (QUADRUPOLE)\n \
    \n \
    src_Mzz = 1.E20\n \
    src_Mxx = 1.E20\n \
    src_Myy = 1.E20\n \
    src_Mxz = 0.E20\n \
    src_Myz = 0.E20\n \
    src_Mxy = 0.E20\n \
    \n \
    source_dp = 100.                         ; source depth [km]\n \
    source_colat = 0.0                       ; source colatitude [degrees]\n \
    source_lon = 0.0                         ; source longitude [ldegrees]\n \
    source_stf = 'dirac_0'                   ; source time function\n \
    #-----------------------------------------------------------------------\n \
    "
    #Options for Event information (CMT solution)
    PyAxiCfg += "# cmtsolut parameters: (not valid for sourceparams)\n \
    cmt_STF = %-33s ; 'dirac_0', 'quheavi', 'gauss_0', 'gauss_1' (1st deriv),\
    'gauss_2' (2nd), 'heavis' \n \
    \n \
    cmt_lat = %-33s ; source latitude [degrees]\n \
    cmt_lon = %-33s ; source longitude [degrees]\n \
    cmt_dp = %-33s ; source depth [km]\n \
    \n \
    cmt_Mrr = %33.1e ; Mrr component\n \
    cmt_Mtt = %33.1e ; Mtt component\n \
    cmt_Mpp = %33.1e ; Mpp component\n \
    cmt_Mrt = %33.1e ; Mrt component\n \
    cmt_Mrp = %33.1e ; Mrp component\n \
    cmt_Mtp = %33.1e ; Mtp component\n \
    \n \
    " % (config.STF, event['latitude'], event['longitude'], event['depth_in_km'],\
    event['m_rr'], event['m_tt'], event['m_pp'], event['m_rt'], \
    event['m_rp'], event['m_tp'])
    #Options for post processing
    PyAxiCfg += "[post_processing]\n \
    post_rotate = .true.                     ; rotate receivers?\n \
    post_components = %-23s ; receiver components: enz,sph,cyl,xyz,src\n \
    post_full_Mij = .true.                   ; sum to full Mij\n \
    post_conv_period = %-22s ; convolve period (0. if not convolved)\n \
    post_STF = %-30s ; source time function type for convolution\n \
    post_dv = %-31s ; disp or velo seismograms\n \
    \n \
    " % ("'" + config.receiver_components + "'", config.postprocessingHalfDuration,\
    config.postprocessingSTF, config.dv)
    #No additional post processing necessary
    PyAxiCfg += "post_Mrr = N                             ; Mrr\n \
    post_Mtt = N                             ; Mtt\n \
    post_Mpp = N                             ; Mpp\n \
    post_Mrt = N                             ; Mrt\n \
    post_Mrp = N                             ; Mrp\n \
    post_Mtp = N                             ; Mtp\n \
    \n \
    post_Scolat = N                          ; Source colatitude\n \
    post_Slon = N                            ; Source longitude\n \
    post_snap = N                            ; plot global snaps?\n \
    post_path = N                            ; Directory for post processed data\n \
    post_negative = N                        ; seismograms at negative time (0 at max. of stf)\n \
    \n \
    "
    # Return Miniseed, no additional processing neccessary
    PyAxiCfg += "[MISC]\n \
    mseed = Y                                ; convert the seismograms into MSEED format\n \
    mseed_all = N                            ; convert the seismograms into MSEED format\n \
    \n \
    convSTF = N                              ; after converting the seismograms into MSEED, convolve with STF\n \
    halfduration = 5.0                       ; halfduration (do not confuse that with halfduration in the test_section)\n \
    \n \
    filter = N                               ; if you want to appy a filter (lowpass and highpass as defined below)\n \
    fmin = 0.012                             ; minimum frequency\n \
    fmax = 0.1                               ; maximum frequency\n \
    \n \
    "
    # Irrelevant options for testing only
    PyAxiCfg += "[testing]\n \
    # Parameters required for the TESING directoryi\n \
    # for running the tests automatically please refer to test_axisem.py\n \
    test = N                                 ; if you want to test AXISEM?\n \
    test_folder = ./automated/test_99        ; address of the test folder\n \
    plot = Y                                 ; to just plot the test waveforms\n \
    save_plots = N                           ; save the plots to a file\n \
    plot_format = png                        ; file format of the plots\n \
    chans = ['Z', 'N']                       ; required channels\n \
    fmin = 0.0005                            ; minimum frequency\n \
    fmax = 0.02                              ; maximum frequency\n \
    halfduration = 20.                       ; halfduration for the Source Time Function\n \
    nstat = 20"
    return PyAxiCfg

def _write_inparam_mesh(config):
    
    """
    Returns input file for Axisem Mesher as string
    """
    
    inparam_mesh = "####### INPUT FILE FOR AXISEM MESHER ##########################################\n \
    %-20s Background model: 'prem','prem_solid' etc (SEE BELOW)\n \
    %-20s DOMINANT period [s]\n \
    %-20s Number of processors to be used\n" % (config.background_model, \
                                            config.dominant_period, \
                                            config.number_of_processors)
    inparam_mesh += "######## Don't change anything below unless you know what you're doing ######## \n \
    .true.              Resolve inner-core shear wave (ignored if not prem)?\n \
    4		    Polynomial order of basis for GLL points inside elements\n \
    1.5d0               Number of elements per DOMINANT wavelength\n \
    0.6d0               Courant number\n \
    6.371e+6            Surface radius [m]\n \
    .false.             Save mesh files (WARNING: large)\n \
    .false.             Write excessive info output (WARNING: for the inclined)\n \
    Diags               Where to dump mesh output files\n \
    3                   Number of expected coarsening levels. Guideline formula:\n \
                             log( vs(cmb)/vs(surf)*r(surf)/r(cmb) ) / log(2)\n \
    ###############################################################################\n \
    #### END OF INPUT #############################################################\n \
    #### BACKGROUND MODELS:\n \
    prem :              Isotropic continental PREM model\n \
    prem_solid:         like 'prem', replace fluid outer core with vs=vp/sqrt(3)\n \
    prem_onecrust:      like 'prem' but extend lower crust to surface\n \
    prem_light:         PREM without crust\n \
    prem_solid_light:   like 'prem_light', but in fluid outer core vs=vp/sqrt(3)\n \
    iasp91:             Isotropic IASP91 model with PREM density\n \
    twolayer_solid:     discontinuity at r0/2: vel./dens. increase by 20%\n \
    homomodel:          homogeneous solid\n \
    ###############################################################################\n"
    return inparam_mesh

def _write_inparam(config):
    
    """
    Returns input file for Axisem solver as string
    """
    
    inparam = "######## INPUT FILE FOR AXISEM (SOLVER): GLOBAL WAVE PROPAGATION IN SNREIs #############\n \
    1                number of simulations. 1: single Mij/f_i; 2: forces; 4: moment tensor \n \
    \n \
    ######## TIME LENGTH/SCHEME CHOICES ####################################################\n \
    %-16s seismogram length [s]\n \
    %-16s time step [s]. Put to 0.0 to use mesher's suggestion (mesh_params.h)\n \
    %-16s time scheme: newmark2,symplec4,ML_SO4m5,ML_SO6m7,KL_O8m17,SS_35o10\n \
    \n \
    " % (config.seismogram_length, config.time_step, config.time_scheme)
    inparam += "######## SOURCE AND RECEIVER SPECIFICATION #############################################\n \
    0.0          source period [s]. Put to 0.0 to use mesh resolution (mesh_params.h)\n \
    cmtsolut     source file type: 'sourceparams','cmtsolut'\n \
    stations     receiver file type: 'colatlon','stations'\n \
    0.0          desired seismogram sampling rate [s] (0.0 for equal to time step (see above)\n \
    \n \
    " % (config.sampling_rate)
    inparam += "######## OUTPUT: PATHS AND SNAPSHOT WAVEFIELDS #########################################\n \
    ./Data           data output path: meshes, wavefields, seismograms, related info\n \
    ./Info           info output path: summary, test and information files\n \
    %-16s save global snapshots for wavefield movies (WARNING: Large!)\n \
    20.0             approximate time intervall between snapshots [s]\n \
    \n \
    " % (config.vtk_output)
    inparam += "######## OUTPUT: SENSITIVITY KERNELS ###################################################\n \
    .false.          save wavefields for sensitivity kernels? (WARNING: Large!)\n \
    8                samples per period between wavefield snapshots for kernels\n \
    mask             source vicinity in wavefields? 'igno','mask','anal'\n \
    1,1              starting (ibeg) and ending (npol-iend) GLL-point index \n \
    \n \
    ######## MISCELLANEOUS PARAMETERS ######################################################\n \
    .false.         save global kinetic/potential energy? Generally not needed\n \
    .false.         overwrite background model with homogeneous parameters? \n \
    10.,5.77,3.    	homogeneous P-vel [km/s]; S-vel [km/s]; rho [g/cm^3]\n \
    .false.         calculate analytical solution around source (for homogeneous model)?\n \
    .false.         add heterogeneity? if true, specify in inparam_hetero (-> get_model)\n \
    .false.         do mesh tests? Suggested to do once per mesh and unedited source code\n \
    .false.         save large test files (valence etc.)? Generally not needed\n \
    binary          Output format for seismograms and wavefields: binary, netcdf\n \
    "
    return inparam

def _write_axisem_station(stations):
    #writes STATION file for input for axisem
    stat = set()
    vals = list()
    #arrange station information
    for station in stations:
        vals.append((station['id'].split('.')[1], station['id'].split('.')[0],\
        station['latitude'], station['longitude'], station['local_depth_in_m'],\
        station['elevation_in_m']))
    stations_input = ''
    for v in vals:
        outputformat = ('%5s %4s %12.5f %17.5f %17.8f %17.8f \n')
        stations_input += outputformat % v
    return stations_input

def _write_cmtsolution(config, event):
    try:
        halfduration=event['half_duration']
    except:
        halfduration=config.dominant_period
    cmtsolution = '%s %s %s %s %s' % (config.STF, event['origin_time'].formatArcLink(),
    event['latitude'], event['longitude'], event['depth_in_km']) + '\n'
    cmtsolution += 'event name: ' + event['origin_time'].formatArcLink() + '\n'
    cmtsolution += 'time shift: ' + '0.0000' + '\n'
    cmtsolution += 'half duration: ' + str(halfduration) + '\n'
    cmtsolution += 'latitude: ' + str(event['latitude']) + '\n'
    cmtsolution += 'longitude: ' + str(event['longitude']) + '\n'
    cmtsolution += 'depth: ' + str(event['depth_in_km']) + '\n'
    cmtsolution += 'Mrr: ' + str(event['m_rr']) + '\n'
    cmtsolution += 'Mtt: ' + str(event['m_tt']) + '\n'
    cmtsolution += 'Mpp: ' + str(event['m_pp']) + '\n'
    cmtsolution += 'Mrt: ' + str(event['m_rt']) + '\n'
    cmtsolution += 'Mrp: ' + str(event['m_rp']) + '\n'
    cmtsolution += 'Mtp: ' + str(event['m_tp']) + '\n'
    return cmtsolution

def _check_configuration(config):
    """
    not yet implemented
    """
    pass

def write(config, events, stations):
    if len(events) > 1:
        raise TypeError, 'Can only simulate one event at a time!'
    event = events[0]
    output_files = {}
    _check_configuration(config)
    output_files['PyAxi.cfg'] = _write_PyAxi_config(config, event)
    output_files['inparam_mesh'] = _write_inparam_mesh(config)
    output_files['inparam'] = _write_inparam_mesh(config)
    output_files['STATIONS'] = _write_axisem_station(stations)
    output_files['CMTSOLUTION'] = _write_cmtsolution(config, event)
    for key in output_files.iterkeys():
        output_files[key] += "\n"
    return output_files

