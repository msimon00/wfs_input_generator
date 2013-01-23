import os

from wfs_input_generator import InputFileGenerator



gen = InputFileGenerator()
gen.add_events("wfs_input_generator/tests/data/event1.xml")
gen.add_stations(["wfs_input_generator/tests/data/dataless.seed.BW_FURT", "wfs_input_generator/tests/data/dataless.seed.BW_RJOB"])
gen.config.mesh={'mesh':"eucrust_small_new"}
gen.config.model={'model':"model_eucrust_small_new.dat"}

gen.config.parameter={'version':18, 'dimension':3, 'advection':0, 'advection_velocity':(1.0, 1.0, 1.0),\
        'anisotropy':0, 'anelasticity':0, 'poroelasticity':0, 
        'adjoint':0, 'material_reference_values':(3600, 9.0e10, 1.11e11),\
        'randomfield':0, 'sourcetype':50, 'source_file':'source.dat', 'sponge':0,\
        'meshgenerator':"Gambit3D-Tetra",'fine_output':0, 'restartfile':0,\
        'DGMethod':1, 'CK':0,'fluxmethod':0, 'DGCycle':1, 'basisfunction_degree':0,\
        'reconstructed_basisfunction_degree':0,'stencil_security_factor':0,\
        'reconstruction_type':0, 'exponent_r':0, 'coefficient_epsilon':0,\
        'linear_weight':0, 'limiter_security_factor':0, 'minspace_order':5,\
        'maxspace_order':5, 'pAdaptivity_file_name':'pAdaptivity_file_name',\
        'material_basis_function_order':1, 'courant_number':0.5, 'min_time_step':10000,\
        'rotational_output':0, 'rotation_components':(1, 1, 1),\
        'variable_output':(0, 0, 0, 0, 0, 0, 1, 1, 1),\
        'material_parameters_output':(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),\
        'output_character':0 , 'output_format':1, 'timestep_output':50,\
        'time_output':25, 'output_index':1, 'sampling':1, 'max_time':2000,\
        'max_iteration':1000,'max_wallclocktime':1e20, 'delay':0}

test_dir="/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/"

if os.path.exists(test_dir+gen.config.mesh['mesh']+".neu") and os.path.exists(test_dir+gen.config.mesh['model']):
    gen.write(format="seissol_1_0", output_dir="/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/")
else:
    msg= 'did not find mesh %s.neu or model %s ' % gen.config.mesh['mesh'], gen.config.mesh['model']
    raise IOError(msg)
    


