#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Input file write for Seissol 1.0

:copyright:
    Marek Simon 2013
:license:
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""
import os

def write(config, events, stations, output_directory, parameter={}, mesh={}, model={}):
    """
    Seissol
    """
    input_files = dict()
    
    #Default Parameters
    PARAMS={'model':"model_eucrust_small_new.dat", 'mesh':"eucrust_small_new",\
        'version':18, 'dimension':3, 'advection':0, 'advection_velocity':(1.0, 1.0, 1.0),\
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
    PARAMS.update(parameter)
    PARAMS.update(model)
    PARAMS.update(mesh)
    
    PARAMS['nbfreesurface'], PARAMS['nbopensurface']\
    =_check_mesh(output_directory+PARAMS['mesh'])
    _check_input_par(PARAMS)
    PARAMS['station_coordinatesXYZ']=_prepare_station_coordinates(stations)
    
    input_files['input.par']=_input_par(PARAMS, stations)
    input_files['source.dat']=_source_dat(event, sampling=PARAMS['sampling'], function='dirac')

    for key in input_files:
        output_file = os.path.join(output_directory, key)
        with open(output_file, "wb") as open_file:
            open_file.write(input_files.get(key))

def _check_mesh(file_path):
    if 'Gambit' in PARAMS['meshgenerator']:
        try:
            f=open(file_path+'.neu','r')
        except:
            msg = "Could not open %s." % file_path
            raise IOError(msg)
        try:
            for i in range(6):
                junk=f.readline()
                tmp=f.readline()
                f.close()
                tmp=[int(s) for s in tmp.split() if s.isdigit()]
                nbsets=tmp[3]
                nbfreesurface=1
                nbopensurface=nbsets-nbfreesurface
                return nbfreesurface, nbopensurface
        except:
            msg "Could not read %s." %file_path
            raise IOError(msg)
    else:
        msg= "Can only handle Gambit *.neu meshes so far!"
        raise TypeError(msg)
        return 1,0
    
def _check_input_par(PARAMS):
    #check seissol options:
    if PARAMS['advection'] is 1:
        PARAMS['advection_velocity']=PARAMS['advection_velocity']+'\n'
    else:
        PARAMS['advection_velocity']=''
    if PARAMS['dimension'] is 2:
        PARAMS['anisotropy']=0
    if PARAMS['dimension'] is 3:
        PARAMS['anelasticity']=0
    #self.poroelasticity=poroelasticity #0=no poroelasticity, 1=Biot's theory to describe 
    #wave propagation through porous materials filled with a fluid. For the
    #high-frequency range, option 1 should be chosen which uses the standard ADERDG
    #For the low-frequency case, option 2 uses the ADER-DG(ST)
    #and option 3 the ADER-DG(FS) method.
    if PARAMS['model'] in 'constant':
        PARAMS['material']=0
    elif PARAMS['model'] in 'PREM':
        PARAMS['material']=2 #0=constant, 1=from file, 2=PREM, 3=linear
    elif PARAMS['model'] in 'linear':
        PARAMS['material']=3
    else:
        PARAMS['material']=1
        print 'Use Material File as Model: ',PARAMS['model']
        PARAMS['materiaLfile_name']=PARAMS['model']
    if PARAMS['anisotropy'] is 1 and PARAMS['material'] is not 1:
        msg= 'To include anisotropy choose material file'
        raise IOError(msg)
    if PARAMS['poroelasticity'] is not 0 and PARAMS['material'] is not 1:
        msg 'To include poroelasticity choose material file'
        raise IOError(msg)
    try:
        PARAMS['material_reference_values']="%-10.3e %-10.3e %-10.3e"\
         % PARAMS['material_reference_values']
    except:
        msg="no material values defined, using default"
        raise UserWarning(msg)
        print "default material reference values!!"
        print "mv_1=3600."
        print "mv_2=9.0e10"
        print "mv_3=1.116e11"
        mv_1=3600.
        mv_2=9.0e10
        mv_3=1.116e11
        PARAMS['material_reference_values']="%-10.3e %-10.3e %-10.3e" % (mv_1, mv_2, mv_3)
    if PARAMS['material'] is 1:
        if len(PARAMS['model'])<32:
            PARAMS['materialfile_name']="%-37s \n" % PARAMS['model']
        else:
            msg= "material file name too long!! <32char!!"
            raise NameError(msg)
    else:
        PARAMS['materialfile_name']=''
    if PARAMS['poroelasticity'] is 3:
        PARAMS['CK']=0
    elif PARAMS['poroelasticity'] is 2:
        PARAMS['CK']=1
    if PARAMS['DGMethod'] is 3:
        PARAMS['DGcycle']=str(PARAMS['DGcycle'])+" \n" # 1: one iteration = one cycle; 2: one iteration = update of all elements
    else:
        PARAMS['DGcycle']=''
    if PARAMS['DGMethod'] is 4:
        PARAMS['basisfunction_degree']="%-37s \n" % PARAMS['basisfunction_degree']
        PARAMS['reconstructed_basisfunction_degree']="%-37s \n" % PARAMS['reconstructed_basisfunction_degree']
        PARAMS['stencil_security_factor']="%-37s \n" % PARAMS['stencil_security_factor']
        PARAMS['reconstruction_type']="%-37s \n" % PARAMS['reconstruction_type'] #0 (linear reconstruction) or ELSE (nonlinear central reconstruction)
        PARAMS['exponent_r']="%-37s \n" % PARAMS['exponent_r']
        PARAMS['coefficient_epsilon']="%-37s \n" % PARAMS['coefficient_epsilon']
        PARAMS['linear_weight']="%-37s \n" % PARAMS['linear_weight']
        PARAMS['limiter_security_factor']="%-37s \n" % PARAMS['limiter_security_factor']
    else:
        PARAMS['basisfunction_degree']=''
        PARAMS['reconstructed_basisfunction_degree']=''
        PARAMS['stencil_security_factor']=''
        PARAMS['reconstruction_type']=''
        PARAMS['exponent_r']=''
        PARAMS['coefficient_epsilon']=''
        PARAMS['linear_weight']=''
        PARAMS['limiter_security_factor']=''
    if PARAMS['minspace_order'] is not PARAMS['maxspace_order']:
        PARAMS['zone_order_flag']=str(0)+' \n' #0 ordering depends on dimension(size of element)
        if str(1) in PARAMS['zone_order_flag']:
            try:
                PARAMS['pAdaptivity_file_name']="%-37s \n" % PARAMS['pAdaptivity_file_name']
            except:
                msg='no pAdaptivity_filename: %s' % PARAMS['pAdaptivity_file_name']
                raise UserWarning(msg)
                PARAMS['pAdaptivity_file_name']=''
        else:
            PARAMS['pAdaptivity_file_name']=''
    else:
        PARAMS['zone_order_flag']=''
        PARAMS['pAdaptivity_file_name']=''
    if PARAMS['rotational_output'] is 1:
        PARAMS['rotation_components']=("%s %s %s \n") % PARAMS['rotation_components']
        #rotation_components=(1, 1, 1) #AXES X, Y, Z
    else:
        PARAMS['rotation_components']=''
    PARAMS['variable_output']=("%s %s %s %s %s %s %s %s %s") % PARAMS['variable_output']
    #variable_output=(0, 0, 0, 0, 0, 0, 1, 1, 1) # up to 9
    PARAMS['material_parameters_output']=("%s %s %s %s %s %s %s %s %s %s %s") % PARAMS['material_parameters_output']
    #material_parameters_output=(1, 1, 1) # up to 11
    PARAMS['station_coordinatesXYZ']=self._prepare_station_coordinates()
    #PARAMS['number_of_pickpoints']=len(self.station_coordinatesXYZ.split('\n'))-1
    #Abort Criteria
    PARAMS['max_time']=PARAMS['max_time']# maximum computing time  \n",
    PARAMS['max_iteration']=PARAMS['max_iteration'] # number of maximum iteration
    PARAMS['max_wallclocktime']=PARAMS['max_wallclocktime'] #Max. wallclock time
    PARAMS['delay']=PARAMS['delay'] # Delay
    #END OF SEISSOL INPUT
        
def _input_par(PARAMS, stations):  
    seissol_input= ''.join([ \
    "<-----------------------------------> \n",\
    "<  S E I S S O L - V E R S I O N    > \n",\
    "<-----------------------------------> \n",\
    "%-37s! seissol-version \n" % PARAMS['version'] ,
    "\n",
    "<-----------------------------------> \n",\
    "<  E Q U A T I O N S                > \n",\
    "<-----------------------------------> \n",\
    "%-37s! number of dimensions ( 2 or 3 ) \n" % PARAMS['dimension'],\
    "%-37s! Advection \n" % PARAMS['advection'],\
    str(PARAMS['advection_velocity']),\
    "%-37s! Anisotropy (only in 3D) \n" % PARAMS['anisotropy'],\
    "%-37s! Anelasticity (only in 2D) \n" % PARAMS['anelasticity'],\
    "%-37s! Poroelasticity \n" % PARAMS['poroelasticity'],\
    "%-37s! adjoint \n" % PARAMS['adjoint'],\
    "%-37s! Material (0=constant;1=from file;2=PREM;3=linear) \n" % PARAMS['material'],\
    "%-37s! Reference values \n" % PARAMS['material_reference_values'],\
    PARAMS['materialfile_name'],\
    "%-37s! random field \n" % PARAMS['randomfield'],\
    " \n",\
    ###########################################################################################ASK STEFAN/CHRISTIAN FOR OPTIONS!!! NOT DOCUMENTED!!!
    "<-----------------------------------> \n",\
    "<  INITIAL CONDITION                > \n",\
    "<-----------------------------------> \n",\
    "Gauss_Puls_Rad                       ! problem-dependent initial values \n",\
    "1                                    ! Not needed for this IC \n",\
    "0.0  0.0  0.0                        ! Center of Gausspulse (x,y,z) \n",\
    "0.                                   ! Amplitudes \n",\
    "5.0e3  5.0e3  5.0e3                  ! Halfwidths (hx>0,hy>0,hz>0) \n",\
    " \n",\
    "<-----------------------------------> \n",\
    "<  BOUNDARIES                       > \n",\
    "<-----------------------------------> \n",\
    #######################check meshfile: 101 and 105 boundaries
    "%-37s! number of free surfaces \n" % PARAMS['nbfreesurface'],\
    "0                                    ! number of no-flux boundaries \n",\
    "0                                    ! number of fault boundaries \n",\
    "0                                    ! number of inflow boundaries \n",\
    "%-37s! number of open (outflow) boundaries \n" % PARAMS['nbopensurface'],\
    "0                                    ! periodic boundaries \n",\
    " \n",\
    ###########################################################################################ASK STEFAN/CHRISTIAN for source options
    "<-----------------------------------> \n",\
    "<  S O U R C E    T E R M S         > \n",\
    "<-----------------------------------> \n",\
    "%-37s! Rupture plane source \n" % PARAMS['sourcetype'],\
    "%-37s! rupture plane input file \n" % PARAMS['source_file'],\
    " \n",\
    "<-----------------------------------> \n",\
    "<  S P O N G E    L A Y E R         > \n",\
    "<-----------------------------------> \n",\
    "%-37s! \n" % PARAMS['sponge'],\
    " \n",\
    "<-----------------------------------> \n",\
    "<  M E S H                          > \n",\
    "<-----------------------------------> \n",\
    "%-37s! name of mesh file \n" % PARAMS['mesh'],\
    "%-37s! name of meshgenerator (format) \n" % PARAMS['meshgenerator'],\
    "0                                    ! Periodic boundaries \n",\
    "0.        0.       0.                ! Mesh displacement   (2D or 3D) \n",\
    "1.        0.       0.                ! Mesh scaling matrix (2D or 3D) \n",\
    "0.        1.       0.                ! Mesh scaling matrix (2D or 3D) \n",\
    "0.        0.       1.                ! Mesh scaling matrix (2D or 3D) \n",\
    " \n",\
    "<-----------------------------------> \n",\
    "<  D I S C R E T I S A T I O N      > \n",\
    "<-----------------------------------> \n",\
    "%-37s! fine output \n" % PARAMS['fine_output'],\
    "%-37s! Restart file \n" % PARAMS['restartfile'],\
    "%-37s! DGmethod \n" % PARAMS['DGMethod'],\
    "%-37s! (0=CK; 1=STDG) \n" % PARAMS['CK'],\
    "%-37s! flux 0-godunov 1-rusanov \n" % PARAMS['fluxmethod'],\
    PARAMS['DGcycle'],\
    PARAMS['basisfunction_degree'],\
    PARAMS['reconstructed_basisfunction_degree'],\
    PARAMS['stencil_security_factor'],\
    PARAMS['reconstruction_type'],\
    PARAMS['exponent_r'],\
    PARAMS['coefficient_epsilon'],\
    PARAMS['linear_weight'],\
    PARAMS['limiter_security_factor'],\
    "%-37s! Min Space accuracy \n" % PARAMS['minspace_order'],\
    "%-37s! Max Space accuracy \n" % PARAMS['maxspace_order'],\
    PARAMS['zone_order_flag'],\
    PARAMS['pAdaptivity_file_name'],\
    "%-37s! Material basis function order \n" % PARAMS['material_basis_function_order'],\
    "%-37s! courant number (<=1.0) \n" % PARAMS['courant_number'],\
    "%-37s! Manualy chosen minimum time \n" % PARAMS['min_time_step'],\
    ' \n',\
    "<-----------------------------------> \n",\
    "<  O U T P U T                      > \n",\
    "<-----------------------------------> \n",\
    'out','\n',\
    "%-37s! rotational output ( 0=no / 1=yes ) \n" % PARAMS['rotational_output'],\
    PARAMS['rotation_components'],\
    "%-37s! variables \n" % PARAMS['variable_output'],\
    "%-37s! material values (rho / mu / lambda) \n" % PARAMS['material_parameters_output'],\
    "%-37s! character of output (0=geo.based, 1=calc. based) \n" % PARAMS['output_character'],\
    "%-37s! format of output (0=IDL, 1=TECPLOT, 2=IBM DX, 4=GiD)) \n" % PARAMS['output_format'],\
    "%-37s! index of printed info at timesteps \n" % PARAMS['timestep_output'],\
    "%-37s! index of printed info at time \n" % PARAMS['time_output'],\
    "%-37s! Criterion for index of printed info: 1=timesteps,2=time,3=timesteps+time \n" % PARAMS['output_index'],\
    "%-37s! Pickpoint Sampling \n" % PARAMS['sampling'],\
    "%-37s! # of coordinates to pick \n" % (len(PARAMS['station_coordinatesXYZ'].split('\n'))-1),\
    PARAMS['station_coordinatesXYZ'],\
    str(0)+' \n',\
    ' \n',\
    "<-----------------------------------> \n",\
    "<  A B O R T - C R I T E R I A      > \n",\
    "<-----------------------------------> \n",\
    "%-37s! maximum computing time \n" % PARAMS['max_time'],\
    "%-37s! number of maximum iteration \n" % PARAMS['max_iteration'],\
    "%-37s! Max. wallclock time \n" % PARAMS['max_wallclocktime'],\
    "%-37s! Delay \n" % PARAMS['delay'],\
    "\n",\
    "<-----------------------------------> \n",\
    "<  A N A L Y S I S  O F   D A T A   > \n",\
    "<-----------------------------------> \n",\
    "0                                    ! Do not analyse data \n",\
    "\n",\
    "<-----------------------------------> \n",\
    "<  D E B U G G I N G     M O D U S  > \n",\
    "<-----------------------------------> \n",\
    "0                                     ! No debugging \n",\
    "\n"])
    return seissol_input

def _prepare_station_coordinates(stations):
    import scipy as sp
    #ugly fix: better check mesh and check if stations are inside!!!:
    deviation_from_sphere=6371000/200 #see also function _source_dat()
    station_coordinates=str()
    coordinate_format=('  %1.7e   %1.7e   %1.7e \n')#      #%s\n')
    receivers=list()
    for station in stations:
        lat = station['latitude'] *np.pi/180
        lon = station['longitude']*np.pi/180
        r = 6371000 - deviation_from_sphere # - dep  # HACK TO MAKE EXAMPLE WORK. SUBTRACTS ESTIMATE OF MAXIMUM DEVIATION FROM SPHERICAL SURFACE
        x = r*sp.cos(lat)*sp.cos(lon)
        y = r*sp.cos(lat)*sp.sin(lon)
        z = r*sp.sin(lat)
        receivers.append(coordinate_format % (x,y,z,))
    station_coordinates=' '.join(list(set(receivers)))
    return station_coordinates

def _source_dat(event, sampling, function):
    import scipy as sp
    #event
    deviation_from_sphere=6371000/200
    lat = event['latitude'] *np.pi/180;
    lon = event['longitude']*np.pi/180;
    dep = event['depth_in_km']    *1000;
    r   = 6371000 - abs(dep); # HACK TO MAKE EXAMPLE WORK. SUBTRACTS ESTIMATE OF MAXIMUM DEVIATION FROM SPHERICAL SURFACE
    x = r*sp.cos(lat)*sp.cos(lon);
    y = r*sp.cos(lat)*sp.sin(lon);
    z = r*sp.sin(lat);
    #Coordinate transformation matrix R_ij:
    R_11= sp.cos(lat)*sp.cos(lon)
    R_12= sp.cos(lat)*sp.sin(lon)
    R_13= sp.sin(lat)
    R_21= sp.sin(lat)*sp.cos(lon)
    R_22= sp.sin(lat)*sp.sin(lon)
    R_23=-sp.cos(lat)
    R_31=-sp.sin(lon)
    R_32= sp.cos(lon)
    R_33= 0.0
    R = np.matrix([[R_11, R_12, R_13],\
                   [R_21, R_22, R_23],\
                   [R_31, R_32, R_33]]);
    #CMT as matrix
    Mrr=event['m_rr']
    Mrt=event['m_rt']
    Mrp=event['m_rp']
    Mtt=event['m_tt']
    Mtp=event['m_tp']
    Mpp=event['m_pp']
    M  = np.matrix([[Mrr, Mrt, Mrp],\
                    [Mrt, Mtt, Mtp],\
                    [Mrp, Mtp, Mpp]])
    #Coordinate Transformation of Momenttensor
    M_rot=R.getI()*M*R
    # SAVE IN SEISSOL FORMAT
    source_dat='Seismic Moment Tensor\n'
    for row in M_rot:
        source_dat+='%2.9e %2.9e %2.9e \n' % (row.tolist()[0][0], row.tolist()[0][1], row.tolist()[0][2])
    source_dat+='Number of subfaults\n'
    source_dat+='1\n'
    source_dat+='               x                y                z         strike        dip             rake            area       Onset time\n'
    source_dat+='%16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e \n' % \
           (     x,     y,     z,   0.0,   0.0,   0.0,   1.0,   0.0     )
    if function in 'dirac':
        source_dat+="source time function \n"
        source_dat+="2.000000000e-02  3  \n"
        source_dat+="samples  \n"
        source_dat+="0.000000000e+00  \n"
        source_dat+="1.000000000e+00  \n"
        source_dat+="0.000000000e+00  \n"
    if function in 'gauss':
        ssampling=float(sampling)/1000
        sd=sampling
        src       = np.arange(0,sd*6,ssampling)
        def gauss(x):
            x=1.0/sd/np.sqrt(2*np.pi)*np.exp(-1/2*(x-6*sd/2)/sd*(x-6*sd/2)/sd)
            return x
        vecgauss=np.vectorize(gauss, otypes=[float])
        src=vecgauss(src)
#        for i,j in fl: #np.nditer(src, op_flags=['readwrite']):
#            x=1.0/hd/np.sqrt(2*np.pi)*np.exp(-1/2*(x-sd/2)/hd*(x-sd/2)/hd)
        src_nsamp = len(src);
        src_dt    = ssampling;
        source_dat+='source time function\n'
        source_dat+='%16.9e %16.9e \n' % (src_dt,src_nsamp)
        source_dat+='samples \n'
        for value in src:
            source_dat+='%16.9e\n' % value
    return source_dat
