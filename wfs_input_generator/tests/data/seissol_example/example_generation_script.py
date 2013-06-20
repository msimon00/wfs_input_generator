from wfs_input_generator import InputFileGenerator
from obspy.core import UTCDateTime
gen = InputFileGenerator()
seissol_example_path = 'wfs_input_generator/tests/data/seissol_example/'
data_dir = 'wfs_input_generator/tests/data/'
import glob
gen.add_stations(glob.glob(data_dir+'*dataless*'))
event = {"latitude": 48.9,"longitude": -2.3,"depth_in_km": 200.0,"origin_time":\
 UTCDateTime(2012, 4, 12, 7, 15, 48, 500000),"m_rr": -2.11e+18,"m_tt": \
 -4.22e+19,"m_pp": 4.43e+19,"m_rt": -9.35e+18,"m_rp": -8.38e+18,"m_tp": -6.44e+18}
gen.add_events(event)
gen.config.mesh = 'most_simple_tet'
gen.config.model = 'PREM'
gen.config.working_directory = seissol_example_path
gen.config.max_time = 1000.0
gen.config.number_of_processors = 16
gen.write(format = 'seissol_1_0', output_dir = seissol_example_path)
