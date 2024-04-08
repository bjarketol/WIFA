from windIO.utils.yml_utils import validate_yaml, load_yaml
from windIO.utils import plant_schemas_path
validate_yaml('examples/cases/windio_4turbines/wind_energy_system/FLOW_toy_study_wind_energy_system.yaml', plant_schemas_path + 'wind_energy_system.yaml')
