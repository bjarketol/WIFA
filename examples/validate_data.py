from windIO.utils.yml_utils import validate_yaml, load_yaml
from windIO.utils import plant_schemas_path

validate_yaml(
    "cases/windio_4turbines/wind_energy_system/system.yaml",
    plant_schemas_path + "wind_energy_system.yaml",
)
