# __init__.py
#
import os

#### User Defined ###
cs_exe_path = ""  # required: add your path to cs exe
salome_exe_path = ""  # required: add your path to salome exe
#
cs_env_command = ""  # optional : add any environment command that must be loaded to run code_sature. Example: "module load Minigorge3 && conda activate my_env"
salome_env_command = (
    ""  # optional : add any environment that must be loaded to run salome
)
#
python_scripts_env_command = ""  # optional : add any environment that must be loaded to run the python prepprocessing and postprocessing scripts
python_scripts_exe = "/usr/bin/python3"  # optional : specify your python exe to run the preprocessing and postprocessing scripts

#### Do not modify
cs_api_path = os.path.dirname(__file__)
