import sys
from os import path, environ, mkdir
import os as os
#from mpi4py import MPI
import string, random
from os import sep, mkdir, walk
from shutil import copy, copytree
import sys
import numpy as np
from os import path, chdir, environ, getcwd
import yaml
#from yaml.loader import SafeLoader
from datetime import datetime, timedelta
from windIO.utils.yml_utils import validate_yaml, Loader, load_yaml
from functools import reduce
import csMeteo.nieuwstadt_stable_profiles_utils as nwstdt

def theta2temp(theta,z_or_dz,P0,Pref=1000.,g=9.81,Rair=287.,Cp=1005.):
    """
    computes the thermodynamic temp(K) as function of the potential temp in a
    hydrostatic atmosphere
    input:
    theta: array of potential temperatures (K)
    z_or_dz: either the altitudes (m) of the profile
             or the differences between levels
    P0: pressure at the lower level (the z[0] level)
        in same units as Pref
    Pref: the 1000 hPa for the conventional definition of potential
          temperature
    g : gravity m/s**2
    Rair: dry air perfect gas constant J/kg/Kelvin
    Cp  : dry air heat capacity J/kg/Kelvin
    """
    RsCp=Rair/Cp
    T0=theta*(P0/Pref)**(RsCp)
    gamma=g/Cp
    psi=1./theta
    try:
        z=np.array(z_or_dz)
        dz=z[...,1:]-z[...,:-1]
        psi[...,1:]=np.cumsum(psi[...,:-1]*dz,axis=-1)
    except:
        dz=np.array(z_or_dz)
        psi[...,1:]=np.cumsum(psi[...,:-1]*dz,axis=-1)
    #
    psi[...,0]=0.
    T=T0-gamma*theta*psi
    return T

def temp2theta(temp,z_or_dz,P0,Pref=1000.,g=9.81,Rair=287.,Cp=1005.):
    """
    computes the potential temperature as function of the thermodynamical
    temperature in a hydrostatic atmosphere
    input:
    temp: array of temperatures (K)
    z_or_dz: either the altitudes (m) of the profile
             or the differences between levels
    P0: pressure at the lower level (the z[0] level)
        in same units as Pref
    Pref: the 1000 hPa for the conventional definition of potential
          temperature
    g : gravity m/s**2
    Rair: dry air perfect gas constant J/kg/Kelvin
    Cp  : dry air heat capacity J/kg/Kelvin
    """
    RsCp=Rair/Cp
    theta0=(Pref/P0)**(RsCp)*temp
    gamma=g/Cp
    #print("g=%f Cp=%f gamma=%f"%(g,Cp,gamma,))
    psi=1./temp
    try:
        z=np.array(z_or_dz)
        dz=z[...,1:]-z[...,:-1] #Note perso : l'usage de l'ellipsis ici ne sert a rien
        psi[...,1:]=np.cumsum(psi[...,:-1]*dz,axis=-1)
    except:  #Note perso : ce try except me parait suspect
        dz=np.array(z_or_dz)
        psi[...,1:]=np.cumsum(psi[...,:-1]*dz,axis=-1)
    psi[...,0]=0.
    #print(psi)
    theta=theta0*np.exp(gamma*psi)
    return theta


def get_keys_for_branch(data, branch, prefix=''):
    result = []
    if isinstance(data, dict):
        for key, value in data.items():
            if key == branch:
                if isinstance(value, dict):
                    for nested_key in value.keys():
                        result.append(f'{prefix}.{key}.{nested_key}')
                else:
                    result.append(f'{prefix}.{key} is not a dictionary.')
                return result
            elif isinstance(value, dict):
                result.extend(get_keys_for_branch(value, branch, f'{prefix}.{key}' if prefix else key))
    return result

def get_child_keys(data, branch):
    keys = branch.split('.')
    current_data = data

    for key in keys:
        if isinstance(current_data, dict) and key in current_data:
            current_data = current_data[key]
        else:
            return []  # Invalid path or key not found

    if isinstance(current_data, dict):
        return list(current_data.keys())
    else:
        return []  # If the final element is not a dictionary

def get_value(data, path):
    keys = path.split('.')
    current_data = data
    for key in keys:
        if key in current_data:
            current_data = current_data[key]
        else:
            return None  # Key not found
    return current_data

class CS_output:
    def __init__(self):
        #
        self.output_folder = None
        #
        self.turbine_output_variables = None
        self.outputs_nc_filename = None
        #
        self.write_fields = None
        self.field_nc_filename = None
        #
        self.field_output_variables = None
        self.zplane_type=None
        self.zmin=None
        self.zmax=None
        self.znumber=None
        self.zcoords = None
        self.xy_sampling_type = None
        self.xbounds=None
        self.ybounds=None
        self.dx=None
        self.dy=None
        #
        self.times = None
        self.run_times = None

class CS_inflow:
    def __init__(self):
        self.data_type = None #'timeseries_hub' or 'timeseries_...' 'sector_occurences'
        #
        #timeseries_hub or timeseries_profiles
        self.wind_velocity = None
        self.wind_dir = None
        #
        #only for timeseries_profiles
        self.u = None
        self.v = None
        self.pottemp = None
        self.tke = []
        self.epsilon = []
        self.heights = None
        #
        #Below only depend on time anyway
        self.times = None
        self.run_times = None #for running particular times only
        self.roughness_height = None
        self.latitude = None
        self.coriolis = True

        #Needed in case of time loop
        self.time_iter = None
        self.damping_lapse = None

        #Run precursor simulation
        self.run_precursor = False
        self.capping_inversion = False
        self.LMO_values = []
        self.ABL_heights = []
        self.lapse_rates = []
        self.dtheta_values = []
        self.dH_values = []
        self.ugeo_values = []

        #Needed for stable profiles
        self.ustar = None
        self.thetastar = None

        #Constants
        self.sea_level_pressure = 101325.
        self.reference_pressure = 1e5

class CS_farm:
    def __init__(self):
        self.turbine_coordinates_xyz = []
        self.turbine_info = []
        self.turbine_types = []
        #
        self.hub_heights = []
        self.rotor_diameters = []
        #
        self.farm_size = 10000.0
        #
        self.cp_curves=[]
        self.ct_curves=[]
        self.cpstar_curves=[]
        self.ctstar_curves=[]

class CS_mesh:
    def __init__(self):
        #
        self.damping_length = 0.0
        self.domain_height = 1000.0
        self.mesh_domain_size = 50000.0
        self.AD_mesh_cell_size = 10.0
        #
        self.remesh = True
        self.mesh_file_name = "mesh.med"

class CS_study:
    def __init__(self,farm_notebook_arg_names=None,farm_notebook_arg_values=None, \
                 prec_notebook_arg_names=None,prec_notebook_arg_values=None, \
                 case_dir=None, result_dir=None, wind_energy_system_file=None, \
                 cs_path=None, salome_path=None, lib_path=None, python_path=None, python_exe=None,\
                 cs_run_folder=None, cs_api_path=None):

        #case
        self.cs_run_folder = cs_run_folder
        self.cs_api_path = cs_api_path
        #
        self.case_dir = case_dir
        self.result_dir = result_dir
        self.case_name = None
        self.cs_path = cs_path
        self.salome_path = salome_path
        self.lib_path = lib_path
        self.python_path = python_path
        self.python_exe = python_exe

        #HPC config
        self.run_node_number = None
        self.run_ntasks_per_node = None
        self.run_wall_time_hours = 10.0
        self.prec_wall_time_hours = 6.0
        self.run_partition = None
        self.mesh_node_number = None
        self.mesh_ntasks_per_node = None
        self.mesh_wall_time_hours = 10.0
        self.mesh_partition = None
        self.wckey = None

        #user input
        self.farm_notebook_arg_names = farm_notebook_arg_names
        self.farm_notebook_arg_values = farm_notebook_arg_values
        self.prec_notebook_arg_names = prec_notebook_arg_names
        self.prec_notebook_arg_values = prec_notebook_arg_values


        #windio configuration
        self.wind_energy_system_file = wind_energy_system_file
        self.wind_system_data = None

        #farm
        self.farm = CS_farm()

        #mesh
        self.mesh = CS_mesh()

        #inflow
        self.inflow = CS_inflow()
        self.wind_origin=270.0

        #output configuration
        self.output = CS_output()

    def set_run_env(self):
        ######CREATE WORK DIRECTORY##########
        #Localize run and copy cs work folders from there to new name
        mkdir(self.cs_run_folder)
        mkdir(self.cs_run_folder + sep + "MESH")
        mkdir(self.cs_run_folder + sep + "logs")
        mkdir(self.cs_run_folder + sep + "meteo_files")
        mkdir(self.cs_run_folder + sep + "precursor_meteo_files")

        #
        os.system("cp -r " + self.cs_api_path + sep + "cs_sources" + sep + "Farm" + " " + self.cs_run_folder+sep+"/.")
        os.system("cp -r " + self.cs_api_path + sep + "cs_sources" + sep + "Precursor" + " " + self.cs_run_folder+sep+"/.")
        ######################################

    def run_case(self, first_case=True, remesh=True, turbine_control=False, damping_layer=False, mesh_file_name="mesh_for_cs_test.med",launch_file_name="launch.sh", job_name="windfarm",log_folder=".", meteo_file_name="",precursor=False, precursor_meteo_file_name=""):
        """
        Function to run code_saturne
        For the moment based on writing a bash file and executing it
        TODO: replace with code_saturne python api
        """
        mesh_nodes=self.mesh_node_number
        mesh_ntasks_per_nodes=self.mesh_ntasks_per_node        
        mesh_wall_time=str(timedelta(hours=self.mesh_wall_time_hours))
        mesh_partition=self.mesh_partition
        cs_nodes=self.run_node_number
        cs_ntasks_per_nodes=self.run_ntasks_per_node
        cs_wall_time=str(timedelta(hours=self.run_wall_time_hours))
        prec_wall_time=str(timedelta(hours=self.prec_wall_time_hours))
        cs_partition=self.run_partition
        wckey=self.wckey

        bash_file_name=self.cs_run_folder + sep + launch_file_name

        if(first_case):
            CS_append='CS_jobids=${CS##* }'
            bash_file=open(bash_file_name,"w")
            #HEADER
            bash_file.write("#!/bin/bash\n"+"#SBATCH --nodes=1\n"+ \
                            "#SBATCH --cpus-per-task=1\n"+ \
                            "#SBATCH --time=00:10:00\n"+ \
                            "#SBATCH --partition="+mesh_partition+"\n"+ \
                            "#SBATCH --wckey="+wckey+"\n"+ \
                            "#SBATCH --output="+log_folder+sep+self.case_name+"_job.out.log\n"+ \
                            "#SBATCH --error="+log_folder+sep+self.case_name+"_job.err.log\n"+ \
                            "#SBATCH --job-name="+self.case_name+"_job"+"\n" + \
                            "\n" + \
                            "module load Miniforge3"+"\n")
        else:
            CS_append='CS_jobids+=":"${CS##* }'
            bash_file=open(bash_file_name,"a")


        #
        split_mesh_folder, split_mesh_file_name= os.path.split(mesh_file_name)
        if(not(remesh)):
            #TODO : raise exception if mesh_file does not exist
            os.system("cp -r " + mesh_file_name+ " " + self.cs_run_folder+sep+"MESH"+sep+".")
        #
        self.mesh.AD_mesh_cell_size=int(np.min(self.farm.rotor_diameters)/8.0)
        #
        self.mesh.mesh_domain_size = np.round(max(self.farm.farm_size*3.2, self.mesh.AD_mesh_cell_size*600),-2) + 2*self.mesh.damping_length


        #
        if('start_rad' in self.farm_notebook_arg_names):
            s_id = self.farm_notebook_arg_names.index('start_rad')
            self.farm_notebook_arg_values[s_id]= 0.5*(self.mesh.mesh_domain_size -2.*self.mesh.damping_length)
        #
        cs_launch_command = self.cs_path +" submit -p setup.xml --case "+self.case_dir+" --notebook-args"
        for k in range(len(self.farm_notebook_arg_names)):
            cs_launch_command += " "+self.farm_notebook_arg_names[k]+"="+str(self.farm_notebook_arg_values[k])
        for k in range(len(self.prec_notebook_arg_names)):
            cs_launch_command += " "+self.prec_notebook_arg_names[k]+"="+str(self.prec_notebook_arg_values[k])

        cs_launch_command += " AD_mesh_cell_size="+str(self.mesh.AD_mesh_cell_size)
        cs_launch_command += " --parametric-args='-m ../"+"MESH"+sep+split_mesh_file_name+"'"

        if(precursor):
            precursor_launch_command = self.cs_path +" submit -p setup.xml --case "+"Precursor"+" --notebook-args"
            for k in range(len(self.prec_notebook_arg_names)):
                precursor_launch_command += " "+self.prec_notebook_arg_names[k]+"="+str(self.prec_notebook_arg_values[k])

            precursor_launch_command += " --kw-args='--meteo ../../../"+precursor_meteo_file_name+"'"
            precursor_launch_command += " --id "+ self.result_dir + " --nodes=1 --partition="+cs_partition+" -J prec_"+job_name+" --time="+prec_wall_time+" --ntasks-per-node=1 --wckey="+wckey+" --exclusive"

        if(meteo_file_name!=""):
            cs_launch_command += " --kw-args='--meteo ../../../"+meteo_file_name+"'"

        cs_launch_command += " --id "+ self.result_dir + " --nodes="+str(cs_nodes)+" --partition="+cs_partition+" -J cs_"+job_name+" --time="+cs_wall_time + " --ntasks-per-node="+str(cs_ntasks_per_nodes) + " --wckey="+wckey+" --exclusive"

        if(remesh and first_case):

            #Original: turbine mesh depending on wind_origin
            #salome_launch_command = self.salome_path + " -t python3 "+self.cs_api_path+sep+"cs_modules"+sep+"csLaunch"+sep+"generate_salome_mesh.py args:--wind_origin="+str(self.wind_origin)+",--disk_mesh_size="+str(self.mesh.AD_mesh_cell_size)+",--domain_size="+str(self.mesh.mesh_domain_size)+",--domain_height="+str(self.mesh.domain_height)+",--output_file='"+mesh_file_name+"'"

            #Non-oriented turbine mesh to avoid multiple remesh in case of timeseries
            #TODO: make dependent on teta variations for single launches. Keyword to force in windio?
            salome_launch_command = self.salome_path + " -t python3 "+self.cs_api_path+sep+"cs_modules"+sep+"csLaunch"+sep+"generate_salome_mesh.py args:--wind_origin="+str("270.0")+",--disk_mesh_size="+str(self.mesh.AD_mesh_cell_size)+",--domain_size="+str(self.mesh.mesh_domain_size)+",--domain_height="+str(self.mesh.domain_height)+",--output_file='MESH"+sep+split_mesh_file_name+"'"
            #
            if(turbine_control):
                salome_launch_command += ",--turbine_control=1.0"
            else:
                salome_launch_command += ",--turbine_control=-1.0"
            #
            if(damping_layer):
                salome_launch_command += ",--damping_layer=1.0"
            else:
                salome_launch_command += ",--damping_layer=-1.0"
            #
            if(precursor):
                precursor_launch_command += " --dependency=afterok:$mesh_jobid\n"
                cs_launch_command += " --dependency=afterok:$meteo_jobid\n"
            else:
                cs_launch_command += " --dependency=afterok:$mesh_jobid\n"

            #
        elif(remesh):
            if(precursor):
                precursor_launch_command += " --dependency=afterok:$mesh_jobid\n"
                cs_launch_command += " --dependency=afterok:$meteo_jobid\n"
            else:
                cs_launch_command += " --dependency=afterok:$mesh_jobid\n"

        elif(precursor):
            cs_launch_command += " --dependency=afterok:$meteo_jobid"

        #============WRITE BASH FILE FOR SLURM================

        #MESH COMMAND
        if(remesh and first_case):
            bash_file.write("#\n"+"# Run Mesh generation with salome and get its jobid\n"+\
                            "mesh_sbatch_output=$(sbatch <<EOF\n"+\
                            "#!/bin/bash\n"+\
                            "#SBATCH --nodes="+str(int(mesh_nodes))+"\n"+\
                            "#SBATCH --cpus-per-task="+str(int(mesh_ntasks_per_nodes))+"\n"+\
                            "#SBATCH --time="+mesh_wall_time+"\n"+\
                            "#SBATCH --partition="+mesh_partition+"\n"+\
                            "#SBATCH --wckey="+wckey+"\n"+\
                            "#SBATCH --output="+log_folder+sep+"mesh_"+job_name+".out.log\n"+\
                            "#SBATCH --error="+log_folder+sep+"mesh_"+job_name+".err.log\n"+\
                            "#SBATCH --job-name=mesh_"+job_name+"\n"+\
                            salome_launch_command+"\n" + "EOF\n"+")\n" + \
                            "# Extract the job ID from the output of sbatch\n" + \
                            "mesh_jobid=${mesh_sbatch_output##* }\n" + \
                            "echo 'meshing job' $mesh_jobid\n")

        #METEO PRECURSOR COMMAND
        Lmoinv_id = self.prec_notebook_arg_names.index('Lmoinv')
        precntmax_id = self.prec_notebook_arg_names.index('precntmax')
        if(precursor):
            bash_file.write("#\n"+"# Run the meteo precursor\n"+\
                            "PRECURSOR=`"+precursor_launch_command+'`\n'+\
                            "precursor_jobid=${PRECURSOR##* }\n"+\
                            "meteo_sbatch=$(sbatch  --parsable --dependency=afterok:$precursor_jobid <<EOF\n"+\
                            "#!/bin/bash"+"\n"+\
                            "#SBATCH --nodes=1"+"\n"+\
                            "#SBATCH --cpus-per-task=1"+"\n"+\
                            "#SBATCH --time=00:10:00"+"\n"+\
                            "#SBATCH --partition=cn"+"\n"+\
                            "#SBATCH --wckey="+wckey+"\n"+\
                            "#SBATCH --output="+log_folder+sep+"meteo_"+job_name+".out.log\n"+ \
                            "#SBATCH --error="+log_folder+sep+"meteo_"+job_name+".err.log\n"+ \
                            "#SBATCH --job-name=meteo_"+job_name+"\n"+\
                            self.python_exe+" "+self.cs_api_path+sep+"cs_modules"+sep+"csMeteo"+sep+"write_meteo_file_from_precursor.py"+ \
                            " --resu_folder='Precursor/RESU/"+self.result_dir + \
                            "' --meteo_filename='"+meteo_file_name + "'"+\
                            " --hub_height="+str(self.farm.hub_heights[0]) + \
                            " --hub_dir="+str(self.wind_origin) + \
                            " --Lmoinv="+str(self.prec_notebook_arg_values[Lmoinv_id]) + \
                            " --precntmax="+str(self.prec_notebook_arg_values[precntmax_id]) + "\n" \
                            "EOF\n"+")\n" + \
                            "# Extract the job ID from the output of sbatch\n" + \
                            "meteo_jobid=${meteo_sbatch##* }\n")

        #CS COMMAND
        #
        bash_file.write("#\n"+"# Run code_saturne after meshing job is finished\n"+\
                        "# with notebooks values:\n"+\
                        "CS=`"+cs_launch_command+'`\n'+\
                        CS_append)

        bash_file.write("\n#\n")
        bash_file.close()
        #===============================================

    def postprocess(self, standalone=False, launch_file_name="launch_postprocess.sh", nodes=1,log_folder="."):
        partition=self.run_partition
        wckey=self.wckey
        ntmax_id = self.farm_notebook_arg_names.index('ntmax')
        ntmax_str=str(self.farm_notebook_arg_values[ntmax_id])
        #
        if not path.exists(self.output.output_folder):
            mkdir(self.output.output_folder)
        os.system("chmod u+rwx "+self.output.output_folder)

        if self.output.zplane_type=="default":
            self.output.zcoords = reduce(lambda re, x: re+[x] if x not in re else re, self.farm.hub_heights, [])
        elif self.output.zplane_type=="grid":
            self.output.zcoords = []
            dz=(self.output.zmax-self.output.zmin)/(self.output.znumber-1)
            for nzi in range(self.output.znumber):
                self.output.zcoords.append(np.round(self.output.zmin + nzi*dz,1))
        postprocess_heights_file = open(self.output.output_folder+"/postprocessing_heights.csv", "w")
        delimiter = ","
        for zcoord in self.output.zcoords[:-1]:
            postprocess_heights_file.write(str(zcoord)+delimiter)
        postprocess_heights_file.write(str(self.output.zcoords[-1]))
        postprocess_heights_file.close()

        postprocess_cases_file = open(self.output.output_folder+"/postprocessing_cases.csv", "w")
        delimiter = ","
        for j in range(len(self.output.run_times)-1):
            case=self.output.run_times[j]
            postprocess_cases_file.write(str(case)+delimiter)
        #
        case=self.output.run_times[-1]
        postprocess_cases_file.write(str(case))
        postprocess_cases_file.close()

        postprocess_fields_file = open(self.output.output_folder+"/postprocessing_fields.csv", "w")
        delimiter = ","
        for field in self.output.field_output_variables[:-1]:
            postprocess_fields_file.write(field+delimiter)
        postprocess_fields_file.write(self.output.field_output_variables[-1])
        postprocess_cases_file.close()

        #
        if(standalone):
            bash_file=open(self.cs_run_folder + sep + launch_file_name,"w")
            bash_file.write("#!/bin/bash\n"+"#SBATCH --nodes=1\n"+ \
                            "#SBATCH --cpus-per-task=1\n"+ \
                            "#SBATCH --time=00:10:00\n"+ \
                            "#SBATCH --partition="+partition+"\n"+ \
                            "#SBATCH --wckey="+wckey+"\n"+ \
                            "#SBATCH --output="+log_folder+sep+"postpro.out.log\n"+ \
                            "#SBATCH --error="+log_folder+sep+"postpro.err.log\n"+ \
                            "#SBATCH --job-name=postpro\n"+ \
                            "#\n")
            bash_file.write("export LD_LIBRARY_PATH="+self.lib_path+" && export PYTHONPATH="+self.python_path+":$PYTHONPATH && "+self.python_exe+" "+self.cs_api_path+sep+"cs_modules"+sep+"csPostpro"+sep+"write_cs_fields.py ../"+self.output.output_folder+" "+self.output.field_nc_filename+" "+self.output.outputs_nc_filename+" "+str(len(self.farm.rotor_diameters))+" "+self.case_name+" "+self.case_dir+" "+ntmax_str)
            bash_file.close()
        else:
            bash_file=open(self.cs_run_folder + sep + launch_file_name,"a")
            bash_file.write("#\n"+"#Postprocessing : writing results in nectdf files\n"+\
                            "postpro_sbatch=$(sbatch --parsable --dependency=afterok:$CS_jobids <<EOF\n"+\
                            "#!/bin/bash\n"+"#SBATCH --nodes=1\n"+ \
                            "#SBATCH --cpus-per-task=1\n"+ \
                            "#SBATCH --time=00:10:00\n"+ \
                            "#SBATCH --partition="+partition+"\n"+ \
                            "#SBATCH --wckey="+wckey+"\n"+ \
                            "#SBATCH --output="+log_folder+sep+"postpro.out.log\n"+ \
                            "#SBATCH --error="+log_folder+sep+"postpro.err.log\n"+ \
                            "#SBATCH --job-name=postpro\n"+ \
                            "export LD_LIBRARY_PATH="+self.lib_path+" && export PYTHONPATH="+self.python_path+":$PYTHONPATH && "+self.python_exe+" "+self.cs_api_path+sep+"cs_modules"+sep+"csPostpro"+sep+"write_cs_fields.py ../"+self.output.output_folder+" "+self.output.field_nc_filename+" "+self.output.outputs_nc_filename+" "+str(len(self.farm.rotor_diameters))+" "+self.case_name+" "+self.case_dir+" "+ntmax_str + "\n" \
                            "EOF\n"+")\n" + \
                            "# Extract the job ID from the output of sbatch\n" + \
                            "postpro_jobid=${postpro_sbatch##* }\n")
            bash_file.close()

    def set_case_dir(self,case_dir):
        """
        change code_saturne case directory folder name
        """
        self.case_dir = case_dir

    def set_result_dir(self,result_dir):
        """
        change code_saturne result directory folder name
        """
        self.result_dir = result_dir

    def set_farm_notebook_arg_names(self,farm_notebook_arg_names):
        """
        change code_saturne notebook_arg_names
        """
        self.farm_notebook_arg_names = farm_notebook_arg_names

    def set_prec_notebook_arg_values(self,prec_notebook_arg_values):
        """
        change code_saturne notebook_arg_values
        """
        self.prec_notebook_arg_values = prec_notebook_arg_values

    def set_notebook_param_from_dictionary(self,farm_notebook_parameters,\
                                           prec_notebook_parameters):
        """
        change code_saturne notebook_arg_values and names from dict
        from both precursor and wind farm simulations
        """
        self.farm_notebook_arg_names = []
        self.farm_notebook_arg_values = []
        for key, value in farm_notebook_parameters.items():
            self.farm_notebook_arg_names.append(key)
            self.farm_notebook_arg_values.append(value)
        #
        self.prec_notebook_arg_names = []
        self.prec_notebook_arg_values = []
        for key, value in prec_notebook_parameters.items():
            self.prec_notebook_arg_names.append(key)
            self.prec_notebook_arg_values.append(value)

    def set_windio(self,wind_energy_system_file):
        """
        change code_saturne notebook_arg_values
        """
        self.wind_energy_system_file = wind_energy_system_file

    def get_windio_data(self):
        #Load the files and store in dictionaries
        self.wind_system_data = load_yaml(self.wind_energy_system_file)
        ####################### LAYOUT and TURBINE DATA ############################
        #
        farm_layout_data = self.wind_system_data['wind_farm']
        turbines_data = farm_layout_data['turbines']

        center_farm=True

        #Turbine coordinates
        turbine_number = len(get_value(self.wind_system_data,'wind_farm.layouts.initial_layout.coordinates.x'))

        #Turbines information. TODO : multiple hub_heights and diameters
        self.farm.rotor_diameters = np.zeros((turbine_number)) + get_value(self.wind_system_data,'wind_farm.turbines.rotor_diameter')
        self.farm.hub_heights =  np.zeros((turbine_number)) + get_value(self.wind_system_data,'wind_farm.turbines.hub_height')

        rotor_center_xy_coordinates = np.zeros((2,turbine_number))
        rotor_center_xy_coordinates[0,:]=get_value(self.wind_system_data,'wind_farm.layouts.initial_layout.coordinates.x')
        rotor_center_xy_coordinates[1,:]=get_value(self.wind_system_data,'wind_farm.layouts.initial_layout.coordinates.y')
        if(center_farm):
            rotor_mean_x=np.mean(rotor_center_xy_coordinates[0,:])
            rotor_center_xy_coordinates[0,:]-=rotor_mean_x
            #
            rotor_mean_y=np.mean(rotor_center_xy_coordinates[1,:])
            rotor_center_xy_coordinates[1,:]-=rotor_mean_y

        farm_xyz = np.zeros((turbine_number,3))
        farm_xyz[:,0] = rotor_center_xy_coordinates[0,:]
        farm_xyz[:,1] = rotor_center_xy_coordinates[1,:]
        farm_xyz[:,2] = self.farm.hub_heights

        self.farm.farm_size = np.sqrt((np.max(farm_xyz[:,0])-np.min(farm_xyz[:,0]))**2+(np.max(farm_xyz[:,1])-np.min(farm_xyz[:,1]))**2) + np.max(self.farm.rotor_diameters)
        self.farm.turbine_coordinates_xyz = farm_xyz
        #
        self.farm.turbine_types = np.zeros((turbine_number))
        self.farm.turbine_types[:] = 1.0 #TODO: multiple turbine_types
        #
        self.farm.turbine_info = np.zeros((farm_xyz.shape[0],5))
        self.farm.turbine_info[:,:3] = farm_xyz
        self.farm.turbine_info[:,3] = self.farm.rotor_diameters
        self.farm.turbine_info[:,4] = self.farm.turbine_types

        #TODO: dev if multiple turbine modes / physical conditions in windIO
        self.farm.ct_curves = []
        self.farm.cp_curves = []
        self.power_curves = []
        performance_info_list = get_child_keys(self.wind_system_data, branch="wind_farm.turbines.performance")
        if('Cp_curve' in performance_info_list):
            cp_curve= np.column_stack((turbines_data['performance']['Cp_curve']['Cp_wind_speeds'],turbines_data['performance']['Cp_curve']['Cp_values']))
            self.farm.cp_curves.append(cp_curve)
        if('power_curve' in performance_info_list):
            power_curve= np.column_stack((turbines_data['performance']['power_curve']['power_wind_speeds'],turbines_data['performance']['power_curve']['power_values']))
            self.power_curves.append(power_curve)
            #
            standard_rho=1.225
            for j in range(len(self.power_curves)):
                cp_curve = np.zeros((len(self.power_curves[j]),2))
                cp_curve[:,0] = self.power_curves[j][:,0] #U0
                for i in range(cp_curve.shape[0]):
                    if(cp_curve[i,0]==0):
                        cp_curve[i,1] = 0
                    else:
                        cp_curve[i,1] = (self.power_curves[j][i,1])/(0.5*standard_rho*(np.pi*((self.farm.rotor_diameters[j]/2.)**2))*(cp_curve[i,0]**3)) #cp

                self.farm.cp_curves.append(cp_curve)

        #
        ct_curve= np.column_stack((turbines_data['performance']['Ct_curve']['Ct_wind_speeds'],turbines_data['performance']['Ct_curve']['Ct_values']))
        self.farm.ct_curves.append(ct_curve)
        #print(get_value(self.wind_system_data, 'wind_farm.turbines.performance.power_curve.power_values'))
        ########################### INFLOW DATA ###################################

        site_data = self.wind_system_data['site']
        resource_data = site_data['energy_resource']
        #====TIMESERIES====
        if('time' in get_child_keys(self.wind_system_data,'site.energy_resource.wind_resource')):
            timeseries_var = get_child_keys(self.wind_system_data,'site.energy_resource.wind_resource')
            
            self.inflow.times = np.array(resource_data['wind_resource']['time'])
            ntimes = len(self.inflow.times)
            
            #TODO: check if we should replace by input from the flow_api in case of binned
            self.inflow.run_times = np.arange(ntimes) #default
            if("cases_run" in get_child_keys(self.wind_system_data, branch="attributes.analysis")):
                if("all_occurences" in get_child_keys(self.wind_system_data, branch="attributes.analysis.cases_run")):
                    all_occurences = get_value(self.wind_system_data,'attributes.analysis.cases_run.all_occurences')
                    if(not(all_occurences)):
                        self.inflow.run_times = get_value(self.wind_system_data,'attributes.analysis.cases_run.occurences_list')
                        self.inflow.times = self.inflow.times[self.inflow.run_times]
                elif("occurences_list" in get_child_keys(self.wind_system_data, branch="attributes.analysis.cases_run")):
                    self.inflow.run_times = get_value(self.wind_system_data,'attributes.analysis.cases_run.occurences_list')
                    self.inflow.times = self.inflow.times[self.inflow.run_times]

            self.inflow.roughness_height = np.array(resource_data['wind_resource']['z0']['data'])
            
            if('lat' in timeseries_var):
                self.inflow.latitude = np.array(resource_data['wind_resource']['lat']['data'])
            else:
                self.inflow.latitude = np.zeros((ntimes)) + 55. #default latitude value

            #====VERTICAL PROFILES====
            if('height' in timeseries_var):
                #PROFILES
                self.inflow.data_type = 'timeseries_profile'
                self.inflow.heights =  np.array(resource_data['wind_resource']['height'])
                nzpoints = len(self.inflow.heights)
                #
                #required
                self.inflow.wind_velocity = np.zeros((nzpoints, ntimes))
                self.inflow.wind_dir = np.zeros((nzpoints, ntimes))
                self.inflow.pottemp = np.zeros((nzpoints, ntimes))
                self.inflow.tke = np.zeros((nzpoints, ntimes))
                self.inflow.epsilon = np.zeros((nzpoints, ntimes))
                
                #Required                
                if('wind_speed' in timeseries_var):
                    self.inflow.wind_velocity[:,:] = np.array(resource_data['wind_resource']['wind_speed']['data']).transpose()
                else:
                    raise ValueError('No "speed" profile in the data')
                if('wind_direction' in timeseries_var):
                    self.inflow.wind_dir[:,:] = np.array(resource_data['wind_resource']['wind_direction']['data']).transpose()
                else:
                    raise ValueError('No "direction" profile in the data')
                if('potential_temperature' in timeseries_var):
                    self.inflow.pottemp[:,:] = np.array(resource_data['wind_resource']['potential_temperature']['data']).transpose()
                else:
                    raise ValueError('No potential_temperature profile in the data')

                #Optional
                if('k' in timeseries_var):
                    self.inflow.tke[:,:] = np.array(resource_data['wind_resource']['k']['data']).transpose()
                else:
                    self.inflow.tke[:,:] = 0.1

                if('epsilon' in timeseries_var):
                    self.inflow.epsilon[:,:] = np.array(resource_data['wind_resource']['epsilon']['data']).transpose()
                else:
                    self.inflow.epsilon[:,:] = 0.003

                #TODO : discriminate stability cases and capping inversion param etc.
                #TODO : list of stabilities for each time occurence
                self.inflow.capping_inversion=True
                self.inflow.run_precursor = True

            #=========================

            #=========HUB HEIGHT======
            else:
                #TIME ONLY
                self.inflow.data_type = 'timeseries_hub'
                self.inflow.run_precursor = False #initialization
                self.inflow.capping_inversion = False #initialization

                #Z profile is generated and tables initialized
                Lz_bl = 2000. ; Nz_bl = 400 ; dz = Lz_bl/Nz_bl; Lz_fa = 26000. ; Nz_fa = 200
                heights_bl = np.linspace(dz/2,Lz_bl-dz/2,Nz_bl,endpoint=True)
                heights_fa = np.linspace(Lz_bl, Lz_fa, Nz_fa, endpoint=True)
                self.inflow.heights = np.concatenate((heights_bl, heights_fa))
                nzpoints = len(self.inflow.heights)
                self.inflow.wind_velocity = np.zeros((nzpoints, ntimes))
                self.inflow.u = np.zeros((nzpoints, ntimes))
                self.inflow.v = np.zeros((nzpoints, ntimes))
                self.inflow.wind_dir = np.zeros((nzpoints, ntimes))
                self.inflow.pottemp = np.zeros((nzpoints, ntimes))
                self.inflow.tke = np.zeros((nzpoints, ntimes))
                self.inflow.epsilon = np.zeros((nzpoints, ntimes))

                #Required
                #TODO: exception if not in file
                self.inflow.wind_velocity =  np.array(resource_data['wind_resource']['wind_speed']['data'])
                self.inflow.wind_dir =  np.array(resource_data['wind_resource']['wind_direction']['data'])
                
                #Optional : stability and capping inversion information
                if('LMO' in timeseries_var):
                    self.inflow.LMO_values = np.array(resource_data['wind_resource']['LMO']['data'])
                else:
                    self.inflow.LMO_values = np.zeros((ntimes)) + 1e10

                #==================Capping inversion========================
                #Activate capping inversion and set default values
                if('ABL_height' in timeseries_var or 'lapse_rate' in timeseries_var \
                   or 'dtheta' in timeseries_var or 'dH' in timeseries_var):                    
                    self.inflow.capping_inversion = True
                    self.inflow.run_precursor = True
                    self.inflow.ABL_heights = np.zeros((ntimes)) + 1500.
                    self.inflow.lapse_rates = np.zeros((ntimes)) + 0.002
                    self.inflow.dtheta_values = np.zeros((ntimes)) + 1.
                    self.inflow.dH_values = np.zeros((ntimes)) + 300.

                #Replace default with user values if any
                if('ABL_height' in timeseries_var):
                    self.inflow.ABL_heights = np.array(resource_data['wind_resource']['ABL_height']['data'])
                if('lapse_rate' in timeseries_var):
                    self.inflow.lapse_rates = np.array(resource_data['wind_resource']['lapse_rate']['data'])
                if('dtheta' in timeseries_var):
                    self.inflow.dtheta_values = np.array(resource_data['wind_resource']['dtheta']['data'])
                if('dH' in timeseries_var):
                    self.inflow.dH_values = np.array(resource_data['wind_resource']['dH']['data'])
                #=========================

            if(not(self.inflow.capping_inversion)):
                self.inflow.coriolis = False
            else:
                self.inflow.coriolis = True

            #print(list(timeseries_dict[1].keys()))
            #print(len(timeseries_dict[1]))
            #print(timeseries_dict[1]['z'])

        #######################  HPC RUNS CONFIG ############################

        self.run_node_number = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.run_node_number')
        self.run_ntasks_per_node = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.run_ntasks_per_node')
        self.run_wall_time_hours = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.run_wall_time_hours')
        self.run_partition = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.run_partition')
        self.mesh_node_number = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.mesh_node_number')
        self.mesh_ntasks_per_node = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.mesh_ntasks_per_node')
        self.mesh_wall_time_hours = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.mesh_wall_time_hours')
        self.mesh_partition = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.mesh_partition')
        self.wckey = get_value(self.wind_system_data,\
                                                      'attributes.analysis.HPC_config.wckey')

        ####################### OUTPUT DATA CONFIG ##########################
        self.output.output_folder = get_value(self.wind_system_data,\
                                                      'attributes.outputs.output_folder')
        self.output.outputs_nc_filename = get_value(self.wind_system_data,\
                                                      'attributes.outputs.turbine_outputs.turbine_nc_filename')
        self.output.turbine_output_variables = get_value(self.wind_system_data,\
                                                      'attributes.outputs.turbine_outputs.output_variables')


        self.output.write_fields = get_value(self.wind_system_data,\
                                                      'attributes.outputs.flow_field.report')
        if(self.output.write_fields):
            self.output.field_nc_filename = get_value(self.wind_system_data,\
                                                      'attributes.outputs.flow_field.flow_nc_filename')
            #
            self.output.field_output_variables = get_value(self.wind_system_data,\
                                                      'attributes.outputs.flow_field.output_variables')
            #
            self.output.zplane_type=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.z_sampling')
            if(self.output.zplane_type == "grid"):
                self.output.zmin=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.z_bounds')[0]
                self.output.zmax=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.z_bounds')[1]
                self.output.znumber=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.z_number')
            elif(self.output.zplane_type == "fixed"):
                self.output.zcoords=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.z_list')
            #
            self.output.xy_sampling_type=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.xy_sampling')
            if(self.output.xy_sampling_type == "bounded"):
                self.output.x_bounds=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.x_bounds')
                self.output.y_bounds=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.y_bounds')
                self.output.dx=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.dx')
                self.output.dy=get_value(self.wind_system_data,'attributes.outputs.flow_field.z_planes.dy')

            ##TODO: process cases where output.run_times are different from inflow.run_times
            self.output.times = self.inflow.times #default
            self.output.run_times = self.inflow.run_times #default
        ########################### CS CONFIG INPUTS ###################################
        self.mesh.remesh = True
        if("remesh" in get_child_keys(self.wind_system_data, branch="attributes.analysis.mesh")):
            self.mesh.remesh=get_value(self.wind_system_data,'attributes.analysis.mesh.remesh')
        if("mesh_file_name" in get_child_keys(self.wind_system_data, branch="attributes.analysis.mesh")):
            self.mesh.mesh_file_name = get_value(self.wind_system_data,'attributes.analysis.mesh.mesh_file_name')

        #TODO: if remesh=False, get mesh domain size and check in case of damping layer

        ################################################################################

    def write_cs_meteo_file(self, cs_meteo_file_name, time_iter, precursor=False):

        sea_level_pressure = self.inflow.sea_level_pressure
        reference_pressure = self.inflow.reference_pressure
        #
        alt = self.inflow.heights
        number_of_alt = len(alt)
        #
        thetaprof = self.inflow.pottemp[:,time_iter]
        T_from_theta = theta2temp(thetaprof, alt, sea_level_pressure, Pref=reference_pressure)
        T = T_from_theta - 273.15
        #
        if precursor:
            uprof = self.inflow.u[:,time_iter]
            vprof = self.inflow.v[:,time_iter]
        else:
            uprof = self.inflow.wind_velocity[:,time_iter]*np.cos(np.pi*(270.0-self.inflow.wind_dir[:,time_iter])/180.0)
            vprof = self.inflow.wind_velocity[:,time_iter]*np.sin(np.pi*(270.0-self.inflow.wind_dir[:,time_iter])/180.0)

        #
        k = self.inflow.tke[:,time_iter]
        epsilon = self.inflow.epsilon[:,time_iter]
        #
        alt = self.inflow.heights
        number_of_alt = len(alt)

        # l'ordre etant: s, x, y, z, epsilon, k, RealTemp, Velocity_x, Velocity_y, Velocity_z
        meteo_file = open(self.cs_run_folder + sep + cs_meteo_file_name, "w")
        meteo_file.write("/ year, quantile, hour, minute, second of the profile:\n")
        meteo_file.write(" 2023,  1, 1,  0,  0\n")
        meteo_file.write("/ location of the meteo profile in the domaine (x,y):\n")
        meteo_file.write(" 0.0  0.0\n")
        meteo_file.write("/ Sea level pressure\n")
        meteo_file.write(" "+str(sea_level_pressure)+"\n")
        meteo_file.write("/ Temperature profile: number of altitudes,(alt.,T in celcius,H in kg/kg ,Nc in N/cm**3)\n")
        meteo_file.write(" "+str(number_of_alt)+"\n")
        delimiter = "  "
        for i in range(number_of_alt):
            meteo_file.write(" "+str(alt[i])+delimiter+str(T[i])+delimiter+str(0.0)+delimiter+str(0.0)+"\n")
        meteo_file.write("/ Wind profile: number of altitudes,(alt.,u,v,ect,eps)\n")
        meteo_file.write(" "+str(number_of_alt)+"\n")
        for i in range(number_of_alt):
            meteo_file.write(" "+str(alt[i])+delimiter+str(uprof[i])+delimiter+str(vprof[i])+delimiter+str(k[i])+delimiter+str(epsilon[i])+"\n")
        #
        meteo_file.close()

    def get_user_data(self, turbine_info, cp_curves, ct_curves):
        center_farm=True
        #Turbine coordinates
        turbine_number = turbine_info.shape[0]
        self.farm.turbine_info = turbine_info
        if(center_farm):
            rotor_mean_x=np.mean(turbine_info[:,0])
            self.farm.turbine_info[:,0]-=rotor_mean_x
            #
            rotor_mean_y=np.mean(turbine_info[:,1])
            self.farm.turbine_info[:,1]-=rotor_mean_y

        self.farm.hub_heights = self.farm.turbine_info[:,2]
        self.farm.rotor_diameters = self.farm.turbine_info[:,3]

        self.farm.farm_size = np.sqrt((np.max(self.farm.turbine_info[:,0])-np.min(self.farm.turbine_info[:,0]))**2+(np.max(self.farm.turbine_info[:,1])-np.min(self.farm.turbine_info[:,1]))**2) + np.max(self.farm.rotor_diameters)
        self.farm.turbine_coordinates_xyz = self.farm.turbine_info[:,:3]
        #
        self.farm.turbine_types[:] = self.farm.turbine_info[:,4]

        #TODO: dev if multiple turbine modes / physical conditions
        self.farm.cp_curves = cp_curves
        self.farm.ct_curves = ct_curves

    def write_cs_input_files(self):
        self.farm.ctstar_curves = []
        self.farm.cpstar_curves = []

        for j in range(len(self.farm.cp_curves)):
            #ct_star
            ct_freestream_speeds = self.farm.ct_curves[j][:,0]
            ct_values = self.farm.ct_curves[j][:,1]
            #
            induct_coeff_a = 0.5*(1.0-np.sqrt(1-ct_values))
            #
            ctstar_values =  ct_values/(1-induct_coeff_a)**2
            ct_disk_speeds = (1-induct_coeff_a)*ct_freestream_speeds

            #cp_star
            cp_freestream_speeds = self.farm.cp_curves[j][:,0]
            cp_values = self.farm.cp_curves[j][:,1]
            ct_values_for_cp_speeds = np.interp(cp_freestream_speeds, ct_freestream_speeds, \
                                                ct_values, left=0, right=0)

            induct_coeff_a_for_cp = 0.5*(1.0-np.sqrt(1-ct_values_for_cp_speeds))
            cpstar_values = cp_values/(1-induct_coeff_a_for_cp)**3
            cp_disk_speeds = (1-induct_coeff_a_for_cp)*cp_freestream_speeds

            cpstar_curve = np.zeros((len(cp_disk_speeds), 2))
            cpstar_curve[:,0] = cp_disk_speeds
            cpstar_curve[:,1] = cpstar_values
            ctstar_curve = np.zeros((len(ct_disk_speeds), 2))
            ctstar_curve[:,0] = ct_disk_speeds
            ctstar_curve[:,1] = ctstar_values

            self.farm.cpstar_curves.append(cpstar_curve)
            self.farm.ctstar_curves.append(ctstar_curve)

        for j in range(len(self.farm.cp_curves)):
            np.savetxt(self.cs_run_folder + sep + self.case_dir + sep +'DATA'+sep+'cp_table_type'+str(j+1)+'.csv', \
                       self.farm.cp_curves[j],delimiter=',')
            np.savetxt(self.cs_run_folder + sep + self.case_dir + sep +'DATA'+sep+'ct_table_type'+str(j+1)+'.csv', \
                       self.farm.ct_curves[j],delimiter=',')
            np.savetxt(self.cs_run_folder + sep + self.case_dir + sep +'DATA'+sep+'cpstar_table_type'+str(j+1)+'.csv', \
                       self.farm.cpstar_curves[j],delimiter=',')
            np.savetxt(self.cs_run_folder + sep + self.case_dir + sep +'DATA'+sep+'ctstar_table_type'+str(j+1)+'.csv', \
                       self.farm.ctstar_curves[j],delimiter=',')

        np.savetxt(self.cs_run_folder + sep + self.case_dir + sep +'DATA'+sep+'turbines_info.csv',\
                   self.farm.turbine_info,delimiter=',',\
                   header="x,y,hub_height,rotor_diameter,turbine_type")

    def set_cs_meteo_file_and_get_wind_dir(self,cs_meteo_file_name, zmeteo):
        os.system("cp "+ cs_meteo_file_name + " " + \
                  self.case_dir+sep+"DATA"+sep+"meteo_file")

        cs_meteo_file = open(self.cs_run_folder + sep + cs_meteo_file_name,"r")
        read_wind_profile = False
        for line in cs_meteo_file:
            if(read_wind_profile and i>-1):
                wind_profile[i,0] = float(line.split()[0])
                wind_profile[i,1] = float(line.split()[1])
                wind_profile[i,2] = float(line.split()[2])
                i+=1
            elif(read_wind_profile and i==-1):
                wind_profile = np.zeros((int(line.split('\n')[0]),3))
                i=i+1
            elif("Wind profile" in line):
                read_wind_profile = True
                i=-1

        windfarm_study.wind_origin = -100.0
        u_hub = np.interp(zmeteo, wind_profile[:,0], wind_profile[:,1])
        v_hub = np.interp(zmeteo, wind_profile[:,0], wind_profile[:,2])
        windfarm_study.wind_origin = np.round(270 - np.arctan2(v_hub,u_hub)*180.0/np.pi,1)
        #print(u_hub,v_hub,windfarm_study.wind_origin)


    def get_wind_dir_from_meteo_file(self,cs_meteo_file_name, zmeteo):
        cs_meteo_file = open(self.cs_run_folder + sep + cs_meteo_file_name,"r")

        read_wind_profile = False
        for line in cs_meteo_file:
            if(read_wind_profile and i>-1):
                wind_profile[i,0] = float(line.split()[0])
                wind_profile[i,1] = float(line.split()[1])
                wind_profile[i,2] = float(line.split()[2])
                i+=1
            elif(read_wind_profile and i==-1):
                wind_profile = np.zeros((int(line.split('\n')[0]),3))
                i=i+1
            elif("Wind profile" in line):
                read_wind_profile = True
                i=-1

        wind_origin = -100.0
        u_hub = np.interp(zmeteo, wind_profile[:,0], wind_profile[:,1])
        v_hub = np.interp(zmeteo, wind_profile[:,0], wind_profile[:,2])
        wind_origin = np.round(270 - np.arctan2(v_hub,u_hub)*180.0/np.pi,1)

        return wind_origin

    def generate_temp_CNBL(self, time_iter):
        ksi = 1.5
        c = 1.0/(2*ksi)
        l = self.inflow.ABL_heights[time_iter] + self.inflow.dH_values[time_iter]/2
        alt = self.inflow.heights
        eta = (alt-l)/(c*self.inflow.dH_values[time_iter])

        f = (np.tanh(eta)+1.0)/2.0
        g = (np.log(2*np.cosh(eta))+eta)/2.0
        for i in np.where(np.isinf(g)):
          g[i] = (np.abs(eta[i])+eta[i])/2.0
        b = self.inflow.lapse_rates[time_iter]*self.inflow.dH_values[time_iter]/3
        a = self.inflow.dtheta_values[time_iter] + b
        th = (293.15 + a*f + b*g)

        return th

    def generate_free_atm_temp_stable(self, time_iter, z, zi, dthdz, th_zi):
        gamma = self.inflow.lapse_rates[time_iter]
        dh = self.inflow.dH_values[time_iter]
        dtheta = self.inflow.dtheta_values[time_iter]
        f = lambda x : (np.tanh(x)+1.0)/2.0
        g = lambda x : (np.abs(x)+x)/2.0 if np.isinf((np.log(2*np.cosh(x))+x)/2.0) else (np.log(2*np.cosh(x))+x)/2.0
        if z<zi:
            return 0
        else:
            ksi = 1.5
            c = 1.0/(2*ksi)
            l = zi + dh/2
            b = gamma*dh/3
            a = dtheta + b
            ztest = np.linspace(zi-100, l, 1000)
            etatest = (ztest-l)/(c*dh)
            derivth = (a/2)*(1-pow(np.tanh(etatest), 2))/(c*dh) + b/2*(1/(c*dh))*(2*np.sinh(etatest))/(2*np.cosh(etatest)) + b/2*(1/(c*dh))
            zoffset = ztest[np.where(derivth<=dthdz)][-1]
            eta_offset =(zoffset-l)/(c*dh)
            l += zi-zoffset
            th_offset = a*f(eta_offset) + b*g(eta_offset)
            eta = (z-l)/(c*dh)
            th = a*f(eta) + b*g(eta)
            return th  + th_zi  - th_offset

    def generate_prof_stable(self, time_iter):
        T0 = 293.15      # Bottom atmospheric temperature (K)
        g = 9.81       # Acceleration of gravity

        # velocity profile parameters
        uref = self.inflow.wind_velocity[time_iter]          # reference wind speed (m/s)
        zref = self.farm.hub_heights[0]

        # Stab Parameters
        Lmo = self.inflow.LMO_values[time_iter]
        lat = self.inflow.latitude[time_iter]
        z0 = self.inflow.roughness_height[time_iter]

        zi = 500.

        alphai = 2

        # vertical grid
        altitude = self.inflow.heights

        kappa = 0.41
        dlmo = 1/Lmo
        omegat = 7.292115e-5
        fcorio = 2*omegat*np.sin(2*np.pi*lat/360)

        if zref>zi:
            ustar = kappa*uref/pow(pow(nwstdt.psim_u_nieuwstadt(zi, zi, dlmo, z0, alphai), 2)\
                                   + pow(nwstdt.psim_v_nieuwstadt(zi, zi, dlmo, z0, alphai), 2.), 0.5)
        else:
            ustar = kappa*uref/pow(pow(nwstdt.psim_u_nieuwstadt(zref, zi, dlmo, z0, alphai), 2)\
                                   + pow(nwstdt.psim_v_nieuwstadt(zref, zi, dlmo, z0, alphai), 2.), 0.5)
        thetastar = pow(ustar, 2)*dlmo/(kappa*g/T0)

        eps = 10
        count = 0
        while eps > 1e-4 and count < 60:
            zi_temp = zi
            zi = 0.4*pow(ustar/(dlmo*fcorio), 0.5)
            if zref>zi:
                ustar = kappa*uref/pow(pow(nwstdt.psim_u_nieuwstadt(zi, zi, dlmo, z0, alphai), 2)\
                                       + pow(nwstdt.psim_v_nieuwstadt(zi, zi, dlmo, z0, alphai), 2.), 0.5)
            else:
                ustar = kappa*uref/pow(pow(nwstdt.psim_u_nieuwstadt(zref, zi, dlmo, z0, alphai), 2)\
                                       + pow(nwstdt.psim_v_nieuwstadt(zref, zi, dlmo, z0, alphai), 2.), 0.5)
            thetastar = pow(ustar, 2)*dlmo/(kappa*9.81/T0)
            eps = np.abs(zi-zi_temp)/zi

        uadim = np.empty(len(altitude))
        vadim = np.empty(len(altitude))
        theta_adim = np.zeros(len(altitude))
        tke_adim = np.zeros(len(altitude))
        epsilon_adim = np.zeros(len(altitude))

        u = np.empty(len(altitude))
        v = np.empty(len(altitude))
        th = np.zeros(len(altitude))

        alpha_prof = np.arctan(nwstdt.psim_v_nieuwstadt(zi, zi, dlmo, z0, alphai)/nwstdt.psim_u_nieuwstadt(zi, zi, dlmo, z0, alphai))
        for i, z in enumerate(altitude):
            if z < zi:
                uadim[i] = nwstdt.psim_u_nieuwstadt(z, zi, dlmo, z0, alphai)/kappa
                vadim[i] = nwstdt.psim_v_nieuwstadt(z, zi, dlmo, z0, alphai)/kappa
                theta_adim[i] = nwstdt.psih_nieuwstadt(z, zi, dlmo, z0, alphai)/kappa
                tke_adim[i] = nwstdt.tke_nieuwstadt(z, zi, dlmo, alphai)
                epsilon_adim[i] = nwstdt.epsilon_nieuwstadt(z, zi, dlmo, alphai)/kappa
            else:
                uadim[i] = nwstdt.psim_u_nieuwstadt(zi, zi, dlmo, z0, alphai)/kappa
                vadim[i] = nwstdt.psim_v_nieuwstadt(zi, zi, dlmo, z0, alphai)/kappa
                # Setting temperature slope based on last points of profiles
                # theta_adim[i] = nwstdt.psih_nieuwstadt(zi, zi, dlmo, z0, alphai)/kappa + (nwstdt.psih_nieuwstadt(zi, zi, dlmo, z0, alphai) - nwstdt.psih_nieuwstadt(zi-1, zi, dlmo, z0, alphai))*(z-zi)/kappa
                tke_adim[i] = 0
                epsilon_adim[i] = 0
            u[i] = ustar*(np.cos(-alpha_prof)*uadim[i]-np.sin(-alpha_prof)*vadim[i])
            v[i] = ustar*(np.sin(-alpha_prof)*uadim[i]+np.cos(-alpha_prof)*vadim[i])
            dthdz_zi = thetastar/kappa*nwstdt.phih_nieuwstadt(zi, zi, dlmo, alphai)
            th_zi = thetastar/kappa*nwstdt.psih_nieuwstadt(zi, zi, dlmo, z0, alphai)
            th[i] = thetastar*theta_adim[i] + T0 + self.generate_free_atm_temp_stable(time_iter, z, zi, dthdz_zi, th_zi)
            tke = pow(ustar, 2)*tke_adim
            tke = np.nan_to_num(tke)
            eps = np.abs(pow(ustar, 3)*epsilon_adim)

        return th, u, v, tke, eps, ustar, thetastar, zi


    def get_lapse_rate_from_meteo_file(self,cs_meteo_file_name):
        # cs_meteo_file = open(cs_meteo_file_name,"r")

        # read_profile = False
        # for line in cs_meteo_file:
        #     if("Wind profile" in line):
        #         read_profile = False
        #     elif(read_profile and i>-1):
        #         temp_profile[i,0] = float(line.split()[0])
        #         temp_profile[i,1] = float(line.split()[1])
        #         i+=1
        #     elif(read_profile and i==-1):
        #         temp_profile = np.zeros((int(line.split('\n')[0]),2))
        #         i=i+1
        #     elif("Temperature profile" in line):
        #         read_profile = True
        #         i=-1

        gamma=-1000
        alt = self.inflow.heights
        #
        theta_prof = self.inflow.pottemp[:, self.inflow.time_iter]
        zlapse = 1500

        slope, intercept = np.polyfit(alt[np.where(alt>zlapse)[0]], theta_prof[np.where(alt>zlapse)[0]], 1)
        gamma=slope
        return gamma
