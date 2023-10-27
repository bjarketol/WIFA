import sys
from os import path, environ
import os as os
#from mpi4py import MPI
import string, random
from os import sep, mkdir, walk
from shutil import copy, copytree
import sys
import numpy as np
from os import path, chdir, environ, getcwd
import yaml
from yaml.loader import SafeLoader
from datetime import datetime

class Loader(yaml.SafeLoader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        super().__init__(stream)

    def include(self, node):

        filename = os.path.join(self._root, self.construct_scalar(node))

        with open(filename, 'r') as f:
            return yaml.load(f, self.__class__)
Loader.add_constructor('!include', Loader.include)

class CS_study:
    def __init__(self, notebook_arg_names=None,notebook_arg_values=None, case_dir=None, result_dir=None, wind_energy_system_file=None, cs_path=None, salome_path=None):

        #case
        self.case_dir = case_dir
        self.result_dir = result_dir
        self.cs_path = cs_path
        self.salome_path = salome_path

        #user input
        self.notebook_arg_names = notebook_arg_names
        self.notebook_arg_values = notebook_arg_values

        #farm configuration
        self.wind_energy_system_file = wind_energy_system_file
        self.turbine_coordinates_xyz = []
        self.hub_height = 100.0
        self.rotor_diameter = 150.0
        self.farm_size = 1000.0
        self.hub_height = 100.0
        self.wind_origin=270.0
        self.cp_curve=[]
        self.ct_curve=[]

    def run_case(self, mesh_nodes=8, mesh_ntasks_per_nodes=1, mesh_wall_time="2-00:00:00", mesh_partition="bm", cs_nodes=2, cs_ntasks_per_nodes=24, cs_wall_time="2-00:00:00", cs_partition="cn", wckey="wckey"):
        """
        Function to run code_saturne
        For the moment based on writing a bash file and executing it
        TODO: replace with code_saturne python api
        """

        bash_file_name="launch_cs_windfarm.sh"
        bash_file=open(bash_file_name,"w")
        #
        cs_launch_command = self.cs_path +" submit -p setup.xml --case "+self.case_dir+" --notebook-args"
        for k in range(len(self.notebook_arg_names)):
            cs_launch_command = cs_launch_command + " "+self.notebook_arg_names[k]+"="+str(round(self.notebook_arg_values[k],2))
        cs_launch_command = cs_launch_command + " AD_mesh_cell_size="+str(int(self.rotor_diameter/27.0))
        cs_launch_command = cs_launch_command + " --id "+ self.result_dir + " --nodes="+str(cs_nodes)+" --partition="+cs_partition+" -J " +  self.result_dir +" --time="+cs_wall_time + " --ntasks-per-node="+str(cs_ntasks_per_nodes) + " --wckey="+wckey+" --exclusive --dependency=afterok:$mesh_jobid\n"


        salome_launch_command = self.salome_path + " -t python3 generate_salome_mesh.py args:--turbine_diameter="+str(self.rotor_diameter)+",--turbine_height="+str(self.hub_height)+",--wind_origin="+str(self.wind_origin)+",--disk_mesh_size="+str(int(self.rotor_diameter/20.0))+",--domain_size="+str(np.round(self.farm_size*5,0))+",--domain_height=1000,--output_file='MESH/mesh_for_cs_test.med'"
        #============WRITE BASH FILE FOR SLURM================
        
        #HEADER 
        bash_file.write("#!/bin/bash\n"+"#SBATCH --nodes=1\n"+ \
                        "#SBATCH --cpus-per-task=1\n"+ \
                        "#SBATCH --time=00:10:00\n"+ \
                        "#SBATCH --partition="+mesh_partition+"\n"+ \
                        "#SBATCH --wckey="+wckey+"\n"+ \
                        "#SBATCH --output=cs_windfarm.out.log\n"+ \
                        "#SBATCH --error=cs_windfarm.err.log\n"+ \
                        "#SBATCH --job-name=cs_windfarm\n")

        #MESH COMMAND
        bash_file.write("#\n"+"# Run Mesh generation with salome and get its jobid\n"+\
                        "mesh_sbatch_output=$(sbatch <<EOF\n"+\
                        "#!/bin/bash\n"+\
                        "#SBATCH --nodes="+str(int(mesh_nodes))+"\n"+\
                        "#SBATCH --cpus-per-task="+str(int(mesh_ntasks_per_nodes))+"\n"+\
                        "#SBATCH --time="+mesh_wall_time+"\n"+\
                        "#SBATCH --partition="+mesh_partition+"\n"+\
                        "#SBATCH --wckey="+wckey+"\n"+\
                        "#SBATCH --output=cs_mesh.out.log\n"+\
                        "#SBATCH --error=cs_mesh.err.log\n"+\
                        "#SBATCH --job-name=cs_mesh\n"+\
                        salome_launch_command+"\n" + "EOF\n"+")\n" + \
                        "# Extract the job ID from the output of sbatch\n" + \
                        "mesh_jobid=${mesh_sbatch_output##* }\n" + \
                        "echo 'meshing job' $mesh_jobid\n")

        #CS COMMAND
        bash_file.write("#\n"+"# Run code_saturne after meshing job is finished\n"+\
                        "# with notebooks values:\n"+\
                        cs_launch_command)
        
        bash_file.close()
        #===============================================

        os.system("sbatch "+ bash_file_name)

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

    def set_notebook_arg_names(self,notebook_arg_names):
        """
        change code_saturne notebook_arg_names
        """
        self.notebook_arg_names = notebook_arg_names

    def set_notebook_arg_values(self,notebook_arg_values):
        """
        change code_saturne notebook_arg_values
        """
        self.notebook_arg_values = notebook_arg_values

    def set_notebook_param_from_dictionary(self,cs_notebook_parameters):
        """
        change code_saturne notebook_arg_values and names from dict
        """
        self.notebook_arg_names = []
        self.notebook_arg_values = []
        for key, value in cs_notebook_parameters.items():
            self.notebook_arg_names.append(key)
            self.notebook_arg_values.append(value)

    def set_windio(self,wind_energy_system_file):
        """
        change code_saturne notebook_arg_values
        """
        self.wind_energy_system_file = wind_energy_system_file

    def get_windio_data(self):
        #Load the files and store in dictionaries
        with open(self.wind_energy_system_file) as f:
            wind_system_data = yaml.load(f, Loader=Loader)
        #
        site_data = wind_system_data['site']
        resource_data = site_data['energy_resource']
        #
        farm_layout_data = wind_system_data['wind_farm']
        turbines_data = farm_layout_data['turbines']

        center_farm=True
        #
        #Turbines information
        self.rotor_diameter = turbines_data['rotor_diameter']
        self.hub_height =turbines_data['hub_height']

        #Turbine coordinates
        turbine_number = len(farm_layout_data['layouts']['initial_layout']['coordinates']['x'])
        rotor_center_xy_coordinates = np.zeros((2,turbine_number))
        rotor_center_xy_coordinates[0,:]=farm_layout_data['layouts']['initial_layout']['coordinates']['x']
        rotor_center_xy_coordinates[1,:]=farm_layout_data['layouts']['initial_layout']['coordinates']['y']
        if(center_farm):
            rotor_mean_x=np.mean(rotor_center_xy_coordinates[0,:])
            rotor_center_xy_coordinates[0,:]-=rotor_mean_x
            #
            rotor_mean_y=np.mean(rotor_center_xy_coordinates[1,:])
            rotor_center_xy_coordinates[1,:]-=rotor_mean_y

        farm_xyz = np.zeros((turbine_number,3))
        farm_xyz[:,0] = rotor_center_xy_coordinates[0,:]
        farm_xyz[:,1] = rotor_center_xy_coordinates[1,:]
        farm_xyz[:,2] = self.hub_height

        self.farm_size = np.sqrt((np.max(farm_xyz[:,0])-np.min(farm_xyz[:,0]))**2+(np.max(farm_xyz[:,1])-np.min(farm_xyz[:,1]))**2) + self.rotor_diameter
        self.turbine_coordinates_xyz = farm_xyz
        
        #TODO: dev if multiple turbine modes / physical conditions
        cp_curve= np.column_stack((turbines_data['performance']['Cp_curve']['Cp_wind_speeds'],turbines_data['performance']['Cp_curve']['Cp_values']))
        self.cp_curve = cp_curve
        #
        ct_curve= np.column_stack((turbines_data['performance']['Ct_curve']['Ct_wind_speeds'],turbines_data['performance']['Ct_curve']['Ct_values']))
        self.ct_curve = ct_curve
        
        #TODO: read from wind resource from windIO
        
    def write_cs_input_files(self):
        np.savetxt(self.case_dir + sep +'DATA'+sep+'cp_table.csv',self.cp_curve,delimiter=',')
        np.savetxt(self.case_dir + sep +'DATA'+sep+'ct_table.csv',self.ct_curve,delimiter=',')
        np.savetxt(self.case_dir + sep +'DATA'+sep+'placement_turbines.csv',\
                   self.turbine_coordinates_xyz,delimiter=',')

        
if __name__ == "__main__":
    """
    Main function of script
    """

    #=================SET STUDY FROM WINDIO =============
    #CS case directory
    work_dir = '.'
    os.chdir(work_dir)

    windfarm_study = CS_study(case_dir=work_dir+sep+"Dynamique", \
                              cs_path="/software/rd/saturne/code_saturne/8.0/arch/cronos_impi/bin/code_saturne", \
                              salome_path="/software/rd/salome/logiciels/salome/V9_10_0/salome")

    windfarm_study.set_windio("windio_inputs" + sep + "wind_energy_system"+sep + \
                              "IEA37_case_study_3_wind_energy_system.yaml")

    windfarm_study.get_windio_data()

    #=================USER MODIFICATIONS===============
    
    #Custom modification, e.g. for sensitivity study
    windfarm_study.wind_origin = 70.0
    windfarm_study.domain_height = 1000.0    
    #windfarm_study.rotor_diameter = 120.0
    #windfarm_study.hub_height = 100
    #windfarm_study.turbine_coordinates_xyz[:,2] = windfarm_study.hub_height
    
    windfarm_study.write_cs_input_files()


    #
    # #Example 1 : inflow is in cs meteo file
    # cs_notebook_parameters = {
    #     'meteo_profile' : 1,  #1 for meteo_file in case_dir/DATA/meteo_example
    #     'teta': windfarm_study.wind_origin,
    #     'WT_d' : windfarm_study.rotor_diameter,
    #     'st_method' : 0, #0 : homogeneous AD
    #     'isol' : 0, #0 if full farm, i>0 if turbine i in isolation
    #     'Coriolis' : 1, #1 if coriolis, 0 if not
    #     'lat' : 55.0, #latitude for coriolis
    #     'energy' : 1 #1 to solve energy equation, 0 for constant density
    # }

    #Example 2 : inflow is in MO profile
    cs_notebook_parameters = {
        'meteo_profile' : 2,
        'Lmoinv': 0.0,
        'z0': 0.1,
        'zref': 10,
        'teta': windfarm_study.wind_origin,
        'ureff' : 10.0,
        't0' : 287.59,
        'WT_d' : windfarm_study.rotor_diameter,
        'st_method' : 0, #0 : homogeneous AD
        'isol' : 0, #0 if full farm, i>0 if turbine i in isolation
        'Coriolis' : 1, #1 if coriolis, 0 if not
        'lat' : 55.0, #latitude for coriolis
        'energy' : 1 #1 to solve energy equation, 0 for constant density
    }
    windfarm_study.set_notebook_param_from_dictionary(cs_notebook_parameters)
    #=================================================

    #Run case
    windfarm_study.set_result_dir("wind_farm_simulation_"+datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    #
    windfarm_study.run_case(mesh_nodes=8, mesh_wall_time="0-00:10:00", mesh_ntasks_per_nodes=1, mesh_partition="bm", cs_nodes=2, cs_ntasks_per_nodes=24, cs_wall_time="0-00:40:00", cs_partition="cn", wckey="P120Y:EFLOW")

    print('Wrote and launched "launch_cs_windfarm.sh". Waiting for job to finish')
    print("After meshing is finished, the mesh will be stored as "+"MESH"+sep+"mesh_for_cs_test.med")
    print("After cs simulation is finished, power output will be stored in "+windfarm_study.case_dir+sep+"RESU"+sep+windfarm_study.result_dir+sep+"power.txt file")
