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
        self.hub_height = 140.0
        self.rotor_diameter = 240.0 

    def run_case(self, nodes=2, ntasks_per_nodes=48, wall_time_in_minutes=1800):
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
        cs_launch_command = cs_launch_command + " --id "+ self.result_dir + " --nodes="+str(nodes)+" --partition=cn -J " +  self.result_dir +" --time="+str(wall_time_in_minutes) + " --ntasks-per-node="+str(ntasks_per_nodes) + " --wckey=P11YD:CODE_SATURNE --exclusive --dependency=afterok:$mesh_jobid\n" 
        #============WRITE BASH FILE FOR SLURM================
        #TODO :  replace with header from user
        bash_file.write(fr"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=bm
#SBATCH --wckey=P11YD:CODE_SATURNE
#SBATCH --output=cs_windfarm.out.log
#SBATCH --error=cs_windfarm.err.log
#SBATCH --job-name=cs_windfarm
#
# Run Mesh generation with salome and get its jobid
mesh_sbatch_output=$(sbatch <<EOF
#!/bin/bash
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=bm
#SBATCH --wckey=P11YD:CODE_SATURNE
#SBATCH --output=cs_mesh.out.log
#SBATCH --error=cs_mesh.err.log
#SBATCH --job-name=cs_mesh
""")
        bash_file.write(self.salome_path + " -t python3 generate_salome_mesh.py")
        bash_file.write("""
EOF
)

# Extract the job ID from the output of sbatch
mesh_jobid=${mesh_sbatch_output##* }

echo "meshing job" $mesh_jobid
#
# Run code_saturne after meshing job is finished
# with notebooks values:
""")
        bash_file.write(cs_launch_command)
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

    def windio_to_cs_files(self):
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
        self.hub_height =turbines_data['hub_height'] #hauteur moyeu

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
        np.savetxt(self.case_dir + sep +'DATA'+sep+'placement_turbines.csv',farm_xyz,delimiter=',')

        #TODO: dev if multiple turbine modes / physical conditions
        cp_curve= np.column_stack((turbines_data['performance']['Cp_curve']['Cp_wind_speeds'],turbines_data['performance']['Cp_curve']['Cp_values']))
        np.savetxt(self.case_dir + sep +'DATA'+sep+'cp_table.csv',cp_curve,delimiter=',')
        #
        ct_curve= np.column_stack((turbines_data['performance']['Ct_curve']['Ct_wind_speeds'],turbines_data['performance']['Ct_curve']['Ct_values']))
        np.savetxt(self.case_dir + sep +'DATA'+sep+'ct_table.csv',ct_curve,delimiter=',')
    
if __name__ == "__main__":
    """
    Main function of script
    """
    
    #CS case directory
    work_dir = '.'
    os.chdir(work_dir)

    windfarm_study = CS_study(case_dir=work_dir+sep+"Dynamique", \
                              cs_path="/software/rd/saturne/code_saturne/8.0/arch/cronos_impi/bin/code_saturne", \
                              salome_path="/software/rd/salome/logiciels/salome/V9_10_0/salome")

    windfarm_study.set_windio("wind_energy_system"+sep + \
                              "IEA37_case_study_3_wind_energy_system.yaml")    


    windfarm_study.windio_to_cs_files()

    #Simulation parameters. TODO: read from wind resource and user
    cs_notebook_parameters = {
        'Lmoinv': 0.0,
        'z0': 0.7,
        'zref': windfarm_study.hub_height,
        'teta': 270.0,
        'ureff' : 10.0,
        't0' : 287.59,
        'st_method' : 0,
        'isol' : 0,
        'Coriolis' : -1
    }
    windfarm_study.set_notebook_param_from_dictionary(cs_notebook_parameters)

    #Run case
    windfarm_study.set_result_dir("wind_farm_simulation")
    windfarm_study.run_case(nodes=2, ntasks_per_nodes=48, wall_time_in_minutes=40)

    print("Waiting for job to finish")
    print("after meshing is finished, the mesh will be stored as "+windfarm_study.case_dir+sep+"DATA"+sep+"mesh_for_cs_test.med")
    print("after cs simulation is finished, power output will be stored in "+windfarm_study.case_dir+sep+"RESU"+sep+windfarm_study.result_dir+sep+"power.txt file")
