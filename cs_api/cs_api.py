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
    def __init__(self, notebook_arg_names=None,notebook_arg_values=None, case_dir=None, result_dir=None, wind_energy_system_file=None, cs_path=None):

        #case
        self.case_dir = case_dir
        self.result_dir = result_dir
        self.cs_path = cs_path
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
      
        bash_file_name="launch_code_saturne.sh"
        bash_file=open(bash_file_name,"w")
        code_saturne_path = self.cs_path
        #============WRITE BASH FILE FOR SLURM================
        #TODO :  replace with header from user
        bash_file.write("#!/bin/bash\n")
        bash_file.write("#SBATCH --nodes=1\n")
        bash_file.write("#SBATCH --cpus-per-task=1\n")
        bash_file.write("#SBATCH --time=00:10:00\n")
        bash_file.write("#SBATCH --partition=bm\n")
        bash_file.write("#SBATCH --wckey=P11YD:CODE_SATURNE\n")
        bash_file.write("#SBATCH --output="+"code_saturne.out.log\n")
        bash_file.write("#SBATCH --error="+"code_saturne.err.log\n")
        bash_file.write("#SBATCH --job-name="+"code_saturne"+"\n")
        bash_file.write("#\n")
        bash_file.write("# Change to submission directory\n")
        bash_file.write('if test -n "$SLURM_SUBMIT_DIR" ; then cd $SLURM_SUBMIT_DIR ; fi\n')
        bash_file.write("unset -v SLURM_SUBMIT_DIR\n")
        bash_file.write("\n")
        bash_file.write("# Ensure the correct command is found:\n")
        bash_file.write("export PATH="+code_saturne_path+":$PATH\n")
        bash_file.write("CASE_DIR="+self.case_dir+sep+"DATA/\n")
        bash_file.write("# Get case directory and go inside\n")
        bash_file.write("cd $CASE_DIR\n")
        bash_file.write("# Run command for dynamic case with notebooks values:\n")
        #
        cs_launch_command = "\code_saturne submit -p setup.xml --notebook-args"
        for k in range(len(self.notebook_arg_names)):
            cs_launch_command = cs_launch_command + " "+self.notebook_arg_names[k]+"="+str(round(self.notebook_arg_values[k],2))
        cs_launch_command = cs_launch_command + " --id "+ self.result_dir + " --nodes="+str(nodes)+" --partition=cn -J " +  self.result_dir +" --time="+str(wall_time_in_minutes) + " --ntasks-per-node="+str(ntasks_per_nodes) + " --wckey=P11YD:CODE_SATURNE --exclusive\n"                    
        bash_file.write(cs_launch_command)
        #
        bash_file.close()
        #===============================================
    
        os.system("sbatch "+"launch_full_simulation.sh")
    
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
                        cs_path="/software/rd/saturne/code_saturne/8.0/arch/cronos_impi/bin")

    windfarm_study.set_windio("wind_energy_system"+sep + \
                              "IEA37_case_study_3_wind_energy_system.yaml")    


    windfarm_study.windio_to_cs_files()

    #Simulation parameters. TODO: read from wind resource and user
    Lmoinv=0.0
    z0=0.7
    zref=windfarm_study.hub_height
    teta=270
    ureff = 10.0
    t0 = 287.59
    notebook_arg_names = ["Lmoinv","z0","zref","teta","ureff","t0"]
    notebook_arg_values = [Lmoinv, z0, zref, teta, ureff, t0]
    windfarm_study.set_notebook_arg_names(notebook_arg_names)
    windfarm_study.set_notebook_arg_values(notebook_arg_values)

    #Run case
    windfarm_study.set_result_dir("wind_farm_simulation")
    windfarm_study.run_case(nodes=2, ntasks_per_nodes=48, wall_time_in_minutes=40)
