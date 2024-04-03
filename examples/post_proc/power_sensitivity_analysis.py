import numpy as np
from os import sys, path, environ
from os import sep, mkdir, walk
from os import chdir, environ, getcwd
import os as os
import sys
import string, random
from shutil import copy, copytree
#
import netCDF4 as nc
from ot_chaos import *
import openturns as ot
import openturns.viewer as viewer
import scipy as scipy
import xarray as xr

stochastic_input_file=""
stochastic_varnames = ["direction", "speed", "z0"]
MC_sample_size=100 ; number_of_turbines=4 #TODO:read
power_file="turbine_data.nc"

#########Load the data############
#stochastic inputs
input_nc_file = nc.Dataset("windio_toy/plant_energy_resource/Stochastic_atHubHeight.nc")
nvar = len(stochastic_varnames)
MC_sample = np.zeros((len(inputnc_file.variables["time"]),nvar))
for j in range(nvar):
    varname=stochastic_varnames[j]
    MC_sample[:,j] = nc_file.variables[varname][:]

#power outputs
power_data = xr.load_dataset(power_file)
power_table = np.zeros((MC_sample_size, number_of_turbines))
for j in range(number_of_turbines):
    power_table[:,j] = power_data['power'].sel(wt=j).values[:MC_sample_size]


#########Compute Sobol indices############
input_variable_array = MC_sample[:MC_sample_size,:]
power_sobol_indices = np.zeros((number_of_turbines,nvar))
power_total_sobol_indices = np.zeros((number_of_turbines,nvar))
sample_std = np.std(power_table,axis=0)
PC_deg=3
#
copula_type = "independent"  #choices: gaussian, independent
for i in range(number_of_turbines):
    if(number_of_turbines>1):
        std_to_test=sample_std[i]
    else:
        std_to_test=sample_std
    if(std_to_test!=0):        
        if(number_of_turbines>1):
            fullSet = StatisticalSet(input_variable_array[:,:],power_table[:,i])
        else:
            fullSet = StatisticalSet(input_variable_array[:,:],power_table[:])
        training_size = fullSet.setSize
        test_size = 0
        otChaosObject = OTChaos(fullSet, training_size=training_size, \
                                test_size=test_size, prediction_size=0, number_of_sets=1)
        otChaosObject.marginals = []

        for v in range(fullSet.inputDimension):
            otChaosObject.marginals.append(uq_var_marginals[v])
            otChaosObject.copula = copula_type

        otChaosObject.polynomialDegree= PC_deg
        otChaosObject.construct_PCE_ot()
        chaosSI = ot.FunctionalChaosSobolIndices(otChaosObject.polynomialChaosResult) 
        for v in range(fullSet.inputDimension):                    
            if(number_of_turbines>1):
                power_sobol_indices[i,v] =  chaosSI.getSobolIndex(v)  
                power_total_sobol_indices[i,v] =  chaosSI.getSobolTotalIndex(v)  
            else:
                power_sobol_indices[i,v] =  chaosSI.getSobolIndex(v)  
                power_total_sobol_indices[i,v] =  chaosSI.getSobolTotalIndex(v)
