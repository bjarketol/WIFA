from netCDF4 import Dataset
import sys
from os import path, sep
import numpy as np
import csPostpro.cs_postprocess_utils as cs_pp


def get_output_at_plane_and_time(ens, output_varname, origin, normal,
                                     time_step):
        '''
        evaluate the model for given variable at specific heigh and time.
        '''
        #
        data = cs_pp.read_ensight_data_at_time(ens, inst=time_step)
        plane_cut = cs_pp.extract_plane_from_ensight(data, origin=origin, normal=normal)
        output_variable = cs_pp.get_field_from_ensight(plane_cut, output_varname)

        return output_variable

def get_plane_info(ens, origin, normal):
        '''
        evaluate the model for given variable at specific heigh and time.
        '''
        #
        data = cs_pp.read_ensight_data_at_time(ens, inst=0)
        coupe = cs_pp.extract_plane_from_ensight(data, origin=origin, normal=normal)
        triang = cs_pp.extract_saturne_triangulation(coupe, normal)
        points = data.GetPoints()

        return triang

postpro_dir = sys.argv[1]
file_name = sys.argv[2]
turbine_file_name = sys.argv[3]
turbine_number = int(sys.argv[4])
case_name = sys.argv[5]
case_dir = sys.argv[6]
ntmax_str = sys.argv[7]

if not path.exists(postpro_dir):
    sys.mkdir(postpro_dir)


postprocess_heights_file = open(postpro_dir+"/postprocessing_heights.csv", "r")
postprocess_cases_file = open(postpro_dir+"/postprocessing_cases.csv", "r")
postprocess_fields_file = open(postpro_dir+"/postprocessing_fields.csv", "r")

cases = np.genfromtxt(postprocess_cases_file, delimiter=",")
zplot = np.array(np.genfromtxt(postprocess_heights_file, delimiter=","), ndmin=1)
fields = np.genfromtxt(postprocess_fields_file, delimiter=",", dtype=str)
################TURBINE DATA###################
# File creation
rootgrp = Dataset(postpro_dir+sep+turbine_file_name, "w", format="NETCDF4")
# Dimensions
turbines = rootgrp.createDimension("turbines", turbine_number)
time = rootgrp.createDimension("time", len(cases))
#one = rootgrp.createDimension("one", 1)
coord = rootgrp.createDimension("coord", 2)
# Variables
#total_power_file = rootgrp.createVariable("total_power", "f8", ("one", "time",))
pos_file = rootgrp.createVariable("WT_positions", "f8", ("turbines", "coord", "time",))
z_hub_file = rootgrp.createVariable("WT_hub_heights", "f8", ("turbines", "time",))
diameters_file = rootgrp.createVariable("WT_diameters", "f8", ("turbines", "time",))
#ux_file = rootgrp.createVariable("ux", "f8", ("turbines", "time",))
#uy_file = rootgrp.createVariable("uy", "f8", ("turbines", "time",))
#uz_file = rootgrp.createVariable("uz", "f8", ("turbines", "time",))
u_file = rootgrp.createVariable("velocity", "f8", ("turbines", "time",))
#u_hub_file = rootgrp.createVariable("u_hub", "f8", ("turbines", "time",))
dir_file = rootgrp.createVariable("direction", "f8", ("turbines", "time",))
#ctstar_file = rootgrp.createVariable("ctstar", "f8", ("turbines", "time",))
#cpstar_file = rootgrp.createVariable("cpstar", "f8", ("turbines", "time",))
#thrust_file = rootgrp.createVariable("thrust", "f8", ("turbines", "time",))
powercp_file = rootgrp.createVariable("power", "f8", ("turbines", "time",))

for j, casei in enumerate(cases):
    print(str(j)+'/'+str(len(cases)),end='\r')
    #TODO: string formatting for id "%05d"
    case_name_id = str(int(1000000+casei+1))[1:]
    result_dir=case_name+"_"+case_name_id
    power_file = case_dir + sep + "RESU" + sep + result_dir + sep + "power_iter"+ntmax_str+".csv"
    total_power = []
    x_coords = []
    y_coords = []
    z_hub = []
    diameters = []
    ux = []
    uy = []
    uz = []
    u = []
    u_hub = []
    dir_table = []
    ctstar = []
    cpstar = []
    thrust = []
    power_table_cpstar = []

    # Get turbine and power output info
    with open(power_file, 'r') as file:
        # Read the first line
        first_line = file.readline()
        total_power.append(float(first_line.split(' ')[-1]))
    power_file_table = np.genfromtxt(power_file, delimiter=',', skip_header=2)
    x_coords.append(power_file_table[:, 1] - np.mean(power_file_table[:, 1]))
    y_coords.append(power_file_table[:, 2] - np.mean(power_file_table[:, 2]))
    z_hub.append(power_file_table[:, 3])
    diameters.append(power_file_table[:, 4])
    ux.append(power_file_table[:, 6])
    uy.append(power_file_table[:, 7])
    uz.append(power_file_table[:, 8])
    u.append(power_file_table[:, 9])
    u_hub.append(power_file_table[:, 10])
    dir_table.append(power_file_table[:, 11])
    ctstar.append(power_file_table[:, 13])
    cpstar.append(power_file_table[:, 14])
    thrust.append(power_file_table[:, 15])
    power_table_cpstar.append(power_file_table[:, 16])

    #total_power_file[:, j] = total_power[0]
    pos_file[:, 0, j] = x_coords[0]
    pos_file[:, 1, j] = y_coords[0]
    z_hub_file[:, j] = z_hub[0]
    diameters_file[:, j] = diameters[0]
    #ux_file[:, j] = ux[0]
    #uy_file[:, j] = uy[0]
    #uz_file[:, j] = uz[0]
    u_file[:, j] = u[0]
    #u_hub_file[:, j] = u_hub[0]
    dir_file[:, j] = dir_table[0]
    #ctstar_file[:, j] = ctstar[0]
    #cpstar_file[:, j] = cpstar[0]
    #thrust_file[:, j] = thrust[0]
    powercp_file[:, j] = power_table_cpstar[0]

rootgrp.close()

###############CS FIELD DATA###################
# File creation
rootgrp = Dataset(postpro_dir+sep+file_name, "w", format="NETCDF4")

for i, casei in enumerate(cases):
    print(str(i)+'/'+str(len(cases)),end='\r')
    #TODO: string formatting for id "%05d"
    case_name_id = str(int(1000000+casei+1))[1:]
    result_dir=case_name+"_"+case_name_id
    ens = cs_pp.load_ensight_data(case_dir + sep + "RESU" + sep + result_dir + sep + "postprocessing", "RESULTS_FLUID_DOMAIN.case")
    if i == 0:
        x, y, z, triang=get_plane_info(ens, (0, 0, zplot[0]), (0, 0, 1))
        npoints=triang.triangles.shape[0]
        xtri=np.zeros((npoints))
        ytri=np.zeros((npoints))
        for k in range(npoints):
            i1, i2, i3 = triang.triangles[k]
            xtri[k] = (x[i1] + x[i2] + x[i3])/3.
            ytri[k] = (y[i1] + y[i2] + y[i3])/3.

        #
        points = rootgrp.createDimension("points", npoints)
        nc_cases = rootgrp.createDimension("time", len(cases))
        altitudes = rootgrp.createDimension("altitudes", len(zplot))
        rootgrp.createVariable("x", "f8", ("points",))
        rootgrp.createVariable("y", "f8", ("points",))
        rootgrp.createVariable("z", "f8", ("altitudes",))
        rootgrp.createVariable("time", "f8", ("time",))
        rootgrp.variables["x"][:] = xtri
        rootgrp.variables["y"][:] = ytri
        rootgrp.variables["z"][:] = zplot
        rootgrp.variables["time"][:] = cases
        for field in fields:
            rootgrp.createVariable(field, "f8", ("points", "altitudes", "time",))
    for j, zj in enumerate(zplot):
        zplane_center = (0, 0, zj)

        if "speed" or ("direction" in fields):
            velocity = get_output_at_plane_and_time(ens, "Velocity", zplane_center, (0, 0, 1), 0)
            if "speed" in fields:
                speed = np.sqrt(pow(velocity[:, 0], 2.0) + pow(velocity[:, 1], 2.0))
                rootgrp.variables["speed"][:, j, i] = speed
            if "direction" in fields:
                direction = np.arctan(velocity[:, 1]/velocity[:, 0])*360/(2*np.pi) + 270
                rootgrp.variables["direction"][:, j, i] = direction
        if "pressure" in fields:
            pressure = get_output_at_plane_and_time(ens, "total_pressure", zplane_center, (0, 0, 1), 0)
            rootgrp.variables["pressure"][:, j, i] = pressure
        if "tke" in fields:
            tke = get_output_at_plane_and_time(ens, "k", zplane_center, (0, 0, 1), 0)
            rootgrp.variables["tke"][:, j, i] = tke

rootgrp.close()
