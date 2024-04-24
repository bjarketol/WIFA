#!/usr/bin/env python3
import sys
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy as vtk_to_np

def load_ensight_data(path,case):
    """
    load the ensight file
    """
    ens = vtk.vtkGenericEnSightReader()
    ens.SetFilePath(path)
    ens.SetCaseFileName(case)
    ens.ReadAllVariablesOn()
    ens.Update()
    return ens

def read_ensight_data_at_time(ens,inst=None):
    """
    read the results for a given time inst
    """
    try:
        times = ens.GetTimeSets().GetItem(0)
        if (inst>times.GetSize()-1 or inst<0 or inst is None):
            print("Time Exceeds max or is negative. Read time %i"%(times.GetSize()-1))
            inst=times.GetSize()-1
        ens.SetTimeValue(times.GetTuple1(inst))
    except AttributeError:
        print("Error reading the required time step")
        pass
    ens.Update()
    data=ens.GetOutput()
    if data.GetClassName()=='vtkMultiBlockDataSet':
        data=data.GetBlock(0)

    return data

def ensight_time_list(ens):
    """
    read the time list
    """
    try:
        times=ens.GetTimeSets().GetItem(0)
        return list(range(times.GetNumberOfTuples())),[times.GetTuple1(i) for i in range(times.GetNumberOfTuples())]
    except AttributeError:
        return []

def ensight_var_list(data):
    """
    read the list of variable names
    """
    fields=getattr(data,"GetCellData")()
    varname_list=[fields.GetArrayName(i) for i in range(fields.GetNumberOfArrays())]
    return varname_list

def get_field_from_ensight(data,varname):
    """
    extract the field 'varname' from ensight 'data'
    """
    field=getattr(data,"GetCellData")().GetArray(varname)
    if field  is not None:
        tab=vtk_to_np(field)
        tab=tab.copy()
        return tab
    else:
        return None

def get_point_coords_from_ensight(data):
    """
    get the coordinates of mesh points
    """
    return np.array([data.GetPoint(i) for i in range(data.GetNumberOfPoints())]).transpose()

def extract_plane_from_ensight(data,origin=(0.,0.,0.),normal=(0.,0.,1.)):
    """
    extract plane of center 'origin' orthogonally to 'normal'
    """
    plane = vtk.vtkPlane()
    plane.SetOrigin(*origin)
    plane.SetNormal(*normal)
    extract = vtk.vtkCutter()
    extract.SetInputData(data)
    extract.SetCutFunction(plane)
    extract.Update()
    plane_cut=extract.GetOutput()

    if not plane_cut.GetNumberOfCells():
        print("Empty plane - PASS")
        return None

    plane_cut.GetCellData().SetActiveScalars("Velocity")

    return plane_cut
