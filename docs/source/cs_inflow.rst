Inflow Boundary Conditions
--------------------------
Inflow refers to the meteorological profiles provided as boundary and initial conditions for the wind farm simulation. It consists of velocity (components :math:`u` and :math:`v`), temperature, turbulent kinetic energy :math:`k` and dissipation rate of turbulent kinetic energy :math:`\varepsilon` profiles. Considering the impact of the methodology used to generate inflow profiles on the model results, the user should be carefull about input data and parameters. Indeed, the inflow condition is determined by the data available as input. A summary of the inflow profile methodology is summarized in the following figure (more details are given in the next sections).

.. figure:: ../img/cs_inflow_Ms4Release.drawio.png

   Schematic representation of the choice of inflow method as a function of input parameters

Input data at hub height or at a reference altitude
+++++++++++++++++++++++++++++++++++++++++++++++++++


Vertical profile as input data
++++++++++++++++++++++++++++++

The Atmospheric Boundary Layer (ABL) is parameterized by setting the system unknowns to a given vertical profile, at the domain boundaries and as initial state. If only the speed and direction at hub height are given, then the atmospheric profiles are extrapolated assuming a neutral stratification, using the Monin Obukhov Similarity theory.
