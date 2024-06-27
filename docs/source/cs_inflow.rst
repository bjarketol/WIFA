Inflow Boundary Conditions
--------------------------
Inflow refers to the meteorological profiles provided as boundary and initial conditions for the wind farm simulation. It consists of velocity (components :math:`u` and :math:`v`), temperature, turbulent kinetic energy :math:`k` and dissipation rate of turbulent kinetic energy :math:`\varepsilon` profiles. Inflow conditions are determined by the data given as an input.

.. warning::
   Considering the impact of the methodology used to generate inflow profiles on the model results, the user should be carefull about input data and parameters

A summary of the inflow profile methodology is summarized in the following figure (more details are given in the next sections).

.. figure:: ../img/cs_inflow_Ms4Release.drawio.png

   Schematic representation of the choice of inflow method as a function of input parameters

Input data at hub height or at a reference altitude
+++++++++++++++++++++++++++++++++++++++++++++++++++

If input data contains only values at hub height or at a reference height (absence of ``heights`` dimension in input data), synthetic meteorological profiles are constructed. In the case of minimal input data available, namely ``wind_speed`` and ``wind_direction``, default values for roughness ``z0``, Monin-Obukhov length ``LMO`` and ground temperature ``T0`` are chosen. Default values are shown below.

+------------+-----------------+
| Variable   | Default value   |
+============+=================+
| ``z0``     | :math:`0.0001`  |
+------------+-----------------+
| ``LMO``    | :math:`\infty`  |
+------------+-----------------+
| ``T0``     | :math:`293.15K` |
+------------+-----------------+

In that case, profiles from Monin-Obukhov similarity theory matching the ``wind_speed`` value at hub height are prescribed in a :math:``800m`` high numerical domain with an angle corresponding to the ``wind_direction``, the enrgy equation is turned off and coriolis force is not accounted for. If input data contains ``z0``, ``LMO`` or ``T0``, default values are replaced by input data.

If input data contains at least one variable describing the upper boundary layer and free atmosphere (``ABL_height``, ``capping_inversion_thickness``, ``capping_inversion_strength``, or ``lapse_rate``), the temperature profile above the boundary layer height is constructed based on the parameters provided according to the following figure.

.. figure:: ../img/temp_profile.png
   :scale: 50%

   Potential temperature profile above the boundary layer height :math:`H` (``ABL_height``) with :math:`\Delta h` corresponding to ``capping_inversion_thickness``, :math:`\Delta \theta` to ``capping_inversion_strength`` and :math`\Gamma` to ``lapse_rate``.


If some of the parameters are not present, default values shown below are used to construct the profile.

+---------------------------------+-----------------+
| Variable                        | Default value   |
+=================================+=================+
| ``ABL_height``                  | :math:`1500m`   |
+---------------------------------+-----------------+
| ``capping_inversion_thickness`` | :math:`300m`    |
+---------------------------------+-----------------+
| ``capping_inversion_strength``  | :math:`2K`      |
+---------------------------------+-----------------+
| ``lapse_rate``                  | :math:`1K/km`   |
+---------------------------------+-----------------+

Below the boundary layer height, depending on the Monin-Obukhov length ``LMO``, potential temperature is constant for neutral boundary layer (either no ``LMO`` in input data or ``LMO`` superior to :math:`1000m` and calculated using a blending between Monin-Obukhov similarity theory in the surface layer and Nieuwstadt (1984) :cite:`nieuwstadt1984` theory toward the top of the boundary layer for positive ``LMO`` lower than :math:`1000m`. Negative ``LMO`` corresponding to unstable atmospheric conditions are not supported with this methodology.



Vertical profile as input data
++++++++++++++++++++++++++++++

