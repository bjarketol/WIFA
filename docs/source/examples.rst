How to simulate
-------------------
For testing, you can clone the repository, add the API to your python path,  copy the  `flow_api/main_api script <https://gitlab.windenergy.dtu.dk/eu-flow/wp4/FLOW_API/-/blob/main/flow_api/main_api.py?ref_type=heads>`_  script in a work directory,  and test it using the `windio_4turbines_* examples <https://gitlab.windenergy.dtu.dk/eu-flow/wp4/FLOW_API/-/tree/main/examples/cases?ref_type=heads>`_.

Then, running the api:
.. code-block:: console
python3  main_api.py    $PATH_TO_examples/windio_4turbines/wind_energy_system/FLOW_toy_study_wind_energy_system.yaml
