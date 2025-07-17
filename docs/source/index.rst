.. WIFA documentation master file, created by
   sphinx-quickstart on Thu Apr 11 11:03:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: ../img/Logo_FLOW.png
    :align: center
    :width: 150

=====================================
Welcome to WIFA's documentation!
=====================================

WIFA integrates different tools (PyWake, foxes, wayve, code_saturne) of different fidelities to model the atmospheric flow of a wind farm, through a Python framework. This includes engineering wake models, an Atmospheric Perturbation Model (APM), and 3D CFD, all of which allow to estimate the power production of wind farms. These tools are interfaced using a common input and output schema: the WindIO ontology.

The objective of WIFA is to make the use of these models harmonized and easier, although they are initially different both in formulation and use. In particular, this would allow chaining to other tools (e.g. load models).

The presented API is a culmination of close collaboration within the FLOW project. Each supported tool is a product of the institutions represented within the project.

Contents
--------
  .. toctree::
      :maxdepth: 2

      installation

  .. toctree::
      :maxdepth: 2

      examples

  .. toctree::
      :maxdepth: 3

      API

  .. toctree::
      :maxdepth: 3

      references


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


License
-------

MIT License

Copyright (c) 2024 EU FLOW Project

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
