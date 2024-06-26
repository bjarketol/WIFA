Installation guidelines
-----------------------


prerequisites
-----------------------
conda environement?



Pywake
~~~~~~~~~~~~~~~~~~~~~~~


FOXES
~~~~~~~~~~~~~~~~~~~~~~~
The installation of *FOXES* is described `here in the documentation <https://fraunhoferiwes.github.io/foxes.docs/installation.html>`_.

For the latest relase, run (for `conda` consult the link above):

.. code-block:: console
  
  pip install foxes

For the latest developments, clone and install the `eu_flow <https://github.com/FraunhoferIWES/foxes/tree/eu_flow>`_
branch from github:

.. code-block:: console
  
  git clone git@github.com:FraunhoferIWES/foxes.git
  cd foxes
  git checkout eu_flow
  pip install -e .


WAYVE
~~~~~~~~~~~~~~~~~~~~~~~

code_saturne and salome
~~~~~~~~~~~~~~~~~~~~~~~
code_saturne and salome should be installed independantly, prior to using the code_saturne FLOW api.

code_saturne, source code and prerequisites for version 8.0 can be found `using this link <https://www.code-saturne.org/cms/web/Download/>`_, namely the `github repository <https://github.com/code-saturne/code_saturne/>`_ with a python script for semi-automated installation with detailed instructions.

For salome, it can be built in two ways:  

  * stand alone `direct download <https://www.salome-platform.org/?page_id=2430/>`_ 
  * building the `salome_cfd extension <https://github.com/code-saturne/salome_cfd_extensions/>`_
