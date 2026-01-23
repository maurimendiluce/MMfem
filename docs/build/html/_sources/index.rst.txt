.. NavierStokes documentation master file, created by
   sphinx-quickstart on Tue Dec 30 09:17:45 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MMFem documentation!
========================================

Libreria basada en NGSolve para resolver el problema de NavierStokes estacionario y no estacionario.

Descripción general
-------------------

El proyecto está organizado en tres componentes principales:

* ``triangulation.py``: definición de dominios geométricos y etiquetas de borde.
* ``FEMSpaces.py``:------
* ``problems.py``: ----
* ``solvers.py``: ----
* ``adapt.py``: ----
* ``run.py``: script principal que ejecuta una simulación completa.
* ``run_adapt.py``: script principal que ejecuta el método adaptativo

Modelos matemáticos
-------------------

.. math::

   \begin{align*}
      - \nu \Delta \mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} + \nabla p &= \mathbf{f} \qquad \Omega\\
      \nabla \cdot \mathbf{u} &= 0 \qquad \Omega\\
      \mathbf{u} &= \mathbf{g} \qquad \partial\Omega
   \end{align*}
   

Ejecución
---------

.. code-block:: bash

   python run.py

.. code-block:: bash

   python run_adapt.py

Referencias
-----------

* Documentación oficial de NGSolve: https://ngsolve.org/documentation.html


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   triangulation
   FEMSpaces
   problems
   solvers
   adapt


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
