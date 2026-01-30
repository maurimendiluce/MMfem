Lid-driven cavity problem
=========================

This example demonstrates the solution of the classical lid-driven cavity
problem for the incompressible Navier–Stokes equations using **MMfem**.

The cavity is the unit square Ω = (0,1)², with a moving lid on the top boundary.

----

Problem setup
--------------

- Domain: unit square
- Boundary conditions:
  - Top: horizontal velocity (lid)
  - Other walls: no-slip
- Discretization: Taylor–Hood elements
- Solver: Picard / Newton iterations

----

Python implementation
---------------------

The complete Python script is shown below.

.. literalinclude:: ../../../examples/lid_driven_cavity.py
   :language: python
   :linenos:
