.. schrodinger documentation master file, created by
   sphinx-quickstart on Wed Sep 12 16:13:05 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2

###########
Schrodinger
###########

************
Introduction
************

A program that solves the 1d schrodinger equation. It reads
a file *schrodinger.inp* that specifies a given potential and calculates
energies, wavefunctions, expectation values and standard deviation. All data
is plotted afterwards and saved to *schrodinger.pdf.*


****************
Input and Output
****************

The name of the input file must be called **schrodinger.inp**. By default 
*schrodinger* assumes it is located in the same directory as itself. 
Other directory can be accessed via the *-d* starting option. The input file 
must have the following structure:

::

	4.0           # Mass
	-5.0 5.0 1999 # xMin xMax nPoint
	1 5           # first and last eigenvalue to include in the output
	polynomial    # interpolation type
	3             # # nr. of interpolation points and xy declarations
	-1.0 0.5
	0.0 0.0
	1.0 0.5

You can use as many interpolation points as you want. The number has to fit 
the number of xy declarations afterwards. Possible interpolation types are 
*linear*, *csplines* and *polynomial*. *nPoints* is the number of discrete 
points that are used for calculating.

All calculations are saved in the following files. They are used to plot
the data and can be used for further work:

* *potential.dat*

  interpolated potential in XY-Format:

  ::

     x1 V(x1)
     x2 V(x2)
     :    :

* *energies.dat*

  calculated eigenvalues

  ::

     E1
     E2
     E3
     : 

* *wavefuncs.dat*

  calculated wavefunctions in NXY-Format

  ::

     x1 wf1(x1) wf2(x1) wf3(x1) ...
     x2 wf1(x2) wf2(x2) wf3(x2) ...
     :

* *expvalues.dat*

  expectation values and standard deviation

  ::

     exp_val1 st_dev1
     exp_val2 st_dev2
     :   

****************
Starting Options
****************
Starting options (optional parameters) for the main module *schrodinger.py*.
Use the long form as -*-name* or use the short version showed below.

-------------------
Optional Parameters
-------------------

* directory: -d [path]

  Used to specify the path of the input file *schrodinger.inp*

* split: -s

  Splitting the wavefunctions, expectation values and standard deviations in
  the plot for a better view.

* stretch: -st [float]

  Multiplies the wavefunctions with a factor for a better view.

* markersize: -m [float]

  Changes the markersize of the expectation values and standard deviation.

*******
Modules
*******
--------------
schrodinger.py
--------------
Main module that can be used with the starting options (optional parameters)
above.

-----------------
schrodinger_io.py
-----------------
.. automodule:: schrodinger_io
   :members:

---------------------
schrodinger_solver.py
---------------------
.. automodule:: schrodinger_solver
   :members:

------------------------
schrodinger_visualize.py
------------------------
.. automodule:: schrodinger_visualize
   :members:

--------------------------
test_schrodinger_solver.py
--------------------------
.. automodule:: test_schrodinger_solver
   :members:

****************
Scientific Notes
****************
All calculations are numeric! The potential defined by discrete reference
points inside the input file is interpolated using numerical algorithms 
(e.g., finite differences and integrals as Riemann sums).
The constructions of a tridiagonal matrix allows solving the time independent
schrodinger equation at discrete points as an eigenvalue problem. This 
results in inaccuracies both due to discretization errors and due to rounding
errors in the floating-point number calculation.

The probability of the particle to be inside the given x-boundaries is treated
as 100%, so the probablity outside has to be 0. Therefore the problem has 
aquivalent additonal infinite high potential walls at 'xmin' and 'xmax'. 
In case of non-decreasing wavefunction amplitudes outside the potential
boundaries, reconsideration of the calculated solutions is recommended.
