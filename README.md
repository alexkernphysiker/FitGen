Genetic algorithm for fitting and solving equations (C++11)
===========================================================

This source can be built with cmake. Just execute two commands

cmake .
make

and you will be able to link the library and run example applications.

paramset.h - set of parameters.

fit_gen.h - base abstract classes used in the whole library.

genetic.h - template classes that provide different types of genetic algorithm. 
Inheriting one template class from another allows to combine several mutation mechanisms.

fitpoints.h - classes that provide optimality criteria for fitting algorithms taking data points into account.
Several template classes allow to extend their functionality.

paramfunc.h - template classes providing parametric functions. 
Mechanism of inheriting template classes allows to construct complicated ones. 
Parameter number is obtained automatically.
There are also template classes that can accept std::function or other function-like objects.

initialcondition.h - classes that provide algorithm of generatins points in fitted parameters space for initialization of population.

filter.h - classes that provide conditions on fitted parameters acting like a filter for new points that appear in the population.

equation.h - template classes and functions providing solving equations using classes declared in other header files.

Example/equation.cpp - example of using genetic algorithm for solving equation.

Example/fit_distribution.cpp - example of fitting of random numbers distribution with using template classes for parametric function and filtering condition that require std::function parameters.

Example/fit_points.cpp - example of fitting data points with x and y errors with using template classes for parametric functions that allow to build complicated ones combining with each other.