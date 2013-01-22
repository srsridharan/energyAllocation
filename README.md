energyAllocation
================

Code (python) for energy allocation for a collection of houses to a set of aggregators.

This makes use of CPLEX called from Python to solve a set of LP and BINLP problems with convex constraints. 
The aim of the code is to understand the optimality gap between the two approaches (since the BINLP reaches 
its limitations sooner). This would guide the development of scalable approaches  to this problem. 
