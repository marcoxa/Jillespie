# Jillespie

A [Julia](https://www.julialang.org) implemention of *Gillespie's
Stochastic Simulation* procedure(s) for biochemical systems.

The implementation is essentially an exploration of Julia's
capabilities; it tries to be "efficient" and "usable", by defining a
*Domain Specific Language* based on macros.


## Installation

Eventually the code will be registered in Julia's ecosystem.  For the
time being, once you have downloaded/cloned the code in a folder of
your choice on your machine, you can just use the normal Julia's
functionalities to use the system.


## Standard Disclaimers about Stochastic Gillespie/Montecarlo Simulations 

As noted in "Numerical Recipes in C++" (reference), stochastic
simulations (hence Gillespie style ones) are a great way to burn
computer cycles.  Of course, Julia is not the Fortran IV used in
Gillespie's original article (Fortran is more efficient), and your
watch has more computing power than the IBM 360 Gillespie's wrote his
code on, and yet you should be careful, especially because the current
implementation does not implement tau-leaping (yet).

Another caveat is that Gillespie was keenly aware of the constraints
under which his simulations produce a decent approximation of a trace
consistent with a solution of the *Chemical Master Equation*.  Using
the algorithm to model intra-cellular phenomena introduces, per se, an
extra layer of ... "approximation".


### A NOTE ON FORKING

Of course you are free to fork the project subject to the current
licensing scheme.  However, before you do so, I ask you to consider
plain old "cooperation" by asking me to become a developer.
It helps keeping the entropy level at an acceptable level.


Enjoy

Marco Antoniotti, Milan, Italy, (C) 2023
