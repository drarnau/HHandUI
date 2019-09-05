[This repository](https://github.com/drarnau/HHandUI) contains all files necessary to replicate the results in _On Households and Unemployment Insurance_ by [Sekyu Choi](http://sekyuchoi.weebly.com/) and [Arnau Valladares-Esteban](http://arnau.eu/).

## Implementation concept
The algorithm receives **parameters** from the outside world, solves the model, and returns **moments of the simulated data** as output.

## Parameters
There are three types of parameters (mutually exclusive categories):
- *Numerical*: Are **not** part of the mathematical representation of the model. Should not play any role in the economic implications of the solution of the model. Should only be relevant for speed and precision. Examples include: grid sizes, tolerance values, etc.
- *Assigned*: Are part of the mathematical representation of the model. Are chosen outside of the calibration exercise.
- *Calibrated*: Are determined by the calibration exercise in which the distance between the moments generated by the model and the moments observed in the data is minimized.

## Mains steps to solve the model
1. Compute **decision rules** using value function iteration.
  - Use Golden search method to find value of assets that solves the Bellman equations.
2. Use **decision rules** to simulate life of a large number of agents and compute targeted moments.

## Compilation instructions
These instructions consider `GNU Fortran (Ubuntu 7.4.0-1ubuntu1~18.04.1) 7.4.0` on a double core machine with 15 GB of RAM memory. Execution time for the above is around 3 minutes. Output of the code is stored in the folder where the executable is saved.

1. Navigate to the folder where the `.f90` files are located.
2. Create an executable file with gfortran:
  - *Normal*: `gfortran -g -fcheck=all -fbacktrace -Wall -mcmodel=large Globals.f90 Utils.f90 Initialisation.f90 VFiteration.f90 Households.f90 Simulation.f90 Main.f90 -o main.out`.
  - Optimised: `gfortran -mcmodel=large Globals.f90 Utils.f90 Initialisation.f90 VFiteration.f90 Households.f90 Simulation.f90 Main.f90 -O3 -o main.out`
3. Execute: `./main`.
