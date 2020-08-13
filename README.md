# Coupled BEM acoustic

This code provides an implementation for the acoustic scattering problem of two
penetrable scatterers under a BEM philosophy, using planar triangular meshes.
Successfully tested under Julia 1.1.
Details from the theoretical formulation are in the work:
https://www.sciencedirect.com/science/article/abs/pii/S0022460X20304405?via%3Dihub
A preprint version of this work exists at https://arxiv.org/abs/1909.11781

## Required packages

The Julia packages `SpecialFunctions`, `Distributed` and `LinearAlgebra` are required,
so make sure you have these insatalled in your Julia environment.

## Using the code

The code is organized in one main file
```
jul_main.jl
```
which must be loaded previously to any script run. This file loads all the necessary
files (which resides in the `src` directory).
In this file the user can configure its parallel environment (i.e. how much processors
he has in the machine running Julia). The line `n_cores` define how many cores are 
added to the main processor. For example, in a machine with 4 cores `n_cores=3` is the
maximum allowable that makes sense. I recommend to use as many cores as you can.

Once loaded the main file, by typing
```
include("jul_main.jl")
```
you are ready to execute the two-spheres example provided (backscattering from two
penetrable spheres).
```
include("Script_Shells_bem_K.jl")
```





## Memory considerations




Thank you for try the code.
Any suggestions and comments will be welcomed at sivasadartantasvueltasnosirve@gmail.com.
