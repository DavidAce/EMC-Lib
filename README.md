# Evolutionary Monte Carlo  (beta)
This program finds the global minimum in a general parameter
space by combining Genetic algorithms and Monte Carlo 
algorithms. 

References:

[On learning strategies for evolutionary Monte Carlo](http://www.people.fas.harvard.edu/~junliu/TechRept/07folder/Goswami%26Liu07.pdf)

[Evolutionary Monte Carlo for protein folding simulations](http://users.phhp.ufl.edu/faliang/papers/2001/JCP2D.pdf)

[Real-parameter evolutionary Monte Carlo with applications to Bayesian mixture models](http://users.phhp.ufl.edu/faliang/papers/2001/RealEMC.pdf)


The program minimizes problems of any dimension and in general performs well despite rough energy landscapes. Problems 1-3 dimensions usually
take under 1 second while 3-6 dimensions may take tens of seconds. Around 10 dimensions takes roughly a minute.

---
## Quick Start

#### From command line
To build project run `build.sh`.

To launch program run `run.sh`

#### From IDE
Some IDE's with CMake support can self-configure from the file CMakeLists.txt found in the project root folder. This
is by far the easiest approach. Recommended: [CLion](https://www.jetbrains.com/clion/download) or [Visual Studio Code](https://code.visualstudio.com/) with C++ and CMake Tools extensions.



#### Requirements
 Please install the following software before building the project.
 * C++ compiler with support for c++14 standard (tested with GNU GCC version >= 5.2).
 * CMake (tested with version >= 3.7).
 * MPI (Preferrably OpenMPI).
 
 The package manager [Hunter](https://github.com/ruslo/hunter) is included to ease the building process.
 During the first build, the dependencies listed in CMakeLists.txt will be downloaded and installed by
 [Hunter](https://github.com/ruslo/hunter) automatically on any platform (Linux/OSX/Win).
 
 The following software is installed by [Hunter](https://github.com/ruslo/hunter):   
 * [Eigen](http://eigen.tuxfamily.org) for tensor support and SVD decomposition.

 The default installation folder in Linux is `~/.hunter`.

---

##### OpenMP
Support for the compiler flag `-fopenmp`, used to parallelize the algorithm.


        
## Usage
See the example_main.cpp file for syntax details. In principle it is enough to include 
the header `source/EMC.h` from your own project to gain access to the class object `objective_function`
which has to be initialized and then passed to the function `minimize()`. The constructor of class object `objective_function`
takes arguments:

        objective_function obj_fun(my_function, lower_bound, upper_bound, tolerance, aux_data1, aux_data2   ... aux_dataN);
 
where `lower_bound` and `upper_bound` are Eigen arrays of type `Eigen::ArrayXd(n)`, where n is the number of parameters you wish to minimize, i.e.,
one lower and upper boundary per parameter. The argument `tolerance` has type `double` and tells the program when
to finish: a small value (say 1e-16) takes longer than a large value (`say 1e-8`) but gives better accuracy. You can
provide as many `aux_data` objects as you wish as long as they are of type `Eigen::ArrayXd` or `Eigen::ArrayXXd`. These
can be used to actually compute your minimizing function and are available as members of the `objective_function` class, like `obj_fun.aux[0]`, `obj_fun.aux[1]`... `obj_fun.aux[N]`.

Lastly, `my_function` is a pointer to a function of return type `long double`, declared as:

        long double my_function(objective_function &obj_fun, Eigen::Array<long double, Dynamic, 1> &inputParameters);

In its definition, you will have to write a function body that maps `inputParameters` (and optionally `obj_fun.aux[]`) into
a real scalar value like an "Fitness" (or "Energy") that you wish to minimize.
        

### The Example program
If you compile and run the program *as is*, it will try to find
the minimum of the paraboloid `f(x,y) = sqrt(x^2 + y^2)` with solution `(0,0)`. 


## Constants
See the file `constants.hpp`. In particular, set the following
variables to your liking: 

`int M` sets the number of parallel subpopulations. Preferrably this
number should be the same as the number cpu-threads available for OpenMP.
(Usually 4 to 8, depending on your cpu).


`int max_generations` controls the maximum duration
of the simulation (default `50000`). 

