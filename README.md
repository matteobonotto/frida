
![](docs/FRIDA_logo_v0.PNG)

# FRIDA
Official repo of the FRIDA code.

FRIDA (FRee-boundary Integro-Differential Axisimmetric) couples Finite Element Method and Boundary Element Method (FEM-BEM) to solve the free-boundary axi-symmetric plasma equilibrium problem (i.e., the Grad-Shafranov equation in an open-boundary formulation).

## Compiling mex functions
A list of compatible compilers can be found at these links ([Windows](https://it.mathworks.com/support/requirements/supported-compilers.html), [Linux](https://it.mathworks.com/support/requirements/supported-compilers-linux.html), [Mac](https://it.mathworks.com/support/requirements/supported-compilers-mac.html)). Some functions are written in Fortran 90 and compiled into Matlab-compatible mex functions, and some native Matlab functions are traislated into mex functions using a C++ compiler, thus both compilers need to be installed to properly compile all the routines. However, the instuctions for Fortran and C++ code are kept separate, to ease things for who want to use only one of them.

To compile Fortran 90 code run `this_script.m`

To compile some native matlab functions into C++ mex funcitons run `that_script.m`

## Examples

## Citations
Please cite the following papers:
```
@article{bonotto2022coupled,
  title={A coupled FEM-BEM approach for the solution of the free-boundary axi-symmetric plasma equilibrium problem},
  author={Bonotto, M and Abate, D and Bettini, P and Villone, F and others},
  journal={Commun. Comput. Phys.},
  volume={31},
  number={1},
  pages={27--59},
  year={2022}
}
```
```
@article{bonotto2023efficient,
  title={Efficient Numerical Solution of Coupled Axisymmetric Plasma Equilibrium and Eddy Current Problems},
  author={Bonotto, Matteo and Abate, Domenico and Bettini, Paolo and Iaiunese, Antonio and Isernia, Nicola and Villone, Fabio},
  journal={IEEE Access},
  volume={11},
  pages={27489--27505},
  year={2023},
  publisher={IEEE}
}
```
