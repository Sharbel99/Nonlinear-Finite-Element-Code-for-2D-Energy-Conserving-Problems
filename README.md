# Nonlinear-Finite-Element-Code-for-2D-Energy-Conserving-Problems
A program showcasing the process of solving a 2D boundary value problem problem in a 4-node quadrilateral element using the Newton-Raphson method, using an energy conserving algorithm
Two cases are presented: tension and shear. Various stress and displacement components are plotted in addition to the residual log.

## Description of the problem
The code solves the problem having the following energy function: 
```math
W\left(\textbf{F}\right) = \mu_1 \left(tr \left[\textbf{F}^T \textbf{F}\right] -3 \right) - \mu ln J + \dfrac{\lambda}{2} \left( ln J\right)^2
```

## Instructions
First run the input file (Input_File_axial_strain.m or Input_File_shear_strain.m) then run FEA_Program.m

To edit the code to run a different energy functional, edit the parameters in the input file, and edit the residual and tangent in Elast2d_Elem.m

To change the plots, edit FEA_Program.m

## Reference
Hughes, T. J. (2012). The finite element method: linear static and dynamic finite element analysis. Courier Corporation.
