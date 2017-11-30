# ICEP TRAJECTORY SCRIPT
By Allan M. Brooks, Syeda Sabrina, and Kyle J. M. Bishop

## DESCRIPTION

This MATLAB script allows the user to input shape tensors C and D and view the resulting trajectory.


## REQUIREMENTS

This script was written in MATLAB R2016a and has not been compatability tested with previous versions.


## INSTRUCTIONS

  1. Run trajectorySolver.m

  2. Input desired particle symmetry when prompted

  3. The script will output the forms of tensors C and D for the chosen symmetry. When prompted, input nonzero terms for C and D. See Section 4 of the Supporting Information for more details on how each term impacts dynamics.

  4. When prompted, input the time (in dimensionless units) you would like to run the trajectory.

  5. When prompted, input the starting angle of the particle [phi, theta, psi]. See Figure 1 of the main text for more information.

  6. The script will output a figure showing the trajectory based on the given inputs. PLEASE NOTE that the particle geometry shown is a representative of the selected symmetry group but does not correspond to the exact values given by the tensor.

  7. To change the trajectory time and starting angle, follow the prompts.
