# Sample Code for Modeling ODEs in MATLAB

Christy Pickering, Mac Gabhann Lab

This repository contains the code for the example ordinary differential equation (ODE) mechanistic model with multiple ligands and receptors.

## File Structure

The general structure used is as follows:

* **Species:** lists the species of the model (ligands, antibodies, receptors, and complexes) and the standard numbers used to refer to them, used as a reference
* **Parameters:** lists the parameters used by the model and the standard numbers used to refer to them, used as a reference
* **Equations**: contains the differential equations that govern production, binding and unbinding, transport, internalization, and degradation of the species of the model
* **Main**: calls the ODE solver on the equations file with the parameters and conditions specified by the driver file, also completes a mole balance on the output
* **Driver**: provides the specifics of a given simulation and gives the input to the main file for the ODE solver; determines the species, parameters, initial conditions, and desired output

## Specific Models

The repository contains the following model:

* **Intro**: an example toy model of two ligands binding with two receptors and one co-receptor, used as an example of building and solving a system of ODEs in MATLAB; all parameter values are "dummy" values and not based on particular ligands or receptors
    * "intro_driver.m"

## Included MATLAB Files

* **intro_species.m**: list of species for reference
* **intro_parameters.m**: list of parameters for reference
* **intro_eqns.m**: equations for the binding of all of the receptors and ligands used in this sample model
* **intro_main.m**: calls the solver on the intro_eqns file and performs mole balance, used for all drivers
* **intro_driver.m**: sets up the simulations for binding of molecules in the defined sample system

## Git LFS

Git LFS is used to store the non-text files in this repository. GitHub stores and updates pointers to the image files, and the actual image files are stored on the Git LFS server.

For more information, see the [Git LFS website](https://git-lfs.github.com/).
