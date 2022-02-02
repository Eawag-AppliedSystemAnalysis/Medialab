# Medialab
Modeling Early DIAgenesis using MATLAB

MEDIALAB (Modeling Early DIAgenesis using MATLAB) is an early diagenesis model, which calculates concentrations and fluxes of chemical species as well as rates of all the biogeochemical pathways at each depth of aquatic sediments for a specific time period. It is based on the model MATSEDLAB, which had been created by B. Shafei as part of his Ph.D. thesis project at the Georgia Institute of Technology under the supervision of P. Van Cappellen (Shafei, 2012). The code was further improved and subsequently applied to simulate early diagenesis in a lake with variable boundary conditions by Steinsberger et al. (2019). 

A system of partial differential equations corresponding to early diagenesis equations are automatically generated through MATLAB’s symbolic programming capabilities and solved using MATLAB’s built-in solver pdepe to evaluate the temporal and spatial distribution of chemical species.

![Output Image](https://github.com/Eawag-AppliedSystemAnalysis/Medialab/blob/master/docs/example_output.png)

## Documentation

Read the docs [here.](https://github.com/Eawag-AppliedSystemAnalysis/Medialab/blob/master/docs/MEDIALAB%20User%20Manual_220201.pdf)

## Install 

- Clone the repository to your local machine using the command: 

 `git clone https://github.com/Eawag-AppliedSystemAnalysis/Medialab.git`
 
 Note that the repository will be copied to your current working directory.
 
 #### Requirements
 
 The execution of MEDIALAB requires an active installation of MATLAB (Version 7.6 release R2008a or later). It is recommended to allow at least 1GB of space on the hard drive for the model output and 1GB of contiguous random-access memory for the initialization routine.

## Quick Start

1. Open the repository in MatLab.
2. Edit the values in the input files (located in the inputs folders) according to your desired scenario. Details can be found in the docs. 
3. Run mainMEDIALAB.m in order to do the computation and produce the output files.
4. Run postprocessMEDIALAB.m in order to produce output figures from the data (for details consult the docs)