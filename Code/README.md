Changelog
==================

This code package will be revised if we find typos and errors. This file keeps track of the changes and the Github version control functionality can be used to identify the exact changes. If you find an error, please contact Emil Bjornson at emil.bjornson@liu.se


## Version 1.05, 2019-04-17

Correcting a few errors in the code comments and adjusting the syntax in a few files to remove Matlab warnings.


## Version 1.04, 2018-10-12

Correcting errors in section7_figure26.m and section7_figure27.m related to the generation of random subarrays. The error changed the horizontal scalings in Figures 7.26 and 7.27, which have been updated in the book, but the qualitative results are the same.


## Version 1.03, 2018-08-01

Revising the definition of the case 20 × 5 × 1 (M = 100) in the case study in Section 7, so that every other antenna has another polarization. This is not a coding error, but makes the simulation setup more practical.


## Version 1.02, 2018-06-08

Correcting "minimize" to "maximize" in functionPowerOptimization_prodSINR.m and cleaning up the code in section7_figure2.m. This error has no impact on the simulation results presented in the book.


## Version 1.01, 2017-11-16

Changing the integration intervals in functionRlocalscattering.m to get accurate values of the small eigenvalues in Figure 2.6.


## Version 1.0, 2017-11-04

First version of the simulation code that reproduces the simulations in:

Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "[Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency](https://www.massivemimobook.com)", Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/2000000093.
