Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency 
==================

This code package is distributed along with the monograph:

Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "[Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency](https://www.massivemimobook.com)", Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/2000000093.

For further information about the monograph, please visit: [https://www.massivemimobook.com](https://www.massivemimobook.com)

The package contains a simulation environment, based on Matlab, that can be used to reproduce all the simulation results in monograph. We hope that the code will support the learning of the Massive MIMO topic and also serve as a baseline for further research endeavors. *We encourage you to also perform reproducible research!*


## Abstract of the Monograph

Massive multiple-input multiple-output (MIMO) is one of the most promising technologies for the next generation of wireless communication networks because it has the potential to provide game-changing improvements in spectral efficiency (SE) and energy efficiency (EE). This monograph summarizes many years of research insights in a clear and self-contained way and provides the reader with the necessary knowledge and mathematical tools to carry out independent research in this area. Starting from a rigorous definition of Massive MIMO, the monograph covers the important aspects of channel estimation, SE, EE, hardware efficiency (HE), and various practical deployment considerations.
From the beginning, a very general, yet tractable, canonical system model with spatial channel correlation is introduced. This model is used to realistically assess the SE and EE, and is later extended to also include the impact of hardware impairments. Owing to this rigorous modeling approach, a lot of classic “wisdom” about Massive MIMO, based on too simplistic system models, is shown to be questionable.
The monograph contains many numerical examples, which can be reproduced using Matlab code that is available online.


## Content of Code Package

This code package contains 74 Matlab scripts, 29 Matlab functions, and 7 binary files with Matlab data.

Each script is used to reproduce a particular simulation figure in the monograph. The scripts are named using the convention sectionX_figureY, which is interpreted as the script that reproduces Figure X.Y. A few scripts are instead named sectionX_figureY_Z and will then generate both Figure X.Y and Figure X.Z.

The functions are used by the scripts to carry out certain tasks, such as initiating a simulation setup, generating channel correlation matrices, generating channel realizations, computing channel estimates, computing SEs, computing the power consumption, etc.

The Matlab data files are of the type .mat and contain measurement results or particular simulation results.

See each script and function for further documentation. Note that some of the functions use CVX and QuaDRiGa, which need to be installed separately; see below.


## Software and Hardware Requirements

The code was written to be used in Matlab and has been tested in Matlab 2015b. Some of the scripts and functions might also work in Octave, but there is no guarantee of compatibility.

A few scripts and functions require additional software packages that are developed independently and are delivered with separate licenses. To generate Figures 7.2, 7.41, and 7.42, you need to solve convex optimization problems using CVX from CVX Research, Inc. (http://cvxr.com/cvx/). The code has been tested with CVX 2.1 using the solver Mosek, version 7.1.0.12. We discourage the use of the solvers SDPT3 and SeDuMi since these crashed during the test. To generate Figures 7.41 and 7.42, you also need to generate channels using QuaDRiGa from the Fraunhofer Heinrich Hertz Institute (http://www.quadriga-channel-model.de). The code has been tested with QuaDRiGa version 1.4.8-571.

Since the running example considers a setup with 16 cells, 100 antennas per BS, and 10 UEs per cell, these simulations require a lot of RAM to store the channel correlation matrices and channel realizations. The code has been tested successfully on a MacBook Pro with 8 GB 1600 MHz DDR3 and 2.6 GHz Intel Core i5, which should be viewed as a minimum requirement. Some of the simulations can take days to run, therefore we recommend that you first set nbrOfSetups = 1 to check how much time it takes for each realization of random UE location and shadow fading.


## Acknowledgements

We would like to thank the editor Robert W. Heath Jr. for organizing the review of this monograph and the anonymous reviewers for their constructive and detailed comments. We are grateful for the feedback provided by our proof-readers Alessio Zappone (University of Cassino and Southern Lazio), Maximilian Arnold (University of Stuttgart), Andrea Pizzo (University of Pisa), Daniel Verenzuela, Hei Victor Cheng, Giovanni Interdonato, Marcus Karlsson, Antzela Kosta, Özgecan Özdogan (Linköping university), and Zahid Aslam (Siradel).
Emil Björnson has been supported by ELLIIT, CENIIT, and the Swedish Foundation for Strategic Research.
Luca Sanguinetti has been supported by the ERC Starting Grant 305123 MORE.


## License and Referencing

This code is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our monograph as described above. We also recommend that you mention the existence of this code package in your publication, to spread the word of its existence and to ensure that you will not be accused of plagiarism. 
