function a = functionSpatialSignature3DLoS(U,varphi,theta,lambda)
%Compute the spatial signature using the 3D LoS model defined in Section
%7.3.1.
%
%INPUT:
%U        = 3 x M matrix containing the positions of the M antennas
%varphi   = Azimuth angle in radians
%theta    = Elevation angle in radians
%
%OUTPUT:
%a        = M x 1 spatial signature vector
%
%
%This Matlab function was developed to generate simulation results to:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), 
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", 
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, 
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.0 (Last edited: 2017-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Define the wave vector in (7.11)
k = -2*pi/lambda * [cos(varphi)*cos(theta); sin(varphi)*cos(theta); sin(theta)];

%Compute the spatial signature as in (7.12)
a = transpose(exp(1i* k'*U));
