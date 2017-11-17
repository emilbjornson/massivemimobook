function R = functionRlocalscatteringApprox(M,theta,ASDdeg,antennaSpacing)
%Generate the spatial correlation matrix for the local scattering model,
%defined in (2.23) with the Gaussian angular distribution. The small-angle
%approximation described in Section 2.6.2 is used to increase efficiency,
%thus this function should only be used for ASDs below 15 degrees.
%
%INPUT:
%M              = Number of antennas
%theta          = Nominal angle
%ASDdeg         = Angular standard deviation around the nominal angle
%                 (measured in degrees)
%antennaSpacing = (Optional) Spacing between antennas (in wavelengths)
%
%OUTPUT:
%R              = M x M spatial correlation matrix
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


%Set the antenna spacing if not specified by input
if  nargin < 4
    
    %Half a wavelength distance
    antennaSpacing = 1/2;
    
end


%Compute the ASD in radians based on input
ASD = ASDdeg*pi/180; 


%The correlation matrix has a Toeplitz structure, so we only need to
%compute the first row of the matrix
firstRow = zeros(M,1);


%Go through all the columns of the first row
for column = 1:M
    
    %Distance from the first antenna
    distance = column-1;
    
    %Compute the approximated integral as in (2.24)
    firstRow(column) = exp(1i*2*pi*antennaSpacing*sin(theta)*distance)*exp(-ASD^2/2 * ( 2*pi*antennaSpacing*cos(theta)*distance )^2);
    
end

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);
