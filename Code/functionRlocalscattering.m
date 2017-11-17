function R = functionRlocalscattering(M,theta,ASDdeg,antennaSpacing,distribution)
%Generate the spatial correlation matrix for the local scattering model,
%defined in (2.23) for different angular distributions.
%
%INPUT:
%M              = Number of antennas
%theta          = Nominal angle
%ASDdeg         = Angular standard deviation around the nominal angle
%                 (measured in degrees)
%antennaSpacing = (Optional) Spacing between antennas (in wavelengths)
%distribution   = (Optional) Choose between 'Gaussian', 'Uniform', and
%                'Laplace' angular distribution. Gaussian is default
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
%This is version 1.1 (Last edited: 2017-11-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Set the antenna spacing if not specified by input
if  nargin < 4
    
    %Half a wavelength distance
    antennaSpacing = 1/2;
    
end

%Set angular distribution to Gaussian if not specified by input
if nargin<5
    distribution = 'Gaussian';
end


%Compute the ASD in radians based on input
ASD = ASDdeg*pi/180;


%The correlation matrix has a Toeplitz structure, so we only need to
%compute the first row of the matrix
firstRow = zeros(M,1);


%Go through all the columns of the first row
for column = 1:M
    
    %Distance from the first antenna
    distance = antennaSpacing*(column-1);
    
    
    %For Gaussian angular distribution
    if strcmp(distribution,'Gaussian')
        
        %Define integrand of (2.23)
        F = @(Delta)exp(1i*2*pi*distance*sin(theta+Delta)).*exp(-Delta.^2/(2*ASD^2))/(sqrt(2*pi)*ASD);
        
        %Compute the integral in (2.23) by including 20 standard deviations
        firstRow(column) = integral(F,-20*ASD,20*ASD);
        
        
    %For uniform angular distribution
    elseif strcmp(distribution,'Uniform')
        
        %Set the upper and lower limit of the uniform distribution
        limits = sqrt(3)*ASD;
        
        %Define integrand of (2.23)
        F = @(Delta)exp(1i*2*pi*distance*sin(theta+Delta))/(2*limits);
        
        %Compute the integral in (2.23) over the entire interval
        firstRow(column) = integral(F,-limits,limits);
        
        
    %For Laplace angular distribution
    elseif strcmp(distribution,'Laplace')
        
        %Set the scale parameter of the Laplace distribution
        LaplaceScale = ASD/sqrt(2);
        
        %Define integrand of (2.23)
        F = @(Delta)exp(1i*2*pi*distance*sin(theta+Delta)).*exp(-abs(Delta)/LaplaceScale)/(2*LaplaceScale);
        
        %Compute the integral in (2.23) by including 20 standard deviations
        firstRow(column) = integral(F,-20*ASD,20*ASD);
        
    end
    
end

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);
