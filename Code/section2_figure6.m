%This Matlab script can be used to reproduce Figure 2.6 in the monograph:
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


%Empty workspace and close figures
close all;
clear;


%Number of BS antennas
M = 100;

%Set the angle of the UE
theta = pi/6;

%Set the ASD
ASD = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Compute spatial correlation matrix with local scattering model and
%different angular distributions
R_Gaussian = functionRlocalscattering(M,theta,ASD,antennaSpacing,'Gaussian');
R_Uniform = functionRlocalscattering(M,theta,ASD,antennaSpacing,'Uniform');
R_Laplace = functionRlocalscattering(M,theta,ASD,antennaSpacing,'Laplace');

%Channel correlation matrix with uncorrelated fading
R_uncorrelated = eye(M);

%Extract the eigenvalues and place them in decreasing order
eigenvalues_Gaussian = flipud(eig(R_Gaussian));
eigenvalues_Uniform = flipud(eig(R_Uniform));
eigenvalues_Laplace = flipud(eig(R_Laplace));
eigenvalues_uncorr = flipud(eig(R_uncorrelated));

%Replace negative eigenvalues with a small positive number 
%(since the correlation matrices should be Hermitian)
eigenvalues_Gaussian(eigenvalues_Gaussian<0) = 1e-16;
eigenvalues_Uniform(eigenvalues_Uniform<0) = 1e-16;
eigenvalues_Laplace(eigenvalues_Laplace<0) = 1e-16;


%% Plot the simulation results
figure;
hold on; box on;

plot(1:M,10*log10(eigenvalues_Laplace),'r--','LineWidth',1);
plot(1:M,10*log10(eigenvalues_Uniform),'b-.','LineWidth',1);
plot(1:M,10*log10(eigenvalues_Gaussian),'k','LineWidth',1);
plot(1:M,10*log10(eigenvalues_uncorr),'k:','LineWidth',1);

xlabel('Eigenvalue number in decreasing order');
ylabel('Normalized eigenvalue [dB]');
legend('Laplace','Uniform','Gaussian','Location','SouthEast');
ylim([-50 10]);
