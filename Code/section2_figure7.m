%This Matlab script can be used to reproduce Figure 2.7 in the monograph:
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


%Define the range of number of BS antennas
M = 1:200;
Mmax = max(M);

%Set the nominal angle of the UE
theta = pi/6;

%Set the range of ASDs
ASDs = [10 30];

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Preallocate matrix for storing the simulation results
variance = zeros(length(M),length(ASDs));


%% Go through the range of ASDs
for n = 1:length(ASDs)

    %Output simulation progress
    disp([num2str(n) ' ASDs out of ' num2str(length(ASDs))]);    
    
    %Compute spatial correlation matrix with the local scattering model for
    %the maximum number of antennas
    R = functionRlocalscattering(Mmax,theta,ASDs(n),antennaSpacing,'Gaussian');
    
    %Go through the range of antennas
    for m = 1:length(M)
        
        %Extract the correlation matrix of right size
        Rm = R(1:M(m),1:M(m));

        %Compute the variance of the channel hardening according to (2.17)
        variance(m,n) = real(trace(Rm*Rm)/(trace(Rm)).^2);
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(M,variance(:,1),'r--','LineWidth',1);
plot(M,variance(:,2),'b-.','LineWidth',1);
plot(M,1./M,'k-','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Variance in (2.17)');
legend('Gaussian, ASD 10^o','Gaussian, ASD 30^o','Uncorrelated','Location','NorthEast');
ylim([0 0.2]);
