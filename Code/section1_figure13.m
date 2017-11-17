%This Matlab script can be used to reproduce Figure 1.13 in the monograph:
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
M = [10 100];

%Select the number of Monte Carlo realizations of the Rayleigh fading
numberOfRealizations = 1000000;

%Generate random UE angles from 0 to 2*pi
varphiDesired = 2*pi*rand(1,numberOfRealizations);
varphiInterfering = 2*pi*rand(1,numberOfRealizations);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Preallocate matrix for storing the simulation results
interferenceGainLoS = zeros(numberOfRealizations,length(M));

%Compute the argument, that appears in (1.28), for different UE angles
argument = (2*pi*antennaSpacing*( sin(varphiDesired(:))  - sin(varphiInterfering(:))) );

%Go through different number of antennas
for mindex = 1:length(M)

    %Compute the g-function in (1.28), but implemented slightly differently
    interferenceGainLoS(:,mindex) = ((1-cos(argument*M(mindex)))./(1-cos(argument)))/M(mindex);

end

CDFvalues_LoS = linspace(0,1,numberOfRealizations);


%Compute the CDF of the relative interference gain for NLoS using the
%Exp(1)-distribution
interferenceGainNLoS = logspace(-6,2,10000);
CDFvalues_NLoS = 1-exp(-interferenceGainNLoS);


%% Plot the simulation results
figure;
hold on; box on;

plot(interferenceGainNLoS,CDFvalues_NLoS,'r--','LineWidth',1);
plot(sort(interferenceGainLoS(:,1)),CDFvalues_LoS,'b-.','LineWidth',1);
plot(sort(interferenceGainLoS(:,2)),CDFvalues_LoS,'k-','LineWidth',1);

set(gca,'Xscale','log');
xlim([1e-6 1e2]),
xlabel('Relative interference gain');
ylabel('CDF');

legend('NLoS, any M','LoS, M=10','LoS, M=100','Location','NorthWest');
