%This Matlab script can be used to reproduce Figure 3.3 in the monograph:
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

%Range of angular standard deviations in the local scattering model (degrees)
ASDs = [1e-2 1:1:9 10:5:50];

%Define the range of nominal angles of arrival
varphiRange = linspace(-pi,+pi,100);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the effective SNR in (3.13), including processing gain
SNRdB = 10;
SNR = 10.^(SNRdB/10);


%Preallocate matrices for storing the simulation results
NMSE_localscattering = zeros(length(ASDs),length(varphiRange));


%% Go through the range of nominal angles
for r = 1:length(varphiRange)
    
    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);    
    
    %Go through all ASDs
    for n = 1:length(ASDs)
        
        %Compute the spatial correlation matrix
        R = functionRlocalscattering(M,varphiRange(r),ASDs(n),antennaSpacing);
        
        %Compute the NMSE according (3.20)
        NMSE_localscattering(n,r) = real(trace(R - SNR*R*((SNR*R+eye(M))\R))/trace(R));
               
    end
    
end

%Compute the NMSE for the uncorrelated fading case
NMSE_uncorrelated = 1/(SNR+1);


%% Plot the simulation results
figure;
hold on; box on;

plot(ASDs,mean(NMSE_localscattering,2),'k-','LineWidth',1);
plot(ASDs,NMSE_uncorrelated*ones(size(ASDs)),'r--','LineWidth',1);

xlabel('ASD $$(\sigma_{\varphi})$$ [degree]','Interpreter','latex');
ylabel('NMSE');
set(gca,'YScale','log');

legend('Local scattering model','Limit: Uncorrelated','Location','SouthEast');
