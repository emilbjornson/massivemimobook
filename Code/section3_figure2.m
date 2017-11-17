%This Matlab script can be used to reproduce Figure 3.2 in the monograph:
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


%Define the range of BS antennas
Mvalues = [1 10 100];

%Angular standard deviation in the local scattering model (in degrees)
ASD = 10;

%Define the range of nominal angles of arrival
varphiRange = linspace(-pi,+pi,100);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the range of effective SNRs in (3.13), including processing gain
SNRdB = -10:1:20;
SNR = 10.^(SNRdB/10);


%Preallocate matrices for storing the simulation results
NMSE = zeros(length(SNR),length(varphiRange),length(Mvalues));


%% Go through the range of nominal angles
for r = 1:length(varphiRange)

    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);    
    
    %Compute the spatial correlation matrix
    R = functionRlocalscattering(max(Mvalues),varphiRange(r),ASD,antennaSpacing);
    
    %Go through all SNRs
    for n = 1:length(SNR)
        
        %Go through all number of antennas
        for mind = 1:length(Mvalues)
            
            %Pick out the correlation matrix for current number of antennas
            Rmatrix = R(1:Mvalues(mind),1:Mvalues(mind));
            
            %Compute the NMSE according (3.20)
            NMSE(n,r,mind) = real(trace(Rmatrix - SNR(n)*Rmatrix*((SNR(n)*Rmatrix+eye(Mvalues(mind)))\Rmatrix)))/trace(Rmatrix);
            
        end
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(SNRdB,mean(NMSE(:,:,1),2),'k-','LineWidth',1);
plot(SNRdB,mean(NMSE(:,:,2),2),'r--','LineWidth',1);
plot(SNRdB,mean(NMSE(:,:,3),2),'b-.','LineWidth',1);

xlabel('Effective SNR [dB]');
ylabel('NMSE');
set(gca,'YScale','log');

legend('M=1','M=10','M=100','Location','SouthWest');
