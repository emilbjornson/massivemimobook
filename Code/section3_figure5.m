%This Matlab script can be used to reproduce Figure 3.5 in the monograph:
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

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Nominal angle of desired UE
varphiDesired = pi/6;

%Range of nominal angles of the interfering UE
varphiInterfererDegrees = -180:1:180;
varphiInterfererRadians = varphiInterfererDegrees*(pi/180);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the effective SNR in (3.13) for the desired UE
SNR1dB = 10;
SNR1 = 10.^(SNR1dB/10);

%Define the effective SNR in (3.13) for the interfering UE
SNR2dB = [10 0 -10];
SNR2 = 10.^(SNR2dB/10);


%Preallocate matrices for storing the simulation results
NMSE_corr = zeros(length(varphiInterfererRadians),length(varphiDesired),length(SNR2dB));
NMSE_uncorr = zeros(2,length(SNR2dB));


%Compute the spatial correlation matrix of the desired UE
R1 = functionRlocalscattering(M,varphiDesired,ASDdeg,antennaSpacing);


%% Go through all angles of interfering UE
for n = 1:length(varphiInterfererRadians)
    
    %Output simulation progress
    disp([num2str(n) ' angles out of ' num2str(length(varphiInterfererRadians))]);
    
    %Compute the spatial correlation matrix of the interfering UE
    R2 = functionRlocalscattering(M,varphiInterfererRadians(n),ASDdeg,antennaSpacing);
    
    %Go through all number of antennas
    for s = 1:length(SNR2dB)
        
        %Compute the NMSE according (3.20) when having spatial correlation
        NMSE_corr(n,s) = 1 - SNR1*abs(trace(R1*((SNR1*R1+SNR2(s)*R2+eye(M))\R1)))/trace(R1);
        
    end
    
end


%Go through all number of antennas
for s = 1:length(SNR2dB)
    
    %Compute the NMSE according (3.20) when having uncorrelated fading
    NMSE_uncorr(:,s) = 1 - SNR1/(SNR1+SNR2(s)+1);
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(varphiInterfererDegrees,NMSE_corr(:,1),'k-','LineWidth',1);
plot(varphiInterfererDegrees,NMSE_corr(:,2),'r--','LineWidth',1);
plot(varphiInterfererDegrees,NMSE_corr(:,3),'b-.','LineWidth',1);

plot([varphiInterfererDegrees(1); varphiInterfererDegrees(end)],NMSE_uncorr(:,1),'k-','LineWidth',1);
plot([varphiInterfererDegrees(1); varphiInterfererDegrees(end)],NMSE_uncorr(:,2),'r--','LineWidth',1);
plot([varphiInterfererDegrees(1); varphiInterfererDegrees(end)],NMSE_uncorr(:,3),'b-.','LineWidth',1);

xlabel('Angle of interfering UE [degree]');
ylabel('NMSE');
xlim([-180 180]);
set(gca,'YScale','log');

legend('Same SNR','10 dB weaker','20 dB weaker','Location','SouthWest');
