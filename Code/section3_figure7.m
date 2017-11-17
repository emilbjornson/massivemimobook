%This Matlab script can be used to reproduce Figure 3.7 in the monograph:
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

%Range of nominal angles of the desired UE
varphiDesiredDegrees = linspace(-180,180,73);
varphiDesiredRadians = varphiDesiredDegrees*(pi/180);

%Range of nominal angles of the interfering UE
varphiInterfererDegrees = linspace(-180,180,73);
varphiInterfererRadians = varphiInterfererDegrees*(pi/180);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the range of the effective SNR in (3.13) for the desired UE
SNR1dB = -10:1:20;
SNR1 = 10.^(SNR1dB/10);

%Define the range of the effective SNR in (3.13) for the interfering UE
SNR2dB = SNR1dB-10;
SNR2 = 10.^(SNR2dB/10);


%Preallocate matrices for storing the simulation results
NMSE_MMSE = zeros(length(SNR1dB),length(varphiInterfererRadians),length(varphiDesiredRadians),length(M));
NMSE_EWMMSE = zeros(length(SNR1dB),length(varphiInterfererRadians),length(varphiDesiredRadians),length(M));
NMSE_LS = zeros(length(SNR1dB),length(varphiInterfererRadians),length(varphiDesiredRadians),length(M));


%% Go through all angles of desired UE
for n1 = 1:length(varphiDesiredRadians)

    %Output simulation progress
    disp([num2str(n1) ' angles out of ' num2str(length(varphiDesiredRadians))]);    
    
    %Compute the spatial correlation matrix of the desired UE
    R1 = functionRlocalscattering(max(M),varphiDesiredRadians(n1),ASDdeg,antennaSpacing);
    
    
    %Go through all angles of interfering UE
    for n2 = 1:length(varphiInterfererRadians)
        
        %Compute the spatial correlation matrix of the interfering UE
        R2 = functionRlocalscattering(M,varphiInterfererRadians(n2),ASDdeg,antennaSpacing);
        
        
        %Go through all number of antennas
        for s = 1:length(SNR1dB)
            
            %Compute the NMSE in (3.20) with the MMSE estimator
            NMSE_MMSE(s,n1,n2) = real(trace(R1 - SNR1(s)*R1*((SNR1(s)*R1+SNR2(s)*R2+eye(M))\R1)))/trace(R1);
            
            
            %Compute the matrix in (3.33) that the received pilot signal is
            %multiplied with in the EW-MMSE estimator
            A_EWMMSE = (sqrt(SNR1(s))/(SNR1(s)+SNR2(s)+1))*eye(M);
            
            %Compute the NMSE in (3.20) with the EW-MMSE estimator
            NMSE_EWMMSE(s,n1,n2) = (trace(R1) + trace(A_EWMMSE*(SNR1(s)*R1+SNR2(s)*R2+eye(M))*A_EWMMSE') - 2*real(trace(A_EWMMSE'*R1))*sqrt(SNR1(s)))/trace(R1);
            

            %Compute the matrix in (3.36) that the received pilot signal is
            %multiplied with in the LS estimator
            A_LS = eye(M)/sqrt(SNR1(s));
            
            %Compute the NMSE in (3.20) with the LS estimator
            NMSE_LS(s,n1,n2) = (trace(R1) + trace(A_LS*(SNR1(s)*R1+SNR2(s)*R2+eye(M))*A_LS') - 2*real(trace(A_LS'*R1))*sqrt(SNR1(s)))/trace(R1);
            
            
        end
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(SNR1dB,mean(mean(NMSE_LS,2),3),'b-.','LineWidth',1);
plot(SNR1dB,mean(mean(NMSE_EWMMSE,2),3),'k-','LineWidth',1);
plot(SNR1dB,mean(mean(NMSE_MMSE,2),3),'r--','LineWidth',1);

xlabel('Effective SNR [dB]');
ylabel('NMSE');
set(gca,'YScale','log');
ylim([1e-2 1e1]);

legend('LS','EW-MMSE','MMSE','Location','SouthWest');
