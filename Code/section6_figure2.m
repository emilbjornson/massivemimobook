%This Matlab script can be used to reproduce Figure 6.2 in the monograph:
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

%Define the range of nominal angles of arrival
varphiRange = linspace(-pi,+pi,100);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the range of effective SNRs in (3.13), including processing gain
SNRdB = -10:1:30;
SNR = 10.^(SNRdB/10);

%Select the range of hardware qualities. The values at the same position in
%the different vectors are considered simultaneously.
kappatRange = [0.95 0.99 1];
kapparRange = [0.95 0.99 1];

%Set length of pilot sequences
tau_p = 10;

%Define an M x M identity matrix
eyeM = eye(M);


%Preallocate matrices for storing the simulation results
NMSE_LMMSE = zeros(length(SNR),length(varphiRange),length(kappatRange));
NMSE_mismatched = zeros(length(SNR),length(varphiRange),length(kapparRange));


%% Go through the range of nominal angles
for r = 1:length(varphiRange)
    
    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);
    
    %Compute the spatial correlation matrix
    R = functionRlocalscattering(M,varphiRange(r),ASDdeg,antennaSpacing);
    
    %Go through all SNRs
    for n = 1:length(SNR)
        
        %Go through all hardware qualities
        for kind = 1:length(kappatRange)
            
            %Extract diagonal of spatial correlation matrix
            Rdiag = diag(diag(R));
            
            %Extract hardware qualities
            kappat = kappatRange(kind);
            kappar = kapparRange(kind);
            
            %Compute the NMSE by using (6.29) to compute tr(C)
            NMSE_LMMSE(n,r,kind) = real(trace(R - SNR(n)*kappat*kappar*tau_p*R*((SNR(n)*(1+kappat*(tau_p-1))*kappar*R + SNR(n)*(1-kappar)*Rdiag + eyeM)\R)))/trace(R);
            
            %Compute deterministic matrix of the MMSE estimator in case of
            %ideal hardware. It will be used as a mismatched estimator
            %under hardware impairments
            A_mismatch = sqrt(SNR(n))*R/(SNR(n)*tau_p*R + eyeM);
            
            %Compute the NMSE with the mismatched estimator
            NMSE_mismatched(n,r,kind) = real( trace(R) - sqrt(kappat*kappar*SNR(n))*tau_p*(trace(A_mismatch*R)+trace(R*A_mismatch')) + tau_p*trace(A_mismatch*(SNR(n)*(1+kappat*(tau_p-1))*kappar*R + SNR(n)*(1-kappar)*Rdiag + eyeM)*A_mismatch') )/M;
            
        end
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(SNRdB,mean(NMSE_mismatched(:,:,1),2),'r--','LineWidth',1);
plot(SNRdB,mean(NMSE_LMMSE(:,:,1),2),'k-','LineWidth',1);
plot(SNRdB,mean(NMSE_LMMSE(:,:,3),2),'b-.','LineWidth',1);
plot(SNRdB,mean(NMSE_mismatched(:,:,2),2),'r--','LineWidth',1);
plot(SNRdB,mean(NMSE_LMMSE(:,:,2),2),'k-','LineWidth',1);

xlabel('Effective SNR [dB]');
ylabel('NMSE');
set(gca,'YScale','log');

legend('Impairments: Mismatched estimator','Impairments: LMMSE estimator','Ideal hardware: MMSE estimator','Location','SouthWest');
