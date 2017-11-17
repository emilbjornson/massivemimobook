%This Matlab script can be used to reproduce Figure 4.8 in the monograph:
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


%Number of BSs
L = 16;

%Number of UEs per BS
K = 10;

%Number of BS antennas
M = 100;

%Define the pilot reuse factor
f = 2;

%Select the number of setups with random UE locations
nbrOfSetups = 50;

%Select the number of channel realizations per setup
nbrOfRealizations = 100;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;

%Use the exact Gaussian local scattering model
accuracy = 1;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Prepare to save simulation results
userSE_MR = zeros(length(ASDdeg)+1,K*L,nbrOfSetups);
userSE_RZF = zeros(length(ASDdeg)+1,K*L,nbrOfSetups);
userSE_MMMSE = zeros(length(ASDdeg)+1,K*L,nbrOfSetups);

%Generate uncorrelated correlation matrices
Runcorr = repmat(eye(M),[1 1 K L L]);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Go through all ASDs, including uncorrelated fading
    for spr = 1:length(ASDdeg)+1
        
        %Output simulation progress
        disp([num2str(spr) ' correlation models out of ' num2str(length(ASDdeg)+1)]);
        
        
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        if spr<=length(ASDdeg) %Consider correlated Rayleigh fading
            
            [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R(:,:,:,:,:,spr),channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
            
        elseif spr==length(ASDdeg)+1 %Consider uncorrelated Rayleigh fading
            
            [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(Runcorr,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
            
        end
        
        
        %Compute SEs using Theorem 4.1
        [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
        
        %Save average sum SE per cell
        userSE_MR(spr,:,n) = SE_MR(:);
        userSE_RZF(spr,:,n) = SE_RZF(:);
        userSE_MMMSE(spr,:,n) = SE_MMMSE(:);
        
        %Delete large matrices
        clear Hhat C Rscaled;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results
CDFnumbers = linspace(0,1,K*L*nbrOfSetups);

figure(1);
hold on; box on;

SE_MR = reshape(userSE_MR(1,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_MR(:)),CDFnumbers,'k-','LineWidth',1);

SE_MR = reshape(userSE_MR(2,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_MR(:)),CDFnumbers,'k--','LineWidth',1);

SE_RZF = reshape(userSE_RZF(1,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_RZF(:)),CDFnumbers,'b-','LineWidth',1);

SE_RZF = reshape(userSE_RZF(2,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_RZF(:)),CDFnumbers,'b--','LineWidth',1);

SE_MMMSE = reshape(userSE_MMMSE(1,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_MMMSE(:)),CDFnumbers,'r-','LineWidth',1);

SE_MMMSE = reshape(userSE_MMMSE(2,:,:),[K*L*nbrOfSetups 1]);
plot(sort(SE_MMMSE(:)),CDFnumbers,'r--','LineWidth',1);

xlabel('SE per UE [bit/s/Hz]');
ylabel('CDF');
xlim([0 12]);
legend('Correlated fading','Uncorrelated fading','Location','SouthEast');
