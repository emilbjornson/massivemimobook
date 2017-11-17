%This Matlab script can be used to reproduce Figure 4.20 in the monograph:
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

%Define the range of pilot reuse factors
fRange = [1 2 4];

%Select the number of setups with random UE locations
nbrOfSetups = 50;

%Select the number of channel realizations per setup
nbrOfRealizations = 500;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per UE (mW)
rho = 100;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
sumSE_MR = zeros(length(fRange),3,nbrOfSetups);
sumSE_RZF = zeros(length(fRange),3,nbrOfSetups);
sumSE_MMMSE = zeros(length(fRange),3,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Go through all pilot reuse factors
    for s = 1:length(fRange)
        
        %Output simulation progress
        disp([num2str(s) ' reuse factors out of ' num2str(length(fRange))]);
        
        %Extract pilot reuse factor
        f = fRange(s);
        
        %Generate channel realizations with estimates and estimation error
        %correlation matrices, using MMSE, EW-MMSE, and LS estimators
        [Hhat_MMSE,C_MMSE,tau_p,Rscaled,H,Hhat_EW_MMSE,C_EW_MMSE,Hhat_LS,C_LS] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
        
        
        %Compute SEs with the estimation bound in Theorem 4.6 using
        %Monte Carlo simulations with the MMSE estimator
        [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE] = functionComputeSE_DL_hardening(H,Hhat_MMSE,C_MMSE,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Compute SEs with the estimation bound in Theorem 4.9 using
        %Monte Carlo simulations with the MMSE estimator
        [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat_MMSE,C_MMSE,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Use the largest of the two bounds
        SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
        SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
        SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);
        
        %Save results
        sumSE_MR(s,1,n) = mean(sum(SE_MR,1));
        sumSE_RZF(s,1,n) = mean(sum(SE_RZF,1));
        sumSE_MMMSE(s,1,n) = mean(sum(SE_MMMSE,1));
        
        
        %Compute SEs with the estimation bound in Theorem 4.6 using
        %Monte Carlo simulations with the EW-MMSE estimator
        [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE] = functionComputeSE_DL_hardening(H,Hhat_EW_MMSE,C_EW_MMSE,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Compute SEs with the estimation bound in Theorem 4.9 using
        %Monte Carlo simulations with the EW-MMSE estimator
        [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat_EW_MMSE,C_EW_MMSE,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Use the largest of the two bounds
        SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
        SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
        SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);
        
        %Save results
        sumSE_MR(s,2,n) = mean(sum(SE_MR,1));
        sumSE_RZF(s,2,n) = mean(sum(SE_RZF,1));
        sumSE_MMMSE(s,2,n) = mean(sum(SE_MMMSE,1));
        
        
        %Compute SEs with the estimation bound in Theorem 4.6 using
        %Monte Carlo simulations with the LS estimator
        [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE] = functionComputeSE_DL_hardening(H,Hhat_LS,C_LS,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Compute SEs with the estimation bound in Theorem 4.9 using
        %Monte Carlo simulations with the LS estimator
        [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat_LS,C_LS,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho);
        
        %Use the largest of the two bounds
        SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
        SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
        SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);
        
        %Save results
        sumSE_MR(s,3,n) = mean(sum(SE_MR,1));
        sumSE_RZF(s,3,n) = mean(sum(SE_RZF,1));
        sumSE_MMMSE(s,3,n) = mean(sum(SE_MMMSE,1));
        
        %Delete large matrices
        clear Hhat_MMSE C_MMSE Rscaled H Hhat_EW_MMSE C_EW_MMSE Hhat_LS C_LS;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results
sumSEs = [max(mean(sumSE_MMMSE,3),[],1)' max(mean(sumSE_RZF,3),[],1)' max(mean(sumSE_MR,3),[],1)'];

bar(sumSEs);
ylabel('Average sum SE [bit/s/Hz/cell]');
ylim([0 60]);
set(gca,'xticklabel',['  MMSE '; 'EW-MMSE'; '   LS  '])
legend('M-MMSE','RZF','MR','Location','NorthEast');
colormap(hot);
