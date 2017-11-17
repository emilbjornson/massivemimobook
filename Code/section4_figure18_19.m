%This Matlab script can be used to reproduce Figures 4.18-19 in the monograph:
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

%Define the range of BS antennas
Mrange = 10:10:100;

%Extract maximum number of BS antennas
Mmax = max(Mrange);

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
sumSE_MR = zeros(length(Mrange),length(fRange),nbrOfSetups);
sumSE_ZF = zeros(length(Mrange),length(fRange),nbrOfSetups);
sumSE_SMMSE = zeros(length(Mrange),length(fRange),nbrOfSetups);
sumSE_RZF = zeros(length(Mrange),length(fRange),nbrOfSetups);
sumSE_MMMSE = zeros(length(Mrange),length(fRange),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,Mmax,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Go through all number of antennas
    for m = 1:length(Mrange)
        
        %Output simulation progress
        disp([num2str(m) ' antennas out of ' num2str(length(Mrange))]);
        
        %Go through all pilot reuse factors
        for s = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(s);
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R(1:Mrange(m),1:Mrange(m),:,:,:),channelGainOverNoise,nbrOfRealizations,Mrange(m),K,L,p,f);
            
            %Compute SEs with the estimation bound in Theorem 4.6 using
            %Monte Carlo simulations 
            [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE,SE_hardening_ZF,SE_hardening_SMMSE] = functionComputeSE_DL_hardening(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p,rho);
            
            %Compute SEs with the estimation bound in Theorem 4.9 using
            %Monte Carlo simulations 
            [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_DL_estimation(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p,rho);
            
            %Use the largest of the two bounds
            SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
            SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
            SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);
            SE_ZF(SE_hardening_ZF>SE_ZF) = SE_hardening_ZF(SE_hardening_ZF>SE_ZF);
            SE_SMMSE(SE_hardening_SMMSE>SE_SMMSE) = SE_hardening_SMMSE(SE_hardening_SMMSE>SE_SMMSE);

            %Save results
            sumSE_MR(m,s,n) = mean(sum(SE_MR,1));
            sumSE_ZF(m,s,n) = mean(sum(SE_ZF,1));
            sumSE_SMMSE(m,s,n) = mean(sum(SE_SMMSE,1));
            sumSE_RZF(m,s,n) = mean(sum(SE_RZF,1));
            sumSE_MMMSE(m,s,n) = mean(sum(SE_MMMSE,1));
            
            %Delete large matrices
            clear H Hhat C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results
for s = 1:length(fRange)
    
    figure(s);
    hold on; box on;
    
    plot(Mrange,mean(sumSE_MMMSE(:,s,:),3),'rd-','LineWidth',1);
    plot(Mrange,mean(sumSE_SMMSE(:,s,:),3),'b:','LineWidth',1);
    plot(Mrange,mean(sumSE_RZF(:,s,:),3),'k-.','LineWidth',1);
    plot(Mrange,mean(sumSE_ZF(:,s,:),3),'r--','LineWidth',1);
    plot(Mrange,mean(sumSE_MR(:,s,:),3),'bs-','LineWidth',1);
    
    xlabel('Number of antennas (M)');
    ylabel('Average sum SE [bit/s/Hz/cell]');
    legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
    ylim([0 60]);
    
end
