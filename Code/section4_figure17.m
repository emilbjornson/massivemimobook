%This Matlab script can be used to reproduce Figure 4.17 in the monograph:
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

%Select range of length of coherence blocks
tau_c_range = round(logspace(1.3,3,50));

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Generate uncorrelated correlation matrices
Runcorr = repmat(eye(M),[1 1 K L L]);


%Prepare to save simulation results
sumSE_hardening_MR = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);
sumSE_hardening_RZF = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);
sumSE_hardening_MMMSE = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);

sumSE_estimation_MR = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);
sumSE_estimation_RZF = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);
sumSE_estimation_MMMSE = zeros(length(tau_c_range),length(fRange),nbrOfSetups,2);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Go through correlated (c=1) and uncorrelated (c=2) channels
    for c = 1:2
        
        %Output simulation progress
        disp([num2str(c) ' correlation models out of ' num2str(2)]);
        
        %Go through all pilot reuse factors
        for s = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(s);
            
            if c == 1 %Correlated fading
                
                %Generate channel realizations with estimates and
                %estimation error correlation matrices
                [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            elseif c == 2 %Uncorrelated fading
                
                %Generate channel realizations with estimates and 
                %estimation error correlation matrices
                [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(Runcorr,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            end
            
            
            %Compute SEs with the estimation bound in Theorem 4.6 using
            %Monte Carlo simulations 
            [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_hardening(H,Hhat,C,Rscaled,tau_c_range,tau_p,nbrOfRealizations,M,K,L,p,rho);

            %Save simulation results
            sumSE_hardening_MR(:,s,n,c) = squeeze(mean(sum(SE_MR,1),2));
            sumSE_hardening_RZF(:,s,n,c) = squeeze(mean(sum(SE_RZF,1),2));
            sumSE_hardening_MMMSE(:,s,n,c) = squeeze(mean(sum(SE_MMMSE,1),2));
            
            
            %Compute SEs with the estimation bound in Theorem 4.9 using
            %Monte Carlo simulations 
            [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat,C,Rscaled,tau_c_range,tau_p,nbrOfRealizations,M,K,L,p,rho);
            
            %Save simulation results
            sumSE_estimation_MR(:,s,n,c) = squeeze(mean(sum(SE_MR,1),2));
            sumSE_estimation_RZF(:,s,n,c) = squeeze(mean(sum(SE_RZF,1),2));
            sumSE_estimation_MMMSE(:,s,n,c) = squeeze(mean(sum(SE_MMMSE,1),2));
            
            %Delete large matrices
            clear H Hhat C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results (c=1 correlated fading, c=2 uncorrelated fading)
for c = 1:2
    
    figure;
    hold on; box on;
    
    plot(tau_c_range,max(mean(sumSE_estimation_RZF(:,:,:,c),3),[],2),'k-.','LineWidth',1);
    plot(tau_c_range,max(mean(sumSE_hardening_RZF(:,:,:,c),3),[],2),'k-','LineWidth',1);
    
    plot(tau_c_range,max(mean(sumSE_estimation_MMMSE(:,:,:,c),3),[],2),'r-.','LineWidth',1);
    plot(tau_c_range,max(mean(sumSE_hardening_MMMSE(:,:,:,c),3),[],2),'r-','LineWidth',1);
    
    plot(tau_c_range,max(mean(sumSE_estimation_MR(:,:,:,c),3),[],2),'b-.','LineWidth',1);
    plot(tau_c_range,max(mean(sumSE_hardening_MR(:,:,:,c),3),[],2),'b-','LineWidth',1);
    
    xlabel('Samples per coherence block (\tau_c)');
    ylabel('Average sum SE [bit/s/Hz/cell]');
    legend('Estimation bound','Hardening bound','Location','NorthWest');
    ylim([0 70]);
    set(gca,'XScale','log');
    
end
