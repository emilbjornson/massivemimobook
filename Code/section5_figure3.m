%This Matlab script can be used to reproduce Figure 5.3 in the monograph:
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

%Number of users per BS
K = 10;

%Range of number of BS antennas
M = logspace(0,6,50);

%Define the pilot reuse factor
f = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 20;

%Select the range of scaling exponents in the power-scaling law
eValues = [0 1/2 1];


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW), before applying scaling law
pOriginal = 100;

%Total downlink transmit power per UE (mW), before applying scaling law
rhoOriginal = 100;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;


%Prepare to save simulation results
sumSE_MR = zeros(length(M),length(eValues),nbrOfSetups);
SE_MR_asymptotic = zeros(K,L,nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup. We will not use the
    %spatial correlation matrices in this simulation, since i.i.d. Rayleigh
    %fading is considered.
    [~,channelGaindB] = functionExampleSetup(L,K,1,2,1);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all scaling exponents
    for r = 1:length(eValues)
        
        %Output simulation progress
        disp([num2str(r) ' scaling exponent out of ' num2str(length(eValues))]);
        
        %Go through all number of BS antennas
        for m = 1:length(M)
            
            %Compute transmit powers by following the power-scaling law
            %in Lemma 5.1
            p = pOriginal/M(m)^(eValues(r));
            rho = rhoOriginal/M(m)^(eValues(r));
            
            %Compute SEs with average-normalized MR precoding using Corollary 4.7
            [SE_MR,SE_MR_asymp] = functionComputeSE_DL_MR_impairments(channelGainOverNoise,tau_c,M(m),K,L,p,rho,f);
            
            %Store the average sum SE
            sumSE_MR(m,r,n) = sum(mean(SE_MR,1),2);
            
        end
        
        %Store the asymptotic limit for constant power case
        if r == 1
            SE_MR_asymptotic(:,:,n) = SE_MR_asymp;
        end
        
    end
    
end


%Compute the average sum SE in the asymptotic limit
sumSE_MR_asymptotic = mean(sum(mean(SE_MR_asymptotic,1),2),3);


%% Plot the simulation results
figure;
hold on; box on;

plot([min(M) max(M)],sumSE_MR_asymptotic*ones(1,2),'k:','LineWidth',1);

plot(M,max(mean(sumSE_MR(:,1,:),3),[],2),'r-','LineWidth',1);

plot(M,max(mean(sumSE_MR(:,2,:),3),[],2),'b-.','LineWidth',1);

plot(M,max(mean(sumSE_MR(:,3,:),3),[],2),'k--','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('Asymptotic limit','Fixed power','Scaling: \epsilon=1/2','Scaling: \epsilon=1','Location','NorthWest');
set(gca,'XScale','log');
ylim([0 110]);
