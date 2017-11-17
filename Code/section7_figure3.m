%This Matlab script can be used to reproduce Figure 7.3 in the monograph:
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
nbrOfSetups = 100;

%Select the number of channel realizations per setup
nbrOfRealizations = 200;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Define noise figure at BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Set the range of Delta values in the power control policy of (7.10)
DeltadB = 0:10:20;

%Generate uncorrelated covariance matrices
Runcorr = repmat(eye(M),[1 1 K L L]);


%Prepare to save simulation results
sumSE_MR = zeros(L*K*nbrOfSetups,length(DeltadB));
sumSE_ZF = zeros(L*K*nbrOfSetups,length(DeltadB));
sumSE_SMMSE = zeros(L*K*nbrOfSetups,length(DeltadB));
sumSE_RZF = zeros(L*K*nbrOfSetups,length(DeltadB));
sumSE_MMMSE = zeros(L*K*nbrOfSetups,length(DeltadB));


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoiseOriginal = channelGaindB - noiseVariancedBm;
    
    
    %Go through all values of Delta 
    for s = 1:length(DeltadB)
        
        %Extract the average channel gains before power control
        channelGainOverNoise = channelGainOverNoiseOriginal;
        
        %Go through all cells
        for j = 1:L
            
            %Compute beta_j,min in the power control policy of (7.10)
            betajMin = min(channelGainOverNoiseOriginal(:,j,j));
            
            %Scale the average channel gains by applying the power control
            %policy of (7.10). Note that we are including this power
            %control here so that we can then view it as if all UEs of 
            %transmitting at maximum power
            differenceSNR = channelGainOverNoiseOriginal(:,j,j)-betajMin;
            backoff = differenceSNR-DeltadB(s);
            backoff(backoff<0) = 0;
            
            channelGainOverNoise(:,j,:) = channelGainOverNoiseOriginal(:,j,:)-repmat(backoff,[1 1 L]);
            
        end

        
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
        
        %Compute SEs using Theorem 4.1
        [SE_MR,SE_ZF,SE_SMMSE,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
        
        %Save results
        sumSE_MR(1+(n-1)*K*L:n*K*L,s) = SE_MR(:);
        sumSE_ZF(1+(n-1)*K*L:n*K*L,s) = SE_ZF(:);
        sumSE_SMMSE(1+(n-1)*K*L:n*K*L,s) = SE_SMMSE(:);
        sumSE_RZF(1+(n-1)*K*L:n*K*L,s) = SE_RZF(:);
        sumSE_MMMSE(1+(n-1)*K*L:n*K*L,s) = SE_MMMSE(:);
        
        %Delete large matrices
        clear Hhat C Rscaled;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results

figure;
hold on; box on;

plot(sort(sumSE_MR(:,1)),linspace(0,1,K*L*nbrOfSetups),'k-','LineWidth',1);
plot(sort(sumSE_MR(:,2)),linspace(0,1,K*L*nbrOfSetups),'b-.','LineWidth',1);
plot(sort(sumSE_MR(:,3)),linspace(0,1,K*L*nbrOfSetups),'r--','LineWidth',1);
xlabel('SE per UE [bit/s/Hz]');
ylabel('CDF');
legend('0 dB','10 dB','20 dB','Location','SouthEast');
xlim([0 8]);

figure;
hold on; box on;

plot(sort(sumSE_RZF(:,1)),linspace(0,1,K*L*nbrOfSetups),'k-','LineWidth',1);
plot(sort(sumSE_RZF(:,2)),linspace(0,1,K*L*nbrOfSetups),'b-.','LineWidth',1);
plot(sort(sumSE_RZF(:,3)),linspace(0,1,K*L*nbrOfSetups),'r--','LineWidth',1);
xlabel('SE per UE [bit/s/Hz]');
ylabel('CDF');
legend('0 dB','10 dB','20 dB','Location','SouthEast');
xlim([0 8]);

figure;
hold on; box on;

plot(sort(sumSE_MMMSE(:,1)),linspace(0,1,K*L*nbrOfSetups),'k-','LineWidth',1);
plot(sort(sumSE_MMMSE(:,2)),linspace(0,1,K*L*nbrOfSetups),'b-.','LineWidth',1);
plot(sort(sumSE_MMMSE(:,3)),linspace(0,1,K*L*nbrOfSetups),'r--','LineWidth',1);
xlabel('SE per UE [bit/s/Hz]');
ylabel('CDF');
legend('0 dB','10 dB','20 dB','Location','SouthEast');
xlim([0 8]);
