%This Matlab script can be used to reproduce Figures 4.21-22 in the monograph:
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
f = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 10;

%Select the number of channel realizations per setup
nbrOfRealizations = 1000;


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

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Prepare to save simulation results
sumchannelGains_MR = zeros(K,L,K,L,nbrOfSetups);
sumchannelGains_RZF = zeros(K,L,K,L,nbrOfSetups);
sumchannelGains_MMMSE = zeros(K,L,K,L,nbrOfSetups);



%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,C,tau_p,~,H] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
    
    %Compute channel gain terms between all users
    [channelGains_MR,channelGains_RZF,channelGains_MMMSE] = functionComputeULDLPowerLevels(H,Hhat,C,nbrOfRealizations,M,K,L,p);
    
    %Delete large matrices
    clear Hhat C H R;
    
    %Save results
    sumchannelGains_MR(:,:,:,:,n) = channelGains_MR;
    sumchannelGains_RZF(:,:,:,:,n) = channelGains_RZF;
    sumchannelGains_MMMSE(:,:,:,:,n) = channelGains_MMMSE;
    
end


%Prepare to compute signal and interference gains in UL and DL
signalGainsUL_MR = zeros(K,L,nbrOfSetups);
interferenceGainsUL_MR = zeros(K,L,nbrOfSetups);

signalGainsDL_MR = zeros(K,L,nbrOfSetups);
interferenceGainsDL_MR = zeros(K,L,nbrOfSetups);

signalGainsUL_MMMSE = zeros(K,L,nbrOfSetups);
interferenceGainsUL_MMMSE = zeros(K,L,nbrOfSetups);

signalGainsDL_MMMSE = zeros(K,L,nbrOfSetups);
interferenceGainsDL_MMMSE = zeros(K,L,nbrOfSetups);


%Go through all cells
for j = 1:L
    
    %Go through all UEs
    for k = 1:K
        
        %Extract the interference gains with MR in the UL
        signalGainsUL_MR(k,j,:) = reshape(p*sumchannelGains_MR(k,j,k,j,:),[1 1 nbrOfSetups]);
        interferenceGainsUL_MR(k,j,:) = reshape(p*sum(sum(sumchannelGains_MR(k,j,:,:,:),3),4),[1 1 nbrOfSetups])-signalGainsUL_MR(k,j,:);
        
        %Extract the interference gains with MR in the DL
        signalGainsDL_MR(k,j,:) = reshape(rho*sumchannelGains_MR(k,j,k,j,:),[1 1 nbrOfSetups]);
        interferenceGainsDL_MR(k,j,:) = reshape(rho*sum(sum(sumchannelGains_MR(:,:,k,j,:),1),2),[1 1 nbrOfSetups])-signalGainsDL_MR(k,j,:);
        
        %Extract the interference gains with M-MMSE in the UL
        signalGainsUL_MMMSE(k,j,:) = reshape(p*sumchannelGains_MMMSE(k,j,k,j,:),[1 1 nbrOfSetups]);
        interferenceGainsUL_MMMSE(k,j,:) = reshape(p*sum(sum(sumchannelGains_MMMSE(k,j,:,:,:),3),4),[1 1 nbrOfSetups])-signalGainsUL_MMMSE(k,j,:);
        
        %Extract the interference gains with M-MMSE in the DL
        signalGainsDL_MMMSE(k,j,:) = reshape(rho*sumchannelGains_MMMSE(k,j,k,j,:),[1 1 nbrOfSetups]);
        interferenceGainsDL_MMMSE(k,j,:) = reshape(rho*sum(sum(sumchannelGains_MMMSE(:,:,k,j,:),1),2),[1 1 nbrOfSetups])-signalGainsDL_MMMSE(k,j,:);
        
        
    end
    
end


%Plot simulation results in the UL with MR
figure; hold on; box on;

plot(10*log10(signalGainsUL_MR(:)),10*log10(interferenceGainsUL_MR(:)),'k*','LineWidth',1);
xlabel('Signal power over noise floor [dB]');
ylabel('Interference power over noise floor [dB]');
axis([-10 70 0 60]);

%Plot simulation results in the UL with M-MMSE
figure; hold on; box on;

plot(10*log10(signalGainsUL_MMMSE(:)),10*log10(interferenceGainsUL_MMMSE(:)),'ks','LineWidth',1);
xlabel('Signal power over noise floor [dB]');
ylabel('Interference power over noise floor [dB]');
axis([-10 70 0 60]);


%Plot simulation results in the DL with MR
figure; hold on; box on;

plot(10*log10(signalGainsDL_MR(:)),10*log10(interferenceGainsDL_MR(:)),'k*','LineWidth',1);
xlabel('Signal power over noise floor [dB]');
ylabel('Interference power over noise floor [dB]');
axis([-10 70 0 60]);


%Plot simulation results in the DL with M-MMSE
figure; hold on; box on;

plot(10*log10(signalGainsDL_MMMSE(:)),10*log10(interferenceGainsDL_MMMSE(:)),'ks','LineWidth',1);
xlabel('Signal power over noise floor [dB]');
ylabel('Interference power over noise floor [dB]');
axis([-10 70 0 60]);
