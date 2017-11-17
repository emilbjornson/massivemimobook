%This Matlab script can be used to reproduce Figure 7.7 in the monograph:
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

%Select range of number of UEs per BS
Krange = fliplr([1 3 5 10:5:20 30:10:80]);

%Number of BS antennas
M = 100;

%Extract maximum number of UEs
Kmax = max(Krange);

%Define the range of pilot reuse factors
fRange = [1 2];

%Select number of setups with random user locations
nbrOfSetups = 20;

%Select the number of setups with random UE locations
nbrOfRealizations = 200;


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

%Select length of coherence block for first plot
tau_c = 200;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
sumSE_MR = zeros(length(Krange),length(fRange),nbrOfSetups);
sumSE_RZF = zeros(length(Krange),length(fRange),nbrOfSetups);
sumSE_MMMSE = zeros(length(Krange),length(fRange),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,Kmax,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all number of UEs
    for m = 1:length(Krange)
        
        %Output simulation progress
        disp([num2str(m) ' UE numbers out of ' num2str(length(Krange))]);
        
        %Go through all pilot reuse factors
        for s = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(s);
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R(:,:,1:Krange(m),:,:),channelGainOverNoise,nbrOfRealizations,M,Krange(m),L,p,f);
            
            %Compute SEs with the estimation bound in Theorem 4.6 using
            %Monte Carlo simulations
            [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE] = functionComputeSE_DL_hardening(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,Krange(m),L,p,rho);
            
            %Compute SEs with the estimation bound in Theorem 4.9 using
            %Monte Carlo simulations
            [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,Krange(m),L,p,rho);
            
            %Use the largest of the two bounds
            SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
            SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
            SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);
            
            %Save results
            sumSE_MR(m,s,n) = mean(sum(SE_MR,1));
            sumSE_RZF(m,s,n) = mean(sum(SE_RZF,1));
            sumSE_MMMSE(m,s,n) = mean(sum(SE_MMMSE,1));
            
            %Delete large matrices
            clear H Hhat C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%Select length of coherence block for second plot
tau_c_2 = 400;

sumSE_MMMSE_2 = zeros(size(sumSE_MMMSE));
sumSE_RZF_2 = zeros(size(sumSE_RZF));
sumSE_MR_2 = zeros(size(sumSE_MR));

for s = 1:length(fRange)
    prelogOld = (1-f*Krange/tau_c);
    prelogNew = (1-f*Krange/tau_c_2);
    ratio = (prelogNew./prelogOld)';
    
    for n = 1:nbrOfSetups
        
        sumSE_MMMSE_2(:,s,n) = ratio.*sumSE_MMMSE(:,s,n);
        sumSE_RZF_2(:,s,n) = ratio.*sumSE_RZF(:,s,n);
        sumSE_MR_2(:,s,n) = ratio.*sumSE_MR(:,s,n);
        
    end
    
end



%% Plot the simulation results

figure;
hold on; box on;

plot(Krange,mean(max(sumSE_MMMSE,[],2),3),'rd-','LineWidth',1);
plot(Krange,mean(max(sumSE_RZF,[],2),3),'k-.','LineWidth',1);
plot(Krange,mean(max(sumSE_MR,[],2),3),'bs-','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('M-MMSE','RZF','MR','Location','NorthWest');
ylim([0 110]);


figure;
hold on; box on;

plot(Krange,mean(max(sumSE_MMMSE_2,[],2),3),'rd-','LineWidth',1);
plot(Krange,mean(max(sumSE_RZF_2,[],2),3),'k-.','LineWidth',1);
plot(Krange,mean(max(sumSE_MR_2,[],2),3),'bs-','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('M-MMSE','RZF','MR','Location','NorthWest');
ylim([0 110]);
