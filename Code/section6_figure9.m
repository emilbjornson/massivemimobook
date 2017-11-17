%This Matlab script can be used to reproduce Figure 6.9 in the monograph:
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

%Select range of number of BS antennas
Mrange = [10:10:200 300:100:2000];

%Define the pilot reuse factor
f = 1;

%Select the number of setups with random UE locations
nbrOfSetups = 100;

%Select the range of hardware qualities. The values at the same position in
%the different vectors are considered simultaneously.
kappatUE = [0.99 0.99 0.99 0.99]; 
kapparBS = [1 0.99 0.95 0.9];
kapparUE = [0.99 0.99 0.99 0.99];
kappatBS = [1 0.99 0.95 0.9];


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per UE (mW)
rho = 100;

%Noise figure at the BS and UE (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Select length of coherence block
tau_c = 200;


%Prepare to save simulation results
sumSE_MR = zeros(length(Mrange),length(kappatUE),nbrOfSetups);
sumSE_MR_asymptotics = zeros(length(Mrange),length(kappatUE),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute the average channel gains for one setup, while we are not
    %using the spatial correlation matrices that the function can produce
    [~,channelGaindB] = functionExampleSetup(L,K,1,2,0);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all hardware qualities
    for r = 1:length(kappatUE)
        
        %Output simulation progress
        disp([num2str(r) ' hardware qualities out of ' num2str(length(kappatUE))]);
        
        
        %Go through all number of BS antennas
        for m = 1:length(Mrange)
            
            %Compute DL SE using Corollary 6.6
            [SE_MR,SE_MR_asymptotic]  = functionComputeSE_DL_MR_impairments(channelGainOverNoise,tau_c,Mrange(m),K,L,p,rho,f,kappatUE(r),kapparBS(r),kappatBS(r),kapparUE(r));
            
            %Save results
            sumSE_MR(m,r,n) = sum(mean(SE_MR,1),2);
            sumSE_MR_asymptotics(m,r,n) = sum(mean(SE_MR_asymptotic,1),2);
            
        end
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(Mrange,max(mean(sumSE_MR(:,1,:),3),[],2),'b-','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR(:,2,:),3),[],2),'k--','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR(:,3,:),3),[],2),'r-.','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR(:,4,:),3),[],2),'k:','LineWidth',1);

plot(Mrange,max(mean(sumSE_MR_asymptotics(:,1,:),3),[],2),'r--','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR_asymptotics(:,2,:),3),[],2),'r--','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR_asymptotics(:,3,:),3),[],2),'r--','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR_asymptotics(:,4,:),3),[],2),'r--','LineWidth',1);
    
xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('\kappa_t^{BS}=\kappa_r^{BS}=1','\kappa_t^{BS}=\kappa_r^{BS}=0.99','\kappa_t^{BS}=\kappa_r^{BS}=0.95','\kappa_t^{BS}=\kappa_r^{BS}=0.9','Asymptotic limit','Location','SouthEast');
ylim([0 60]);
