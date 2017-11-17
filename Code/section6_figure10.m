%This Matlab script can be used to reproduce Figure 6.10 in the monograph:
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
Mrange = logspace(0,4,50);

%Define the pilot reuse factor
f = 2;

%Select the number of setups with random UE locations
nbrOfSetups = 100;

%Hardware quality of the UEs' transmitters
kappatUE = 0.997;

%Hardware quality of the BSs' receivers before applying the
%hardware-quality scaling law in Corollary 6.8
kapparBS_original = 0.997;

%Select range of scaling exponents in the hardware-quality scaling law in
%Corollary 6.8
eValues = [0 1/243 1 1/81 1/27 1/9 1/3];


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


%Prepare to save simulation results
sumSE_MR = zeros(length(Mrange),length(eValues),nbrOfSetups);
sumSE_MR_asymptotic = zeros(nbrOfSetups,1);


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
    
    
    %Go through all scaling exponents
    for r = 1:length(eValues)
        
        %Output simulation progress
        disp([num2str(r) ' scaling exponent out of ' num2str(length(eValues))]);
        
        
        %Go through all number of BS antennas
        for m = 1:length(Mrange)
            
            %Compute the BSs' receivers by following the hardware-quality
            %scaling law in Corollary 6.8
            kapparBS = kapparBS_original/Mrange(m)^eValues(r);
            
            %Compute UL SE using Corollary 6.4
            [SE_MR,SE_MR_asymptotic] = functionComputeSE_UL_MR_impairments(channelGainOverNoise,tau_c,Mrange(m),K,L,p,f,kappatUE,kapparBS);
            
            %Save results
            sumSE_MR(m,r,n) = sum(mean(SE_MR,1),2);
            
        end
        
    end
    
    %Compute the average sum SE in the asymptotic limit
    sumSE_MR_asymptotic(n) = sum(mean(SE_MR_asymptotic,1),2);
end


%% Plot the simulation results
figure;
hold on; box on;

plot([min(Mrange) max(Mrange)],mean(sumSE_MR_asymptotic)*ones(1,2),'k:','LineWidth',1);

plot(Mrange,max(mean(sumSE_MR(:,1,:),3),[],2),'r-','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR(:,2,:),3),[],2),'b-.','LineWidth',1);
plot(Mrange,max(mean(sumSE_MR(:,3,:),3),[],2),'k--','LineWidth',1);

for r = 4:length(eValues)
    
    plot(Mrange,max(mean(sumSE_MR(:,r,:),3),[],2),'b-.','LineWidth',1);
    
end

xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('Asymptotic limit','Fixed','Scaling law','Faster than scaling law','Location','NorthWest');
set(gca,'XScale','log');


figure;
hold on; box on;

plot(Mrange,sqrt(1-kapparBS_original./(Mrange.^eValues(1))),'r-','LineWidth',1);
plot(Mrange,sqrt(1-kapparBS_original./(Mrange.^eValues(2))),'b-.','LineWidth',1);
plot(Mrange,sqrt(1-kapparBS_original./(Mrange.^eValues(3))),'k--','LineWidth',1);

for r = 4:length(eValues)
    
    plot(Mrange,sqrt(1-kapparBS_original./(Mrange.^eValues(r))),'b-.','LineWidth',1);
    
end

xlabel('Number of antennas (M)');
ylabel('EVM');
legend('Fixed','Scaling law','Faster than scaling law','Location','NorthWest');
set(gca,'XScale','log');
