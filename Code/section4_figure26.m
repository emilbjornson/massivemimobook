%This Matlab script can be used to reproduce Figure 4.26 in the monograph:
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
L = 2;

%Number of UEs per BS
K = 2;

%Number of BS antennas
M = 200;

%Select the number of setups of large-scale fading realizations
nbrOfSetups = 100;

%Select the number of channel realizations per setup
nbrOfRealizations = 100;

%Define the range of standard deviations of the large-scale fading variations
largeFadingStds = 0:0.5:4;


%% Scenario setup

%Define BS positions as in Figure 4.24 using complex coordinates
BSpositions = [0 200+1i*200];

%Define UE positions as in Figure 4.24 using complex coordinates
UEdistance = 140; %Distance from serving BS
angleUE1 = 1.2*pi/5;
angleUE2 = pi/2-angleUE1;
UEpositions = [UEdistance*exp(1i*angleUE1) BSpositions(2)-UEdistance*exp(1i*angleUE1); UEdistance*exp(1i*angleUE2) BSpositions(2)-UEdistance*exp(1i*angleUE2)];


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

%Pathloss exponent
alpha = 3.76;

%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

%Define the pilot reuse factor
f = 1;

%Prepare to save simulation results
sumSE_MR = zeros(length(largeFadingStds),nbrOfSetups);
sumSE_ZF = zeros(length(largeFadingStds),nbrOfSetups);
sumSE_SMMSE = zeros(length(largeFadingStds),nbrOfSetups);
sumSE_RZF = zeros(length(largeFadingStds),nbrOfSetups);
sumSE_MMMSE = zeros(length(largeFadingStds),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Prepare to store normalized spatial correlation matrices
    R = zeros(M,M,K,L,L,length(largeFadingStds));
    
    %Prepare to store average channel gain numbers (in dB)
    channelGaindB = zeros(K,L,L);
    
    %Go through all the cells
    for l = 1:L
        
        %Go through all BSs
        for j = 1:L
            
            %Compute distance between UEs in cell l and BS j
            distancesBSj = abs(UEpositions(:,l)-BSpositions(j));
            
            %Compute average channel gain using the large-scale fading model in
            %(2.3), while neglecting the shadow fading
            channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
            
            
            %Go through all UEs in cell l
            for k = 1:K
                
                %Generate large-scale fading realizations over the array
                variationsOverArray = randn(M,1);
                
                for c = 1:length(largeFadingStds)
                    
                    %Compute spatial correlation matrix for a particular
                    %standard deviation of the large-scale fading variation
                    R(:,:,k,l,j,c) = diag(10.^(largeFadingStds(c)*variationsOverArray/10));
                    
                end
                
            end
            
        end
        
    end
    
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all standard deviation of large-scale fading variations
    for c = 1:length(largeFadingStds)
        
        %Generate channel realizations with estimates and estimation
        [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R(:,:,:,:,:,c),channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
        
        %Compute SE using Monte-Carlo realizations
        [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
        
        %Save results
        sumSE_MR(c,n) = mean(sum(SE_MR,1));
        sumSE_ZF(c,n) = mean(sum(SE_ZF,1));
        sumSE_SMMSE(c,n) = mean(sum(SE_SMMSE,1));
        sumSE_RZF(c,n) = mean(sum(SE_RZF,1));
        sumSE_MMMSE(c,n) = mean(sum(SE_MMMSE,1));
        
        %Delete large matrices
        clear Hhat C Rscaled;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot simulation results
figure;
hold on; box on;

plot(largeFadingStds,mean(sumSE_MMMSE,2),'rd-','LineWidth',1);
plot(largeFadingStds,mean(sumSE_SMMSE,2),'b:','LineWidth',1);
plot(largeFadingStds,mean(sumSE_RZF,2),'k-.','LineWidth',1);
plot(largeFadingStds,mean(sumSE_ZF,2),'r--','LineWidth',1);
plot(largeFadingStds,mean(sumSE_MR,2),'bs-','LineWidth',1);

xlabel('Standard deviation of large-scale fading variations over array');
ylabel('Average sum SE [bit/s/Hz/cell]');

legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
