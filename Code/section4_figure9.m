%This Matlab script can be used to reproduce Figure 4.9 in the monograph:
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
nbrOfSetups = 10;

%Select the number of channel realizations per setup
nbrOfRealizations = 1000;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Total uplink transmit power per UE (mW)
p = 100;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
channelVariations_MR = zeros(K,L,nbrOfSetups,2);
channelVariations_RZF = zeros(K,L,nbrOfSetups,2);
channelVariations_MMMSE = zeros(K,L,nbrOfSetups,2);

%Generate uncorrelated correlation matrices
Runcorr = repmat(eye(M),[1 1 K L L]);


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
        
        if c == 1 %Correlated fading
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
            
        elseif c == 2 %Uncorrelated fading
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(Runcorr,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
            
        end
        
        %Compute channel gain terms between all users
        [variation_MR,variation_RZF,variation_MMMSE] = functionComputeGainVariations(H,Hhat,C,nbrOfRealizations,M,K,L,p);
        
        %Save results
        channelVariations_MR(:,:,n,c) = variation_MR;
        channelVariations_RZF(:,:,n,c) = variation_RZF;
        channelVariations_MMMSE(:,:,n,c) = variation_MMMSE;
        
        %Delete large matrices
        clear Hhat C Rscaled;
        
    end
    
    %Delete large matrices
    clear R;
    
end



%% Go through correlated and uncorrelated scenarios
for c = 1:2
    
    variations_MR = real(channelVariations_MR(:,:,:,c));
    variations_RZF = real(channelVariations_RZF(:,:,:,c));
    variations_MMMSE = real(channelVariations_MMMSE(:,:,:,c));
    
    y = linspace(0,1,length(variations_MR(:)));
    
    
    %Plot simulation results
    figure(1); hold on; box on;
    plot(sort(variations_MMMSE(:)),y,'r-','LineWidth',1);
    plot(sort(variations_RZF(:)),y,'k-.','LineWidth',1);
    plot(sort(variations_MR(:)),y,'b--','LineWidth',1);
    
    xlabel('Variance in (4.12)');
    ylabel('CDF');
    legend('M-MMSE','RZF','MR','Location','SouthEast');
    set(gca,'XScale','log');
    
end

%Plot reference curves
Mplot = [1 30 100];

for n = 1:length(Mplot)
    
    plot((Mplot(n)*(gamma(Mplot(n))/gamma(Mplot(n)+0.5))^2-1)*[1; 1],[0; 1],'k:','LineWidth',1);
    
end
