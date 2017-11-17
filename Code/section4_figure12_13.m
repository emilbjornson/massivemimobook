%This Matlab script can be used to reproduce Figures 4.12-13 in the monograph:
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
fRange = [1 2 4 16];

%Select the number of setups with random UE locations
nbrOfSetups = 100;

%Select the number of channel realizations per setup
nbrOfRealizations = 100;


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
strongPowerlevels_MR = zeros(4,length(fRange),nbrOfSetups,2);
strongPowerlevels_RZF = zeros(4,length(fRange),nbrOfSetups,2);
strongPowerlevels_MMMSE = zeros(4,length(fRange),nbrOfSetups,2);
weakPowerlevels_MR = zeros(4,length(fRange),nbrOfSetups,2);
weakPowerlevels_RZF = zeros(4,length(fRange),nbrOfSetups,2);
weakPowerlevels_MMMSE = zeros(4,length(fRange),nbrOfSetups,2);


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
        
        %Go through all pilot reuse factors
        for s = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(s);
            
            if c == 1 %Correlated fading
                
                %Generate channel realizations with estimates and estimation
                %error correlation matrices
                [Hhat,C] = functionChannelEstimates(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            elseif c == 2 %Uncorrelated fading
                
                %Generate channel realizations with estimates and estimation
                %error correlation matrices
                [Hhat,C] = functionChannelEstimates(Runcorr,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            end
            
            %Compute signal and interference powers with MR, RZF, and M-MMSE
            [powerterms_MR,powerterms_RZF,powerterms_MMMSE] = functionComputeULPowerLevels(Hhat,C,nbrOfRealizations,M,K,L,p,f);
            
            %Delete large matrices
            clear Hhat C;
            
            %Go through all cells and save simulation results
            for j = 1:L
                
                %Determine the strongest UE in each cell and average over
                %the signal and interference powers in different cells
                [~,UEind] = max(powerterms_MR(1,:,j)); %Consider MR
                strongPowerlevels_MR(:,s,n,c) = strongPowerlevels_MR(:,s,n,c) + powerterms_MR(:,UEind,j)/L;
                
                [~,UEind] = max(powerterms_RZF(1,:,j)); %Consider RZF
                strongPowerlevels_RZF(:,s,n,c) = strongPowerlevels_RZF(:,s,n,c) + powerterms_RZF(:,UEind,j)/L;
                
                [~,UEind] = max(powerterms_MMMSE(1,:,j)); %Consider M-MMSE
                strongPowerlevels_MMMSE(:,s,n,c) = strongPowerlevels_MMMSE(:,s,n,c) + powerterms_MMMSE(:,UEind,j)/L;
                
                %Determine the weakest UE in each cell and average over
                %the signal and interference powers in different cells
                [~,UEind] = min(powerterms_MR(1,:,j)); %Consider MR
                weakPowerlevels_MR(:,s,n,c) = weakPowerlevels_MR(:,s,n,c) + powerterms_MR(:,UEind,j)/L;
                
                [~,UEind] = min(powerterms_RZF(1,:,j)); %Consider RZF
                weakPowerlevels_RZF(:,s,n,c) = weakPowerlevels_RZF(:,s,n,c) + powerterms_RZF(:,UEind,j)/L;
                
                [~,UEind] = min(powerterms_MMMSE(1,:,j)); %Consider M-MMSE
                weakPowerlevels_MMMSE(:,s,n,c) = weakPowerlevels_MMMSE(:,s,n,c) + powerterms_MMMSE(:,UEind,j)/L;
                
            end
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%Compute the average signal and interference powers with respect to
%different setups with random UE locations
avgStrongPowerlevels_MR = mean(strongPowerlevels_MR,3);
avgStrongPowerlevels_RZF = mean(strongPowerlevels_RZF,3);
avgStrongPowerlevels_MMMSE = mean(strongPowerlevels_MMMSE,3);

avgWeakPowerlevels_MR = mean(weakPowerlevels_MR,3);
avgWeakPowerlevels_RZF = mean(weakPowerlevels_RZF,3);
avgWeakPowerlevels_MMMSE = mean(weakPowerlevels_MMMSE,3);


%Go through correlated and uncorrelated scenarios
for c = 1:2
    
    %Extract signal power
    signalLevelMR = avgStrongPowerlevels_MR(1,:,:,c);
    signalLevelRZF = avgStrongPowerlevels_RZF(1,:,:,c);
    signalLevelMMMSE = avgStrongPowerlevels_MMMSE(1,:,:,c);
    
    %Compute average coherent interference by taking the interference from
    %UEs that cause pilot contamination to a particular UE and subtract the
    %average interference that the other UEs in the same cell are causing
    %(which is an estimate of the non-coherent interference).
    coherentMR = avgStrongPowerlevels_MR(4,:,:,c) - (avgStrongPowerlevels_MR(2,:,:,c) - avgStrongPowerlevels_MR(4,:,:,c))/(K-1);
    coherentMR(4) = 0;
    
    coherentRZF = avgStrongPowerlevels_RZF(4,:,:,c) - (avgStrongPowerlevels_RZF(2,:,:,c) - avgStrongPowerlevels_RZF(4,:,:,c))/(K-1);
    coherentRZF(4) = 0;
    
    coherentMMMSE = avgStrongPowerlevels_MMMSE(4,:,:,c) - (avgStrongPowerlevels_MMMSE(2,:,:,c) - avgStrongPowerlevels_MMMSE(4,:,:,c))/(K-1);
    coherentMMMSE(4) = 0;
    
    %Subtract the estimated coherent interference from the total
    %interference power to get the sum of non-coherent interference
    noncoherentMR = avgStrongPowerlevels_MR(3,:,:,c) - coherentMR;
    noncoherentRZF = avgStrongPowerlevels_RZF(3,:,:,c) - coherentRZF;
    noncoherentMMMSE = avgStrongPowerlevels_MMMSE(3,:,:,c) - coherentMMMSE;
    
    
    %Plot simulation results for Figure 4.12a and Figure 4.13a
    figure; hold on; box on;
    
    plot(1:4,10*log10(signalLevelMR),'kd-','LineWidth',1);
    plot(1:4,10*log10(signalLevelRZF),'k*-','LineWidth',1);
    plot(1:4,10*log10(signalLevelMMMSE),'ks-','LineWidth',1);
    
    plot(1:4,10*log10(noncoherentMR),'bd--','LineWidth',1);
    plot(1:4,10*log10(noncoherentRZF),'b*--','LineWidth',1);
    plot(1:4,10*log10(noncoherentMMMSE),'bs--','LineWidth',1);
    
    plot(1:4,10*log10(coherentMR),'rd-.','LineWidth',1);
    plot(1:4,10*log10(coherentRZF),'r*-.','LineWidth',1);
    plot(1:4,10*log10(coherentMMMSE),'rs-.','LineWidth',1);
    
    legend('MR','RZF','M-MMSE','Location','SouthEast');
    
    set(gca,'xTick',1:4);
    set(gca,'xTickLabel',{'1','2','4','16'});
    xlabel('Pilot reuse factor (f)');
    ylabel('Power over noise floor [dB]');
    ylim([-10 50]);
    title('Strongest UE in the cell');
    
end



%Go through correlated and uncorrelated scenarios
for c = 1:2
    
    %Extract signal power
    signalLevelMR = avgWeakPowerlevels_MR(1,:,:,c);
    signalLevelRZF = avgWeakPowerlevels_RZF(1,:,:,c);
    signalLevelMMMSE = avgWeakPowerlevels_MMMSE(1,:,:,c);
    
    %Compute average coherent interference by taking the interference from
    %UEs that cause pilot contamination to a particular UE and subtract the
    %average interference that the other UEs in the same cell are causing
    %(which is an estimate of the non-coherent interference).
    coherentMR = avgWeakPowerlevels_MR(4,:,:,c) - (avgWeakPowerlevels_MR(2,:,:,c) - avgWeakPowerlevels_MR(4,:,:,c))/(K-1);
    coherentMR(4) = 0;
    
    coherentRZF = avgWeakPowerlevels_RZF(4,:,:,c) - (avgWeakPowerlevels_RZF(2,:,:,c) - avgWeakPowerlevels_RZF(4,:,:,c))/(K-1);
    coherentRZF(4) = 0;
    
    coherentMMMSE = avgWeakPowerlevels_MMMSE(4,:,:,c) - (avgWeakPowerlevels_MMMSE(2,:,:,c) - avgWeakPowerlevels_MMMSE(4,:,:,c))/(K-1);
    coherentMMMSE(4) = 0;
    
    %Subtract the estimated coherent interference from the total
    %interference power to get the sum of non-coherent interference
    noncoherentMR = avgWeakPowerlevels_MR(3,:,:,c) - coherentMR;
    noncoherentRZF = avgWeakPowerlevels_RZF(3,:,:,c) - coherentRZF;
    noncoherentMMMSE = avgWeakPowerlevels_MMMSE(3,:,:,c) - coherentMMMSE;
    
    
    %Plot simulation results
    figure; hold on; box on;
    
    plot(1:4,10*log10(signalLevelMR),'kd-','LineWidth',1);
    plot(1:4,10*log10(signalLevelRZF),'k*-','LineWidth',1);
    plot(1:4,10*log10(signalLevelMMMSE),'ks-','LineWidth',1);
    
    plot(1:4,10*log10(noncoherentMR),'bd--','LineWidth',1);
    plot(1:4,10*log10(noncoherentRZF),'b*--','LineWidth',1);
    plot(1:4,10*log10(noncoherentMMMSE),'bs--','LineWidth',1);
    
    plot(1:4,10*log10(coherentMR),'rd-.','LineWidth',1);
    plot(1:4,10*log10(coherentRZF),'r*-.','LineWidth',1);
    plot(1:4,10*log10(coherentMMMSE),'rs-.','LineWidth',1);
    
    legend('MR','RZF','M-MMSE','Location','SouthEast');
    
    set(gca,'xTick',1:4);
    set(gca,'xTickLabel',{'1','2','4','16'});
    xlabel('Pilot reuse factor (f)');
    ylabel('Power over noise floor [dB]');
    ylim([-10 50]);
    title('Weakest UE in the cell');
    
end
