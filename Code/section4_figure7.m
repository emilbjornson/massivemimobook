%This Matlab script can be used to reproduce Figure 4.7 in the monograph:
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
nbrOfSetups = 10;

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

%Select length of coherence block
tau_c = 200;

%Use the exact Gaussian local scattering model
accuracy = 1;

%Range of angular standard deviation in the local scattering model (in degrees)
ASDdegRange = [1e-4 5:5:50];

%Generate uncorrelated correlation matrices
Runcorr = repmat(eye(M),[1 1 K L L]);

%Prepare to save simulation results
sumSE_MR = zeros(length(ASDdegRange)+1,length(fRange),nbrOfSetups);
sumSE_RZF = zeros(length(ASDdegRange)+1,length(fRange),nbrOfSetups);
sumSE_MMMSE = zeros(length(ASDdegRange)+1,length(fRange),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdegRange);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    %Go through all ASDs, including uncorrelated fading
    for spr = 1:length(ASDdegRange)+1
        
        %Output simulation progress
        disp([num2str(spr) ' correlation models out of ' num2str(length(ASDdegRange)+1)]);
        
        %Go through all pilot reuse factors
        for s = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(s);
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            if spr<=length(ASDdegRange) %Consider correlated Rayleigh fading
                
                [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R(:,:,:,:,:,spr),channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            elseif spr==length(ASDdegRange)+1 %Consider uncorrelated Rayleigh fading
                
                [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(Runcorr,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
                
            end
            
            
            %Compute SEs using Theorem 4.1
            [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
            
            %Save average sum SE per cell
            sumSE_MR(spr,s,n) = mean(sum(SE_MR,1));
            sumSE_RZF(spr,s,n) = mean(sum(SE_RZF,1));
            sumSE_MMMSE(spr,s,n) = mean(sum(SE_MMMSE,1));
            
            %Delete large matrices
            clear Hhat C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%Compute the average sum SE over setups and maximize with respect to the
%pilot reuse factors
sumSE_MMSE_opt = max(mean(sumSE_MMMSE,3),[],2);
sumSE_RZF_opt = max(mean(sumSE_RZF,3),[],2);
sumSE_MR_opt = max(mean(sumSE_MR,3),[],2);


%% Plot the simulation results
figure(1);
hold on; box on;

plot(ASDdegRange,sumSE_MMSE_opt(1:end-1),'rd-','LineWidth',1);
plot(ASDdegRange,sumSE_RZF_opt(1:end-1),'k-.','LineWidth',1);
plot(ASDdegRange,sumSE_MR_opt(1:end-1),'bs-','LineWidth',1);

plot([min(ASDdegRange) max(ASDdegRange)],sumSE_MMSE_opt(end)*ones(1,2),'k:','LineWidth',1);
plot([min(ASDdegRange) max(ASDdegRange)],sumSE_RZF_opt(end)*ones(1,2),'k:','LineWidth',1);
plot([min(ASDdegRange) max(ASDdegRange)],sumSE_MR_opt(end)*ones(1,2),'k:','LineWidth',1);

xlabel('ASD $$(\sigma_{\varphi})$$ [degree]','Interpreter','latex');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('M-MMSE','RZF','MR','Location','SouthWest');
ylim([0 70]);
