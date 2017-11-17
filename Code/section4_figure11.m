%This Matlab script can be used to reproduce Figure 4.11 in the monograph:
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

%Select length of coherence block
tau_c = 200;

%Use the exact Gaussian local scattering model
accuracy = 1;

%Range of angular standard deviation in the local scattering model (in degrees)
ASDdegRange = [1e-4 5:5:50];


%Prepare to save simulation results
sumSE_MR_theorem = zeros(length(ASDdegRange),nbrOfSetups);
sumSE_MR_UatF_corollary = zeros(length(ASDdegRange),nbrOfSetups);
sumSE_MR_UatF_alt1 = zeros(length(ASDdegRange),nbrOfSetups);
sumSE_MR_UatF_alt2 = zeros(length(ASDdegRange),nbrOfSetups);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup and many ASDs
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdegRange);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all ASDs
    for spr = 1:length(ASDdegRange)
        
        %Output simulation progress
        disp([num2str(spr) ' ASDs out of ' num2str(length(ASDdegRange))]);
        
        %Generate channel realizations with estimates and estimation error
        %correlation matrices
        [Hhat,C,tau_p,~,H] = functionChannelEstimates(R(:,:,:,:,:,spr),channelGainOverNoise,nbrOfRealizations,M,K,L,p,f);
        
        %Compute SE using Monte Carlo simulations
        [signal,interference,prelogFactor,SE_MR] = functionComputeSE_UL_MR_MC(Hhat,H,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p);
        
        
        %Save results
        sumSE_MR_theorem(spr,n) = mean(sum(SE_MR,1));
        
        sumSE_MR_UatF_corollary(spr,n) = prelogFactor*mean(sum(log2(1+signal(:,:,1)./interference(:,:,1)),1));
        
        sumSE_MR_UatF_alt1(spr,n) = prelogFactor*mean(sum(log2(1+signal(:,:,2)./interference(:,:,2)),1));
        
        sumSE_MR_UatF_alt2(spr,n) = prelogFactor*mean(sum(log2(1+signal(:,:,3)./interference(:,:,3)),1));
        
        %Delete large matrices
        clear Hhat H C;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%Compute the average sum SE over setups
sumSE_MR_theorem_opt = mean(sumSE_MR_theorem,2);
sumSE_MR_UatF_corollary_opt = mean(sumSE_MR_UatF_corollary,2);
sumSE_MR_UatF_alt1_opt = mean(sumSE_MR_UatF_alt1,2);
sumSE_MR_UatF_alt2_opt = mean(sumSE_MR_UatF_alt2,2);


%% Plot the simulation results
figure(1);
hold on; box on;

plot(ASDdegRange,sumSE_MR_theorem_opt,'rs--','LineWidth',1);
plot(ASDdegRange,sumSE_MR_UatF_alt2_opt,'b*-','LineWidth',1);
plot(ASDdegRange,sumSE_MR_UatF_alt1_opt,'kd:','LineWidth',1);
plot(ASDdegRange,sumSE_MR_UatF_corollary_opt,'kp--','LineWidth',1);

xlabel('ASD $$(\sigma_{\varphi})$$ [degree]','Interpreter','latex');
ylabel('Average sum SE [bit/s/Hz/cell]');
legend('Original bound (Th. 4.1)','UatF: h/||h||^2','UatF: h/||h||','UatF: h (Corr. 4.5)','Location','NorthEast');
ylim([0 50]);
