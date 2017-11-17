%This Matlab script can be used to reproduce Figure 7.5 in the monograph:
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
fRange = [1 2];

%Only one snapshot with random UE locations
nbrOfSetups = 1;

%Select the number of random pilot assignments
nbrOfRandomAssignments = 1000;

%Select the number of channel realizations per setup
nbrOfChannelRealizations = 500;


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

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
allSE_MR = zeros(L*K*nbrOfSetups,nbrOfRandomAssignments,length(fRange));
allSE_RZF = zeros(L*K*nbrOfSetups,nbrOfRandomAssignments,length(fRange));
allSE_MMMSE = zeros(L*K*nbrOfSetups,nbrOfRandomAssignments,length(fRange));


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [Roriginal,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoiseOriginal = channelGaindB - noiseVariancedBm;
    
    
    %Go through all random pilot assignments
    for s = 1:nbrOfRandomAssignments
        
        %Output simulation progress
        disp([num2str(s) ' pilot assignment out of ' num2str(nbrOfRandomAssignments)]);
        
        %Make a copy of the UEs' channel statistics
        R = Roriginal;
        channelGainOverNoise = channelGainOverNoiseOriginal;
        
        %Go through all cells except the first one
        for j = 2:L
            
            %Randomize the UE indicies in the cell j
            order = randperm(K);
            
            %Move around the UEs in cell j
            for l = 1:L
                
                R(:,:,:,j,l) = Roriginal(:,:,order,j,l);
                channelGainOverNoise(:,j,l) = channelGainOverNoiseOriginal(order,j,l);
                
            end
            
        end
        
        
        %Go through all pilot reuse factors
        for m = 1:length(fRange)
            
            %Extract pilot reuse factor
            f = fRange(m);
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R,channelGainOverNoise,nbrOfChannelRealizations,M,K,L,p,f);
            
            %Compute SEs using Theorem 4.1
            [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfChannelRealizations,M,K,L,p);
            
            %Save results
            allSE_MR(1+(n-1)*K*L:n*K*L,s,m) = SE_MR(:);
            allSE_RZF(1+(n-1)*K*L:n*K*L,s,m) = SE_RZF(:);
            allSE_MMMSE(1+(n-1)*K*L:n*K*L,s,m) = SE_MMMSE(:);
            
            %Delete large matrices
            clear Hhat C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end

%Prepare to plot simulation results for sum SE
sumSE_MR = reshape(mean(allSE_MR,1),[nbrOfRandomAssignments,length(fRange)]);
sumSE_RZF = reshape(mean(allSE_RZF,1),[nbrOfRandomAssignments,length(fRange)]);
sumSE_MMMSE = reshape(mean(allSE_MMMSE,1),[nbrOfRandomAssignments,length(fRange)]);

%Prepare to plot simulation reults for the weakest user's SE in cell 1 (an
%arbitary cell, since all cells are 
SE_MR = reshape(allSE_MR,[K L nbrOfRandomAssignments,length(fRange)]);
SE_RZF = reshape(allSE_RZF,[K L nbrOfRandomAssignments,length(fRange)]);
SE_MMMSE = reshape(allSE_MMMSE,[K L nbrOfRandomAssignments,length(fRange)]);

[~,weakestUE] = min(channelGainOverNoiseOriginal(:,1,1));

minSE_MR = reshape(SE_MR(weakestUE,1,:,:),[nbrOfRandomAssignments,length(fRange)]);
minSE_RZF = reshape(SE_RZF(weakestUE,1,:,:),[nbrOfRandomAssignments,length(fRange)]);
minSE_MMMSE = reshape(SE_MMMSE(weakestUE,1,:,:),[nbrOfRandomAssignments,length(fRange)]);


%% Plot simulation results
figure;
hold on; box on;

plot(K*sort(sumSE_MR(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(K*sort(sumSE_MR(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

plot(K*sort(sumSE_RZF(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(K*sort(sumSE_RZF(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

plot(K*sort(sumSE_MMMSE(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(K*sort(sumSE_MMMSE(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

xlabel('Sum SE [bit/s/Hz/cell]');
ylabel('CDF');
legend('f=1','f=2','Location','SouthWest');
xlim([0 60]);


figure;
hold on; box on;

plot(sort(minSE_MR(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(sort(minSE_MR(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

plot(sort(minSE_RZF(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(sort(minSE_RZF(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

plot(sort(minSE_MMMSE(:,1)),linspace(0,1,nbrOfRandomAssignments),'k-','LineWidth',1);
plot(sort(minSE_MMMSE(:,2)),linspace(0,1,nbrOfRandomAssignments),'b-.','LineWidth',1);

xlabel('SE of weakest UE [bit/s/Hz]');
ylabel('CDF');
legend('f=1','f=2','Location','SouthEast');
