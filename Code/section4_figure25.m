%This Matlab script can be used to reproduce Figure 4.25 in the monograph:
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

%Define the range of BS antennas
Mrange = round(logspace(1,3,10));

%Extract maximum number of BS antennas
Mmax = max(Mrange);

%Select the number of channel realizations
nbrOfRealizations = 5000;


%% Scenario setup

%Define BS positions using complex coordinates
BSpositions = [0 200+1i*200];

%Define UE positions using complex coordinates
userDistance = 140;
anglediff = 1.2*pi/5;
anglediff2 = pi/2-anglediff;
UEpositions = [userDistance*exp(1i*anglediff) BSpositions(2)-userDistance*exp(1i*anglediff); userDistance*exp(1i*anglediff2) BSpositions(2)-userDistance*exp(1i*anglediff2)];


%% Propagation model

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

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Prepare to store normalized spatial correlation matrices
R = zeros(Mmax,Mmax,K,L,L);

%Prepare to store average channel gain numbers (in dB)
channelGaindB = zeros(K,L,L);

%Define the pilot reuse factor
f = 1;


%% Go through all the cells
for l = 1:L
    
    %Go through all BSs
    for j = 1:L
        
        %Compute distance between UEs in cell l and BS j
        distancesBSj = abs(UEpositions(:,l)-BSpositions(j));
        
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading
        channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
        
        %Compute nominal angles between UE k in cell l and BS j, and
        %generate spatial correlation matrices for the channels using the
        %local scattering model
        for k = 1:K
            
            angleBSj = angle(UEpositions(k,l)-BSpositions(j));
            
            R(:,:,k,l,j) = functionRlocalscatteringApprox(Mmax,angleBSj,ASDdeg,antennaSpacing);
            
        end
        
    end
    
end


%Prepare to save simulation results
sumSE_MR = zeros(length(Mrange),1);
sumSE_ZF = zeros(length(Mrange),1);
sumSE_SMMSE = zeros(length(Mrange),1);
sumSE_RZF = zeros(length(Mrange),1);
sumSE_MMMSE = zeros(length(Mrange),1);

%Compute the normalized average channel gain, where the normalization
%is based on the noise power
channelGainOverNoise = channelGaindB - noiseVariancedBm;


%% Go through all number of antennas
for m = 1:length(Mrange)
    
    %Output simulation progress
    disp([num2str(m) ' antennas out of ' num2str(length(Mrange))]);
    
    %Generate channel realizations with estimates and estimation
    [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R(1:Mrange(m),1:Mrange(m),:,:,:),channelGainOverNoise,nbrOfRealizations,Mrange(m),K,L,p,f);
    
    %Compute SE using Monte Carlo simulations
    [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p);
    
    %Save results
    sumSE_MR(m) = mean(sum(SE_MR,1));
    sumSE_ZF(m) = mean(sum(SE_ZF,1));
    sumSE_SMMSE(m) = mean(sum(SE_SMMSE,1));
    sumSE_RZF(m) = mean(sum(SE_RZF,1));
    sumSE_MMMSE(m) = mean(sum(SE_MMMSE,1));
    
    %Delete large matrices
    clear Hhat C Rscaled;
    
end

%Delete large matrices
clear R;


%% Plot simulation results
figure;
hold on; box on;

plot(Mrange,sumSE_MMMSE,'rd-','LineWidth',1);
plot(Mrange,sumSE_SMMSE,'b:','LineWidth',1);
plot(Mrange,sumSE_RZF,'k-.','LineWidth',1);
plot(Mrange,sumSE_ZF,'r--','LineWidth',1);
plot(Mrange,sumSE_MR,'bs-','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average sum SE [bit/s/Hz/cell]');

legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
set(gca,'XScale','log');
