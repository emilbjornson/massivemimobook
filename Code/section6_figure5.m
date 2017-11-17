%This Matlab script can be used to reproduce Figure 6.5 in the monograph:
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

%Define pilot reuse factor
f = 2;

%Select the number of setups with random UE locations
nbrOfSetups = 10;

%Select the number of channel realizations per setup
nbrOfRealizations = 1000;

%Select the range of hardware qualities. The values at the same position in
%the different vectors are considered simultaneously.
kappatUE = [0.95 0.99 1];
kapparBS = kappatUE;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Define total uplink transmit power per UE (mW)
p = 100;

%Define noise figure at BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
powerterms_MR_total = zeros(length(kappatUE),6);
powerterms_RZF_total = zeros(length(kappatUE),6);
powerterms_MMMSE_total = zeros(length(kappatUE),6);


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all hardware quality values
    for r = 1:length(kappatUE)
        
        %Generate channel realizations with estimates and estimation
        %error correlation matrices
        [Hhat,C,~,~,H] = functionChannelEstimates_impairments(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f,kappatUE(r),kapparBS(r));
        
        %Compute received power divided into six categories
        [powerterms_MR,powerterms_RZF,powerterms_MMMSE] = functionComputeULPowerLevels_impairments(H,Hhat,C,nbrOfRealizations,M,K,L,p,f,kappatUE(r),kapparBS(r));

        %Average over different setups
        powerterms_MR_total(r,:) = powerterms_MR_total(r,:) + powerterms_MR/nbrOfSetups;
        powerterms_RZF_total(r,:) = powerterms_RZF_total(r,:) + powerterms_RZF/nbrOfSetups;
        powerterms_MMMSE_total(r,:) = powerterms_MMMSE_total(r,:) + powerterms_MMMSE/nbrOfSetups;
        
        %Delete large matrices
        clear Hhat C H;
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot simulation results
figure(1);
bar(10*log10(powerterms_MR_total));
legend('Desired signal','Interf: Same pilot','Interf: Other pilots','Transmit dist.','Receiver dist.','Self-dist./interf.');
ylabel('Power over noise floor [dB]');
set(gca,'xTickLabel',{'\kappa=0.95','\kappa=0.99',' \kappa=1 '});
colormap('hot');
ylim([-10 60]);

figure(2);
bar(10*log10(powerterms_RZF_total));
legend('Desired signal','Interf: Same pilot','Interf: Other pilots','Transmit dist.','Receiver dist.','Self-dist./interf.');
ylabel('Power over noise floor [dB]');
set(gca,'xTickLabel',{'\kappa=0.95','\kappa=0.99',' \kappa=1 '});
colormap('hot');
ylim([-10 60]);

figure(3);
bar(10*log10(powerterms_MMMSE_total));
legend('Desired signal','Interf: Same pilot','Interf: Other pilots','Transmit dist.','Receiver dist.','Self-dist./interf.');
ylabel('Power over noise floor [dB]');
set(gca,'xTickLabel',{'\kappa=0.95','\kappa=0.99',' \kappa=1 '});
colormap('hot');
ylim([-10 60]);
