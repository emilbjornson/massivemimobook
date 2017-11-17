%This Matlab script can be used to reproduce Figure 1.9 in the monograph:
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


%Define the SNR range for analytical curves
SNRdB = -10:0.1:30;
SNR = 10.^(SNRdB/10);

%Define the SNR range for Monte Carlo simulations
SNRdB_montecarlo = -10:5:30;
SNR_montecarlo = 10.^(SNRdB_montecarlo/10);

%Define the different beta_bar values (strength of inter-cell interference)
betabar = [1e-1 1e-3]';

%Preallocate matrices for storing the simulation results
SE_LoS = zeros(length(betabar),length(SNR));
SE_NLoS = zeros(length(betabar),length(SNR));
SE_NLoS_montecarlo = zeros(length(betabar),length(SNR_montecarlo));

%Select number of Monte Carlo realizations of the Rayleigh fading
numberOfFadingRealizations = 100000;


%% Go through different strengths of the interference
for b = 1:length(betabar)
    
    %Compute SE under line-of-sight (LoS) propagation as in (1.17)
    SE_LoS(b,:) = log2(1+SNR./(betabar(b)*SNR+1));
    
    
    %Generate uncorrelated Rayleigh fading channel realizations
    fadingRealizationsDesired = (randn(numberOfFadingRealizations,1)+1i*randn(numberOfFadingRealizations,1))/sqrt(2);
    fadingRealizationsInterference = (randn(numberOfFadingRealizations,1)+1i*randn(numberOfFadingRealizations,1))/sqrt(2);
    
    %Compute SE under non-line-of-sight (NLoS) propagation from the first
    %line in (1.18), using Monte Carlo simulations for the channel realizations
    SE_NLoS_montecarlo(b,:) = mean(log2(1+abs(fadingRealizationsDesired).^2*SNR_montecarlo ./ ( abs(fadingRealizationsInterference).^2*SNR_montecarlo*betabar(b) +1)),1);
    
    
    %Compute SE under non-line-of-sight (NLoS) propagation as in (1.18)
    SE_NLoS(b,:) = (exp(1./SNR).*expint(1./SNR) - exp(1./(betabar(b)*SNR)).*expint(1./(betabar(b)*SNR)))/((1-betabar(b))*log(2));
    
end


%% Plot the simulation results
figure;
hold on; box on;

for b = 1:length(betabar)
    
    plot(SNRdB,SE_LoS(b,:),'k-','LineWidth',1);
    plot(SNRdB_montecarlo(1),SE_NLoS_montecarlo(1),'bd-.','LineWidth',1);
    plot(SNRdB_montecarlo,SE_NLoS_montecarlo,'bd','LineWidth',1);
    plot(SNRdB,SE_NLoS(b,:),'b-.','LineWidth',1);
    
end

xlabel('SNR [dB]');
ylabel('Average SE [bit/s/Hz]');

legend('LoS','NLoS','Location','SouthEast');
ylim([0 10]);
