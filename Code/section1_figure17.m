%This Matlab script can be used to reproduce Figure 1.17 in the monograph:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), 
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", 
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4, 
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.01 (Last edited: 2019-04-17)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Empty workspace and close figures
close all;
clear;


%Define the SNR
SNR = 1;

%Define betabar (strength of inter-cell interference)
betabar = 1e-1;

%Define the range of number of UEs
K = 1:20;

%Define range of antenna-UE ratios
c = [1 2 4 8];

%Extract the maximal number of UEs and BS antennas
Kmax = max(K);
Mmax = Kmax*max(c);

%Select number of Monte Carlo realizations for the line-of-sight (LoS)
%angles and of the non-line-of-sight (NLoS) Rayleigh fading
numberOfRealizations = 10000;


%Generate NLoS channels using uncorrelated Rayleigh fading
H_NLoS_desired = sqrt(1/2)*(randn(Mmax,Kmax,numberOfRealizations)+1i*randn(Mmax,Kmax,numberOfRealizations));
H_NLoS_interfering = sqrt(betabar/2)*(randn(Mmax,Kmax,numberOfRealizations)+1i*randn(Mmax,Kmax,numberOfRealizations));


%Generate random UE angles from 0 to 2*pi
varphiDesired = 2*pi*rand(1,Kmax,numberOfRealizations);
varphiInterfering = 2*pi*rand(1,Kmax,numberOfRealizations);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2;

%Generate LoS channels with different random UE angles
H_LoS_desired =  exp( repmat((0:Mmax-1)',[1 Kmax numberOfRealizations]) .* repmat(-2i*pi*antennaSpacing*sin(varphiDesired),[Mmax 1 1]) );
H_LoS_interfering =  sqrt(betabar)*exp( repmat((0:Mmax-1)',[1 Kmax numberOfRealizations]) .* repmat(-2i*pi*antennaSpacing*sin(varphiInterfering),[Mmax 1 1]) );


%Preallocate matrices for storing the simulation results
SE_MMSE_NLoS_montecarlo = zeros(length(K),length(c));
SE_MMSE_LoS_montecarlo = zeros(length(K),length(c));


%% Go through all Monte Carlo realizations
for n = 1:numberOfRealizations
    
    %Output simulation progress
    disp([num2str(n) ' realizations out of ' num2str(numberOfRealizations)]);
    
    %Go through the range of number of UEs
    for kindex = 1:length(K)
        
        %Go through the range of antenna-UE ratios
        for cindex = 1:length(c)
            
            %Compute the number of antennas
            M = K(kindex)*c(cindex);
            

            %Compute the SE with M-MMSE under LoS propagation for one
            %realization of the UE angles
            
            %Compute the M-MMSE combining filter
            MMMSEfilter = ( SNR* H_LoS_desired(1:M,1:K(kindex),n)*H_LoS_desired(1:M,1:K(kindex),n)' + SNR*H_LoS_interfering(1:M,1:K(kindex),n)*H_LoS_interfering(1:M,1:K(kindex),n)' + eye(M) ) \ (SNR*H_LoS_desired(1:M,1:K(kindex),n));
            
            %Compute the intra-cell effective channel gains after M-MMSE filtering
            channelgainsIntracell = abs(MMMSEfilter'*H_LoS_desired(1:M,1:K(kindex),n)).^2;
            
            %Extract desired signal power for each UE
            signalpowers = diag(channelgainsIntracell);
            
            %Extract and compute interference powers for each UE
            interferencepowers = sum(channelgainsIntracell,2) - signalpowers + sum(abs(MMMSEfilter'*H_LoS_interfering(1:M,1:K(kindex),n)).^2,2);
            
            %Compute the effective 1/SNR after noise amplification
            scalednoisepower = (1/SNR)*sum(abs(MMMSEfilter').^2,2);
            
            %Compute the uplink SE with M-MMSE combining
            SE_MMSE_LoS_montecarlo(kindex,cindex) = SE_MMSE_LoS_montecarlo(kindex,cindex) + sum(log2(1 + signalpowers./(interferencepowers+scalednoisepower)))/numberOfRealizations;
            
            
            %Compute the SE with M-MMSE under NLoS propagation for one
            %realization of the Rayleigh fading

            %Compute the M-MMSE combining filter
            MMMSEfilter = ( SNR* H_NLoS_desired(1:M,1:K(kindex),n)*H_NLoS_desired(1:M,1:K(kindex),n)' + SNR* H_NLoS_interfering(1:M,1:K(kindex),n)*H_NLoS_interfering(1:M,1:K(kindex),n)' + eye(M) ) \ (SNR*H_NLoS_desired(1:M,1:K(kindex),n));
            
            %Compute the intra-cell effective channel gains after MSE filtering
            channelgainsIntracell = abs(MMMSEfilter'*H_NLoS_desired(1:M,1:K(kindex),n)).^2;
            
            %Extract desired signal power for each UE
            signalpowers = diag(channelgainsIntracell);
            
            %Extract and compute interference powers for each UE
            interferencepowers = sum(channelgainsIntracell,2) - signalpowers + sum(abs(MMMSEfilter'*H_NLoS_interfering(1:M,1:K(kindex),n)).^2,2);
            
            %Compute the effective 1/SNR after noise amplification
            scalednoisepower = (1/SNR)*sum(abs(MMMSEfilter').^2,2);
            
            %Compute the uplink SE with M-MMSE combining
            SE_MMSE_NLoS_montecarlo(kindex,cindex) = SE_MMSE_NLoS_montecarlo(kindex,cindex) + sum(log2(1 + signalpowers./(interferencepowers+scalednoisepower)))/numberOfRealizations;
            
        end
        
    end
    
end


%% Plot the simulation results
figure(1);
hold on; box on;

for cindex = 1:length(c)
    
    plot(K,SE_MMSE_LoS_montecarlo(:,cindex),'k-','LineWidth',1);
    plot(K,SE_MMSE_NLoS_montecarlo(:,cindex),'b-.','LineWidth',1);
    
end

xlabel('Number of UEs (K)');
ylabel('Average sum SE [bit/s/Hz/cell]');

legend('LoS','NLoS','Location','NorthWest');
ylim([0 120]);
