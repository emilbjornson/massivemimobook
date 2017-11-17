%This Matlab script can be used to reproduce Figure 1.16 in the monograph:
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


%Define the SNR
SNR = 1;

%Define betabar (strength of inter-cell interference)
betabar = 1e-1;

%Define the range of number of UEs
K = 1:20;

%Define range of number of BS antennas
M = [10 100];

%Extract the maximal number of UEs and BS antennas
Kmax = max(K);
Mmax = max(M);

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
SE_MR_NLoS_montecarlo = zeros(length(K),length(M));
SE_MR_LoS = zeros(length(K),length(M));
SE_MMSE_NLoS_montecarlo = zeros(length(K),length(M));
SE_MMSE_LoS_montecarlo = zeros(length(K),length(M));


%% Go through all Monte Carlo realizations
for n = 1:numberOfRealizations
    
    %Output simulation progress
    disp([num2str(n) ' realizations out of ' num2str(numberOfRealizations)]);
    
    %Go through the range of number of UEs
    for kindex = 1:length(K);
        
        %Go through the range of number of BS antennas
        for mindex = 1:length(M);
            
            
            %Compute the SE with MR under LoS propagation using (1.43) for
            %one realization of the UE angles
            
            %Compute the uplink SE with MR combining
            argumentsDesired = 2*pi*antennaSpacing*( repmat(sin(varphiDesired(1,1:K(kindex),n)),[K(kindex) 1])  -  repmat(sin(varphiDesired(1,1:K(kindex),n))',[1 K(kindex)]) );
            argumentsInterfering = 2*pi*antennaSpacing*( repmat(sin(varphiDesired(1,1:K(kindex),n)),[K(kindex) 1])  -  repmat(sin(varphiInterfering(1,1:K(kindex),n))',[1 K(kindex)]) );
            oneminuscos = (1-cos(argumentsDesired)) + eye(K(kindex));
            
            %Compute the uplink SE with MR combining
            SE_MR_LoS(kindex,mindex) = SE_MR_LoS(kindex,mindex) + sum(log2(1 + SNR*M(mindex)*ones(1,K(kindex))  ./ ( (SNR/M(mindex))*sum( (1-cos(M(mindex)*argumentsDesired)) ./ oneminuscos,1) + (SNR/M(mindex))*betabar*sum( (1-cos(M(mindex)*argumentsInterfering)) ./ (1-cos(argumentsInterfering)), 1) + 1)))/numberOfRealizations;
            
            
            %Compute the SE with M-MMSE under LoS propagation using (1.43)
            
            %Compute the M-MMSE combining vectors
            MMMSEvectors = ( SNR* H_LoS_desired(1:M(mindex),1:K(kindex),n)*H_LoS_desired(1:M(mindex),1:K(kindex),n)' + SNR*H_LoS_interfering(1:M(mindex),1:K(kindex),n)*H_LoS_interfering(1:M(mindex),1:K(kindex),n)' + eye(M(mindex)) ) \ (SNR*H_LoS_desired(1:M(mindex),1:K(kindex),n));
            
            %Compute the intra-cell channel powers after M-MMSE combining
            channelgainsIntracell = abs(MMMSEvectors'*H_LoS_desired(1:M(mindex),1:K(kindex),n)).^2;
            
            %Extract desired signal power for each UE
            signalpowers = diag(channelgainsIntracell);
            
            %Extract and compute interference powers for each UE
            interferencepowers = sum(channelgainsIntracell,2) - signalpowers + sum(abs(MMMSEvectors'*H_LoS_interfering(1:M(mindex),1:K(kindex),n)).^2,2);
            
            %Compute the effective 1/SNR after noise amplification
            scalednoisepower = (1/SNR)*sum(abs(MMMSEvectors').^2,2);
            
            %Compute the uplink SE with M-MMSE combining
            SE_MMSE_LoS_montecarlo(kindex,mindex) = SE_MMSE_LoS_montecarlo(kindex,mindex) + sum(log2(1 + signalpowers./(interferencepowers+scalednoisepower)))/numberOfRealizations;
            
            
            
            %Compute the SE with MR under NLoS propagation using the first line
            %in (1.44) for one realization of the Rayleigh fading
            
            %Compute the MR combining vectors
            MRvectors = H_NLoS_desired(1:M(mindex),1:K(kindex),n);
            
            %Compute the intra-cell effective channel gains after MR filtering
            channelgainsIntracell = abs(MRvectors'*H_NLoS_desired(1:M(mindex),1:K(kindex),n)).^2;
            
            %Extract desired signal power for each UE
            signalpowers = diag(channelgainsIntracell);
            
            %Extract and compute interference powers for each UE
            interferencepowers = sum(channelgainsIntracell,2) - signalpowers + sum(abs(MRvectors'*H_NLoS_interfering(1:M(mindex),1:K(kindex),n)).^2,2);
            
            %Compute the effective 1/SNR after noise amplification
            scalednoisepower = (1/SNR)*sum(abs(MRvectors').^2,2);
            
            %Compute the uplink SE with MR combining
            SE_MR_NLoS_montecarlo(kindex,mindex) = SE_MR_NLoS_montecarlo(kindex,mindex) + sum(log2(1 + signalpowers./(interferencepowers+scalednoisepower)))/numberOfRealizations;
            
            
            %Compute the SE with M-MMSE under NLoS propagation for one
            %realization of the Rayleigh fading
            
            %Compute the M-MMSE combining filter
            MMMSEvectors = ( SNR* H_NLoS_desired(1:M(mindex),1:K(kindex),n)*H_NLoS_desired(1:M(mindex),1:K(kindex),n)' + SNR* H_NLoS_interfering(1:M(mindex),1:K(kindex),n)*H_NLoS_interfering(1:M(mindex),1:K(kindex),n)' + eye(M(mindex)) ) \ (SNR*H_NLoS_desired(1:M(mindex),1:K(kindex),n));
            
            %Compute the intra-cell effective channel gains after MSE filtering
            channelgainsIntracell = abs(MMMSEvectors'*H_NLoS_desired(1:M(mindex),1:K(kindex),n)).^2;
            
            %Extract desired signal power for each UE
            signalpowers = diag(channelgainsIntracell);
            
            %Extract and compute interference powers for each UE
            interferencepowers = sum(channelgainsIntracell,2) - signalpowers + sum(abs(MMMSEvectors'*H_NLoS_interfering(1:M(mindex),1:K(kindex),n)).^2,2);
            
            %Compute the effective 1/SNR after noise amplification
            scalednoisepower = (1/SNR)*sum(abs(MMMSEvectors').^2,2);
            
            %Compute the uplink SE with M-MMSE combining
            SE_MMSE_NLoS_montecarlo(kindex,mindex) = SE_MMSE_NLoS_montecarlo(kindex,mindex) + sum(log2(1 + signalpowers./(interferencepowers+scalednoisepower)))/numberOfRealizations;
            
            
        end
        
    end
    
end


%Compute the lower bound in (1.44) on the SE under NLoS propagation
SE_MR_NLoS_lower = zeros(length(K),length(M));

for mindex = 1:length(M);
    SE_MR_NLoS_lower(:,mindex) = K .* log2(1 + SNR*(M(mindex)-1)  ./ ( (K-1)*SNR + K*SNR*betabar + 1 ) );
end


%% Plot the simulation results for MR combining
figure(1);
hold on; box on;

for mindex = 1:length(M)
    
    plot(K,SE_MR_LoS(:,mindex),'k-','LineWidth',1);
    plot(K(1),SE_MR_NLoS_montecarlo(1,mindex),'bd-.','LineWidth',1);
    plot(K,SE_MR_NLoS_lower(:,mindex)','r--','LineWidth',1);
    plot(K(5:5:end),SE_MR_NLoS_montecarlo(5:5:end,mindex),'bd','LineWidth',1);
    plot(K,SE_MR_NLoS_montecarlo,'b-.','LineWidth',1);
end

xlabel('Number of UEs (K)');
ylabel('Average sum SE [bit/s/Hz/cell]');

legend('LoS','NLoS','NLoS (lower bound)','Location','NorthWest');


%% Plot the simulation results for M-MMSE combining
figure(2);
hold on; box on;

for mindex = 1:length(M)
    
    plot(K,SE_MMSE_LoS_montecarlo(:,mindex),'k-','LineWidth',1);
    plot(K,SE_MMSE_NLoS_montecarlo(:,mindex),'b-.','LineWidth',1);
    
end

xlabel('Number of UEs (K)');
ylabel('Average sum SE [bit/s/Hz/cell]');

legend('LoS','NLoS','Location','NorthWest');
ylim([0 120]);
