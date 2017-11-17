%This Matlab script can be used to reproduce Figure 1.14 in the monograph:
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

%Define range of number of BS antennas
M = 1:100;

%Extract the maximal number of BS antennas
Mmax = max(M);

%Select number of Monte Carlo realizations of the Rayleigh fading
numberOfRealizations = 100000;

%Generate random UE angles from 0 to 2*pi
varphiDesired = 2*pi*rand(1,numberOfRealizations);
varphiInterfering = 2*pi*rand(1,numberOfRealizations);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Preallocate matrices for storing the simulation results
SE_LoS = zeros(Mmax,numberOfRealizations);
SE_NLoS_montecarlo = zeros(Mmax,numberOfRealizations);


%Generate uncorrelated Rayleigh fading realizations
fadingRealizationsDesired = (randn(Mmax,numberOfRealizations)+1i*randn(Mmax,numberOfRealizations))/sqrt(2);
fadingRealizationsInterfering = (randn(Mmax,numberOfRealizations)+1i*randn(Mmax,numberOfRealizations))/sqrt(2);


%Go through different channel realizations
for k = 1:numberOfRealizations
    
    %Compute the argument, that appears in (1.28), for different UE angles
    argument = (2*pi*antennaSpacing*( sin(varphiDesired(k))  - sin(varphiInterfering(k))) );
    
    %Compute the SE under line-of-sight (LoS) propagation, for a particular
    %set of UE angles, using (1.27). Note that the g-function is
    %implemented slightly differently from (1.28)
    SE_LoS(:,k) = log2(1 + SNR*M  ./ (  (SNR*betabar/(1-cos(argument)))*(1-cos(argument*M)) ./ M + 1 ) );
    
    %Compute the SE under NLoS propagation for one fading realization
    SE_NLoS_montecarlo(:,k) = log2(1 + SNR*cumsum(abs(fadingRealizationsDesired(:,k)).^2) ./ (  betabar*SNR*(abs(cumsum(conj(fadingRealizationsDesired(:,k)).*fadingRealizationsInterfering(:,k))).^2)./cumsum(abs(fadingRealizationsDesired(:,k)).^2) +1 ) );
    
end

%Compute the lower bound on SE under NLoS propagation in (1.32)
SE_NLoS_lower = log2(1 + SNR*(M-1)  ./ ( SNR*betabar + 1 ) );



%Compute the SE under NLoS propagation using (1.29)
SE_NLoS =  (exp(1/(SNR*betabar))*expint(1/(SNR*betabar))/log(2))*(1./(1-1/betabar).^M -1 );

%First summation in (1.29) by considering all the M values at once
for m = 1:Mmax
    
    disp([num2str(m) ' antenna numbers out of ' num2str(length(Mmax))]);
    
    %First term that depend only on m
    term1 = 1/betabar/(1-1/betabar)^m;
    
    %Second summation
    for l = 0:Mmax-m
        
        %Second term that depend only on m and l
        term2 = term1 * (-1).^(M-m-l+1) ./ (SNR.^(M-m-l)  .* factorial(abs(M-m-l)) * log(2));
        
        %Third and fourth summation
        term3 = exp(1/SNR)*expint(1/SNR);
        
        for n = 1:l
            
            for j = 0:n-1
                
                term3 = term3 + 1/n/factorial(j)/SNR^j;
                
            end
            
        end
        
        %Determine which of the considered M values that contain the current term
        feasible = (m<=M) & (l<=(M-m));
        
        %Add the terms to the SE computation
        SE_NLoS(feasible) = SE_NLoS(feasible) + term2(feasible)*term3;
        
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(M,mean(SE_LoS,2),'k-','LineWidth',1);
plot(M(1),mean(SE_NLoS_montecarlo(1,:),2),'bd-.','LineWidth',1);
plot(M,SE_NLoS_lower,'r--','LineWidth',1);
plot(M(10:10:end),mean(SE_NLoS_montecarlo(10:10:end,:),2),'bd','LineWidth',1);
plot(M,SE_NLoS,'b-.','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Average SE [bit/s/Hz]');

legend('LoS','NLoS: Rayleigh fading','NLoS: Rayleigh fading (lower bound)','Location','SouthEast');
