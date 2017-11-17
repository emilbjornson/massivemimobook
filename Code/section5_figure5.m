%This Matlab script can be used to reproduce Figure 5.5 in the monograph:
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


%% Propagation and hardware parameters

%Communication bandwidth
B = 0.1*10^(6); %100 kHz

%PA efficiency
mu = 0.4;

%Range of number of BS antennas
Mrange = [1000 100 10 2];

%Fixed circuit power per BS (in Watt)
P_FIX = 10;

%Range of SE values
SE = (0:0.0001:20)';

%Select ratio between noise power and beta_0^0
sigma2_beta = 10^(-6*0.1);

%Compute nu_0 in (5.14)
nu_0 = sigma2_beta/mu;


%Prepare to save simulation results
EE = zeros(length(SE),length(Mrange));
maxEE = zeros(length(Mrange),1);
maxSE = zeros(length(Mrange),1);
maxEE_theory = zeros(length(Mrange),1);
maxSE_theory = zeros(length(Mrange),1);


%% Go through range of number of antennas
for index1 = 1:length(Mrange)
    
    %Extract number of antennas
    M = Mrange(index1);
    
    
    %Go through range of SE values
    for index2 = 1:length(SE)
        
        %Compute transmit power using (5.12)
        p = (2^SE(index2) - 1)/(M-1)*sigma2_beta;
        
        %Compute the EE using (5.11)
        EE(index2,index1) = B*SE(index2)/(p/mu + P_FIX);
        
    end
    
    %Find the EE-maximizing point on each curve
    [max_value, index_max ] = max(EE(:,index1));
    maxEE(index1) = max_value;
    maxSE(index1) = SE(index_max);
    
    %Find the EE-maximizing pair of SE and EE values, using (5.18) and (5.19)
    argument  = (M-1)*P_FIX/(nu_0*exp(1)) - 1/exp(1);
    maxSE_theory(index1) = (lambertw(argument) + 1)/log(2);
    maxEE_theory(index1) = (M-1)*B*2^(-maxSE_theory(index1))/(nu_0*log(2));
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(SE,EE(:,1),'k','LineWidth',1);
plot(SE,EE(:,2),'k-.','LineWidth',1);
plot(SE,EE(:,3),'k--','LineWidth',1);
plot(SE,EE(:,4),'k:','LineWidth',1);

set(gca,'YScale','log');

legend('M = 1000','M = 100','M = 10','M = 2','Location','NorthEast');
xlabel('SE [bit/s/Hz]');
ylabel('EE [bit/Joule]');

axis([0 max(SE) 10^2 10^7])

plot(maxSE,maxEE,'ko','LineWidth',1);
plot(maxSE_theory,maxEE_theory,'ro-','LineWidth',1,'MarkerFaceColor','r');
