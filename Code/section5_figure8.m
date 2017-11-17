%This Matlab script can be used to reproduce Figure 5.8 in the monograph:
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

%Range of antenna-UE ratios
cRange = [1 2 4 8];

%Number of UEs
K = 10;

%Fixed circuit power per BS (in Watt)
P_FIX = 10;

%Circuit power per BS antenna (in Watt)
P_BS = 1;

%Circuit power UE (in Watt)
P_UE = 0.5;

%Range of SE values
sumSE = (0:0.00002:40)';

%Select ratio between noise power and beta_0^0
sigma2_beta = 10^(-6*0.1);

%Compute nu_0 in (5.14)
nu_0 = sigma2_beta/mu;

%Strength of inter-cell interference  = -15 dB
barbeta = 10^(-15/10);


%Prepare to save simulation results
EE = zeros(length(sumSE),length(cRange));
maxEE = zeros(length(cRange),1);
maxSE = zeros(length(cRange),1);


%% Go through range of antenna-UE ratios
for index1 = 1:length(cRange)
    
    %Compute number of BS antennas
    M = cRange(index1)*K;
    
    
    %Go through range of SE values
    for index2 = 1:length(sumSE)
        
        p = sigma2_beta*((M-1)/(2^sumSE(index2) - 1) - K*barbeta +1 - K)^-1;
        
        if p < 0
            
            break;
            
        end
        
        EE(index2,index1) = B*K*sumSE(index2)/(K*p/mu + P_FIX + M*P_BS + K*P_UE);
        
    end
    
    [max_value, index_max ] = max(EE(:,index1));
    
    maxEE(index1) = max_value;
    
    maxSE(index1) = K*sumSE(index_max);
    
end


%% Plot the simulation results
figure(1);
hold on; box on;

plot(K*sumSE,EE(:,1),'k','LineWidth',1);
plot(K*sumSE,EE(:,2),'k-.','LineWidth',1);
plot(K*sumSE,EE(:,3),'k--','LineWidth',1);
plot(K*sumSE,EE(:,4),'k:','LineWidth',1);

set(gca,'YScale','log');
xlabel('Sum SE [bit/s/Hz/cell]');
ylabel('EE [bit/Joule]');
legend('M/K = 1','M/K = 2','M/K = 4','M/K = 8','Location','NorthEast');

axis([0 max(sumSE) 10^2 10^6])

semilogy(maxSE,maxEE,'ro','LineWidth',1,'MarkerFaceColor','r');
