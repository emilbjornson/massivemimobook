%This Matlab script can be used to reproduce Figure 5.7 in the monograph:
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

%Number of BS antennas
M = 10;

%Define range of number of UEs
Krange = [5 10 30];

%Fixed circuit power per BS (in Watt)
P_FIX = 10;

%Circuit power per BS antenna (in Watt)
P_BS = 1;

%Circuit power UE (in Watt)
P_UE = 0.5;

%Range of sum SE values
sumSE = (0:0.00002:14)';

%Select ratio between noise power and beta_0^0
sigma2_beta = 10^(-6*0.1);

%Compute nu_0 in (5.14)
nu_0 = sigma2_beta/mu;

%Define range of strengths of inter-cell interference = [-3 dB, -15 dB]
barbetaRange = [10^(-3/10) 10^(-15/10)];


%% Go through range of strengths of inter-cell interference
for index3 = 1:length(barbetaRange)
    
    %Prepare to save simulation results
    txPower = zeros(length(sumSE),length(Krange));
    EE = zeros(length(sumSE),length(Krange));
    maxEE = zeros(length(Krange),1);
    maxSumSE = zeros(length(Krange),1);
    
    
    %Go through range of number of UEs
    for index1 = 1:length(Krange)
        
        %Extract strength of inter-cell interference
        barbeta = barbetaRange(index3);
        
        %Extract number of UEs
        K = Krange(index1);
        
        
        %Go through range of SE values
        for index2 = 1:length(sumSE)
            
            %Compute transmit power using (5.27)
            p = sigma2_beta*((M-1)/(2^sumSE(index2) - 1) - K*barbeta +1 - K)^(-1);
            
            %Skip SEs for which there are no feasible transmit powers
            if p < 0
                
                break;
                
            end
            
            %Compute the EE using (5.28)
            EE(index2,index1) = B*K*sumSE(index2)/(K*p/mu + P_FIX + M*P_BS + K*P_UE);
            
            %Compute transmit power
            txPower(index2,index1) = K*p/mu;
            
        end
        
        %Find EE-maximizing point on each curve
        [max_value, index_max ] = max(EE(:,index1));
        maxEE(index1) = max_value;
        maxSumSE(index1) = K*sumSE(index_max);
        
    end
    
    %Plot the simulation reults
    figure(1); hold on; box on;
    plot(Krange(1)*sumSE,EE(:,1),'k','LineWidth',1);
    plot(Krange(2)*sumSE,EE(:,2),'k-.','LineWidth',1);
    plot(Krange(3)*sumSE,EE(:,3),'k:','LineWidth',1);
    plot(maxSumSE,maxEE,'ro','LineWidth',1,'MarkerFaceColor','r');
    
end

set(gca,'YScale','log');
axis([0 14 10^2 10^6])

xlabel('Sum SE [bit/s/Hz/cell]');
ylabel('EE [bit/Joule]');
legend('K = 5','K = 10', 'K = 30','Location','NorthEast');
