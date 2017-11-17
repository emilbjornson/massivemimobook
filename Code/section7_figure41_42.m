%This Matlab script can be used to reproduce Figures 7.41-42 in the monograph:
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
%
%Note: This script require additional software packages to be used, which
%need to be downloaded and installed separately. These packages are
%developed independently and are delivered with separate licenses.
%
%The channels are generated using QuaDRiGa from the Fraunhofer Heinrich
%Hertz Institute (http://www.quadriga-channel-model.de). This script has
%been tested with QuaDRiGa version 1.4.8-571.
%
%The downlink power allocation is optimized using CVX from CVX Research,
%Inc. (http://cvxr.com/cvx/). This script has been tested with CVX 2.1,
%using the solver Mosek, version 7.1.0.12. We discourage the use of the
%solvers SDPT3 and SeDuMi since these crashed during the test.


%Empty workspace and close figures
close all;
clear;


%Number of BSs
L = 16;

%Number of UEs to drop within the area closest to a BS
Kdrop = 10;

%Maximum number of UEs to be served per BS
Kmax = 15;

%Pilot reuse factor
f = 2;

%Select the number of setups with random UE locations
nbrOfSetups = 50;

%Select the number of channel realizations per setup
nbrOfRealizations = 400;

%Fractions of UL and DL data
ULfraction = 1/3;
DLfraction = 2/3;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Effective bandwidth after removing 5% overhead from cyclic prefix
B_effective = B*0.95;

%Noise figure at the BS and UE (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Define total uplink transmit power per UE (mW)
p = 100;

%Maximum downlink transmit power per BS (mW)
Pmax = 1000;

%Select length of coherence block
tau_c = 190;

%Compute pilot length
tau_p = Kmax*f;


%Prepare to save simulation results
SE_MR_DL_maxmin = zeros(tau_p,L,nbrOfSetups,3);
SE_RZF_DL_maxmin = zeros(tau_p,L,nbrOfSetups,3);

SE_MR_DL_prodSINR = zeros(tau_p,L,nbrOfSetups,3);
SE_RZF_DL_prodSINR = zeros(tau_p,L,nbrOfSetups,3);

SE_MR_UL_20control = zeros(tau_p,L,nbrOfSetups,3);
SE_RZF_UL_20control = zeros(tau_p,L,nbrOfSetups,3);

SE_MR_UL_0control = zeros(tau_p,L,nbrOfSetups,3);
SE_RZF_UL_0control = zeros(tau_p,L,nbrOfSetups,3);

nbrOfUEs = zeros(nbrOfSetups,3);


%% Go through all antenna configurations
for config = 1:3
    
    %Output simulation progress
    disp(['Antenna configuration ' num2str(config) ' out of 3']);
    
    if config == 1 %10x5x2
        
        %Select range of BS antennas
        M = 100;
        
        %Select number of polarizations
        polarizations = 2;
        
    elseif config == 2 %20x5x1
        
        %Select range of BS antennas
        M = 100;
        
        %Select number of polarizations
        polarizations = 1;
        
    elseif config == 3 %20x5x2
        
        %Select range of BS antennas
        M = 200;
        
        %Select number of polarizations
        polarizations = 2;
        
    end
    
    
    %Go through all setups
    for n = 1:nbrOfSetups
        
        %Output simulation progress
        disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        
        %Generate channel realizations for current setup
        [H,Rest,activeUEs] = functionExampleSetup_Quadriga(L,Kdrop,B,noiseVariancedBm,Kmax,f,M,polarizations);
        
        %Update how many UEs that have been scheduled
        nbrOfUEs(n,config) = sum(activeUEs(:));
        
        
        %Compute DL results
        
        %Compute the prelog factor for the DL
        prelogFactorDL = DLfraction*(1-tau_p/tau_c);
        
        %Compute LS channel estimates
        [Hhat_LS,C_LS] = functionChannelEstimates_LS(H,Rest,nbrOfRealizations,M,tau_p,L,p);
        
        %Compute the signal and interference terms of the DL SEs using the
        %hardening bound in Theorem 4.6
        [signal_MR,interf_MR,signal_RZF,interf_RZF] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c,tau_p,nbrOfRealizations,M,tau_p,L,p);
        
        
        %Delete large matrices
        clear Hhat_LS C_LS;
        
        %Optimize transmit powers for DL for all scenarios
        disp('DL max-min optimization: MR');
        SE_MR_DL_maxmin(:,:,n,config) = functionPowerOptimization_maxmin(signal_MR,interf_MR,Pmax,prelogFactorDL);
        
        disp('DL max-min optimization: RZF');
        SE_RZF_DL_maxmin(:,:,n,config) = functionPowerOptimization_maxmin(signal_RZF,interf_RZF,Pmax,prelogFactorDL);
        
        disp('Max product optimization: MR');
        SE_MR_DL_prodSINR(:,:,n,config) = functionPowerOptimization_prodSINR(signal_MR,interf_MR,Pmax,prelogFactorDL);
        
        disp('Max product optimization: RZF');
        SE_RZF_DL_prodSINR(:,:,n,config) = functionPowerOptimization_prodSINR(signal_RZF,interf_RZF,Pmax,prelogFactorDL);
        
        
        
        %Compute UL results
        
        %20 dB power control
        powerDiffdB = 20;
        
        Hscaled = H;
        RestScaled = Rest;
        
        %Apply the power control policy in (7.11)
        for j = 1:L
            
            betaValues = 10*log10(squeeze(sum(sum(Rest(:,:,activeUEs(:,j)==1,j,j),1),2)/M));
            
            betaMin = min(betaValues);
            
            differenceSNR = betaValues-betaMin;
            backoff = differenceSNR-powerDiffdB;
            backoff(backoff<0) = 0;
            
            activeIndices = find(activeUEs(:,j));
            
            for k = 1:length(activeIndices)
                
                Hscaled(:,:,activeIndices(k),j,:) = H(:,:,activeIndices(k),j,:)/10^(backoff(k)/20);
                
                RestScaled(:,:,activeIndices(k),j,:) = Rest(:,:,activeIndices(k),j,:)/10^(backoff(k)/10);
                
            end
            
        end
        
        %Compute LS channel estimates
        [Hhat_LS,C_LS] = functionChannelEstimates_LS(Hscaled,RestScaled,nbrOfRealizations,M,tau_p,L,p);
        
        %Compute UL SEs using the hardening bound
        [SE_MR,SE_RZF] = functionComputeSE_UL_impairments(Hscaled,Hhat_LS,C_LS,tau_c,tau_p,nbrOfRealizations,M,tau_p,L,p,1,1);
        
        %Store simulation results
        SE_MR_UL_20control(:,:,n,config) = ULfraction*SE_MR;
        SE_RZF_UL_20control(:,:,n,config) = ULfraction*SE_RZF;
        
        %Delete large matrices
        clear Hhat_LS C_LS;
        
        
        %0 dB power control
        powerDiffdB = 0;
        
        Hscaled = H;
        RestScaled = Rest;
        
        %Apply the power control policy in (7.11)
        for j = 1:L
            
            betaValues = 10*log10(squeeze(sum(sum(Rest(:,:,activeUEs(:,j)==1,j,j),1),2)/M));
            
            betaMin = min(betaValues);
            
            differenceSNR = betaValues-betaMin;
            backoff = differenceSNR-powerDiffdB;
            backoff(backoff<0) = 0;
            
            activeIndices = find(activeUEs(:,j));
            
            for k = 1:length(activeIndices)
                
                Hscaled(:,:,activeIndices(k),j,:) = H(:,:,activeIndices(k),j,:)/10^(backoff(k)/20);
                
                RestScaled(:,:,activeIndices(k),j,:) = Rest(:,:,activeIndices(k),j,:)/10^(backoff(k)/10);
                
            end
            
        end
        
        %Compute LS channel estimates
        [Hhat_LS,C_LS] = functionChannelEstimates_LS(Hscaled,RestScaled,nbrOfRealizations,M,tau_p,L,p);
        
        %Compute UL SEs using the hardening bound
        [SE_MR,SE_RZF] = functionComputeSE_UL_impairments(Hscaled,Hhat_LS,C_LS,tau_c,tau_p,nbrOfRealizations,M,tau_p,L,p,1,1);
        
        %Store simulation results
        SE_MR_UL_0control(:,:,n,config) = ULfraction*SE_MR;
        SE_RZF_UL_0control(:,:,n,config) = ULfraction*SE_RZF;
        
        %Delete large matrices
        clear Hhat_LS C_LS;
        
        %Delete large matrices
        clear H Rest;
        
    end
    
end


%% Plot the simulation results
for config = [3 2 1]
    
    %Only keep the SE results for the active UEs
    SE_MR_DL_maxmin_active = SE_MR_DL_maxmin(:,:,:,config);
    SE_MR_DL_maxmin_active = sort(SE_MR_DL_maxmin_active(find(SE_MR_DL_maxmin_active>0)));
    
    SE_RZF_DL_maxmin_active = SE_RZF_DL_maxmin(:,:,:,config);
    SE_RZF_DL_maxmin_active = sort(SE_RZF_DL_maxmin_active(find(SE_RZF_DL_maxmin_active>0)));
    
    SE_MR_DL_prodSINR_active = SE_MR_DL_prodSINR(:,:,:,config);
    SE_MR_DL_prodSINR_active = sort(SE_MR_DL_prodSINR_active(find(SE_MR_DL_prodSINR_active>0)));
    
    SE_RZF_DL_prodSINR_active = SE_RZF_DL_prodSINR(:,:,:,config);
    SE_RZF_DL_prodSINR_active = sort(SE_RZF_DL_prodSINR_active(find(SE_RZF_DL_prodSINR_active>0)));
    
    SE_MR_UL_20control_active = SE_MR_UL_20control(:,:,:,config);
    SE_MR_UL_20control_active = sort(SE_MR_UL_20control_active(find(SE_MR_UL_20control_active>0)));
    
    SE_RZF_UL_20control_active = SE_RZF_UL_20control(:,:,:,config);
    SE_RZF_UL_20control_active = sort(SE_RZF_UL_20control_active(find(SE_RZF_UL_20control_active>0)));
    
    SE_MR_UL_0control_active = SE_MR_UL_0control(:,:,:,config);
    SE_MR_UL_0control_active = sort(SE_MR_UL_0control_active(find(SE_MR_UL_0control_active>0)));
    
    SE_RZF_UL_0control_active = SE_RZF_UL_0control(:,:,:,config);
    SE_RZF_UL_0control_active = sort(SE_RZF_UL_0control_active(find(SE_RZF_UL_0control_active>0)));
    
    
    %Extract total number of UEs
    Ktotal = sum(nbrOfUEs(:,config));
    
    %Set line styles
    if config == 1
        
        MR = 'k-.';
        RZF = 'b-.';
        
    elseif config == 2
        
        MR = 'k--';
        RZF = 'b--';
        
    elseif config == 3
        
        MR = 'k-';
        RZF = 'b-';
        
    end
    
    
    figure(1);
    hold on; box on;
    
    plot((B_effective/1e6)*SE_MR_DL_maxmin_active,linspace(0,1,Ktotal),MR,'LineWidth',1);
    
    if config == 3
        plot((B_effective/1e6)*SE_MR_DL_maxmin_active,linspace(0,1,Ktotal),'k--','LineWidth',1);
        plot((B_effective/1e6)*SE_MR_DL_maxmin_active,linspace(0,1,Ktotal),'k-.','LineWidth',1);
        legend('M=200, 20 x 5 x 2','M=100, 20 x 5 x 1','M=100, 10 x 5 x 2','Location','SouthEast');
        xlabel('DL throughput per UE [Mbit/s]');
        ylabel('CDF');
        xlim([0 100]);
    end
    
    plot((B_effective/1e6)*SE_RZF_DL_maxmin_active,linspace(0,1,Ktotal),RZF,'LineWidth',1);
    
    
    figure(2);
    hold on; box on;
    
    plot((B_effective/1e6)*SE_MR_DL_prodSINR_active,linspace(0,1,Ktotal),MR,'LineWidth',1);
    
    if config == 3
        plot((B_effective/1e6)*SE_MR_DL_prodSINR_active,linspace(0,1,Ktotal),'k--','LineWidth',1);
        plot((B_effective/1e6)*SE_MR_DL_prodSINR_active,linspace(0,1,Ktotal),'k-.','LineWidth',1);
        legend('M=200, 20 x 5 x 2','M=100, 20 x 5 x 1','M=100, 10 x 5 x 2','Location','SouthEast');
        xlabel('DL throughput per UE [Mbit/s]');
        ylabel('CDF');
        xlim([0 100]);
    end
    
    plot((B_effective/1e6)*SE_RZF_DL_prodSINR_active,linspace(0,1,Ktotal),RZF,'LineWidth',1);
    
    
    figure(3);
    hold on; box on;
    
    plot((B_effective/1e6)*SE_MR_UL_20control_active,linspace(0,1,Ktotal),MR,'LineWidth',1);
    
    if config == 3
        plot((B_effective/1e6)*SE_MR_UL_20control_active,linspace(0,1,Ktotal),'k--','LineWidth',1);
        plot((B_effective/1e6)*SE_MR_UL_20control_active,linspace(0,1,Ktotal),'k-.','LineWidth',1);
        legend('M=200, 20 x 5 x 2','M=100, 20 x 5 x 1','M=100, 10 x 5 x 2','Location','SouthEast');
        xlabel('UL throughput per UE [Mbit/s]');
        ylabel('CDF');
        xlim([0 50]);
    end
    
    plot((B_effective/1e6)*SE_RZF_UL_20control_active,linspace(0,1,Ktotal),RZF,'LineWidth',1);
    
    
    figure(4);
    hold on; box on;
    
    plot((B_effective/1e6)*SE_MR_UL_0control_active,linspace(0,1,Ktotal),MR,'LineWidth',1);
    
    if config == 3
        plot((B_effective/1e6)*SE_MR_UL_0control_active,linspace(0,1,Ktotal),'k--','LineWidth',1);
        plot((B_effective/1e6)*SE_MR_UL_0control_active,linspace(0,1,Ktotal),'k-.','LineWidth',1);
        legend('M=200, 20 x 5 x 2','M=100, 20 x 5 x 1','M=100, 10 x 5 x 2','Location','SouthEast');
        xlabel('UL throughput per UE [Mbit/s]');
        ylabel('CDF');
        xlim([0 50]);
    end
    
    plot((B_effective/1e6)*SE_RZF_UL_0control_active,linspace(0,1,Ktotal),RZF,'LineWidth',1);
    
end
