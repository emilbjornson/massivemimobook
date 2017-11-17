%This Matlab script can be used to reproduce Figure 5.13 in the monograph:
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


%Load SE simulation data, generated using the code from Section 4
load section5_Mvarying_K10_20;

%Number of UEs per BS
k_index = 2; %Selecting K = 10 from the loaded SE results
K = Krange(k_index);

%Fractions of data samples used for UL and DL
ULfraction = 1/3;
DLfraction = 2/3;

%Compute joint UL/DL sum SE using the fractions of UL/DL data
sumSE_MMMSE = ULfraction*sumSE_MMMSE_UL(:,k_index) + DLfraction*sumSE_MMMSE_DL(:,k_index);
sumSE_SMMSE = ULfraction*sumSE_SMMSE_UL(:,k_index) + DLfraction*sumSE_SMMSE_DL(:,k_index);
sumSE_RZF = ULfraction*sumSE_RZF_UL(:,k_index) + DLfraction*sumSE_RZF_DL(:,k_index);
sumSE_ZF = ULfraction*sumSE_ZF_UL(:,k_index) + DLfraction*sumSE_ZF_DL(:,k_index);
sumSE_MR = ULfraction*sumSE_MR_UL(:,k_index) + DLfraction*sumSE_MR_DL(:,k_index);

%Number of BSs
L = 16;

%Communication bandwidth
B = 20e6;

%PA efficiency UEs and BSs
mu_UE = 0.4;
mu_BS = 0.5;

%Define the pilot reuse factor
f = 1;

%Select length of coherence block
tau_c = 200;

%Compute length of pilot sequences
tau_p = f*K;

%Transmit power per UE in W
p = 0.1;

%Compute total effective transmit power
ETP_total =  K*p*(tau_p/mu_UE + (tau_c-tau_p)*(ULfraction/mu_UE + DLfraction/mu_BS))/tau_c;


%% Go through the two value sets of the CP model
for valueset = 1:2
    
    %Compute the total CP with different schemes
    [P_MR,P_RZF,P_MMMSE,P_ZF,P_SMMSE] = functionCPcomputation(Mrange,K,L,B,tau_c,tau_p,valueset,sumSE_MR,sumSE_RZF,sumSE_MMMSE,sumSE_ZF,sumSE_SMMSE);
    
    
    %Compute EE with M-MMSE
    EE_MMMSE = (B*sumSE_MMMSE)./(ETP_total + P_MMMSE);
    
    %Compute EE with S-MMSE
    EE_SMMSE = (B*sumSE_SMMSE)./(ETP_total + P_SMMSE);
    
    %Compute EE with RZF
    EE_RZF = (B*sumSE_RZF)./(ETP_total + P_RZF);
    
    %Compute EE with ZF
    EE_ZF = (B*sumSE_ZF)./(ETP_total + P_ZF);
    
    %Compute EE with MR
    EE_MR = (B*sumSE_MR)./(ETP_total + P_MR);
    
    
    %Plot simulation results
    figure;
    
    plot(B*sumSE_MMMSE/10^6,EE_MMMSE/10^6, 'rd-','LineWidth',1);hold on;
    plot(B*sumSE_SMMSE/10^6,EE_SMMSE/10^6,'b:','LineWidth',1);
    plot(B*sumSE_RZF/10^6,EE_RZF/10^6,'k-.','LineWidth',1);
    plot(B*sumSE_ZF/10^6,EE_ZF/10^6,'r--','LineWidth',1);
    plot(B*sumSE_MR/10^6,EE_MR/10^6,'bs-','LineWidth',1);
    
    xlabel('Throughput [Mbit/s/cell]');
    ylabel('EE [Mbit/Joule/cell]');
    legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
    
end
