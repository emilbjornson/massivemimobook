%This Matlab script can be used to reproduce Figure 5.9a in the monograph:
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
k_index = 1; %Selecting K = 10 from the loaded SE results
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

%Select length of coherence block
tau_c = 200;

%Select length of pilot sequences
tau_p = K;

%Prepare to store simulation results
P_MMMSE = zeros(length(Mrange),2);
P_SMMSE = zeros(length(Mrange),2);
P_RZF = zeros(length(Mrange),2);
P_ZF = zeros(length(Mrange),2);
P_MR = zeros(length(Mrange),2);


%% Compute CP values for different value sets
for valueset = 1:2
    
    [P_MR(:,valueset),P_RZF(:,valueset),P_MMMSE(:,valueset),P_ZF(:,valueset),P_SMMSE(:,valueset)] = functionCPcomputation(Mrange,K,L,B,tau_c,tau_p,valueset,sumSE_MR,sumSE_RZF,sumSE_MMMSE,sumSE_ZF,sumSE_SMMSE);
    
end


%% Plot the simulation results
figure;
hold on; box on;

for valueset = 1:2
    
    plot(Mrange,10*log10(P_MMMSE(:,valueset)/0.001),'rd-','LineWidth',1);
    plot(Mrange,10*log10(P_SMMSE(:,valueset)/0.001),'b:','LineWidth',1);
    plot(Mrange,10*log10(P_RZF(:,valueset)/0.001),'k-.','LineWidth',1);
    plot(Mrange,10*log10(P_ZF(:,valueset)/0.001),'r--','LineWidth',1);
    plot(Mrange,10*log10(P_MR(:,valueset)/0.001),'bs-','LineWidth',1);
    
end

xlabel('Number of antennas (M)');
ylabel('Total CP [dBm]');
legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','SouthEast');
