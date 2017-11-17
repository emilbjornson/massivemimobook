%This Matlab script can be used to reproduce Figure 5.9b in the monograph:
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
load section5_M100_Kvarying;

%Fractions of data samples used for UL and DL
ULfraction = 1/3;
DLfraction = 2/3;

%Compute joint UL/DL sum SE using the fractions of UL/DL data
sumSE_MMMSE = ULfraction*sumSE_MMMSE_UL + DLfraction*sumSE_MMMSE_DL;
sumSE_SMMSE = ULfraction*sumSE_SMMSE_UL + DLfraction*sumSE_SMMSE_DL;
sumSE_RZF = ULfraction*sumSE_RZF_UL + DLfraction*sumSE_RZF_DL;
sumSE_ZF = ULfraction*sumSE_ZF_UL + DLfraction*sumSE_ZF_DL;
sumSE_MR = ULfraction*sumSE_MR_UL + DLfraction*sumSE_MR_DL;

%Number of BSs
L = 16;

%Number of BS antennas
M = 100;

%Communication bandwidth
B = 20e6;

%Select length of coherence block
tau_c = 200;

%Select length of pilot sequences
tau_p = Krange;

%Prepare to store simulation results
P_MMMSE = zeros(length(Krange),2);
P_SMMSE = zeros(length(Krange),2);
P_RZF = zeros(length(Krange),2);
P_ZF = zeros(length(Krange),2);
P_MR = zeros(length(Krange),2);


%% Compute CP values for different value sets
for valueset = 1:2
    
    for k = 1:length(Krange)
        
        %Compute the total CP with different schemes
        [P_MR(k,valueset),P_RZF(k,valueset),P_MMMSE(k,valueset),P_ZF(k,valueset),P_SMMSE(k,valueset)] = functionCPcomputation(M,Krange(k),L,B,tau_c,tau_p(k),valueset,sumSE_MR,sumSE_RZF,sumSE_MMMSE,sumSE_ZF,sumSE_SMMSE);

    end

end


%% Plot the simulation results
figure;
hold on; box on;

for valueset = 1:2
    
    plot(Krange,10*log10(P_MMMSE(:,valueset)/0.001),'rd-','LineWidth',1);hold on;
    plot(Krange,10*log10(P_SMMSE(:,valueset)/0.001),'b:','LineWidth',1);
    plot(Krange,10*log10(P_RZF(:,valueset)/0.001),'k-.','LineWidth',1);
    plot(Krange,10*log10(P_ZF(:,valueset)/0.001),'r--','LineWidth',1);
    plot(Krange,10*log10(P_MR(:,valueset)/0.001),'bs-','LineWidth',1);
    
end

xlabel('Number of UEs (K)');
ylabel('Total CP [dBm]');
legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
