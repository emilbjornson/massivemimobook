%This Matlab script can be used to reproduce Figure 4.3 in the monograph:
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


%Select length of coherence block
tau_c = 200;

%Define range of number of BS antennas
Mrange = 10:10:100;

%Define range of number of UEs
Krange = 1:40;

%Set number of cells considered in the M-MMSE scheme
L = 9;


%% Consider M=100 and varying K
K = Krange;
M = max(Mrange);

%Compute number of samples for uplink data
tau_u = (tau_c-K);


%Compute complexity of receive combining
receiverProcessing = tau_u.*K*M;

%Add complexity of computing combining matrix
complexity_MMMSE = receiverProcessing + L*K*(M^2+M)/2 + M^2*K + (M^3-M)/3; %M-MMSE
complexity_SMMSE = receiverProcessing + 3*M^2*K/2 + M*K/2 + (M^3-M)/3; %S-MMSE
complexity_RZF = receiverProcessing + 3*K.^2*M/2 + 3*M*K/2 + (K.^3-K)/3; %S-MMSE2
complexity_ZF = receiverProcessing + 3*K.^2*M/2 + M*K/2 + (K.^3-K)/3; %ZF
complexity_MR = receiverProcessing; %MR



%Plot the simulation results for M=100 and varying K
figure(1);
hold on; box on;

plot(K(1),complexity_MMMSE(1),'rd-','LineWidth',1);
plot(K(1),complexity_SMMSE(1),'b:','LineWidth',1);
plot(K(1),complexity_RZF(1),'k-.','LineWidth',1);
plot(K(1),complexity_ZF(1),'r--','LineWidth',1);
plot(K(1),complexity_MR(1),'bs-','LineWidth',1);

plot(K,complexity_MMMSE,'r-','LineWidth',1);
plot(K,complexity_SMMSE,'b:','LineWidth',1);
plot(K,complexity_RZF,'k-.','LineWidth',1);
plot(K,complexity_ZF,'r--','LineWidth',1);
plot(K,complexity_MR,'b-','LineWidth',1);

plot(K([1 5:5:40]),complexity_MMMSE([1 5:5:40]),'rd','LineWidth',1);
plot(K([1 5:5:40]),complexity_MR([1 5:5:40]),'bs','LineWidth',1);

xlabel('Number of UEs (K)');
ylabel('Number of complex multiplications');
set(gca,'YScale','log');

legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','SouthEast');



%% Consider K=10 and varying M
K = Krange(10);
M = Mrange;

%Compute number of samples for uplink data
tau_u = (tau_c-K);


%Compute complexity of receive combining
receiverProcessing = tau_u*K*M;


%Add complexity of computing combining matrix
complexity_MMMSE = receiverProcessing + L*K*(M.^2+M)/2 + M.^2*K + (M.^3-M)/3; %M-MMSE
complexity_SMMSE = receiverProcessing + 3*M.^2*K/2 + M*K/2 + (M.^3-M)/3; %S-MMSE
complexity_RZF = receiverProcessing + 3*K.^2*M/2 + 3*M*K/2 + (K^3-K)/3; %RZF
complexity_ZF = receiverProcessing + 3*K^2*M/2 + M*K/2 + (K^3-K)/3; %ZF
complexity_MR = receiverProcessing; %MR



%Plot the simulation results for K=10 and varying M
figure(2);
hold on; box on;

plot(M(1),complexity_MMMSE(1),'rd-','LineWidth',1);
plot(M(1),complexity_SMMSE(1),'b:','LineWidth',1);
plot(M(1),complexity_RZF(1),'k-.','LineWidth',1);
plot(M(1),complexity_ZF(1),'r--','LineWidth',1);
plot(M(1),complexity_MR(1),'bs-','LineWidth',1);


plot(M,complexity_MMMSE,'r-','LineWidth',1);
plot(M,complexity_SMMSE,'b:','LineWidth',1);
plot(M,complexity_RZF,'k-.','LineWidth',1);
plot(M,complexity_ZF,'r--','LineWidth',1);
plot(M,complexity_MR,'b-','LineWidth',1);

plot(M,complexity_MMMSE,'rd','LineWidth',1);
plot(M,complexity_MR,'bs','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Number of complex multiplications');
set(gca,'YScale','log');

legend('M-MMSE','S-MMSE','RZF','ZF','MR','Location','NorthWest');
