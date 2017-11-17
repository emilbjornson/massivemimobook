%This Matlab script can be used to reproduce Figure 5.11 in the monograph:
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

%Obtain CP model coefficients for Value set 2 in Table 5.3
[P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(2);

%Only consider M=100 antennas
index = 11;

%Extract the number of antennas
M = Mrange(index);


%% Compute CP values

%Compute CP for coding and decoding using (5.35)
P_CD_MMMSE = (P_COD + P_DEC)*B*sumSE_MMMSE(index);
P_CD_SMMSE = (P_COD + P_DEC)*B*sumSE_SMMSE(index);
P_CD_RZF = (P_COD + P_DEC)*B*sumSE_RZF(index);
P_CD_ZF = (P_COD + P_DEC)*B*sumSE_ZF(index);
P_CD_MR = (P_COD + P_DEC)*B*sumSE_MR(index);

%Compute CP for backhaul traffic using (5.36)
P_BH_MMMSE = P_BT*B*sumSE_MMMSE(index);
P_BH_SMMSE = P_BT*B*sumSE_SMMSE(index);
P_BH_RZF = P_BT*B*sumSE_RZF(index);
P_BH_ZF = P_BT*B*sumSE_ZF(index);
P_BH_MR = P_BT*B*sumSE_MR(index);

%Compute CP for transceiver chains using (5.34)
P_TC = M*P_BS + P_LO + K*P_UE;

%Compute CP for channel estimation with all other schemes, where only
%the channels to UEs in other cells are estimated, using (5.37)
P_CE = 3*K*B/(tau_c*L_BS)*(M*tau_p + M^2);

%Compute CP for UL reception and DL transmission
P_SP_RT = (tau_c - tau_p)*3*B/(tau_c*L_BS)*M*K;

%Compute CP for computation of precoding vectors
P_SP_DL = 4*M*K*B/(tau_c*L_BS);

%Sum up the power terms that are independent of the processing scheme
P_SAME = P_FIX + P_TC + P_SP_RT + P_SP_DL;


%Compute CP for computation of the combining vectors with different
%schemes, based on Table 5.2
P_SP_UL_MMMSE = 3*B*(L*(3*M^2 + M)*K/2 + M^3/3 + 2*M + M*tau_p*(tau_p-K))/(tau_c*L_BS);
P_SP_UL_SMMSE = 3*B*(3*M^2*K/2 + M*K/2 + (M^3 - M)/3 + (7/3)*M)/(tau_c*L_BS);
P_SP_UL_RZF = 3*B*(3*K^2*M/2 + 3*M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
P_SP_UL_ZF = 3*B*(3*K^2*M/2 + M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
P_SP_UL_MR = 7*B*K/(tau_c*L_BS);


%% Plot the simulation results
figure;
power = 10*log10([P_FIX/0.001; P_TC/0.001; P_SP_RT/0.001; P_SP_DL/0.001]);
power = [power , power , power ];

bar(power);
ylabel('Power [dBm]');
set(gca,'xticklabel',['             Fixed power            ';'          Transceiver chain         '; 'Signal processing: Reception/transm.'; '        Computing precoding         ']);
legend('M-MMSE','RZF','MR','Location','NorthEast');
colormap(hot);
ylim([-10 50]);


figure;
P_MMMSE = [ P_CE, P_SP_UL_MMMSE, P_BH_MMMSE, P_CD_MMMSE];
P_RZF = [  P_CE, P_SP_UL_RZF, P_BH_RZF, P_CD_RZF];
P_MR = [ P_CE, P_SP_UL_MR, P_BH_MR, P_CD_MR];

power = 10*log10([P_MMMSE'/0.001, P_RZF'/0.001, P_MR'/0.001]);

bar(power);
ylabel('Power [dBm]');
set(gca,'xticklabel',['    Channel estimation     '; 'Computing receive combining'; '          Backhaul         ';'      Coding/decoding      ']);
legend('M-MMSE','RZF','MR','Location','NorthEast');
colormap(hot);
ylim([-10 50]);
