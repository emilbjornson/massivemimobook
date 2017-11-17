%This Matlab script can be used to reproduce Figure 3.8 in the monograph:
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


%Define number of UEs per cell
K = 10;

%Define the range of BS antennas
Mrange = 10:10:100;

%Define the pilot length
tau_p = K;

%Complexity of MMSE channel estimation, based on Table 3.1
complexity_MMSE = Mrange*tau_p + Mrange.^2;

%Complexity of EW-MMSE channel estimation, based on Table 3.1
complexity_EW_MMSE = Mrange*(tau_p + 1);

%Complexity of LS channel estimation, based on Table 3.1
complexity_LS = Mrange*tau_p;


%% Plot the simulation results
figure;
hold on; box on;

plot(Mrange,complexity_MMSE,'r--','LineWidth',1);
plot(Mrange,complexity_EW_MMSE,'k-','LineWidth',1);
plot(Mrange,complexity_LS,'b-.','LineWidth',1);

xlabel('Number of antennas (M)');
ylabel('Number of complex multiplications');
set(gca,'YScale','log');

legend('MMSE','EW-MMSE','LS','Location','NorthWest');
