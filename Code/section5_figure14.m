%This Matlab script can be used to reproduce Figure 5.14 in the monograph:
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


%Select the value set from Table 5.3 that is considered
valueset = 1;

%Load SE simulation data, generated using the code from Section 4
load section5_Mvarying_Kvarying.mat;

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

%Fractions of data samples used for UL and DL
ULfraction = 1/3;
DLfraction = 2/3;

%Transmit power per UE in W
p = 0.1;

%Compute joint UL/DL sum SE using the fractions of UL/DL data
sumSE_MMMSE = ULfraction*sumSE_MMMSE_UL + DLfraction*sumSE_MMMSE_DL;
sumSE_RZF = ULfraction*sumSE_RZF_UL + DLfraction*sumSE_RZF_DL;
sumSE_MR = ULfraction*sumSE_MR_UL + DLfraction*sumSE_MR_DL;

%Prepare to save simulation results
EE_MR = zeros(length(Mrange),length(Krange));
EE_RZF = zeros(length(Mrange),length(Krange));
EE_MMMSE = zeros(length(Mrange),length(Krange));


%% Go through all number of BS antennas
for m = 1:length(Mrange)
    
    %Go through all number of UEs
    for k = 1:length(Krange)
        
        %Compute length of pilot sequences
        tau_p = f*Krange(k);
        
        %Compute the total CP with different schemes
        [P_MR,P_RZF,P_MMMSE] = functionCPcomputation(Mrange(m),Krange(k),L,B,tau_c,tau_p,valueset,sumSE_MR,sumSE_RZF,sumSE_MMMSE);
        
        
        %Compute total effective transmit power
        ETP_total =  Krange(k)*p*(tau_p/mu_UE + (tau_c-tau_p)*(ULfraction/mu_UE + DLfraction/mu_BS))/tau_c;
        
        
        %Compute EE with MR
        EE_MR(m,k) = (B*sumSE_MR(m,k))./(ETP_total + P_MR);
        
        %Compute EE with RZF
        EE_RZF(m,k) = (B*sumSE_RZF(m,k))./(ETP_total + P_RZF);
        
        %Compute EE with M-MMSE
        EE_MMMSE(m,k) = (B*sumSE_MMMSE(m,k))./(ETP_total + P_MMMSE);
        
    end
    
end


%% Plot simulation results

%Plot Figure 5.14a
figure;
hold on; box on; grid on;

surfc(Krange,Mrange,EE_MMMSE/10^6,'EdgeColor','none');
colormap(autumn);
shading interp
hold on
contour3(Krange,Mrange,EE_MMMSE/1e6,10,'k')
[~,I] = max(EE_MMMSE(:));
[row,col] = ind2sub(size(EE_MMMSE),I); %2D maximizer
Kopt = Krange(col);  %Optimal number of UEs
Mopt = Mrange(row);  %Optimal number of BS antennas
hold on
plot3(Kopt,Mopt,EE_MMMSE(row,col)/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(Kopt,Mopt,min(min(EE_MMMSE))/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');

view([-17 32]);

xlabel('Number of UEs (K)')
ylabel('Number of BS antennas (M)');
zlabel('EE [Mbit/Joule]');


%Plot Figure 5.14b
figure;
hold on; box on; grid on;

surfc(Krange,Mrange,EE_RZF/10^6,'EdgeColor','none');
colormap(autumn);
shading interp
hold on
contour3(Krange,Mrange,EE_RZF/1e6,10,'k')
[~,I] = max(EE_RZF(:));
[row,col] = ind2sub(size(EE_RZF),I); %2D maximizer
Kopt = Krange(col);  %Optimal number of UEs
Mopt = Mrange(row);  %Optimal number of BS antennas
hold on
plot3(Kopt,Mopt,EE_RZF(row,col)/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(Kopt,Mopt,min(min(EE_RZF))/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');

view([-17 32]);

xlabel('Number of UEs (K)')
ylabel('Number of BS antennas (M)');
zlabel('EE [Mbit/Joule]');



%Plot Figure 5.14b
figure;
hold on; box on; grid on;

surfc(Krange,Mrange,EE_MR/10^6,'EdgeColor','none');
colormap(autumn);
shading interp
hold on
contour3(Krange,Mrange,EE_MR/1e6,10,'k')
[EE_max,I] = max(EE_MR(:));
[row,col] = ind2sub(size(EE_MR),I); %2D maximizer
Kopt = Krange(col);  %Optimal number of UEs
Mopt = Mrange(row);  %Optimal number of BS antennas
hold on
plot3(Kopt,Mopt,EE_MR(row,col)/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(Kopt,Mopt,min(min(EE_MR))/1e6,'k*','MarkerSize',16,'MarkerFaceColor','black');

view([-17 32]);

xlabel('Number of UEs (K)')
ylabel('Number of BS antennas (M)');
zlabel('EE [Mbit/Joule]');
