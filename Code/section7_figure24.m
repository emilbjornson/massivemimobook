%This Matlab script can be used to reproduce Figure 7.24 in the monograph:
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


%% Define parameters
thetas = [75 65 45 0]; %Values of the elevation angle
costhetas = cos(thetas/180*pi);
Psid = -1:0.001:1; %Range of values of Phi'
T = @(M, d, x) abs(1/M * sin(pi*M*d.*x) ./ sin(pi*d.*x)); %Define T(Psi)
M_H = 16; %Number of horizontal antennas
d_H = 1/2; %Horizontal antenna spacing


%% Plot the simulation results
figure
hold on
box on

plot(Psid, T(M_H,d_H,costhetas(1)*Psid), 'k:', 'Linewidth', 1)
plot(Psid, T(M_H,d_H,costhetas(2)*Psid), 'b-.', 'Linewidth', 1)
plot(Psid, T(M_H,d_H,costhetas(3)*Psid), 'r--', 'Linewidth', 1)
plot(Psid, T(M_H,d_H,costhetas(4)*Psid), 'k-', 'Linewidth', 1)

xlabel('\Psi''');
ylabel('T''(\Psi'')');
legend('\theta=75^o','\theta=65^o','\theta=45^o','\theta=0^o')
