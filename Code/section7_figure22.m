%This Matlab script can be used to reproduce Figure 7.22 in the monograph:
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
M_V = 16; %Number of antennas per vertical column
d_V = [2 1 0.5]; %Vertical antenna spacings (in number of wavelength)
f_c = 1:0.25:5; %Carrier frequency (in GHz)
c = 3e8; %Speed of light
lambda_c = c./(f_c*1e9); %Wavelength at carrier frequency
h = M_V*d_V'*lambda_c; %Compute antenna height


%% Plot the simulation results
figure;
hold on; box on;
plot(f_c,h(1,:), '-ko', 'Linewidth', 1)
plot(f_c,h(2,:), '-rs', 'Linewidth', 1)
plot(f_c,h(3,:), '-bd', 'Linewidth', 1)

xlabel('Carrier frequency [GHz]')
ylabel('Array height [m]')
legend('d_V=2', 'd_V=1', 'd_V=1/2')
