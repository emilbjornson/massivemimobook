%This Matlab script can be used to reproduce Figure 7.15 in the monograph:
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
M_H = 8; %Number of antenna per horizontal row
M_V = 8; %Number of rows
d_H = 1/2; %Horizontal antenna spacing
d_V = 1/2; %Vertical antenna spacing
h_BS = 25; %BS height
h_UE = 1.5; %UE height

%Range of scatter radii
r = [50,100,200];

%UE position
d = 200; %UE horizontal distance
varphi = pi/8; %UE azimuth angle
POS = repmat(d,[1,2]).* [cos(varphi); sin(varphi)]'; % Resulting UE position


%% Compute correlation matrix for each scatter radius
h = h_BS - h_UE; %Difference in height
M = M_H*M_V; %Number of antennas
EV = zeros(length(r),M);

for i = 1:length(r)
    dvarphi = atan(r(i)/d);
    theta_max = atan(h/(max(d-r(i),0)));
    theta_min = atan(h/(d+r(i)));
    theta = (theta_max + theta_min)/2;
    dtheta = (theta_max - theta_min)/2;    
    R = functionRlocalscattering3D(M_H, M_V, d_H, d_V, varphi, dvarphi, theta, dtheta, 'Uniform');
    EV(i,:) = sort(abs(eig(R)),'descend');   
end


%% Plot the simulation results
figure;
hold on; box on;

plot(1:M,10*log10(EV(1,:)),'k','LineWidth',1);
plot(1:M,10*log10(EV(2,:)),'b--','LineWidth',1);
plot(1:M,10*log10(EV(3,:)),'r-.','LineWidth',1);

xlabel('Eigenvalue number in decreasing order');
ylabel('Normalized eigenvalue [dB]');
ylim([-40 20]);
legend({'r = 50 m','r = 100 m','r = 200 m'},'Location','NorthEast','Interpreter','latex');
