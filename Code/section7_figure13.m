%This Matlab script can be used to reproduce Figure 7.13 in the monograph:
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
lambda = 1; %Wavelength (normalized)
M_H = 8; %Number of antenna per horizontal row
M_V = 4; %Number of rows
d_H = 0.5*lambda; %Horizontal antenna spacing
d_V = 0.5*lambda; %Vertical antenna spacing

%Define the antenna geometry
M = M_H*M_V; %Total number of antennas
U = zeros(3,M); %Matrix containing the position of the antennas

i = @(m) mod(m-1,M_H); %Horizontal index
j = @(m) floor((m-1)/M_H); %Vertical index

for m = 1:M
    U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the mth element
end


%% Compute the spatial signature for various directions
varphi = linspace(-pi/2,pi/2,1000);
theta= [0,-pi/3];
P = zeros(length(varphi),length(theta));
v = ones(M,1);

for i = 1:length(varphi)
    
    for j = 1:length(theta)
        
        P(i,j) = abs(v'*functionSpatialSignature3DLoS(U,varphi(i),theta(j),lambda))/M;
        
    end
    
end

P = P/max(max(P));
P = 10*log10(P);


%% Plot the simulation results
figure;
hold on; box on;
plot(varphi/pi,P(:,1)','-k','LineWidth',1)
plot(varphi/pi,P(:,2)','r--','LineWidth',1)

ylim([-30,0])
ax = gca;
set(ax, 'XTick', [-0.5, -0.25, 0, 0.25, 0.5]);
xlabel('Azimuth angle $\varphi$ in multiples of $\pi$','Interpreter','latex');
ylabel('Normalized array response [dB]');
legend({'$\theta = 0$', '$\theta = -\pi/3$   '},'Interpreter','latex');
