%This Matlab script can be used to reproduce Figure 7.16 in the monograph:
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
M_max = 12; %Largest number of antennas in one dimension
d_H = 1/2; %Horizontal antenna spacing
d_V = 1/2; %Vertical antenna spacing
h_BS = 25; %BS height
h_UE = 1.5; %UE height
r = 50; %Scatter radius
p = 6; %Select dimension of p-dominant eigenspace

%UE position
d = 200; %UE horizontal distance
varphiRange = [0,pi/180*15,pi/180*30,pi/180*45,pi/180*60]; %UE azimuth angles
POS = d*[cos(varphiRange); sin(varphiRange)]'; %Resulting UE position


%% Compute correlation matrix for all parameter combinations
h = h_BS - h_UE;
R = cell(size(POS,1),M_max);
Up = cell(size(POS,1),M_max);

for i = 1:size(POS,1)
    
    %Compute correlation matrix for largest antenna geometry
    M_H = M_max; %Number of antenna per horizontal row
    M_V = M_max; %Number of rows
    M = M_H*M_V; %Number of antennas
    varphi = varphiRange(i); %Azimuth angle
    dvarphi = atan(r/d); %Horizontal angular spread
    theta_max = atan(h/(max(d-r,0))); %Maximum elevation angle
    theta_min = atan(h/(d+r)); %Minimum elevation angle
    theta = (theta_max + theta_min)/2;  %Mean elevation angle
    dtheta = (theta_max - theta_min)/2; %Vertical spread
    R{i,M_max} = functionRlocalscattering3D(M_H, M_V, d_H, d_V, varphi, dvarphi, theta, dtheta, 'Uniform');
    [Up{i,M_max},~,~] = svds(R{i,M_max},min(p,M));
    
    %Compute correlation matrices for smaller antenna geometries by
    %selecting the corresponding antenna indices
    for m = 1:M_max-1
        M_H = m; %Number of antenna per horizontal row
        M_V = m; %Number of rows
        M = M_H*M_V; %Number of antennas
        indices = repmat([0:M_H-1]*M_max + 1, [M_V,1]) + repmat([0:M_V-1]', [1,M_H]);
        indices = indices(:);
        R{i,m} = R{i,M_max}(indices,indices);
        [U,~,~] = svd(R{i,m});
        [Up{i,m},~,~] = svds(R{i,m},min(p,M));
    end
    
end

%Compute the chordal distance of the p-dominant eigenspaces using
%Definition 7.2
dc = @(X,Y) norm(X*X' - Y*Y','fro')^2;
Dc = zeros(size(POS,1)-1,M_max);
for m=1:M_max
    for i=1:size(POS,1)-1
        Dc(i,m) = dc(Up{1,m},Up{i+1,m});
    end
end


%% Plot the simulation results
figure;
hold on; box on;

TMP = Dc;
TMP(1,:) = Dc(end,:);
TMP(2,:) = Dc(end-1,:);
TMP(3,:) = Dc(end-2,:);
TMP(4,:) = Dc(end-3,:);


plot((3:M_max).^2, TMP(1,3:M_max)','bs-','LineWidth', 1);
plot((3:M_max).^2, TMP(2,3:M_max)','r^--','LineWidth', 1);
plot((3:M_max).^2, TMP(3,3:M_max)','ko:','LineWidth', 1);
plot((3:M_max).^2, TMP(4,3:M_max)','bd-.','LineWidth', 1);

xlabel('Number of antennas (M)')
ylabel('Chordal distance')
xlim([9,M_max^2])
ylim([0, 2*p])
legend({'$\varphi=60^\circ$','$\varphi=45^\circ$','$\varphi=30^\circ$','$\varphi=15^\circ$'}, 'Interpreter', 'latex', 'Location', 'SouthEast')
