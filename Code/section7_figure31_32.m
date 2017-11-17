%This Matlab script can be used to reproduce Figures 7.31-32 in the monograph:
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

%The total number of antennas is twice the actual value, since we assume
%that every antenna is dual-polarized and later select a subset that
%corresponds to the two different cases illustrated in Figure 7.30
Mtotal = 512; 

%Effective number of antennas that are used in the simulation
Meff = [2 4 8 16 32 64 128 256];


%% Compute covariance matrices for both UEs

%Define the angles of the three UEs
varphi1 = pi/180*30;
varphi2 = pi/180*25;
varphi3 = -pi/180*10;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 5;

%Antenna spacing
antennaSpacing = 1/2;

%Compute spatial correlation matrices using the local scattering model
R1 = functionRlocalscattering(Mtotal/2,varphi1,ASDdeg,antennaSpacing);
R2 = functionRlocalscattering(Mtotal/2,varphi2,ASDdeg,antennaSpacing);
R3 = functionRlocalscattering(Mtotal/2,varphi3,ASDdeg,antennaSpacing);


%% Precompute polarization correlation matrices
F_BS = eye(2); %Value of (7.35) with perfect XPI
m_UE = [1; 0]; %Value of (7.36) with perfect XPI and no rotation
XPDdB = 7; %Set the XPD in dB
q_XPD = 1/(1 + 10^(XPDdB/10)); %Compute q_XPD as in (7.38)
Sigma = [sqrt(1-q_XPD), sqrt(q_XPD); sqrt(q_XPD), sqrt(1-q_XPD)]; %(7.41)
S = repmat(Sigma,[Mtotal/2,1]);
r_p = 0.2; %Set polarization parameters in (7.42)
t_p = 0.2; %Set polarization parameters in (7.42)
C_r = [1, r_p; r_p', 1]; %Defined in (7.42)
C_t = [1, t_p; t_p', 1]; %Defined in (7.42)
RC = C_t^(1/2); %Compute the third term in (7.44)
M_BSL = kron(eye(Mtotal/2),F_BS); %Compute the last term in (7.44)

THETA = @(x) [cos(x), -sin(x); sin(x), cos(x)]; %Define a rotation matrix

%Compute the first term in (7.37)
LC1 = sqrtm(kron(R1,C_r));
LC2 = sqrtm(kron(R2,C_r));
LC3 = sqrtm(kron(R3,C_r));


%% Monte Carlo simulations
numberOfRealizations = 10000;

%Prepare to save simulation results
Cor_unipl2 = zeros(numel(Meff),numberOfRealizations);
Cor_dpl2 = zeros(numel(Meff),numberOfRealizations);
Cor_unipl3 = zeros(numel(Meff),numberOfRealizations);
Cor_dpl3 = zeros(numel(Meff),numberOfRealizations);

Restup1 = cell(numel(Meff),1);
Restdp1 = cell(numel(Meff),1);
Restup3 = cell(numel(Meff),1);
Restdp3 = cell(numel(Meff),1);

for m = 1:numel(Meff)
    Restup1{m} = zeros(Meff(m));
    Restdp1{m} = zeros(Meff(m));
    Restup3{m} = zeros(Meff(m));
    Restdp3{m} = zeros(Meff(m));
end

%Go through all channel realizations
for it = 1:numberOfRealizations
    
    % Generate UE polarization vectors in (7.36) with random phase offsets
    theta_r1 = 2*pi*rand(1);
    theta_r2 = 2*pi*rand(1);
    theta_r3 = 2*pi*rand(1);
    m_UE1 = THETA(theta_r1) * m_UE;
    m_UE2 = THETA(theta_r2) * m_UE;
    m_UE3 = THETA(theta_r3) * m_UE;
    
    % Generate channel vectors for the three UEs using (7.38)
    W1 = 1/sqrt(2)*(randn(Mtotal,2) + 1i*randn(Mtotal,2));
    W2 = 1/sqrt(2)*(randn(Mtotal,2) + 1i*randn(Mtotal,2));
    W3 = 1/sqrt(2)*(randn(Mtotal,2) + 1i*randn(Mtotal,2));
    Z1 = (LC1*W1*RC) .* S;
    Z2 = (LC2*W2*RC) .* S;
    Z3 = (LC3*W3*RC) .* S;
    h1 = M_BSL*Z1*m_UE1;
    h2 = M_BSL*Z2*m_UE2;
    h3 = M_BSL*Z3*m_UE3;
    
    % Compute correlation coefficients
    for m = 1:numel(Meff)
        
        %By selecting every other antenna, we get a uni-polarized array
        uniplv_ind = 1:2:2*Meff(m);
        
        %By selecting the antennas in a particular pattern, we get a
        %dual-polarized array without having co-located antennas
        dpl_ind = ones(1,Meff(m));
        
        for l = 2:Meff(m)
            
            if (~mod(l,2))
                dpl_ind(l) = dpl_ind(l-1) + 3;
            else
                dpl_ind(l) = dpl_ind(l-1) + 1;
            end
            
        end
        
        %Compute the correlation between channels between UE 1 and UE 2
        Cor_unipl2(m,it) = (abs(h1(uniplv_ind)'*h2(uniplv_ind))/norm(h1(uniplv_ind))/norm(h2(uniplv_ind)))^2;
        Cor_dpl2(m,it) = (abs(h1(dpl_ind)'*h2(dpl_ind))/norm(h1(dpl_ind))/norm(h2(dpl_ind)))^2;
        
        %Compute the correlation between channels between UE 1 and UE 3
        Cor_unipl3(m,it) = (abs(h1(uniplv_ind)'*h3(uniplv_ind))/norm(h1(uniplv_ind))/norm(h3(uniplv_ind)))^2;
        Cor_dpl3(m,it) = (abs(h1(dpl_ind)'*h3(dpl_ind))/norm(h1(dpl_ind))/norm(h3(dpl_ind)))^2;
        
        %Compute spatial correlation matrices of UE 1
        Restup1{m} = Restup1{m} + h1(uniplv_ind)*h1(uniplv_ind)'/numberOfRealizations;
        Restdp1{m} = Restdp1{m} + h1(dpl_ind)*h1(dpl_ind)'/numberOfRealizations;
        
        %Compute spatial correlation matrices of UE 1
        Restup3{m} = Restup3{m} + h3(uniplv_ind)*h3(uniplv_ind)'/numberOfRealizations;
        Restdp3{m} = Restdp3{m} + h3(dpl_ind)*h3(dpl_ind)'/numberOfRealizations;
        
    end
    
end

%Compute the average correlation between UE channels, to be plotted in
%Figure 7.31
Avg_Cor_unipl1 = mean(Cor_unipl2,2);
Avg_Cor_dpl1 = mean(Cor_dpl2,2);
Avg_Cor_unipl3 = mean(Cor_unipl3,2);
Avg_Cor_dpl3 = mean(Cor_dpl3,2);

rankup = zeros(numel(Meff),1);
rankdp = zeros(numel(Meff),1);
chdup = zeros(numel(Meff),1);
chddp = zeros(numel(Meff),1);

for m = 1:numel(Meff)
    
    %Compute rank of correlation matrices to be plotted in Figure 7.32a
    rankup(m) = rank(Restup1{m});
    rankdp(m) = rank(Restdp1{m});
    
    %Compute chordal distance between correlation matrices of UE 1 and UE 3
    %to be plotted in Figure 7.32b
    [U1up,~] = svd(Restup1{m});
    U1up = U1up(:,1:rank(Restup1{m}));
    [U3up,~] = svd(Restup3{m}, 'econ');
    U3up = U3up(:,1:rank(Restup3{m}));
    
    chdup(m) = norm(U1up*U1up' - U3up*U3up','fro')^2;
    [U1dp,~] = svd(Restdp1{m}, 'econ');
    U1dp = U1dp(:,1:rank(Restdp1{m}));
    [U3dp,~] = svd(Restdp3{m}, 'econ');
    U3dp = U3dp(:,1:rank(Restdp3{m}));
    chddp(m) = norm(U1dp*U1dp' - U3dp*U3dp','fro')^2;
    
end


%% Plot the simulation results
figure;
hold on; box on;
plot(Meff,Avg_Cor_unipl1, 'k-', 'LineWidth', 1)
plot(Meff,Avg_Cor_dpl1, 'r--', 'LineWidth', 1)
plot(Meff,Avg_Cor_unipl3, 'k-', 'LineWidth', 1)
plot(Meff,Avg_Cor_dpl3, 'r--', 'LineWidth', 1)
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Number of antennas (M)')
ylabel('Average UE correlation')
legend('Uni-polarized','Dual-polarized', 'Location', 'SouthWest')
xlim([2,256])
set(gca,'xtick',[2,4,8,16,32,64,128,256])

figure;
hold on; box on;
plot(Meff,rankdp, 'r--', 'LineWidth', 1)
plot(Meff,rankup, 'k-', 'LineWidth', 1)
xlim([2,256])
ylim([0,180])
legend('Dual-polarized','Uni-polarized', 'Location', 'NorthWest')
set(gca,'xtick',[2,32,64,128,256])
xlabel('Number of antennas (M)')
ylabel('Rank of correlation matrix')

figure;
hold on; box on;
plot(Meff,chdup, 'k-', 'LineWidth', 1)
plot(Meff,chddp, 'r--', 'LineWidth', 1)
legend('Uni-polarized','Dual-polarized', 'Location', 'NorthWest')
xlim([2,256])
ylim([0,180])
set(gca,'xtick',[2,32,64,128,256])
xlabel('Number of antennas (M)')
ylabel('Chordal distance')
