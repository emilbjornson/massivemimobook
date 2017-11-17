%This Matlab script can be used to reproduce Figure 7.18 in the monograph:
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
M_max = 15;          % number of antenna columns/rows
d_H = 1/2;          % horizontal antenna spacing
h_BS = 25;          % height of the antenna array
h_UE = 1.5;         % height of the UE
CellRadius = 500;   % cell radius
K = 500;            % number of UEs
fc = 2.6e9;         % center frequency
bw = 20e6;          % bandwidth
N = 1024;           % number of subcarriers

% Compute derived parameters
h = h_BS - h_UE;    % effective height
d_V = 1/2;          % vertical antenna spacing
d = d_V;

%Simulate new channel reponses or load existing results (the latter is the
%preset)
simulate = false;


%% Simulate channels
if (simulate)
    % Draw random positions
    A = rand(1,K)*CellRadius^2;
    B = rand(1,K)*pi*2/3 - pi/3;
    POS = repmat(sqrt(A)',[1,2]).* [cos(B); sin(B)]'; %UE positions uniformly distributed on a disc of radius R
    
    % Create a new layout
    s = simulation_parameters; % Create object with general parameters
    s.center_frequency = fc;
    s.sample_density = 1;
    s.use_absolute_delays = 1;
    
    l = layout(s);  % Create new layout from general parameters
    l.no_rx = K;  % Generate KUEs
    l.rx_array.generate('omni'); % use omnidirectional antennas at the UEs
    
    l.randomize_rx_positions(CellRadius, h_UE, h_UE, 1); % randomly distribute UEs within a disc of radius R and height 1.5m
    for k=1:K
        l.rx_position(:,k) = [POS(k,1); POS(k,2);  h_UE]; % use the same positions as for the other simulations
    end
    
    for i=1:l.no_rx % for each receiver
        l.track(i).generate('linear',0,0) % define a linear track consiting of only one position
        l.track(i).scenario = '3GPP_3D_UMa_LOS'; % select the Urban Macrocell NLOS scenario for the link
    end
    
    % Simulate channel for planar array
    l.tx_array.generate('3gpp-3d', 1, M_max, M_max, fc, 1, 0, d);
    [h_channel, h_parset, h_cb] = l.get_channels();
    h_freq = h_channel.fr(bw,N);
    H_plan = zeros(K,M_max^2, N);
    for k=1:K
        H_plan(k,:,:) = squeeze(h_freq{k});
    end
    
    % Simulate channel for horizontal array
    l.tx_array.generate('3gpp-3d', 1, 1, M_max^2, fc, 1, 0, d);
    [h_channel, h_parset, h_cb] = l.get_channels();
    h_freq = h_channel.fr(bw,N);
    H_hori = zeros(K,M_max^2, N);
    for k=1:K
        H_hori(k,:,:) = squeeze(h_freq{k});
    end
    
    % Simulate channel for vertical array
    l.tx_array.generate('3gpp-3d', 1, M_max^2, 1, fc, 1, 0, d);
    [h_channel, h_parset, h_cb] = l.get_channels();
    h_freq = h_channel.fr(bw,N);
    H_vert = zeros(K,M_max^2, N);
    for k=1:K
        H_vert(k,:,:) = squeeze(h_freq{k});
    end
end

%% Compute channel orthogonality
if (~simulate)
    load('section7_orthogonality_3GPP.mat')
else
    
    IT = 100000;
    COR_3GPP_PLAN = zeros(1,M_max);
    COR_3GPP_HORI = zeros(1,M_max);
    COR_3GPP_VERT = zeros(1,M_max);
    
    for it=1:IT
        
        %Output simulation progress
        disp([num2str(it) ' iterations out of ' num2str(IT)]);
        
        % Draw random position indices
        posInd = randperm(K,2);
        
        % For all antenna size
        for m = 1:M_max
            M_H = m;
            M_V = m;
            M = M_H*M_V;
            indices = repmat((0:M_H-1)*M_max + 1, [M_V,1]) + repmat((0:M_V-1)', [1,M_H]);
            indices = indices(:);
            for n=1:N
                % Planar array
                h1 = H_plan(posInd(1),indices,n);
                h2 = H_plan(posInd(2),indices,n);
                COR_3GPP_PLAN(m) = abs(h1*h2')^2 / (norm(h1)^2*norm(h2)^2)/IT/N + COR_3GPP_PLAN(m);
                
                % Horizontal array
                h1 = H_hori(posInd(1),1:m^2,n);
                h2 = H_hori(posInd(2),1:m^2,n);
                COR_3GPP_HORI(m) = abs(h1*h2')^2 / (norm(h1)^2*norm(h2)^2)/IT/N + COR_3GPP_HORI(m);
                
                % Vertical array
                h1 = H_vert(posInd(1),1:m^2,n);
                h2 = H_vert(posInd(2),1:m^2,n);
                COR_3GPP_VERT(m) = abs(h1*h2')^2 / (norm(h1)^2*norm(h2)^2)/IT/N + COR_3GPP_VERT(m);
            end
        end
    end
end


%% Load channel measurements
hor_out=load('section7_measurement_horizontal.mat');
vert_out=load('section7_measurement_vertical.mat');
planar_square = load('section7_measurement_planar.mat');

x_measured=(1:8).^2;

%Compute average UE correlation
COR_MEAS_HORI = mean(hor_out.crossco).^2;
COR_MEAS_PLAN = mean(planar_square.crossco).^2;
COR_MEAS_VERT = mean(vert_out.crossco).^2;


%% Plot the simulation results
figure;
hold on; box on;

xvals = (1:M_max).^2;

plot(xvals, COR_3GPP_VERT, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-', 'Marker', 'x');
plot(x_measured, COR_MEAS_VERT(x_measured), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--', 'Marker', 'x');

plot(xvals, COR_3GPP_PLAN, 'LineWidth', 1, 'Color', 'b', 'LineStyle', '-', 'Marker', 'o');
plot(x_measured, COR_MEAS_PLAN, 'LineWidth', 1, 'Color', 'b', 'LineStyle', '--', 'Marker', 'o');

plot(x_measured, COR_MEAS_HORI(x_measured), 'LineWidth', 1, 'Color', 'r', 'LineStyle', '--','Marker', 'square');
plot(xvals, COR_3GPP_HORI, 'LineWidth', 1, 'Color', 'r', 'LineStyle', '-','Marker', 'square');

COR_IID = 1./xvals;
plot(xvals, COR_IID, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-.');

legend({'3GPP Vertical', 'Meas. Vertical', '3GPP Planar', 'Meas. Planar', 'Meas. Horizontal', '3GPP Horizontal', 'i.i.d.'}, 'location', 'southwest');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Number of antennas (M)');
ylabel('Average UE correlation');
set(gca,'xtick',[1,4,9,16,25,36,49,64,121,225]);
xlim([1,M_max^2]);
ylim([COR_IID(end), 1]);
