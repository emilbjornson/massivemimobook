%This Matlab script can be used to reproduce Figures 7.38-39 in the monograph:
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


%% Parameters
targetRates = [1,2,5,10,15,20:10:100];  %Target rates in Mbit/s
S = [20,40,81];                         %Number of randomly selected SBSs out of 81
L = 16;                                 %Number of BSs
B = 20;                                 %Bandwidth in MHz
noiseFigure = 7;                        %Noise figure at the SBSs and BSs
DLfraction = 2/3;                       %Fraction of time used for DL
nbrOfSetups = 100;                      %Number of random SBS selections


%Compute target SINRs
f_SINR = @(x) 2.^(x/B/DLfraction)-1; % Rate = B*DLfraction*log2(1+targetSINR)
targetSINR = f_SINR(targetRates);

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B*1e6) + noiseFigure;

%Maximum transmit power per BS (in dBm)
Pmax = 46;


%% Monte Carlo simulations
Pmin = zeros(numel(S),numel(targetSINR));
Mmin = zeros(numel(S),numel(targetSINR));

%Go through all numbers of SBSs
for s = 1:numel(S)
    
    %Output simulation progress
    disp([num2str(s) ' small cell numbers out of ' num2str(length(S))]);
    
    %Go through all setups with random SBS selections
    for n = 1:nbrOfSetups
        
        %Output simulation progress
        disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        
        %Compute average channel gains from all SBSs to all BSs
        channelGaindB = functionExampleSetup_backhaul(S(s),L);
        
        %Compute the normalized average channel gain, where the 
        %normalization is based on the noise power
        channelGainOverNoise = channelGaindB - noiseVariancedBm;
        
        %Compute average channel gains in linear scale
        betas = 10.^(channelGainOverNoise/10);
        
        Mstart = 1; % number of antennas to start with
        
        %Go through all target rates
        for i = 1:numel(targetRates)
            
            %Obtain the number of antennas and minimal power
            [M, P] = functionMinAntennasPower(targetSINR(i),betas,Mstart,Pmax);
            Mmin(s,i) = Mmin(s,i) + M/nbrOfSetups;
            Pmin(s,i) = Pmin(s,i) + P/nbrOfSetups;
            Mstart = M; %Increase Mstart because we need at least M for a higher SINR value
            
        end
        
    end
    
end


%% Plot the simulation results

%Minimum number of BS-antennas plot
figure;
hold on; box on;
plot(targetRates, Mmin(3,:), '-r^', 'LineWidth', 1)
plot(targetRates, Mmin(2,:), '-bsquare', 'LineWidth', 1)
plot(targetRates, Mmin(1,:), '-ko', 'LineWidth', 1)
legend({'S=81', 'S=40', 'S=20'}, 'Location', 'NorthWest')
ylabel('Required number of BS antennas')
xlabel('Downlink backhaul rate [Mbit/s]')

first_axis = gca;
sqz = 0.12; %Distance to squeeze the first plot
set(first_axis, 'Position', get(first_axis, 'Position') + [0 sqz 0 -sqz ]);
ax2 = axes('Position', get(first_axis, 'Position') .* [1 1 1 0.001] - [0 sqz 0 0],'Color','none');
xlim(get(first_axis, 'XLim') * 0.5);
set(ax2, 'XScale', get(first_axis, 'XScale'));
xlabel(ax2,'Uplink backhaul rate [Mbit/s]')



%Minimum transmit power plot
figure;
hold on; box on;
plot(targetRates, Pmin(3,:), '-r^', 'LineWidth', 1)
plot(targetRates, Pmin(2,:), '-bsquare', 'LineWidth', 1)
plot(targetRates, Pmin(1,:), '-ko', 'LineWidth', 1)
plot(targetRates, 46*ones(1,numel(targetRates)), '--k', 'LineWidth', 1)
legend({'S=81', 'S=40', 'S=20'}, 'Location', 'SouthEast')
ylabel('Minimum transmit power per BS [dBm]')
xlabel('Downlink backhaul rate [Mbit/s]')

first_axis = gca;
sqz = 0.12; %Distance to squeeze the first plot
set(first_axis, 'Position', get(first_axis, 'Position') + [0 sqz 0 -sqz ]);
ax2 = axes('Position', get(first_axis, 'Position') .* [1 1 1 0.001] - [0 sqz 0 0],'Color','none');
xlim(get(first_axis, 'XLim') * 0.5);
set(ax2, 'XScale', get(first_axis, 'XScale'));
xlabel(ax2,'Uplink backhaul rate [Mbit/s]')
