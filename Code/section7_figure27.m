%This Matlab script can be used to reproduce Figure 7.27 in the monograph:
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
S = [1, 2, 4, 8, 16, 32, 64, 128, 256]; %Number of subarrays
M = 256./S; %Number of antennas per subarray
d0 = 10; %Guard distance to the closest subarray
W = 350; %Cell radii
numberOfSetups = 5000; %Number of random subarray drops

%Propagation model parameters as in (2.3)
alpha = 3.76;
Upsilon = -148.1;
sigma_sf = 8;
pathlossdB = @(d) Upsilon-10*alpha*log10(d/1000);
addShadowing = @(pldB) 10.^((pldB +  sigma_sf*randn(size(pldB)))/10);


%% Monte Carlo simulations
betaValues = zeros(numel(S), numberOfSetups);

%Go through the number of subarrays
for s = 1:numel(S)
    
    %Go through all drops of subarrays
    for it1 = 1:numberOfSetups
        
        d = sort(sqrt(rand(S(s),1)*(W-sqrt(d0)) + sqrt(d0))); %The squared distances of randomly distributed point on a circle are uniformly distributed
        pldB = pathlossdB(d); %Compute pathloss without shadow fading
        pl = addShadowing(pldB); %Add random shadow fading
        
        %Compute average channel gain
        betaValues(s,it1) = sum(pl)/S(s);
        
    end
    
end


%% Plot the simulation results
figure;
hold on; box on;
plot(sort(10*log10(betaValues(1,:))),linspace(0,1,numberOfSetups), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
plot(sort(10*log10(betaValues(5,:))),linspace(0,1,numberOfSetups), 'LineStyle', '--', 'Color', 'r',  'LineWidth', 1);
plot(sort(10*log10(betaValues(9,:))),linspace(0,1,numberOfSetups), 'LineStyle', '-.', 'Color', 'b',  'LineWidth', 1);

xlabel('Average channel gain [dB]')
ylabel('CDF')
legend({['S=', num2str(S(1))],['S=', num2str(S(5))],['S=', num2str(S(9))]}, 'Location', 'NorthWest')
