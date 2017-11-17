%This Matlab script can be used to reproduce Figure 7.26 in the monograph:
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
W = [350, 100]; %Cell radii
numberOfSetups = 1000; %Number of random subarray drops
numberOfRealizations = 200; %Number of channel realizations per subarray drop

%Propagation model parameters as in (2.3)
alpha = 3.76;
Upsilon = -148.1;
sigma_sf = 8;
pathlossdB = @(d) Upsilon-10*alpha*log10(d/1000);
addShadowing = @(pldB) 10.^((pldB +  sigma_sf*randn(size(pldB)))/10);


%% Monte Carlo simulations
metric = zeros(numel(S), numberOfSetups*numberOfRealizations, numel(W));

%Go through all cell radii
for w = 1:numel(W)
    
    %Output simulation progress
    disp([num2str(w) ' cell radii out of ' num2str(length(W))]);
    
    %Go through the number of subarrays
    for s = 1:numel(S)
        
        %Output simulation progress
        disp([num2str(s) ' subarray size out of ' num2str(length(S))]);
        
        %Go through all drops of subarrays
        for it1 = 1:numberOfSetups
            
            d = sqrt(rand(S(s),1)*(W(w)-sqrt(d0)) + sqrt(d0)); %The squared distances of randomly distributed point on a circle are uniformly distributed
            pldB = pathlossdB(d); %Compute pathloss without shadow fading
            pl = addShadowing(pldB); %Add random shadow fading
            R = kron(diag(pl),eye(M(s)));
            R12 = sqrt(R);
            
            %Generate channel realizations and compute the metric
            h = R12 * 1/sqrt(2)*(randn(M(1),numberOfRealizations) +1i*randn(M(1),numberOfRealizations));
            metric(s, (it1-1)*numberOfRealizations+1:it1*numberOfRealizations,w) = sum(abs(h).^2,1)/trace(R);
            
        end
        
    end
    
end

%Compute the variance in (2.17)
variance = zeros(numel(W),numel(S));
for w = 1:numel(W)
    variance(w,:) = var(metric(:,:,w),0,2);
end


%% Plot the simulation results
figure;
hold on; box on;
set(gca,'xtick', S, 'XScale', 'log');
xlim([1,M(1)]);

plot(S, variance(1,:),  'LineWidth', 1, 'Color', 'k', 'LineStyle', '-', 'Marker', 'o')
plot(S, variance(2,:),  'LineWidth', 1, 'Color', 'r', 'LineStyle', '-', 'Marker', '^')
xlabel('Number of subarrays');
ylabel('Variance in (2.17)');

legend({'Cell radius: 350 m','Cell radius: 100 m'}, 'Location', 'NorthWest');
