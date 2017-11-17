%This Matlab script can be used to reproduce Figure 7.8 in the monograph:
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


%Maximum number of pilots
pilotMax = 40;

%Define range of average number of UEs
Kaverages = [1 10 20 40];

%Define range of effective number of UEs
Krange = 0:pilotMax;

%Prepare to store simulation results
Kdistributions = zeros(length(Krange),length(Kaverages));


%% Go through all average numbers of UEs
for n = 1:length(Kaverages)
    
    %Compute the probabilities of having different numbers of active UEs, 
    %based on a Poission distribution
    Kdistributions(:,n) = poisspdf(Krange,Kaverages(n));
    
    %Enforce scheduling to make sure that there cannot be more than 
    %pilotMax UEs active per cell
    Kdistributions(end,n) = 1 - sum(Kdistributions(1:end-1,n));
    
end


%% Plot simulation results
figure;
hold on; box on;

plot(Krange,Kdistributions(:,1),'r*-','LineWidth',1);
plot(Krange,Kdistributions(:,2),'ko-','LineWidth',1);
plot(Krange,Kdistributions(:,3),'bs-','LineWidth',1);
plot(Krange,Kdistributions(:,4),'rd--','LineWidth',1);

xlabel('Number of active UEs');
ylabel('Probability');
legend('K=1  ','K=10  ','K=20  ','K=40  ','Location','NorthWest');
ylim([0 0.6]);
