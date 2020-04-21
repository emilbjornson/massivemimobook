%This Matlab script can be used to reproduce Figure 1.15 in the monograph:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.01 (Last edited: 2020-04-21)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Empty workspace and close figures
close all;
clear;


%Set maximum number of antennas
Mmax = 400;

%Set the range of number of antennas
M = (1:Mmax)';


%Prepare to save simulation results
minval = zeros(Mmax,1);
maxval = zeros(Mmax,1);


%Compute percentiles using the inverse chi-squared distribution
for m = 1:Mmax
    
    minval(m) = chi2inv(0.1,2*m);
    maxval(m) = chi2inv(0.9,2*m);
    
end


%Generate one random channel realization
h = (randn(Mmax,1)+1i*randn(Mmax,1))/sqrt(2);

%Compute the normalized squared norm
h2 = cumsum(abs(h(:,1)).^2)./M;


%% Plot the simulation results
figure;
hold on; box on;
plot(M,ones(Mmax,1),'k-','LineWidth',1);
plot(M,minval./(2*M),'b--','LineWidth',1);
plot(M,h2,'r-','LineWidth',1);
plot(M,maxval./(2*M),'b--','LineWidth',1);

legend('Mean value','Percentiles','One realization');
xlabel('Number of antennas (M)');
ylabel('Normalized channel gain');
ylim([0 3]);
