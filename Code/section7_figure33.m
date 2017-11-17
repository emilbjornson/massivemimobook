%This Matlab script can be used to reproduce Figure 7.33 in the monograph:
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
fc0 = 2e9; %Initial carrier frequency
c = 3e8; %Speed of light
lambda0 = c/fc0; %Initial wavelength
A = @(lambda)  (1/8)*lambda.^2; %Function for computing effective area

%Define the function in (7.48)
N = @(lambda,lambda0) lambda./lambda0 .* A(lambda0)./A(lambda);


%% Plot the simulation results
figure;
hold on; box on;

f = fc0:1e9:100e9;
lambda = c./f;
semilogx(f/1e9,N(lambda,lambda0), '-k', 'Linewidth', 1);
xlabel('Carrier frequency $f_c$ [GHz]','interpreter', 'latex');
ylabel('$N(f_c,f_{c_0})$','interpreter', 'latex');
xlim([1,100])
