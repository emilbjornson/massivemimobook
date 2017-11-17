%This Matlab script can be used to reproduce Figure 7.23 in the monograph:
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


%% Define parameters and plot results
Omegas = -2:0.001:2; %Range of Omega values
S = @(M, d, x) abs(1/M * sin(pi*M*d.*x) ./ sin(pi*d.*x)); %Define S(Omega)
M_V = [8 32]; %Number of antennas per vertical column
d = [1/4 1/2 1]; %Vertical antenna spacings (in number of wavelength)
dstring = {'\frac14', '\frac12', '1'};

figure;
hold on; box on;

for i = 1:numel(d)
    
    for m = 1:numel(M_V)
        
        subplot(numel(d),numel(M_V),(i-1)*numel(M_V)+m);
        plot(Omegas, S(M_V(m),d(i),Omegas), '-k', 'Linewidth', 1);
        xlabel('$\Omega$', 'fontsize', 14, 'interpreter', 'latex');
        ylabel('$S(\Omega)$', 'fontsize', 14, 'interpreter', 'latex');
        set(gca, 'fontsize', 14);
        name = strcat('$M_V=', num2str(M_V(m)), ',\,d_V=', dstring(i), '$');
        title(name, 'fontsize', 14, 'interpreter', 'latex');
        set(gca, 'fontsize', 12);
        
    end
    
end
