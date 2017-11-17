%This Matlab script can be used to reproduce Figure 7.20 in the monograph:
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
M = 64; %Number of antennas
Mprime = 30; %Number of 
beta = [0,10^(-5/10),1];
IT = 10000;

%Derived parameters
betad = (M-beta*Mprime)/(M-Mprime);


%% Simulate channels
correlationAA = zeros(numel(beta),IT);
correlationAB = zeros(numel(beta),IT);

%Go through all beta values
for b = 1:numel(beta)
    R_A = diag([beta(b)*ones(1,Mprime), betad(b)*ones(1,M-Mprime)]);
    R_B = diag([betad(b)*ones(1,M-Mprime), beta(b)*ones(1,Mprime)]);
    R_A2 = sqrt(R_A);
    R_B2 = sqrt(R_B);
    
    %Compute UE correlation metric in (7.22)
    for it = 1:IT
        
        hA1 = R_A2*1/sqrt(2)*(randn(M,1) + 1i*randn(M,1));
        hA2 = R_A2*1/sqrt(2)*(randn(M,1) + 1i*randn(M,1));
        hB1 = R_B2*1/sqrt(2)*(randn(M,1) + 1i*randn(M,1));
        correlationAA(b,it) = 10*log10((abs(hA1'*hA2)/(norm(hA1)*norm(hA2)))^2);
        correlationAB(b,it) = 10*log10((abs(hA1'*hB1)/(norm(hA1)*norm(hB1)))^2);
        
    end
end


%% Plot the simulation results
figure;
hold on; box on;
nth = 1000;
Colors = {'r','b'};
Markers = {'square', 'o'};

plot(-100,-100,'Color', 'k', 'LineStyle', '-');
plot(-100,-100,'Color', 'k', 'LineStyle', '--');
plot(-100,-100,'Color','k', 'LineStyle', 'none', 'Marker', '*');
plot(-100,-100,'Color', Colors{2}, 'LineStyle', 'none', 'Marker', Markers{2});
plot(-100,-100,'Color', Colors{1}, 'LineStyle', 'none', 'Marker', Markers{1});

[y,x] = ecdf(correlationAA(end,:));
plot(x,y,'Color','k', 'LineWidth', 1);
plot(x(1:nth:end),y(1:nth:end),'Color','k', 'LineStyle', 'none', 'Marker', '*');

for b = 1:numel(beta)-1
    [y,x] = ecdf(correlationAA(b,:));
    plot(x,y,'Color', Colors{b}, 'LineWidth', 1);
end

for b = 1:numel(beta)-1
    [y,x] = ecdf(correlationAB(b,:));
    plot(x,y,'--','Color', Colors{b}, 'LineWidth', 1)
end


for b = 1:numel(beta)-1
    [y,x] = ecdf(correlationAA(b,:));
    plot(x(1:nth:end),y(1:nth:end),'LineStyle', 'none','Color', Colors{b}, 'LineWidth', 1, 'Marker', Markers{b});
end

for b = 1:numel(beta)-1
    [y,x] = ecdf(correlationAB(b,:));
    plot(x(1:nth:end),y(1:nth:end),'LineStyle', 'none','Color', Colors{b}, 'LineWidth', 1, 'Marker', Markers{b})
end

xlim([-40,-10])
legend('Same region', 'Different regions','\beta=1', '\beta=-5 dB', '\beta=0','Location','NorthWest')
xlabel('Average UE correlation [dB]')
ylabel('CDF')
