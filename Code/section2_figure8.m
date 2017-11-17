%This Matlab script can be used to reproduce Figure 2.8 in the monograph:
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


%Number of BS antennas
M = 100;

%Angular standard deviation in the local scattering model (in degrees)
ASDs = [10 30];

%Set the nominal angle of the desired UE
varphiDesired = pi/6;

%Set range of nominal angles of the interfering UE
varphiInterfererDegrees = -180:1:180;
varphiInterfererRadians = varphiInterfererDegrees*(pi/180);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Preallocate matrix for storing the simulation results
variance = zeros(length(varphiInterfererRadians),length(ASDs));


%% Go through the range of ASDs
for n = 1:length(ASDs)

    %Output simulation progress
    disp([num2str(n) ' ASDs out of ' num2str(length(ASDs))]);    
    
    %Compute spatial correlation matrix of the desired UE
    R1 = functionRlocalscattering(M,varphiDesired,ASDs(n),antennaSpacing);
    
    %Go through all angles of the interfering UE
    for r = 1:length(varphiInterfererRadians)
        
        %Compute spatial correlation matrix of the interfering UE
        R2 = functionRlocalscattering(M,varphiInterfererRadians(r),ASDs(n),antennaSpacing);
        
        %Compute variance of favorable propagation according to (2.19)
        variance(r,n) = real(trace(R1*R2)/(trace(R1)*trace(R2)));

    end
    
end


%% Plot the simulation results
figure;
hold on; box on;

plot(varphiInterfererDegrees,variance(:,1),'r--','LineWidth',1);
plot(varphiInterfererDegrees,variance(:,2),'b-.','LineWidth',1);
plot(varphiInterfererDegrees,1/M*ones(length(varphiInterfererDegrees),1),'k-','LineWidth',1);

xlabel('Angle of interfering UE [degree]');
ylabel('Variance in (2.19)');
xlim([-180 180]);

legend('Gaussian, ASD 10^o','Gaussian, ASD 30^o','Uncorrelated','Location','NorthWest');
