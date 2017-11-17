%This Matlab script can be used to reproduce Figure 6.3 in the monograph:
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


%Number of BSs
L = 16;

%Number of UEs per BS
K = 10;

%Number of BS antennas
M = 100;

%Define the range of pilot reuse factors
fRange = [1 2 4];

%Select the number of setups with random UE locations
nbrOfSetups = 100;

%Select the number of channel realizations per setup
%(no realizations since only the error correlation matrices are used)
nbrOfRealizations = 0;

%Select the range of hardware qualities. The values at the same position in
%the different vectors are considered simultaneously.
kappatUE = [0.95 1];
kapparBS = kappatUE;


%% Propagation parameters

%Communication bandwidth
B = 20e6;

%Define total uplink transmit power per UE (mW)
p = 100;

%Define noise figure at BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;


%Prepare to save simulation results
NMSE_nonideal = zeros(L,K,nbrOfSetups,length(fRange),length(kappatUE));


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    %Compute channel statistics for one setup
    [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg);
    
    %Compute the normalized average channel gain, where the normalization
    %is based on the noise power
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    
    
    %Go through all pilot reuse factors
    for s = 1:length(fRange)
        
        %Output simulation progress
        disp([num2str(s) ' reuse factors out of ' num2str(length(fRange))]);
        
        %Extract pilot reuse factor
        f = fRange(s);
        
        
        %Go through all hardware quality values
        for r = 1:length(kappatUE)
            
            %Generate channel realizations with estimates and estimation
            %error correlation matrices
            [~,C,~,Rscaled] = functionChannelEstimates_impairments(R,channelGainOverNoise,nbrOfRealizations,M,K,L,p,f,kappatUE(r),kapparBS(r));
            
            %Go through all cells
            for j = 1:L
                
                %Go through all UEs
                for k = 1:K
                    
                    %Save the NMSE
                    NMSE_nonideal(j,k,n,s,r) = real(trace(C(:,:,k,j,j)))/trace(Rscaled(:,:,k,j,j));
                    
                end
                
            end
            
            %Delete large matrices
            clear C Rscaled;
            
        end
        
    end
    
    %Delete large matrices
    clear R;
    
end


%% Plot the simulation results
figure;
hold on; box on;

for r = 1:length(kappatUE)
    
    NMSEcdf = NMSE_nonideal(:,:,:,1,r);
    plot(sort(NMSEcdf(:)),linspace(0,1,K*L*nbrOfSetups),'b--','LineWidth',1);
    
    NMSEcdf = NMSE_nonideal(:,:,:,2,r);
    plot(sort(NMSEcdf(:)),linspace(0,1,K*L*nbrOfSetups),'k-.','LineWidth',1);
    
    NMSEcdf = NMSE_nonideal(:,:,:,3,r);
    plot(sort(NMSEcdf(:)),linspace(0,1,K*L*nbrOfSetups),'r','LineWidth',1);
    
end

set(gca,'XScale','log');
xlim([1e-5 1]);
xlabel('NMSE');
ylabel('CDF')
legend('f=1','f=2','f=4','Location','NorthWest');
