%This Matlab script can be used to reproduce Figure 3.6 in the monograph:
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

%Define range of samples in the sample correlation matrices
sampleRange = [1 10 15 20:10:90 100:50:500];

%Number of realizations in the Monte Carlo simulation (per angle)
numberOfRealizations = 50;

%Angular standard deviation in the local scattering model (in degrees)
ASD = 10;

%Define the range of nominal angles of arrival
varphiRange = linspace(-pi,+pi,50);

%Range of coefficients in the linear combination between correlation matrix
%estimates
linearCombination = linspace(0,1,100);

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Define the effective SNR in (3.13), including processing gain
SNRdB = 10;
SNR = 10.^(SNRdB/10);


%Preallocate matrices for storing the simulation results
NMSE_ideal = zeros(1,length(varphiRange));
NMSE_samples = zeros(length(sampleRange),length(linearCombination),length(varphiRange));



%% Go through the range of nominal angles
for r = 1:length(varphiRange)
    
    %Output simulation progress
    disp([num2str(r) ' angles out of ' num2str(length(varphiRange))]);    
    
    %Compute the spatial correlation matrix
    R = functionRlocalscattering(M,varphiRange(r),ASD,antennaSpacing);
    
    %Compute square root of the spatial correlation matrix
    Rsqrt = sqrtm(R);
    
    %Generate uncorrelated channel realizations
    uncorrelatedRealizations = (randn(M,max(sampleRange),numberOfRealizations)+1i*randn(M,max(sampleRange),numberOfRealizations))/sqrt(2);
    
    
    %Go through all Monte Carlo realizations
    for s = 1:numberOfRealizations
        
        %Generate correlated channel realizations
        channelRealizations = Rsqrt*uncorrelatedRealizations(:,:,s);
        
        %Compute the NMSE in (3.20) if the correlation matrix is known
        NMSE_ideal(1,r) = real(trace(R - SNR*R*((SNR*R+eye(M))\R)))/trace(R);
        
        
        %Go through the different number of samples in the correlation
        %estimates
        for n = 1:length(sampleRange)
            
            %Compute sample correlation matrix
            Rsample = (channelRealizations(:,1:sampleRange(n))*channelRealizations(:,1:sampleRange(n))')/sampleRange(n);
            Rsample = (Rsample + Rsample')/2; %Make sure that it is Hermitian
            
            %Extract only the main diagonal of the sample correlation matrix
            RsampleDiag = diag(diag(Rsample));
            
            
            %Go through all different linear combinations of the sample
            %correlation matrix and the diagonalized version of it
            for c = 1:length(linearCombination)
                
                %Compute the linear combination in (3.27)
                Rest = Rsample*linearCombination(c) + RsampleDiag*(1-linearCombination(c));
                
                %Compute the matrix that the received pilot signal is
                %multiplied with in the approximate MMSE estimator
                A = sqrt(SNR)*Rest/(SNR*Rest+eye(M));
                
                %Compute the NMSE in (3.20) with this approximate estimator
                NMSE_samples(n,c,r) = NMSE_samples(n,c,r) + (trace(R) + trace(A*(SNR*R+eye(M))*A') - 2*real(trace(A'*R))*sqrt(SNR))/(trace(R)*numberOfRealizations);
                
            end
            
        end
        
    end
    
end


%Compute the average over the different nominal angles and then find the
%linar combination that minimizes this average NMSE
NMSE_samples_opt = min(mean(NMSE_samples,3),[],2);

%Copmute the NMSE for uncorrelated channels
NMSE_uncorrelated = 1/(SNR+1);


%% Plot the simulation results
figure;
hold on; box on;

plot(sampleRange,real(NMSE_samples_opt),'k-','LineWidth',1);
plot([min(sampleRange) max(sampleRange)],mean(NMSE_ideal)*ones(1,2),'r--','LineWidth',1);
plot([min(sampleRange) max(sampleRange)],mean(NMSE_uncorrelated)*ones(1,2),'b-.','LineWidth',1);

xlabel('Number of samples (N)');
ylabel('NMSE');
set(gca,'YScale','log');
ylim([1e-2 1e-1]);

legend('Estimated correlation matrix','Location','SouthWest');
