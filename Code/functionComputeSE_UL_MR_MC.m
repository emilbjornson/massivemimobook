function [signal,interference,prelogFactor,SE_MR] = functionComputeSE_UL_MR_MC(Hhat,H,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p)
%Use Monte Carlo methods to compute the numerator and denominator of the
%uplink effective SINR in (4.11) with MR, using three different
%normalizations. The SE with Theorem 4.1 is also computed as a reference.
%
%INPUT:
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates 
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact realizations of all channels
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%signal       = K x L x 3 matrix with numerator in (4.14) for all UEs in
%               the network. The result with MR from (4.11) is given by
%               signal(:,:,1), MR as in (4.23) is given by signal(:,:,2),
%               and MR as in (4.24) is given by signal(:,:,3)
%interference = K x L x 3 matrix with denominator in (4.14) for all UEs in
%               the network. The result with MR from (4.11) is given by
%               signal(:,:,1), MR as in (4.23) is given by signal(:,:,2),
%               and MR as in (4.24) is given by signal(:,:,3)
%prelogFactor = Pre-log factor of the SE expression
%SE_MR        = K x L matrix where element (k,l) is the uplink SE of UE k 
%               in cell l achieved with MR combining using Theorem 4.1
%
%
%This Matlab function was developed to generate simulation results to:
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


%Store M x M identity matrix
eyeM = eye(M);

%Compute sum of all estimation error correlation matrices at every BS
C_totM = reshape(p*sum(sum(C(:,:,:,:,:),3),4),[M M L]);

%Compute the pre-log factor assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/tau_c;

%Prepare to store simulation results
SE_MR = zeros(K,L);
combinedChannel = zeros(K,L,nbrOfRealizations,3);
sumpower = zeros(K,L,nbrOfRealizations,3);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining with normalization in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        %Go through all UEs in cell j
        for k = 1:K
        
            %MR combining with normalization in (4.11)
            v = V_MR(:,k);
            
            %Compute signal and interference+noise terms using in (4.3)
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            interfnoise = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - signal;
            
            %Compute instantaneous achievable SE for one channel
            %realization, using lower bound in Theorem 4.1
            SE_MR(k,j) = SE_MR(k,j) + prelogFactor*real(log2(1+signal/interfnoise))/nbrOfRealizations;

            
            %Compute a realization of channel after receive combining when
            %using the original MR normalization in (4.11)
            combinedChannel(k,j,n,1) = v'*H(:,n,k,j,j);
            
            %Compute a realization of the sum power of the received signal
            %after receive combining, corresponding to all terms in the
            %denominator of (4.14), except for the subtraction of the
            %desired signal
            sumpower(k,j,n,1) = p*sum(abs(v'*Hallj).^2) + norm(v)^2;

            
            %MR combining with normalization in (4.23)
            v = V_MR(:,k)/norm(V_MR(:,k));
            
            %Compute a realization of channel after receive combining when
            %using the alternative MR normalization in (4.23)
            combinedChannel(k,j,n,2) = v'*H(:,n,k,j,j);
            
            %Compute a realization of the sum power of the received signal
            %after receive combining, corresponding to all terms in the
            %denominator of (4.14), except for the subtraction of the
            %desired signal
            sumpower(k,j,n,2) = p*sum(abs(v'*Hallj).^2) + norm(v)^2;
            

            %MR combining with normalization in (4.24)
            v = V_MR(:,k)/norm(V_MR(:,k))^2;
            
            %Compute a realization of channel after receive combining when
            %using the alternative MR normalization in (4.24)
            combinedChannel(k,j,n,3) = v'*H(:,n,k,j,j);
            
            %Compute a realization of the sum power of the received signal
            %after receive combining, corresponding to all terms in the
            %denominator of (4.14), except for the subtraction of the
            %desired signal
            sumpower(k,j,n,3) = p*sum(abs(v'*Hallj).^2)  + norm(v)^2;
            
            
        end
        
    end
    
end


%Compute the numerator in (4.14) using Monte Carlo simulations for all UEs
%in the network
signal = squeeze(p*abs(mean(combinedChannel,3)).^2);

%Compute the denominator in (4.14) using Monte Carlo simulations for all
%UEs in the network
interference = squeeze(mean(sumpower,3)) - signal;
