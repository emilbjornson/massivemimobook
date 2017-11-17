function [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL_impairments(H,Hhat,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p,kappatUE,kapparBS)
%Compute UL SE for different receive combining schemes with hardware
%impairments, based on Theorem 6.2.
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact channel realizations
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%kappatUE          = Hardware quality of the UEs' transmitters
%kapparBS          = Hardware quality of the BSs' receivers
%
%OUTPUT:
%SE_MR    = K x L matrix where element (k,l) is the uplink SE of UE k in
%           cell l achieved with MR combining
%SE_RZF   = Same as SE_MR but with RZF combining
%SE_MMMSE = Same as SE_MR but with M-MMSE combining
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


%Store identity matrices of different sizes
eyeK = eye(K);
eyeM = eye(M);

%Compute sum of all estimation error correlation matrices at every BS
C_totM = reshape(p*sum(sum(C,3),4),[M M L]);

%Compute the pre-log factor assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L);
signal_RZF = zeros(K,L);
signal_MMMSE = zeros(K,L);

%Prepare to store simulation results for norms of combining vectors
combiningNorm_MR = zeros(K,L);
combiningNorm_RZF = zeros(K,L);
combiningNorm_MMMSE = zeros(K,L);

%Prepare to store simulation results for sum interference powers
interf_MR = zeros(K,L);
interf_RZF = zeros(K,L);
interf_MMMSE = zeros(K,L);

%Prepare to store simulation results for sum impairment-caused interference
impair_MR = zeros(K,L);
impair_RZF = zeros(K,L);
impair_MMMSE = zeros(K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all users to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        if nargout > 1 %Compute RZF combining in (4.9)
            V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        end
        
        if nargout > 2 %Compute M-MMSE combining in (4.7)
            V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        end
        
        
        %Go through all users in cell j
        for k = 1:K
            
            %Check if the UE is active
            if norm(V_MR(:,k))>0
                
                %%MR combining
                w = V_MR(:,k)/norm(V_MR(:,k))^2; %Extract combining vector
                wrep = repmat(w,[1 K*L]);
                
                %Compute realizations of the expectations in signal and
                %interference terms in (6.32)
                signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                combiningNorm_MR(k,j) = combiningNorm_MR(k,j) + norm(w).^2/nbrOfRealizations;
                interf_MR(k,j) = interf_MR(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
                impair_MR(k,j) = impair_MR(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
                
                
                %%RZF combining
                if nargout >=2
                    
                    w = V_RZF(:,k); %Extract combining vector
                    wrep = repmat(w,[1 K*L]);
                    
                    %Compute realizations of the expectations in signal and
                    %interference terms in (6.32)
                    signal_RZF(k,j) = signal_RZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    combiningNorm_RZF(k,j) = combiningNorm_RZF(k,j) + norm(w).^2/nbrOfRealizations;
                    interf_RZF(k,j) = interf_RZF(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
                    impair_RZF(k,j) = impair_RZF(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
                    
                end
                
                
                %%M-MMSE combining
                if nargout >=3
                    
                    w = V_MMMSE(:,k); %Extract combining vector
                    wrep = repmat(w,[1 K*L]);
                    
                    %Compute realizations of the expectations in signal and
                    %interference terms in (6.32)
                    signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    combiningNorm_MMMSE(k,j) = combiningNorm_MMMSE(k,j) + norm(w).^2/nbrOfRealizations;
                    interf_MMMSE(k,j) = interf_MMMSE(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
                    impair_MMMSE(k,j) = impair_MMMSE(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
                    
                end
                
                %Take care of special cases when the UE is not active
            else
                
                
                %%MR combining
                combiningNorm_MR(k,j) = 1;
                
                %%RZF combining
                if nargout >=2
                    
                    combiningNorm_RZF(k,j) = 1;
                    
                end
                
                %%MMMSE combining
                if nargout >=3
                    
                    combiningNorm_MMMSE(k,j) = 1;
                    
                end
                
                
            end
            
        end
        
    end
    
end

%Precompute terms that appear multiple times
sigma2 = 1/(kappatUE*kapparBS);
factor1 = 1/kappatUE;
factor2 = (1-kapparBS)/(kappatUE*kapparBS);

%Compute SEs according to Theorem 6.2 and (6.32)
SE_MR = prelogFactor*real(log2(1+(p*abs(signal_MR).^2) ./ (factor1*interf_MR + factor2*impair_MR - p*abs(signal_MR).^2 + sigma2*combiningNorm_MR)));

if nargout >=2
    SE_RZF = prelogFactor*real(log2(1+(p*abs(signal_RZF).^2) ./ (factor1*interf_RZF + factor2*impair_RZF - p*abs(signal_RZF).^2 + sigma2*combiningNorm_RZF)));
end

if nargout >=3
    SE_MMMSE = prelogFactor*real(log2(1+(p*abs(signal_MMMSE).^2) ./ (factor1*interf_MMMSE + factor2*impair_MMMSE - p*abs(signal_MMMSE).^2 +sigma2*combiningNorm_MMMSE)));
end
