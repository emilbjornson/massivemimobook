function [signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE,interf_MMMSE,prelogFactor] = functionComputeSINR_DL(H,Hhat,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p)
%Compute DL SE for different transmit precoding schemes using Theorem 4.6.
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact channel realizations
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation.
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%signal_MR    = K x L matrix where element (k,j) is a_jk in (7.2) with MR
%               precoding
%interf_MR    = K x L x K x L matrix where (l,i,j,k) is b_lijk in (7.3)
%               with MR precoding
%signal_RZF   = Same as signal_MR but with RZF precoding
%interf_RZF   = Same as interf_MR but with RZF precoding
%signal_MMMSE = Same as signal_MR but with M-MMSE precoding
%interf_MMMSE = Same as interf_MR but with M-MMSE precoding
%prelogFactor = Prelog factor
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

%Compute the prelog factor assuming only downlink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);


%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L);
signal_RZF = zeros(K,L);
signal_MMMSE = zeros(K,L);

%Prepare to store simulation results for interference powers
interf_MR = zeros(K,L,K,L);
interf_RZF = zeros(K,L,K,L);
interf_MMMSE = zeros(K,L,K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Extract channel realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        %Compute RZF combining in (4.9)
        V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        
        %Compute M-MMSE combining in (4.7)
        V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            if norm(V_MR(:,k))>0
                
                %%MR precoding
                w = V_MR(:,k)/norm(V_MR(:,k)); %Extract precoding vector
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms of (7.2) and (7.3)
                signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                interf_MR(k,j,:,:) = interf_MR(k,j,:,:) + reshape(abs(w'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
                
                
                %%RZF precoding
                w = V_RZF(:,k)/norm(V_RZF(:,k)); %Extract precoding vector
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms of (7.2) and (7.3)
                signal_RZF(k,j) = signal_RZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                interf_RZF(k,j,:,:) = interf_RZF(k,j,:,:) + reshape(abs(w'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
                
                
                %%M-MMSE precoding
                w = V_MMMSE(:,k)/norm(V_MMMSE(:,k)); %Extract precoding vector
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms of (7.2) and (7.3)
                signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                interf_MMMSE(k,j,:,:) = interf_MMMSE(k,j,:,:) + reshape(abs(w'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
                
            end
            
        end
        
    end
    
end

%Compute the terms in (7.2)
signal_MR = abs(signal_MR).^2;
signal_RZF = abs(signal_RZF).^2;
signal_MMMSE = abs(signal_MMMSE).^2;

%Compute the terms in (7.3)
for j = 1:L
    
    for k = 1:K
        
        interf_MR(k,j,k,j) = interf_MR(k,j,k,j) - signal_MR(k,j);
        interf_RZF(k,j,k,j) = interf_RZF(k,j,k,j) - signal_RZF(k,j);
        interf_MMMSE(k,j,k,j) = interf_MMMSE(k,j,k,j) - signal_MMMSE(k,j);
        
    end
    
end
