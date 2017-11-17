function [variation_MR,variation_RZF,variation_MMMSE] = functionComputeGainVariations(H,Hhat,C,nbrOfRealizations,M,K,L,p)
%Compute normalized variance of channel after receive combining, using
%MR, RZF, and M-MMSE receive combining. 
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact realizations of all channels
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates 
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%variation_MR    = K x L matrix with the normalized variance of each UE's
%                  channel after MR combining
%variation_RZF   = Same as variation_MR but with RZF combining
%variation_MMMSE = Same as variation_MR but with M-MMSE combining
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

%Prepare to store simulation results
combinedChannels_MR = zeros(K,L,nbrOfRealizations);
combinedChannels_RZF = zeros(K,L,nbrOfRealizations);
combinedChannels_MMMSE = zeros(K,L,nbrOfRealizations);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        %Compute RZF combining in (4.9)
        V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        
        %Compute M-MMSE combining in (4.7)
        V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            
            %MR combining
            v = V_MR(:,k);
            v = v/norm(v);
            

            %Compute UL channels after receive combining for the
            %interfering channels from all UEs to UE k in cell j.
            combinedChannels_MR(k,j,n) = (v'*H(:,n,k,j,j));
            
            
            %RZF combining
            v = V_RZF(:,k);
            v = v/norm(v);
            
            %Compute UL channels after receive combining for the
            %interfering channels from all UEs to UE k in cell j.
            combinedChannels_RZF(k,j,n) = (v'*H(:,n,k,j,j));
            

            %M-MMSE combining
            v = V_MMMSE(:,k);
            v = v/norm(v);
            
            %Compute UL channels after receive combining for the
            %interfering channels from all UEs to UE k in cell j.
            combinedChannels_MMMSE(k,j,n) = (v'*H(:,n,k,j,j));
            
        end
        
    end
    
end

%Compute normalized variance of channels after receive combining
variation_MR = var(combinedChannels_MR,[],3)./(mean(combinedChannels_MR,3).^2);
variation_RZF = var(combinedChannels_RZF,[],3)./(mean(combinedChannels_RZF,3).^2);
variation_MMMSE = var(combinedChannels_MMMSE,[],3)./(mean(combinedChannels_MMMSE,3).^2);
