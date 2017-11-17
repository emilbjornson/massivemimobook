function [channelGains_MR,channelGains_RZF,channelGains_MMMSE] = functionComputeULDLPowerLevels(H,Hhat,C,nbrOfRealizations,M,K,L,p)
%Compute the average squared inner product between combining/precoding
%vectors and the channels of different UEs. These represent the power of
%the desired and interfering signals.
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact channel realizations
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
%channelGains_MR    = K x L x K x L matrix where (k,j,i,l) is the average
%                     squared inner product between the MR
%                     combining/precoding vector of UE k in cell j and the
%                     channel from BS j to UE i in cell l
%channelGains_RZF   = Same as channelGains_MR but with RZF combining/precoding
%channelGains_MMMSE = Same as channelGains_MR but with M-MMSE combining/precoding
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
C_totM = reshape(p*sum(sum(C(:,:,:,:,:),3),4),[M M L]);


%Prepare to store simulation results
channelGains_MR = zeros(K,L,K,L);
channelGains_RZF = zeros(K,L,K,L);
channelGains_MMMSE = zeros(K,L,K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Extract channel estimate realizations from UEs in cell j to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        
        %Compute three different combining/precoding schemes
        V_MR = Hhatallj(:,K*(j-1)+1:K*j); %MR in Eq. (XX)
        V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK); %RZF in Eq. (XX)
        V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR; %M-MMSE in Eq. (XX)
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            
            %%MR combining/precoding
            v = V_MR(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute UL powers for the interfering channels from all UEs to
            %UE k in cell j. These are also the interfering DL gains
            %from UE k in cell j to all other UEs.
            channelGains_MR(k,j,:,:) = channelGains_MR(k,j,:,:) + reshape(abs(v'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
            
            
            %%RZF combining/precoding
            v = V_RZF(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute UL powers for the interfering channels from all UEs to
            %UE k in cell j. These are also the interfering DL gains
            %from UE k in cell j to all other UEs.
            channelGains_RZF(k,j,:,:) = channelGains_RZF(k,j,:,:) + reshape(abs(v'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
            
            
            %%M-MMSE combining/precoding
            v = V_MMMSE(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute UL powers for the interfering channels from all UEs to
            %UE k in cell j. These are also the interfering DL gains
            %from UE k in cell j to all other UEs.
            channelGains_MMMSE(k,j,:,:) = channelGains_MMMSE(k,j,:,:) + reshape(abs(v'*Hallj).^2,[1 1 K L])/nbrOfRealizations;
            
        end
        
    end
    
end
