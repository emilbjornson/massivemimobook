function [powerterms_MR,powerterms_RZF,powerterms_MMMSE] = functionComputeULPowerLevels(Hhat,C,nbrOfRealizations,M,K,L,p,f)
%Compute and categorize the UL signal power for different receive combining
%schemes.
%
%INPUT:
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates 
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%f                 = Pilot reuse factor
%
%OUTPUT:
%powerterms_MR    = 4 x K x L matrix with signal and interference powers
%                   for each of the UEs with MR combining.
%                   powerterms_MR(:,k,l) is the vector for UE k in cell l,
%                   where the first element is the average signal power,
%                   the second element is the average sum interference from
%                   cells that use the same pilot, the third element is the
%                   average sum interference from pilot-sharing UEs, and
%                   the fourth term is the total interference from all UEs
%powerterms_RZF   = Same as powerterms_MR but with RZF combining
%powerterms_MMMSE = Same as powerterms_MR but with M-MMSE combining
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


%Generate pilot pattern
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
    
elseif f == 4 %Only works in the running example with its 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]);
    
elseif f == 16 %Only works in the running example with its 16 BSs
    
    pilotPattern = (1:L)';
    
end


%Compute sum of all estimation error correlation matrices at every BS
C_totM = reshape(p*sum(sum(C(:,:,:,:,:),3),4),[M M L]);

%Compute sum of single-cell estimation error correlation matrix for all cells
C_totMk = zeros(M,M,K,L);

for j = 1:L
    
    %Extract cells that use same pilots as cell j
    groupMembers = find(pilotPattern==pilotPattern(j))';
    
    for k = 1:K
        
        C_totMk(:,:,k,j) = p*sum(C(:,:,k,groupMembers,j),4);
        
    end
    
end


%Prepare to store simulation results
powerterms_MR = zeros(4,K,L);
powerterms_RZF = zeros(4,K,L);
powerterms_MMMSE = zeros(4,K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Extract channel estimate realizations from UEs in cell j to BS j
        Hintraj = reshape(Hhat(:,n,:,j,j),[M K]);
        
        %Extract cells that use the same pilots as cell j
        groupMembers = find(pilotPattern==pilotPattern(j))';
        

        %Compute three different combining schemes
        V_MR = Hallj(:,K*(j-1)+1:K*j); %MR combining in (4.11)
        V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK); %RZF combining in (4.9)
        V_MMMSE = p*(p*(Hallj*Hallj')+C_totM(:,:,j)+eyeM)\V_MR; %M-MMSE combining in (4.7)
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            %Extract channel realizations from UEs that cause pilot
            %contamination to UE k in cell j
            HjkPC = reshape(Hhat(:,n,k,groupMembers,j),[M length(groupMembers)]);
            
            %Extract channel realizations from UEs that cause pilot
            %contamination to any of the UEs in cell j
            HjPC = reshape(Hhat(:,n,:,groupMembers,j),[M K*length(groupMembers)]);
            
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute signal and interference+noise terms using (4.3).
            %The first row is the signal term, the second row is the
            %interference from all UEs in the cells that cause pilot
            %contamination, the third row is the interference only from the
            %UE that causes pilot contamination to UE k in cell j, and the
            %last row is the total interference.
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            membercells = p*sum(abs(v'*HjPC).^2) + v'*sum(C_totMk(:,:,:,j)-C(:,:,:,j,j),3)*v - p*sum(abs(v'*Hintraj).^2);
            pilotcont = p*sum(abs(v'*HjkPC).^2) + v'*(C_totMk(:,:,k,j)-C(:,:,k,j,j))*v - signal;
            totalinterference = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j))*v - signal;
            
            %Average over the realizations of the power terms
            powerterms_MR(:,k,j) = powerterms_MR(:,k,j) + real([signal membercells totalinterference pilotcont])'/nbrOfRealizations;
            
            
            
            %%RZF combining
            v = V_RZF(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute signal and interference+noise terms using (4.3),
            %see the details provided above
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            membercells = p*sum(abs(v'*HjPC).^2) + v'*sum(C_totMk(:,:,:,j)-C(:,:,:,j,j),3)*v - p*sum(abs(v'*Hintraj).^2);
            pilotcont = p*sum(abs(v'*HjkPC).^2) + v'*(C_totMk(:,:,k,j)-C(:,:,k,j,j))*v - signal;
            totalinterference = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j))*v - signal;
            
            %Average over the realizations of the power terms
            powerterms_RZF(:,k,j) = powerterms_RZF(:,k,j) + real([signal membercells totalinterference pilotcont])'/nbrOfRealizations;
            
            
            %%M-MMSE combining
            v = V_MMMSE(:,k); %Extract combining vector
            v = v/norm(v);
            
            %Compute signal and interference+noise terms using (4.3),
            %see the details provided above
            signal = p*abs(v'*Hhat(:,n,k,j,j))^2;
            membercells = p*sum(abs(v'*HjPC).^2) + v'*sum(C_totMk(:,:,:,j)-C(:,:,:,j,j),3)*v - p*sum(abs(v'*Hintraj).^2);
            pilotcont = p*sum(abs(v'*HjkPC).^2) + v'*(C_totMk(:,:,k,j)-C(:,:,k,j,j))*v - signal;
            totalinterference = p*sum(abs(v'*Hallj).^2) + v'*(C_totM(:,:,j))*v - signal;
            
            %Average over the realizations of the power terms
            powerterms_MMMSE(:,k,j) = powerterms_MMMSE(:,k,j) + real([signal membercells totalinterference pilotcont])'/nbrOfRealizations;
            
        end
        
    end
    
end
