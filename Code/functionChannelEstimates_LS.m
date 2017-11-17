function [Hhat_LS,C_LS,tau_p] = functionChannelEstimates_LS(H,R,nbrOfRealizations,M,K,L,p)
%Generate LS channel estimates of all UEs in the entire network, given a
%set of channel realizations of arbitary distribution.
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    channel realizations.
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices for all UEs in the network. R(:,:,k,j,l) is
%                    the correlation matrix for the channel between UE k
%                    in cell j and the BS in cell l.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat_LS      = M x nbrOfRealizations x K x L x L matrix with the LS
%               channel estimates. The matrix Hhat_MMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_LS         = M x M x K x L x L matrix with estimation error correlation
%               matrices. The matrix is organized in the same way as R.
%tau_p        = Length of pilot sequences
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


%Length of pilot sequences
tau_p = K;

%Store identity matrix of size M x M
eyeM = eye(M);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L) + 1i*randn(M,nbrOfRealizations,K,L));

%Prepare to store LS channel estimates
Hhat_LS = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error correlation matrices
C_LS = zeros(M,M,K,L,L);


%% Go through all cells
for j = 1:L
    
    %Compute processed pilot signal for all UEs that use these pilots, according to (3.5)
    yp = sqrt(p)*tau_p*sum(H(:,:,:,:,j),4) + sqrt(tau_p)*Np(:,:,:,j);
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the matrix that is inverted in the MMSE estimator
        PsiInv = (p*tau_p*sum(R(:,:,k,:,j),4) + eyeM);
        
        %Go through all cells
        for l = 1:L
            
            %Check if the UE is active (inactive UEs have zero matrices)
            if trace(R(:,:,k,l,j))>0
                
                %Compute LS estimate of channel between BS l and UE k
                %in cell j using (3.35) and (3.36)
                A_LS = 1/(sqrt(p)*tau_p);
                Hhat_LS(:,:,k,l,j) = A_LS*yp(:,:,k);
                
                %Compute corresponding estimation error correlation
                %matrix, using the principle from (3.29)
                productAR = A_LS * R(:,:,k,l,j);
                
                C_LS(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_LS*PsiInv*A_LS';
                
            end
            
        end
        
    end
    
end
