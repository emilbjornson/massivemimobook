function [Hhat_MMSE,C_MMSE,tau_p,R,H,Hhat_EW_MMSE,C_EW_MMSE,Hhat_LS,C_LS] = functionChannelEstimates(R,channelGaindB,nbrOfRealizations,M,K,L,p,f)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are assumed to be correlated
%Rayleigh fading. The MMSE estimator, EW-MMSE estimator, and LS estimator
%are used. The latter two estimators are only computed if their outputs are
%requested when calling the function.
%
%INPUT:
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices for all UEs in the network. R(:,:,k,j,l) is
%                    the correlation matrix for the channel between UE k 
%                    in cell j and the BS in cell l. This such matrix can
%                    either include the average channel gain or can be
%                    normalized arbitrarily.
%channelGaindB     = K x L x L matrix containing the average channel gains
%                    in dB of all the channels, if these are not already
%                    included in the spatial channel correlation matrices. 
%                    The product R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10) 
%                    is the full spatial channel correlation matrix.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%f                 = Pilot reuse factor
%
%OUTPUT:
%Hhat_MMSE    = M x nbrOfRealizations x K x L x L matrix with the MMSE
%               channel estimates. The matrix Hhat_MMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_MMSE       = M x M x K x L x L matrix with estimation error correlation
%               matrices when using MMSE estimation. The matrix is
%               organized in the same way as R.
%tau_p        = Length of pilot sequences
%R            = Scaled version of the input spatial correlation matrices R,
%               where the channel gains from channelGaindB are included
%H            = M x nbrOfRealizations x K x L x L matrix with the true
%               channel realizations. The matrix is organized in the same
%               way as Hhat_MMSE.
%Hhat_EW_MMSE = Same as Hhat_MMSE, but using the EW-MMSE estimator
%C_EW_MMSE    = Same as C_MMSE, but using the EW-MMSE estimator
%Hhat_LS      = Same as Hhat_MMSE, but using the LS estimator
%C_LS         = Same as C_MMSE, but using the LS estimator
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


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

%Prepare a matrix to save the channel gains per UE
betas = zeros(K,L,L);


%Go through all channels and apply the channel gains to the spatial 
%correlation matrices
for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            
            if channelGaindB(k,j,l)>-Inf
                
                %Extract channel gain in linear scale
                betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
                
                %Apply channel gain to correlation matrix
                R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
                
                %Apply correlation to the uncorrelated channel realizations
                Rsqrt = sqrtm(R(:,:,k,j,l));
                H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
                
            else
                
                betas(k,j,l) = 0;
                R(:,:,k,j,l) = 0;
                H(:,:,k,j,l) = 0;
                
            end
            
        end
        
    end
    
end



%% Perform channel estimation


%Length of pilot sequences
tau_p = f*K;


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


%Store identity matrix of size M x M
eyeM = eye(M);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));



%Prepare for MMSE estimation

%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error correlation matrices
C_MMSE = zeros(M,M,K,L,L);


%Prepare for EW-MMSE estimation
if nargout >= 5
    
    %Prepare to store EW-MMSE channel estimates
    Hhat_EW_MMSE = zeros(M,nbrOfRealizations,K,L,L);
    
    %Prepare to store estimation error correlation matrices
    C_EW_MMSE = zeros(M,M,K,L,L);
    
end

%Prepare for LS estimation
if nargout >= 7
    
    %Prepare to store EW-MMSE channel estimates
    Hhat_LS = zeros(M,nbrOfRealizations,K,L,L);
    
    %Prepare to store estimation error correlation matrices
    C_LS = zeros(M,M,K,L,L);
    
end


%% Go through all cells
for j = 1:L
    
    %Go through all f pilot groups
    for g = 1:f
        
        %Extract the cells that belong to pilot group g
        groupMembers = find(g==pilotPattern)';
        
        %Compute processed pilot signal for all UEs that use these pilots, according to (3.5)
        yp = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
        
        %Go through all UEs
        for k = 1:K
            
            %Compute the matrix that is inverted in the MMSE estimator
            PsiInv = (p*tau_p*sum(R(:,:,k,g==pilotPattern,j),4) + eyeM);
            
            %If EW-MMSE estimation is to be computed
            if nargout >= 5
                %Compute a vector with elements that are inverted in the EW-MMSE estimator
                PsiInvDiag = diag(PsiInv);
            end
            
            %Go through the cells in pilot group g
            for l = groupMembers
                
                %Compute MMSE estimate of channel between BS l and UE k in
                %cell j using (3.9) in Theorem 3.1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = sqrt(p)*RPsi*yp(:,:,k);
                
                %Compute corresponding estimation error correlation matrix, using (3.11)
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);
                
                
                %If EW-MMSE estimation is to be computed
                if nargout >= 5
                    
                    %Compute EW-MMSE estimate of channel between BS l and
                    %UE k in cell j using (3.33)
                    A_EW_MMSE = diag(sqrt(p)*diag(R(:,:,k,l,j)) ./ PsiInvDiag);
                    Hhat_EW_MMSE(:,:,k,l,j) = A_EW_MMSE*yp(:,:,k);
                    
                    %Compute corresponding estimation error correlation
                    %matrix, using the principle from (3.29)
                    productAR = A_EW_MMSE * R(:,:,k,l,j);
                    
                    C_EW_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_EW_MMSE*PsiInv*A_EW_MMSE';
                
                end
                
                
                %If LS estimation is to be computed
                if nargout >= 7
                    
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
    
end
