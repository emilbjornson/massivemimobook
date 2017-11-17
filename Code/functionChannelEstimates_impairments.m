function [Hhat_LMMSE,C_LMMSE,tau_p,R,H] = functionChannelEstimates_impairments(R,channelGaindB,nbrOfRealizations,M,K,L,p,f,kappatUE,kapparBS)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are assumed to be correlated
%Rayleigh fading. The LMMSE estimator that takes transceiver hardware
%impairments into account is used.
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
%kappatUE          = Hardware quality of the UEs' transmitters
%kapparBS          = Hardware quality of the BSs' receivers
%
%OUTPUT:
%Hhat_LMMSE   = M x nbrOfRealizations x K x L x L matrix with the LMMSE
%               channel estimates. The matrix Hhat_LMMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_LMMSE      = M x M x K x L x L matrix with estimation error correlation
%               matrices when using LMMSE estimation. The matrix is
%               organized in the same way as R.
%tau_p        = Length of pilot sequences
%R            = Scaled version of the input spatial correlation matrices R,
%               where the channel gains from channelGaindB are included
%H            = M x nbrOfRealizations x K x L x L matrix with the true
%               channel realizations. The matrix is organized in the same
%               way as Hhat_LMMSE.
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
            
            %Extract channel gain in linear scale
            betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
            
            %Apply channel gain to correlation matrix
            R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
            
            %Apply correlation to the uncorrelated channel realizations
            Rsqrt = sqrtm(R(:,:,k,j,l));
            H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
            
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


%Prepare to store LMMSE channel estimates
Hhat_LMMSE = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error correlation matrices
C_LMMSE = zeros(M,M,K,L,L);


%% Go through all f pilot groups
for g = 1:f
    
    %Create transmit distortion term
    etaUE = sqrt(0.5*p*tau_p*kapparBS*(1-kappatUE))*(randn(1,nbrOfRealizations,K,L,K) + 1i*randn(1,nbrOfRealizations,K,L,K));
    
    
    %Go through all cells
    for j = 1:L
        
        %Extract the cells that belong to pilot group g
        groupMembers = find(g==pilotPattern)';
        
        
        %Go through all UEs
        for k = 1:K
            
            %Create receive distortion term
            etabarBS = sqrt(0.5*p*tau_p*(1-kapparBS))*(randn(M,nbrOfRealizations,K,L) + 1i*randn(M,nbrOfRealizations,K,L));
            
            %Sum up transmit and receive distortion terms
            etaSum = repmat(etaUE(:,:,:,:,k),[M 1 1 1]) + etabarBS;
            
            %Compute processed pilot signal for all UEs that use these pilots, according to (6.22)
            ypk = sqrt(p*kappatUE*kapparBS)*tau_p*sum(H(:,:,k,g==pilotPattern,j),4) + squeeze(sum(sum(H(:,:,:,:,j).*etaSum,4),3)) + sqrt(tau_p)*Np(:,:,k,j,g);
            
            %Compute sum of all UEs' correlation matrices to BS j
            RjAll = sum(sum(R(:,:,:,:,j),3),4);
            
            %Compute the matrix that is inverted in the LMMSE estimator
            PsiInv = (p*kappatUE*kapparBS*tau_p*sum(R(:,:,k,g==pilotPattern,j),4) + p*(1-kappatUE)*kapparBS*RjAll + p*(1-kapparBS)*diag(diag(RjAll)) + eyeM);
            
            %Go through the cells in pilot group g
            for l = groupMembers
                
                %Compute LMMSE estimate of channel between BS l and UE k in
                %cell j using (6.23) in Theorem 6.1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_LMMSE(:,:,k,l,j) = sqrt(p*kappatUE*kapparBS)*RPsi*ypk;
                
                %Compute corresponding estimation error correlation matrix, using (6.26)
                C_LMMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*kappatUE*kapparBS*tau_p*RPsi*R(:,:,k,l,j);
                
            end
            
        end
        
    end
    
end
