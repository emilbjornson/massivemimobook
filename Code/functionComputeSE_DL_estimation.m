function [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_DL_estimation(H,Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho)
%Compute DL SE for different transmit precoding schemes using Theorem 4.9.
%
%INPUT:
%H                 = M x nbrOfRealizations x K x L x L matrix with the
%                    exact channel realizations
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%rho               = Downlink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MR    = K x L x length(tau_c) matrix where element (k,l,:) is the
%           downlink SE of UE k in cell l achieved with MR precoding for
%           different lengths of the coherence block
%SE_RZF   = Same as SE_MR but with RZF precoding
%SE_MMMSE = Same as SE_MR but with M-MMSE precoding
%SE_ZF    = Same as SE_MR but with ZF precoding
%SE_SMMSE = Same as SE_MR but with S-MMSE precoding
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
%This is version 1.01 (Last edited: 2018-11-21)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Store identity matrices of different sizes
eyeK = eye(K);
eyeM = eye(M);

if nargout > 2
    
    %Compute sum of all estimation error correlation matrices at every BS
    C_totM = reshape(p*sum(sum(C,3),4),[M M L]);
    
end

%Compute sum of intra-cell estimation error correlation matrices at every BS
if nargout>4
    
    CR_totS = zeros(M,M,L);
    
    for j = 1:L
        CR_totS(:,:,j) = p*(sum(C(:,:,:,j,j),3)+sum(sum(R(:,:,:,[1:j-1 j+1:end],j),3),4));
    end
    
end


%Compute the prelog factor assuming only downlink transmission
prelogFactor = (tau_c-tau_p)./tau_c;


%Prepare to store simulation results for intra-cell channel gains
gain_MR = zeros(K,L,K);

%Prepare to store simulation results for squared intra-cell channel gains
gain2_MR = zeros(K,L,K);

%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L,nbrOfRealizations);

%Prepare to store simulation results for intra-cell interference powers
intraInterf_MR = zeros(K,L,nbrOfRealizations);

%Prepare to store simulation results for inter-cell interference powers
interInterf_MR = zeros(K,L);


if nargout > 1
    
    %Prepare to store simulation results for intra-cell channel gains
    gain_RZF = zeros(K,L,K);
    
    %Prepare to store simulation results for squared intra-cell channel gains
    gain2_RZF = zeros(K,L,K);
    
    %Prepare to store simulation results for signal gains
    signal_RZF = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for intra-cell interference powers
    intraInterf_RZF = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for inter-cell interference powers
    interInterf_RZF = zeros(K,L);
    
end


if nargout > 2
    
    %Prepare to store simulation results for intra-cell channel gains
    gain_MMMSE = zeros(K,L,K);
    
    %Prepare to store simulation results for squared intra-cell channel gains
    gain2_MMMSE = zeros(K,L,K);
    
    %Prepare to store simulation results for signal gains
    signal_MMMSE = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for intra-cell interference powers
    intraInterf_MMMSE = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for inter-cell interference powers
    interInterf_MMMSE = zeros(K,L);
    
end


if nargout > 3
    
    %Prepare to store simulation results for intra-cell channel gains
    gain_ZF = zeros(K,L,K);
    
    %Prepare to store simulation results for squared intra-cell channel gains
    gain2_ZF = zeros(K,L,K);
    
    %Prepare to store simulation results for signal gains
    signal_ZF = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for intra-cell interference powers
    intraInterf_ZF = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for inter-cell interference powers
    interInterf_ZF = zeros(K,L);
    
end


if nargout > 4
    
    %Prepare to store simulation results for intra-cell channel gains
    gain_SMMSE = zeros(K,L,K);
    
    %Prepare to store simulation results for squared intra-cell channel gains
    gain2_SMMSE = zeros(K,L,K);
    
    %Prepare to store simulation results for signal gains
    signal_SMMSE = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for intra-cell interference powers
    intraInterf_SMMSE = zeros(K,L,nbrOfRealizations);
    
    %Prepare to store simulation results for inter-cell interference powers
    interInterf_SMMSE = zeros(K,L);
    
end



%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j); %MR in Eq. (XX)
        
        if nargout > 1 %Compute RZF combining in (4.9)
            V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        end
        
        if nargout > 2 %Compute M-MMSE combining in (4.7)
            V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        end
        
        if nargout > 3 %Compute ZF combining in (4.10), with the small regularization term 1e-12 for numerical stability
            V_ZF = V_MR/(V_MR'*V_MR+1e-12*eyeK);
        end
        
        if nargout > 4 %Compute S-MMSE combining in (4.8)
            V_SMMSE = p*(p*(V_MR*V_MR')+CR_totS(:,:,j)+eyeM)\V_MR;
        end
        
        
        %Go through all users in cell j
        for k = 1:K
            
            %Consider non-zero channels
            if norm(V_MR(:,k))>0
                
                %%MR precoding in (4.37) computed based on combining vector
                w = V_MR(:,k)/norm(V_MR(:,k)); %Extract precoding vector
                
                %Compute signal term in (4.39)
                signal_MR(k,j,n) = abs(w'*H(:,n,k,j,j))^2;
                
                %Compute inner product with all channels
                innerproductsAll = reshape(w'*Hallj,[K L]);
                
                %Compute inner product with intra-cell channels
                interferenceIntra = Hallj(:,K*(j-1)+1:K*j)'*w;
                
                %Add interference caused to all UEs, which will be averaged
                %with respect to small-scale fading realizations
                interInterf_MR = interInterf_MR + abs(innerproductsAll).^2/nbrOfRealizations;
                
                %Remove the interference caused to intra-cell UEs
                interInterf_MR(:,j) = interInterf_MR(:,j) - abs(interferenceIntra).^2/nbrOfRealizations;
                
                %Compute the interference caused to all intra-cell UEs
                intraInterf_MR(:,j,n) = intraInterf_MR(:,j,n) + abs(interferenceIntra).^2;
                
                %Remove the interference caused to the UE itself
                intraInterf_MR(k,j,n) = intraInterf_MR(k,j,n) - signal_MR(k,j,n);
                
                %Compute mean and quadratic mean of the precoded channel by
                %Monte Carlo simulations
                gain_MR(:,j,k) = gain_MR(:,j,k) + interferenceIntra/nbrOfRealizations;
                gain2_MR(:,j,k) = gain2_MR(:,j,k) + abs(interferenceIntra).^2/nbrOfRealizations;
                
                
                %%RZF precoding
                if nargout > 1
                    
                    %RZF precoding in (4.37) computed based on combining vector
                    w = V_RZF(:,k)/norm(V_RZF(:,k));
                    
                    %Compute signal and interference+noise terms using Eq. (XX)
                    signal_RZF(k,j,n) = abs(w'*H(:,n,k,j,j))^2;
                    
                    %Compute inner product with all channels
                    innerproductsAll = reshape(w'*Hallj,[K L]);
                    
                    %Compute inner product with intra-cell channels
                    interferenceIntra = Hallj(:,K*(j-1)+1:K*j)'*w;
                    
                    %Add interference caused to all UEs, which will be averaged
                    %with respect to small-scale fading realizations
                    interInterf_RZF = interInterf_RZF + abs(innerproductsAll).^2/nbrOfRealizations;
                    
                    %Remove the interference caused to intra-cell UEs
                    interInterf_RZF(:,j) = interInterf_RZF(:,j) - abs(interferenceIntra).^2/nbrOfRealizations;
                    
                    %Compute the interference caused to all intra-cell UEs
                    intraInterf_RZF(:,j,n) = intraInterf_RZF(:,j,n) + abs(interferenceIntra).^2;
                    
                    %Remove the interference caused to the UE itself
                    intraInterf_RZF(k,j,n) = intraInterf_RZF(k,j,n) - signal_RZF(k,j,n);
                    
                    %Compute mean and quadratic mean of the precoded channel by
                    %Monte Carlo simulations
                    gain_RZF(:,j,k) = gain_RZF(:,j,k) + interferenceIntra/nbrOfRealizations;
                    gain2_RZF(:,j,k) = gain2_RZF(:,j,k) + abs(interferenceIntra).^2/nbrOfRealizations;
                    
                end
                
                
                %%M-MMSE precoding
                if nargout > 2
                    
                    %M-MMSE precoding in (4.37) computed based on combining vector
                    w = V_MMMSE(:,k)/norm(V_MMMSE(:,k));
                    
                    %Compute signal and interference+noise terms using Eq. (XX)
                    signal_MMMSE(k,j,n) = abs(w'*H(:,n,k,j,j))^2;
                    
                    %Compute inner product with all channels
                    innerproductsAll = reshape(w'*Hallj,[K L]);
                    
                    %Compute inner product with intra-cell channels
                    interferenceIntra = Hallj(:,K*(j-1)+1:K*j)'*w;
                    
                    %Add interference caused to all UEs, which will be averaged
                    %with respect to small-scale fading realizations
                    interInterf_MMMSE = interInterf_MMMSE + abs(innerproductsAll).^2/nbrOfRealizations;
                    
                    %Remove the interference caused to intra-cell UEs
                    interInterf_MMMSE(:,j) = interInterf_MMMSE(:,j) - abs(interferenceIntra).^2/nbrOfRealizations;
                    
                    %Compute the interference caused to all intra-cell UEs
                    intraInterf_MMMSE(:,j,n) = intraInterf_MMMSE(:,j,n) + abs(interferenceIntra).^2;
                    
                    %Remove the interference caused to the UE itself
                    intraInterf_MMMSE(k,j,n) = intraInterf_MMMSE(k,j,n) - signal_MMMSE(k,j,n);
                    
                    %Compute mean and quadratic mean of the precoded channel by
                    %Monte Carlo simulations
                    gain_MMMSE(:,j,k) = gain_MMMSE(:,j,k) + interferenceIntra/nbrOfRealizations;
                    gain2_MMMSE(:,j,k) = gain2_MMMSE(:,j,k) + abs(interferenceIntra).^2/nbrOfRealizations;
                    
                end
                
                
                %%ZF precoding
                if nargout > 3
                    
                    %ZF precoding in (4.37) computed based on combining vector
                    w = V_ZF(:,k)/norm(V_ZF(:,k));
                    
                    %Compute signal and interference+noise terms using Eq. (XX)
                    signal_ZF(k,j,n) = abs(w'*H(:,n,k,j,j))^2;
                    
                    %Compute inner product with all channels
                    innerproductsAll = reshape(w'*Hallj,[K L]);
                    
                    %Compute inner product with intra-cell channels
                    interferenceIntra = Hallj(:,K*(j-1)+1:K*j)'*w;
                    
                    %Add interference caused to all UEs, which will be averaged
                    %with respect to small-scale fading realizations
                    interInterf_ZF = interInterf_ZF + abs(innerproductsAll).^2/nbrOfRealizations;
                    
                    %Remove the interference caused to intra-cell UEs
                    interInterf_ZF(:,j) = interInterf_ZF(:,j) - abs(interferenceIntra).^2/nbrOfRealizations;
                    
                    %Compute the interference caused to all intra-cell UEs
                    intraInterf_ZF(:,j,n) = intraInterf_ZF(:,j,n) + abs(interferenceIntra).^2;
                    
                    %Remove the interference caused to the UE itself
                    intraInterf_ZF(k,j,n) = intraInterf_ZF(k,j,n) - signal_ZF(k,j,n);
                    
                    %Compute mean and quadratic mean of the precoded channel by
                    %Monte Carlo simulations
                    gain_ZF(:,j,k) = gain_ZF(:,j,k) + interferenceIntra/nbrOfRealizations;
                    gain2_ZF(:,j,k) = gain2_ZF(:,j,k) + abs(interferenceIntra).^2/nbrOfRealizations;
                    
                end
                
                
                %%S-MMSE precoding
                if nargout > 4
                    
                    %S-MMSE precoding in (4.37) computed based on combining vector
                    w = V_SMMSE(:,k)/norm(V_SMMSE(:,k));
                    
                    %Compute signal and interference+noise terms using Eq. (XX)
                    signal_SMMSE(k,j,n) = abs(w'*H(:,n,k,j,j))^2;
                    
                    %Compute inner product with all channels
                    innerproductsAll = reshape(w'*Hallj,[K L]);
                    
                    %Compute inner product with intra-cell channels
                    interferenceIntra = Hallj(:,K*(j-1)+1:K*j)'*w;
                    
                    %Add interference caused to all UEs, which will be averaged
                    %with respect to small-scale fading realizations
                    interInterf_SMMSE = interInterf_SMMSE + abs(innerproductsAll).^2/nbrOfRealizations;
                    
                    %Remove the interference caused to intra-cell UEs
                    interInterf_SMMSE(:,j) = interInterf_SMMSE(:,j) - abs(interferenceIntra).^2/nbrOfRealizations;
                    
                    %Compute the interference caused to all intra-cell UEs
                    intraInterf_SMMSE(:,j,n) = intraInterf_SMMSE(:,j,n) + abs(interferenceIntra).^2;
                    
                    %Remove the interference caused to the UE itself
                    intraInterf_SMMSE(k,j,n) = intraInterf_SMMSE(k,j,n) - signal_SMMSE(k,j,n);
                    
                    %Compute mean and quadratic mean of the precoded channel by
                    %Monte Carlo simulations
                    gain_SMMSE(:,j,k) = gain_SMMSE(:,j,k) + interferenceIntra/nbrOfRealizations;
                    gain2_SMMSE(:,j,k) = gain2_SMMSE(:,j,k) + abs(interferenceIntra).^2/nbrOfRealizations;
                    
                end
                
            end
            
        end
        
    end
    
end


%% Prepare to compute the SEs with different pilot reuse factors
SE_MR = zeros(K,L,length(tau_c));
SE_ZF = zeros(K,L,length(tau_c));
SE_SMMSE = zeros(K,L,length(tau_c));
SE_RZF = zeros(K,L,length(tau_c));
SE_MMMSE = zeros(K,L,length(tau_c));


%Compute the first term in (4.38), except for prelog factor, for different
%precoding schemes. These terms represent the SE with perfect intra-cell CSI
SE_MR_perfect = mean(log2(1+(rho*signal_MR) ./ (rho*intraInterf_MR + repmat(rho*interInterf_MR,[1 1 nbrOfRealizations]) + 1)),3);

if nargout > 1
    SE_RZF_perfect = mean(log2(1+(rho*signal_RZF) ./ (rho*intraInterf_RZF + repmat(rho*interInterf_RZF,[1 1 nbrOfRealizations]) + 1)),3);
end

if nargout > 2
    SE_MMMSE_perfect = mean(log2(1+(rho*signal_MMMSE) ./ (rho*intraInterf_MMMSE + repmat(rho*interInterf_MMMSE,[1 1 nbrOfRealizations]) + 1)),3);
end

if nargout > 3
    SE_ZF_perfect = mean(log2(1+(rho*signal_ZF) ./ (rho*intraInterf_ZF + repmat(rho*interInterf_ZF,[1 1 nbrOfRealizations]) + 1)),3);
end

if nargout > 4
    SE_SMMSE_perfect = mean(log2(1+(rho*signal_SMMSE) ./ (rho*intraInterf_SMMSE + repmat(rho*interInterf_SMMSE,[1 1 nbrOfRealizations]) + 1)),3);
end


%Go through all lengths of the pilot reuse factors
for n = 1:length(tau_c)
    
    SE_MR(:,:,n) = prelogFactor(n)*SE_MR_perfect - log2(prod(1+rho*(tau_c(n)-tau_p)*(gain2_MR-abs(gain_MR).^2),3))/tau_c(n);
    
    if nargout > 1
        SE_RZF(:,:,n) = prelogFactor(n)*SE_RZF_perfect - log2(prod(prod(1+rho*(tau_c(n)-tau_p)*(gain2_RZF-abs(gain_RZF).^2),3),4))/tau_c(n);
    end
    
    if nargout > 2
        SE_MMMSE(:,:,n) = prelogFactor(n)*SE_MMMSE_perfect - log2(prod(prod(1+rho*(tau_c(n)-tau_p)*(gain2_MMMSE-abs(gain_MMMSE).^2),3),4))/tau_c(n);
    end
    
    if nargout > 3
        SE_ZF(:,:,n) = prelogFactor(n)*SE_ZF_perfect - log2(prod(prod(1+rho*(tau_c(n)-tau_p)*(gain2_ZF-abs(gain_ZF).^2),3),4))/tau_c(n);
    end
    
    if nargout > 4
        SE_SMMSE(:,:,n) = prelogFactor(n)*SE_SMMSE_perfect - log2(prod(prod(1+rho*(tau_c(n)-tau_p)*(gain2_SMMSE-abs(gain_SMMSE).^2),3),4))/tau_c(n);
    end
    
end

%Turn negative SE values into zero, since the SE cannot be negative
SE_MR(SE_MR<0) = 0;

if nargout > 1
    SE_RZF(SE_RZF<0) = 0;
end

if nargout > 2
    SE_MMMSE(SE_MMMSE<0) = 0;
end

if nargout > 3
    SE_ZF(SE_ZF<0) = 0;
end

if nargout > 4
    SE_SMMSE(SE_SMMSE<0) = 0;
end
