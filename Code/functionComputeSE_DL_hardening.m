function [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_DL_hardening(H,Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho)
%Compute DL SE for different transmit precoding schemes using Theorem 4.6.
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
prelogFactor = (tau_c-tau_p)./(tau_c);

%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L);

if nargout > 1
    signal_RZF = zeros(K,L);
end

if nargout > 2
    signal_MMMSE = zeros(K,L);
end

if nargout > 3
    signal_ZF = zeros(K,L);
end

if nargout > 4
    signal_SMMSE = zeros(K,L);
end


%Prepare to store simulation results for sum interference powers
interf_MR = zeros(K,L);

if nargout > 1
    interf_RZF = zeros(K,L);
end

if nargout > 2
    interf_MMMSE = zeros(K,L);
end

if nargout > 3
    interf_ZF = zeros(K,L);
end

if nargout > 4
    interf_SMMSE = zeros(K,L);
end



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
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            if norm(V_MR(:,k))>0
                
                %%MR precoding in (4.37) computed based on combining vector
                w = V_MR(:,k)/norm(V_MR(:,k));
                
                %Compute realizations of the terms inside the expectations
                %of the signal and interference terms of (4.26)
                signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                interf_MR = interf_MR + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                
                
                %%RZF precoding
                if nargout > 1
                    
                    %RZF precoding in (4.37) computed based on combining vector
                    w = V_RZF(:,k)/norm(V_RZF(:,k));
                    
                    %Compute realizations of the terms inside the expectations
                    %of the signal and interference terms of (4.26)
                    signal_RZF(k,j) = signal_RZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_RZF = interf_RZF + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
                
                %%M-MMSE precoding
                if nargout > 2
                    
                    %M-MMSE precoding in (4.37) computed based on combining vector
                    w = V_MMMSE(:,k)/norm(V_MMMSE(:,k));
                    
                    %Compute realizations of the terms inside the expectations
                    %of the signal and interference terms of (4.26)
                    signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_MMMSE = interf_MMMSE + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
                
                %%ZF precoding
                if nargout > 3
                    
                    %ZF precoding in (4.37) computed based on combining vector
                    w = V_ZF(:,k)/norm(V_ZF(:,k));
                    
                    %Compute realizations of the terms inside the expectations
                    %of the signal and interference terms of (4.26)
                    signal_ZF(k,j) = signal_ZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_ZF = interf_ZF + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
                
                %%S-MMSE precoding
                if nargout > 4
                    
                    %S-MMSE precoding in (4.37) computed based on combining vector
                    w = V_SMMSE(:,k)/norm(V_SMMSE(:,k));
                    
                    %Compute realizations of the terms inside the expectations
                    %of the signal and interference terms of (4.26)
                    signal_SMMSE(k,j) = signal_SMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                    interf_SMMSE = interf_SMMSE + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                    
                end
                
            end
            
        end
        
    end
    
end


%% Prepare to compute the SEs with different pilot reuse factors
SE_MR = zeros(K,L,length(tau_c));

if nargout > 1
    SE_RZF = zeros(K,L,length(tau_c));
end

if nargout > 2
    SE_MMMSE = zeros(K,L,length(tau_c));
end

if nargout > 3
    SE_ZF = zeros(K,L,length(tau_c));
end

if nargout > 4
    SE_SMMSE = zeros(K,L,length(tau_c));
end

%Go through all lengths of the pilot reuse factors
for n = 1:length(tau_c)
    
    %Compute SEs according to Theorem 4.6
    SE_MR(:,:,n) = prelogFactor(n)*real(log2(1+(rho*abs(signal_MR).^2) ./ (interf_MR - rho*abs(signal_MR).^2 + 1)));
    
    if nargout > 1
        SE_RZF(:,:,n) = prelogFactor(n)*real(log2(1+(rho*abs(signal_RZF).^2) ./ (interf_RZF - rho*abs(signal_RZF).^2 + 1)));
    end
    
    if nargout > 2
        SE_MMMSE(:,:,n) = prelogFactor(n)*real(log2(1+(rho*abs(signal_MMMSE).^2) ./ (interf_MMMSE - rho*abs(signal_MMMSE).^2 +1)));
    end
    
    if nargout > 3
        SE_ZF(:,:,n) = prelogFactor(n)*real(log2(1+(rho*abs(signal_ZF).^2) ./ (interf_ZF - rho*abs(signal_ZF).^2 + 1)));
    end
    
    if nargout > 4
        SE_SMMSE(:,:,n) = prelogFactor(n)*real(log2(1+(rho*abs(signal_SMMSE).^2) ./ (interf_SMMSE - rho*abs(signal_SMMSE).^2 + 1)));
    end
    
end
