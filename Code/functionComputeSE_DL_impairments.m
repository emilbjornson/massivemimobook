function [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_impairments(H,Hhat,C,tau_c,tau_p,nbrOfRealizations,M,K,L,p,rho,kappatBS,kapparUE)
%Compute DL SE for different transmit precoding schemes with hardware
%impairments, based on Theorem 6.5.
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
%rho               = Downlink transmit power per UE (same for everyone)
%kappatUE          = Hardware quality of the UEs' transmitters
%kapparBS          = Hardware quality of the BSs' receivers
%
%OUTPUT:
%SE_MR    = K x L matrix where element (k,l) is the downlink SE of UE k in
%           cell l achieved with MR precoding
%SE_RZF   = Same as SE_MR but with RZF precoding
%SE_MMMSE = Same as SE_MR but with M-MMSE precoding
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

%Compute the pre-log factor assuming only downlink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store simulation results for signal gains
signal_MR = zeros(K,L);
signal_RZF = zeros(K,L);
signal_MMMSE = zeros(K,L);

%Prepare to store simulation results for norms of precoding vectors
precodingNorm_MR = zeros(K,L);
precodingNorm_RZF = zeros(K,L);
precodingNorm_MMMSE = zeros(K,L);

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
            
            %%MR precoding
            w = V_MR(:,k)/norm(V_MR(:,k)); %Extract precoding vector
            wrep = repmat(w,[1 K*L]);
            
            %Compute realizations of the expectations in signal and
            %interference terms in (6.46)
            signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
            precodingNorm_MR(k,j) = precodingNorm_MR(k,j) + norm(w).^2/nbrOfRealizations;
            interf_MR = interf_MR + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
            impair_MR = impair_MR + rho*reshape(sum(abs(wrep.*Hallj).^2,1),[K L])/nbrOfRealizations;
            
            
            %%RZF precoding
            if nargout >=2
                
                w = V_RZF(:,k)/norm(V_RZF(:,k)); %Extract precoding vector
                wrep = repmat(w,[1 K*L]);
                
                %Compute realizations of the expectations in signal and
                %interference terms in (6.46)
                signal_RZF(k,j) = signal_RZF(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                precodingNorm_RZF(k,j) = precodingNorm_RZF(k,j) + norm(w).^2/nbrOfRealizations;
                interf_RZF = interf_RZF + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                impair_RZF = impair_RZF + rho*reshape(sum(abs(wrep.*Hallj).^2,1),[K L])/nbrOfRealizations;
                
            end
            
            %%M-MMSE precoding
            if nargout >=3
                
                w = V_MMMSE(:,k)/norm(V_MMMSE(:,k)); %Extract precoding vector
                wrep = repmat(w,[1 K*L]);
                
                %Compute realizations of the expectations in signal and
                %interference terms in (6.46)
                signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/nbrOfRealizations;
                precodingNorm_MMMSE(k,j) = precodingNorm_MMMSE(k,j) + norm(w).^2/nbrOfRealizations;
                interf_MMMSE = interf_MMMSE + rho*reshape(abs(w'*Hallj).^2,[K L])/nbrOfRealizations;
                impair_MMMSE = impair_MMMSE + rho*reshape(sum(abs(wrep.*Hallj).^2,1),[K L])/nbrOfRealizations;
                
            end
            
        end
        
    end
    
end

%Precompute terms that appear multiple times
sigma2 = 1/(kapparUE*kappatBS);
factor1 = 1/kapparUE;
factor2 = (1-kappatBS)/(kapparUE*kappatBS);

%Compute SEs according to Theorem 6.5 and (6.46)
SE_MR = prelogFactor*real(log2(1+(rho*abs(signal_MR).^2) ./ (factor1*interf_MR + factor2*impair_MR - rho*abs(signal_MR).^2 + sigma2)));

if nargout >=2
    SE_RZF = prelogFactor*real(log2(1+(rho*abs(signal_RZF).^2) ./ (factor1*interf_RZF + factor2*impair_RZF - rho*abs(signal_RZF).^2 + sigma2)));
end

if nargout >=3
    SE_MMMSE = prelogFactor*real(log2(1+(rho*abs(signal_MMMSE).^2) ./ (factor1*interf_MMMSE + factor2*impair_MMMSE - rho*abs(signal_MMMSE).^2 +sigma2)));
end
