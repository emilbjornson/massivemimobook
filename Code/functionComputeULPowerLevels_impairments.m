function [powerterms_MR,powerterms_RZF,powerterms_MMMSE] = functionComputeULPowerLevels_impairments(H,Hhat,C,nbrOfRealizations,M,K,L,p,f,kappatUE,kapparBS)
%Compute and categorize the UL signal power for different receive combining
%schemes, under hardware impairments.
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
%f                 = Pilot reuse factor
%kappatUE          = Hardware quality of the UEs' transmitters
%kapparBS          = Hardware quality of the BSs' receivers
%
%OUTPUT:
%powerterms_MR    = 6 x K x L matrix with signal and interference powers
%                   for each of the UEs with MR combining.
%                   powerterms_MR(:,k,l) is the vector for UE k in cell l,
%                   where the first element is the average desired signal
%                   power, the second element is the interference from UEs
%                   having the same pilot, the third term is the
%                   interference from UEs having different pilots, the
%                   fourth term is the transmitter distortion from other
%                   UEs, the fifth term is the receiver distortion from all
%                   UEs, and the last term is the self-distortion/interference
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
C_totM = reshape(p*sum(sum(C,3),4),[M M L]);


%Prepare to store simulation results for the signal gains
signal_MR = zeros(K,L);
signal_RZF = zeros(K,L);
signal_MMMSE = zeros(K,L);

%Prepare to store simulation results for the norms of combining vectors
vectornorm_MR = zeros(K,L);
vectornorm_RZF = zeros(K,L);
vectornorm_MMMSE = zeros(K,L);

%Prepare to store simulation results for the sum of interference powers
interference_MR = zeros(K,L);
interference_RZF = zeros(K,L);
interference_MMMSE = zeros(K,L);

%Prepare to store simulation results for the sum of interference powers
%from UEs having the same pilot signals, including the user itself
reuseinterf_MR = zeros(K,L);
reuseinterf_RZF = zeros(K,L);
reuseinterf_MMMSE = zeros(K,L);

%Prepare to store simulation results for the self-interference power
selfinterf_MR = zeros(K,L);
selfinterf_RZF = zeros(K,L);
selfinterf_MMMSE = zeros(K,L);

%Prepare to store simulation results for the sum of receiver distortion
receiverdist_MR = zeros(K,L);
receiverdist_RZF = zeros(K,L);
receiverdist_MMMSE = zeros(K,L);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L

        %Extract channel realizations from all UEs to BS j
        Hallj = reshape(H(:,n,:,:,j),[M K*L]);
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Extract cells that use same pilots as cell j
        groupMembers = find(pilotPattern==pilotPattern(j))';
        
        
        %Compute three different combining schemes
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        V_RZF = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        V_MMMSE = p*(p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\V_MR;
        
        
        %Go through all UEs in cell j
        for k = 1:K
            
            
            %Extract channel realizations from UEs that cause pilot
            %contamination to UE k in cell j
            HjkPC = reshape(H(:,n,k,groupMembers,j),[M length(groupMembers)]);
            
            %Extract channel from BS j to its k:th UE
            Hjk = H(:,n,k,j,j);
            
        
            %%MR combining
            w = V_MR(:,k)/norm(V_MR(:,k))^2; %Extract combining vector
            wrep = repmat(w,[1 K*L]);
            
            %Use Monte Carlo simulations to compute the expectations that
            %appear in the six different components in Section 6.3.3
            signal_MR(k,j) = signal_MR(k,j) + (w'*Hjk)/nbrOfRealizations;
            vectornorm_MR(k,j) = vectornorm_MR(k,j) + norm(w).^2/nbrOfRealizations;
            interference_MR(k,j) = interference_MR(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
            receiverdist_MR(k,j) = receiverdist_MR(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
            reuseinterf_MR(k,j) = reuseinterf_MR(k,j) + p*sum(abs(w'*HjkPC).^2)/nbrOfRealizations;
            selfinterf_MR(k,j) = selfinterf_MR(k,j) + p*sum(abs(w'*Hjk).^2)/nbrOfRealizations;
            
           
            
            %%RZF combining
            w = V_RZF(:,k); %Extract combining vector
            wrep = repmat(w,[1 K*L]);
            
            %Use Monte Carlo simulations to compute the expectations that
            %appear in the six different components in Section 6.3.3
            signal_RZF(k,j) = signal_RZF(k,j) + (w'*Hjk)/nbrOfRealizations;
            vectornorm_RZF(k,j) = vectornorm_RZF(k,j) + norm(w).^2/nbrOfRealizations;
            interference_RZF(k,j) = interference_RZF(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
            receiverdist_RZF(k,j) = receiverdist_RZF(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
            reuseinterf_RZF(k,j) = reuseinterf_RZF(k,j) + p*sum(abs(w'*HjkPC).^2)/nbrOfRealizations;
            selfinterf_RZF(k,j) = selfinterf_RZF(k,j) + p*sum(abs(w'*Hjk).^2)/nbrOfRealizations;
            
            
            %%MMMSE combining
            w = V_MMMSE(:,k); %Extract combining vector
            wrep = repmat(w,[1 K*L]);
            
            %Use Monte Carlo simulations to compute the expectations that
            %appear in the six different components in Section 6.3.3
            signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*Hjk)/nbrOfRealizations;
            vectornorm_MMMSE(k,j) = vectornorm_MMMSE(k,j) + norm(w).^2/nbrOfRealizations;
            interference_MMMSE(k,j) = interference_MMMSE(k,j) + p*sum(abs(w'*Hallj).^2)/nbrOfRealizations;
            receiverdist_MMMSE(k,j) = receiverdist_MMMSE(k,j) + p*sum(sum(abs(wrep.*Hallj).^2,1))/nbrOfRealizations;
            reuseinterf_MMMSE(k,j) = reuseinterf_MMMSE(k,j) + p*sum(abs(w'*HjkPC).^2)/nbrOfRealizations;
            selfinterf_MMMSE(k,j) = selfinterf_MMMSE(k,j) + p*sum(abs(w'*Hjk).^2)/nbrOfRealizations;
            
        end
        
    end
    
end


%Set the normalized noise variance
sigma2 = 1;


%Compute the desired signal term, as described in Section 6.3.3
signalpower_MR = kappatUE*kapparBS*mean(mean(p*abs(signal_MR).^2./vectornorm_MR))/sigma2;
signalpower_RZF = kappatUE*kapparBS*mean(mean(p*abs(signal_RZF).^2./vectornorm_RZF))/sigma2;
signalpower_MMMSE = kappatUE*kapparBS*mean(mean(p*abs(signal_MMMSE).^2./vectornorm_MMMSE))/sigma2;

%Compute the interference from UEs having the same pilot, as described in
%Section 6.3.3
reuseinterference_MR = kappatUE*kapparBS*mean(mean((reuseinterf_MR-selfinterf_MR)./vectornorm_MR))/sigma2; 
reuseinterference_RZF = kappatUE*kapparBS*mean(mean((reuseinterf_RZF-selfinterf_RZF)./vectornorm_RZF))/sigma2; 
reuseinterference_MMMSE = kappatUE*kapparBS*mean(mean((reuseinterf_MMMSE-selfinterf_MMMSE)./vectornorm_MMMSE))/sigma2; 

%Compute the interference from UEs having different pilots, as described in
%Section 6.3.3
noreuseinterference_MR = kappatUE*kapparBS*mean(mean((interference_MR-reuseinterf_MR)./vectornorm_MR))/sigma2;
noreuseinterference_RZF = kappatUE*kapparBS*mean(mean((interference_RZF-reuseinterf_RZF)./vectornorm_RZF))/sigma2;
noreuseinterference_MMMSE = kappatUE*kapparBS*mean(mean((interference_MMMSE-reuseinterf_MMMSE)./vectornorm_MMMSE))/sigma2;

%Compute the transmitter distortion from other UEs, as described in
%Section 6.3.3
transmitDistortion_MR = (1-kappatUE)*kapparBS*mean(mean((interference_MR-selfinterf_MR)./vectornorm_MR))/sigma2;
transmitDistortion_RZF = (1-kappatUE)*kapparBS*mean(mean((interference_RZF-selfinterf_RZF)./vectornorm_RZF))/sigma2;
transmitDistortion_MMMSE = (1-kappatUE)*kapparBS*mean(mean((interference_MMMSE-selfinterf_MMMSE)./vectornorm_MMMSE))/sigma2;

%Compute the receiver distortion from all UEs, as described in
%Section 6.3.3
receiverDistortion_MR = (1-kapparBS)*mean(mean(receiverdist_MR./vectornorm_MR))/sigma2;
receiverDistortion_RZF = (1-kapparBS)*mean(mean(receiverdist_RZF./vectornorm_RZF))/sigma2;
receiverDistortion_MMMSE = (1-kapparBS)*mean(mean(receiverdist_MMMSE./vectornorm_MMMSE))/sigma2;

%Compute the self-distortion and self-interference, as described in
%Section 6.3.3
selfinterference_MR = kapparBS*mean(mean((selfinterf_MR-kappatUE*p*abs(signal_MR).^2)./vectornorm_MR))/sigma2; 
selfinterference_RZF = kapparBS*mean(mean((selfinterf_RZF-kappatUE*p*abs(signal_RZF).^2)./vectornorm_RZF))/sigma2; 
selfinterference_MMMSE = kapparBS*mean(mean((selfinterf_MMMSE-kappatUE*p*abs(signal_MMMSE).^2)./vectornorm_MMMSE))/sigma2; 


%Prepare to output the power values
powerterms_MR = [signalpower_MR reuseinterference_MR noreuseinterference_MR transmitDistortion_MR receiverDistortion_MR selfinterference_MR];
powerterms_RZF = [signalpower_RZF reuseinterference_RZF noreuseinterference_RZF transmitDistortion_RZF receiverDistortion_RZF selfinterference_RZF];
powerterms_MMMSE = [signalpower_MMMSE reuseinterference_MMMSE noreuseinterference_MMMSE transmitDistortion_MMMSE receiverDistortion_MMMSE selfinterference_MMMSE];

