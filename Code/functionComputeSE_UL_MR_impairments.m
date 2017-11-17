function [SE_MR,SE_MR_asymptotic] = functionComputeSE_UL_MR_impairments(channelGainOverNoise,tau_c,M,K,L,p,f,kappatUE,kapparBS)
%Compute UL SE for different receive combining schemes using Corollary 6.4.
%
%INPUT:
%channelGainOverNoise = K x L x L matrix containing the average channel
%                       gains normalized by the noise variance of al
%                       channels (in dB). channelGainOverNoise(k,j,l) is 
%                       the channel gain between between UE k in cell j 
%                       and BS l
%tau_c                = Length of coherence block
%M                    = Number of antennas per BS
%K                    = Number of UEs per cell
%L                    = Number of BSs and cells
%p                    = Uplink transmit power per UE (same for everyone)
%f                    = Pilot reuse factor
%kappatUE             = Hardware quality of the UEs' transmitters
%kapparBS             = Hardware quality of the BSs' receivers
%
%OUTPUT:
%SE_MR            = K x L matrix where element (k,l) is the downlink SE of
%                   UE k in cell l achieved with MR precoding
%SE_MR_asymptotic = K x L matrix where element (k,l) is the asymptotic
%                   downlink SE of UE k in cell l achieved with MR precoding
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


%Initiate the hardware qualities as perfect if these are not given as input
if nargin<8
    kappatUE = 1;
end

if nargin<9
    kapparBS = 1;
end


%Transform the (normalized) average channel gains from dB to linear scale
betaValues = 10.^(channelGainOverNoise/10);

%Compute length of pilot sequences
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


%Prepare to store signal and interference terms for all UEs
signalTerm = zeros(K,L);
interferenceTerm = zeros(K,L);

%If the asymptotic SE should be computed
if nargout>1
    
    signalTerm_asymptotic = zeros(K,L);
    interferenceTerm_asymptotic = zeros(K,L);
    
end


%% Go through all cells
for j = 1:L
    
    %Compute the parameter in (6.50)
    G = (1 + kapparBS*(M-1))/(M*kappatUE*kapparBS);
    
    
    %Go through all UEs
    for k = 1:K
        
        %Extract cells that use same pilots as cell j
        groupMembers = find(pilotPattern==pilotPattern(j))';
        
        %Compute the parameter in (6.38)
        psi_jk = 1/( p*kappatUE*kapparBS*tau_p*sum(betaValues(k,groupMembers,j)) + p*(1-kappatUE*kapparBS)*sum(sum(betaValues(:,:,j))) + 1 );
        
        %Compute the numerator of (6.39)
        signalTerm(k,j) = (p*betaValues(k,j,j))^2*tau_p*psi_jk*M;
        
        %Compute parameter in (6.40)
        Flijk = (1+p*betaValues(:,:,j)*psi_jk*( 1 - kappatUE*kapparBS + (1-kappatUE)*kapparBS^2*(M-1)  ))/(kappatUE*kapparBS)^2;
        
        %Compute the denominator of (6.39)
        interferenceTerm(k,j) = sum(sum(p*betaValues(:,:,j).*Flijk)) + p^2*sum(betaValues(k,groupMembers,j).^2)*tau_p*psi_jk*M*G - signalTerm(k,j) + 1/(kappatUE*kapparBS)^2;
        
        
        %If the asymptotic SE in (6.42) should be computed
        if nargout>1
            
            %Compute the numerator of (6.42)
            signalTerm_asymptotic(k,j) = (p*betaValues(k,j,j))^2;
            
            %Compute the denominator of (6.42)
            interferenceTerm_asymptotic(k,j) = p^2*sum(sum(betaValues(:,:,j).^2))*(1-kappatUE)/kappatUE^2/tau_p + p^2*sum(betaValues(k,groupMembers,j).^2)/kappatUE - (p*betaValues(k,j,j))^2;
            
        end
        
    end
    
end


%Compute the prelog factor assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/tau_c;

%Compute the final SE expression based on Corollary 6.3
SE_MR = prelogFactor * log2(1+signalTerm ./ real(interferenceTerm));

if nargout>1
    
    %Compute the final asymptotic SE expression based on Corollary 6.4
    SE_MR_asymptotic = prelogFactor * log2(1+signalTerm_asymptotic ./ real(interferenceTerm_asymptotic));
    
end
