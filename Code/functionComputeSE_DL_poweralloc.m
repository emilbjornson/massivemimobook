function SE = functionComputeSE_DL_poweralloc(rho,signal,interference,prelogFactor)
%Compute the SE in Theorem 4.6 using the formulation in (7.1) for a given
%power allocation scheme.
%
%INPUT:
%rho          = K x L matrix where element (k,j) is the downlink transmit
%               power allocated to UE k in cell j
%signal       = K x L matrix where element (k,j,n) is a_jk in (7.2)
%interference = K x L x K x L matrix where (l,i,jk,n) is b_lijk in (7.3)
%prelogFactor = Prelog factor
%
%OUTPUT:
%SE = K x L matrix where element (k,j) is the downlink SE of UE k in cell j
%     using the power allocation given as input
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


%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);

%Prepare to save results
SE = zeros(K,L);


%% Go through all cells
for j = 1:L
    
    %Go through all UEs in cell j
    for k = 1:K
        
        %Compute the SE in Theorem 4.6 using the formulation in (7.1)
        SE(k,j) = prelogFactor*log2(1+(rho(k,j)*signal(k,j)) / (sum(sum(rho.*interference(:,:,k,j))) + 1));
        
    end
    
end
