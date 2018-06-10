function SE = functionPowerOptimization_prodSINR(signal,interference,Pmax,prelogFactor)
%Compute DL power allocation that solves the max product SINR problem,
%using the algorithm in Theorem 7.2.
%
%This function require additional software packages to be used, which
%need to be downloaded and installed separately. These packages are
%developed independently and are delivered with separate licenses.
%The implementation uses CVX (http://cvxr.com/cvx) and has been tested
%using CVX version 2.1. We recommend the use of the Mosek solver (we
%have tested using version 7.1.0.12).
%
%INPUT:
%signal       = K x L matrix where element (k,j) is a_jk in (7.2)
%interference = K x L x K x L matrix where (l,i,j,k) is b_lijk in (7.3)
%Pmax         = Maximum transmit power per BS
%prelogFactor = Prelog factor
%
%OUTPUT:
%SE = K x L matrix where element (k,j) is the downlink SE of UE k in cell j
%     using the max product power allocation solution
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
%This is version 1.1 (Last edited: 2018-06-08)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);


%% Solve geometric program in (7.8) using CVX
cvx_begin gp
cvx_quiet(true); % This suppresses screen output from the solver

variable rho(K,L);
variable c(K,L);

maximize prod(prod(c))

subject to

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j)>0
            %SINR constraints of UE k in cell j
            c(k,j)*(sum(sum(rho.*interference(:,:,k,j))) + 1) <= (rho(k,j)*signal(k,j));
            
            rho(k,j) >= 0;
            
        else
            %This applies if UE k in cell j is inactive
            c(k,j) == 1;
            rho(k,j) >= 0;
            
        end
        
    end
    
    sum(rho(:,j)) <= Pmax;
    
end

cvx_end

%% Analyze the CVX output and prepare the output variables
if isempty(strfind(cvx_status,'Solved')) %The problem was not solved by CVX, for some reason, and we then consider equal power allocation
    rhoSolution = (Pmax/K)*ones(K,L);
else %The problem was solved by CVX
    rhoSolution = rho;
end

%Refine the solution obtained from CVX using the Matlab command fmincon.
%This is particularly important in case CVX fails to solve the problem
A = kron(eye(L),ones(1,K));
B = Pmax*ones(L,1);
options = optimoptions('fmincon','Algorithm','interior-point','MaxFunEvals',50000,'MaxIter',5000);
xend = fmincon(@(x) -functionComputesumSE_DL_poweralloc(x,signal,interference,prelogFactor),rhoSolution(:),A,B,[],[],zeros(K*L,1),[],[],options);
rhoBest = reshape(xend,[K L]);

%Compute the SEs using Theorem 4.6
SE = functionComputeSE_DL_poweralloc(rhoBest,signal,interference,prelogFactor);


function sumSE = functionComputesumSE_DL_poweralloc(rho,signal,interference,prelogFactor)

%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);

%Reshape power variable since fmincon optimizes vectors
rho = reshape(rho,[K L]);

%Prepare to save results
SE = zeros(K,L);

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j) > 0 %Check if the UE k in cell j is active
            
            %Compute the SE in Theorem 4.6 using the formulation in (7.1)
            SE(k,j) = prelogFactor*log2((rho(k,j)*signal(k,j)) / (sum(sum(rho.*interference(:,:,k,j))) + 1));
            
        else %If the UE is inactive
            
            SE(k,j) = 0;
            
        end
        
    end
    
end

%Compute the sum SE of all cells
sumSE = sum(SE(:));
