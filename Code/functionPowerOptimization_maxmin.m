function SE = functionPowerOptimization_maxmin(signal,interference,Pmax,prelogFactor)
%Compute DL power allocation that solves the max-min fairness problem,
%using the algorithm in Theorem 7.1.
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
%     using the max-min power allocation solution
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
%This is version 1.01 (Last edited: 2019-03-16)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%Extract number of UEs
K = size(signal,1);

%Extract number of cells
L = size(signal,2);

%Check which UEs that have non-zero channels, because these ones are
%excluded (they are considered inactive)
nonzero = reshape(signal,[K*L 1]);
nonzero = nonzero(nonzero>0);

%Initalize the gamma-variables in Algorithm 1
rateLower = 0;
rateUpper = log2(1+Pmax*min(nonzero));

%Set the accuracy of the bisection
delta = 0.01;

%Prepare to save the power solution
rhoBest = zeros(K,L);

%Solve the max-min problem by bisection - iterate until the difference
%between the lower and upper points in the interval is smaller than delta
while norm(rateUpper - rateLower) > delta
    
    %Compute the midpoint of the line. Note that we are performing the
    %bisection in the SE domain instead of the SINR domain as in Algorithm
    %1, since we can then specify delta as the SE difference
    rateCandidate = (rateLower+rateUpper)/2; 
    
    %Transform the midpoints into SINR requirements
    gammaCandidate = 2.^(rateCandidate)-1;
    
    %Solve the feasibility problem using CVX
    [feasible,rhoSolution] = functionFeasibilityProblem_cvx(signal,interference,Pmax,K,L,gammaCandidate);
    
    
    %If the problem was feasible, then replace rateLower with
    %gammaCandidate and store rhoSolution as the new best solution.
    if feasible
        rateLower = rateCandidate;
        rhoBest = rhoSolution;
    else
        %If the problem was not feasible, then replace ratePoint with
        %gammaCandidate
        rateUpper = rateCandidate;
    end
    
end

%Compute the SEs using Theorem 4.6
SE = functionComputeSE_DL_poweralloc(rhoBest,signal,interference,prelogFactor);



function [feasible,rhoSolution] = functionFeasibilityProblem_cvx(signal,interference,Pmax,K,L,SINR)
%Solve the linear feasibility problem in Algorithm 1 using CVX, by adding
%an extra variable so that it becomes a minimization problem with better
%properties.

cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable rho(K,L);  %Variable for the K x L power allocation matrix
variable scaling    %Scaling parameter for power constraints

minimize scaling %Minimize the power indirectly by scaling the power constraints

subject to

for j = 1:L
    
    for k = 1:K
        
        if signal(k,j)>0
            
            %SINR constraints
            SINR*(sum(sum(rho.*interference(:,:,k,j))) + 1) - (rho(k,j)*signal(k,j)) <= 0
            
        end
        
        rho(k,j)>=0
        
    end
    
    sum(rho(:,j)) <= scaling*Pmax;
    
end

scaling >= 0; %Power constraints must be positive

cvx_end


%% Analyze the CVX output and prepare the output variables
if isempty(strfind(cvx_status,'Solved')) %Both the power minimization problem and the feasibility problem are infeasible
    feasible = false;
    rhoSolution = [];
elseif scaling>1 %Only the power minimization problem is feasible
    feasible = false;
    rhoSolution = rho;
else %Both the power minimization problem and feasibility problem are feasible
    feasible = true;
    rhoSolution = rho;
end
