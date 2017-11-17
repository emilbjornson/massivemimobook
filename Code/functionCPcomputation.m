function [P_MR,P_RZF,P_MMMSE,P_ZF,P_SMMSE] = functionCPcomputation(Mrange,K,L,B,tau_c,tau_p,valueset,sumSE_MR,sumSE_RZF,sumSE_MMMSE,sumSE_ZF,sumSE_SMMSE)
%This function returns the total CP with different processing schemes,
%using the model defined Section 5.4 and one of the value sets in Table
%5.3.
%
%INPUT:
%Mrange      = Vector with different number of BS antennas
%K           = Number of UEs per cell
%L           = Number of cells
%sumSE_MR    = Average sum SE per cell with MR (same size as Mrange)
%sumSE_RZF   = Average sum SE per cell with RZF (same size as Mrange)
%sumSE_MMMSE = Average sum SE per cell with M-MMSE (same size as Mrange)
%sumSE_ZF    = Average sum SE per cell with ZF (same size as Mrange)
%sumSE_SMMSE = Average sum SE per cell with S-MMSE (same size as Mrange)
%B           = Bandwidth
%tau_c       = Length of coherence block
%tau_p       = Length of pilot sequence
%valueset    = Select which one of the value sets in Table 5.3 that is
%              considered. Either valueset=1 or valueset=2 are possible.
%
%OUTPUT:
%P_MR    = Vector of same size as Mrange with total CP for MR processing
%P_RZF   = Same as P_MR but with RZF processing
%P_MMMSE = Same as P_MMMSE but with M-MMSE processing
%P_ZF    = Same as P_ZF but with XF processing
%P_SMMSE = Same as P_SMMSE but with S-MMSE processing
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


%Obtain CP model coefficients for one of the value sets in Table 5.3
[P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(valueset);

%Prepare to store simulation results
P_MR = zeros(length(Mrange),1);
P_TC = zeros(length(Mrange),1);
P_CE = zeros(length(Mrange),1);
P_SP_RT = zeros(length(Mrange),1);
P_SP_DL = zeros(length(Mrange),1);
P_SP_UL_MR = zeros(length(Mrange),1);

%Compute CP for coding and decoding using (5.35)
P_CD_MR = (P_COD + P_DEC)*B*sumSE_MR;

%Compute CP for backhaul traffic using (5.36)
P_BH_MR = P_BT*B*sumSE_MR;


%Repeat computations for RZF
if nargin>8
    P_RZF = zeros(length(Mrange),1);
    P_CD_RZF = (P_COD + P_DEC)*B*sumSE_RZF;
    P_BH_RZF = P_BT*B*sumSE_RZF;
    P_SP_UL_RZF = zeros(length(Mrange),1);
end

%Repeat computations for M-MMSE
if nargin>9
    P_MMMSE = zeros(length(Mrange),1);
    P_CD_MMMSE = (P_COD + P_DEC)*B*sumSE_MMMSE;
    P_BH_MMMSE = P_BT*B*sumSE_MMMSE;
    P_SP_UL_MMMSE = zeros(length(Mrange),1);
end

%Repeat computations for ZF
if nargin>10
    P_ZF = zeros(length(Mrange),1);
    P_CD_ZF = (P_COD + P_DEC)*B*sumSE_ZF;
    P_BH_ZF = P_BT*B*sumSE_ZF;
    P_SP_UL_ZF = zeros(length(Mrange),1);
end

%Repeat computations for S-MMSE
if nargin>11
    P_SMMSE = zeros(length(Mrange),1);
    P_CD_SMMSE = (P_COD + P_DEC)*B*sumSE_SMMSE;
    P_BH_SMMSE = P_BT*B*sumSE_SMMSE;
    P_SP_UL_SMMSE = zeros(length(Mrange),1);
end


%Go through all number of antennas
for index = 1:length(Mrange)
    
    %Extract current number of antennas
    M = Mrange(index);
    
    %Compute CP for transceiver chains using (5.34)
    P_TC(index) = M*P_BS + P_LO + K*P_UE;
    
    %Compute CP for channel estimation with all other schemes, where only
    %the channels to UEs in other cells are estimated, using (5.37)
    P_CE(index) = 3*K*B/(tau_c*L_BS)*(M*tau_p + M^2);
    
    %Compute CP for UL reception and DL transmission
    P_SP_RT(index) = (tau_c - tau_p)*3*B/(tau_c*L_BS)*M*K;
    
    %Compute CP for computation of precoding vectors
    P_SP_DL(index) = 4*M*K*B/(tau_c*L_BS);
    
    %Sum up the power terms that are independent of the processing scheme
    P_SAME = P_FIX + P_TC(index) + P_SP_RT(index) + P_SP_DL(index);
    
    
    %Compute CP for computation of the combining vectors with different
    %schemes, based on Table 5.2
    P_SP_UL_MR(index) = 7*B*K/(tau_c*L_BS);
    
    %Compute the final CP values
    P_MR(index) = P_SAME + P_CE(index) + P_CD_MR(index) + P_BH_MR(index) + P_SP_UL_MR(index);
    
    %Repeat same computations for RZF
    if nargin>8
        P_SP_UL_RZF(index) = 3*B*(3*K^2*M/2 + 3*M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
        P_RZF(index) = P_SAME + P_CE(index) + P_CD_RZF(index) + P_BH_RZF(index) + P_SP_UL_RZF(index);
    end
    
    %Repeat same computations for M-MMSE
    if nargin>9
        P_SP_UL_MMMSE(index) = 3*B*(L*(3*M^2 + M)*K/2 + M^3/3 + 2*M + M*tau_p*(tau_p-K))/(tau_c*L_BS);
        P_MMMSE(index) = P_SAME + P_CE(index) + P_CD_MMMSE(index) + P_BH_MMMSE(index) + P_SP_UL_MMMSE(index);
    end
    
    %Repeat same computations for ZF
    if nargin>10
        P_SP_UL_ZF(index) = 3*B*(3*K^2*M/2 + M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
        P_ZF(index) = P_SAME + P_CE(index) + P_CD_ZF(index) + P_BH_ZF(index) + P_SP_UL_ZF(index);
    end
    
    %Repeat same computations for S-MMSE
    if nargin>11
        P_SP_UL_SMMSE(index) = 3*B*(3*M^2*K/2 + M*K/2 + (M^3 - M)/3 + (7/3)*M)/(tau_c*L_BS);
        P_SMMSE(index) = P_SAME + P_CE(index) + P_CD_SMMSE(index) + P_BH_SMMSE(index) + P_SP_UL_SMMSE(index);
    end
    
end
