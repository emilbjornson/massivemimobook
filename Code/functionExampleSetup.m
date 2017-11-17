function [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg)
%This function generates the channel statistics between UEs at random
%locations and the BSs in the running example, defined in Section 4.1.3.
%
%Note that all distances in this code are measured in meters, while
%distances is often measured in kilometer in the monograph.
%
%INPUT:
%L        = Number of BSs and cells
%K        = Number of UEs per cell
%M        = Number of antennas per BS
%accuracy = Compute exact correlation matrices from the local scattering
%           model if approx=1. Compute a small-angle approximation of the
%           model if approx=2
%ASDdeg   = Angular standard deviation around the nominal angle
%           (measured in degrees)
%
%OUTPUT:
%R             = M x M x K x L x L matrix with spatial correlation matrices
%                for all UEs in the network. R(:,:,k,j,l) is the correlation
%                matrix for the channel between UE k in cell j and the BS
%                in cell l. This matrix is normalized such that trace(R)=M.
%channelGaindB = K x L x L matrix containing the average channel gains in
%                dB of all the channels. The product 
%                R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10) is the full
%                spatial channel correlation matrix.
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


%% Model parameters

%Set the length in meters of the total square area
squareLength = 1000;

%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);

%Pathloss exponent
alpha = 3.76;

%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

%Standard deviation of shadow fading
sigma_sf = 10;

%Minimum distance between BSs and UEs
minDistance = 35;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Prepare to put out UEs in the cells
UEpositions = zeros(K,L);
perBS = zeros(L,1);

%Prepare to store normalized spatial correlation matrices
R = zeros(M,M,K,L,L,length(ASDdeg));

%Prepare to store average channel gain numbers (in dB)
channelGaindB = zeros(K,L,L);


%% Go through all the cells
for l = 1:L
    
    %Put out K UEs in the cell, uniformly at random. The procedure is
    %iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs
    while perBS(l)<K
        
        %Put out new UEs
        UEremaining = K-perBS(l);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;
        
        %Keep those that satisfy the minimum distance
        posXY = posXY(abs(posXY)>=minDistance);
        
        %Store new UEs
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
        
    end
    
    
    %Go through all BSs
    for j = 1:L
        
        %Compute the distance from the UEs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a UE and the
        %nine different locations of a BS is considered 
        [distancesBSj,whichpos] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading 
        channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
        
        %Compute nominal angles between UE k in cell l and BS j, and
        %generate spatial correlation matrices for the channels using the
        %local scattering model
        for k = 1:K
            
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,whichpos(k)));
            
            if accuracy == 1 %Use the exact implementation of the local scattering model
                
                for spr = 1:length(ASDdeg)
                    
                    R(:,:,k,l,j,spr) = functionRlocalscattering(M,angleBSj,ASDdeg(spr),antennaSpacing);
                    
                end
                
            elseif accuracy == 2 %Use the approximate implementation of the local scattering model
                
                for spr = 1:length(ASDdeg)
                    
                    R(:,:,k,l,j,spr) = functionRlocalscatteringApprox(M,angleBSj,ASDdeg(spr),antennaSpacing);
                    
                end
                
            end
            
        end
        
    end
    
    
    %Go through all UEs in cell l and generate shadow fading realizations
    for k = 1:K
        
        %Generate shadow fading realizations
        shadowing = sigma_sf*randn(1,1,L);
        channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        
        %Check if another BS has a larger average channel gain to the UE
        %than BS l
        while channelGainShadowing(l) < max(channelGainShadowing)
            
            %Generate new shadow fading realizations (until all UEs in cell
            %l has its largest average channel gain from BS l)
            shadowing = sigma_sf*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            
        end
        
        %Store average channel gains with shadowing fading
        channelGaindB(k,l,:) = channelGainShadowing;
        
    end
    
end
