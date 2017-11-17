function channelGaindB = functionExampleSetup_backhaul(S,L)
%This function generates the channel gains between SBSs at random
%locations on a 9x9 grid in L cells and the Massive MIMO BSs, for the
%scenario described in Section 7.6.
%
%INPUT:
%S             = Number of randomly selected SBSs per cell out of 81
%L             = Total number of cells and BSs
%
%OUTPUT:
%channelGaindB = L x L x S matrix containing the average channel gains in
%                dB of all the channels between S. channelGaindB(l,j,s)
%                corresponds to the average channel gain from SBS s in cell
%                j to BS l
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


%Deploy 81 SBSs on a regular grid in each cell
SBSpositions = zeros(81,L);
for b = 1:L
    
    counter = 0;
    
    for i = 1:9
        
        for j = 1:9
            
            counter = counter + 1;
            SBSpositions(counter,b) = interBSDistance*(-0.5 + (i-1)*1/9 + 1i*(-0.5 + (j-1)*1/9))  + BSpositions(b);

        end
        
    end
    
end


%Randomly select K SBSs within each cell
SELECTION = zeros(S,L);
for b=1:L
    tmp = randperm(81);
    SELECTION(:,b) = SBSpositions(tmp(1:S),b);
end

SBSpositions = SELECTION;


%Prepare to store average channel gains (in dB)
channelGaindB = zeros(S,L,L);


%% Go through all the cells
for l = 1:L
    
    
    %Go through all BSs
    for j = 1:L
        
        %Compute the distance from the SBSs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a SBS and the
        %nine different locations of a BS is considered 
        distancesBSj = min(abs( repmat(SBSpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[S 1]) ),[],2);
        
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading
        channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
        
    end
    
    %Go through all SBSs in cell l and generate shadow fading
    for k = 1:S
        
        %Generate shadow fading realization
        shadowing = sigma_sf*randn(1,1,L);
        channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        
        %Check if another BS has a larger average channel gain to the SBS
        %than BS l
        while channelGainShadowing(l) < max(channelGainShadowing)
            
            %Generate new shadow fading realizations (until all SBSs in 
            %cell l has its largest average channel gain from BS l)
            shadowing = sigma_sf*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            
        end
        
        %Store average channel gains with shadowing fading
        channelGaindB(k,l,:) = channelGainShadowing;
        
    end
    
end

%Permute dimensions to obtain the designated output
channelGaindB = permute(channelGaindB,[3,2,1]);
