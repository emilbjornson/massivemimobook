function [H,Rest,activeUEs] = functionExampleSetup_Quadriga(L,Kdrop,B,noiseVariancedBm,Kmax,f,M,polarizations)
%This function generates the channel realizations between UEs at random
%locations and the BSs in the running example, defined in Section 4.1.3.
%BSs with cylindrical arrays and channel are generated using QuaDRiGa from
%the Fraunhofer Heinrich Hertz Institute. The Urban Microcell NLOS scenario
%is used for channel modeling. QuaDRiGa needs to be installed separately
%(http://www.quadriga-channel-model.de) and is delivered with a separate
%license. This function has been tested using QuaDRiGa version 1.4.8-571.
%
%INPUT:
%L                = Number of BSs and cells
%Kdrop            = Number of UEs to be dropped in the square around a BS
%B                = Bandwidth in Hz
%noiseVariancedBm = Noise variance in dBm
%Kmax             = Maximum number of UEs served by a BS
%f                = Pilot reuse factor, giving pilot length Kmax*f
%M                = Number of BS antennas
%polarizations    = Select number of antenna polarizations (1 or 2)
%
%OUTPUT:
%H         = M x 400 x K x L x L matrix with the channel realizations over
%            400 subcarriers at one time instance
%Rest      = M x M x K x L x L matrix with estimates of the spatial
%            correlation matrices for all UEs in the network.
%            Rest(:,:,k,j,l) is the correlation matrix for the channel
%            between UE k in cell j and the BS in cell l.
%activeUEs = Kmax x L with zeros and ones. activeUEs(k,l)==1 means that
%            pilot k is used by a UE in cell l
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


%% Create new Quadriga layout

%Set irrelevant parameters
s = simulation_parameters;
s.sample_density = 1;
s.use_absolute_delays = 1;

%Set center frequency
center_frequency = 2e9;
s.center_frequency = center_frequency;

%Number of subcarriers
nbrOfSubcarriers = 400;

%Create new layout from general parameters
lay = layout(s);

%Generate BSs
lay.no_tx = L;

%Set BS heights
lay.tx_position(3,:) = 25;

%Set the length in meters of the total square area
squareLength = 1000;

%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);


%% Deploy BSs

%Minimum distance between BSs and UEs
minDistance = 35;

%Distance between BSs in vertical/horizontal direction
interSiteDistance = squareLength/nbrBSsPerDim;

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);

for j = 1:length(BSpositions)
    
    lay.tx_position(1:2,j) = [real(BSpositions(j)); imag(BSpositions(j))];
    
end


%% Create a circular antenna array for each BS
M_V = 5; %Number of vertical antennas
M_H = M/M_V; %Number of antennas on each horizontal circle

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

if polarizations == 1
    PolarizationIndicator = 1; %Vertical polarization
elseif polarizations == 2
    PolarizationIndicator = 3; %Dual +/-45deg polarization
end

%Compute height of array
arrayHeight = (M_V-1)*antennaSpacing*3e8/center_frequency;

if (PolarizationIndicator==1) %Single polarization
    
    circumference = M_H*antennaSpacing*3e8/center_frequency;
    radius = circumference/(2*pi);
    delta_angle = 2*pi/M_H;
    
    %Go through all BSs
    for b = 1:L
        
        %Create rectangular array of size M_V x M_H
        lay.tx_array(b).generate('3gpp-3d', 1, M_V, M_H, center_frequency, PolarizationIndicator, 0, antennaSpacing);
        
        %Place antennas on a circle and rotate radiation patters
        for i = 1:M_V
            for j = 1:M_H
                indices = (i-1)*M_H + j;
                angle = (j-1)*delta_angle;
                lay.tx_array(b).element_position(1, indices) = radius*cos(angle);
                lay.tx_array(b).element_position(2, indices) = radius*sin(angle);
                lay.tx_array(b).element_position(3, indices) = (i-1)*antennaSpacing*3e8/center_frequency - arrayHeight/2;
                lay.tx_array(b).rotate_pattern(rad2deg(angle), 'z', indices, 0);
            end
        end
    end
    
elseif (PolarizationIndicator==3) %Dual polarization
    
    circumference = M_H/2*antennaSpacing*3e8/center_frequency;
    radius = circumference/(2*pi);
    delta_angle = 2*pi/(M_H/2);
    
    %Go through all BSs
    for b = 1:L
        
        %Create rectangular array of size M_V x M_H/2
        lay.tx_array(b).generate('3gpp-3d', 1, M_V, M_H/2, center_frequency, PolarizationIndicator, 0, antennaSpacing);
        
        %Place antennas on a circle and rotate radiation patters (while keeping co-located antennas together)
        for i = 1:M_V
            for j = 1:M_H/2
                indices = (i-1)*M_H + 2*j-1 : (i-1)*M_H + 2*j;
                angle = (j-1)*delta_angle;
                lay.tx_array(b).element_position(1, indices) = radius*cos(angle);
                lay.tx_array(b).element_position(2, indices) = radius*sin(angle);
                lay.tx_array(b).element_position(3, indices) = (i-1)*antennaSpacing*3e8/center_frequency - arrayHeight/2;
                lay.tx_array(b).rotate_pattern(rad2deg(angle), 'z', indices, 0);
            end
        end
    end
end

%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';

%Compute the exact dimension of the square where the users are located
maxDistance = interSiteDistance;



%Prepare to put out UEs in the cells
UEpositions = zeros(Kdrop,L);
UEpositionsWrapped = zeros(Kdrop,L,length(wrapLocations));
perBS = zeros(L,1);

%Go through all the cells
for l = 1:L
    
    %Put out K UEs in the cell, uniformly at random. The procedure is
    %iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs
    while perBS(l)<Kdrop
        
        %Put out users
        UEremaining = Kdrop-perBS(l);
        posX = rand(UEremaining,1)*maxDistance - maxDistance/2;
        posY = rand(UEremaining,1)*maxDistance - maxDistance/2;
        posXY = posX + 1i*posY;
        
        %Keep those that satisfy the minimum distance
        posXY = posXY(abs(posXY)>=minDistance);
        
        %Store new UEs
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
        
    end
    
    %Create alternative UE positions using wrap around
    for k = 1:Kdrop
        
        UEpositionsWrapped(k,l,:) = UEpositions(k,l) + wrapLocations;
        
    end
    
end



%% Configure UEs

Ktotal = Kdrop*L*length(wrapLocations); %Total number of UEs

%Define UE heights
UE_heights = 1.5*ones(Kdrop,L);
UE_heightsWrapped = repmat(UE_heights,[1 1 length(wrapLocations)]);

%Generate UEs
lay.no_rx = Ktotal;

%Define UE antennas
lay.rx_array.generate('omni');


%% Simulate channels

% Randomly distribute UEs
lay.rx_position =  [real(UEpositionsWrapped(:))'; imag(UEpositionsWrapped(:))'; UE_heightsWrapped(:)'];

%Define tracks for each UE, assuming a fixed UE location
for k=1:Ktotal
    lay.track(k).generate('linear',0,0) %Define a linear track consisting of only one position
    lay.track(k).scenario = '3GPP_3D_UMi_NLOS'; %Select the Urban Microcell NLOS scenario
end

%Generate pilot patterns
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2 %Only works for 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
    
end


%Randomize pilot allocation in each cell
randOrder = zeros(Kmax*f,L);

for j = 1:L
    
    randOrder(1+(pilotPattern(j)-1)*Kmax:pilotPattern(j)*Kmax,j) = randperm(Kmax)+(pilotPattern(j)-1)*Kmax;
    
end


%Compute variance and standard deviation of the noise
noiseVar = 10^(noiseVariancedBm/10);
noiseStd = sqrt(noiseVar);


%Prepare to store channel realizations
H = zeros(M,nbrOfSubcarriers,Kmax,L,L);
Rest = zeros(M,M,Kmax,L,L);
perBS = zeros(L,1);
activeUEs = zeros(Kmax,L);


%% Go through all cells
for j = 1:L
    
    %Output simulation progress
    disp([num2str(j) ' cells generated out of ' num2str(L)]);
    
    %Go through all UEs
    for k = 1:Kdrop
        
        Huser = zeros(M,nbrOfSubcarriers,1,1,L);
        Ruser = zeros(M,M,1,1,L);
        
        %Extract the channels to all BSs
        for l = 1:L
            
            [~,minr] = min(abs(UEpositionsWrapped(k,j,:)-BSpositions(l)));
            
            userind = k+(j-1)*Kdrop+(minr-1)*Kdrop*L;
            
            [ h_channel, h_cb ] = lay.get_channels_seg(l, userind);
            Hextract = h_channel.fr(B, nbrOfSubcarriers);
            
            Huser(:,:,1,1,l) = reshape(Hextract,[M nbrOfSubcarriers])/noiseStd;
            Ruser(:,:,1,1,l) = diag(mean(abs(Huser(:,:,1,1,l)).^2,2)/noiseVar);
            
        end
        
        %Determine which BS should serve the UE
        [~,bestBS] = max(mean(sum(abs(Huser(:,:,1,1,:)).^2,1),2));
        
        %Check if the selected BS has pilots available
        if perBS(bestBS)<Kmax
            
            %Add the UE to the cell of the selected BS
            perBS(bestBS) = perBS(bestBS) + 1;
            H(:,:,randOrder(perBS(bestBS)+(pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS,:) = Huser;
            Rest(:,:,randOrder(perBS(bestBS)+(pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS,:) = Ruser;
            activeUEs(randOrder(perBS(bestBS)+(pilotPattern(bestBS)-1)*Kmax,bestBS),bestBS) = 1;
            
        end
        
    end
    
end
