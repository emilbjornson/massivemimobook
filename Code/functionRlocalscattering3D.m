function R = functionRlocalscattering3D(M_H, M_V, d_H, d_V, varphi, stdphi, theta, stdtheta, dist)
%Generate the spatial correlation matrix for the 3D local scattering model,
%defined in (7.17) for different angular distributions.
%
%INPUT:
% M_H           = Number of antennas per horizontal row, i.e., # of columns
% M_V           = Number of horizontal rows
% d_H           = Horizontal antenna spacing in multiples of wavelength
% d_V           = Vertical antenna spacing in multiples of wavelength
% varphi        = Azimuth angle in radians
% stdphi        = Azimuth angular spread (standard deviation or limit)
% theta         = Elevation angle in radians
% stdtheta      = Elevation angular spread (standard deviation or limit)
% dist          = 'Gaussian' or 'Laplace' or 'Uniform'
%
%OUTPUT:
% R             = M x M channel covariance matrix, where M = M_H x M_V
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


% Index mapping functions
i = @(m) mod(m-1, M_H);
j = @(m) floor((m-1)/M_H);

%Create the correlation matrix
M = M_H*M_V;
R = zeros(M);

% Define angle distributions as functions
if (strcmp(dist,'Laplace'))
    f1 = @(x) exp(-sqrt(2)*abs(x-varphi)/stdphi)/(sqrt(2)*stdphi);
    f2 = @(x) exp(-sqrt(2)*abs(x-theta)/stdtheta)/(sqrt(2)*stdtheta);
elseif (strcmp(dist,'Gaussian'))
    f1 = @(x) exp(-(x-varphi).^2/(2*stdphi^2))/(sqrt(2*pi)*stdphi);
    f2 = @(x) exp(-(x-theta).^2/(2*stdtheta^2))/(sqrt(2*pi)*stdtheta);
elseif (strcmp(dist,'Uniform'))
    f1 =  @(x)1/(2*stdphi);
    f2 =  @(x) 1/(2*stdtheta);
else
    error('Please provide a valid angle distribution')
end

%% Compute correlation matrix
LOOKUP = zeros(2*M_H-1,M_V); %Define a lookup table to speedup the computation process
for m = 1:M
    for l = 1:m
        indi = i(l)-i(m) + M_H;
        indj = j(l)-j(m) + M_V;
        if (LOOKUP(indi,indj)==0)
            integrand = @(p,t) exp(1j*2*pi*d_V*(j(l)-j(m)).*sin(t)) .* exp(1j*2*pi*d_H*(i(l)-i(m)).*cos(t).*sin(p)) .* f1(p) .* f2(t);
            if (strcmp(dist,'Uniform'))
                LOOKUP(indi,indj) = integral2(integrand, varphi-stdphi, varphi+stdphi, theta-stdtheta, theta+stdtheta);
            else
                LOOKUP(indi,indj) = integral2(integrand, varphi-20*stdphi, varphi+20*stdphi, theta-20*stdtheta, theta+20*stdtheta);
            end
        end
        R(m,l) = LOOKUP(indi,indj);
        R(l,m) = R(m,l)';
    end
end
