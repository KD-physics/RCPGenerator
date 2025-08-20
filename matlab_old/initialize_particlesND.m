function [x, D] = initialize_particlesND(phi, N, Box, distribution)
%INITIALIZE_PARTICLESND Initializes particle positions and diameters in ND.
%
% INPUTS:
%   phi         : Target packing fraction (0 < phi < 1).
%   N           : Number of particles.
%   Box         : [1 x Ndim] array of box dimensions in each dimension.
%   distribution: Structure specifying diameter distribution parameters.
%                 If empty, prints available distributions and their tags.
%
% OUTPUTS:
%   x           : [N x Ndim] array of particle positions in ND.
%   D           : [N x 1] array of particle diameters scaled to the target phi.
%
% Example distribution structure:
%   distribution.type = 'Gaussian'; % Distribution type
%   distribution.mu = 1;            % Mean diameter (for Gaussian)
%   distribution.sigma = 0.2;       % Standard deviation (for Gaussian)

% Check if distribution is empty and print options
if nargin < 4 || isempty(distribution)
    fprintf('Available distributions and required parameters:\n');
    fprintf('- Mono: d (single diameter value)\n');
    fprintf('- Bidisperse: d1, d2, p (fraction of d1 particles)\n');
    fprintf('- Gaussian: mu (mean), sigma (standard deviation)\n');
    fprintf('- BiGaussian: mu1, sigma1, mu2, sigma2, p (fraction for first Gaussian)\n');
    fprintf('- Lognormal: mu (mean), sigma (log-scale std deviation)\n');
    fprintf('- Flat: d_min (min diameter), d_max (max diameter)\n');
    fprintf('- PowerLaw: d_min (min diameter), d_max (max diameter), exponent\n');
    fprintf('- Exponential: d_min (min diameter), d_max (max diameter)\n');
    fprintf('- Weibull: scale, shape\n');
    fprintf('- Custom: Provide a vector D of size [1, N]\n');
    return;
end

% Extract dimensionality and box dimensions
Ndim = length(Box);

% Generate diameters based on the specified distribution
D = generate_diameter_distribution(N, distribution);
D = D(:);

% Scale diameters to achieve the target packing fraction
D = scale_diametersND(D, phi, Box);

% Initialize positions
x = rand(N, Ndim) .* Box; % Random initial positions
count = 1;

% Ensure particles are placed without overlap
while count < N
    % Generate a new random position
    x_tmp = rand(1, Ndim) .* Box;

    % Compute distances to existing particles
    d = sqrt(sum((x(1:count, :) - x_tmp).^2, 2));
    r = (D(1:count) + D(count + 1)) / 2;

    % Check if the new position is valid
    if all(d > r)
        count = count + 1;
        x(count, :) = x_tmp;
    end
end

end

function D = generate_diameter_distribution(N, distribution)
%GENERATE_DIAMETER_DISTRIBUTION Generates particle diameters based on a specified distribution.
%
% INPUTS:
%   N           : Number of particles.
%   distribution: Structure specifying distribution type and parameters.
%
% OUTPUTS:
%   D           : [1 x N] array of particle diameters.

% Parse distribution type
dist_type = lower(distribution.type);

switch dist_type
    case 'mono'
        % Monodisperse: single diameter
        D = ones(1, N) * distribution.d;

    case 'gaussian'
        % Gaussian distribution: mu (mean), sigma (std dev)
        D = abs(normrnd(distribution.mu, distribution.sigma, 1, N));

    case 'bigaussian'
        % Bi-Gaussian distribution: two Gaussians with weights
        mu1 = distribution.mu1;
        sigma1 = distribution.sigma1;
        mu2 = distribution.mu2;
        sigma2 = distribution.sigma2;
        p = distribution.p; % Fraction for first Gaussian
        N1 = round(p * N);
        N2 = N - N1;
        D1 = abs(normrnd(mu1, sigma1, 1, N1));
        D2 = abs(normrnd(mu2, sigma2, 1, N2));
        D = [D1, D2];
        D = D(randperm(N)); % Shuffle

    case 'bidisperse'
        % Bi-disperse distribution: two fixed diameters and a fraction
        d1 = distribution.d1; % Smaller diameter
        d2 = distribution.d2; % Larger diameter
        p = distribution.p;  % Fraction of particles with d1
        N1 = round(p * N);   % Number of particles with d1
        N2 = N - N1;         % Number of particles with d2
        D = [repmat(d1, 1, N1), repmat(d2, 1, N2)];
        D = D(randperm(N)); % Shuffle to randomize positions

    case 'lognormal'
        % Lognormal distribution: mu (mean), sigma (std dev in log space)
        D = lognrnd(distribution.mu, distribution.sigma, 1, N);

    case 'flat'
        % Flat distribution: d_min, d_max
        D = distribution.d_min + (distribution.d_max - distribution.d_min) * rand(1, N);

    case 'powerlaw'
        % Power-law distribution: d_min, d_max, exponent
        d_min = distribution.d_min;
        d_max = distribution.d_max;
        exponent = distribution.exponent;
        u = rand(1, N);
        %D = ((d_max^(exponent + 1) - d_min^(exponent + 1)) * u + d_min^(exponent + 1)).^(1 / (exponent + 1));
        if exponent == -1
            % Special case for exponent = -1
            D = d_min * exp(u * log(d_max / d_min));
        else
            % General power-law case
            D = ((d_max^(exponent + 1) - d_min^(exponent + 1)) * u + d_min^(exponent + 1)).^(1 / (exponent + 1));
        end

    case 'exponential'
        % Exponential distribution: d_min, d_max
        d_min = distribution.d_min;
        d_max = distribution.d_max;
        lambda = -log(0.5) / (d_max - d_min); % Decay rate for the range
        D = d_min + (-1 / lambda) * log(1 - rand(1, N));

    case 'weibull'
        % Weibull distribution: scale, shape
        scale = distribution.scale;
        shape = distribution.shape;
        D = wblrnd(scale, shape, 1, N);

    case 'custom'
        % Custom distribution: user-provided array D
        if isfield(distribution, 'D') && length(distribution.D) == N
            D = distribution.D;
        else
            error('Custom distribution requires a "D" field with N diameters.');
        end

    otherwise
        error('Unsupported distribution type: %s', dist_type);
end

end

function [D_scaled, scale_factor] = scale_diametersND(D, phi, L)
%SCALE_DIAMETERSND Scales particle diameters to achieve the desired packing fraction in ND.
%
% INPUTS:
%   D   : [N x 1] Array of particle diameters.
%   phi : Desired packing fraction (0 < phi < 1).
%   L   : [1 x Ndim] Array of box dimensions in each dimension.
%
% OUTPUT:
%   D_scaled : [N x 1] Array of scaled particle diameters.
%
% NOTES:
% - The packing fraction \( \phi \) is defined as the ratio of the total particle volume
%   to the total box volume.

% Compute the dimensionality (Ndim)
Ndim = length(L);

% Compute the current total volume of particles (ND hypersphere volume)
current_volume = sum((pi^(Ndim/2) / gamma(Ndim/2 + 1)) * (D / 2).^Ndim);

% Compute the total box volume
box_volume = prod(L);

% Compute the desired total particle volume based on the packing fraction
desired_volume = phi * box_volume;

% Compute the scaling factor
scale_factor = (desired_volume / current_volume)^(1 / Ndim);

% Scale the diameters
D_scaled = D * scale_factor;

end