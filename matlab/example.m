
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% List Available Particle Distribution Types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialize_particlesND();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 2D packing with periodic boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndim = 2; % Number of dimensions
Box = ones(1,Ndim); % Set container size to 1 x 1
verbose = 0; % Print packing update as it goes
N = 500; % Number of particles

phi = 0.1; % Initial Packing Fraction for seeding particles

% Particle size distribution
% See function initialize_particlesND for other distribution types
exponent = -3;
distribution.type = 'PowerLaw';
distribution.d_min = 1; % Smaller diameter
distribution.d_max = 3; % Larger diameter
distribution.exponent = exponent;  % exponent of power law

walls = [0,0]; % periodic boundaries (0 for periodic)
fix_height = 0; % do not constrain container height to particle diameter

% Initialize Particles
[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

plot_particles_periodic(x0, D0, Box)

% Create Packing
[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, walls, fix_height, verbose);

% Plot packing
plot_particles_periodic(x, D, Box)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 2D packing with 1 hard boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndim = 2; % Number of dimensions
Box = ones(1,Ndim); % Set container size to 1 x 1
verbose = 0; % Print packing update as it goes
N = 500; % Number of particles

phi = 0.1; % Initial Packing Fraction for seeding particles

% Particle size distribution
distribution.type = 'BiGaussian';
distribution.mu1 = 1; % Smaller diameter
distribution.sigma1 = 0.2; % Smaller std diameter
distribution.mu2 = 5; % Larger diameter
distribution.sigma2 = 1; % Larger std diameter
distribution.p = 0.6;  % 60% of particles with d1

walls = [0,1]; % periodic boundary for x and hard for y (1 for hard wall)
fix_height = 0; % do not constrain container height to particle diameter

% Initialize Particles
[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

plot_particles_periodic(x0, D0, Box)

% Create Packing
[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, walls, fix_height, verbose);

% Plot packing
plot_particles_periodic(x, D, Box)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 2D packing with 1 circular boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndim = 2; % Number of dimensions
Box = ones(1,Ndim); % Set container size to 1 x 1
verbose = 0; % Print packing update as it goes
N = 500; % Number of particles

phi = 0.1; % Initial Packing Fraction for seeding particles

% Particle size distribution
distribution.type = 'BiGaussian';
distribution.mu1 = 1; % Smaller diameter
distribution.sigma1 = 0.2; % Smaller std diameter
distribution.mu2 = 5; % Larger diameter
distribution.sigma2 = 1; % Larger std diameter
distribution.p = 0.6;  % 60% of particles with d1

walls = [-2,1]; % negative integer -N means first N dimensions are circular with diameter of Box(1). So [-2,0,1] in 3D would be a cylinder because first two dimensions are circular while last is straight
fix_height = 0; % do not constrain container height to particle diameter

% Initialize Particles
[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

plot_particles_periodic(x0, D0, Box)

% Create Packing
[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, walls, fix_height, verbose);

% Plot packing
plot_particles_periodic(x, D, Box)
% Draw boundary
hold on
theta = linspace(0,2*pi,1000);
plot(cos(theta)/2+0.5, sin(theta)/2+0.5, '-k', 'LineWidth',2)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 3D packing within sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndim = 3; % Number of dimensions
Box = ones(1,Ndim); % Set container size to 1 x 1
verbose = 0; % Print packing update as it goes
N = 500; % Number of particles

phi = 0.1; % Initial Packing Fraction for seeding particles

% Particle size distribution
distribution.type = 'BiGaussian';
distribution.mu1 = 1; % Smaller diameter
distribution.sigma1 = 0.2; % Smaller std diameter
distribution.mu2 = 5; % Larger diameter
distribution.sigma2 = 1; % Larger std diameter
distribution.p = 0.6;  % 60% of particles with d1

walls = [-3,1,1]; % Make first 3 dimensions circular boundaries
fix_height = 0; % do not constrain container height to particle diameter

% Initialize Particles
[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

% Create Packing
[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, walls, fix_height, verbose);

% Plot packing
plot_particles_3D(x, D, [1,1,1],5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating 2D packing where final container height is fixed multiple of
%%% fist particle diameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndim = 2; % Number of dimensions
Box = [1,0.2]; % Set container size to 1 x 1
verbose = 0; % Print packing update as it goes
N = 500; % Number of particles

phi = 0.1; % Initial Packing Fraction for seeding particles

% Particle size distribution
exponent = -3;
distribution.type = 'PowerLaw';
distribution.d_min = 1; % Smaller diameter
distribution.d_max = 3; % Larger diameter
distribution.exponent = exponent;  % exponent of power law

walls = [0,1]; % negative integer -N means first N dimensions are circular with diameter of Box(1). So [-2,0,1] in 3D would be a cylinder because first two dimensions are circular while last is straight
fix_height = 1; % Constrain container height to first particle diameter as packing generates
ContainerHeight = 6; % Set container height to 6 times the first particle's diameter

% Initialize Particles
[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

%Resize container and particle positions
ratio = (Box(2)/D0(1))/ContainerHeight;
x0(:,2) = x0(:,2)/ratio; % Make sure particles are in the container still
Box(2) = Box(2)/ratio;  % Whatever the ratio of Box(end)/D0(1) that is supplied to CreatePacking is maintained during packing generation when fix_height=1

plot_particles_periodic(x0, D0, [Box(1), Box(2)])
axis([0.3, 0.6, -0.001, Box(2)])

% Create Packing
[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, walls, fix_height, verbose);
Box(2) = ContainerHeight*D(1);

% Plot packing
plot_particles_periodic(x, D, Box)
axis([0.3, 0.6, -0.001, Box(2)+0.001])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
