

Ndim = 2;
x0 = rand(500, 2);
D0 = rand(500, 1)+0.1;
Box = ones(1,Ndim);
verbose = 0;

[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, 0);

plot_particles_periodic(x, D, Box)


Ndim = 2;
Box = ones(1,Ndim);
verbose = 0;
N = 500;

phi = 0.1;
exponent = -3;
distribution.type = 'PowerLaw';
distribution.d_min = 1; % Smaller diameter
distribution.d_max = 3; % Larger diameter
distribution.exponent = exponent;  % exponent of power law

[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

plot_particles_periodic(x0, D0, Box)

[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, 0);

plot_particles_periodic(x, D, Box)


Ndim = 2;
Box = ones(1,Ndim);
verbose = 0;
N = 500;

phi = 0.1;
exponent = -1;
distribution.type = 'BiGaussian';
distribution.mu1 = 1; % Smaller diameter
distribution.sigma1 = 0.2; % Smaller std diameter
distribution.mu2 = 5; % Larger diameter
distribution.sigma2 = 1; % Larger std diameter
distribution.p = 0.6;  % 60% of particles with d1

[x0, D0] = initialize_particlesND(phi, N, Box, distribution);

plot_particles_periodic(x0, D0, Box)

[x, D, U_history, phi_history, Fx] = CreatePacking(x0, D0, Box, 0);

plot_particles_periodic(x, D, Box)


