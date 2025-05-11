function plot_particles_3D(x, D, Box, palette_choice, sphere_resolution)
% plot_particles_3D  Render 3D spheres at specified positions with smooth, flat-ish shading
%   x: N-by-3 matrix of sphere centers [x, y, z]
%   D: N-by-1 vector of sphere diameters
%   Box: 1-by-3 vector [Lx, Ly, Lz] specifying axis limits
%   palette_choice: integer selecting color palette (1-12). Optional, default = 1
%   sphere_resolution: resolution of sphere mesh. Optional, default = 50

if nargin < 4
    palette_choice = 1;
end
if nargin < 5
    sphere_resolution = 50;
end

% Unpack Box for axis limits
Lx = Box(1);
Ly = Box(2);
Lz = Box(3);

% Define color palettes (normalized)
palette1 = [204, 85, 85; 255, 179, 71; 255, 255, 128; 119, 221, 119; 127, 106, 217; 0, 225, 225] / 255;
palette2 = [0, 128, 128; 236, 224, 200; 112, 128, 144; 240, 128, 128] / 255;
palette3 = [170, 189, 145; 222, 210, 186; 198, 134, 103; 176, 182, 186] / 255;
palette4 = [143, 170, 180; 244, 221, 220; 196, 177, 162; 60, 76, 89] / 255;
palette5 = [255, 99, 71; 255, 165, 0; 255, 223, 0; 192, 241, 192; 255, 250, 205; 50, 205, 50] / 255;
palette6 = [0, 191, 255; 30, 124, 245; 255, 105, 180; 235, 215, 230; 220, 235, 235; 255, 69, 0] / 255;
palette7 = [255, 110, 255; 0, 255, 127; 173, 255, 47; 255, 235, 250; 155, 80, 210] / 255;
palette8 = [255, 20, 147; 0, 255, 255; 238, 130, 238; 124, 252, 0] / 255;
palette9 = [176, 224, 230; 173, 216, 230; 255, 250, 205; 234, 154, 86; 240, 255, 255; 70, 130, 180] / 255;
palette10 = [255, 162, 173; 142, 251, 142; 212, 251, 212; 173, 216, 230; 255, 235, 250; 90, 150, 200] / 255;
palette11 = [255, 165, 0; 255, 99, 71; 186, 184, 108; 255, 250, 205; 135, 206, 235] / 255;
palette12 = [205, 92, 92; 213, 181, 110; 156, 111, 68; 210, 105, 30; 234, 154, 86; 245, 218, 171] / 255;

switch palette_choice
    case 1, colors = palette1;
    case 2, colors = palette2;
    case 3, colors = palette3;
    case 4, colors = palette4;
    case 5, colors = palette5;
    case 6, colors = palette6;
    case 7, colors = palette7;
    case 8, colors = palette8;
    case 9, colors = palette9;
    case 10, colors = palette10;
    case 11, colors = palette11;
    case 12, colors = palette12;
    otherwise, colors = palette1;
end

% Set up figure
%figure;
clf
hold on;
axis equal;
xlim([0, Lx]);
ylim([0, Ly]);
zlim([0, Lz]);
view(3);

% Lighting setup: softer highlights
lightangle(-45, 30);
lighting phong;
material dull;

% Generate unit sphere mesh
[sx, sy, sz] = sphere(sphere_resolution);

% Number of particles and random color assignment
N = size(x, 1);
rng(42);
particle_colors = colors(randi(size(colors,1), [N,1]), :);

% Plot each sphere with reduced specular and exponent
for i = 1:N
    Xi = x(i,1);
    Yi = x(i,2);
    Zi = x(i,3);
    Ri = D(i)/2;
    h = surf(Ri*sx + Xi, Ri*sy + Yi, Ri*sz + Zi, ...
        'FaceColor', particle_colors(i,:), 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', 'AmbientStrength', 0.07, ...
        'DiffuseStrength', 0.8, 'SpecularStrength', 0.007, ...
        'SpecularExponent', 1);
    set(h, 'LineStyle', 'none');
end

% Final rendering adjustments
camlight headlight;
axis off;
hold off;
