function plot_particles_periodic(x, D, Box, palette_choice)

if nargin < 4
    palette_choice = 1;
end

Lx = Box(1);
Ly = Box(2);

y = x(:,2);
x = x(:,1);

% Define colors
colors = [
    204, 85, 85;     % Dull apple red
    255, 179, 71;    % 70's pastel orange
    255, 255, 128;   % 70's pastel yellow
    119, 221, 119;   % 70's pastel green
    177, 156 - 50, 217;  % Deeper purple
    0, 225, 225      % Faded cyan
] / 255; % Normalize to [0, 1] for MATLAB

% Define colors
palette1 = [
    204, 85, 85;     % Dull apple red
    255, 179, 71;    % 70's pastel orange
    255, 255, 128;   % 70's pastel yellow
    119, 221, 119;   % 70's pastel green
    177, 156 - 50, 217;  % Deeper purple
    0, 225, 225      % Faded cyan
] / 255; % Normalize to [0, 1] for MATLAB

% % Define color palettes
% palette1 = [
%     204, 85, 85;     % Dull apple red
%     255, 179, 71;    % 70's pastel orange
%     255, 255, 128;   % 70's pastel yellow
%     119, 221, 119];  % 70's pastel green

palette2 = [
    0, 128, 128;    % Teal
    236, 224, 200;  % Pale Sand
    112, 128, 144;  % Slate Blue
    240, 128, 128]/255; % Muted Coral

palette3 = [
    170, 189, 145;  % Soft Olive Green
    222, 210, 186;  % Warm Beige
    198, 134, 103;  % Muted Terracotta
    176, 182, 186]/255; % Cool Gray

palette4 = [
    143, 170, 180;  % Soft Blue-Gray
    244, 221, 220;  % Light Blush
    196, 177, 162;  % Muted Taupe
    60, 76, 89]/255;    % Cool Navy

% Bright and Vibrant Color Palettes (Normalized by 255)

palette5 = [ 
    255, 99, 71;     % Tomato Red
    255, 165, 0;     % Vibrant Orange
    255, 223, 0;     % Bright Yellow
    192, 241, 192;    % Pale Green
    255, 250, 205;   % Lemon Chiffon
    50, 205, 50] / 255;    % Lime Green

palette6 = [ 
    0, 191, 255;     % Deep Sky Blue
    30, 124, 245;    % Dodger Blue
    255, 105, 180;   % Hot Pink
    235, 215, 230;    % Lavender Blush    
    220, 235, 235;    % Azure White
    255, 69, 0] / 255;     % Red-Orange

palette7 = [
    255, 110, 255;     % Magenta
    0, 255, 127;     % Spring Green
    173, 255, 47;    % Bright Green-Yellow
    255, 235, 250;    % Lavender Blush    
    155, 80, 210] / 255;     % Indigo

palette8 = [
    255, 20, 147;    % Deep Pink
    0, 255, 255;     % Cyan
    238, 130, 238;   % Violet
    124, 252, 0] / 255;    % Lawn Green

palette9 = [  % Winter - Cool and Crisp
    176, 224, 230;    % Powder Blue
    173, 216, 230;    % Light Blue
    255, 250, 205;   % Lemon Chiffon
    234, 154, 86;    % Sandy Brown
    240, 255, 255;    % Azure White
    70, 130, 180] / 255;   % Steel Blue

palette10 = [  % Spring - Fresh and Blooming
    255, 162, 173;    % Light Pink
    142, 251, 142;    % Pale Green
    212, 251, 212;    % Pale Green
    173, 216, 230;    % Light Blue
    255, 235, 250;    % Lavender Blush
    90, 150, 200] / 255;   % Light Green

palette11 = [  % Summer - Bright and Playful
    255, 165, 0;     % Orange
    255, 99, 71;     % Tomato Red
    186, 184, 108;   % Olive
    255, 250, 205;   % Lemon Chiffon
    135, 206, 235] / 255;  % Sky Blue

palette12 = [  % Fall - Warm and Earthy
    205, 92, 92;     % Indian Red
    213, 181, 110;     % Gold
    156, 111, 68;     % Coffee
    210, 105, 30;    % Chocolate
    234, 154, 86;    % Sandy Brown
    245, 218, 171] / 255;  % Moccasin


switch palette_choice
    case 1
        colors = palette1;
    case 2
        colors = palette2;
    case 3
        colors = palette3;
    case 4
        colors = palette4;
    case 5
        colors = palette5;
    case 6
        colors = palette6;
    case 7
        colors = palette7;
    case 8
        colors = palette8;        
    case 9
        colors = palette9;
    case 10
        colors = palette10;
    case 11
        colors = palette11;
    case 12
        colors = palette12;
    otherwise
        colors = palette1;
end
rng(42); 

% Number of particles
N = length(x);

% Randomly assign colors to particles
particle_colors = colors(randi(size(colors, 1), [1, N]), :);

%fig = figure('Units', 'inches', 'Position', [0, 0, 20, 10]); % Wider and taller
% Create an invisible figure
% filename = 'particles.svg';
% fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0, 0, 20, 10]);

% Create figure
cla
hold on;
% Offsets for periodic boundaries
offsets = [-1, 0, 1]; % To create periodic replicas in x and y directions

% Draw particles
for i = 1:N
    for dx = offsets
        for dy = offsets
            % Compute periodic positions
            x_shifted = x(i) + dx * Lx;
            y_shifted = y(i) + dy * Ly;

            if x_shifted > -D(i) && x_shifted < Box(1) + D(i) && y_shifted > -D(i) && y_shifted < Box(2) + D(i)
                % Draw particle as a filled circle with black boundary
                rectangle('Position', [x_shifted - D(i)/2, y_shifted - D(i)/2, D(i), D(i)], ...
                          'Curvature', [1, 1], ... % Make it a circle
                          'FaceColor', particle_colors(i, :), ...
                          'EdgeColor', 'k', ...
                          'LineWidth', 1); % Black boundary width
            end
        end
    end
end

% Set axis limits
axis equal;
xlim([0, Lx]);
ylim([0, Ly]);
%title('Particles with Periodic Boundaries');
%xlabel('x');
%ylabel('y');
axis off
hold off;

%saveas(fig, 'particles_plot.svg');
% Save to vector format
% print(fig, filename, '-dsvg', '-painters'); % SVG format

% Close figure to save memory
%close(fig);
