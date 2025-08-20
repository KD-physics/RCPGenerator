function [x, D, U_history, phi_history, F_history] = CreatePacking2(x0, D, Box, walls, fix_height, verbose)
% Optimizes particle packing using the Adam algorithm in ND.
%
% INPUTS:
%   x0       : [N x Ndim] matrix of initial particle positions.
%   D        : [N x 1] array of particle diameters.
%   Box      : [1 x Ndim] array of periodic box dimensions in each dimension.
%   verbose  : Boolean flag for progress display (default: false).
%
% OUTPUTS:
%   x         : [N x Ndim] final particle positions.
%   D         : [N x 1] final particle diameters (scaled to target phi).
%   U_history : Array of energy values over time.
%   phi_history: Array of packing fractions over time.
%   F_history : Array of mean forces over time.
%
% NOTES:
% - Utilizes an adaptive adjustment to the packing fraction (phi).
% - The optimization terminates when phi changes by a negligible amount.
% 
% Last Update
%   - December 14, 2024
%   - Kenneth Desmond
%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% --- Configure System
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MaxNeighbors = [300, 750, 5500, 5500];
    
    % phi0 = 0.025; 
    % D = D/sqrt(sum(pi*(D/2).^2)/phi0);
    
    % Extract box dimensions and particle counts
    N = size(x0, 1); % Number of particles
    Ndim = size(x0, 2); % Number of dimensions
    L = Box; % Box dimensions in each direction
    L0 = L;

    % Default verbosity
    if nargin < 4
        walls = zeros(1,Ndim);
        fix_height = 0;
        verbose = false;
    end
    if nargin < 5
        fix_height = 0;
        verbose = false;
    end
    if nargin < 6
        verbose = false;
    end
    
    %%%%%%%%%%%%%%%%%%
    %%% Settings
    %%%%%%%%%%%%%%%%%%
    
    % Fix algorithm coefficients and settings
    alpha_max = 0.0025;
    alpha = alpha_max; % Adam learning rate
    beta1 = 0.9; % Exponential decay rate for first moment
    beta2 = 0.999; % Exponential decay rate for second moment
    epsilon = 1e-8; % Small value to prevent division by zero
    N_steps = 80000; % Max steps before giving up and terminating early
    dt = 0.1; % Time step for verlet
    method = 'ADAM';
    
    %%%%%%%%%%%%%%%%%%
    %%% Initialize Variables
    %%%%%%%%%%%%%%%%%%
    
    % Initialize particle positions and moment estimates
    x = x0; % Initial positions
    F = zeros(N, Ndim); % Initialize force matrix
    m = zeros(N, Ndim); % First moment for positions
    v = zeros(N, Ndim); % Second moment for positions
    v_max = zeros(N, Ndim); % max first moment for positions
    t = 0; % Time step for bias correction
    
    % Initialize history tracking
    U_history = zeros(1, N_steps); % Store energy at each step
    phi_history = zeros(1, N_steps); % Store phi at each step
    F_history = zeros(1, N_steps); % Store mean force at each step
    
    v_verlet = zeros(N, Ndim); 
    a = zeros(N, Ndim);
    a_old = zeros(N, Ndim);
    verlet_drag = 2;
    
    v_update = zeros(N, Ndim); 
    m_hat = zeros(N, Ndim); 
    v_hat = zeros(N, Ndim); 
    
    kappa = 1;
    dkappa = 0;
    m_kappa = 0; 
    v_kappa = 0; 
    m_hat_kappa = 0; 
    v_hat_kappa = 0; 
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Initialize Particle Diameters
    %%%%%%%%%%%%%%%%%%
    phi_modifier = 1;
    if walls(1) < 0
        if walls(1) == -1
            walls(1) = -2;
        end
        for k=2:(-walls(1))
            walls(k) = -1;
        end
        phi_modifier = (pi^(abs(walls(1))/2) / gamma(abs(walls(1))/2 + 1)) * (1 / 2).^abs(walls(1)); % provide correction to container volume now its a hypersphere and not a box
    end

    phi = 0.15; % Initial packing fraction
    if Ndim == 2; phi = 0.76; end % Initial packing fraction
    if Ndim == 3; phi = 0.45; end % Initial packing fraction
    if Ndim == 4; phi = 0.15; end % Initial packing fraction
    if Ndim == 5; phi = 0.1; end % Initial packing fraction
    if Ndim > 5 ; phi = 0.1/sqrt(Ndim); end % Initial packing fraction
    [D, scale] = scale_diametersND(D, phi*phi_modifier, L, fix_height); % Set Initial Diameters
    if fix_height
        L(end) = L(end)*scale;
    end
    D = D(:); % Make sure diameters are in a column vector
    
    %%%%%%%%%%%%%%%%%%
    %%% Set Phi Stepping and Termination Conditions
    %%%%%%%%%%%%%%%%%%
    delta_phi0 = 1.5*1.5E-3; % Initial phi step size
    delta_phi = delta_phi0;
    dphi = 0; % Current step size
    count = 0; % Counter for successful phi adjustments
    direction_flag = 0; % Tracks direction of phi updates
    U_threshold = (2.5E-4)^2; % Energy threshold for minimal overlap
    F_tol = sqrt(U_threshold)/50; % Force threshold for concluding phi too large and needs to be lowered
    
    
    %%%%%%%%%%%%%%%%%%
    %%% Initialize pairs and pairs updating parameters
    %%%%%%%%%%%%%%%%%%
    
    Nneighbors = MaxNeighbors(min([4,Ndim-1]));
    pairs = zeros(N, Nneighbors, 'uint16');
    
    refresh = ones(N,1);
    pairs = GetPairsND_3(N, x, D, L, pairs, refresh); % pair list
    Dmin = min(D); % Track the minimum diameter. Used to determine updating pair list
    x_last = x; % tracking position of particle during last update to trigger new update
    update_flag = false; % flag to trigger update
    Dmin_last = Dmin; % Dmin during last update. If Dmin grows too much, trigger entire pairs rebuild
    LastAlphaUpdate = 1;
    LastPhiUpdate = 0;
    
    phi_min = 0.8;
    if Ndim == 2; phi_min = 0.36; end
    if Ndim == 3; phi_min = 0.61; end
    if Ndim == 4; phi_min = 0.4; end
    if Ndim == 5; phi_min = 0.18; end
        
    R = 1/1.0002;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% --- Main Function ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DD = D;
    N_steps = 150000;
    mu = 5E-4;
    mu_flag = 1;
    alpha = 0.003;
    alpha_max = 0.003;
    mu_change = 1E9;
    phi_max = [0,0];

    % Diameter Growth and Optimization Loop
    for step = 1:N_steps
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Optimizer Given Current Conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % if dphi ~= 0
        %     if verbose; fprintf("Step %d: Adjusting phi %.3e\n",step, dphi); end        
        %     if ~strcmp(method, 'ADAM') && verbose
        %         fprintf('Swithcing to Method: %s\n', 'Adam')
        %     end
        %     method = 'ADAM';
        %     LastPhiUpdate = step;
        %     if dphi > 0
        %         alpha = min(alpha*1.25, alpha_max);
        %     end
        % end
        % 
        % 
        % if step - LastPhiUpdate > 2500 && strcmp(method, 'ADAM')
        %     if verbose
        %         fprintf('Swithcing to Method: %s\n', 'AMSGrad')
        %     end
        %     method = 'AMSGrad';
        % end
        % 
        % if step - LastPhiUpdate > 4000 && strcmp(method, 'AMSGrad')
        %     if verbose
        %         fprintf('Swithcing to Method: %s\n', 'Verlet')
        %     end
        %     method = 'Verlet';
        % 
        %     v_verlet(:) = 0;
        %     a(:) = 0;
        %     a_old(:) = 0;
        % end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Learning Rate Management
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Update timestep for bias correction
        t = t + 1;
    
        % % Decay learning rate
        % alpha = alpha*R;
        % 
        % if dphi > 0
        %     LastAlphaUpdate = step;
        % end
        if step - LastAlphaUpdate > 500 && strcmp(method, 'ADAM') && (mu_flag < 1 || phi - phi_history(max([step-250,1])) < 5E-4) && step > 4000
            trend = mean(F_history(step-450:step-250))/mean(F_history(step-75:step-1));        
            if trend < 0.85
                % if verbose; fprintf('Step %d: Lowering Alpha by 2x\n', step); end
                LastAlphaUpdate = step;
                alpha = alpha/1.5;
            elseif step - LastAlphaUpdate > 1500
                % if verbose; fprintf('Step %d: Raising Alpha by 1.5x\n', step); end
                trend = mean(F_history(step-200:step-1))/mean(F_history(step-3000:step-2500));
                if trend > 0.5
                    LastAlphaUpdate = step;
                    alpha = alpha*1.15;
                end

            end
            if mu_flag == 1
                alpha = max([alpha, alpha_max/10]);
            end
        end

        alpha = min([alpha*R, alpha_max]);        
        if mod(step, 500) == 0
            alpha = min([alpha_max,1.1*alpha]);
        end

    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Updating Phi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Update Particle Diameters
        D = DD*kappa;
        Dmin_latest = min(D);

        % Compute the current total volume of particles (ND hypersphere volume)
        current_volume = sum((pi^(Ndim/2) / gamma(Ndim/2 + 1)) * (D / 2).^Ndim);
        
        % Compute the total box volume
        box_volume = prod(L);
        if fix_height
            L(end) = L0(end)*kappa;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% mu Management
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        phi = current_volume/box_volume;
        phi_history(step) = phi/phi_modifier;


        if phi > phi_max(1)
            phi_max = [phi, step];
        end
        if step > 500 && mu_flag == 1
            XX = phi_history(step-250:step);
            delta_XX = max(XX)-min(XX);
            if delta_XX < 1E-5 || step - phi_max(2) > max([2000, N/3])
                mu = mu/10;
                alpha = alpha/2;
                mu_flag = 0;
                mu_change = step;
            end
        end

        
        if step > 500 && mu_flag == 0 && step > mu_change+250
            XX = phi_history(step-250:step);
            delta_XX = max(XX)-min(XX);
            if delta_XX < 1E-6
                mu = mu/10;
                mu_flag = -1;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Managing Pairs List
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Full rebuild ever 750 steps or Diameters grew too much
        % if mod(step, 750) == 0 || Dmin_last/Dmin > 1.05
        if Dmin_latest/Dmin > 1.05
            % reset refresh to 1
            for k = 1:N
                refresh(k) = 1;
            end
            pairs = GetPairsND_3(N, x, D, L, pairs, refresh); 
            x_last = x;
            Dmin = min(D);
        else    
            % reset refresh to 0
            for k = 1:N
                refresh(k) = 0;
            end
            % For particles that move more than threshold, update their pair list    
            for k = 1:N
                if norm(x(k,:)-x_last(k,:)) > Dmin/4
                    update_flag = true;
                    x_last(k,:) = x(k,:);
                    refresh(k) = 1;
                end
            end
        
            % If update_flag indicate pair update trigger then update
            if update_flag 
                pairs = GetPairsND_3(N, x, D, L, pairs, refresh);
            end
        end
        update_flag = false;
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Forces
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [F, U, min_dist, Lc, z, dkappa, Fmean] = GetForcesND_3(pairs, x, D, F, L, walls, mu);
        dkappa = dkappa/kappa;
        %F_magnitude = mean(sqrt(sum(F.^2,2))sqrt(sum(F.^2,2)))/Lc/mean(z)*Ndim/sqrt(Ndim); % Compute mean force magnitude for update condition
        F_magnitude = sqrt(sum(F.^2,2));        
        F_magnitude = mean(F_magnitude(F_magnitude > 1E-16))/Fmean/sqrt(Ndim);
        F_history(step) = F_magnitude;

        % Store energy
        U_history(step) = U;
        
        % if F_magnitude < 1E-6 && phi > 0.76
        % % if max(min_dist(:))/F_magnitude < 1E-2 && phi > 0.76
        %     mu = mu*0.75;
        %     % fprintf('step = %d, mu = %.2e \n', step, mu)
        % end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Optimizer updates: moments, acceleration, and velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Update Adam moment estimates
        m = beta1 * m - (1 - beta1) * F;
        v = beta2 * v + (1 - beta2) * (F.^2);
        v_max = max(v_max, v);                     % Element-wise max
    
        m_kappa = beta1 * m_kappa - (1 - beta1) * dkappa;
        v_kappa = beta2 * v_kappa + (1 - beta2) * (dkappa.^2);
        

        % AMSGrad: Update maximum second moment estimate
        if strcmp(method, 'AMSGrad')
            v_update = v_max;                            % Use v_max in updates
            v_update_kappa = v_kappa;
        else
            v_update = v;                                % Use standard v in ADAM
            v_update_kappa = v_kappa;
        end
    
        % Compute adam bias-corrected moment estimates
        m_hat = m / (1 - beta1^t);
        v_hat = v_update / (1 - beta2^t);
    
        m_hat_kappa = m_kappa / (1 - beta1^t);
        v_hat_kappa = v_update_kappa / (1 - beta2^t);

        kappa = kappa - alpha * m_hat_kappa ./ (sqrt(v_hat_kappa) + epsilon);

        % Update verlet accelerations
        a = F - verlet_drag * v_verlet;
    
        % Update verlet velocities
        v_verlet = v_verlet + 0.5 * (a_old + a) * dt;
    
        % Store current verlet accelerations
        a_old = a;        
        % [m, v, v_max, v_update, m_hat, v_hat, a, v_verlet, a_old] = AdamStep(method, N, Ndim,  F, beta1, beta2, t, dt, verlet_drag, m, v, v_max, v_update, m_hat, v_hat, a, v_verlet, a_old);     
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Termination and phi updating conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        if mu_flag == -1 && F_history(step) < 5E-3 && max(min_dist(:)) > 1E-16 && std(phi_history(step-1000:step))/mean(phi_history(step-1000:step)) < 1E-6
            if verbose
                fprintf('Success: Packing achieved.\n');
                fprintf('Step %d: Phi = %.6f, mu = %.2e, min_dist = %.2e, |F|/<F> = %.2e, mu = %.2e \n', ...
                        step, phi_history(step), mu, max(min_dist(:)), F_history(step), mu);

            end
            break;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Positions updates
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Update positions
        x_old = x;
        if ~strcmp(method, 'Verlet')
            x = x - alpha * m_hat ./ (sqrt(v_hat) + epsilon);
        else
            x = x + v_verlet * dt + 0.5 * a_old * dt^2;
        end
        delta_x = mean(sqrt(sum((x_old-x).^2,2))./D);
        % Apply periodic boundary conditions
        x = mod(x, L);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        % Optional progress display
        if (mod(step, 500) == 0) && verbose
            fprintf('Step %d: φ = %.6f, max overlap = %.2e, mu = %.2e, |F|/<F> = %.2e, Neighbors/changes = %d / %d, α = %.2e, max d = %.2e, κ = %.3e \n', ...
                    step, phi_history(step), max(min_dist(:)), mu, F_history(step), max(pairs(:,1)), sum(refresh), alpha, delta_x, kappa);
        end
        

    end

    % Trim histories to actual steps
    U_history = U_history(1:step);
    phi_history = phi_history(1:step);
    F_history = F_history(1:step);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- Sub Functions ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m, v, v_max, v_update, m_hat, v_hat, a, v_verlet, a_old] = AdamStep(method, N, Ndim,  F, beta1, beta2, t, dt, verlet_drag, m, v, v_max, v_update, m_hat, v_hat, a, v_verlet, a_old)
% Make an Adam Step updating all parameters 
%
% INPUTS:
%   method    : Update method; Adam, AMSGrad, Verlet
%   N         :
%   Ndim      :
%   F         : [N x Ndim] Array of particle forces.
%   beta1 [2] : Fixed Adam parameters.
%   t         : time stamp
%   dt        : time step for when switching to verlet
%   remaining : Adam and verlet update parameters
%
% OUTPUT:
%   Updated parameters.
%
% NOTES:
% - The packing fraction \( \phi \) is defined as the ratio of the total particle volume
%   to the total box volume.

    for k = 1:N
        for kk = 1:Ndim
            % Update Adam moment estimates
            m(k,kk) = beta1 * m(k,kk) - (1 - beta1) * F(k,kk);
            v(k,kk) = beta2 * v(k,kk) + (1 - beta2) * (F(k,kk)^2);
            v_max(k,kk) = max(v_max(k,kk), v(k,kk));                     % Element-wise max
        
            % AMSGrad: Update maximum second moment estimate
            if strcmp(method, 'AMSGrad')
                v_update(k,kk) = v_max(k,kk);                            % Use v_max in updates
            else
                v_update(k,kk) = v(k,kk);                                % Use standard v in ADAM
            end
        
            % Compute adam bias-corrected moment estimates
            m_hat(k,kk) = m(k,kk) / (1 - beta1^t);
            v_hat(k,kk) = v_update(k,kk) / (1 - beta2^t);
        
            % Update verlet accelerations
            a(k,kk) = F(k,kk) - verlet_drag * v_verlet(k,kk);
        
            % Update verlet velocities
            v_verlet(k,kk) = v_verlet(k,kk) + 0.5 * (a_old(k,kk) + a(k,kk)) * dt;
        
            % Store current verlet accelerations
            a_old(k,kk) = a(k,kk);
        end
    end
end

function [D_scaled, scale_factor] = scale_diametersND(D, phi, L, fix_height)
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
    if fix_height
        scale_factor = (desired_volume / current_volume)^(1 / (Ndim - 1));
    else
        scale_factor = (desired_volume / current_volume)^(1 / Ndim);
    end

    % Scale the diameters
    D_scaled = D * scale_factor;

end

function pairs = GetPairsND_3(N, x, D, L, pairs, refresh)
%GETPAIRSND Identifies pairs of particles within the cutoff distance in ND.
%
% INPUTS:
%   N    : Number of particles.
%   x    : [N x Ndim] matrix of particle positions in Ndim dimensions.
%   D    : [N x 1] vector of particle diameters.
%   L    : [1 x Ndim] array of box dimensions in each dimension.
%   refresh    : Particle indices to update because displacement exceeded threshold
%
% OUTPUT:
%   pairs : [M x 2] array of interacting particle pairs (indices),
%           where M is the number of pairs satisfying the cutoff distance.
%
% NOTES:
% - Periodic boundary conditions are applied to account for box wrapping.
% - Pairs are symmetric (e.g., if (i, j) is in the list, (j, i) is omitted).
% - The function is optimized for performance and avoids unnecessary calculations.
    
    
    Ndim = size(x, 2); % Number of dimensions
    
    Dmax = max(D);
    t = min([max([median(D)*1.05, min(D)*3]), max(D)]);
    [~, sort_idx] = sort(x(:,1));
    
    % Create the inverse mapping
    sort_location = zeros(size(sort_idx));
    sort_location(sort_idx) = 1:length(sort_idx);
    %'starting'
    j_list = zeros(N+1,1);
    for i = 1:N
        %fprintf("%d \n",i)
        if refresh(i) == 1
            
            pairs_old = pairs(i,1:pairs(i,1)+1);
            pairs_old(2:end)=sort(pairs_old(2:end));
        
            pairs(i,1) = 0;
            
            %Reduce to pair search to nearby indices is sorted list
            j_list(1) = 0;
            %Get lower limit
            r_c_max = (D(i)+1.01*Dmax)/2+t;
            for direction = -1:2:1
                j_list_flag = 1;
                jdx = 0;
                while j_list_flag
                    jdx = jdx+direction;
                    if direction*jdx > N/2
                        j_list_flag = 0;
                    end
                    j = sort_idx(mod(sort_location(i)+jdx-1,N)+1);
                    dx = x(j,1) - x(i,1);
                    dx = dx - round(dx / L(1)) * L(1);
                    r_c = (D(i)+D(j))/2+t;
                    if abs(dx) > r_c_max  
                        j_list_flag = 0;
                    elseif j_list_flag && abs(dx) < r_c
                        % Check if within cutoff in x-direction
                        if true %abs(dx) < r_c
                            proceed = true;
                            if proceed
                                for k = 2:Ndim
                                    dz = x(j,k) - x(i,k);
                                    dz = dz - round(dz / L(k)) * L(k);
                                    if abs(dz) > r_c
                                        proceed = false;
                                        break
                                    end
                                end
                            end
            
                            if proceed
                                pairs(i,1) = pairs(i,1) + 1;
                                pairs(i,pairs(i,1)+1) = j;
            
                                j_flag = false;
                                if ~j_flag && pairs_old(1) > 0
                                    %Check for duplicate
                                    not_a_duplicate = true;
                                    for jdx_check = 2:pairs(j,1)+1
                                        if i == pairs(j,jdx_check)
                                            not_a_duplicate = false;
                                        end
                                    end
                                    if not_a_duplicate
                                        pairs(j,1) = pairs(j,1) + 1;
                                        pairs(j,pairs(j,1)+1) = i;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %"finished"
end

function [F, U, min_dist, Lc, z, dkappa, Fmean] = GetForcesND_3(pairs, x, D, F, L, walls, mu)
%GETFORCESND Computes the pairwise forces, energy, and contact count for particles in ND.
%
% INPUTS:
%   pairs : [Npairs x 2] array of interacting particle pairs (indices).
%   x     : [N x Ndim] matrix of particle positions in Ndim dimensions.
%   D     : [N x 1] array of particle diameters.
%   F     : [N x Ndim] matrix for storing the forces in Ndim dimensions.
%   L     : [1 x Ndim] array of box dimensions in each dimension.
%
% OUTPUTS:
%   F     : Updated [N x Ndim] force matrix with computed pairwise forces.
%   U     : Total potential energy of the system.
%   z     : Coordination number or average number of particle contacts (overlapping pairs).
%
% NOTES:
% - Periodic boundary conditions are applied to ensure particles interact
%   across box boundaries.
% - Pairwise forces are computed using a linear spring model:
%     F = -K*(r_ij / d_ij - 1), with K set to 1.
%   r_ij is the sum of particle radii, and d_ij is the distance.
% - Potential energy for each overlapping pair:
%     E_pair = K/2 * (r_ij - d_ij)^2
% - \( K \) is the spring constant.
    
    
    % Parameters
    N = size(pairs, 1); % Number of particles
    Ndim = size(x, 2); % Number of dimensions
    if Ndim == 1
        Ndim = length(x)/N;
    end
    min_dist = zeros(size(pairs));
    
    
    % Initialize outputs
    z = x(:,1)*0; % Contact count
    U = 0; % Total potential energy
    Lc = 0;
    count = 1;
    %min_dist = 0;
    dx = zeros(1,Ndim);
    F(:) = 0;
    Fmean = 0;

    wall_flag = max(walls);
    circle_flag = walls(1) < 0;

    dkappa = 0;
    % Loop over all interacting pairs
    for i = 1:N
        % Particle indices
        for jdx = 2:pairs(i, 1)+1
            j = pairs(i, jdx);
    
            if j > i% && not_a_duplicate
            
                % Sum of particle radii
                r_ij = (D(i) + D(j)) / 2;
                
                flag = 1;
                d_ij = 0;
                for k=1:Ndim
                    if flag == 1
                        dx(k) = x(j, k) - x(i, k);
                        if (walls(k)==0); dx(k) = dx(k) - round(dx(k) ./ L(k)) .* L(k);end
                        if abs(dx(k)) > r_ij
                            flag = 0;
                        end
                        d_ij = d_ij + dx(k)^2;
                    end
                end
        
                if flag == 1                
                    % Compute Euclidean distance between particles
                    d_ij = sqrt(d_ij);
            
                    % Check for overlap (contact)
                    if d_ij < r_ij

                        % Increment contact count
                        z(i) = z(i) + 1;
                        z(j) = z(j) + 1;
            
                        K = 1;
                        % Compute pairwise force magnitude
                        F_mag = -K * (r_ij / d_ij - 1);
            
                        % Compute potential energy contribution
                        % U = U + (1 - d_ij/r_ij);
                        U = U + K*(d_ij - r_ij)^2/2;
                        count = count + 1;
                        min_dist(i,jdx) = (1 - d_ij/r_ij);
                        Lc = Lc + r_ij;
            
                        % Update forces in all dimensions
                        F(i, :) = F(i, :) + F_mag * dx; % Force on particle i
                        F(j, :) = F(j, :) - F_mag * dx; % Reaction force on particle j                    
        
                        dkappa = dkappa + r_ij*(d_ij - r_ij);%2*F_mag*r_ij;

                        Fmean = Fmean + F_mag;

                    end
                end
            end
        end

        if wall_flag
            r_ij = D(i)/2;            
            for M=1:Ndim
                if walls(M)
                    for k=0:1
                        dx(M) = k*L(M) - x(i, M);
                        d_ij = abs(dx(M));
        
                        if d_ij < r_ij
                            K = 2.;
                            % Compute pairwise force magnitude
                            F_mag = -K * (r_ij / d_ij - 1);
        
                            % Compute potential energy contribution
                            U = U + K*(d_ij - r_ij)^2/2;
                            count = count + 1;
        
                            % Update forces in all dimensions
                            F(i, M) = F(i, M) + F_mag * dx(M); % Force on particle i
                            
                            dkappa = dkappa + r_ij*(d_ij - r_ij);
                            Fmean = Fmean + F_mag;
                            
                        end                
                    end
                end
            end
        end

        if circle_flag
            r_ij = D(i)/2;
            d_ij = 0;
            R = L(1)/2;
            for M=1:(-walls(1))
                d_ij = d_ij + (x(i,M)-R)^2;
            end
            d_ij = sqrt(d_ij);
            delta = R - d_ij;
            if delta < r_ij
                delta = d_ij;
                d_ij = 0;
                for M=1:(-walls(1))
                    dx(M) = abs(R/delta - 1)*(x(i, M)-R);
                    d_ij = d_ij + dx(M)^2;
                end
                d_ij = sqrt(d_ij);

                K = 2.;
                % Compute pairwise force magnitude
                F_mag = -K * abs(r_ij / d_ij - 1);

                % Compute potential energy contribution
                U = U + K*(d_ij - r_ij)^2/2;
                count = count + 1;

                % Update forces in all dimensions
                for M=1:(-walls(1))
                    F(i, M) = F(i, M) + F_mag * dx(M); % Force on particle i
                end
                
                dkappa = dkappa + K*r_ij*(d_ij - r_ij);
                Fmean = Fmean + F_mag;

            end
        end

    end
    U = U^2 - mu*sum(D);
    Lc = Lc/count;
    dkappa = dkappa + mu*sum(D);
    Fmean = -Fmean/count;


end




