% Molecular Clutches Model for 2D Cases with Viscoelasticity
% Simulates clutch dynamics in 2D using Monte Carlo method with fixed time steps
% and stochastic updates, inspired by the Gillespie algorithm.

clc
clear all

%% Simulation Parameters
% Physical and simulation constants
params.ks1 = 2;          % Short-term stiffness (pN/nm)
params.ks2 = 0.1;        % Long-term stiffness (pN/nm)
params.eta = 1;        % Viscosity parameter
params.kc = 5;           % Clutch stiffness (pN/nm)
params.Fb = 2.5;         % Characteristic force for unbinding (pN)
params.Fm = 2;           % Motor force (pN)
params.vu = 120;         % Unloaded motor velocity (nm/s)
params.vp = 120;         % Protrusion velocity (nm/s)
params.ron = 1;          % Binding rate (s^-1)
params.roff = 0.09;      % Baseline unbinding rate (s^-1)
params.tau = 1;          % Timescale for unbinding adjustment (s)
params.a = 1.5;          % Exponent for unbinding adjustment
params.thres = 0.85;     % Threshold for fraction of bound clutches
params.kcy = 0.001;       % Cytoskeletal stiffness (pN/nm) 

% System sizes
params.Nc = 200;         % Number of clutches
params.Nm = 200;         % Number of motors

% Time parameters
params.Tfinal = 600;     % Total simulation time (s)
params.dt = 0.01;        % Time step (s)
num_steps = ceil(params.Tfinal / params.dt);  % Total number of steps (60000)

%% Initialize System State
% Define variables for each side: left, right, top, bottom
sides = {'left', 'right', 'top', 'bottom'};
for side = sides
    s = side{1};
    % Clutch-specific variables
    eval([s '.Pbi = zeros(params.Nc, 1);']);        % Binding probability
    eval([s '.ronin = params.ron * ones(params.Nc, 1);']);  % Effective on-rate
    eval([s '.Fci = zeros(params.Nc, 1);']);        % Clutch forces (pN)
    eval([s '.xci = zeros(params.Nc, 1);']);        % Clutch positions (nm)
    eval([s '.Ncbound = 0;']);                      % Number of bound clutches
    eval([s '.Frec = 0;']);                         % Accumulated force record
    % Time series arrays
    eval([s '.xs2 = zeros(num_steps, 1);']);        % Long-term substrate position
    eval([s '.xst = zeros(num_steps, 1);']);        % Total substrate position
    eval([s '.Fmy = zeros(num_steps, 1);']);        % Myosin force
    eval([s '.vf = zeros(num_steps, 1);']);         % Retrograde flow velocity
    eval([s '.Pbt = zeros(num_steps, 1);']);        % Fraction of bound clutches
end

% Initialize positions and migration distances
system.D_left = zeros(num_steps + 1, 1);
system.D_right = zeros(num_steps + 1, 1);
system.D_bottom = zeros(num_steps + 1, 1);
system.D_top = zeros(num_steps + 1, 1);
system.Ds_h = zeros(num_steps + 1, 1);  % Horizontal migration distance
system.Ds_v = zeros(num_steps + 1, 1);  % Vertical migration distance

% Set initial positions
system.D_left(1) = -3000;
system.D_right(1) = 3000;
system.D_bottom(1) = -3000;
system.D_top(1) = 3000;
system.Ds_h(1) = 0;
system.Ds_v(1) = 0;

% Initialize coupling forces
Fc_left = 0;
Fc_right = 0;
Fc_top = 0;
Fc_bottom = 0;

% Initialize unbinding rate
roff = params.roff;

%% Main Simulation Loop
for j = 1:num_steps
    % Update bonds for all sides
    [left.Pbi, left.Ncbound, left.Pbt(j), left.roffi] = update_bonds(...
        left.Pbi, left.ronin, roff, left.Fci, params.Fb, params.dt, params.Nc);
    [right.Pbi, right.Ncbound, right.Pbt(j), right.roffi] = update_bonds(...
        right.Pbi, right.ronin, roff, right.Fci, params.Fb, params.dt, params.Nc);
    [top.Pbi, top.Ncbound, top.Pbt(j), top.roffi] = update_bonds(...
        top.Pbi, top.ronin, roff, top.Fci, params.Fb, params.dt, params.Nc);
    [bottom.Pbi, bottom.Ncbound, bottom.Pbt(j), bottom.roffi] = update_bonds(...
        bottom.Pbi, bottom.ronin, roff, bottom.Fci, params.Fb, params.dt, params.Nc);

    % Solve constitutive equations
    [left.xst(j), left.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        left.xci, left.Pbi, left.xs2, left.xst, left.Ncbound);
    [right.xst(j), right.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        right.xci, right.Pbi, right.xs2, right.xst, right.Ncbound);
    [top.xst(j), top.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        top.xci, top.Pbi, top.xs2, top.xst, top.Ncbound);
    [bottom.xst(j), bottom.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        bottom.xci, bottom.Pbi, bottom.xs2, bottom.xst, bottom.Ncbound);

    % Compute previous myosin forces
    Fmy_left_old = sum(params.kc * (left.xci - left.xst(j)) .* left.Pbi);
    Fmy_right_old = sum(params.kc * (right.xci - right.xst(j)) .* right.Pbi);
    Fmy_top_old = sum(params.kc * (top.xci - top.xst(j)) .* top.Pbi);
    Fmy_bottom_old = sum(params.kc * (bottom.xci - bottom.xst(j)) .* bottom.Pbi);

    % Adjust unbinding rate based on conditions
    if any([left.Ncbound, right.Ncbound, top.Ncbound, bottom.Ncbound] == 0) || ...
       any([left.Pbt(j), right.Pbt(j), top.Pbt(j), bottom.Pbt(j)] > params.thres)
        ru = rand;
        toff = params.tau * (1 - ru)^(1 / (1 - params.a));
        roff = 1 / toff;
    end

    % Compute retrograde flow for horizontal (left-right)
    [dxc_left, dxc_right, left.vf(j), right.vf(j)] = compute_retrograde_flow(...
        Fmy_left_old, Fmy_right_old, left.Ncbound, right.Ncbound, Fc_left, params);
    % Compute retrograde flow for vertical (bottom-top)
    [dxc_bottom, dxc_top, bottom.vf(j), top.vf(j)] = compute_retrograde_flow(...
        Fmy_bottom_old, Fmy_top_old, bottom.Ncbound, top.Ncbound, Fc_bottom, params);

    % Update clutch positions and forces
    left.xci = (1 - left.Pbi) .* left.xst(j) + left.Pbi .* (dxc_left + left.xci);
    right.xci = (1 - right.Pbi) .* right.xst(j) + right.Pbi .* (dxc_right + right.xci);
    top.xci = (1 - top.Pbi) .* top.xst(j) + top.Pbi .* (dxc_top + top.xci);
    bottom.xci = (1 - bottom.Pbi) .* bottom.xst(j) + bottom.Pbi .* (dxc_bottom + bottom.xci);

    left.Fci = params.kc * (left.xci - left.xst(j)) .* left.Pbi;
    right.Fci = params.kc * (right.xci - right.xst(j)) .* right.Pbi;
    top.Fci = params.kc * (top.xci - top.xst(j)) .* top.Pbi;
    bottom.Fci = params.kc * (bottom.xci - bottom.xst(j)) .* bottom.Pbi;

    left.Fmy(j) = Fmy_left_old + dxc_left * left.Ncbound * params.kc + Fc_left;
    right.Fmy(j) = Fmy_right_old + dxc_right * right.Ncbound * params.kc + Fc_right;
    top.Fmy(j) = Fmy_top_old + dxc_top * top.Ncbound * params.kc + Fc_top;
    bottom.Fmy(j) = Fmy_bottom_old + dxc_bottom * bottom.Ncbound * params.kc + Fc_bottom;

    % Update positions for the next step
    ra = [1, 1];  % Fixed asymmetry factor
    system.D_left(j+1) = system.D_left(j) - params.dt * (params.vp * ra(1) - left.vf(j));
    system.D_right(j+1) = system.D_right(j) + params.dt * (params.vp * ra(2) - right.vf(j));
    system.D_bottom(j+1) = system.D_bottom(j) - params.dt * (params.vp * ra(1) - bottom.vf(j));
    system.D_top(j+1) = system.D_top(j) + params.dt * (params.vp * ra(2) - top.vf(j));

    % Update migration distances
    system.Ds_h(j+1) = system.Ds_h(j) + params.dt * (left.vf(j) - right.vf(j)) / 2;
    system.Ds_v(j+1) = system.Ds_v(j) + params.dt * (bottom.vf(j) - top.vf(j)) / 2;

    % Compute coupling forces for the next iteration
    delta = (system.D_right(j+1) - system.D_left(j+1)) - (system.D_top(j+1) - system.D_bottom(j+1));
    Fc_left = params.kcy * delta;
    Fc_right = -params.kcy * delta;
    Fc_top = -params.kcy * delta;
    Fc_bottom = params.kcy * delta;
end

%% Post-Processing and Plotting
tt = params.dt * (0:num_steps-1);  % Time vector for plotting

% **Horizontal Components (Left-Right)**
figure('Position', [100, 100, 1200, 800]);
subplot(2, 3, 1)
plot(tt, left.Pbt, 'LineWidth', 2)
hold on
plot(tt, right.Pbt, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fraction probability P_b')
legend('Left', 'Right')
title(sprintf('Horizontal: η = %.0f', params.eta))

subplot(2, 3, 2)
plot(tt, (left.vf - right.vf)/2, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Migration speed (nm/s)')
title('Horizontal Migration Speed')

subplot(2, 3, 4)
plot(tt, sum(left.Fci)/max(left.Ncbound, 0.001), 'LineWidth', 2)
hold on
plot(tt, sum(right.Fci)/max(right.Ncbound, 0.001), 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Average clutch force (pN)')
legend('Left', 'Right')
title('Horizontal Clutch Forces')

subplot(2, 3, 3)
plot(tt, system.Ds_h(1:end-1), 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Migration distance (nm)')
title('Horizontal Migration Distance')

subplot(2, 3, 6)
plot(tt, left.xst, 'LineWidth', 2)
hold on
plot(tt, right.xst, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Substrate displacement (nm)')
legend('Left', 'Right')
title('Horizontal Substrate Displacement')

% **Vertical Components (Bottom-Top)**
figure('Position', [100, 100, 1200, 800]);
subplot(2, 3, 1)
plot(tt, bottom.Pbt, 'LineWidth', 2)
hold on
plot(tt, top.Pbt, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fraction probability P_b')
legend('Bottom', 'Top')
title(sprintf('Vertical: η = %.0f', params.eta))

subplot(2, 3, 2)
plot(tt, (bottom.vf - top.vf)/2, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Migration speed (nm/s)')
title('Vertical Migration Speed')

subplot(2, 3, 4)
plot(tt, sum(bottom.Fci)/max(bottom.Ncbound, 0.001), 'LineWidth', 2)
hold on
plot(tt, sum(top.Fci)/max(top.Ncbound, 0.001), 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Average clutch force (pN)')
legend('Bottom', 'Top')
title('Vertical Clutch Forces')

subplot(2, 3, 3)
plot(tt, system.Ds_v(1:end-1), 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Migration distance (nm)')
title('Vertical Migration Distance')

subplot(2, 3, 6)
plot(tt, bottom.xst, 'LineWidth', 2)
hold on
plot(tt, top.xst, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Substrate displacement (nm)')
legend('Bottom', 'Top')
title('Vertical Substrate Displacement')

%% Helper Functions
% Update Bonds
function [Pbi, Ncbound, Pbt, roffi] = update_bonds(Pbi, ronin, roff, Fci, Fb, dt, Nc)
    roffi = roff * exp(abs(Fci) ./ Fb);  % Force-dependent unbinding rate
    Pbi = Pbi + dt * ((1 - Pbi) .* ronin - Pbi .* roffi);  % Euler step
    Pbi = double(Pbi > rand(Nc, 1));     % Stochastic update
    Ncbound = sum(Pbi);                  % Count bound clutches
    Pbt = Ncbound / Nc;                  % Fraction bound
end

% Motor Clutch Step
function [xst_j, xs2_j] = motor_clutch_step(dt, j, ks1, ks2, eta, kc, xci, Pbi, xs2, xst, Ncbound)
    a1 = ks1 / (ks1 + ks2 + kc * Ncbound);
    a2 = eta / dt;
    xcsum = sum(xci .* Pbi);
    if j == 1
        xs2_j = a1 * kc * xcsum / (ks1 - a1 * ks1 + a2);
    else
        xs2_j = (a2 * xs2(j-1) + a1 * kc * xcsum) / (ks1 - a1 * ks1 + a2);
    end
    if Ncbound == 0 && j > 10
        xs2_j = 0;
    end
    xst_j = (ks1 * xs2_j + kc * xcsum) / (ks2 + ks1 + kc * Ncbound);
end

% Compute Retrograde Flow
function [dxc_A, dxc_B, vf_A, vf_B] = compute_retrograde_flow(Fmy_A_old, Fmy_B_old, Ncbound_A, Ncbound_B, Fc_A, params)
    dFmy = Fmy_A_old - Fmy_B_old;
    Kmy = params.Fm * params.Nm / params.dt / params.vu;
    Ncbound_A_safe = max(Ncbound_A, 0.001);  % Avoid division by zero
    Ncbound_B_safe = max(Ncbound_B, 0.001);
    dxc_A = (Fc_A + params.Fm * params.Nm - Fmy_A_old - Kmy * dFmy / (params.kc * Ncbound_B_safe)) / ...
            (Kmy * (1 + Ncbound_A_safe / Ncbound_B_safe) + Ncbound_A_safe * params.kc);
    dxc_B = (dFmy + params.kc * Ncbound_A_safe * dxc_A) / (params.kc * Ncbound_B_safe);
    vf_A = dxc_A / params.dt;
    vf_B = dxc_B / params.dt;
end