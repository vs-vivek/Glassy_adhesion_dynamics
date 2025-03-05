% Molecular Clutches Model for 1D Cases with Viscoelasticity
% Simulates clutch dynamics using Monte Carlo method with fixed time steps
% and stochastic updates, inspired by the Gillespie algorithm.
% Refer to readme file for instructions on how to use this code.

clc
clear all

%% Simulation Parameters
% Define physical and simulation constants
params.ks1 = 1;          % Short-term stiffness 
params.ks2 = 0.3;        % Long-term stiffness 
params.eta = 1000;          % Viscosity parameter CHANGE THIS PARAMETER TO CHANGE RELAXATION TIMESCALE
                         % Fast relaxing --> eta = 1, Slow relaxing --> eta = 1000   
params.kc = 5;           % Clutch stiffness 
params.Fb = 2.5;         % Characteristic force for unbinding 
params.Fm = 2;           % Motor force 
params.vu = 120;         % Unloaded motor velocity 
params.vp = 120;         % Protrusion velocity 
params.ron = 1;          % Binding rate 
params.roff = 0.1;       % Baseline unbinding rate 
params.tau = 1;          % Timescale for unbinding adjustment 
params.a = 1.5;          % Exponent for unbinding adjustment
params.thres = 0.86;     % Threshold for fraction of bound clutches
params.Ftalin = 0.2;     % Talin force (unused in dynamics but kept)
params.kcy = 0;

% System sizes
params.Nc = 200;         % Number of clutches
params.Nm = 200;         % Number of motors

% Time parameters
params.Tfinal = 2000;     % Total simulation time
params.dt = 0.01;        % Time step
num_steps = ceil(params.Tfinal / params.dt);  % Total number of steps 

%% Initialize System State
% Left side
left.Pbi = zeros(params.Nc, 1);        % Clutch binding probability
left.ronin = params.ron * ones(params.Nc, 1);  % Effective on-rate
left.Fci = zeros(params.Nc, 1);        % Clutch forces 
left.xci = zeros(params.Nc, 1);        % Clutch positions 
left.Ncbound = 0;                      % Number of bound clutches
left.Frec = 0;                         % Accumulated force record

% Right side
right.Pbi = zeros(params.Nc, 1);
right.ronin = params.ron * ones(params.Nc, 1);
right.Fci = zeros(params.Nc, 1);
right.xci = zeros(params.Nc, 1);
right.Ncbound = 0;
right.Frec = 0;

% Preallocate time series arrays
left.xs2 = zeros(num_steps, 1);    % Long-term substrate position (nm)
left.xst = zeros(num_steps, 1);    % Total substrate position (nm)
left.Fmy = zeros(num_steps, 1);    % Myosin force (pN)
left.vf = zeros(num_steps, 1);     % Retrograde flow velocity (nm/s)
left.Pbt = zeros(num_steps, 1);    % Fraction of bound clutches
right.xs2 = zeros(num_steps, 1);
right.xst = zeros(num_steps, 1);
right.Fmy = zeros(num_steps, 1);
right.vf = zeros(num_steps, 1);
right.Pbt = zeros(num_steps, 1);

% Global variables
system.Vm = zeros(num_steps, 1);   % Migration velocity 
system.Ds = zeros(num_steps + 1, 1);   % Substrate displacement 
system.D_left = zeros(num_steps + 1, 1);  % Left edge position 
system.D_right = zeros(num_steps + 1, 1); % Right edge position 
system.roff_record = zeros(num_steps, 1);  % Adjusted unbinding rate
system.Fci_record_left = zeros(num_steps, 1);  % Total clutch force (left)
system.Fci_record_right = zeros(num_steps, 1); % Total clutch force (right)
system.roff_record_left = zeros(params.Nc, num_steps); % Per-clutch unbinding rates (left)
system.roff_record_right = zeros(params.Nc, num_steps); % Per-clutch unbinding rates (right)
system.Pbidoc_1 = zeros(num_steps, 2);  % First clutch binding state

% Initial conditions
system.Ds(1) = 0;         % Initial substrate displacement
system.D_left(1) = -3000; % Initial left edge 
system.D_right(1) = 3000; % Initial right edge 
roff = params.roff;       % Initial unbinding rate

%% Main Simulation Loop
for j = 1:num_steps
    % Update bond states
    [left.Pbi, left.Ncbound, left.Pbt(j), roffi_left] = update_bonds(...
        left.Pbi, left.ronin, roff, left.Fci, params.Fb, params.dt, params.Nc);
    [right.Pbi, right.Ncbound, right.Pbt(j), roffi_right] = update_bonds(...
        right.Pbi, right.ronin, roff, right.Fci, params.Fb, params.dt, params.Nc);

    % Record additional variables
    system.roff_record_left(:, j) = roffi_left;
    system.roff_record_right(:, j) = roffi_right;
    system.Pbidoc_1(j, :) = [left.Pbi(1), right.Pbi(1)];

    % Accumulate forces
    left.Frec = left.Frec + left.Fmy(j) / (left.Ncbound + 1e-10);
    right.Frec = right.Frec + right.Fmy(j) / (right.Ncbound + 1e-10);

    % Update constitutive equations
    [left.xst(j), left.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        left.xci, left.Pbi, left.xs2, left.xst, left.Ncbound);
    Fmy_left_old = sum(params.kc * (left.xci - left.xst(j)) .* left.Pbi);

    [right.xst(j), right.xs2(j)] = motor_clutch_step(...
        params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
        right.xci, right.Pbi, right.xs2, right.xst, right.Ncbound);
    Fmy_right_old = sum(params.kc * (right.xci - right.xst(j)) .* right.Pbi);

    % Adjust unbinding rate
    if left.Ncbound == 0 || right.Ncbound == 0 || ...
       left.Pbt(j) > params.thres || right.Pbt(j) > params.thres
        ru = rand;
        toff = params.tau * (1 - ru)^(1 / (1 - params.a));
        roff = 1 / toff;
    end
    system.roff_record(j) = roff;

    % Compute retrograde flow
    dFmy = Fmy_left_old - Fmy_right_old;
    Kmy = params.Fm * params.Nm / params.dt / params.vu;
    Ncbound_left_safe = max(left.Ncbound, 0.001);
    Ncbound_right_safe = max(right.Ncbound, 0.001);
    Fcyto = params.kcy * ((system.D_right(j) - system.D_left(j)) - ...
                          (system.D_right(1) - system.D_left(1)));
    dxc_left = (params.Fm * params.Nm - Fmy_left_old - ...
                Kmy * dFmy / (params.kc * Ncbound_right_safe) + Fcyto) / ...
               (Kmy * (1 + Ncbound_left_safe / Ncbound_right_safe) + ...
                Ncbound_left_safe * params.kc);
    dxc_right = (dFmy + params.kc * Ncbound_left_safe * dxc_left) / ...
                (params.kc * Ncbound_right_safe);
    left.vf(j) = dxc_left / params.dt;
    right.vf(j) = dxc_right / params.dt;

    % Update positions and forces
    left.xci = (1 - left.Pbi) .* left.xst(j) + left.Pbi .* (dxc_left + left.xci);
    right.xci = (1 - right.Pbi) .* right.xst(j) + right.Pbi .* (dxc_right + right.xci);
    left.Fci = params.kc * (left.xci - left.xst(j)) .* left.Pbi;
    right.Fci = params.kc * (right.xci - right.xst(j)) .* right.Pbi;
    left.Fmy(j) = Fmy_left_old + dxc_left * left.Ncbound * params.kc;
    right.Fmy(j) = Fmy_right_old + dxc_right * right.Ncbound * params.kc;

    % Record total clutch forces
    system.Fci_record_left(j) = sum(left.Fci);
    system.Fci_record_right(j) = sum(right.Fci);

    % Update cell migration for the next step
    if j < num_steps
        system.Vm(j) = (left.vf(j) - right.vf(j)) / 2;
        ra = [1, 1];  % Fixed asymmetry factor
        system.D_left(j+1) = system.D_left(j) - params.dt * (params.vp * ra(1) - left.vf(j));
        system.D_right(j+1) = system.D_right(j) + params.dt * (params.vp * ra(2) - right.vf(j));
        system.Ds(j+1) = system.Ds(j) + params.dt * system.Vm(j);
    end
end

%% Post-Processing
tt = params.dt * (1:num_steps);  % Time vector for plotting
Vsa = params.vp - mean(left.vf + right.vf) / 2;  % Average speed
MSDfinal = system.Ds(end)^2;  % Final mean squared displacement
Vma = mean(abs(system.Vm(100:end)));  % Mean absolute migration velocity

%% Plotting Results
figure('Position', [100, 100, 1200, 800]);
set(gcf, 'Color', 'w');

% Fraction of bound clutches
subplot(2, 3, 1)
plot(tt, left.Pbt, 'LineWidth', 2);
hold on
plot(tt, right.Pbt, 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Fraction probability P_b')
legend('Pbt_Left', 'Pbt_Right')
title(sprintf(' η = %0.f', params.eta))
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

% Migration speed
subplot(2, 3, 2)
plot(tt, system.Vm, 'LineWidth', 2)
xlabel('Time')
ylabel('Migration speed')
title(sprintf(' η = %0.f', params.eta))
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

% Total clutch force
subplot(2, 3, 4)
plot(tt, system.Fci_record_left, 'LineWidth', 2)
hold on
plot(tt, system.Fci_record_right, 'LineWidth', 2)
xlabel('Time')
ylabel('Total clutch force')
legend('Left', 'Right')
title(sprintf(' η = %0.f', params.eta))
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

% Migration distance
subplot(2, 3, 3)
plot(tt, system.Ds(1:num_steps), 'LineWidth', 2)  % Match length with tt
xlabel('Time')
ylabel('Migration distance')
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

% Adjusted unbinding rate
subplot(2, 3, 5)
plot(tt, system.roff_record, 'LineWidth', 2)
xlabel('Time')
ylabel('Koff')
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

% Substrate displacement
subplot(2, 3, 6)
plot(tt, left.xst, 'LineWidth', 2)
hold on
plot(tt, right.xst, 'LineWidth', 2)
xlabel('Time')
ylabel('Substrate displacement')
legend('Left', 'Right')
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);

%% Helper Functions
function [Pbi, Ncbound, Pbt, roffi] = update_bonds(Pbi, ronin, roff, Fci, Fb, dt, Nc)
    % Update clutch binding states and return roffi
    roffi = roff * exp(abs(Fci) ./ Fb);  % Force-dependent unbinding rate
    Pbi = Pbi + dt * ((1 - Pbi) .* ronin - Pbi .* roffi);  % Euler step
    Pbi = double(Pbi > rand(Nc, 1));     % Stochastic update
    Ncbound = sum(Pbi);                  % Count bound clutches
    Pbt = Ncbound / Nc;                  % Fraction bound
end

function [xst, xs2] = motor_clutch_step(dt, j, ks1, ks2, eta, kc, xci, Pbi, xs2, xst, Ncbound)
    % Compute substrate positions
    a1 = ks1 / (ks1 + ks2 + kc * Ncbound);
    a2 = eta / dt;
    xcsum = sum(xci .* Pbi);
    if j == 1
        xs2 = a1 * kc * xcsum / (ks1 - a1 * ks1 + a2);
    else
        xs2 = (a2 * xs2(j-1) + a1 * kc * xcsum) / (ks1 - a1 * ks1 + a2);
    end
    if Ncbound == 0 && j > 10
        xs2 = 0;
    end
    xst = (ks1 * xs2 + kc * xcsum) / (ks2 + ks1 + kc * Ncbound);
end