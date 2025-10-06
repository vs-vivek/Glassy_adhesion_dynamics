clc;
clear all;
close all;

%% Parameters
params.ks1 = 1;          
params.ks2 = 0.3;         
params.eta = 1000;        
params.kc = 5;            
params.Fb = 2.5;          
params.Fm = 2;            
params.vu = 120;          
params.vp = 120;          
params.ron = 1;           
params.roff = 0.2;        
params.Ftalin = 0.2;      
params.kcy = 0;           
params.Nc = 200;          
params.Nm = 200;          
params.Tfinal = 200;      
params.dt = 0.01;         
params.num_steps = ceil(params.Tfinal / params.dt);

%% Run 50 trajectories
num_trajectories = 50;
all_trap_times = [];

for i = 1:num_trajectories
    fprintf("Running trajectory %d/%d...\n", i, num_trajectories);
    [system, params] = run_conventional_model(params);

    % Compute total clutch force and time
    total_force = system.Fci_record_left + system.Fci_record_right;
    tt1 = (1:length(total_force)) * params.dt;

    % Threshold to detect collapse (can tune this if needed)
    force_threshold = 1e-2;

    % Find collapse indices
    collapse_idx = find(total_force < force_threshold);

    % Remove collapses too close to each other (less than 1s apart)
    min_separation = round(1 / params.dt);  % 1 second in steps
    filtered_idx = collapse_idx([true; diff(collapse_idx) > min_separation]);

    % Get trapping times
    collapse_times = tt1(filtered_idx);
    trap_times = diff(collapse_times);

    all_trap_times = [all_trap_times, trap_times];
end

%% Plot histogram
figure;
histogram(all_trap_times, 'Normalization', 'probability');
xlabel('Trapping time duration (s)');
ylabel('Probability of occurrence');
title('Trapping Time Distribution (50 Conventional Model Trajectories)');
set(gca, 'FontSize', 14, 'LineWidth', 2);

%% Conventional Model Function (Wrapped Original Code)
function [system, params] = run_conventional_model(params)
    num_steps = params.num_steps;

    % Initialize left side
    left.Pbi = zeros(params.Nc, 1);
    left.ronin = params.ron * ones(params.Nc, 1);
    left.Fci = zeros(params.Nc, 1);
    left.xci = zeros(params.Nc, 1);
    left.Ncbound = 0;
    left.Frec = 0;

    % Initialize right side
    right.Pbi = zeros(params.Nc, 1);
    right.ronin = params.ron * ones(params.Nc, 1);
    right.Fci = zeros(params.Nc, 1);
    right.xci = zeros(params.Nc, 1);
    right.Ncbound = 0;
    right.Frec = 0;

    % Preallocate
    left.xs2 = zeros(num_steps, 1);
    left.xst = zeros(num_steps, 1);
    left.Fmy = zeros(num_steps, 1);
    left.vf = zeros(num_steps, 1);
    left.Pbt = zeros(num_steps, 1);
    right.xs2 = zeros(num_steps, 1);
    right.xst = zeros(num_steps, 1);
    right.Fmy = zeros(num_steps, 1);
    right.vf = zeros(num_steps, 1);
    right.Pbt = zeros(num_steps, 1);

    system.Vm = zeros(num_steps, 1);
    system.Ds = zeros(num_steps + 1, 1);
    system.D_left = zeros(num_steps + 1, 1);
    system.D_right = zeros(num_steps + 1, 1);
    system.roff_record = zeros(num_steps, 1);
    system.Fci_record_left = zeros(num_steps, 1);
    system.Fci_record_right = zeros(num_steps, 1);
    system.roff_record_left = zeros(params.Nc, num_steps);
    system.roff_record_right = zeros(params.Nc, num_steps);
    system.Pbidoc_1 = zeros(num_steps, 2);

    system.Ds(1) = 0;
    system.D_left(1) = -3000;
    system.D_right(1) = 3000;

    roff = params.roff;

    % Main simulation loop
    for j = 1:num_steps
        [left.Pbi, left.Ncbound, left.Pbt(j), roffi_left] = update_bonds(...
            left.Pbi, left.ronin, roff, left.Fci, params.Fb, params.dt, params.Nc);
        [right.Pbi, right.Ncbound, right.Pbt(j), roffi_right] = update_bonds(...
            right.Pbi, right.ronin, roff, right.Fci, params.Fb, params.dt, params.Nc);

        system.roff_record_left(:, j) = roffi_left;
        system.roff_record_right(:, j) = roffi_right;
        system.Pbidoc_1(j, :) = [left.Pbi(1), right.Pbi(1)];

        left.Frec = left.Frec + left.Fmy(j) / (left.Ncbound + 1e-10);
        right.Frec = right.Frec + right.Fmy(j) / (right.Ncbound + 1e-10);

        [left.xst(j), left.xs2(j)] = motor_clutch_step(...
            params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
            left.xci, left.Pbi, left.xs2, left.xst, left.Ncbound);
        Fmy_left_old = sum(params.kc * (left.xci - left.xst(j)) .* left.Pbi);

        [right.xst(j), right.xs2(j)] = motor_clutch_step(...
            params.dt, j, params.ks1, params.ks2, params.eta, params.kc, ...
            right.xci, right.Pbi, right.xs2, right.xst, right.Ncbound);
        Fmy_right_old = sum(params.kc * (right.xci - right.xst(j)) .* right.Pbi);

        system.roff_record(j) = roff;

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

        left.xci = (1 - left.Pbi) .* left.xst(j) + left.Pbi .* (dxc_left + left.xci);
        right.xci = (1 - right.Pbi) .* right.xst(j) + right.Pbi .* (dxc_right + right.xci);
        left.Fci = params.kc * (left.xci - left.xst(j)) .* left.Pbi;
        right.Fci = params.kc * (right.xci - right.xst(j)) .* right.Pbi;
        left.Fmy(j) = Fmy_left_old + dxc_left * left.Ncbound * params.kc;
        right.Fmy(j) = Fmy_right_old + dxc_right * right.Ncbound * params.kc;

        system.Fci_record_left(j) = sum(left.Fci);
        system.Fci_record_right(j) = sum(right.Fci);

        if j < num_steps
            system.Vm(j) = (left.vf(j) - right.vf(j)) / 2;
            ra = [1, 1];
            system.D_left(j+1) = system.D_left(j) - params.dt * (params.vp * ra(1) - left.vf(j));
            system.D_right(j+1) = system.D_right(j) + params.dt * (params.vp * ra(2) - right.vf(j));
            system.Ds(j+1) = system.Ds(j) + params.dt * system.Vm(j);
        end
    end
end

%% Helper Functions
function [Pbi, Ncbound, Pbt, roffi] = update_bonds(Pbi, ronin, roff, Fci, Fb, dt, Nc)
    roffi = roff * exp(abs(Fci) ./ Fb);
    Pbi = Pbi + dt * ((1 - Pbi) .* ronin - Pbi .* roffi);
    Pbi = double(Pbi > rand(Nc, 1));
    Ncbound = sum(Pbi);
    Pbt = Ncbound / Nc;
end

function [xst, xs2] = motor_clutch_step(dt, j, ks1, ks2, eta, kc, xci, Pbi, xs2, xst, Ncbound)
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
