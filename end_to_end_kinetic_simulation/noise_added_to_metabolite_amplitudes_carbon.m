%% Hyperpolarized Metabolite Kinetic Analysis
% This script implements metabolite kinetic modeling for hyperpolarized 13C MRS
% focusing on pyruvate to bicarbonate conversion with various noise characteristics

% Usage: (1) download the UCSF hyperpolarized MRS toolkit and add it to your matlab path 
% (link: https://zenodo.org/records/10915690 -- we used version 1.7) 
%clear all;
close all;
base_path = fileparts(mfilename('fullpath'));            % end_to_end_kinetic_simulation/
addpath(genpath(base_path));                              % this folder + hyperpolarized-mri-toolbox
addpath(genpath(fullfile(fileparts(base_path),'OXSA'))); % OXSA (AMARES)

%% Simulation parameters for 13C 
k_true = 0.05;         % True conversion rate (1/s)
r1_pyr = 1/30;         % Pyruvate T1 relaxation rate (1/s)
r1_bic = 1/20;         % Bicarbonate T1 relaxation rate (1/s)
pyr_0 = 1.0;           % Initial pyruvate magnetization
TR = 3.0;              % Repetition time (s)
n_timepoints = 60;     % Number of time points
t = (0:n_timepoints-1)*TR;  % Time points

global rate_V
rate_V=2.94;            % Non-central noise parameters found semi-empirically
global scale_delta
scale_delta=0.041;      % for a chi^2 distribution

% Parameters for new noise types
global lognormal_sigma
lognormal_sigma = 0.6098475;  % Shape parameter for lognormal distribution – from experimental data 
global lognormal_meanlog
lognormal_meanlog = 4.1264617; 

global gamma_shape
gamma_shape = 3.3255825;      % Shape parameter for gamma distribution
global gamma_rate
gamma_rate = 0.0431788;       % Rate parameter for gamma distribution

% Flip angles (could be variable in practice)
flip_angle_deg_pyr = 15;  % degrees
flip_angle_deg_bic = 15;  % degrees
flip_angles_pyr = deg2rad(flip_angle_deg_pyr) * ones(1, n_timepoints);
flip_angles_bic = deg2rad(flip_angle_deg_bic) * ones(1, n_timepoints);
flip_angles = [flip_angles_pyr; flip_angles_bic];

% Noise parameters
noise_level = 0.1;      % Very roughly noise level relative to max signal
noise_types = {'gaussian', 'chi-squared', 'censored', 'lognormal', 'gamma'};
censoring_threshold = 0.05;  % For censored noise type

%% Generate true metabolite curves and add noise
rng(42);  % Set random seed for reproducibility

figure('Position', [100 100 1200 1000]);
set(gcf, 'color', 'w');
for i = 1:length(noise_types)
    % Generate simulated data
    [pyr_signal, bic_signal, pyr_true, bic_true] = generate_metabolite_signals(...
        k_true, r1_pyr, r1_bic, pyr_0, TR, n_timepoints, ...
        flip_angles_pyr, flip_angles_bic, noise_level, ...
        noise_types{i}, censoring_threshold);
    
    % Stack signals for fitting
    signals = zeros(1, 2, n_timepoints);
    signals(1, 1, :) = pyr_signal;
    signals(1, 2, :) = bic_signal;
    
    % Fixed parameters (fixing R1P as recommended)
    params_fixed = struct('R1P', r1_pyr);
    
    % Initial parameter estimates
    params_est = struct('kPB', 0.02, 'R1B', 1/25, 'B0_start', 0.0);
    
    % Fit using least squares method
    [params_fit, bic_fit, ~, err_metrics] = fit_pyr_kinetics(signals, TR, flip_angles, params_fixed, params_est);
    
    % Plot the data and fit
    subplot(5, 2, 2*i-1);
    plot(t, pyr_signal, 'b.', 'MarkerSize', 10); hold on;
    plot(t, bic_signal, 'r.', 'MarkerSize', 10);
    plot(t, pyr_true, 'b-', 'LineWidth', 2);
    plot(t, bic_true, 'r-', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Signal', 'FontSize', 12);
    title(['Simulated Data - ' noise_types{i} ' Noise'], 'FontSize', 14);
    legend('Pyruvate (noisy)', 'Bicarbonate (noisy)', 'Pyruvate (true)', 'Bicarbonate (true)');
    grid on;
    
    subplot(5, 2, 2*i);
    plot(t, bic_signal, 'r.', 'MarkerSize', 10); hold on;
    plot(t, bic_true, 'r-', 'LineWidth', 2);
    plot(t, reshape(bic_fit, [1, n_timepoints]), 'g--', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Signal', 'FontSize', 12);
    title(['Fitting Results - kPB=' num2str(params_fit.kPL, '%.4f') ' (true=' num2str(k_true, '%.4f') ')'], 'FontSize', 14);
    legend('Bicarbonate (noisy)', 'Bicarbonate (true)', 'Bicarbonate (fit)');
    grid on;
    
    % Display results
    disp(['Noise Type: ' noise_types{i}]);
    disp(['Estimated kPB: ' num2str(params_fit.kPL)]);
    disp(['Estimated R1B: ' num2str(params_fit.R1L)]);
    disp(['R^2: ' num2str(err_metrics.lactate.Rsq)]);
    disp('-------------------');
end

sgtitle('Metabolite Kinetic Fitting with Different Noise Types', 'FontSize', 16);
drawnow;
%% Monte Carlo simulation for error analysis
n_mc = 250;
k_estimates = zeros(length(noise_types), n_mc);

for i = 1:length(noise_types)
    for j = 1:n_mc
        % Generate data with noise
        [pyr_signal, bic_signal, ~, ~] = generate_metabolite_signals(...
            k_true, r1_pyr, r1_bic, pyr_0, TR, n_timepoints, ...
            flip_angles_pyr, flip_angles_bic, noise_level, ...
            noise_types{i}, censoring_threshold);
        
        % Stack signals for fitting
        signals = zeros(1, 2, n_timepoints);
        signals(1, 1, :) = pyr_signal;
        signals(1, 2, :) = bic_signal;
        
        % Fit the data
        [params_fit, ~, ~, ~] = fit_pyr_kinetics(signals, TR, flip_angles, params_fixed, params_est);
        
        % Store the estimated rate constant
        k_estimates(i, j) = params_fit.kPL;
    end
end

%% Plot Monte Carlo results
figure('Position', [100 100 1500 400]);
set(gcf, 'color', 'w');
for i = 1:length(noise_types)
    subplot(1, 5, i);
    histogram(k_estimates(i, :), 20, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    xline(k_true, 'r--', 'LineWidth', 2);
    xline(mean(k_estimates(i, :)), 'g-', 'LineWidth', 2);
    
    bias = mean(k_estimates(i, :)) - k_true;
    bias_percent = 100 * bias / k_true;
    std_dev = std(k_estimates(i, :));
    
    xlabel('Estimated kPB (1/s)', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title({[noise_types{i} ' Noise'], ...
        ['Bias: ' num2str(bias_percent, '%.1f') '%'], ...
        ['StdDev: ' num2str(std_dev, '%.4f')]}, 'FontSize', 14);
    grid on;
    
    % Display statistics
    disp(['Monte Carlo Results for ' noise_types{i} ' Noise:']);
    disp(['Mean Estimated kPB: ' num2str(mean(k_estimates(i, :)))]);
    disp(['Standard Deviation: ' num2str(std_dev)]);
    disp(['Bias: ' num2str(bias) ' (' num2str(bias_percent) '%)']);
    disp(['95% Confidence Interval: [' num2str(prctile(k_estimates(i, :), 2.5)) ', ' ...
        num2str(prctile(k_estimates(i, :), 97.5)) ']']);
    disp('-------------------');
    drawnow;
end

sgtitle(sprintf('Distribution of kPB Estimates (n=%d)',size(k_estimates,2)), 'FontSize', 16);
drawnow;
%% Effect of noise level on rate constant estimation with multiple simulations -- 13C 
%noise_levels = logspace(log10(0.05), log10(50), 20); % Noise levels from 0.05% to 50%
noise_levels = logspace(1,-3,200);
n_noise = length(noise_levels);
n_mc_per_noise_level = 100; % Number of Monte Carlo iterations per noise level

% Set up arrays to store results
k_vs_noise_mean = zeros(length(noise_types), n_noise);
k_vs_noise_std = zeros(length(noise_types), n_noise);
k_vs_noise_ci_lower = zeros(length(noise_types), n_noise);
k_vs_noise_ci_upper = zeros(length(noise_types), n_noise);
k_vs_noise_all = zeros(length(noise_types), n_noise, n_mc_per_noise_level);

% For each noise type and level, run multiple Monte Carlo simulations
for i = 1:length(noise_types)
    disp(['Processing noise type: ' noise_types{i}]);
    for j = 1:n_noise
        disp(['  Noise level: ' num2str(noise_levels(j)*100) '%']);
        for k = 1:n_mc_per_noise_level  % Using parfor for parallel execution if available
            % Generate data with different noise levels
            [pyr_signal, bic_signal, ~, ~] = generate_metabolite_signals(...
                k_true, r1_pyr, r1_bic, pyr_0, TR, n_timepoints, ...
                flip_angles_pyr, flip_angles_bic, noise_levels(j), ...
                noise_types{i}, censoring_threshold);
            
            % Stack signals for fitting
            signals = zeros(1, 2, n_timepoints);
            signals(1, 1, :) = pyr_signal;
            signals(1, 2, :) = bic_signal;
            
            % Fit the data
            [params_fit, ~, ~, ~] = fit_pyr_kinetics(signals, TR, flip_angles, params_fixed, params_est);
            
            % Store the estimated rate constant
            k_vs_noise_all(i, j, k) = params_fit.kPL;
        end
        
        % Calculate statistics for this noise type and level
        k_estimates = squeeze(k_vs_noise_all(i, j, :));
        k_vs_noise_mean(i, j) = mean(k_estimates);
        k_vs_noise_std(i, j) = std(k_estimates);
        
        % Calculate empirical 95% confidence intervals
        k_vs_noise_ci_lower(i, j) = prctile(k_estimates, 2.5);
        k_vs_noise_ci_upper(i, j) = prctile(k_estimates, 97.5);
    end
end

%% Plot results with shaded error bars
figure('Position', [100 100 800 600]);
set(gcf, 'color', 'w');
colors = {'b', 'g', 'r', 'm', 'c'};  % Added colors for new noise types
markers = {'o', 's', 'd', '^', 'v'};  % Added markers for new noise types
line_styles = {'-', '-', '-', '-', '-'};

for i = 1:length(noise_types)
    % Add shaded error region (95% confidence interval)
    x_region = [noise_levels*100, fliplr(noise_levels*100)];
    y_region = [k_vs_noise_ci_upper(i, :), fliplr(k_vs_noise_ci_lower(i, :))];
    
    % Create color with alpha for shaded region
    face_color = colors{i};
    face_alpha = 0.2;
    
    h = fill(x_region, y_region, face_color, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
    set(h, 'HandleVisibility', 'off'); % Don't show in legend
    hold on;
    
    % Plot mean estimates
    plot(noise_levels*100, k_vs_noise_mean(i, :), [colors{i} line_styles{i} markers{i}], 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', noise_types{i});
end

% Add true value line
yline(k_true, 'k--', 'LineWidth', 2, 'DisplayName', 'True kPB');

xlabel('Noise Level (%)', 'FontSize', 14);
ylabel('Estimated kPB (1/s)', 'FontSize', 14);
title('Effect of Noise Level on Rate Constant Estimation (with 95% CIs)', 'FontSize', 16);
legend('Location', 'northeast');
grid on;
set(gca, 'XScale', 'log');
xlim([min(noise_levels)*100, max(noise_levels)*100]);
ylim([0, 0.07]);

%% Helper functions
function [pyr_signal, bic_signal, pyr_true, bic_true] = generate_metabolite_signals(...
    k_true, r1_pyr, r1_bic, pyr_0, TR, n_timepoints, ...
    flip_angles_pyr, flip_angles_bic, noise_level, noise_type, censoring_threshold)
% Generate simulated pyruvate and bicarbonate signals with realistic noise.
% This is taken from the UCSF HP Toolkit. 

% Initialize arrays
pyr_mz = zeros(1, n_timepoints);
bic_mz = zeros(1, n_timepoints);
pyr_signal = zeros(1, n_timepoints);
bic_signal = zeros(1, n_timepoints);

% Initial conditions
pyr_mz(1) = pyr_0;
bic_mz(1) = 0.0;

% Apply first RF pulse
pyr_signal(1) = pyr_mz(1) * sin(flip_angles_pyr(1));
bic_signal(1) = bic_mz(1) * sin(flip_angles_bic(1));

% Update magnetization after first pulse
pyr_mz(1) = pyr_mz(1) * cos(flip_angles_pyr(1));
bic_mz(1) = bic_mz(1) * cos(flip_angles_bic(1));

% Propagate the system
for i = 2:n_timepoints
    % Magnetization evolution during TR
    pyr_next = pyr_mz(i-1) * exp(-(k_true + r1_pyr) * TR);
    bic_next = bic_mz(i-1) * exp(-r1_bic * TR) + ...
        k_true * pyr_mz(i-1) * (exp(-r1_bic * TR) - exp(-(k_true + r1_pyr) * TR)) / ...
        (k_true + r1_pyr - r1_bic + eps);  % Avoid division by zero
    
    % Store the magnetization
    pyr_mz(i) = pyr_next;
    bic_mz(i) = bic_next;
    
    % Apply RF pulse and get signal
    pyr_signal(i) = pyr_mz(i) * sin(flip_angles_pyr(i));
    bic_signal(i) = bic_mz(i) * sin(flip_angles_bic(i));
    
    % Update magnetization after RF pulse
    pyr_mz(i) = pyr_mz(i) * cos(flip_angles_pyr(i));
    bic_mz(i) = bic_mz(i) * cos(flip_angles_bic(i));
end

% Store the true signals before adding noise
pyr_true = pyr_signal;
bic_true = bic_signal;

% Add noise
max_pyr = max(pyr_signal);
max_bic = max(bic_signal);
global rate_V
global scale_delta
global lognormal_sigma
global lognormal_meanlog
global gamma_shape
global gamma_rate

%Parameter to normalise to unit SD...I don't think this is required IF you
%take the parameters for the distributions from data. 

normalise_sd = 1; 

if strcmp(noise_type, 'gaussian')
    % Gaussian noise of mean zero and standard deviation (noise_level * max_pyr)
    pyr_noise = noise_level * max_pyr * randn(1, n_timepoints);
    bic_noise = noise_level * max_bic * randn(1, n_timepoints);
    pyr_signal = pyr_signal + pyr_noise;
    bic_signal = bic_signal + bic_noise;
    
elseif strcmp(noise_type, 'chi-squared') 
    % Non-central chi-squared noise (approximating AMARES output)
    % Generate random numbers
    pn = ncx2rnd(rate_V, scale_delta, 1, n_timepoints);
    ln = ncx2rnd(rate_V, scale_delta, 1, n_timepoints);
    
    % Calculate standard deviations
    if normalise_sd
    std_pn = sqrt(2*(rate_V + 2*scale_delta));
    std_ln = sqrt(2*(rate_V + 2*scale_delta));
    else 
        std_pn = 1; 
        std_ln = 1; 
    end 
    % Normalize to be of variance one, then apply noise level scaling as
    % above
    pyr_noise = noise_level * max_pyr * (pn / std_pn);
    bic_noise = noise_level * max_bic * (ln / std_ln);
    
    pyr_signal = pyr_signal + pyr_noise;
    bic_signal = bic_signal + bic_noise;
    
elseif strcmp(noise_type, 'lognormal')
    % Lognormal noise with specified sigma parameter
    % Generate random numbers from lognormal distribution
    % mean and variance of underlying normal distribution
    %mu_pyr =  lognormal_meanlog - 0.5 * lognormal_sigma^2;  % Set mean of log(X) to have E[X]=1
    %mu_bic =  lognormal_meanlog - 0.5 * lognormal_sigma^2;
    
    %pyr_noise_raw = exp(mu_pyr + lognormal_sigma * randn(1, n_timepoints));
    %bic_noise_raw = exp(mu_bic + lognormal_sigma * randn(1, n_timepoints));
    pyr_noise_raw = lognrnd( lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);
    bic_noise_raw = lognrnd( lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);
    
    % Normalize to have unit standard deviation
    if normalise_sd
    std_pyr_noise = std(pyr_noise_raw);
    std_bic_noise = std(bic_noise_raw);
    else 
        std_pyr_noise=1; 
        std_bic_noise=1; 
    end
        
    pyr_noise = noise_level * max_pyr * (pyr_noise_raw / std_pyr_noise );
    bic_noise = noise_level * max_bic * (bic_noise_raw / std_bic_noise );
    
    pyr_signal = pyr_signal + pyr_noise;
    bic_signal = bic_signal + bic_noise;
    
elseif strcmp(noise_type, 'gamma')|| strcmp(noise_type, 'censored')
    % Gamma noise with specified shape and rate parameters
    % Generate random numbers from gamma distribution
    pyr_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);
    bic_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);
    
    % Normalize to have unit standard deviation
     if normalise_sd
    std_pyr_noise = std(pyr_noise_raw);
    std_bic_noise = std(bic_noise_raw);
     else 
           std_pyr_noise=1; 
        std_bic_noise=1; 
     end 
    pyr_noise = noise_level * max_pyr * (pyr_noise_raw / std_pyr_noise );
    bic_noise = noise_level * max_bic * (bic_noise_raw / std_bic_noise );
    
    pyr_signal = pyr_signal + pyr_noise;
    bic_signal = bic_signal + bic_noise;
end

% Ensure non-negative signals (as in real MR data)
pyr_signal = max(pyr_signal, 0);
bic_signal = max(bic_signal, 0);

is_censored=0; %Can flip on to see the effects on other distributions 

% Apply censoring if requested
if strcmp(noise_type, 'censored') || is_censored
    pyr_threshold = censoring_threshold * max_pyr;
    bic_threshold = censoring_threshold * max_bic;
    
    % Censor
    pyr_signal(pyr_signal < pyr_threshold) = 0;
    bic_signal(bic_signal < bic_threshold) = 0;
end
end


%% 