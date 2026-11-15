%% Deuterium Metabolic Imaging (DMI) Kinetic Analysis
% This script implements metabolite kinetic modeling for 2H DMI
% focusing on glucose to lactate conversion with various noise characteristics

% Usage: Requires fit_dmi_kinetics.m and solve_dmi_model.m in the path
%clear all;
close all;
base_path = fileparts(mfilename('fullpath'));            % end_to_end_kinetic_simulation/
addpath(genpath(base_path));                              % this folder + hyperpolarized-mri-toolbox
addpath(genpath(fullfile(fileparts(base_path),'OXSA'))); % OXSA (AMARES)

%% Simulation parameters for 2H DMI
% True kinetic parameters (based on Kreis et al. Radiology 2019, https://pubs.rsna.org/doi/10.1148/radiol.2019191242 -- see the SI and table E2)
a1_true = 17.0;          % Initial plasma glucose (mM)
kp_true = 0.042;         % Plasma glucose clearance (min^-1)
v_true = 0.72;           % Extracellular volume fraction
kg_true = 0.141;         % Glucose transfer rate (min^-1)
Vmax_true = 0.985;       % Maximum glucose consumption (mM/min)
kl_true = 0.124;         % Lactate clearance rate (min^-1)

% Fixed parameters
f = 0.6;                 % Fractional deuterium enrichment
Km = 2.8;                % Glucose transport Km (mM)

% Time parameters
TR = 2.0;                % Repetition time (minutes)
n_timepoints = 40;       % Number of time points
t = (0:n_timepoints-1)*TR;  % Time points (minutes)

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
    [glc_signal, lac_signal, glc_true, lac_true] = generate_dmi_signals(...
        a1_true, kp_true, v_true, kg_true, Vmax_true, kl_true, f, Km, ...
        t, noise_level, noise_types{i}, censoring_threshold);

    % Fixed parameters for fitting
    params_fixed = struct('f', f, 'Km', Km);

    % Initial parameter estimates
    params_est = struct('a1', 15.0, 'kp', 0.05, 'v', 0.6, 'kg', 0.1, 'Vmax', 0.8, 'kl', 0.1);

    % Fit using least squares method
    [params_fit, lac_fit, glc_fit, err_metrics] = fit_dmi_kinetics(...
        glc_signal, lac_signal, t, params_fixed, params_est);

    % Plot the data and fit
    subplot(5, 2, 2*i-1);
    plot(t, glc_signal, 'b.', 'MarkerSize', 10); hold on;
    plot(t, lac_signal, 'r.', 'MarkerSize', 10);
    plot(t, glc_true, 'b-', 'LineWidth', 2);
    plot(t, lac_true, 'r-', 'LineWidth', 2);
    xlabel('Time (min)', 'FontSize', 12);
    ylabel('Concentration (mM)', 'FontSize', 12);
    title(['Simulated Data - ' noise_types{i} ' Noise'], 'FontSize', 14);
    legend('Glucose (noisy)', 'Lactate (noisy)', 'Glucose (true)', 'Lactate (true)');
    grid on;

    subplot(5, 2, 2*i);
    plot(t, lac_signal, 'r.', 'MarkerSize', 10); hold on;
    plot(t, lac_true, 'r-', 'LineWidth', 2);
    plot(t, lac_fit, 'g--', 'LineWidth', 2);
    xlabel('Time (min)', 'FontSize', 12);
    ylabel('Concentration (mM)', 'FontSize', 12);
    title(['Fitting Results - Vmax=' num2str(params_fit.Vmax, '%.4f') ' (true=' num2str(Vmax_true, '%.4f') ')'], 'FontSize', 14);
    legend('Lactate (noisy)', 'Lactate (true)', 'Lactate (fit)');
    grid on;

    % Display results
    disp(['Noise Type: ' noise_types{i}]);
    disp(['Estimated Vmax: ' num2str(params_fit.Vmax)]);
    disp(['Estimated kg: ' num2str(params_fit.kg)]);
    disp(['Estimated kl: ' num2str(params_fit.kl)]);
    disp(['R^2 (Lactate): ' num2str(err_metrics.lactate.Rsq)]);
    disp(['R^2 (Glucose): ' num2str(err_metrics.glucose.Rsq)]);
    disp('-------------------');
end

sgtitle('DMI Kinetic Fitting with Different Noise Types', 'FontSize', 16);
drawnow;

%% Monte Carlo simulation for error analysis
n_mc = 250;
Vmax_estimates = zeros(length(noise_types), n_mc);

for i = 1:length(noise_types)
    for j = 1:n_mc
        % Generate data with noise
        [glc_signal, lac_signal, ~, ~] = generate_dmi_signals(...
            a1_true, kp_true, v_true, kg_true, Vmax_true, kl_true, f, Km, ...
            t, noise_level, noise_types{i}, censoring_threshold);

        % Fit the data
        [params_fit, ~, ~, ~] = fit_dmi_kinetics(...
            glc_signal, lac_signal, t, params_fixed, params_est);

        % Store the estimated Vmax
        Vmax_estimates(i, j) = params_fit.Vmax;
    end
end

%% Plot Monte Carlo results
figure('Position', [100 100 1500 400]);
set(gcf, 'color', 'w');
for i = 1:length(noise_types)
    subplot(1, 5, i);
    histogram(Vmax_estimates(i, :), 20, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    xline(Vmax_true, 'r--', 'LineWidth', 2);
    xline(mean(Vmax_estimates(i, :)), 'g-', 'LineWidth', 2);

    bias = mean(Vmax_estimates(i, :)) - Vmax_true;
    bias_percent = 100 * bias / Vmax_true;
    std_dev = std(Vmax_estimates(i, :));

    xlabel('Estimated Vmax (mM/min)', 'FontSize', 12);
    ylabel('Frequency', 'FontSize', 12);
    title({[noise_types{i} ' Noise'], ...
        ['Bias: ' num2str(bias_percent, '%.1f') '%'], ...
        ['StdDev: ' num2str(std_dev, '%.4f')]}, 'FontSize', 14);
    grid on;

    % Display statistics
    disp(['Monte Carlo Results for ' noise_types{i} ' Noise:']);
    disp(['Mean Estimated Vmax: ' num2str(mean(Vmax_estimates(i, :)))]);
    disp(['Standard Deviation: ' num2str(std_dev)]);
    disp(['Bias: ' num2str(bias) ' (' num2str(bias_percent) '%)']);
    disp(['95% Confidence Interval: [' num2str(prctile(Vmax_estimates(i, :), 2.5)) ', ' ...
        num2str(prctile(Vmax_estimates(i, :), 97.5)) ']']);
    disp('-------------------');
    drawnow;
end

sgtitle(sprintf('Distribution of Vmax Estimates (n=%d) - 2H DMI',size(Vmax_estimates,2)), 'FontSize', 16);
drawnow;

%% Effect of noise level on Vmax estimation with multiple simulations -- 2H DMI
noise_levels = logspace(1,-3,200);
n_noise = length(noise_levels);
n_mc_per_noise_level = 100; % Number of Monte Carlo iterations per noise level

% Set up arrays to store results
Vmax_vs_noise_mean = zeros(length(noise_types), n_noise);
Vmax_vs_noise_std = zeros(length(noise_types), n_noise);
Vmax_vs_noise_ci_lower = zeros(length(noise_types), n_noise);
Vmax_vs_noise_ci_upper = zeros(length(noise_types), n_noise);
Vmax_vs_noise_all = zeros(length(noise_types), n_noise, n_mc_per_noise_level);

% For each noise type and level, run multiple Monte Carlo simulations
for i = 1:length(noise_types)
    disp(['Processing noise type: ' noise_types{i}]);
    for j = 1:n_noise
        disp(['  Noise level: ' num2str(noise_levels(j)*100) '%']);
        for k = 1:n_mc_per_noise_level
            % Generate data with different noise levels
            [glc_signal, lac_signal, ~, ~] = generate_dmi_signals(...
                a1_true, kp_true, v_true, kg_true, Vmax_true, kl_true, f, Km, ...
                t, noise_levels(j), noise_types{i}, censoring_threshold);

            % Fit the data
            [params_fit, ~, ~, ~] = fit_dmi_kinetics(...
                glc_signal, lac_signal, t, params_fixed, params_est);

            % Store the estimated Vmax
            Vmax_vs_noise_all(i, j, k) = params_fit.Vmax;
        end

        % Calculate statistics for this noise type and level
        Vmax_estimates_temp = squeeze(Vmax_vs_noise_all(i, j, :));
        Vmax_vs_noise_mean(i, j) = mean(Vmax_estimates_temp);
        Vmax_vs_noise_std(i, j) = std(Vmax_estimates_temp);

        % Calculate empirical 95% confidence intervals
        Vmax_vs_noise_ci_lower(i, j) = prctile(Vmax_estimates_temp, 2.5);
        Vmax_vs_noise_ci_upper(i, j) = prctile(Vmax_estimates_temp, 97.5);
    end
end

%% Plot results with shaded error bars
figure('Position', [100 100 800 600]);
set(gcf, 'color', 'w');
colors = {'b', 'g', 'r', 'm', 'c'};
markers = {'o', 's', 'd', '^', 'v'};
line_styles = {'-', '-', '-', '-', '-'};

for i = 1:length(noise_types)
    % Add shaded error region (95% confidence interval)
    x_region = [noise_levels*100, fliplr(noise_levels*100)];
    y_region = [Vmax_vs_noise_ci_upper(i, :), fliplr(Vmax_vs_noise_ci_lower(i, :))];

    % Create color with alpha for shaded region
    face_color = colors{i};
    face_alpha = 0.2;

    h = fill(x_region, y_region, face_color, 'FaceAlpha', face_alpha, 'EdgeColor', 'none');
    set(h, 'HandleVisibility', 'off'); % Don't show in legend
    hold on;

    % Plot mean estimates
    plot(noise_levels*100, Vmax_vs_noise_mean(i, :), [colors{i} line_styles{i} markers{i}], 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', noise_types{i});
end

% Add true value line
yline(Vmax_true, 'k--', 'LineWidth', 2, 'DisplayName', 'True Vmax');

xlabel('Noise Level (%)', 'FontSize', 14);
ylabel('Estimated Vmax (mM/min)', 'FontSize', 14);
title('Effect of Noise Level on Vmax Estimation (with 95% CIs) - 2H DMI', 'FontSize', 16);
legend('Location', 'northeast');
grid on;
set(gca, 'XScale', 'log');
xlim([min(noise_levels)*100, max(noise_levels)*100]);
ylim([0, max(Vmax_true * 1.5, 1.5)]);

%% Helper functions
function [glc_signal, lac_signal, glc_true, lac_true] = generate_dmi_signals(...
    a1, kp, v, kg, Vmax, kl, f, Km, time_points, noise_level, noise_type, censoring_threshold)
% Generate simulated glucose and lactate signals with realistic noise.
% See the SI of Kreis et al, Radiology 2019 -- 10.1148/radiol.2019191242/

% Create parameter vector for solve_dmi_model
params = [a1, kp, v, kg, Vmax, kl];
params_fixed = struct('f', f, 'Km', Km);

% Generate true signals using the DMI model
[glc_true, lac_true] = solve_dmi_model(params, time_points, params_fixed);

% Make sure they are row vectors
glc_true = glc_true(:)';
lac_true = lac_true(:)';

% Add noise
max_glc = max(glc_true);
max_lac = max(lac_true);
n_timepoints = length(time_points);

global rate_V
global scale_delta
global lognormal_sigma
global lognormal_meanlog
global gamma_shape
global gamma_rate

% Parameter to normalise to unit SD
normalise_sd = 1;

if strcmp(noise_type, 'gaussian')
    % Gaussian noise of mean zero and standard deviation (noise_level * max_glc)
    glc_noise = noise_level * max_glc * randn(1, n_timepoints);
    lac_noise = noise_level * max_lac * randn(1, n_timepoints);
    glc_signal = glc_true + glc_noise;
    lac_signal = lac_true + lac_noise;

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
    % Normalize to be of variance one, then apply noise level scaling
    glc_noise = noise_level * max_glc * (pn / std_pn);
    lac_noise = noise_level * max_lac * (ln / std_ln);

    glc_signal = glc_true + glc_noise;
    lac_signal = lac_true + lac_noise;

elseif strcmp(noise_type, 'lognormal')
    % Lognormal noise with specified sigma parameter
    glc_noise_raw = lognrnd(lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);
    lac_noise_raw = lognrnd(lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);

    % Normalize to have unit standard deviation
    if normalise_sd
        std_glc_noise = std(glc_noise_raw);
        std_lac_noise = std(lac_noise_raw);
    else
        std_glc_noise = 1;
        std_lac_noise = 1;
    end

    glc_noise = noise_level * max_glc * (glc_noise_raw / std_glc_noise);
    lac_noise = noise_level * max_lac * (lac_noise_raw / std_lac_noise);

    glc_signal = glc_true + glc_noise;
    lac_signal = lac_true + lac_noise;

elseif strcmp(noise_type, 'gamma') || strcmp(noise_type, 'censored')
    % Gamma noise with specified shape and rate parameters
    glc_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);
    lac_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);

    % Normalize to have unit standard deviation
    if normalise_sd
        std_glc_noise = std(glc_noise_raw);
        std_lac_noise = std(lac_noise_raw);
    else
        std_glc_noise = 1;
        std_lac_noise = 1;
    end
    glc_noise = noise_level * max_glc * (glc_noise_raw / std_glc_noise);
    lac_noise = noise_level * max_lac * (lac_noise_raw / std_lac_noise);

    glc_signal = glc_true + glc_noise;
    lac_signal = lac_true + lac_noise;
end

% Ensure non-negative signals (as in real MR data)
glc_signal = max(glc_signal, 0);
lac_signal = max(lac_signal, 0);

is_censored = 0; % Can flip on to see the effects on other distributions

% Apply censoring if requested
if strcmp(noise_type, 'censored') || is_censored
    glc_threshold = censoring_threshold * max_glc;
    lac_threshold = censoring_threshold * max_lac;

    % Censor
    glc_signal(glc_signal < glc_threshold) = 0;
    lac_signal(lac_signal < lac_threshold) = 0;
end
end

%%
