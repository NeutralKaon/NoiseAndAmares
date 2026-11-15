function [params_fit, lactate_fit, glucose_fit, err_metrics] = fit_dmi_kinetics(glucose_data, lactate_data, time_points, params_fixed, params_est)
%FIT_DMI_KINETICS Fit kinetic model for deuterium metabolic imaging (DMI)
%
% This function implements the kinetic model described in the Radiology paper
% for [6,6'-2H2]glucose to [3,3'-2H2]lactate conversion in tumor tissue.
%
% USAGE:
%   [params_fit, lactate_fit, glucose_fit, err_metrics] = fit_dmi_kinetics(...
%       glucose_data, lactate_data, time_points, params_fixed, params_est)
%
% INPUTS:
%   glucose_data - Vector of glucose concentrations over time (mM)
%   lactate_data - Vector of lactate concentrations over time (mM)
%   time_points  - Vector of time points (minutes)
%   params_fixed - Struct with fixed parameters (optional)
%                  Fields: f (fractional enrichment), Km (glucose transport Km)
%   params_est   - Struct with initial parameter estimates
%                  Fields: a1, kp, v, kg, Vmax, kl
%
% OUTPUTS:
%   params_fit   - Fitted parameters structure
%   lactate_fit  - Fitted lactate time course
%   glucose_fit  - Fitted glucose time course
%   err_metrics  - Error metrics (R^2, RMSE, etc.)
%
% MODEL EQUATIONS (from Radiology supplementary material):
%   Plasma glucose: G*p(t) = a1 * exp(-kp * t)
%   
%   Tissue glucose: dG*/dt = kg * (G*p - G*e) - f * Vglc
%   where Vglc = Vmax * G*e / (Km + G*e)
%   
%   Tissue lactate: dL*/dt = f * Vglc - kl * L*
%
% REFERENCE: 
%   Kreis et al, https://pubs.rsna.org/doi/10.1148/radiol.2019191242, Radiology 2020 (Supplementary Material)

%% Input validation and default parameters
if nargin < 4 || isempty(params_fixed)
    params_fixed = struct();
end

if nargin < 5 || isempty(params_est)
    params_est = struct();
end

% Set default fixed parameters
if ~isfield(params_fixed, 'f')
    params_fixed.f = 0.6;  % Fractional deuterium enrichment (from paper)
end
if ~isfield(params_fixed, 'Km')
    params_fixed.Km = 2.8;  % mM, glucose transport Km (from paper)
end

% Set default initial estimates
default_est = struct(...
    'a1', 17.0, ...      % Initial plasma glucose (mM)
    'kp', 0.042, ...     % Plasma glucose clearance (min^-1)
    'v', 0.7, ...        % Extracellular volume fraction
    'kg', 0.14, ...      % Glucose transfer rate (min^-1)
    'Vmax', 1.0, ...     % Maximum glucose consumption (mM/min)
    'kl', 0.12);         % Lactate clearance rate (min^-1)

field_names = fieldnames(default_est);
for i = 1:length(field_names)
    if ~isfield(params_est, field_names{i})
        params_est.(field_names{i}) = default_est.(field_names{i});
    end
end

%% Prepare data
time_points = time_points(:);
glucose_data = glucose_data(:);
lactate_data = lactate_data(:);

% Remove any NaN or infinite values
valid_idx = isfinite(glucose_data) & isfinite(lactate_data);
time_points = time_points(valid_idx);
glucose_data = glucose_data(valid_idx);
lactate_data = lactate_data(valid_idx);

%% Set up optimization
% Parameter vector: [a1, kp, v, kg, Vmax, kl]
x0 = [params_est.a1, params_est.kp, params_est.v, params_est.kg, params_est.Vmax, params_est.kl];

% Parameter bounds (all parameters should be positive)
lb = [0.1, 0.001, 0.01, 0.001, 0.001, 0.001];
ub = [100, 1.0, 1.0, 1.0, 10.0, 1.0];

% Optimization options
options = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 10000, ...
    'TolFun', 1e-8, ...
    'TolX', 1e-8);

% Define objective function
objective_fun = @(x) dmi_objective(x, time_points, glucose_data, lactate_data, params_fixed);

% Perform optimization
try
    [x_fit, resnorm, residual, exitflag] = lsqnonlin(objective_fun, x0, lb, ub, options);
    
    if exitflag <= 0
        warning('Optimization did not converge properly (exitflag = %d)', exitflag);
    end
catch ME
    warning('Optimization failed: %s', ME.message);
    x_fit = x0;  % Use initial guess if optimization fails
    resnorm = inf;
    residual = ones(size([glucose_data; lactate_data])) * inf;
end

%% Extract fitted parameters
params_fit = struct();
params_fit.a1 = x_fit(1);
params_fit.kp = x_fit(2);
params_fit.v = x_fit(3);
params_fit.kg = x_fit(4);
params_fit.Vmax = x_fit(5);
params_fit.kl = x_fit(6);

% Include fixed parameters
params_fit.f = params_fixed.f;
params_fit.Km = params_fixed.Km;

%% Generate fitted curves
[glucose_fit, lactate_fit] = solve_dmi_model(x_fit, time_points, params_fixed);

%% Calculate error metrics
% Split residuals
n_points = length(glucose_data);
glucose_residuals = residual(1:n_points);
lactate_residuals = residual((n_points+1):end);

% Calculate R-squared
glucose_ss_res = sum(glucose_residuals.^2);
glucose_ss_tot = sum((glucose_data - mean(glucose_data)).^2);
glucose_r2 = 1 - glucose_ss_res / (glucose_ss_tot + eps);

lactate_ss_res = sum(lactate_residuals.^2);
lactate_ss_tot = sum((lactate_data - mean(lactate_data)).^2);
lactate_r2 = 1 - lactate_ss_res / (lactate_ss_tot + eps);

% Store error metrics
err_metrics = struct();
err_metrics.glucose = struct();
err_metrics.glucose.Rsq = glucose_r2;
err_metrics.glucose.RMSE = sqrt(glucose_ss_res / n_points);
err_metrics.glucose.residuals = glucose_residuals;

err_metrics.lactate = struct();
err_metrics.lactate.Rsq = lactate_r2;
err_metrics.lactate.RMSE = sqrt(lactate_ss_res / n_points);
err_metrics.lactate.residuals = lactate_residuals;

err_metrics.overall_resnorm = resnorm;

end

function residuals = dmi_objective(x, time_points, glucose_data, lactate_data, params_fixed)
%DMI_OBJECTIVE Objective function for DMI kinetic fitting
    
    % Solve the model with current parameters
    [glucose_model, lactate_model] = solve_dmi_model(x, time_points, params_fixed);
    
    % Calculate residuals (weighted equally for glucose and lactate)
    glucose_residuals = (glucose_data - glucose_model) ./ (max(glucose_data) + eps);
    lactate_residuals = (lactate_data - lactate_model) ./ (max(lactate_data) + eps);
    
    % Combine residuals
    residuals = [glucose_residuals; lactate_residuals];
end

%% Example usage and demonstration
function demo_dmi_kinetics()
%DEMO_DMI_KINETICS Demonstrate the DMI kinetic fitting function

    fprintf('Running DMI kinetics demonstration...\n');
    
    % Time points (minutes)
    time_points = 0:2:80;  % Every 2 minutes for 80 minutes
    
    % True parameters for simulation
    true_params = struct(...
        'a1', 17.0, ...      % mM
        'kp', 0.042, ...     % min^-1  
        'v', 0.72, ...       % dimensionless
        'kg', 0.141, ...     % min^-1
        'Vmax', 0.985, ...   % mM/min
        'kl', 0.124, ...     % min^-1
        'f', 0.6, ...        % dimensionless
        'Km', 2.8);          % mM
    
    % Generate true curves
    params_vec = [true_params.a1, true_params.kp, true_params.v, ...
                  true_params.kg, true_params.Vmax, true_params.kl];
    params_fixed = struct('f', true_params.f, 'Km', true_params.Km);
    
    [glucose_true, lactate_true] = solve_dmi_model(params_vec, time_points, params_fixed);
    
    % Add realistic noise
    noise_level = 0.1;  % 10% noise
    glucose_noisy = glucose_true + noise_level * max(glucose_true) * randn(size(glucose_true));
    lactate_noisy = lactate_true + noise_level * max(lactate_true) * randn(size(lactate_true));
    
    % Ensure non-negative (as would be the case with AMARES)
    glucose_noisy = max(glucose_noisy, 0);
    lactate_noisy = max(lactate_noisy, 0);
    
    % Initial parameter estimates (perturbed from true values)
    params_est = struct(...
        'a1', 15.0, ...
        'kp', 0.05, ...
        'v', 0.6, ...
        'kg', 0.1, ...
        'Vmax', 0.8, ...
        'kl', 0.1);
    
    % Fit the model
    fprintf('Fitting DMI kinetic model...\n');
    [params_fit, lactate_fit, glucose_fit, err_metrics] = fit_dmi_kinetics(...
        glucose_noisy, lactate_noisy, time_points, params_fixed, params_est);
    
    % Display results
    fprintf('\nFitting Results:\n');
    fprintf('================\n');
    param_names = {'a1', 'kp', 'v', 'kg', 'Vmax', 'kl'};
    for i = 1:length(param_names)
        pname = param_names{i};
        fprintf('%6s: True = %8.4f, Fitted = %8.4f, Error = %6.2f%%\n', ...
            pname, true_params.(pname), params_fit.(pname), ...
            100 * abs(params_fit.(pname) - true_params.(pname)) / true_params.(pname));
    end
    
    fprintf('\nGoodness of Fit:\n');
    fprintf('Glucose R² = %.4f, RMSE = %.4f mM\n', err_metrics.glucose.Rsq, err_metrics.glucose.RMSE);
    fprintf('Lactate R² = %.4f, RMSE = %.4f mM\n', err_metrics.lactate.Rsq, err_metrics.lactate.RMSE);
    
    % Plot results
    figure('Position', [100, 100, 1200, 500]);
    
    subplot(1, 2, 1);
    plot(time_points, glucose_noisy, 'g.', 'MarkerSize', 10); hold on;
    plot(time_points, glucose_true, 'g-', 'LineWidth', 2);
    plot(time_points, glucose_fit, 'g--', 'LineWidth', 2);
    xlabel('Time (min)');
    ylabel('Concentration (mM)');
    title('[6,6''-²H₂]glucose');
    legend('Noisy data', 'True', 'Fitted', 'Location', 'best');
    grid on;
    
    subplot(1, 2, 2);
    plot(time_points, lactate_noisy, 'b.', 'MarkerSize', 10); hold on;
    plot(time_points, lactate_true, 'b-', 'LineWidth', 2);
    plot(time_points, lactate_fit, 'b--', 'LineWidth', 2);
    xlabel('Time (min)');
    ylabel('Concentration (mM)');
    title('[3,3''-²H₂]lactate');
    legend('Noisy data', 'True', 'Fitted', 'Location', 'best');
    grid on;
    
    sgtitle('DMI Kinetic Model Fitting Results');
    
    fprintf('\nDemo completed successfully!\n');
end

% Uncomment the line below to run the demonstration
% demo_dmi_kinetics();
