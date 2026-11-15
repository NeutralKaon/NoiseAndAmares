%% Test Pure Gaussian Noise Through Quantification Pipeline
% This script generates 100 spectra containing ONLY white Gaussian noise
% (no signal), quantifies them via both AMARES and integration, and
% saves the results for statistical analysis in R.
%
% Purpose: Verify what distribution pure Gaussian noise follows after
% passing through the quantification pipeline.

clear; close all;

% Add paths (resolved relative to this script's location)
base_path = fileparts(mfilename('fullpath'));      % StatisticalAnalysis/
pkg_root  = fileparts(base_path);                  % package root
addpath(genpath(fullfile(pkg_root,'end_to_end_kinetic_simulation')));  % PK_3T_*, quantify_spectra_pipline, hyperpolarized-mri-toolbox
addpath(genpath(fullfile(pkg_root,'OXSA')));                            % AMARES

%% Parameters
n_samples = 100;  % Number of pure noise spectra
n_points = 1024;  % FID length (reduced from 2048 to avoid exact zero averaging)
noise_std = 1.0;  % Standard deviation of noise

% Spectral parameters
dwell_time = 1/2000;  % 2 kHz spectral width
time_axis = (0:n_points-1)' * dwell_time;

%% Initialize storage
% For 13C: We expect 5 peaks (Pyruvate, Lactate, Bicarbonate, Alanine, Pyruvate-hydrate)
% For 2H: We expect 5 peaks (Water, Glucose, Glx, Lactate, + one more)
carbon_amares_amplitudes = zeros(n_samples, 5);
carbon_integration_amplitudes = zeros(n_samples, 5);
deuterium_amares_amplitudes = zeros(n_samples, 5);
deuterium_integration_amplitudes = zeros(n_samples, 5);

fprintf('Generating %d pure noise spectra...\n', n_samples);

%% Generate and quantify 13C noise spectra
fprintf('\n=== 13C Noise Quantification ===\n');
base_freq_13C = 128.3e6;  % 13C at 3T (Hz)

% Get prior knowledge structure for 13C
pk_struct_13C = PK_3T_1C_13C_5Peaks(0, 1);

% Create a compatible prior_knowledge structure for integration
% The integration method needs values.boundsCellArray
prior_knowledge_13C = struct();
prior_knowledge_13C.values = struct();
prior_knowledge_13C.values.boundsCellArray = {
    '1-13C Pyruvate',       [170, 172], [0, inf], [0, inf], [0, 360], [], [];
    '1-13C Lactate',        [183, 185], [0, inf], [0, inf], [0, 360], [], [];
    '1-13C Bicarbonate',    [160, 162], [0, inf], [0, inf], [0, 360], [], [];
    '1-13C Alanine',        [176, 178], [0, inf], [0, inf], [0, 360], [], [];
    '1-13C Pyruvate-Hydrate', [179, 181], [0, inf], [0, inf], [0, 360], [], [];
};

for i = 1:n_samples
    if mod(i, 10) == 0
        fprintf('  13C sample %d/%d\n', i, n_samples);
    end

    % Generate pure complex Gaussian noise
    noise_real = noise_std * randn(n_points, 1);
    noise_imag = noise_std * randn(n_points, 1);
    noise_fid = noise_real + 1i * noise_imag;

    % Package as fake_spectra structure for quantification
    fake_spectra_13C = struct();
    fake_spectra_13C.nucleus = '13C';
    fake_spectra_13C.noisy_fids = noise_fid;
    fake_spectra_13C.time_axis = time_axis;
    fake_spectra_13C.dwell_time = dwell_time;
    fake_spectra_13C.base_freq = base_freq_13C;
    fake_spectra_13C.prior_knowledge = prior_knowledge_13C;
    fake_spectra_13C.n_samples = n_points;
    fake_spectra_13C.pk = pk_struct_13C;  % For AMARES
    fake_spectra_13C.clip_negative = false;  % keep signed noise integrals
    fake_spectra_13C.metabolite_names = {'1-13C Pyruvate', '1-13C Lactate', ...
        '1-13C Bicarbonate', '1-13C Alanine', '1-13C Pyruvate-Hydrate'};

    % Add dummy true_amplitudes and time_points (required by quantify_spectra_pipline)
    fake_spectra_13C.true_amplitudes = zeros(5, 1);  % No signal, all zeros
    fake_spectra_13C.time_points = 0;

    % Calculate frequency axis
    fake_spectra_13C.freqAxis = ((-n_points/2):(n_points/2-1))' / (dwell_time * n_points);
    fake_spectra_13C.ppmAxis = fake_spectra_13C.freqAxis / (base_freq_13C * 1e-6);

    % Quantify via AMARES
    try
        quantified_amares = quantify_spectra_pipline(fake_spectra_13C, 'amares');
        if isfield(quantified_amares, 'quantified_amplitudes')
            carbon_amares_amplitudes(i, :) = quantified_amares.quantified_amplitudes(:, 1)';
        end
    catch ME
        fprintf('    AMARES failed for 13C sample %d: %s\n', i, ME.message);
        carbon_amares_amplitudes(i, :) = NaN;
    end

    % Quantify via Integration
    try
        quantified_integration = quantify_spectra_pipline(fake_spectra_13C, 'integration');
        if isfield(quantified_integration, 'quantified_amplitudes')
            carbon_integration_amplitudes(i, :) = quantified_integration.quantified_amplitudes(:, 1)';
            if all(carbon_integration_amplitudes(i,:) == 0)
                warning('All carbon amplitudes zero, try decreasing spectral resolution or bins...'); 
            end 
        end
    catch ME
        fprintf('    Integration failed for 13C sample %d: %s\n', i, ME.message);
        carbon_integration_amplitudes(i, :) = NaN;
    end
end

%% Generate and quantify 2H noise spectra
fprintf('\n=== 2H Noise Quantification ===\n');
base_freq_2H = 19.8e6;  % 2H at 3T (Hz)

% Get prior knowledge structure for 2H
pk_struct_2H = PK_3T_2H_5Peaks(0);

% Create a compatible prior_knowledge structure for integration
prior_knowledge_2H = struct();
prior_knowledge_2H.values = struct();
prior_knowledge_2H.values.boundsCellArray = {
    'Water',    [4.5, 4.9],   [0, inf], [0, inf], [0, 360], [], [];
    'Glucose',  [3.6, 4.0],   [0, inf], [0, inf], [0, 360], [], [];
    'Glx',      [2.1, 2.5],   [0, inf], [0, inf], [0, 360], [], [];
    'Lactate',  [0.7, 1.35],  [0, inf], [0, inf], [0, 360], [], [];
    'Other',    [0.0, 0.5],   [0, inf], [0, inf], [0, 360], [], [];
};

for i = 1:n_samples
    if mod(i, 10) == 0
        fprintf('  2H sample %d/%d\n', i, n_samples);
    end

    % Generate pure complex Gaussian noise
    noise_real = noise_std * randn(n_points, 1);
    noise_imag = noise_std * randn(n_points, 1);
    noise_fid = noise_real + 1i * noise_imag;

    % Package as fake_spectra structure for quantification
    fake_spectra_2H = struct();
    fake_spectra_2H.nucleus = '2H';
    fake_spectra_2H.noisy_fids = noise_fid;
    fake_spectra_2H.time_axis = time_axis;
    fake_spectra_2H.dwell_time = dwell_time;
    fake_spectra_2H.base_freq = base_freq_2H;
    fake_spectra_2H.prior_knowledge = prior_knowledge_2H;
    fake_spectra_2H.n_samples = n_points;
    fake_spectra_2H.pk = pk_struct_2H;  % For AMARES
    fake_spectra_2H.clip_negative = false;  % keep signed noise integrals
    fake_spectra_2H.metabolite_names = {'Water', 'Glucose', 'Glx', 'Lactate', 'Other'};

    % Add dummy true_amplitudes and time_points (required by quantify_spectra_pipline)
    fake_spectra_2H.true_amplitudes = zeros(5, 1);  % No signal, all zeros
    fake_spectra_2H.time_points = 0;

    % Calculate frequency axis
    fake_spectra_2H.freqAxis = ((-n_points/2):(n_points/2-1))' / (dwell_time * n_points);
    fake_spectra_2H.ppmAxis = fake_spectra_2H.freqAxis / (base_freq_2H * 1e-6);

    % Quantify via AMARES
    try
        quantified_amares = quantify_spectra_pipline(fake_spectra_2H, 'amares');
        if isfield(quantified_amares, 'quantified_amplitudes')
            deuterium_amares_amplitudes(i, :) = quantified_amares.quantified_amplitudes(:, 1)';
        end
    catch ME
        fprintf('    AMARES failed for 2H sample %d: %s\n', i, ME.message);
        deuterium_amares_amplitudes(i, :) = NaN;
    end

    % Quantify via Integration
    try
        quantified_integration = quantify_spectra_pipline(fake_spectra_2H, 'integration');
        if isfield(quantified_integration, 'quantified_amplitudes')
            deuterium_integration_amplitudes(i, :) = quantified_integration.quantified_amplitudes(:, 1)';
        end
    catch ME
        fprintf('    Integration failed for 2H sample %d: %s\n', i, ME.message);
        deuterium_integration_amplitudes(i, :) = NaN;
    end
end

%% Save results
output_file = 'noise_statistics_test_results.mat';
save(output_file, 'carbon_amares_amplitudes', 'carbon_integration_amplitudes', ...
    'deuterium_amares_amplitudes', 'deuterium_integration_amplitudes', ...
    'n_samples', 'noise_std');

fprintf('\nResults saved to %s\n', output_file);

%% Display summary statistics
fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('\n13C AMARES:\n');
fprintf('  Mean: %.4f (should be ~0 for pure noise)\n', mean(carbon_amares_amplitudes(:), 'omitnan'));
fprintf('  Std:  %.4f\n', std(carbon_amares_amplitudes(:), 'omitnan'));
fprintf('  Min:  %.4f (should be >=0 for magnitudes)\n', min(carbon_amares_amplitudes(:), [], 'omitnan'));

fprintf('\n13C Integration:\n');
fprintf('  Mean: %.4f\n', mean(carbon_integration_amplitudes(:), 'omitnan'));
fprintf('  Std:  %.4f\n', std(carbon_integration_amplitudes(:), 'omitnan'));
fprintf('  Min:  %.4f\n', min(carbon_integration_amplitudes(:), [], 'omitnan'));

fprintf('\n2H AMARES:\n');
fprintf('  Mean: %.4f\n', mean(deuterium_amares_amplitudes(:), 'omitnan'));
fprintf('  Std:  %.4f\n', std(deuterium_amares_amplitudes(:), 'omitnan'));
fprintf('  Min:  %.4f\n', min(deuterium_amares_amplitudes(:), [], 'omitnan'));

fprintf('\n2H Integration:\n');
fprintf('  Mean: %.4f\n', mean(deuterium_integration_amplitudes(:), 'omitnan'));
fprintf('  Std:  %.4f\n', std(deuterium_integration_amplitudes(:), 'omitnan'));
fprintf('  Min:  %.4f\n', min(deuterium_integration_amplitudes(:), [], 'omitnan'));

fprintf('\nNext step: Run noise_statistics_analysis.R to fit distributions\n');
