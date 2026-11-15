%% Simulate the effect of noise on spectral quantification and variously used metrics
%% Enhanced Multi-Nuclear MRS AMARES vs Integration Noise Analysis Script
% Demonstrates effect of non-Gaussian noise on parameter estimation
% Compares AMARES vs integration quantification methods

%clear; close all; clc;

% Add necessary paths (resolved relative to this script's location)
base_path = fileparts(mfilename('fullpath'));   % end_to_end_kinetic_simulation/
addpath(genpath(fullfile(base_path,'hyperpolarized-mri-toolbox')));
addpath(genpath(fullfile(fileparts(base_path),'OXSA')));

% Set global noise parameters (from experimental data)
global rate_V scale_delta lognormal_sigma lognormal_meanlog gamma_shape gamma_rate
rate_V = 2.94;
scale_delta = 0.041;
lognormal_sigma = 0.6098475;
lognormal_meanlog = 4.1264617;
gamma_shape = 3.3255825;
gamma_rate = 0.0431788;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2H DMI AMARES vs Integration Noise Analysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting 2H DMI noise analysis (AMARES vs Integration)...\n');

%% DMI Simulation parameters
noise_levels_dmi = logspace(-3, 1, 30);  % 30 noise levels from 0.1% to 1000%
% Note that "percent" is perhaps misleading - rather than noise in the
% quantified spectral domain this is now noise in the FID directly.
n_mc_dmi = 50;  % Monte Carlo iterations per noise level

% True DMI parameters
dmi_params_true = struct(...
    'a1', 17.0, ...
    'kp', 0.0416, ...
    'v', 0.72, ...
    'kg', 0.141, ...
    'Vmax', 0.985, ...
    'kl', 0.124);

% Initialize DMI results storage
n_noise_dmi = length(noise_levels_dmi);
vmax_estimates_dmi_amares = zeros(n_noise_dmi, n_mc_dmi);
kg_estimates_dmi_amares = zeros(n_noise_dmi, n_mc_dmi);
vmax_estimates_dmi_integration = zeros(n_noise_dmi, n_mc_dmi);
kg_estimates_dmi_integration = zeros(n_noise_dmi, n_mc_dmi);

fprintf('Testing %d noise levels with %d Monte Carlo iterations each\n', n_noise_dmi, n_mc_dmi);


dmi_a1 = dmi_params_true.a1;
dmi_kp = dmi_params_true.kp;
dmi_v = dmi_params_true.v;
dmi_kg = dmi_params_true.kg;
dmi_Vmax = dmi_params_true.Vmax;
dmi_kl = dmi_params_true.kl;

%% DMI Main simulation loop
for noise_idx = 1:n_noise_dmi
    current_noise = noise_levels_dmi(noise_idx);
    fprintf('********************************************\n')
    fprintf('DMI: Processing Noise Power %.1f%% (%d/%d)\n', current_noise*100, noise_idx, n_noise_dmi);
    fprintf('********************************************\n')

    % Pre-allocate temporary arrays for this noise level
    vmax_estimates_dmi_amares_temp = zeros(1, n_mc_dmi);
    kg_estimates_dmi_amares_temp = zeros(1, n_mc_dmi);
    vmax_estimates_dmi_integration_temp = zeros(1, n_mc_dmi);
    kg_estimates_dmi_integration_temp = zeros(1, n_mc_dmi);
    fprintf('Starting MC loops...\n')
    for mc_idx = 1:n_mc_dmi
        try
            fprintf('MC iteration: %d, noise iteration %d\n',mc_idx, noise_idx)
            % Generate synthetic 2H spectra with current noise level
            fake_spectra_dmi_temp = generate_fake_fids('2H', 'noise_level', current_noise, ...
                'a1', dmi_a1, 'kp', dmi_kp, ...
                'v', dmi_v, 'kg', dmi_kg, ...
                'Vmax', dmi_Vmax, 'kl', dmi_kl);

            fake_spectra_dmi_temp.noise_level = current_noise;
            fake_spectra_dmi_temp.mc_id = mc_idx;

            % Quantify using AMARES
            quantified_results_dmi_amares_temp = quantify_spectra_pipline(fake_spectra_dmi_temp, 'amares');
            % Quantify using Integration
            quantified_results_dmi_integration_temp = quantify_spectra_pipline(fake_spectra_dmi_temp, 'integration');

            % %             A semi-empirically found "concentration" scale, obtained on a
            % %             noiseless spectrum from above, ii.e. run "generate_fake_fids"
            % %             with 0 as a noise level, and then look at the ratio of the
            % %             (constant!) D2O peak in amares and integration for the
            % %             chemical shift regions defined. This yields:
            % %             normFacDmi=mean(quantified_results_dmi_amares.quantified_amplitudes(1,:)./quantified_results_dmi_integration.quantified_amplitudes(1,:));
            % %             normFacDmi = 2.160012799913669e+02;

            % Extract reference peak for proper concentration normalization
            ref_idx = find(strcmp(fake_spectra_dmi_temp.metabolite_names, 'CDCl3'));
            if ~isempty(ref_idx)
                % Calculate normalization from constant reference peak
                ref_amares = quantified_results_dmi_amares_temp.quantified_amplitudes(ref_idx, :);
                ref_integration = quantified_results_dmi_integration_temp.quantified_amplitudes(ref_idx, :);

                % Use mean ratio (should be very consistent since reference is constant)
                normFacDmi = mean(ref_amares ./ ref_integration);
            else
                % Fail loudly: otherwise normFacDmi silently retains a stale
                % value from the previous MC iteration (or is undefined on the
                % first), corrupting the integration scaling.
                error('end_to_end:refPeakMissing', ...
                    'DMI reference peak ''CDCl3'' not found in metabolite_names; cannot normalise integration to AMARES.');
            end



            % Extract glucose and lactate timecourses for AMARES
            glucose_idx = find(strcmp(fake_spectra_dmi_temp.metabolite_names, 'Glucose'));
            lactate_idx = find(strcmp(fake_spectra_dmi_temp.metabolite_names, 'Lactate'));

            if ~isempty(glucose_idx) && ~isempty(lactate_idx)
                % AMARES results
                glucose_amares = quantified_results_dmi_amares_temp.quantified_amplitudes(glucose_idx, :);
                lactate_amares = quantified_results_dmi_amares_temp.quantified_amplitudes(lactate_idx, :);

                % Integration results
                glucose_integration = quantified_results_dmi_integration_temp.quantified_amplitudes(glucose_idx, :)*normFacDmi;
                lactate_integration = quantified_results_dmi_integration_temp.quantified_amplitudes(lactate_idx, :)*normFacDmi;
                %Used in case the kinetic estimates care about scale --
                %which they probably shouldn't given that absolute
                %quantification in MR is hard...

                % Fit DMI kinetics to AMARES-quantified data
                dmi_params_fixed = struct('f', 0.6, 'Km', 2.8);
                dmi_params_est = struct('a1', 15.0, 'kp', 0.05, 'v', 0.6, ...
                    'kg', 0.1, 'Vmax', 0.8, 'kl', 0.1);

                [dmi_params_fit_amares, ~, ~, ~] = fit_dmi_kinetics(...
                    glucose_amares', lactate_amares', fake_spectra_dmi_temp.time_points/60, ...
                    dmi_params_fixed, dmi_params_est);

                % Fit DMI kinetics to Integration-quantified data
                [dmi_params_fit_integration, ~, ~, ~] = fit_dmi_kinetics(...
                    glucose_integration', lactate_integration', fake_spectra_dmi_temp.time_points/60, ...
                    dmi_params_fixed, dmi_params_est);

                % Store AMARES results
                vmax_estimates_dmi_amares_temp(mc_idx) = dmi_params_fit_amares.Vmax;
                kg_estimates_dmi_amares_temp(mc_idx) = dmi_params_fit_amares.kg;

                % Store Integration results
                vmax_estimates_dmi_integration_temp(mc_idx) = dmi_params_fit_integration.Vmax;
                kg_estimates_dmi_integration_temp(mc_idx) = dmi_params_fit_integration.kg;
            else
                % If metabolites not found, store NaN
                vmax_estimates_dmi_amares_temp(mc_idx) = NaN;
                kg_estimates_dmi_amares_temp(mc_idx) = NaN;
                vmax_estimates_dmi_integration_temp(mc_idx) = NaN;
                kg_estimates_dmi_integration_temp(mc_idx) = NaN;
            end

        catch ME
            fprintf('DMI Error at Noise Power %.1f%%, iteration %d: %s\n', ...
                current_noise*100, mc_idx, ME.message);
            % keyboard;  % uncomment to debug
            vmax_estimates_dmi_amares_temp(mc_idx) = NaN;
            kg_estimates_dmi_amares_temp(mc_idx) = NaN;
            vmax_estimates_dmi_integration_temp(mc_idx) = NaN;
            kg_estimates_dmi_integration_temp(mc_idx) = NaN;
        end
    end

    % Store results from parfor loop
    vmax_estimates_dmi_amares(noise_idx, :) = vmax_estimates_dmi_amares_temp;
    kg_estimates_dmi_amares(noise_idx, :) = kg_estimates_dmi_amares_temp;
    vmax_estimates_dmi_integration(noise_idx, :) = vmax_estimates_dmi_integration_temp;
    kg_estimates_dmi_integration(noise_idx, :) = kg_estimates_dmi_integration_temp;
end

%% Calculate DMI statistics
vmax_mean_dmi_amares = nanmean(vmax_estimates_dmi_amares, 2);
vmax_std_dmi_amares = nanstd(vmax_estimates_dmi_amares, [], 2);
vmax_ci_lower_dmi_amares = prctile(vmax_estimates_dmi_amares, 2.5, 2);
vmax_ci_upper_dmi_amares = prctile(vmax_estimates_dmi_amares, 97.5, 2);

vmax_mean_dmi_integration = nanmean(vmax_estimates_dmi_integration, 2);
vmax_std_dmi_integration = nanstd(vmax_estimates_dmi_integration, [], 2);
vmax_ci_lower_dmi_integration = prctile(vmax_estimates_dmi_integration, 2.5, 2);
vmax_ci_upper_dmi_integration = prctile(vmax_estimates_dmi_integration, 97.5, 2);

kg_mean_dmi_amares = nanmean(kg_estimates_dmi_amares, 2);
kg_std_dmi_amares = nanstd(kg_estimates_dmi_amares, [], 2);
kg_ci_lower_dmi_amares = prctile(kg_estimates_dmi_amares, 2.5, 2);
kg_ci_upper_dmi_amares = prctile(kg_estimates_dmi_amares, 97.5, 2);


kg_mean_dmi_integration = nanmean(kg_estimates_dmi_integration, 2);
kg_std_dmi_integration = nanstd(kg_estimates_dmi_integration, [], 2);
kg_ci_lower_dmi_integration = prctile(kg_estimates_dmi_integration, 2.5,2);
kg_ci_upper_dmi_integration = prctile(kg_estimates_dmi_integration, 97.5,2);

%% Plot DMI results
% Example FIDs
figure('Position', [200 200 400 400]);
set(gcf,'color','w');
plot(fake_spectra_dmi_temp.ppmAxis, sum(real(specFft(fake_spectra_dmi_temp.clean_fids)),2),...
    'DisplayName','Noiseless Spectra','linewidth',2,'color','k');
hold on;
plot(fake_spectra_dmi_temp.ppmAxis, sum(real(specFft(fake_spectra_dmi_temp.noisy_fids)),2),...
    'DisplayName','With Gaussian Noise','linewidth',2,'color',[0.000 0.447 0.741 0.4]);
xlim([-5, 10]); set(gca,'XDir','reverse')
legend('Location', 'best'); grid on;
xlabel('Frequency (ppm)');
ylabel('Intensity (au)');
sgtitle('2H DMI: Generated Spectra and Noise', 'FontSize', 16);

% MC panel
figure('Position', [100 100 1200 800]); set(gcf, 'color', 'w');

% Vmax results with confidence intervals - AMARES
subplot(2, 2, 1);
semilogx(noise_levels_dmi, vmax_mean_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES Estimated Vmax');
yline(dmi_params_true.Vmax, 'k--', 'LineWidth', 2, 'DisplayName', 'True Vmax');
hold on;
x_region = [noise_levels_dmi, fliplr(noise_levels_dmi)];
y_region = [vmax_ci_upper_dmi_amares', fliplr(vmax_ci_lower_dmi_amares')];
fill(x_region, y_region, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none','DisplayName','95% CI');

% Add integration results
semilogx(noise_levels_dmi, vmax_mean_dmi_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration Estimated Vmax');
y_region_int = [vmax_ci_upper_dmi_integration', fliplr(vmax_ci_lower_dmi_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none','DisplayName','Integration 95% CI');

xlabel('Noise Power'); ylabel('Estimated Vmax (mM/min)');
title('2H DMI: Vmax Estimation vs Noise Power (AMARES vs Integration)');
legend('Location', 'best'); grid on;

% Vmax bias comparison
subplot(2, 2, 2);
bias_percent_amares = 100 * (vmax_mean_dmi_amares - dmi_params_true.Vmax) / dmi_params_true.Vmax;
bias_percent_integration = 100 * (vmax_mean_dmi_integration - dmi_params_true.Vmax) / dmi_params_true.Vmax;
semilogx(noise_levels_dmi, bias_percent_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
semilogx(noise_levels_dmi, bias_percent_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
yline(0, 'k--', 'LineWidth', 1);
xlabel('Noise Power'); ylabel('Vmax Bias (%)');
title('2H DMI: Vmax Bias vs Noise Level');
legend('Location', 'best'); grid on;

% kg results comparison
subplot(2, 2, 3);
semilogx(noise_levels_dmi, kg_mean_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES K_g');
hold on;
semilogx(noise_levels_dmi, kg_mean_dmi_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration K_g');
yline(dmi_params_true.kg, 'k--', 'LineWidth', 2, 'DisplayName', 'True K_g');
xlabel('Noise Power'); ylabel('Estimated K_g (min^{-1})');
title('2H DMI: K_g Estimation vs Noise Level');
legend('Location', 'best'); grid on;

x_region = [noise_levels_dmi, fliplr(noise_levels_dmi)];
y_region = [kg_ci_upper_dmi_amares', fliplr(kg_ci_lower_dmi_amares')];
fill(x_region, y_region, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none','DisplayName','95% CI');
y_region_int = [kg_ci_upper_dmi_integration', fliplr(kg_ci_lower_dmi_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none','DisplayName','Integration 95% CI');


% Distribution at 10% noise level
subplot(2, 2, 4);
noise_10_idx_dmi = find(abs(noise_levels_dmi - 0.1) < 0.05, 1);
if ~isempty(noise_10_idx_dmi)
    vmax_at_10_amares = vmax_estimates_dmi_amares(noise_10_idx_dmi, :);
    vmax_at_10_amares = vmax_at_10_amares(~isnan(vmax_at_10_amares));
    vmax_at_10_integration = vmax_estimates_dmi_integration(noise_10_idx_dmi, :);
    vmax_at_10_integration = vmax_at_10_integration(~isnan(vmax_at_10_integration));

    if ~isempty(vmax_at_10_amares) && ~isempty(vmax_at_10_integration)
        histogram(vmax_at_10_amares, 10, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'Normalization', 'probability','DisplayName','AMARES');
        hold on;
        histogram(vmax_at_10_integration, 10, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'Normalization', 'probability','DisplayName','Integration');
        xline(dmi_params_true.Vmax, 'k--', 'LineWidth', 2, 'DisplayName', 'True Vmax');
        xlabel('Estimated Vmax (mM/min)'); ylabel('Probability');
        title(['2H DMI: Vmax Distribution at ' sprintf('%.1f',noise_levels_dmi(noise_10_idx_dmi)*100) '% Noise']);
        legend('Location', 'best');
    end
end

sgtitle('2H DMI: AMARES vs Integration Quantification Comparison', 'FontSize', 16);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 13C Hyperpolarized MRS AMARES vs Integration Noise Analysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting 13C HP-MRS noise analysis (AMARES vs Integration)...\n');

%% 13C Simulation parameters
noise_levels_hp = logspace(-3, 1, 30);  % 30 noise levels from 0.1% to 1000 "%"
n_mc_hp = 50;  % Monte Carlo iterations per noise level

% True HP-13C kinetic parameters
k_true_hp = 0.05;         % True conversion rate (1/s)
r1_pyr_hp = 1/30;         % Pyruvate T1 relaxation rate (1/s)
r1_bic_hp = 1/20;         % Bicarbonate T1 relaxation rate (1/s)
pyr_0_hp = 1.0;           % Initial pyruvate magnetization
FlipAngle_hp = 15;        % Flip angle (degrees)

% Initialize 13C results storage
n_noise_hp = length(noise_levels_hp);
kpb_estimates_hp_amares = zeros(n_noise_hp, n_mc_hp);
r1b_estimates_hp_amares = zeros(n_noise_hp, n_mc_hp);
kpb_estimates_hp_integration = zeros(n_noise_hp, n_mc_hp);
r1b_estimates_hp_integration = zeros(n_noise_hp, n_mc_hp);

fprintf('Testing %d noise levels with %d Monte Carlo iterations each\n', n_noise_hp, n_mc_hp);

%% 13C Main simulation loop
for noise_idx = 1:n_noise_hp
    current_noise = noise_levels_hp(noise_idx);
    fprintf('********************************************\n')
    fprintf('13C: Processing Noise Power %.1f%% (%d/%d)\n', current_noise*100, noise_idx, n_noise_hp);
    fprintf('********************************************\n')

    % Pre-allocate temporary arrays for this noise level
    kpb_estimates_hp_amares_temp = zeros(1, n_mc_hp);
    r1b_estimates_hp_amares_temp = zeros(1, n_mc_hp);
    kpb_estimates_hp_integration_temp = zeros(1, n_mc_hp);
    r1b_estimates_hp_integration_temp = zeros(1, n_mc_hp);
    fprintf('Starting MC loops...\n');
    for mc_idx = 1:n_mc_hp
        try
            fprintf('MC iteration: %d, noise iteration %d\n',mc_idx, noise_idx)
            % Generate synthetic 13C spectra with current noise level
            fake_spectra_hp_temp = generate_fake_fids('13C', 'noise_level', current_noise, ...
                'k_true', k_true_hp, 'r1_pyr', r1_pyr_hp, 'r1_bic', r1_bic_hp, ...
                'pyr_0', pyr_0_hp, 'FlipAngle', FlipAngle_hp);

            % Quantify using AMARES
            quantified_results_hp_amares_temp = quantify_spectra_pipline(fake_spectra_hp_temp, 'amares');

            % Quantify using Integration
            quantified_results_hp_integration_temp = quantify_spectra_pipline(fake_spectra_hp_temp, 'integration');
% 
%             %Again, same argument as before but with pyruvate as a
%             %"reference":
%             normFacCarbon = 68.059772137922565;
%             %(found by generating noiseless fids and then)
%             %mean(quantified_results_hp_amares.quantified_amplitudes(1,:)./quantified_results_hp_integration.quantified_amplitudes(1,:))

 % Extract reference peak for proper concentration normalization
  ref_idx = find(strcmp(fake_spectra_hp_temp.metabolite_names, '1-13C UreaRef'));
  if ~isempty(ref_idx)
      % Calculate normalization from constant reference peak
      ref_amares = quantified_results_hp_amares_temp.quantified_amplitudes(ref_idx, :);
      ref_integration = quantified_results_hp_integration_temp.quantified_amplitudes(ref_idx, :);

      % Use mean ratio (should be very consistent since reference is constant)
      normFacCarbon = mean(ref_amares ./ ref_integration);
  else
      % Fail loudly: otherwise normFacCarbon silently retains a stale value
      % from the previous MC iteration (or is undefined on the first),
      % corrupting the integration scaling.
      error('end_to_end:refPeakMissing', ...
          '13C reference peak ''1-13C UreaRef'' not found in metabolite_names; cannot normalise integration to AMARES.');
  end



            % Extract pyruvate and bicarbonate timecourses
            pyr_idx = find(strcmp(fake_spectra_hp_temp.metabolite_names, '1-13C Pyruvate'));
            bic_idx = find(strcmp(fake_spectra_hp_temp.metabolite_names, '1-13C Bicarbonate'));

            if ~isempty(pyr_idx) && ~isempty(bic_idx)
                % AMARES results
                pyr_amares = quantified_results_hp_amares_temp.quantified_amplitudes(pyr_idx, :);
                bic_amares = quantified_results_hp_amares_temp.quantified_amplitudes(bic_idx, :);

                % Integration results
                pyr_integration = quantified_results_hp_integration_temp.quantified_amplitudes(pyr_idx, :)*normFacCarbon;
                bic_integration = quantified_results_hp_integration_temp.quantified_amplitudes(bic_idx, :)*normFacCarbon;

                % Prepare data for HP kinetic fitting - AMARES
                n_timepoints = length(fake_spectra_hp_temp.time_points);
                TR = fake_spectra_hp_temp.time_points(2) - fake_spectra_hp_temp.time_points(1);

                signals_amares = zeros(1, 2, n_timepoints);
                signals_amares(1, 1, :) = pyr_amares*normFacCarbon;
                signals_amares(1, 2, :) = bic_amares*normFacCarbon;

                signals_integration = zeros(1, 2, n_timepoints);
                signals_integration(1, 1, :) = pyr_integration*normFacCarbon;
                signals_integration(1, 2, :) = bic_integration*normFacCarbon;

                flip_angles_pyr = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                flip_angles_bic = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                flip_angles = [flip_angles_pyr; flip_angles_bic];

                % Fit HP kinetics to AMARES-quantified data
                params_fixed = struct('R1P', r1_pyr_hp);
                params_est = struct('kPL', 0.02, 'R1L', 1/25, 'B0_start', 0.0,'kPL_lb',0,'R1L_lb',0,'R1P_lb',0);

                [params_fit_amares, ~, ~, ~] = fit_pyr_kinetics(signals_amares, TR, ...
                    flip_angles, params_fixed, params_est);

                % Fit HP kinetics to Integration-quantified data
                [params_fit_integration, ~, ~, ~] = fit_pyr_kinetics(signals_integration, TR, ...
                    flip_angles, params_fixed, params_est);

                % Store AMARES results
                kpb_estimates_hp_amares_temp(mc_idx) = params_fit_amares.kPL;
                r1b_estimates_hp_amares_temp(mc_idx) = params_fit_amares.R1L;

                % Store Integration results
                kpb_estimates_hp_integration_temp(mc_idx) = params_fit_integration.kPL;
                r1b_estimates_hp_integration_temp(mc_idx) = params_fit_integration.R1L;
            else
                % If metabolites not found, store NaN
                kpb_estimates_hp_amares_temp(mc_idx) = NaN;
                r1b_estimates_hp_amares_temp(mc_idx) = NaN;
                kpb_estimates_hp_integration_temp(mc_idx) = NaN;
                r1b_estimates_hp_integration_temp(mc_idx) = NaN;
            end

        catch ME
            fprintf('13C Error at Noise Power %.1f%%, iteration %d: %s\n', ...
                current_noise*100, mc_idx, ME.message);
            % keyboard;  % uncomment to debug
            kpb_estimates_hp_amares_temp(mc_idx) = NaN;
            r1b_estimates_hp_amares_temp(mc_idx) = NaN;
            kpb_estimates_hp_integration_temp(mc_idx) = NaN;
            r1b_estimates_hp_integration_temp(mc_idx) = NaN;
        end
    end

    % Store results from parfor loop
    kpb_estimates_hp_amares(noise_idx, :) = kpb_estimates_hp_amares_temp;
    r1b_estimates_hp_amares(noise_idx, :) = r1b_estimates_hp_amares_temp;
    kpb_estimates_hp_integration(noise_idx, :) = kpb_estimates_hp_integration_temp;
    r1b_estimates_hp_integration(noise_idx, :) = r1b_estimates_hp_integration_temp;
end

%% Calculate 13C statistics
kpb_mean_hp_amares = nanmean(kpb_estimates_hp_amares, 2);
kpb_std_hp_amares = nanstd(kpb_estimates_hp_amares, [], 2);
kpb_ci_lower_hp_amares = prctile(kpb_estimates_hp_amares, 2.5, 2);
kpb_ci_upper_hp_amares = prctile(kpb_estimates_hp_amares, 97.5, 2);

kpb_mean_hp_integration = nanmean(kpb_estimates_hp_integration, 2);
kpb_std_hp_integration = nanstd(kpb_estimates_hp_integration, [], 2);
kpb_ci_lower_hp_integration = prctile(kpb_estimates_hp_integration, 2.5, 2);
kpb_ci_upper_hp_integration = prctile(kpb_estimates_hp_integration, 97.5, 2);

r1b_mean_hp_amares = nanmean(r1b_estimates_hp_amares, 2);
r1b_std_hp_amares = nanstd(r1b_estimates_hp_amares, [], 2);
r1b_ci_lower_hp_amares =  prctile(r1b_estimates_hp_amares, 2.5, 2);
r1b_ci_upper_hp_amares =  prctile(r1b_estimates_hp_amares, 97.5, 2);


r1b_mean_hp_integration = nanmean(r1b_estimates_hp_integration, 2);
r1b_std_hp_integration = nanstd(r1b_estimates_hp_integration, [], 2);
r1b_ci_hp_lower_integration = prctile(r1b_estimates_hp_integration, 2.5, 2);
r1b_ci_hp_upper_integration = prctile(r1b_estimates_hp_integration, 97.5, 2);

%% Plot 13C results
figure('Position', [200 200 400 400]);
set(gcf,'color','w');
plot(fake_spectra_hp_temp.ppmAxis, sum(real(specFft(fake_spectra_hp_temp.clean_fids)),2),...
    'DisplayName','Noiseless Spectra','linewidth',2,'color','k');
hold on;
plot(fake_spectra_hp_temp.ppmAxis, sum(real(specFft(fake_spectra_hp_temp.noisy_fids)),2),...
    'DisplayName','With Gaussian Noise','linewidth',2,'color',[0.000 0.447 0.741 0.4]);
set(gca,'XDir','reverse')
legend('Location', 'best'); grid on;
xlabel('Frequency (ppm)');
ylabel('Intensity (au)');
sgtitle('13C HP-MRS: Generated Spectra and Noise', 'FontSize', 16);

% MC panel
figure('Position', [100 100 1200 800]); set(gcf, 'color', 'w');

% kPB results with confidence intervals
subplot(2, 2, 1);
semilogx(noise_levels_hp, kpb_mean_hp_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES Estimated kPB');
yline(k_true_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True kPB');
hold on;
x_region = [noise_levels_hp, fliplr(noise_levels_hp)];
y_region = [kpb_ci_upper_hp_amares', fliplr(kpb_ci_lower_hp_amares')];
fill(x_region, y_region, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName','AMARES 95% CI');

% Add integration results
semilogx(noise_levels_hp, kpb_mean_hp_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration Estimated kPB');
y_region_int = [kpb_ci_upper_hp_integration', fliplr(kpb_ci_lower_hp_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName','Integration 95% CI');




xlabel('Noise Power'); ylabel('Estimated kPB (1/s)');
title('13C HP-MRS: kPB Estimation vs Noise Power (AMARES vs Integration)');
legend('Location', 'best'); grid on;

% kPB bias comparison
subplot(2, 2, 2);
bias_percent_hp_amares = 100 * (kpb_mean_hp_amares - k_true_hp) / k_true_hp;
bias_percent_hp_integration = 100 * (kpb_mean_hp_integration - k_true_hp) / k_true_hp;
yline(0, 'k--', 'LineWidth', 1);
semilogx(noise_levels_hp, bias_percent_hp_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
semilogx(noise_levels_hp, bias_percent_hp_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
xlabel('Noise Power'); ylabel('kPB Bias (%)');
title('13C HP-MRS: kPB Bias vs Noise Level');
legend('Location', 'best'); grid on;

% R1B results comparison
subplot(2, 2, 3);
semilogx(noise_levels_hp, r1b_mean_hp_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES R1B');
hold on;
semilogx(noise_levels_hp, r1b_mean_hp_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration R1B');
yline(r1_bic_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True R1B');
xlabel('Noise Power'); ylabel('Estimated R1B (1/s)');
title('13C HP-MRS: R1B Estimation vs Noise Level');
legend('Location', 'best'); grid on;
x_region = [noise_levels_hp, fliplr(noise_levels_hp)];
y_region = [r1b_ci_upper_hp_amares', fliplr(r1b_ci_lower_hp_amares')];
fill(x_region, y_region, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName','AMARES 95% CI');
y_region_int = [r1b_ci_hp_upper_integration', fliplr(r1b_ci_hp_lower_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName','Integration 95% CI');




% Distribution at 10% noise level
subplot(2, 2, 4);
noise_10_idx_hp = find(abs(noise_levels_hp - 0.1) < 0.05, 1);
if ~isempty(noise_10_idx_hp)
    kpb_at_10_amares = kpb_estimates_hp_amares(noise_10_idx_hp, :);
    kpb_at_10_amares = kpb_at_10_amares(~isnan(kpb_at_10_amares));
    kpb_at_10_integration = kpb_estimates_hp_integration(noise_10_idx_hp, :);
    kpb_at_10_integration = kpb_at_10_integration(~isnan(kpb_at_10_integration));

    if ~isempty(kpb_at_10_amares) && ~isempty(kpb_at_10_integration)
        histogram(kpb_at_10_amares, 10, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'Normalization', 'probability','DisplayName','AMARES');
        hold on;
        histogram(kpb_at_10_integration, 10, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'Normalization', 'probability','DisplayName','Integration');
        xline(k_true_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True kPB');
        xlabel('Estimated kPB (1/s)'); ylabel('Probability');
        title(sprintf('13C HP-MRS: kPB Distribution at %.1f%% Noise',noise_levels_hp(noise_10_idx_hp)*100));
        legend('Location', 'best');
    end
end

sgtitle('13C Hyperpolarized MRS: AMARES vs Integration Quantification Comparison', 'FontSize', 16);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Comparative Analysis
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nGenerating comprehensive comparative analysis...\n');

%% 1. AMARES vs Integration Comparison for 2H DMI
figure('Position', [50 50 1600 1200]); set(gcf, 'color', 'w');

% DMI Vmax: Method comparison
subplot(3, 4, 1);
loglog(noise_levels_dmi, vmax_mean_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_dmi, vmax_mean_dmi_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
yline(dmi_params_true.Vmax, 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
xlabel('Noise Power'); ylabel('Estimated Vmax (mM/min)');
title('2H DMI: Vmax Estimation');
legend('Location', 'best'); grid on;

% DMI Vmax: Bias comparison
subplot(3, 4, 2);
bias_dmi_amares = 100 * (vmax_mean_dmi_amares - dmi_params_true.Vmax) / dmi_params_true.Vmax;
bias_dmi_integration = 100 * (vmax_mean_dmi_integration - dmi_params_true.Vmax) / dmi_params_true.Vmax;
loglog(noise_levels_dmi, abs(bias_dmi_amares), 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_dmi, abs(bias_dmi_integration), 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
xlabel('Noise Power'); ylabel('|Bias| (%)');
title('2H DMI: Absolute Bias Comparison');
legend('Location', 'best'); grid on;

% DMI Vmax: Precision comparison (CV)
subplot(3, 4, 3); hold off;
cv_dmi_amares = 100 * vmax_std_dmi_amares ./ vmax_mean_dmi_amares;
cv_dmi_integration = 100 * vmax_std_dmi_integration ./ vmax_mean_dmi_integration;
loglog(noise_levels_dmi, cv_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_dmi, cv_dmi_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
xlabel('Noise Power'); ylabel('Coefficient of Variation (%)');
title('2H DMI: Precision Comparison');
legend('Location', 'northwest'); grid on;

% DMI: Bland altman plot
subplot(3, 4, 4); hold off;
x1 = vmax_mean_dmi_amares;
x2 = vmax_mean_dmi_integration;
% Calculate mean and difference
mean_vals = (x1 + x2) / 2;
diff_vals = x1 - x2; n = length(diff_vals);
mean_diff = mean(diff_vals); std_diff = std(diff_vals);
t_crit = tinv(1 - 0.05/2, n-1);  % t-critical value for 95% CI
upper_loa = mean_diff + t_crit * std_diff;
lower_loa = mean_diff - t_crit * std_diff;
scatter(mean_vals, diff_vals, 50, log10(noise_levels_dmi), 'filled');
hold on;
cb = colorbar;
cb.Label.String = 'Log_{10}(noise)';
%colormap(flipud(jet));  % Dark = high values, light = low values
yline(mean_diff, 'r--', 'LineWidth', 2);
yline(upper_loa, 'g--', 'LineWidth', 1.5);
yline(lower_loa, 'g--', 'LineWidth', 1.5);

% Add labels and formatting
xlabel('Mean of two methods');
ylabel('Difference (Amares - Integration)');
title('Bland-Altman Plot');
legend('Data points', sprintf('Mean diff = %.3f', mean_diff), ...
    sprintf('Upper LoA = %.3f', upper_loa), sprintf('Lower LoA = %.3f', lower_loa), ...
    'Location', 'best');
grid on;


% selected_noise_levels = [0.01, 0.05, 0.1, 0.2]; % 1%, 5%, 10%, 20%
% colors_noise = {'g', 'b', 'r', 'm'};
% for i = 1:length(selected_noise_levels)
%     noise_idx = find(abs(noise_levels_dmi - selected_noise_levels(i)) < 0.005, 1);
%     if ~isempty(noise_idx)
%         vmax_amares = vmax_estimates_dmi_amares(noise_idx, :);
%         vmax_integration = vmax_estimates_dmi_integration(noise_idx, :);
%         valid_idx = ~isnan(vmax_amares) & ~isnan(vmax_integration);
%         scatter(vmax_amares(valid_idx), vmax_integration(valid_idx), 50, colors_noise{i}, 'filled', ...
%             'DisplayName', sprintf('%.0f%% Noise', selected_noise_levels(i)*100));
%         hold on;
%     end
% end
% plot([0, 2], [0, 2], 'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
% xlabel('AMARES Vmax (mM/min)'); ylabel('Integration Vmax (mM/min)');
% title('2H DMI: Method Agreement');
legend('Location', 'best'); grid on; axis equal;

%% 2. AMARES vs Integration Comparison for 13C HP-MRS
% 13C kPB: Method comparison
subplot(3, 4, 5);
loglog(noise_levels_hp, kpb_mean_hp_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_hp, kpb_mean_hp_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
yline(k_true_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True Value');
xlabel('Noise Power'); ylabel('Estimated kPB (1/s)');
title('13C HP-MRS: kPB Estimation');
legend('Location', 'best'); grid on;

% 13C kPB: Bias comparison
subplot(3, 4, 6);
bias_hp_amares = 100 * (kpb_mean_hp_amares - k_true_hp) / k_true_hp;
bias_hp_integration = 100 * (kpb_mean_hp_integration - k_true_hp) / k_true_hp;
loglog(noise_levels_hp, abs(bias_hp_amares), 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_hp, abs(bias_hp_integration), 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
xlabel('Noise Power'); ylabel('|Bias| (%)');
title('13C HP-MRS: Absolute Bias Comparison');
legend('Location', 'best'); grid on;

% 13C kPB: Precision comparison (CV)
subplot(3, 4, 7);
cv_hp_amares = 100 * kpb_std_hp_amares ./ kpb_mean_hp_amares;
cv_hp_integration = 100 * kpb_std_hp_integration ./ kpb_mean_hp_integration;
loglog(noise_levels_hp, cv_hp_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
hold on;
loglog(noise_levels_hp, cv_hp_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
xlabel('Noise Power'); ylabel('Coefficient of Variation (%)');
title('13C HP-MRS: Precision Comparison');
legend('Location', 'northwest'); grid on;

% 13C: Method agreement at different noise levels
subplot(3, 4, 8);
% for i = 1:length(selected_noise_levels)
%     noise_idx = find(abs(noise_levels_hp - selected_noise_levels(i)) < 0.005, 1);
%     if ~isempty(noise_idx)
%         kpb_amares = kpb_estimates_hp_amares(noise_idx, :);
%         kpb_integration = kpb_estimates_hp_integration(noise_idx, :);
%         valid_idx = ~isnan(kpb_amares) & ~isnan(kpb_integration);
%         scatter(kpb_amares(valid_idx), kpb_integration(valid_idx), 50, colors_noise{i}, 'filled', ...
%             'DisplayName', sprintf('%.0f%% Noise', selected_noise_levels(i)*100));
%         hold on;
%     end
% end
% plot([0, 0.1], [0, 0.1], 'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
% xlabel('AMARES kPB (1/s)'); ylabel('Integration kPB (1/s)');
% title('13C HP-MRS: Method Agreement');
% legend('Location', 'best'); grid on; axis equal; xlim([0, 0.1]); ylim([0, 0.1]);
x1 = kpb_mean_hp_amares;
x2 = kpb_mean_hp_integration;
% Calculate mean and difference
mean_vals = (x1 + x2) / 2;
diff_vals = x1 - x2; n = length(diff_vals);
mean_diff = mean(diff_vals); std_diff = std(diff_vals);
t_crit = tinv(1 - 0.05/2, n-1);  % t-critical value for 95% CI
upper_loa = mean_diff + t_crit * std_diff;
lower_loa = mean_diff - t_crit * std_diff;
scatter(mean_vals, diff_vals, 50, log10(noise_levels_dmi), 'filled');
hold on;
cb = colorbar;
cb.Label.String = 'Log_{10}(noise)';
%colormap(flipud(jet));  % Dark = high values, light = low values
yline(mean_diff, 'r--', 'LineWidth', 2);
yline(upper_loa, 'g--', 'LineWidth', 1.5);
yline(lower_loa, 'g--', 'LineWidth', 1.5);

% Add labels and formatting
xlabel('Mean of two methods');
ylabel('Difference (Amares - Integration)');
title('Bland-Altman Plot');
legend('Data points', sprintf('Mean diff = %.3f', mean_diff), ...
    sprintf('Upper LoA = %.3f', upper_loa), sprintf('Lower LoA = %.3f', lower_loa), ...
    'Location', 'best');
grid on;
%% 3. Cross-Nuclear Comparison
% Method effectiveness across nuclei - Bias
subplot(3, 4, 9);
% Normalize biases to compare across different parameter scales
norm_bias_dmi_amares = abs(bias_dmi_amares);
norm_bias_dmi_integration = abs(bias_dmi_integration);
norm_bias_hp_amares = abs(bias_hp_amares);
norm_bias_hp_integration = abs(bias_hp_integration);

loglog(noise_levels_dmi, norm_bias_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
hold on;
loglog(noise_levels_dmi, norm_bias_dmi_integration, 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
loglog(noise_levels_hp, norm_bias_hp_amares, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
loglog(noise_levels_hp, norm_bias_hp_integration, 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');

%semilogx(noise_levels_dmi, norm_bias_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
%hold on;
%semilogx(noise_levels_dmi, norm_bias_dmi_integration, 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
%semilogx(noise_levels_hp, norm_bias_hp_amares, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
%semilogx(noise_levels_hp, norm_bias_hp_integration, 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
xlabel('Noise Power'); ylabel('Normalized |Bias| (%)');
title('Cross-Nuclear: Method Bias Comparison');
legend('Location', 'best'); grid on;

% Method effectiveness across nuclei - Precision
subplot(3, 4, 10);
loglog(noise_levels_dmi, cv_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
hold on;
loglog(noise_levels_dmi, cv_dmi_integration, 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
loglog(noise_levels_hp, cv_hp_amares, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
loglog(noise_levels_hp, cv_hp_integration, 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');

%semilogx(noise_levels_dmi, cv_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
%hold on;
%semilogx(noise_levels_dmi, cv_dmi_integration, 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
%semilogx(noise_levels_hp, cv_hp_amares, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
%semilogx(noise_levels_hp, cv_hp_integration, 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
xlabel('Noise Power'); ylabel('Coefficient of Variation (%)');
title('Cross-Nuclear: Precision Comparison');
legend('Location', 'northwest'); grid on;

% Relative method performance (Integration/AMARES ratio)
subplot(3, 4, 11);
ratio_bias_dmi = norm_bias_dmi_integration ./ (norm_bias_dmi_amares + eps);
ratio_bias_hp = norm_bias_hp_integration ./ (norm_bias_hp_amares + eps);
loglog(noise_levels_dmi, ratio_bias_dmi, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI');
hold on;
loglog(noise_levels_hp, ratio_bias_hp, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP-MRS');
yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
xlabel('Noise Power'); ylabel('Integration/AMARES Bias Ratio');
title('Relative Method Performance (Bias)');
legend('Location', 'best'); grid on;

% Root Mean Square Error comparison
subplot(3, 4, 12);
% Calculate RMSE for each method
rmse_dmi_amares = sqrt(mean((vmax_estimates_dmi_amares - dmi_params_true.Vmax).^2, 2, 'omitnan'));
rmse_dmi_integration = sqrt(mean((vmax_estimates_dmi_integration - dmi_params_true.Vmax).^2, 2, 'omitnan'));
rmse_hp_amares = sqrt(mean((kpb_estimates_hp_amares - k_true_hp).^2, 2, 'omitnan'));
rmse_hp_integration = sqrt(mean((kpb_estimates_hp_integration - k_true_hp).^2, 2, 'omitnan'));

% Normalize RMSE by true values for comparison
norm_rmse_dmi_amares = 100 * rmse_dmi_amares / dmi_params_true.Vmax;
norm_rmse_dmi_integration = 100 * rmse_dmi_integration / dmi_params_true.Vmax;
norm_rmse_hp_amares = 100 * rmse_hp_amares / k_true_hp;
norm_rmse_hp_integration = 100 * rmse_hp_integration / k_true_hp;

loglog(noise_levels_dmi, norm_rmse_dmi_amares, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
hold on;
loglog(noise_levels_dmi, norm_rmse_dmi_integration, 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
loglog(noise_levels_hp, norm_rmse_hp_amares, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
loglog(noise_levels_hp, norm_rmse_hp_integration, 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
xlabel('Noise Power'); ylabel('Normalized RMSE (%)');
title('Root Mean Square Error Comparison');
legend('Location', 'northwest'); grid on;

sgtitle('Comprehensive AMARES vs Integration Analysis: Cross-Nuclear Comparison', 'FontSize', 16);

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Rate Constants Correlation Analysis - NEW SECTION
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fprintf('\nStarting correlation analysis with varied rate constants...\n');
%
% %% Rate Constants Correlation Analysis Parameters
% noise_levels_corr = logspace(-2, 0, 15);  % 15 noise levels from 1% to 100%
% n_mc_corr = 20;  % Fewer Monte Carlo iterations for correlation analysis
%
% % Define parameter variation ranges (±50% from true values)
% dmi_param_variations = struct(...
%     'a1_range', [dmi_params_true.a1 * 0.5, dmi_params_true.a1 * 1.5], ...
%     'kp_range', [dmi_params_true.kp * 0.5, dmi_params_true.kp * 1.5], ...
%     'v_range', [dmi_params_true.v * 0.5, dmi_params_true.v * 1.5], ...
%     'kg_range', [dmi_params_true.kg * 0.5, dmi_params_true.kg * 1.5], ...
%     'Vmax_range', [dmi_params_true.Vmax * 0.5, dmi_params_true.Vmax * 1.5], ...
%     'kl_range', [dmi_params_true.kl * 0.5, dmi_params_true.kl * 1.5]);
%
% hp_param_variations = struct(...
%     'k_true_range', [k_true_hp * 0.5, k_true_hp * 1.5], ...
%     'r1_pyr_range', [r1_pyr_hp * 0.5, r1_pyr_hp * 1.5], ...
%     'r1_bic_range', [r1_bic_hp * 0.5, r1_bic_hp * 1.5]);
%
% % Initialize correlation results storage
% n_noise_corr = length(noise_levels_corr);
% n_params_dmi = 6;
% n_params_hp = 3;
%
% % Storage for DMI correlation analysis
% dmi_corr_vmax_amares = zeros(n_noise_corr, n_params_dmi);
% dmi_corr_kg_amares = zeros(n_noise_corr, n_params_dmi);
% dmi_corr_vmax_integration = zeros(n_noise_corr, n_params_dmi);
% dmi_corr_kg_integration = zeros(n_noise_corr, n_params_dmi);
%
% % Storage for HP correlation analysis
% hp_corr_kpb_amares = zeros(n_noise_corr, n_params_hp);
% hp_corr_r1b_amares = zeros(n_noise_corr, n_params_hp);
% hp_corr_kpb_integration = zeros(n_noise_corr, n_params_hp);
% hp_corr_r1b_integration = zeros(n_noise_corr, n_params_hp);
%
% fprintf('Testing %d noise levels with %d Monte Carlo iterations each for correlation analysis\n', n_noise_corr, n_mc_corr);
%
% %% Correlation Analysis
%
% % Initialize storage for correlation results
% correlation_results_dmi = struct();
% correlation_results_hp = struct();
%
% % Storage for correlation coefficients
% dmi_corr_vmax_ratios_amares = zeros(n_noise_corr, 1);
% dmi_corr_vmax_ratios_integration = zeros(n_noise_corr, 1);
% dmi_corr_kg_ratios_amares = zeros(n_noise_corr, 1);
% dmi_corr_kg_ratios_integration = zeros(n_noise_corr, 1);
% hp_corr_kpb_ratios_amares = zeros(n_noise_corr, 1);
% hp_corr_kpb_ratios_integration = zeros(n_noise_corr, 1);
%
% % Storage for ratio and parameter data for each noise level
% dmi_ratio_data_amares = cell(n_noise_corr, 1);
% dmi_ratio_data_integration = cell(n_noise_corr, 1);
% dmi_vmax_data_amares = cell(n_noise_corr, 1);
% dmi_vmax_data_integration = cell(n_noise_corr, 1);
% dmi_kg_data_amares = cell(n_noise_corr, 1);
% dmi_kg_data_integration = cell(n_noise_corr, 1);
%
% hp_ratio_data_amares = cell(n_noise_corr, 1);
% hp_ratio_data_integration = cell(n_noise_corr, 1);
% hp_kpb_data_amares = cell(n_noise_corr, 1);
% hp_kpb_data_integration = cell(n_noise_corr, 1);
%
% fprintf('Starting correlation analysis between spectral ratios and kinetic parameters...\n');
%
% %% DMI Correlation Analysis
% fprintf('Analyzing DMI correlations...\n');
% for noise_idx = 1:n_noise_corr
%     current_noise = noise_levels_corr(noise_idx);
%     fprintf('DMI Correlation: Processing Noise Level %.1f%% (%d/%d)\n', current_noise*100, noise_idx, n_noise_corr);
%
%     % Temporary storage for this noise level
%     dmi_ratios_amares_temp = [];
%     dmi_ratios_integration_temp = [];
%     dmi_vmax_fitted_amares_temp = [];
%     dmi_vmax_fitted_integration_temp = [];
%     dmi_kg_fitted_amares_temp = [];
%     dmi_kg_fitted_integration_temp = [];
%
%     for mc_idx = 1:n_mc_corr
%         try
%             % Generate synthetic 2H spectra with current noise level
%             fake_spectra_dmi_corr = generate_fake_fids('2H', 'noise_level', current_noise, ...
%                 'a1', dmi_a1, 'kp', dmi_kp, ...
%                 'v', dmi_v, 'kg', dmi_kg, ...
%                 'Vmax', dmi_Vmax, 'kl', dmi_kl);
%
%             % Quantify using AMARES
%             quantified_results_dmi_amares_corr = quantify_spectra_pipline(fake_spectra_dmi_corr, 'amares');
%
%             % Quantify using Integration
%             quantified_results_dmi_integration_corr = quantify_spectra_pipline(fake_spectra_dmi_corr, 'integration');
%
%             % Extract glucose and lactate for kinetic fitting
%             glucose_idx = find(strcmp(fake_spectra_dmi_corr.metabolite_names, 'Glucose'));
%             lactate_idx = find(strcmp(fake_spectra_dmi_corr.metabolite_names, 'Lactate'));
%
%             if ~isempty(glucose_idx) && ~isempty(lactate_idx)
%                 % AMARES kinetic fitting
%                 glucose_amares = quantified_results_dmi_amares_corr.quantified_amplitudes(glucose_idx, :);
%                 lactate_amares = quantified_results_dmi_amares_corr.quantified_amplitudes(lactate_idx, :);
%
%                 dmi_params_fixed = struct('f', 0.6, 'Km', 2.8);
%                 dmi_params_est = struct('a1', 15.0, 'kp', 0.05, 'v', 0.6, ...
%                     'kg', 0.1, 'Vmax', 0.8, 'kl', 0.1);
%
%                 [dmi_params_fit_amares_corr, ~, ~, ~] = fit_dmi_kinetics(...
%                     glucose_amares', lactate_amares', fake_spectra_dmi_corr.time_points/60, ...
%                     dmi_params_fixed, dmi_params_est);
%
%                 % Integration kinetic fitting
%                 glucose_integration = quantified_results_dmi_integration_corr.quantified_amplitudes(glucose_idx, :)*normFacDmi;
%                 lactate_integration = quantified_results_dmi_integration_corr.quantified_amplitudes(lactate_idx, :)*normFacDmi;
%
%                 [dmi_params_fit_integration_corr, ~, ~, ~] = fit_dmi_kinetics(...
%                     glucose_integration', lactate_integration', fake_spectra_dmi_corr.time_points/60, ...
%                     dmi_params_fixed, dmi_params_est);
%
%                 % Extract ratios of interest (Lac/Total 2H signal ratio)
%                 if isfield(quantified_results_dmi_amares_corr, 'ratios_of_interest') && ...
%                    isfield(quantified_results_dmi_integration_corr, 'ratios_of_interest')
%
%                     % Use Lac/Total ratio (ratios_of_interest(2))
%                     dmi_ratios_amares_temp(end+1) = quantified_results_dmi_amares_corr.ratios_of_interest(2);
%                     dmi_ratios_integration_temp(end+1) = quantified_results_dmi_integration_corr.ratios_of_interest(2);
%                     dmi_vmax_fitted_amares_temp(end+1) = dmi_params_fit_amares_corr.Vmax;
%                     dmi_vmax_fitted_integration_temp(end+1) = dmi_params_fit_integration_corr.Vmax;
%                     dmi_kg_fitted_amares_temp(end+1) = dmi_params_fit_amares_corr.kg;
%                     dmi_kg_fitted_integration_temp(end+1) = dmi_params_fit_integration_corr.kg;
%                 else
%                     keyboard;
%                 end
%             end
%
%         catch ME
%             fprintf('DMI Correlation Error at Noise %.1f%%, MC %d: %s\n', current_noise*100, mc_idx, ME.message);
%             keyboard;
%             continue;
%         end
%     end
%
%     % Store data for this noise level
%     dmi_ratio_data_amares{noise_idx} = dmi_ratios_amares_temp;
%     dmi_ratio_data_integration{noise_idx} = dmi_ratios_integration_temp;
%     dmi_vmax_data_amares{noise_idx} = dmi_vmax_fitted_amares_temp;
%     dmi_vmax_data_integration{noise_idx} = dmi_vmax_fitted_integration_temp;
%     dmi_kg_data_amares{noise_idx} = dmi_kg_fitted_amares_temp;
%     dmi_kg_data_integration{noise_idx} = dmi_kg_fitted_integration_temp;
%
%     % Calculate correlations if we have enough data
%     % NOTE: These correlations are EXPECTED to be low (~0.2) because:
%     %   - The TRUE Vmax is CONSTANT across all MC iterations (same biological state)
%     %   - Both ratio and fitted Vmax vary only due to MEASUREMENT NOISE
%     %   - Low correlation indicates that noise affects them independently
%     %   - This tells us: ratios are NOT good for recovering kinetic parameters
%     %     from a single noisy measurement at a given biological state
%     % For comparison, see Parameter Variation section below where TRUE Vmax varies
%     % and correlations are high (~1.0), showing ratios DO track real biological changes
%     if length(dmi_ratios_amares_temp) > 2 && length(dmi_vmax_fitted_amares_temp) > 2
%         corr_matrix_amares = corrcoef(dmi_ratios_amares_temp, dmi_vmax_fitted_amares_temp, 'rows', 'complete');
%         if size(corr_matrix_amares, 1) == 2
%             dmi_corr_vmax_ratios_amares(noise_idx) = corr_matrix_amares(1, 2);
%         end
%     end
%
%     if length(dmi_ratios_integration_temp) > 2 && length(dmi_vmax_fitted_integration_temp) > 2
%         corr_matrix_integration = corrcoef(dmi_ratios_integration_temp, dmi_vmax_fitted_integration_temp, 'rows', 'complete');
%         if size(corr_matrix_integration, 1) == 2
%             dmi_corr_vmax_ratios_integration(noise_idx) = corr_matrix_integration(1, 2);
%         end
%     end
%
%     % Calculate K_g correlations if we have enough data
%     if length(dmi_ratios_amares_temp) > 2 && length(dmi_kg_fitted_amares_temp) > 2
%         corr_matrix_kg_amares = corrcoef(dmi_ratios_amares_temp, dmi_kg_fitted_amares_temp, 'rows', 'complete');
%         if size(corr_matrix_kg_amares, 1) == 2
%             dmi_corr_kg_ratios_amares(noise_idx) = corr_matrix_kg_amares(1, 2);
%         end
%     end
%
%     if length(dmi_ratios_integration_temp) > 2 && length(dmi_kg_fitted_integration_temp) > 2
%         corr_matrix_kg_integration = corrcoef(dmi_ratios_integration_temp, dmi_kg_fitted_integration_temp, 'rows', 'complete');
%         if size(corr_matrix_kg_integration, 1) == 2
%             dmi_corr_kg_ratios_integration(noise_idx) = corr_matrix_kg_integration(1, 2);
%         end
%     end
% end
%
% %% 13C HP Correlation Analysis
% fprintf('Analyzing 13C HP correlations...\n');
% for noise_idx = 1:n_noise_corr
%     current_noise = noise_levels_corr(noise_idx);
%     fprintf('13C HP Correlation: Processing Noise Level %.1f%% (%d/%d)\n', current_noise*100, noise_idx, n_noise_corr);
%
%     % Temporary storage for this noise level
%     hp_ratios_amares_temp = [];
%     hp_ratios_integration_temp = [];
%     hp_kpb_fitted_amares_temp = [];
%     hp_kpb_fitted_integration_temp = [];
%
%     for mc_idx = 1:n_mc_corr
%         try
%             % Generate synthetic 13C spectra with current noise level
%             fake_spectra_hp_corr = generate_fake_fids('13C', 'noise_level', current_noise, ...
%                 'k_true', k_true_hp, 'r1_pyr', r1_pyr_hp, 'r1_bic', r1_bic_hp, ...
%                 'pyr_0', pyr_0_hp, 'FlipAngle', FlipAngle_hp);
%
%             % Quantify using AMARES
%             quantified_results_hp_amares_corr = quantify_spectra_pipline(fake_spectra_hp_corr, 'amares');
%
%             % Quantify using Integration
%             quantified_results_hp_integration_corr = quantify_spectra_pipline(fake_spectra_hp_corr, 'integration');
%
%             % Extract pyruvate and bicarbonate for kinetic fitting
%             pyr_idx = find(strcmp(fake_spectra_hp_corr.metabolite_names, '1-13C Pyruvate'));
%             bic_idx = find(strcmp(fake_spectra_hp_corr.metabolite_names, '1-13C Bicarbonate'));
%
%             if ~isempty(pyr_idx) && ~isempty(bic_idx)
%                 % AMARES kinetic fitting
%                 pyr_amares = quantified_results_hp_amares_corr.quantified_amplitudes(pyr_idx, :);
%                 bic_amares = quantified_results_hp_amares_corr.quantified_amplitudes(bic_idx, :);
%
%                 n_timepoints = length(fake_spectra_hp_corr.time_points);
%                 TR = fake_spectra_hp_corr.time_points(2) - fake_spectra_hp_corr.time_points(1);
%
%                 signals_amares = zeros(1, 2, n_timepoints);
%                 signals_amares(1, 1, :) = pyr_amares*normFacCarbon;
%                 signals_amares(1, 2, :) = bic_amares*normFacCarbon;
%
%                 flip_angles_pyr = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
%                 flip_angles_bic = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
%                 flip_angles = [flip_angles_pyr; flip_angles_bic];
%
%                 params_fixed = struct('R1P', r1_pyr_hp);
%                 params_est = struct('kPL', 0.02, 'R1L', 1/25, 'B0_start', 0.0,'kPL_lb',0,'R1L_lb',0,'R1P_lb',0);
%
%                 [params_fit_amares_corr, ~, ~, ~] = fit_pyr_kinetics(signals_amares, TR, ...
%                     flip_angles, params_fixed, params_est);
%
%                 % Integration kinetic fitting
%                 pyr_integration = quantified_results_hp_integration_corr.quantified_amplitudes(pyr_idx, :)*normFacCarbon;
%                 bic_integration = quantified_results_hp_integration_corr.quantified_amplitudes(bic_idx, :)*normFacCarbon;
%
%                 signals_integration = zeros(1, 2, n_timepoints);
%                 signals_integration(1, 1, :) = pyr_integration*normFacCarbon;
%                 signals_integration(1, 2, :) = bic_integration*normFacCarbon;
%
%                 [params_fit_integration_corr, ~, ~, ~] = fit_pyr_kinetics(signals_integration, TR, ...
%                     flip_angles, params_fixed, params_est);
%
%                 % Extract ratios of interest (Bic/Pyr ratio)
%                 if isfield(quantified_results_hp_amares_corr, 'ratios_of_interest') && ...
%                    isfield(quantified_results_hp_integration_corr, 'ratios_of_interest')
%
%                     % Use Bic/Pyr ratio (ratios_of_interest(2))
%                     hp_ratios_amares_temp(end+1) = quantified_results_hp_amares_corr.ratios_of_interest(2);
%                     hp_ratios_integration_temp(end+1) = quantified_results_hp_integration_corr.ratios_of_interest(2);
%                     hp_kpb_fitted_amares_temp(end+1) = params_fit_amares_corr.kPL;
%                     hp_kpb_fitted_integration_temp(end+1) = params_fit_integration_corr.kPL;
%                 end
%             end
%
%         catch ME
%             fprintf('13C HP Correlation Error at Noise %.1f%%, MC %d: %s\n', current_noise*100, mc_idx, ME.message);
%             keyboard
%             continue;
%         end
%     end
%
%     % Store data for this noise level
%     hp_ratio_data_amares{noise_idx} = hp_ratios_amares_temp;
%     hp_ratio_data_integration{noise_idx} = hp_ratios_integration_temp;
%     hp_kpb_data_amares{noise_idx} = hp_kpb_fitted_amares_temp;
%     hp_kpb_data_integration{noise_idx} = hp_kpb_fitted_integration_temp;
%
%     % Calculate correlations if we have enough data
%     if length(hp_ratios_amares_temp) > 2 && length(hp_kpb_fitted_amares_temp) > 2
%         corr_matrix_amares = corrcoef(hp_ratios_amares_temp, hp_kpb_fitted_amares_temp, 'rows', 'complete');
%         if size(corr_matrix_amares, 1) == 2
%             hp_corr_kpb_ratios_amares(noise_idx) = corr_matrix_amares(1, 2);
%         end
%     end
%
%     if length(hp_ratios_integration_temp) > 2 && length(hp_kpb_fitted_integration_temp) > 2
%         corr_matrix_integration = corrcoef(hp_ratios_integration_temp, hp_kpb_fitted_integration_temp, 'rows', 'complete');
%         if size(corr_matrix_integration, 1) == 2
%             hp_corr_kpb_ratios_integration(noise_idx) = corr_matrix_integration(1, 2);
%         end
%     end
% end
%
% %% Plot Correlation Results
% fprintf('Generating correlation plots...\n');
%
% % Create comprehensive correlation figure
% figure('Position', [50 50 1600 1000]); set(gcf, 'color', 'w');
%
% % DMI Correlation plots
% subplot(2, 3, 1);
% semilogx(noise_levels_corr, dmi_corr_vmax_ratios_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
% hold on;
% semilogx(noise_levels_corr, dmi_corr_vmax_ratios_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
% yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'No Correlation');
% xlabel('Noise Power'); ylabel('Correlation Coefficient');
% title('2H DMI: Lac/Total Signal Ratio vs Vmax Correlation');
% legend('Location', 'best'); grid on;
% ylim([-1, 1]);
%
% % 13C HP Correlation plots
% subplot(2, 3, 2);
% semilogx(noise_levels_corr, hp_corr_kpb_ratios_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
% hold on;
% semilogx(noise_levels_corr, hp_corr_kpb_ratios_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
% yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'No Correlation');
% xlabel('Noise Power'); ylabel('Correlation Coefficient');
% title('13C HP: Bic/Pyr Ratio vs kPB Correlation');
% legend('Location', 'best'); grid on;
% ylim([-1, 1]);
%
% % Scatter plot examples at low noise (1%)
% low_noise_idx = find(abs(noise_levels_corr - 0.01) < 0.005, 1);
% if ~isempty(low_noise_idx) && ~isempty(dmi_ratio_data_amares{low_noise_idx})
%     subplot(2, 3, 3);
%     scatter(dmi_ratio_data_amares{low_noise_idx}, dmi_vmax_data_amares{low_noise_idx}, 50, 'b', 'filled', 'DisplayName', 'AMARES');
%     hold on;
%     scatter(dmi_ratio_data_integration{low_noise_idx}, dmi_vmax_data_integration{low_noise_idx}, 50, 'r', 'filled', 'DisplayName', 'Integration');
%     xlabel('Lac/Total Signal Ratio'); ylabel('Fitted Vmax (mM/min)');
%     title('2H DMI: Scatter at 1% Noise');
%     legend('Location', 'best'); grid on;
% end
%
% % Scatter plot examples at low noise (1%) for 13C
% if ~isempty(low_noise_idx) && ~isempty(hp_ratio_data_amares{low_noise_idx})
%     subplot(2, 3, 4);
%     scatter(hp_ratio_data_amares{low_noise_idx}, hp_kpb_data_amares{low_noise_idx}, 50, 'b', 'filled', 'DisplayName', 'AMARES');
%     hold on;
%     scatter(hp_ratio_data_integration{low_noise_idx}, hp_kpb_data_integration{low_noise_idx}, 50, 'r', 'filled', 'DisplayName', 'Integration');
%     xlabel('Bic/Pyr Ratio'); ylabel('Fitted kPB (1/s)');
%     title('13C HP: Scatter at 1% Noise');
%     legend('Location', 'best'); grid on;
% end
%
% % Comparison of correlation strength between methods
% subplot(2, 3, 5);
% semilogx(noise_levels_corr, abs(dmi_corr_vmax_ratios_amares), 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI AMARES');
% hold on;
% semilogx(noise_levels_corr, abs(dmi_corr_vmax_ratios_integration), 'b--', 'LineWidth', 2, 'DisplayName', '2H DMI Integration');
% semilogx(noise_levels_corr, abs(hp_corr_kpb_ratios_amares), 'r-', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
% semilogx(noise_levels_corr, abs(hp_corr_kpb_ratios_integration), 'r--', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
% xlabel('Noise Power'); ylabel('|Correlation Coefficient|');
% title('Correlation Strength Comparison');
% legend('Location', 'best'); grid on;
%
% % Cross-method correlation comparison
% subplot(2, 3, 6);
% dmi_corr_ratio = abs(dmi_corr_vmax_ratios_integration) ./ (abs(dmi_corr_vmax_ratios_amares) + eps);
% hp_corr_ratio = abs(hp_corr_kpb_ratios_integration) ./ (abs(hp_corr_kpb_ratios_amares) + eps);
% semilogx(noise_levels_corr, dmi_corr_ratio, 'b-', 'LineWidth', 2, 'DisplayName', '2H DMI (Int/AMARES)');
% hold on;
% semilogx(noise_levels_corr, hp_corr_ratio, 'r-', 'LineWidth', 2, 'DisplayName', '13C HP (Int/AMARES)');
% yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
% xlabel('Noise Power'); ylabel('Integration/AMARES Correlation Ratio');
% title('Relative Correlation Performance');
% legend('Location', 'best'); grid on;
%
% sgtitle('Spectral Ratio vs Kinetic Parameter Correlation Analysis', 'FontSize', 16);
%
% % Store correlation results
% correlation_results_dmi.noise_levels = noise_levels_corr;
% correlation_results_dmi.corr_vmax_amares = dmi_corr_vmax_ratios_amares;
% correlation_results_dmi.corr_vmax_integration = dmi_corr_vmax_ratios_integration;
% correlation_results_dmi.corr_kg_amares = dmi_corr_kg_ratios_amares;
% correlation_results_dmi.corr_kg_integration = dmi_corr_kg_ratios_integration;
% correlation_results_dmi.ratio_data_amares = dmi_ratio_data_amares;
% correlation_results_dmi.ratio_data_integration = dmi_ratio_data_integration;
% correlation_results_dmi.vmax_data_amares = dmi_vmax_data_amares;
% correlation_results_dmi.vmax_data_integration = dmi_vmax_data_integration;
% correlation_results_dmi.kg_data_amares = dmi_kg_data_amares;
% correlation_results_dmi.kg_data_integration = dmi_kg_data_integration;
%
% correlation_results_hp.noise_levels = noise_levels_corr;
% correlation_results_hp.corr_kpb_amares = hp_corr_kpb_ratios_amares;
% correlation_results_hp.corr_kpb_integration = hp_corr_kpb_ratios_integration;
% correlation_results_hp.ratio_data_amares = hp_ratio_data_amares;
% correlation_results_hp.ratio_data_integration = hp_ratio_data_integration;
% correlation_results_hp.kpb_data_amares = hp_kpb_data_amares;
% correlation_results_hp.kpb_data_integration = hp_kpb_data_integration;
%
% % Create K_g correlation figure
% figure('Position', [50 50 1200 800]); set(gcf, 'color', 'w');
%
% % DMI K_g Correlation plot
% subplot(2, 2, 1);
% semilogx(noise_levels_corr, dmi_corr_kg_ratios_amares, 'b-', 'LineWidth', 2, 'DisplayName', 'AMARES');
% hold on;
% semilogx(noise_levels_corr, dmi_corr_kg_ratios_integration, 'r-', 'LineWidth', 2, 'DisplayName', 'Integration');
% yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'No Correlation');
% xlabel('Noise Power'); ylabel('Correlation Coefficient');
% title('2H DMI: Lac/Total Signal Ratio vs K_g Correlation');
% legend('Location', 'best'); grid on;
% ylim([-1, 1]);
%
% % Scatter plot examples at low noise (1%) for K_g
% subplot(2, 2, 2);
% if ~isempty(low_noise_idx) && ~isempty(dmi_ratio_data_amares{low_noise_idx})
%     scatter(dmi_ratio_data_amares{low_noise_idx}, dmi_kg_data_amares{low_noise_idx}, 50, 'b', 'filled', 'DisplayName', 'AMARES');
%     hold on;
%     scatter(dmi_ratio_data_integration{low_noise_idx}, dmi_kg_data_integration{low_noise_idx}, 50, 'r', 'filled', 'DisplayName', 'Integration');
%     xlabel('Lac/Total Signal Ratio'); ylabel('Fitted K_g (min^{-1})');
%     title('2H DMI: Scatter at 1% Noise');
%     legend('Location', 'best'); grid on;
% end
%
% % Comparison of K_g correlation strength
% subplot(2, 2, 3);
% semilogx(noise_levels_corr, abs(dmi_corr_kg_ratios_amares), 'b-', 'LineWidth', 2, 'DisplayName', 'K_g AMARES');
% hold on;
% semilogx(noise_levels_corr, abs(dmi_corr_kg_ratios_integration), 'b--', 'LineWidth', 2, 'DisplayName', 'K_g Integration');
% semilogx(noise_levels_corr, abs(dmi_corr_vmax_ratios_amares), 'r-', 'LineWidth', 2, 'DisplayName', 'Vmax AMARES');
% semilogx(noise_levels_corr, abs(dmi_corr_vmax_ratios_integration), 'r--', 'LineWidth', 2, 'DisplayName', 'Vmax Integration');
% xlabel('Noise Power'); ylabel('|Correlation Coefficient|');
% title('K_g vs Vmax Correlation Strength');
% legend('Location', 'best'); grid on;
%
% % K_g correlation method comparison
% subplot(2, 2, 4);
% kg_corr_ratio = abs(dmi_corr_kg_ratios_integration) ./ (abs(dmi_corr_kg_ratios_amares) + eps);
% semilogx(noise_levels_corr, kg_corr_ratio, 'b-', 'LineWidth', 2, 'DisplayName', 'K_g (Int/AMARES)');
% hold on;
% semilogx(noise_levels_corr, dmi_corr_ratio, 'r-', 'LineWidth', 2, 'DisplayName', 'Vmax (Int/AMARES)');
% yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
% xlabel('Noise Power'); ylabel('Integration/AMARES Correlation Ratio');
% title('K_g vs Vmax: Relative Correlation Performance');
% legend('Location', 'best'); grid on;
%
% sgtitle('2H DMI: K_g Correlation Analysis', 'FontSize', 16);
%
% fprintf('Correlation analysis complete!\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rate Constant Variation Analysis at Fixed Noise Levels
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting rate constant variation analysis at fixed noise levels...\n');

% Fixed noise levels (not percentages!)
fixed_noise_levels = [1e-3, 0.01, 0.1, 1, 10];
n_fixed_noise = length(fixed_noise_levels);

% Parameter variation ranges
n_param_variations = 20;  % Number of different parameter values to test
vmax_range = linspace(dmi_params_true.Vmax * 0.3, dmi_params_true.Vmax * 1.7, n_param_variations);
kpb_range = linspace(k_true_hp * 0.3, k_true_hp * 1.7, n_param_variations);

% Number of Monte Carlo iterations per parameter combination
n_mc_param = 10;

% Initialize storage for parameter variation results
param_var_results = struct();

% Storage for correlations at each noise level
param_corr_dmi_amares = zeros(n_fixed_noise, 1);
param_corr_dmi_integration = zeros(n_fixed_noise, 1);
param_corr_hp_amares = zeros(n_fixed_noise, 1);
param_corr_hp_integration = zeros(n_fixed_noise, 1);

% Storage for all data points
param_var_vmax_true_all = cell(n_fixed_noise, 1);
param_var_vmax_fitted_amares_all = cell(n_fixed_noise, 1);
param_var_vmax_fitted_integration_all = cell(n_fixed_noise, 1);
param_var_kpb_true_all = cell(n_fixed_noise, 1);
param_var_kpb_fitted_amares_all = cell(n_fixed_noise, 1);
param_var_kpb_fitted_integration_all = cell(n_fixed_noise, 1);

% Storage for ratios of interest
param_var_dmi_ratios_amares_all = cell(n_fixed_noise, 1);
param_var_dmi_ratios_integration_all = cell(n_fixed_noise, 1);
param_var_hp_ratios_amares_all = cell(n_fixed_noise, 1);
param_var_hp_ratios_integration_all = cell(n_fixed_noise, 1);

fprintf('Testing %d noise levels with %d parameter variations each\n', n_fixed_noise, n_param_variations);

%% Parameter Variation Analysis Loop
for noise_idx = 1:n_fixed_noise
    current_noise = fixed_noise_levels(noise_idx);
    fprintf('\n********************************************\n');
    fprintf('Parameter Variation: Processing Noise Level %.3f (%d/%d)\n', current_noise, noise_idx, n_fixed_noise);
    fprintf('********************************************\n');

    % Temporary storage for this noise level
    vmax_true_temp = [];
    vmax_fitted_amares_temp = [];
    vmax_fitted_integration_temp = [];
    kpb_true_temp = [];
    kpb_fitted_amares_temp = [];
    kpb_fitted_integration_temp = [];

    % Temporary storage for ratios
    dmi_ratios_amares_temp = [];
    dmi_ratios_integration_temp = [];
    hp_ratios_amares_temp = [];
    hp_ratios_integration_temp = [];

    %% DMI Parameter Variation
    fprintf('DMI parameter variation...\n');
    for vmax_idx = 1:n_param_variations
        current_vmax = vmax_range(vmax_idx);

        for mc_idx = 1:n_mc_param
            try
                % Generate synthetic 2H spectra with varied Vmax
                fake_spectra_dmi_param = generate_fake_fids('2H', 'noise_level', current_noise, ...
                    'a1', dmi_a1, 'kp', dmi_kp, ...
                    'v', dmi_v, 'kg', dmi_kg, ...
                    'Vmax', current_vmax, 'kl', dmi_kl);

                % Quantify using AMARES
                quantified_results_dmi_amares_param = quantify_spectra_pipline(fake_spectra_dmi_param, 'amares');

                % Quantify using Integration
                quantified_results_dmi_integration_param = quantify_spectra_pipline(fake_spectra_dmi_param, 'integration');

                % Extract glucose and lactate for kinetic fitting
                glucose_idx = find(strcmp(fake_spectra_dmi_param.metabolite_names, 'Glucose'));
                lactate_idx = find(strcmp(fake_spectra_dmi_param.metabolite_names, 'Lactate'));

                if ~isempty(glucose_idx) && ~isempty(lactate_idx)
                    % AMARES kinetic fitting
                    glucose_amares = quantified_results_dmi_amares_param.quantified_amplitudes(glucose_idx, :);
                    lactate_amares = quantified_results_dmi_amares_param.quantified_amplitudes(lactate_idx, :);

                    dmi_params_fixed = struct('f', 0.6, 'Km', 2.8);
                    dmi_params_est = struct('a1', 15.0, 'kp', 0.05, 'v', 0.6, ...
                        'kg', 0.1, 'Vmax', 0.8, 'kl', 0.1);

                    [dmi_params_fit_amares_param, ~, ~, ~] = fit_dmi_kinetics(...
                        glucose_amares', lactate_amares', fake_spectra_dmi_param.time_points/60, ...
                        dmi_params_fixed, dmi_params_est);

                    % Integration kinetic fitting
                    glucose_integration = quantified_results_dmi_integration_param.quantified_amplitudes(glucose_idx, :)*normFacDmi;
                    lactate_integration = quantified_results_dmi_integration_param.quantified_amplitudes(lactate_idx, :)*normFacDmi;

                    [dmi_params_fit_integration_param, ~, ~, ~] = fit_dmi_kinetics(...
                        glucose_integration', lactate_integration', fake_spectra_dmi_param.time_points/60, ...
                        dmi_params_fixed, dmi_params_est);

                    % Store results
                    vmax_true_temp(end+1) = current_vmax;
                    vmax_fitted_amares_temp(end+1) = dmi_params_fit_amares_param.Vmax;
                    vmax_fitted_integration_temp(end+1) = dmi_params_fit_integration_param.Vmax;

                    % Extract and store ratios of interest (Lac/Total 2H)
                    if isfield(quantified_results_dmi_amares_param, 'ratios_of_interest') && ...
                            isfield(quantified_results_dmi_integration_param, 'ratios_of_interest')
                        dmi_ratios_amares_temp(end+1) = quantified_results_dmi_amares_param.ratios_of_interest(2);
                        dmi_ratios_integration_temp(end+1) = quantified_results_dmi_integration_param.ratios_of_interest(2);
                    end
                end

            catch ME
                fprintf('DMI Param Var Error at Noise %.3f, Vmax %.3f, MC %d: %s\n', current_noise, current_vmax, mc_idx, ME.message);
                continue;
            end
        end
    end

    %% 13C HP Parameter Variation
    fprintf('13C HP parameter variation...\n');
    for kpb_idx = 1:n_param_variations
        current_kpb = kpb_range(kpb_idx);

        for mc_idx = 1:n_mc_param
            try
                % Generate synthetic 13C spectra with varied kPB
                fake_spectra_hp_param = generate_fake_fids('13C', 'noise_level', current_noise, ...
                    'k_true', current_kpb, 'r1_pyr', r1_pyr_hp, 'r1_bic', r1_bic_hp, ...
                    'pyr_0', pyr_0_hp, 'FlipAngle', FlipAngle_hp);

                % Quantify using AMARES
                quantified_results_hp_amares_param = quantify_spectra_pipline(fake_spectra_hp_param, 'amares');

                % Quantify using Integration
                quantified_results_hp_integration_param = quantify_spectra_pipline(fake_spectra_hp_param, 'integration');

                % Extract pyruvate and bicarbonate for kinetic fitting
                pyr_idx = find(strcmp(fake_spectra_hp_param.metabolite_names, '1-13C Pyruvate'));
                bic_idx = find(strcmp(fake_spectra_hp_param.metabolite_names, '1-13C Bicarbonate'));

                if ~isempty(pyr_idx) && ~isempty(bic_idx)
                    % AMARES kinetic fitting
                    pyr_amares = quantified_results_hp_amares_param.quantified_amplitudes(pyr_idx, :);
                    bic_amares = quantified_results_hp_amares_param.quantified_amplitudes(bic_idx, :);

                    n_timepoints = length(fake_spectra_hp_param.time_points);
                    TR = fake_spectra_hp_param.time_points(2) - fake_spectra_hp_param.time_points(1);

                    signals_amares = zeros(1, 2, n_timepoints);
                    signals_amares(1, 1, :) = pyr_amares*normFacCarbon;
                    signals_amares(1, 2, :) = bic_amares*normFacCarbon;

                    flip_angles_pyr = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                    flip_angles_bic = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                    flip_angles = [flip_angles_pyr; flip_angles_bic];

                    params_fixed = struct('R1P', r1_pyr_hp);
                    params_est = struct('kPL', 0.02, 'R1L', 1/25, 'B0_start', 0.0,'kPL_lb',0,'R1L_lb',0,'R1P_lb',0);

                    [params_fit_amares_param, ~, ~, ~] = fit_pyr_kinetics(signals_amares, TR, ...
                        flip_angles, params_fixed, params_est);

                    % Integration kinetic fitting
                    pyr_integration = quantified_results_hp_integration_param.quantified_amplitudes(pyr_idx, :)*normFacCarbon;
                    bic_integration = quantified_results_hp_integration_param.quantified_amplitudes(bic_idx, :)*normFacCarbon;

                    signals_integration = zeros(1, 2, n_timepoints);
                    signals_integration(1, 1, :) = pyr_integration*normFacCarbon;
                    signals_integration(1, 2, :) = bic_integration*normFacCarbon;

                    [params_fit_integration_param, ~, ~, ~] = fit_pyr_kinetics(signals_integration, TR, ...
                        flip_angles, params_fixed, params_est);

                    % Store results
                    kpb_true_temp(end+1) = current_kpb;
                    kpb_fitted_amares_temp(end+1) = params_fit_amares_param.kPL;
                    kpb_fitted_integration_temp(end+1) = params_fit_integration_param.kPL;

                    % Extract and store ratios of interest (Bic/Pyr)
                    if isfield(quantified_results_hp_amares_param, 'ratios_of_interest') && ...
                            isfield(quantified_results_hp_integration_param, 'ratios_of_interest')
                        hp_ratios_amares_temp(end+1) = quantified_results_hp_amares_param.ratios_of_interest(2);
                        hp_ratios_integration_temp(end+1) = quantified_results_hp_integration_param.ratios_of_interest(2);
                    end
                end

            catch ME
                fprintf('13C HP Param Var Error at Noise %.3f, kPB %.3f, MC %d: %s\n', current_noise, current_kpb, mc_idx, ME.message);
                % keyboard;  % uncomment to debug
                continue;
            end
        end
    end

    % Store data for this noise level
    param_var_vmax_true_all{noise_idx} = vmax_true_temp;
    param_var_vmax_fitted_amares_all{noise_idx} = vmax_fitted_amares_temp;
    param_var_vmax_fitted_integration_all{noise_idx} = vmax_fitted_integration_temp;
    param_var_kpb_true_all{noise_idx} = kpb_true_temp;
    param_var_kpb_fitted_amares_all{noise_idx} = kpb_fitted_amares_temp;
    param_var_kpb_fitted_integration_all{noise_idx} = kpb_fitted_integration_temp;

    % Store ratio data for this noise level
    param_var_dmi_ratios_amares_all{noise_idx} = dmi_ratios_amares_temp;
    param_var_dmi_ratios_integration_all{noise_idx} = dmi_ratios_integration_temp;
    param_var_hp_ratios_amares_all{noise_idx} = hp_ratios_amares_temp;
    param_var_hp_ratios_integration_all{noise_idx} = hp_ratios_integration_temp;

    % Calculate correlations for this noise level
    if length(vmax_true_temp) > 5 && length(vmax_fitted_amares_temp) > 5
        corr_matrix = corrcoef(vmax_true_temp, vmax_fitted_amares_temp, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            param_corr_dmi_amares(noise_idx) = corr_matrix(1, 2);
        end

        corr_matrix = corrcoef(vmax_true_temp, vmax_fitted_integration_temp, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            param_corr_dmi_integration(noise_idx) = corr_matrix(1, 2);
        end
    end

    if length(kpb_true_temp) > 5 && length(kpb_fitted_amares_temp) > 5
        corr_matrix = corrcoef(kpb_true_temp, kpb_fitted_amares_temp, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            param_corr_hp_amares(noise_idx) = corr_matrix(1, 2);
        end

        corr_matrix = corrcoef(kpb_true_temp, kpb_fitted_integration_temp, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            param_corr_hp_integration(noise_idx) = corr_matrix(1, 2);
        end
    end

    fprintf('Correlations at noise %.3f: DMI AMARES=%.3f, DMI Int=%.3f, HP AMARES=%.3f, HP Int=%.3f\n', ...
        current_noise, param_corr_dmi_amares(noise_idx), param_corr_dmi_integration(noise_idx), ...
        param_corr_hp_amares(noise_idx), param_corr_hp_integration(noise_idx));
end

%% Plot Parameter Variation Results
fprintf('\nGenerating parameter variation plots...\n');

% Create individual correlation plots for each noise level
for noise_idx = 1:n_fixed_noise
    current_noise = fixed_noise_levels(noise_idx);

    % Create figure for this noise level
    figure('Position', [100 + noise_idx*50, 100 + noise_idx*30, 1400, 700]);
    set(gcf, 'color', 'w');

    % DMI scatter plots
    subplot(2, 4, 1);
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_amares_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_vmax_fitted_amares_all{noise_idx}, 50, 'b', 'filled', 'DisplayName', 'AMARES');
        hold on;
        plot([min(param_var_vmax_true_all{noise_idx}), max(param_var_vmax_true_all{noise_idx})], ...
            [min(param_var_vmax_true_all{noise_idx}), max(param_var_vmax_true_all{noise_idx})], ...
            'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
        xlabel('True Vmax (mM/min)'); ylabel('Fitted Vmax (mM/min)');
        title(sprintf('2H DMI AMARES: Noise = %.3f', current_noise));
        legend('Location', 'best'); grid on; axis equal;
        text(0.05, 0.95, sprintf('r = %.3f', param_corr_dmi_amares(noise_idx)), ...
            'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');
    end

    subplot(2, 4, 2);
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_integration_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_vmax_fitted_integration_all{noise_idx}, 50, 'r', 'filled', 'DisplayName', 'Integration');
        hold on;
        plot([min(param_var_vmax_true_all{noise_idx}), max(param_var_vmax_true_all{noise_idx})], ...
            [min(param_var_vmax_true_all{noise_idx}), max(param_var_vmax_true_all{noise_idx})], ...
            'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
        xlabel('True Vmax (mM/min)'); ylabel('Fitted Vmax (mM/min)');
        title(sprintf('2H DMI Integration: Noise = %.3f', current_noise));
        legend('Location', 'best'); grid on; axis equal;
        text(0.05, 0.95, sprintf('r = %.3f', param_corr_dmi_integration(noise_idx)), ...
            'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');
    end

    % 13C HP scatter plots
    subplot(2, 4, 3);
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_amares_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_kpb_fitted_amares_all{noise_idx}, 50, 'b', 'filled', 'DisplayName', 'AMARES');
        hold on;
        plot([min(param_var_kpb_true_all{noise_idx}), max(param_var_kpb_true_all{noise_idx})], ...
            [min(param_var_kpb_true_all{noise_idx}), max(param_var_kpb_true_all{noise_idx})], ...
            'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
        xlabel('True kPB (1/s)'); ylabel('Fitted kPB (1/s)');
        title(sprintf('13C HP AMARES: Noise = %.3f', current_noise));
        legend('Location', 'best'); grid on; axis equal;
        text(0.05, 0.95, sprintf('r = %.3f', param_corr_hp_amares(noise_idx)), ...
            'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');
    end

    subplot(2, 4, 4);
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_integration_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_kpb_fitted_integration_all{noise_idx}, 50, 'r', 'filled', 'DisplayName', 'Integration');
        hold on;
        plot([min(param_var_kpb_true_all{noise_idx}), max(param_var_kpb_true_all{noise_idx})], ...
            [min(param_var_kpb_true_all{noise_idx}), max(param_var_kpb_true_all{noise_idx})], ...
            'k--', 'LineWidth', 1, 'DisplayName', 'Perfect Agreement');
        xlabel('True kPB (1/s)'); ylabel('Fitted kPB (1/s)');
        title(sprintf('13C HP Integration: Noise = %.3f', current_noise));
        legend('Location', 'best'); grid on; axis equal;
        text(0.05, 0.95, sprintf('r = %.3f', param_corr_hp_integration(noise_idx)), ...
            'Units', 'normalized', 'FontSize', 12, 'BackgroundColor', 'white');
    end

    % Bland-Altman plots for DMI
    subplot(2, 4, 5);
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_amares_all{noise_idx})
        mean_vals = (param_var_vmax_true_all{noise_idx} + param_var_vmax_fitted_amares_all{noise_idx}) / 2;
        diff_vals = param_var_vmax_fitted_amares_all{noise_idx} - param_var_vmax_true_all{noise_idx};
        scatter(mean_vals, diff_vals, 50, 'b', 'filled');
        hold on;
        yline(mean(diff_vals), 'r--', 'LineWidth', 2);
        yline(mean(diff_vals) + 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        yline(mean(diff_vals) - 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        xlabel('Mean Vmax'); ylabel('Difference (Fitted - True)');
        title('DMI AMARES: Bland-Altman'); grid on;
    end

    subplot(2, 4, 6);
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_integration_all{noise_idx})
        mean_vals = (param_var_vmax_true_all{noise_idx} + param_var_vmax_fitted_integration_all{noise_idx}) / 2;
        diff_vals = param_var_vmax_fitted_integration_all{noise_idx} - param_var_vmax_true_all{noise_idx};
        scatter(mean_vals, diff_vals, 50, 'r', 'filled');
        hold on;
        yline(mean(diff_vals), 'r--', 'LineWidth', 2);
        yline(mean(diff_vals) + 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        yline(mean(diff_vals) - 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        xlabel('Mean Vmax'); ylabel('Difference (Fitted - True)');
        title('DMI Integration: Bland-Altman'); grid on;
    end

    % Bland-Altman plots for 13C HP
    subplot(2, 4, 7);
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_amares_all{noise_idx})
        mean_vals = (param_var_kpb_true_all{noise_idx} + param_var_kpb_fitted_amares_all{noise_idx}) / 2;
        diff_vals = param_var_kpb_fitted_amares_all{noise_idx} - param_var_kpb_true_all{noise_idx};
        scatter(mean_vals, diff_vals, 50, 'b', 'filled');
        hold on;
        yline(mean(diff_vals), 'r--', 'LineWidth', 2);
        yline(mean(diff_vals) + 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        yline(mean(diff_vals) - 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        xlabel('Mean kPB'); ylabel('Difference (Fitted - True)');
        title('13C HP AMARES: Bland-Altman'); grid on;
    end

    subplot(2, 4, 8);
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_integration_all{noise_idx})
        mean_vals = (param_var_kpb_true_all{noise_idx} + param_var_kpb_fitted_integration_all{noise_idx}) / 2;
        diff_vals = param_var_kpb_fitted_integration_all{noise_idx} - param_var_kpb_true_all{noise_idx};
        scatter(mean_vals, diff_vals, 50, 'r', 'filled');
        hold on;
        yline(mean(diff_vals), 'r--', 'LineWidth', 2);
        yline(mean(diff_vals) + 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        yline(mean(diff_vals) - 1.96*std(diff_vals), 'g--', 'LineWidth', 1.5);
        xlabel('Mean kPB'); ylabel('Difference (Fitted - True)');
        title('13C HP Integration: Bland-Altman'); grid on;
    end

    sgtitle(sprintf('Parameter Recovery Analysis: Noise Level = %.3f', current_noise), 'FontSize', 16);
end

% Create comprehensive heatmap figure
figure('Position', [50 50 1400 1000]); set(gcf, 'color', 'w');

% Correlation heatmap
subplot(2, 3, 1);
heatmap_data = [param_corr_dmi_amares'; param_corr_dmi_integration'; param_corr_hp_amares'; param_corr_hp_integration'];
imagesc(heatmap_data);
colorbar;
colormap(gca, 'parula');
caxis([-1, 1]);
set(gca, 'XTick', 1:n_fixed_noise, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), fixed_noise_levels, 'UniformOutput', false));
set(gca, 'YTick', 1:4, 'YTickLabel', {'DMI AMARES', 'DMI Integration', '13C HP AMARES', '13C HP Integration'});
xlabel('Noise Power'); ylabel('Method');
title('True vs Fitted Parameter Correlations');

% % Add correlation values as text
% for i = 1:4
%     for j = 1:n_fixed_noise
%         text(j, i, sprintf('%.2f', heatmap_data(i, j)), 'HorizontalAlignment', 'center', ...
%              'Color', abs(heatmap_data(i, j)) > 0.5 ? 'white' : 'black', 'FontSize', 10);
%     end
% end

% Correlation comparison across noise levels
subplot(2, 3, 2);
plot(fixed_noise_levels, param_corr_dmi_amares, 'b-o', 'LineWidth', 2, 'DisplayName', 'DMI AMARES');
hold on;
plot(fixed_noise_levels, param_corr_dmi_integration, 'b--s', 'LineWidth', 2, 'DisplayName', 'DMI Integration');
plot(fixed_noise_levels, param_corr_hp_amares, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
plot(fixed_noise_levels, param_corr_hp_integration, 'r--s', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Correlation Coefficient');
title('Parameter Recovery vs Noise Level');
legend('Location', 'best'); grid on;
ylim([-1, 1]);

% Method comparison ratios
subplot(2, 3, 3);
dmi_ratio = param_corr_dmi_integration ./ (param_corr_dmi_amares + eps);
hp_ratio = param_corr_hp_integration ./ (param_corr_hp_amares + eps);
plot(fixed_noise_levels, dmi_ratio, 'b-o', 'LineWidth', 2, 'DisplayName', 'DMI (Int/AMARES)');
hold on;
plot(fixed_noise_levels, hp_ratio, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP (Int/AMARES)');
yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Integration/AMARES Ratio');
title('Relative Method Performance');
legend('Location', 'best'); grid on;

% Summary statistics heatmap - mean absolute correlation
subplot(2, 3, 4);
mean_corr_by_method = [mean(abs(param_corr_dmi_amares)); mean(abs(param_corr_dmi_integration));
    mean(abs(param_corr_hp_amares)); mean(abs(param_corr_hp_integration))];
bar(mean_corr_by_method, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', {'DMI AMARES', 'DMI Int', '13C HP AMARES', '13C HP Int'});
ylabel('Mean |Correlation|');
title('Overall Method Performance');
grid on;

% Standard deviation of correlations
subplot(2, 3, 5);
std_corr_by_method = [std(param_corr_dmi_amares); std(param_corr_dmi_integration);
    std(param_corr_hp_amares); std(param_corr_hp_integration)];
bar(std_corr_by_method, 'FaceColor', [0.8 0.4 0.2]);
set(gca, 'XTickLabel', {'DMI AMARES', 'DMI Int', '13C HP AMARES', '13C HP Int'});
ylabel('Std Dev of Correlations');
title('Method Consistency Across Noise Levels');
grid on;

% Nuclear comparison
subplot(2, 3, 6);
dmi_mean_corr = mean([param_corr_dmi_amares, param_corr_dmi_integration], 2);
hp_mean_corr = mean([param_corr_hp_amares, param_corr_hp_integration], 2);
plot(fixed_noise_levels, dmi_mean_corr, 'b-o', 'LineWidth', 2, 'DisplayName', '2H DMI (mean)');
hold on;
plot(fixed_noise_levels, hp_mean_corr, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP (mean)');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Mean Correlation');
title('Nuclear Method Comparison');
legend('Location', 'best'); grid on;

sgtitle('Comprehensive Parameter Recovery Analysis: Rate Constant Variation', 'FontSize', 16);

% Store parameter variation results
param_var_results.fixed_noise_levels = fixed_noise_levels;
param_var_results.vmax_range = vmax_range;
param_var_results.kpb_range = kpb_range;
param_var_results.param_corr_dmi_amares = param_corr_dmi_amares;
param_var_results.param_corr_dmi_integration = param_corr_dmi_integration;
param_var_results.param_corr_hp_amares = param_corr_hp_amares;
param_var_results.param_corr_hp_integration = param_corr_hp_integration;
param_var_results.vmax_true_all = param_var_vmax_true_all;
param_var_results.vmax_fitted_amares_all = param_var_vmax_fitted_amares_all;
param_var_results.vmax_fitted_integration_all = param_var_vmax_fitted_integration_all;
param_var_results.kpb_true_all = param_var_kpb_true_all;
param_var_results.kpb_fitted_amares_all = param_var_kpb_fitted_amares_all;
param_var_results.kpb_fitted_integration_all = param_var_kpb_fitted_integration_all;

fprintf('Parameter variation analysis complete!\n');

%% Plot Ratio Correlation Analysis
fprintf('\nGenerating ratio correlation plots...\n');

markerSizeForScatter = 30;
alphaForScatter = 0.2;
% Create separate figure for SNR scatter plots (all noise levels combined)
figure('Position', [50 50 1200 900]); set(gcf, 'color', 'w');
colors = lines(n_fixed_noise);

% DMI: Vmax vs Lac/Total Signal ratio correlation plots for each noise level
subplot(2, 2, 1);
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_amares_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_amares_all{noise_idx}, ...
            markerSizeForScatter,colors(noise_idx,:), 'filled', 'DisplayName', sprintf('Noise=%.3f', fixed_noise_levels(noise_idx)));
        hold on;
    end
end
xlabel('True Vmax (mM/min)'); ylabel('Lac/Total Signal Ratio (AMARES)');
title('2H DMI: True Vmax vs Lac/Total Signal Ratio (AMARES)');
legend('Location', 'best'); grid on;

subplot(2, 2, 2);
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_integration_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_integration_all{noise_idx}, ...
            markerSizeForScatter,colors(noise_idx,:), 'filled', 'DisplayName', sprintf('Noise=%.3f', fixed_noise_levels(noise_idx)));
        hold on;
    end
end
xlabel('True Vmax (mM/min)'); ylabel('Lac/Total Signal Ratio (Integration)');
title('2H DMI: True Vmax vs Lac/Total Signal Ratio (Integration)');
legend('Location', 'best'); grid on;

% 13C HP: kPB vs Bic/Pyr ratio correlation plots for each noise level
subplot(2, 2, 3);
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_amares_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_amares_all{noise_idx}, ...
            markerSizeForScatter,colors(noise_idx,:), 'filled', 'DisplayName', sprintf('Noise=%.3f', fixed_noise_levels(noise_idx)));
        hold on;
    end
end
xlabel('True kPB (1/s)'); ylabel('Bic/Pyr Ratio (AMARES)');
title('13C HP: True kPB vs Bic/Pyr Ratio (AMARES)');
legend('Location', 'best'); grid on;

subplot(2, 2, 4);
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_integration_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_integration_all{noise_idx}, ...
            markerSizeForScatter,colors(noise_idx,:), 'filled', 'DisplayName', sprintf('Noise=%.3f', fixed_noise_levels(noise_idx)));
        hold on;
    end
end
xlabel('True kPB (1/s)'); ylabel('Bic/Pyr Ratio (Integration)');
title('13C HP: True kPB vs Bic/Pyr Ratio (Integration)');
legend('Location', 'best'); grid on;

sgtitle('Spectral Ratios vs True Kinetic Parameters: All Noise Levels Combined', 'FontSize', 16);

% Create figure for individual noise level ratio correlations (4x2 layout)
figure('Position', [100 50 1400 1400]); set(gcf, 'color', 'w');

% Individual noise level plots - 4x2 layout (DMI left column, 13C HP right column)
for noise_idx = 1:min(4, n_fixed_noise)  % Show first 4 noise levels
    current_noise = fixed_noise_levels(noise_idx);

    % DMI Vmax vs ratio (left column)
    subplot(4, 2, 2*noise_idx - 1);

    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_amares_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_amares_all{noise_idx}, ...
            markerSizeForScatter,'b', 'filled', 'DisplayName', 'AMARES');
        hold on;

        % Calculate and display correlation (coloured blue for AMARES)
        if length(param_var_vmax_true_all{noise_idx}) > 2
            corr_coef = corrcoef(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_amares_all{noise_idx}, 'rows', 'complete');
            if size(corr_coef, 1) == 2
                text(0.05, 0.95, sprintf('r = %.3f', corr_coef(1,2)), 'Units', 'normalized', ...
                    'FontSize', 12, 'BackgroundColor', 'white', 'Color', 'blue');
            end
        end
    end

    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_integration_all{noise_idx})
        scatter(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_integration_all{noise_idx}, ...
            markerSizeForScatter,'r', 'filled', 'DisplayName', 'Integration');

        % Calculate and display correlation (coloured red for integration)
        if length(param_var_vmax_true_all{noise_idx}) > 2
            corr_coef = corrcoef(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_integration_all{noise_idx}, 'rows', 'complete');
            if size(corr_coef, 1) == 2
                text(0.05, 0.85, sprintf('r = %.3f', corr_coef(1,2)), 'Units', 'normalized', ...
                    'FontSize', 12, 'BackgroundColor', 'white', 'Color', 'red');
            end
        end
    end

    xlabel('True Vmax (mM/min)'); ylabel('Lac/Total Signal Ratio');
    title(sprintf('2H DMI: Noise = %.3f', current_noise));
    legend('Location', 'best'); grid on;

    % 13C HP kPB vs ratio (right column)
    subplot(4, 2, 2*noise_idx);

    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_amares_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_amares_all{noise_idx}, ...
            markerSizeForScatter,'b', 'filled', 'DisplayName', 'AMARES');
        hold on;

        % Calculate and display correlation (coloured blue for AMARES)
        if length(param_var_kpb_true_all{noise_idx}) > 2
            corr_coef = corrcoef(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_amares_all{noise_idx}, 'rows', 'complete');
            if size(corr_coef, 1) == 2
                text(0.05, 0.95, sprintf('r = %.3f', corr_coef(1,2)), 'Units', 'normalized', ...
                    'FontSize', 12, 'BackgroundColor', 'white', 'Color', 'blue');
            end
        end
    end

    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_integration_all{noise_idx})
        scatter(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_integration_all{noise_idx}, ...
            markerSizeForScatter,'r', 'filled', 'DisplayName', 'Integration');

        % Calculate and display correlation (coloured red for integration)
        if length(param_var_kpb_true_all{noise_idx}) > 2
            corr_coef = corrcoef(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_integration_all{noise_idx}, 'rows', 'complete');
            if size(corr_coef, 1) == 2
                text(0.05, 0.85, sprintf('r = %.3f', corr_coef(1,2)), 'Units', 'normalized', ...
                    'FontSize', 12, 'BackgroundColor', 'white', 'Color', 'red');
            end
        end
    end

    xlabel('True kPB (1/s)'); ylabel('Bic/Pyr Ratio');
    title(sprintf('13C HP: Noise = %.3f', current_noise));
    legend('Location', 'best'); grid on;
end

sgtitle('Spectral Ratios vs True Kinetic Parameters: Individual Noise Levels', 'FontSize', 16);
%% Create comprehensive fitted vs true parameter figure across noise levels (3D visualization)
fprintf('\nGenerating fitted vs true parameter plots across noise levels (3D)...\n');

figure('Position', [150 50 1600 1000]); set(gcf, 'color', 'w');

% Prepare data for 3D scatter plots
all_vmax_true = [param_var_vmax_true_all{:}];
all_kpb_true = [param_var_kpb_true_all{:}];

% DMI: Fitted Vmax vs True Vmax (AMARES) - 3D scatter
subplot(2, 2, 1);
hold on;
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_amares_all{noise_idx})
        % Create vectors with noise level repeated
        noise_vec = fixed_noise_levels(noise_idx) * ones(size(param_var_vmax_true_all{noise_idx}));
        scatter3(param_var_vmax_true_all{noise_idx}, param_var_vmax_fitted_amares_all{noise_idx}, noise_vec, ...
            40, noise_vec, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end
% Add identity plane
if ~isempty(all_vmax_true)
    vmax_range = linspace(min(all_vmax_true), max(all_vmax_true), 20);
    noise_range = linspace(min(fixed_noise_levels), max(fixed_noise_levels), 20);
    [V_grid, N_grid] = meshgrid(vmax_range, noise_range);
    surf(V_grid, V_grid, N_grid, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('True Vmax (mM/min)'); ylabel('Fitted Vmax (mM/min)'); zlabel('Noise Power');
title('2H DMI: AMARES Quantification');
grid on; view(45, 30);
set(gca, 'ZScale', 'log');
colormap(gca, 'parula'); colorbar;
caxis([min(fixed_noise_levels), max(fixed_noise_levels)]);

% DMI: Fitted Vmax vs True Vmax (Integration) - 3D scatter
subplot(2, 2, 2);
hold on;
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_integration_all{noise_idx})
        noise_vec = fixed_noise_levels(noise_idx) * ones(size(param_var_vmax_true_all{noise_idx}));
        scatter3(param_var_vmax_true_all{noise_idx}, param_var_vmax_fitted_integration_all{noise_idx}, noise_vec, ...
            40, noise_vec, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end
% Add identity plane
if ~isempty(all_vmax_true)
    vmax_range = linspace(min(all_vmax_true), max(all_vmax_true), 20);
    noise_range = linspace(min(fixed_noise_levels), max(fixed_noise_levels), 20);
    [V_grid, N_grid] = meshgrid(vmax_range, noise_range);
    surf(V_grid, V_grid, N_grid, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('True Vmax (mM/min)'); ylabel('Fitted Vmax (mM/min)'); zlabel('Noise Power');
title('2H DMI: Integration Quantification');
grid on; view(45, 30);
set(gca, 'ZScale', 'log');
colormap(gca, 'parula'); colorbar;
caxis([min(fixed_noise_levels), max(fixed_noise_levels)]);

% 13C HP: Fitted kPB vs True kPB (AMARES) - 3D scatter
subplot(2, 2, 3);
hold on;
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_amares_all{noise_idx})
        noise_vec = fixed_noise_levels(noise_idx) * ones(size(param_var_kpb_true_all{noise_idx}));
        scatter3(param_var_kpb_true_all{noise_idx}, param_var_kpb_fitted_amares_all{noise_idx}, noise_vec, ...
            40, noise_vec, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end
% Add identity plane
if ~isempty(all_kpb_true)
    kpb_range = linspace(min(all_kpb_true), max(all_kpb_true), 20);
    noise_range = linspace(min(fixed_noise_levels), max(fixed_noise_levels), 20);
    [K_grid, N_grid] = meshgrid(kpb_range, noise_range);
    surf(K_grid, K_grid, N_grid, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('True kPB (1/s)'); ylabel('Fitted kPB (1/s)'); zlabel('Noise Power');
title('13C HP: AMARES Quantification');
grid on; view(45, 30);
set(gca, 'ZScale', 'log');
colormap(gca, 'parula'); colorbar;
caxis([min(fixed_noise_levels), max(fixed_noise_levels)]);

% 13C HP: Fitted kPB vs True kPB (Integration) - 3D scatter
subplot(2, 2, 4);
hold on;
for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_integration_all{noise_idx})
        noise_vec = fixed_noise_levels(noise_idx) * ones(size(param_var_kpb_true_all{noise_idx}));
        scatter3(param_var_kpb_true_all{noise_idx}, param_var_kpb_fitted_integration_all{noise_idx}, noise_vec, ...
            40, noise_vec, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end
% Add identity plane
if ~isempty(all_kpb_true)
    kpb_range = linspace(min(all_kpb_true), max(all_kpb_true), 20);
    noise_range = linspace(min(fixed_noise_levels), max(fixed_noise_levels), 20);
    [K_grid, N_grid] = meshgrid(kpb_range, noise_range);
    surf(K_grid, K_grid, N_grid, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end
xlabel('True kPB (1/s)'); ylabel('Fitted kPB (1/s)'); zlabel('Noise Power');
title('13C HP: Integration Quantification');
grid on; view(45, 30);
set(gca, 'ZScale', 'log');
colormap(gca, 'parula'); colorbar;
caxis([min(fixed_noise_levels), max(fixed_noise_levels)]);

sgtitle('Parameter Estimation Accuracy: Effect of Noise on Quantification Methods (3D)', 'FontSize', 16);

% Create additional figure showing error metrics as a function of noise
figure('Position', [200 50 1600 900]); set(gcf, 'color', 'w');

% Calculate RMSE and bias for each noise level
rmse_dmi_amares = zeros(n_fixed_noise, 1);
rmse_dmi_integration = zeros(n_fixed_noise, 1);
rmse_hp_amares = zeros(n_fixed_noise, 1);
rmse_hp_integration = zeros(n_fixed_noise, 1);
bias_dmi_amares = zeros(n_fixed_noise, 1);
bias_dmi_integration = zeros(n_fixed_noise, 1);
bias_hp_amares = zeros(n_fixed_noise, 1);
bias_hp_integration = zeros(n_fixed_noise, 1);

for noise_idx = 1:n_fixed_noise
    % DMI RMSE and bias
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_amares_all{noise_idx})
        errors = param_var_vmax_fitted_amares_all{noise_idx} - param_var_vmax_true_all{noise_idx};
        rmse_dmi_amares(noise_idx) = sqrt(mean(errors.^2, 'omitnan'));
        bias_dmi_amares(noise_idx) = mean(errors, 'omitnan');
    end
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_vmax_fitted_integration_all{noise_idx})
        errors = param_var_vmax_fitted_integration_all{noise_idx} - param_var_vmax_true_all{noise_idx};
        rmse_dmi_integration(noise_idx) = sqrt(mean(errors.^2, 'omitnan'));
        bias_dmi_integration(noise_idx) = mean(errors, 'omitnan');
    end

    % 13C HP RMSE and bias
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_amares_all{noise_idx})
        errors = param_var_kpb_fitted_amares_all{noise_idx} - param_var_kpb_true_all{noise_idx};
        rmse_hp_amares(noise_idx) = sqrt(mean(errors.^2, 'omitnan'));
        bias_hp_amares(noise_idx) = mean(errors, 'omitnan');
    end
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_kpb_fitted_integration_all{noise_idx})
        errors = param_var_kpb_fitted_integration_all{noise_idx} - param_var_kpb_true_all{noise_idx};
        rmse_hp_integration(noise_idx) = sqrt(mean(errors.^2, 'omitnan'));
        bias_hp_integration(noise_idx) = mean(errors, 'omitnan');
    end
end

% Plot RMSE vs noise level
subplot(2, 3, 1);
loglog(fixed_noise_levels, rmse_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI AMARES');
hold on;
loglog(fixed_noise_levels, rmse_dmi_integration, 'b--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI Integration');
loglog(fixed_noise_levels, rmse_hp_amares, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP AMARES');
loglog(fixed_noise_levels, rmse_hp_integration, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP Integration');
xlabel('Noise Power'); ylabel('RMSE');
title('Root Mean Square Error vs Noise');
legend('Location', 'best'); grid on;

% Plot bias vs noise level
subplot(2, 3, 2);
semilogx(fixed_noise_levels, bias_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI AMARES');
hold on;
semilogx(fixed_noise_levels, bias_dmi_integration, 'b--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI Integration');
semilogx(fixed_noise_levels, bias_hp_amares, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP AMARES');
semilogx(fixed_noise_levels, bias_hp_integration, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP Integration');
yline(0, 'k--', 'LineWidth', 1);
xlabel('Noise Power'); ylabel('Bias (Fitted - True)');
title('Estimation Bias vs Noise');
legend('Location', 'best'); grid on;

% Plot correlation vs noise level
subplot(2, 3, 3);
semilogx(fixed_noise_levels, param_corr_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI AMARES');
hold on;
semilogx(fixed_noise_levels, param_corr_dmi_integration, 'b--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI Integration');
semilogx(fixed_noise_levels, param_corr_hp_amares, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP AMARES');
semilogx(fixed_noise_levels, param_corr_hp_integration, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP Integration');
yline(1, 'k--', 'LineWidth', 1);
xlabel('Noise Power'); ylabel('Correlation (r)');
title('True vs Fitted Correlation vs Noise');
legend('Location', 'best'); grid on;
ylim([0, 1.1]);

% Plot coefficient of variation vs noise
subplot(2, 3, 4);
cv_dmi_amares = zeros(n_fixed_noise, 1);
cv_dmi_integration = zeros(n_fixed_noise, 1);
cv_hp_amares = zeros(n_fixed_noise, 1);
cv_hp_integration = zeros(n_fixed_noise, 1);

for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_fitted_amares_all{noise_idx})
        cv_dmi_amares(noise_idx) = std(param_var_vmax_fitted_amares_all{noise_idx}, 'omitnan') / ...
            mean(param_var_vmax_fitted_amares_all{noise_idx}, 'omitnan');
    end
    if ~isempty(param_var_vmax_fitted_integration_all{noise_idx})
        cv_dmi_integration(noise_idx) = std(param_var_vmax_fitted_integration_all{noise_idx}, 'omitnan') / ...
            mean(param_var_vmax_fitted_integration_all{noise_idx}, 'omitnan');
    end
    if ~isempty(param_var_kpb_fitted_amares_all{noise_idx})
        cv_hp_amares(noise_idx) = std(param_var_kpb_fitted_amares_all{noise_idx}, 'omitnan') / ...
            mean(param_var_kpb_fitted_amares_all{noise_idx}, 'omitnan');
    end
    if ~isempty(param_var_kpb_fitted_integration_all{noise_idx})
        cv_hp_integration(noise_idx) = std(param_var_kpb_fitted_integration_all{noise_idx}, 'omitnan') / ...
            mean(param_var_kpb_fitted_integration_all{noise_idx}, 'omitnan');
    end
end

loglog(fixed_noise_levels, cv_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI AMARES');
hold on;
loglog(fixed_noise_levels, cv_dmi_integration, 'b--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI Integration');
loglog(fixed_noise_levels, cv_hp_amares, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP AMARES');
loglog(fixed_noise_levels, cv_hp_integration, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP Integration');
xlabel('Noise Power'); ylabel('Coefficient of Variation');
title('Precision (CV) vs Noise');
legend('Location', 'best'); grid on;

% Plot relative RMSE vs noise
subplot(2, 3, 5);
rel_rmse_dmi_amares = zeros(n_fixed_noise, 1);
rel_rmse_dmi_integration = zeros(n_fixed_noise, 1);
rel_rmse_hp_amares = zeros(n_fixed_noise, 1);
rel_rmse_hp_integration = zeros(n_fixed_noise, 1);

for noise_idx = 1:n_fixed_noise
    if ~isempty(param_var_vmax_true_all{noise_idx})
        mean_true_vmax = mean(param_var_vmax_true_all{noise_idx}, 'omitnan');
        rel_rmse_dmi_amares(noise_idx) = rmse_dmi_amares(noise_idx) / mean_true_vmax * 100;
        rel_rmse_dmi_integration(noise_idx) = rmse_dmi_integration(noise_idx) / mean_true_vmax * 100;
    end
    if ~isempty(param_var_kpb_true_all{noise_idx})
        mean_true_kpb = mean(param_var_kpb_true_all{noise_idx}, 'omitnan');
        rel_rmse_hp_amares(noise_idx) = rmse_hp_amares(noise_idx) / mean_true_kpb * 100;
        rel_rmse_hp_integration(noise_idx) = rmse_hp_integration(noise_idx) / mean_true_kpb * 100;
    end
end

loglog(fixed_noise_levels, rel_rmse_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI AMARES');
hold on;
loglog(fixed_noise_levels, rel_rmse_dmi_integration, 'b--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI Integration');
loglog(fixed_noise_levels, rel_rmse_hp_amares, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP AMARES');
loglog(fixed_noise_levels, rel_rmse_hp_integration, 'r--s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP Integration');
xlabel('Noise Power'); ylabel('Relative RMSE (%)');
title('Relative Error vs Noise');
legend('Location', 'best'); grid on;

% Method performance comparison
subplot(2, 3, 6);
rmse_ratio_dmi = rmse_dmi_integration ./ (rmse_dmi_amares + eps);
rmse_ratio_hp = rmse_hp_integration ./ (rmse_hp_amares + eps);
semilogx(fixed_noise_levels, rmse_ratio_dmi, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'DMI (Int/AMARES)');
hold on;
semilogx(fixed_noise_levels, rmse_ratio_hp, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'HP (Int/AMARES)');
yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
xlabel('Noise Power'); ylabel('RMSE Ratio');
title('Relative Method Performance: Integration vs AMARES');
legend('Location', 'best'); grid on;

sgtitle('Quantification Accuracy Metrics: Effect of Noise on Low-SNR Data', 'FontSize', 16);

fprintf('Fitted vs true parameter analysis complete!\n');

% Create correlation summary heatmap for ratios
figure('Position', [100 100 1200 800]); set(gcf, 'color', 'w');

% Calculate correlation matrices for each noise level
dmi_ratio_corr_amares = zeros(n_fixed_noise, 1);
dmi_ratio_corr_integration = zeros(n_fixed_noise, 1);
hp_ratio_corr_amares = zeros(n_fixed_noise, 1);
hp_ratio_corr_integration = zeros(n_fixed_noise, 1);

for noise_idx = 1:n_fixed_noise
    % DMI correlations
    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_amares_all{noise_idx}) && ...
            length(param_var_vmax_true_all{noise_idx}) > 2
        corr_matrix = corrcoef(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_amares_all{noise_idx}, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            dmi_ratio_corr_amares(noise_idx) = corr_matrix(1, 2);
        end
    end

    if ~isempty(param_var_vmax_true_all{noise_idx}) && ~isempty(param_var_dmi_ratios_integration_all{noise_idx}) && ...
            length(param_var_vmax_true_all{noise_idx}) > 2
        corr_matrix = corrcoef(param_var_vmax_true_all{noise_idx}, param_var_dmi_ratios_integration_all{noise_idx}, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            dmi_ratio_corr_integration(noise_idx) = corr_matrix(1, 2);
        end
    end

    % 13C HP correlations
    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_amares_all{noise_idx}) && ...
            length(param_var_kpb_true_all{noise_idx}) > 2
        corr_matrix = corrcoef(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_amares_all{noise_idx}, 'rows', 'complete');
        corr_matrix = corrcoef(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_amares_all{noise_idx}, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            hp_ratio_corr_amares(noise_idx) = corr_matrix(1, 2);
        end
    end

    if ~isempty(param_var_kpb_true_all{noise_idx}) && ~isempty(param_var_hp_ratios_integration_all{noise_idx}) && ...
            length(param_var_kpb_true_all{noise_idx}) > 2
        corr_matrix = corrcoef(param_var_kpb_true_all{noise_idx}, param_var_hp_ratios_integration_all{noise_idx}, 'rows', 'complete');
        if size(corr_matrix, 1) == 2
            hp_ratio_corr_integration(noise_idx) = corr_matrix(1, 2);
        end
    end
end

% Ratio correlation heatmap
subplot(2, 3, 1);
ratio_heatmap_data = [dmi_ratio_corr_amares'; dmi_ratio_corr_integration'; hp_ratio_corr_amares'; hp_ratio_corr_integration'];
imagesc(ratio_heatmap_data);
colorbar;
colormap(gca, 'parula');
caxis([-1, 1]);
set(gca, 'XTick', 1:n_fixed_noise, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), fixed_noise_levels, 'UniformOutput', false));
set(gca, 'YTick', 1:4, 'YTickLabel', {'DMI AMARES', 'DMI Integration', '13C HP AMARES', '13C HP Integration'});
xlabel('Noise Power'); ylabel('Method');
title('True Parameter vs Ratio Correlations');



% Correlation trends
subplot(2, 3, 2);
plot(fixed_noise_levels, dmi_ratio_corr_amares, 'b-o', 'LineWidth', 2, 'DisplayName', 'DMI AMARES');
hold on;
plot(fixed_noise_levels, dmi_ratio_corr_integration, 'b--s', 'LineWidth', 2, 'DisplayName', 'DMI Integration');
plot(fixed_noise_levels, hp_ratio_corr_amares, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP AMARES');
plot(fixed_noise_levels, hp_ratio_corr_integration, 'r--s', 'LineWidth', 2, 'DisplayName', '13C HP Integration');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Correlation Coefficient');
title('Ratio-Parameter Correlation vs Noise');
legend('Location', 'best'); grid on;
ylim([-1, 1]);

% Method comparison for ratios
subplot(2, 3, 3);
dmi_ratio_performance = abs(dmi_ratio_corr_integration) ./ (abs(dmi_ratio_corr_amares) + eps);
hp_ratio_performance = abs(hp_ratio_corr_integration) ./ (abs(hp_ratio_corr_amares) + eps);
plot(fixed_noise_levels, dmi_ratio_performance, 'b-o', 'LineWidth', 2, 'DisplayName', 'DMI (Int/AMARES)');
hold on;
plot(fixed_noise_levels, hp_ratio_performance, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP (Int/AMARES)');
yline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Equal Performance');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Integration/AMARES Ratio');
title('Relative Ratio Correlation Performance');
legend('Location', 'best'); grid on;

% Overall ratio correlation performance
subplot(2, 3, 4);
mean_ratio_corr = [mean(abs(dmi_ratio_corr_amares)); mean(abs(dmi_ratio_corr_integration));
    mean(abs(hp_ratio_corr_amares)); mean(abs(hp_ratio_corr_integration))];
bar(mean_ratio_corr, 'FaceColor', [0.6 0.2 0.8]);
set(gca, 'XTickLabel', {'DMI AMARES', 'DMI Int', '13C HP AMARES', '13C HP Int'});
ylabel('Mean |Correlation|');
title('Overall Ratio Correlation Performance');
grid on;

% Ratio correlation consistency
subplot(2, 3, 5);
std_ratio_corr = [std(dmi_ratio_corr_amares); std(dmi_ratio_corr_integration);
    std(hp_ratio_corr_amares); std(hp_ratio_corr_integration)];
bar(std_ratio_corr, 'FaceColor', [0.8 0.6 0.2]);
set(gca, 'XTickLabel', {'DMI AMARES', 'DMI Int', '13C HP AMARES', '13C HP Int'});
ylabel('Std Dev of Correlations');
title('Ratio Correlation Consistency');
grid on;

% Nuclear comparison for ratios
subplot(2, 3, 6);
dmi_mean_ratio_corr = mean([dmi_ratio_corr_amares, dmi_ratio_corr_integration], 2);
hp_mean_ratio_corr = mean([hp_ratio_corr_amares, hp_ratio_corr_integration], 2);
plot(fixed_noise_levels, dmi_mean_ratio_corr, 'b-o', 'LineWidth', 2, 'DisplayName', '2H DMI (mean)');
hold on;
plot(fixed_noise_levels, hp_mean_ratio_corr, 'r-o', 'LineWidth', 2, 'DisplayName', '13C HP (mean)');
set(gca, 'XScale', 'log');
xlabel('Noise Power'); ylabel('Mean Correlation');
title('Nuclear Ratio Correlation Comparison');
legend('Location', 'best'); grid on;

sgtitle('Ratio Correlation Analysis: Parameter Variation Study', 'FontSize', 16);

% Store ratio correlation results in param_var_results
param_var_results.dmi_ratios_amares_all = param_var_dmi_ratios_amares_all;
param_var_results.dmi_ratios_integration_all = param_var_dmi_ratios_integration_all;
param_var_results.hp_ratios_amares_all = param_var_hp_ratios_amares_all;
param_var_results.hp_ratios_integration_all = param_var_hp_ratios_integration_all;
param_var_results.dmi_ratio_corr_amares = dmi_ratio_corr_amares;
param_var_results.dmi_ratio_corr_integration = dmi_ratio_corr_integration;
param_var_results.hp_ratio_corr_amares = hp_ratio_corr_amares;
param_var_results.hp_ratio_corr_integration = hp_ratio_corr_integration;

fprintf('Ratio correlation analysis complete!\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LINEWIDTH VARIATION ANALYSIS (High SNR)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyze how linewidth affects parameter recovery at high SNR
% This shows whether integration vs AMARES differ systematically with peak width

fprintf('\n========================================\n');
fprintf('Starting Linewidth Variation Analysis...\n');
fprintf('========================================\n');

%% Linewidth analysis parameters
linewidth_levels = linspace(5,15,10); %[5, 8, 10, 15, 20, 30, 40, 50]; % Hz
%NOTE! The previously published 13C PK has bounds of around 0-20 in most
%resonances, 0-120 Hz for pyruvate, and 5-15 Hz for a reference. This is
%kinda a crude way of doing it -- we set a global linewidth. So for now
%it's 5:15 in 10 steps. 
n_linewidths = length(linewidth_levels);
n_mc_linewidth = 20; % Fewer iterations since high SNR
fixed_noise_linewidth = 0.001; % Very low noise (high SNR)

fprintf('Testing %d linewidth levels at noise=%.4f with %d MC iterations each\n', ...
    n_linewidths, fixed_noise_linewidth, n_mc_linewidth);

%% 2H DMI Linewidth Analysis
fprintf('\n--- 2H DMI Linewidth Analysis ---\n');

% Initialize storage
vmax_lw_dmi_amares = zeros(n_linewidths, n_mc_linewidth);
kg_lw_dmi_amares = zeros(n_linewidths, n_mc_linewidth);
vmax_lw_dmi_integration = zeros(n_linewidths, n_mc_linewidth);
kg_lw_dmi_integration = zeros(n_linewidths, n_mc_linewidth);

for lw_idx = 1:n_linewidths
    current_lw = linewidth_levels(lw_idx);
    fprintf('DMI: Processing Linewidth %.1f Hz (%d/%d)\n', current_lw, lw_idx, n_linewidths);

    vmax_lw_dmi_amares_temp = zeros(1, n_mc_linewidth);
    kg_lw_dmi_amares_temp = zeros(1, n_mc_linewidth);
    vmax_lw_dmi_integration_temp = zeros(1, n_mc_linewidth);
    kg_lw_dmi_integration_temp = zeros(1, n_mc_linewidth);

    for mc_idx = 1:n_mc_linewidth
        try
            % Generate synthetic 2H spectra with current linewidth
            fake_spectra_dmi_lw = generate_fake_fids('2H', ...
                'noise_level', fixed_noise_linewidth, ...
                'fixed_linewidth', current_lw, ...
                'flag_adjust_linewidths', true, ...
                'a1', dmi_a1, 'kp', dmi_kp, 'v', dmi_v, ...
                'kg', dmi_kg, 'Vmax', dmi_Vmax, 'kl', dmi_kl);

            % Quantify using AMARES
            quantified_lw_dmi_amares = quantify_spectra_pipline(fake_spectra_dmi_lw, 'amares');
            % Quantify using Integration
            quantified_lw_dmi_integration = quantify_spectra_pipline(fake_spectra_dmi_lw, 'integration');

            % Extract reference peak for normalization
            ref_idx = find(strcmp(fake_spectra_dmi_lw.metabolite_names, 'CDCl3'));
            if ~isempty(ref_idx)
                ref_amares = quantified_lw_dmi_amares.quantified_amplitudes(ref_idx, :);
                ref_integration = quantified_lw_dmi_integration.quantified_amplitudes(ref_idx, :);
                normFacDmi_lw = mean(ref_amares ./ ref_integration);
            else
                normFacDmi_lw = 1;
            end

            % Extract glucose and lactate
            glucose_idx = find(strcmp(fake_spectra_dmi_lw.metabolite_names, 'Glucose'));
            lactate_idx = find(strcmp(fake_spectra_dmi_lw.metabolite_names, 'Lactate'));

            if ~isempty(glucose_idx) && ~isempty(lactate_idx)
                % AMARES results
                glucose_amares = quantified_lw_dmi_amares.quantified_amplitudes(glucose_idx, :);
                lactate_amares = quantified_lw_dmi_amares.quantified_amplitudes(lactate_idx, :);

                % Integration results
                glucose_integration = quantified_lw_dmi_integration.quantified_amplitudes(glucose_idx, :) * normFacDmi_lw;
                lactate_integration = quantified_lw_dmi_integration.quantified_amplitudes(lactate_idx, :) * normFacDmi_lw;

                % Fit DMI kinetics
                dmi_params_fixed = struct('f', 0.6, 'Km', 2.8);
                dmi_params_est = struct('a1', 15.0, 'kp', 0.05, 'v', 0.6, ...
                    'kg', 0.1, 'Vmax', 0.8, 'kl', 0.1);

                [dmi_fit_amares, ~, ~, ~] = fit_dmi_kinetics(...
                    glucose_amares', lactate_amares', fake_spectra_dmi_lw.time_points/60, ...
                    dmi_params_fixed, dmi_params_est);

                [dmi_fit_integration, ~, ~, ~] = fit_dmi_kinetics(...
                    glucose_integration', lactate_integration', fake_spectra_dmi_lw.time_points/60, ...
                    dmi_params_fixed, dmi_params_est);

                % Store results
                vmax_lw_dmi_amares_temp(mc_idx) = dmi_fit_amares.Vmax;
                kg_lw_dmi_amares_temp(mc_idx) = dmi_fit_amares.kg;
                vmax_lw_dmi_integration_temp(mc_idx) = dmi_fit_integration.Vmax;
                kg_lw_dmi_integration_temp(mc_idx) = dmi_fit_integration.kg;
            else
                vmax_lw_dmi_amares_temp(mc_idx) = NaN;
                kg_lw_dmi_amares_temp(mc_idx) = NaN;
                vmax_lw_dmi_integration_temp(mc_idx) = NaN;
                kg_lw_dmi_integration_temp(mc_idx) = NaN;
            end

        catch ME
            % keyboard;  % uncomment to debug
            fprintf('DMI Linewidth Error at LW=%.1f Hz, iter %d: %s\n', current_lw, mc_idx, ME.message);
            vmax_lw_dmi_amares_temp(mc_idx) = NaN;
            kg_lw_dmi_amares_temp(mc_idx) = NaN;
            vmax_lw_dmi_integration_temp(mc_idx) = NaN;
            kg_lw_dmi_integration_temp(mc_idx) = NaN;
        end
    end

    % Store results
    vmax_lw_dmi_amares(lw_idx, :) = vmax_lw_dmi_amares_temp;
    kg_lw_dmi_amares(lw_idx, :) = kg_lw_dmi_amares_temp;
    vmax_lw_dmi_integration(lw_idx, :) = vmax_lw_dmi_integration_temp;
    kg_lw_dmi_integration(lw_idx, :) = kg_lw_dmi_integration_temp;
end

%% Calculate DMI linewidth statistics
vmax_lw_mean_dmi_amares = nanmean(vmax_lw_dmi_amares, 2);
vmax_lw_std_dmi_amares = nanstd(vmax_lw_dmi_amares, [], 2);
vmax_lw_ci_lower_dmi_amares = prctile(vmax_lw_dmi_amares, 2.5, 2);
vmax_lw_ci_upper_dmi_amares = prctile(vmax_lw_dmi_amares, 97.5, 2);

vmax_lw_mean_dmi_integration = nanmean(vmax_lw_dmi_integration, 2);
vmax_lw_std_dmi_integration = nanstd(vmax_lw_dmi_integration, [], 2);
vmax_lw_ci_lower_dmi_integration = prctile(vmax_lw_dmi_integration, 2.5, 2);
vmax_lw_ci_upper_dmi_integration = prctile(vmax_lw_dmi_integration, 97.5, 2);

kg_lw_mean_dmi_amares = nanmean(kg_lw_dmi_amares, 2);
kg_lw_mean_dmi_integration = nanmean(kg_lw_dmi_integration, 2);

%% 13C HP Linewidth Analysis
fprintf('\n--- 13C HP Linewidth Analysis ---\n');

% Initialize storage
kpb_lw_hp_amares = zeros(n_linewidths, n_mc_linewidth);
r1b_lw_hp_amares = zeros(n_linewidths, n_mc_linewidth);
kpb_lw_hp_integration = zeros(n_linewidths, n_mc_linewidth);
r1b_lw_hp_integration = zeros(n_linewidths, n_mc_linewidth);

for lw_idx = 1:n_linewidths
    current_lw = linewidth_levels(lw_idx);
    fprintf('13C: Processing Linewidth %.1f Hz (%d/%d)\n', current_lw, lw_idx, n_linewidths);

    kpb_lw_hp_amares_temp = zeros(1, n_mc_linewidth);
    r1b_lw_hp_amares_temp = zeros(1, n_mc_linewidth);
    kpb_lw_hp_integration_temp = zeros(1, n_mc_linewidth);
    r1b_lw_hp_integration_temp = zeros(1, n_mc_linewidth);

    for mc_idx = 1:n_mc_linewidth
        try
            % Generate synthetic 13C spectra with current linewidth
            fake_spectra_hp_lw = generate_fake_fids('13C', ...
                'noise_level', fixed_noise_linewidth, ...
                'fixed_linewidth', current_lw, ...
                'flag_adjust_linewidths', true, ...
                'k_true', k_true_hp, 'r1_pyr', r1_pyr_hp, ...
                'r1_bic', r1_bic_hp, 'pyr_0', pyr_0_hp, 'FlipAngle', FlipAngle_hp);

            % Quantify using AMARES
            quantified_lw_hp_amares = quantify_spectra_pipline(fake_spectra_hp_lw, 'amares');
            % Quantify using Integration
            quantified_lw_hp_integration = quantify_spectra_pipline(fake_spectra_hp_lw, 'integration');

            % Extract reference peak for normalization
            ref_idx = find(strcmp(fake_spectra_hp_lw.metabolite_names, '1-13C UreaRef'));
            if ~isempty(ref_idx)
                ref_amares = quantified_lw_hp_amares.quantified_amplitudes(ref_idx, :);
                ref_integration = quantified_lw_hp_integration.quantified_amplitudes(ref_idx, :);
                normFacHP_lw = mean(ref_amares ./ ref_integration);
            else
                normFacHP_lw = 1;
            end

            % Extract pyruvate and bicarbonate
            pyr_idx = find(strcmp(fake_spectra_hp_lw.metabolite_names, '1-13C Pyruvate'));
            bic_idx = find(strcmp(fake_spectra_hp_lw.metabolite_names, '1-13C Bicarbonate'));

            if ~isempty(pyr_idx) && ~isempty(bic_idx)
                % AMARES results
                pyr_amares = quantified_lw_hp_amares.quantified_amplitudes(pyr_idx, :);
                bic_amares = quantified_lw_hp_amares.quantified_amplitudes(bic_idx, :);

                % Integration results
                pyr_integration = quantified_lw_hp_integration.quantified_amplitudes(pyr_idx, :) * normFacHP_lw;
                bic_integration = quantified_lw_hp_integration.quantified_amplitudes(bic_idx, :) * normFacHP_lw;

                % Prepare signals for fitting
                n_timepoints = length(fake_spectra_hp_lw.time_points);
                signals_amares = zeros(1, 2, n_timepoints);
                signals_amares(1, 1, :) = pyr_amares;
                signals_amares(1, 2, :) = bic_amares;

                signals_integration = zeros(1, 2, n_timepoints);
                signals_integration(1, 1, :) = pyr_integration;
                signals_integration(1, 2, :) = bic_integration;

                flip_angles_pyr = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                flip_angles_bic = deg2rad(FlipAngle_hp) * ones(1, n_timepoints);
                flip_angles = [flip_angles_pyr; flip_angles_bic];

                TR_hp = fake_spectra_hp_lw.time_points(2)-fake_spectra_hp_lw.time_points(1);

                % Fit HP kinetics to AMARES-quantified data
                params_fixed = struct('R1P', r1_pyr_hp);
                params_est = struct('kPL', 0.05, 'R1L', 0.05, 'B0_start', 0.0,'kPL_lb',0,'R1L_lb',0,'R1P_lb',0);

                [params_fit_amares, ~, ~, ~] = fit_pyr_kinetics(signals_amares, TR_hp, ...
                    flip_angles, params_fixed, params_est);

                % Fit HP kinetics to Integration-quantified data
                [params_fit_integration, ~, ~, ~] = fit_pyr_kinetics(signals_integration, TR_hp, ...
                    flip_angles, params_fixed, params_est);

                % Store results
                kpb_lw_hp_amares_temp(mc_idx) = params_fit_amares.kPL;
                r1b_lw_hp_amares_temp(mc_idx) = params_fit_amares.R1L;
                kpb_lw_hp_integration_temp(mc_idx) = params_fit_integration.kPL;
                r1b_lw_hp_integration_temp(mc_idx) = params_fit_integration.R1L;
            else
                kpb_lw_hp_amares_temp(mc_idx) = NaN;
                r1b_lw_hp_amares_temp(mc_idx) = NaN;
                kpb_lw_hp_integration_temp(mc_idx) = NaN;
                r1b_lw_hp_integration_temp(mc_idx) = NaN;
            end

        catch ME
            % keyboard;  % uncomment to debug
            fprintf('13C Linewidth Error at LW=%.1f Hz, iter %d: %s\n', current_lw, mc_idx, ME.message);
            kpb_lw_hp_amares_temp(mc_idx) = NaN;
            r1b_lw_hp_amares_temp(mc_idx) = NaN;
            kpb_lw_hp_integration_temp(mc_idx) = NaN;
            r1b_lw_hp_integration_temp(mc_idx) = NaN;
        end
    end

    % Store results
    kpb_lw_hp_amares(lw_idx, :) = kpb_lw_hp_amares_temp;
    r1b_lw_hp_amares(lw_idx, :) = r1b_lw_hp_amares_temp;
    kpb_lw_hp_integration(lw_idx, :) = kpb_lw_hp_integration_temp;
    r1b_lw_hp_integration(lw_idx, :) = r1b_lw_hp_integration_temp;
end

%% Calculate 13C linewidth statistics
kpb_lw_mean_hp_amares = nanmean(kpb_lw_hp_amares, 2);
kpb_lw_std_hp_amares = nanstd(kpb_lw_hp_amares, [], 2);
kpb_lw_ci_lower_hp_amares = prctile(kpb_lw_hp_amares, 2.5, 2);
kpb_lw_ci_upper_hp_amares = prctile(kpb_lw_hp_amares, 97.5, 2);

kpb_lw_mean_hp_integration = nanmean(kpb_lw_hp_integration, 2);
kpb_lw_std_hp_integration = nanstd(kpb_lw_hp_integration, [], 2);
kpb_lw_ci_lower_hp_integration = prctile(kpb_lw_hp_integration, 2.5, 2);
kpb_lw_ci_upper_hp_integration = prctile(kpb_lw_hp_integration, 97.5, 2);

r1b_lw_mean_hp_amares = nanmean(r1b_lw_hp_amares, 2);
r1b_lw_mean_hp_integration = nanmean(r1b_lw_hp_integration, 2);

%% Plot Linewidth Analysis Results

% 2H DMI Linewidth Analysis Figure
figure('Position', [100 100 1400 900]); set(gcf, 'color', 'w');
sgtitle(sprintf('2H DMI: Parameter Recovery vs Linewidth (Noise Level = %.4f)', fixed_noise_linewidth), ...
    'FontSize', 16, 'FontWeight', 'bold');

% Vmax vs Linewidth
subplot(2, 3, 1);
plot(linewidth_levels, vmax_lw_mean_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'AMARES Estimated Vmax');
hold on;
plot(linewidth_levels, vmax_lw_mean_dmi_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Integration Estimated Vmax');
yline(dmi_params_true.Vmax, 'k--', 'LineWidth', 2, 'DisplayName', 'True Vmax');

% Add confidence intervals
x_region = [linewidth_levels, fliplr(linewidth_levels)];
y_region_amares = [vmax_lw_ci_upper_dmi_amares', fliplr(vmax_lw_ci_lower_dmi_amares')];
fill(x_region, y_region_amares, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
y_region_int = [vmax_lw_ci_upper_dmi_integration', fliplr(vmax_lw_ci_lower_dmi_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xlabel('Linewidth (Hz)'); ylabel('Estimated Vmax (mM/min)');
title('Vmax Estimation vs Linewidth');
legend('Location', 'best'); grid on;

% Vmax Bias vs Linewidth
subplot(2, 3, 2);
bias_vmax_amares = 100 * (vmax_lw_mean_dmi_amares - dmi_params_true.Vmax) / dmi_params_true.Vmax;
bias_vmax_integration = 100 * (vmax_lw_mean_dmi_integration - dmi_params_true.Vmax) / dmi_params_true.Vmax;
plot(linewidth_levels, bias_vmax_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES');
hold on;
plot(linewidth_levels, bias_vmax_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration');
yline(0, 'k--', 'LineWidth', 1);
xlabel('Linewidth (Hz)'); ylabel('Vmax Bias (%)');
title('Vmax Bias vs Linewidth');
legend('Location', 'best'); grid on;

% Vmax Variance vs Linewidth
subplot(2, 3, 3);
plot(linewidth_levels, vmax_lw_std_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES Std');
hold on;
plot(linewidth_levels, vmax_lw_std_dmi_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration Std');
xlabel('Linewidth (Hz)'); ylabel('Vmax Std Dev');
title('Vmax Variance vs Linewidth');
legend('Location', 'best'); grid on;

% kg vs Linewidth
subplot(2, 3, 4);
plot(linewidth_levels, kg_lw_mean_dmi_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES K_g');
hold on;
plot(linewidth_levels, kg_lw_mean_dmi_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration K_g');
yline(dmi_params_true.kg, 'k--', 'LineWidth', 2, 'DisplayName', 'True K_g');
xlabel('Linewidth (Hz)'); ylabel('Estimated K_g (min^{-1})');
title('K_g Estimation vs Linewidth');
legend('Location', 'best'); grid on;

% Method comparison scatter
subplot(2, 3, 5);
scatter(vmax_lw_mean_dmi_amares, vmax_lw_mean_dmi_integration, 100, linewidth_levels, 'filled');
hold on;
plot([min([vmax_lw_mean_dmi_amares; vmax_lw_mean_dmi_integration]), ...
    max([vmax_lw_mean_dmi_amares; vmax_lw_mean_dmi_integration])], ...
    [min([vmax_lw_mean_dmi_amares; vmax_lw_mean_dmi_integration]), ...
    max([vmax_lw_mean_dmi_amares; vmax_lw_mean_dmi_integration])], 'k--');
colorbar; colormap(jet);
xlabel('AMARES Vmax'); ylabel('Integration Vmax');
title('Method Comparison (color = linewidth)');
axis equal; grid on;

% Summary statistics
subplot(2, 3, 6);
rmse_amares = sqrt(mean((vmax_lw_mean_dmi_amares - dmi_params_true.Vmax).^2));
rmse_integration = sqrt(mean((vmax_lw_mean_dmi_integration - dmi_params_true.Vmax).^2));
mean_bias_amares = mean(bias_vmax_amares);
mean_bias_integration = mean(bias_vmax_integration);

bar([rmse_amares, rmse_integration; mean_bias_amares, mean_bias_integration]);
set(gca, 'XTickLabel', {'RMSE', 'Mean Bias (%)'});
ylabel('Value');
title('Overall Performance Summary');
legend('AMARES', 'Integration', 'Location', 'best');
grid on;

% 13C HP Linewidth Analysis Figure
figure('Position', [150 50 1400 900]); set(gcf, 'color', 'w');
sgtitle(sprintf('13C HP: Parameter Recovery vs Linewidth (Noise Level = %.4f)', fixed_noise_linewidth), ...
    'FontSize', 16, 'FontWeight', 'bold');

% kPB vs Linewidth
subplot(2, 3, 1);
plot(linewidth_levels, kpb_lw_mean_hp_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'AMARES Estimated kPB');
hold on;
plot(linewidth_levels, kpb_lw_mean_hp_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Integration Estimated kPB');
yline(k_true_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True kPB');

% Add confidence intervals
x_region = [linewidth_levels, fliplr(linewidth_levels)];
y_region_amares = [kpb_lw_ci_upper_hp_amares', fliplr(kpb_lw_ci_lower_hp_amares')];
fill(x_region, y_region_amares, 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
y_region_int = [kpb_lw_ci_upper_hp_integration', fliplr(kpb_lw_ci_lower_hp_integration')];
fill(x_region, y_region_int, 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xlabel('Linewidth (Hz)'); ylabel('Estimated kPB (1/s)');
title('kPB Estimation vs Linewidth');
legend('Location', 'best'); grid on;

% kPB Bias vs Linewidth
subplot(2, 3, 2);
bias_kpb_amares = 100 * (kpb_lw_mean_hp_amares - k_true_hp) / k_true_hp;
bias_kpb_integration = 100 * (kpb_lw_mean_hp_integration - k_true_hp) / k_true_hp;
plot(linewidth_levels, bias_kpb_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES');
hold on;
plot(linewidth_levels, bias_kpb_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration');
yline(0, 'k--', 'LineWidth', 1);
xlabel('Linewidth (Hz)'); ylabel('kPB Bias (%)');
title('kPB Bias vs Linewidth');
legend('Location', 'best'); grid on;

% kPB Variance vs Linewidth
subplot(2, 3, 3);
plot(linewidth_levels, kpb_lw_std_hp_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES Std');
hold on;
plot(linewidth_levels, kpb_lw_std_hp_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration Std');
xlabel('Linewidth (Hz)'); ylabel('kPB Std Dev');
title('kPB Variance vs Linewidth');
legend('Location', 'best'); grid on;

% R1B vs Linewidth
subplot(2, 3, 4);
plot(linewidth_levels, r1b_lw_mean_hp_amares, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'AMARES R1B');
hold on;
plot(linewidth_levels, r1b_lw_mean_hp_integration, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Integration R1B');
yline(r1_bic_hp, 'k--', 'LineWidth', 2, 'DisplayName', 'True R1B');
xlabel('Linewidth (Hz)'); ylabel('Estimated R1B (1/s)');
title('R1B Estimation vs Linewidth');
legend('Location', 'best'); grid on;

% Method comparison scatter
subplot(2, 3, 5);
scatter(kpb_lw_mean_hp_amares, kpb_lw_mean_hp_integration, 100, linewidth_levels, 'filled');
hold on;
plot([min([kpb_lw_mean_hp_amares; kpb_lw_mean_hp_integration]), ...
    max([kpb_lw_mean_hp_amares; kpb_lw_mean_hp_integration])], ...
    [min([kpb_lw_mean_hp_amares; kpb_lw_mean_hp_integration]), ...
    max([kpb_lw_mean_hp_amares; kpb_lw_mean_hp_integration])], 'k--');
colorbar; colormap(jet);
xlabel('AMARES kPB'); ylabel('Integration kPB');
title('Method Comparison (color = linewidth)');
axis equal; grid on;

% Summary statistics
subplot(2, 3, 6);
rmse_amares = sqrt(mean((kpb_lw_mean_hp_amares - k_true_hp).^2));
rmse_integration = sqrt(mean((kpb_lw_mean_hp_integration - k_true_hp).^2));
mean_bias_amares = mean(bias_kpb_amares);
mean_bias_integration = mean(bias_kpb_integration);

bar([rmse_amares, rmse_integration; mean_bias_amares, mean_bias_integration]);
set(gca, 'XTickLabel', {'RMSE', 'Mean Bias (%)'});
ylabel('Value');
title('Overall Performance Summary');
legend('AMARES', 'Integration', 'Location', 'best');
grid on;

% Store linewidth results
linewidth_results = struct();
linewidth_results.linewidth_levels = linewidth_levels;
linewidth_results.fixed_noise = fixed_noise_linewidth;
linewidth_results.n_mc = n_mc_linewidth;

% DMI results
linewidth_results.dmi.vmax_amares = vmax_lw_dmi_amares;
linewidth_results.dmi.kg_amares = kg_lw_dmi_amares;
linewidth_results.dmi.vmax_integration = vmax_lw_dmi_integration;
linewidth_results.dmi.kg_integration = kg_lw_dmi_integration;
linewidth_results.dmi.vmax_mean_amares = vmax_lw_mean_dmi_amares;
linewidth_results.dmi.vmax_mean_integration = vmax_lw_mean_dmi_integration;
linewidth_results.dmi.kg_mean_amares = kg_lw_mean_dmi_amares;
linewidth_results.dmi.kg_mean_integration = kg_lw_mean_dmi_integration;

% HP results
linewidth_results.hp.kpb_amares = kpb_lw_hp_amares;
linewidth_results.hp.r1b_amares = r1b_lw_hp_amares;
linewidth_results.hp.kpb_integration = kpb_lw_hp_integration;
linewidth_results.hp.r1b_integration = r1b_lw_hp_integration;
linewidth_results.hp.kpb_mean_amares = kpb_lw_mean_hp_amares;
linewidth_results.hp.kpb_mean_integration = kpb_lw_mean_hp_integration;
linewidth_results.hp.r1b_mean_amares = r1b_lw_mean_hp_amares;
linewidth_results.hp.r1b_mean_integration = r1b_lw_mean_hp_integration;

fprintf('\nLinewidth variation analysis complete!\n');
fprintf('Key findings at high SNR (noise=%.4f):\n', fixed_noise_linewidth);
fprintf('  2H DMI Vmax RMSE: AMARES=%.4f, Integration=%.4f\n', ...
    sqrt(mean((vmax_lw_mean_dmi_amares - dmi_params_true.Vmax).^2)), ...
    sqrt(mean((vmax_lw_mean_dmi_integration - dmi_params_true.Vmax).^2)));
fprintf('  13C HP kPB RMSE: AMARES=%.4f, Integration=%.4f\n', ...
    sqrt(mean((kpb_lw_mean_hp_amares - k_true_hp).^2)), ...
    sqrt(mean((kpb_lw_mean_hp_integration - k_true_hp).^2)));

%% Save all results
save('comprehensive_noise_analysis_results.mat', 'linewidth_results', '-append');

fprintf('\nAnalysis complete! Results saved to comprehensive_noise_analysis_results.mat\n');
