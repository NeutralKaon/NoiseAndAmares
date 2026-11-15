function quantified_results = quantify_spectra_pipline(fake_spectra, method, varargin)
%QUANTIFY_SPECTRA Quantify synthetic spectra using AMARES or integration
%
% USAGE:
%   results = quantify_spectra(fake_spectra, 'amares');
%   results = quantify_spectra(fake_spectra, 'integration');
%   results = quantify_spectra(fake_spectra, 'integration', 'integration_width_factor', 2.0);
%
% INPUTS:
%   fake_spectra - Output from generate_fake_spectra
%   method       - 'amares' or 'integration'
%
% OPTIONAL PARAMETERS:
%   plot_fits                - Show fitting plots (default: false)
%   integration_width_factor - Widen integration bounds by this factor (default: 1.5)
%                              Only applies to 'integration' method
%   baseline_poly_order      - Polynomial order for baseline correction (default: 1)
%                              0=no baseline, 1=linear, 2=quadratic, etc.
%                              Only applies to 'integration' method

%% Parse inputs
p = inputParser;
addRequired(p, 'fake_spectra', @isstruct);
addRequired(p, 'method', @(x) ismember(x, {'amares', 'integration'}));
addParameter(p, 'plot_fits', false, @islogical);
addParameter(p, 'integration_width_factor', 1.5, @(x) x > 0); % Widen integration by this factor (bounds mode)
addParameter(p, 'baseline_poly_order', 1, @(x) x >= 0 && x <= 3); % Polynomial order for baseline
addParameter(p, 'window_mode', 'bounds', @(x) ismember(x, {'bounds','linewidth'})); % integration window definition
addParameter(p, 'linewidth_factor', 5, @(x) x > 0); % half-width in units of HWHM (linewidth mode)

parse(p, fake_spectra, method, varargin{:});
params = p.Results;

%% Extract data
n_timepoints = size(fake_spectra.noisy_fids, 2);
n_metabolites = length(fake_spectra.metabolite_names);
quantified_amplitudes = zeros(n_metabolites, n_timepoints);

fprintf('Quantifying %d timepoints using %s method...\n', n_timepoints, method);

%% Quantify each timepoint
for t = 1:n_timepoints
    if mod(t, 10) == 0
        fprintf('  Processing timepoint %d/%d\n', t, n_timepoints);
    end
    
    current_fid = fake_spectra.noisy_fids(:, t);
    
    if strcmp(method, 'amares')
        % AMARES quantification using OXSA
        if isfield(fake_spectra,'pk') && ~isempty(fake_spectra.pk)
            % caller supplied its own prior knowledge (noise_statistics_test.m)
        elseif strcmp(fake_spectra.nucleus,'13C')
            % 5 peaks has no reference 
            %fake_spectra.pk=AMARES.priorKnowledge.PK_3T_1C_13C_5Peaks(0,1); 
            % 6 peaks has a "concentration reference" that is constnat
            % throughout the generated simulations 
            fake_spectra.pk=AMARES.priorKnowledge.PK_3T_1C_13C_6Peaks(0,1); 
        elseif strcmp(fake_spectra.nucleus,'2H')
            %fake_spectra.pk=AMARES.priorKnowledge.PK_3T_2H_5Peaks(0); 
            fake_spectra.pk=AMARES.priorKnowledge.PK_3T_2H_6Peaks(0); 
        else
            error('Check input!')
        end
        amps = quantify_with_oxsa_amares(current_fid, fake_spectra, params.plot_fits);
  
    else
        % Frequency domain integration
        amps = quantify_with_integration(current_fid, fake_spectra, params.plot_fits, ...
            params.integration_width_factor, params.baseline_poly_order, ...
            params.window_mode, params.linewidth_factor);
    end
    
    quantified_amplitudes(:, t) = amps;
end

%% Calculate the time-averaged spectrum, quantify it, package it separately for later ratios 

if strcmp(method, 'amares')
    summed_amps = quantify_with_oxsa_amares(mean(fake_spectra.noisy_fids,2), fake_spectra, params.plot_fits);
    
else
    % Frequency domain integration
    summed_amps = quantify_with_integration(mean(fake_spectra.noisy_fids,2), fake_spectra, ...
        params.plot_fits, params.integration_width_factor, params.baseline_poly_order, ...
        params.window_mode, params.linewidth_factor);
end


%% Calculate errors and package results
quantified_results = struct();
quantified_results.method = method;
quantified_results.quantified_amplitudes = quantified_amplitudes;
quantified_results.true_amplitudes = fake_spectra.true_amplitudes;
quantified_results.time_points = fake_spectra.time_points;

if strcmp(method, 'amares')
    %Order isn't necessarily guaranteed...
    quantified_results.metabolite_names = {fake_spectra.pk.bounds.peakName};
end 
% Calculate error metrics
relative_errors = abs(quantified_amplitudes - fake_spectra.true_amplitudes) ./ ...
                 (fake_spectra.true_amplitudes + eps);
quantified_results.absolute_errors = abs(quantified_amplitudes - fake_spectra.true_amplitudes); 
quantified_results.scaled_to_unity_absolute_errors = abs(quantified_amplitudes./quantified_amplitudes(1,1) - fake_spectra.true_amplitudes./fake_spectra.true_amplitudes(1,1));

quantified_results.relative_errors = relative_errors;
quantified_results.mean_relative_error = mean(relative_errors(:));
quantified_results.mean_absolute_error = mean(quantified_results.absolute_errors(:)); 
quantified_results.scaled_to_unity_mean_absolute_error = mean(quantified_results.scaled_to_unity_absolute_errors(:)); 
quantified_results.summed_amps = summed_amps; 

if strcmp(fake_spectra.nucleus,'13C')
    %Lac/pyr and bic/pyr – both the same here...
        quantified_results.ratios_of_interest = ...
        [quantified_results.summed_amps(3)/(quantified_results.summed_amps(1)); ...
        quantified_results.summed_amps(2)/(quantified_results.summed_amps(1))];
elseif strcmp(fake_spectra.nucleus,'2H')
    %Glc/total 2H and lac/total 2H 
    quantified_results.ratios_of_interest = ...
        [quantified_results.summed_amps(2)/sum(quantified_results.summed_amps); ...
        quantified_results.summed_amps(4)/sum(quantified_results.summed_amps)];
else
    error('Check input!')
end
% Calculate correlations for each metabolite
correlations = zeros(n_metabolites, 1);
for i = 1:n_metabolites
    true_vals = fake_spectra.true_amplitudes(i, :);
    quant_vals = quantified_amplitudes(i, :);
    if std(true_vals) > eps && std(quant_vals) > eps
        corr_matrix = corrcoef(true_vals, quant_vals);
        correlations(i) = corr_matrix(1, 2);
    else
        correlations(i) = NaN;
    end
end
quantified_results.correlations = correlations;

fprintf('Quantification completed.\n'); 
fprintf('Mean relative error:  %.2f%%\n',quantified_results.mean_relative_error * 100);
fprintf('Mean absolute error:  %.2f%%\n',quantified_results.mean_absolute_error * 100);
fprintf('Mean unity-scaled AE: %.2f%%\n',quantified_results.scaled_to_unity_mean_absolute_error * 100);

end

function amplitudes = quantify_with_oxsa_amares(fid, fake_spectra, show_plot)
%QUANTIFY_WITH_OXSA_AMARES Quantify FID using OXSA AMARES
%
% Fid: fids; fake_spectra: struct containing useful stuff; show_plot:
% logical. 

n_metabolites = length(fake_spectra.metabolite_names);
amplitudes = zeros(n_metabolites, 1);

try
    % Set up experiment parameters for OXSA
    exptParams = struct();
    exptParams.dwellTime = fake_spectra.dwell_time;
    exptParams.beginTime = 0; % ms
    exptParams.centreFreq = fake_spectra.base_freq;
    exptParams.acqVoxSize = [1, 1, 1]; % Dummy voxel size
    
    % Extract prior knowledge
    pk = fake_spectra.pk; 
    
    %Fake more things
    exptParams.samples=size(fid,1);
    exptParams.imagingFrequency=fake_spectra.base_freq*1e-6; 
    exptParams.timeAxis=fake_spectra.time_axis;
    
    exptParams.freqAxis = ((-exptParams.samples/2):(exptParams.samples/2-1)).'/(fake_spectra.dwell_time*exptParams.samples);
    exptParams.ppmAxis=exptParams.freqAxis /(exptParams.imagingFrequency); 
    %10.7084e6

    
    % Call OXSA AMARES fitting
    [fitResults, fitStatus, ~, ~] = AMARES.amaresFit(fid, exptParams, pk, show_plot);
    
    if fitStatus.EXITFLAG >= 0
        % Extract amplitudes from fit results
        amplitudes = fitResults.amplitude; 
    else
        warning('AMARES fitting failed for timepoint');
        amplitudes(:) = 0;
    end
    
catch ME
    warning('OXSA AMARES error: %s', ME.message);
    amplitudes(:) = 0;
end

end

function amplitudes = quantify_with_integration(fid, fake_spectra, show_plot, integration_width_factor, baseline_poly_order, window_mode, linewidth_factor)
%QUANTIFY_WITH_INTEGRATION Quantify using frequency domain integration
% with wider integration bounds and polynomial baseline correction
%
% INPUTS:
%   fid                      - FID to quantify
%   fake_spectra             - Spectrum structure with prior knowledge
%   show_plot                - Boolean to show plots
%   integration_width_factor - Factor to widen integration bounds (default: 1.5)
%   baseline_poly_order      - Polynomial order for baseline (default: 1, 0=none)

if nargin < 4, integration_width_factor = 1.5; end
if nargin < 5, baseline_poly_order = 1; end
if nargin < 6 || isempty(window_mode), window_mode = 'bounds'; end
if nargin < 7 || isempty(linewidth_factor), linewidth_factor = 5; end

n_metabolites = length(fake_spectra.metabolite_names);
amplitudes = zeros(n_metabolites, 1);

% Convert to frequency domain
spectrum = specFft(fid); %This is an OXSA dependency and puts a 1/sqrt(2) in the front point of the FID.
%Or %fftshift(fft(fid));


%Fake more things
exptParams.samples=size(fid,1);
exptParams.imagingFrequency=fake_spectra.base_freq*1e-6;
exptParams.timeAxis=fake_spectra.time_axis;

exptParams.freqAxis = ((-exptParams.samples/2):(exptParams.samples/2-1)).'/(fake_spectra.dwell_time*exptParams.samples);
exptParams.ppmAxis=exptParams.freqAxis /(exptParams.imagingFrequency);

freq_axis_ppm =exptParams.ppmAxis;
spectrum_real = real(spectrum);

if(show_plot)
    figure(9152); close(9152); figure(9152); set(gcf,'color','w'); hold on;
end

% Integrate each metabolite peak
for i = 1:n_metabolites
    metabolite_name = fake_spectra.metabolite_names{i};

    % Find metabolite in prior knowledge bounds
    pk_idx = find(strcmp(fake_spectra.prior_knowledge.values.boundsCellArray(:,1), metabolite_name), 1);

    if ~isempty(pk_idx)
        if strcmpi(window_mode, 'linewidth')
            % Scale the integration window to each peak's linewidth, so that the
            % captured Lorentzian-area fraction is matched across peaks of
            % different width (and therefore cancels against the constant
            % reference peak). This removes the per-peak proportionality offset
            % between integration and AMARES seen at high SNR when one peak
            % (e.g. pyruvate, 50 Hz) is much broader than the others (10 Hz).
            % Centre and linewidth are taken from the generated prior-knowledge
            % initial values (IVCellArray), which carry the true generated width.
            iv = fake_spectra.prior_knowledge.values.IVCellArray;
            iv_idx = find(strcmp(iv(:,1), metabolite_name), 1);
            ppm_center = iv{iv_idx, 2};                       % ppm
            lw_hz      = iv{iv_idx, 3};                       % Hz (FWHM)
            hwhm_ppm   = (lw_hz/2) / exptParams.imagingFrequency; % imagingFrequency in MHz
            w_half = linewidth_factor*hwhm_ppm;
            lo = ppm_center - w_half;
            hi = ppm_center + w_half;
            % Non-overlapping windows: where two linewidth-scaled windows would
            % overlap, split the overlap evenly so adjacent windows just meet and
            % never double-count shared signal. Well-separated peaks (e.g. 13C) keep
            % their full linewidth window; only the crowded peaks (e.g. 2H DMI) are
            % trimmed, and only by the amount they actually overlapped.
            all_c  = cell2mat(iv(:,2));
            all_hw = linewidth_factor*(cell2mat(iv(:,3))/2)/exptParams.imagingFrequency;
            bl = find(all_c < ppm_center - 1e-6);   % nearest resonance below
            if ~isempty(bl)
                [cb,ib]=max(all_c(bl)); br = cb + all_hw(bl(ib));
                if br > lo, lo = (br + lo)/2; end
            end
            ab = find(all_c > ppm_center + 1e-6);   % nearest resonance above
            if ~isempty(ab)
                [ca,ia]=min(all_c(ab)); al = ca - all_hw(ab(ia));
                if al < hi, hi = (al + hi)/2; end
            end
            ppm_range  = [lo, hi];
        else
            % Original behaviour: window = AMARES bounds widened by a factor.
            ppm_range_amares = fake_spectra.prior_knowledge.values.boundsCellArray{pk_idx, 2};
            ppm_center = mean(ppm_range_amares);
            ppm_width_amares = diff(ppm_range_amares);
            ppm_width_integration = ppm_width_amares * integration_width_factor;
            ppm_range = [ppm_center - ppm_width_integration/2, ...
                         ppm_center + ppm_width_integration/2];
        end

        % Find indices for integration
        freq_indices = freq_axis_ppm >= min(ppm_range) & freq_axis_ppm <= max(ppm_range);

        if any(freq_indices)
            % Extract region
            ppm_region = freq_axis_ppm(freq_indices);
            spectrum_region = spectrum_real(freq_indices);

            % Polynomial baseline correction
            % Fit polynomial to edges of the region for baseline estimation
            n_edge_points = min(10, floor(length(ppm_region)/4)); % Use 10 points or 1/4 of region

            if length(ppm_region) > 2*n_edge_points && baseline_poly_order >= 1
                % Combine left and right edge points
                edge_ppm = [ppm_region(1:n_edge_points); ppm_region(end-n_edge_points+1:end)];
                edge_spectrum = [spectrum_region(1:n_edge_points); spectrum_region(end-n_edge_points+1:end)];

                % Fit polynomial baseline
                baseline_coeffs = polyfit(edge_ppm, edge_spectrum, baseline_poly_order);
                baseline = polyval(baseline_coeffs, ppm_region);

                % Subtract baseline
                spectrum_corrected = spectrum_region - baseline;
            else
                % If region too small, no baseline correction
                spectrum_corrected = spectrum_region;
            end

            % Integrate baseline-corrected spectrum. Clip to be non-negative:
            % metabolite amounts (and AMARES amplitudes) are positive
            % semidefinite, so a negative integral -- which can arise at low
            % SNR or under peak overlap -- is unphysical and would otherwise
            % drive metabolite ratios negative.
            integration_value = trapz(ppm_region, spectrum_corrected);
            % (noise_statistics_test.m turns this off, as in Fig. 8.)
            if isfield(fake_spectra,'clip_negative') && ~fake_spectra.clip_negative
                amplitudes(i) = integration_value;
            else
                amplitudes(i) = max(integration_value, 0);
            end
        else
            amplitudes(i) = 0;
        end

        if(show_plot)
           plot(freq_axis_ppm, freq_indices, 'DisplayName', metabolite_name);
           hold on;
        end
    else
        amplitudes(i) = 0;
    end
end

if(show_plot)
    plot(freq_axis_ppm, spectrum_real./max(spectrum_real(:)),'DisplayName',['Norm Spec']);
    legend('-DynamicLegend');
            end
        end
