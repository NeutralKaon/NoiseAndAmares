function fake_spectra = generate_fake_fids(nucleus, varargin)
%GENERATE_FAKE_SPECTRA Generate synthetic MRS spectra with kinetic timecourses
%
% USAGE:
%   fake_spectra = generate_fake_spectra('13C');
%   fake_spectra = generate_fake_spectra('2H', 'noise_level', 0.1, 'magnetic_field', 7.0);
%
% INPUTS:
%   nucleus      - '13C' or '2H'
%
% OPTIONAL PARAMETERS:
%   noise_level  - Gaussian noise level (default: 0.1)
%   magnetic_field - Field strength in Tesla (default: 3.0)
%   n_timepoints - Number of dynamic timepoints (default: 120 for 13C, 40 for 2H)
%   n_points     - FID length in points (default: 2048)
%   ppm_offset   - Chemical shift offset in ppm (default: 0)
%   default_TR   - TR for each spectrum, defaults to 1s (carbon) and 2 mins (2H)
%   ConcStandardAmp - a numerical "concentration" for a reference peak that
%   doesn't significantly spectrally overlap -- default 1 (DMI) or 1e-2 (13C)
%   fixed_linewidth - ADDITIONAL linewidth in Hz added to all peaks (default: 0.1 Hz = 1e-1)
%                     Applies to both 13C and 2H when flag_adjust_linewidths=true
%                     This is ADDED to the PK structure default linewidths (not set)
%   flag_adjust_linewidths - Boolean to add fixed_linewidth to baseline (default: true)
%                            If false, uses PK structure default linewidths only
%
% Model specific parameters and defaults:
%
% 13C: (all per second or seconds, save FA which is degrees)
% -- k_true = 0.05; r1_pyr = 1/30; r1_bic = 1/20; pyr_0 = 1.0;  FlipAngle=15;
% The model generates "lactate" which is identical to bicarb, "alanine",
% which is
%
%
% 2H: (units more complex, see paper)
% --  a1 =  17.0; kp =  0.0416; v = 0.72; kg=0.141; Vmax= 0.985; kl= 0.124

%% Parse inputs
p = inputParser;
addRequired(p, 'nucleus', @(x) ismember(x, {'13C', '2H'}));
addParameter(p, 'noise_level', 0.1, @(x) x >= 0);
% noise_reference selects what noise_level is a fraction of:
%   'peakamp' - the maximum substrate PEAK AMPLITUDE (pyruvate for 13C, glucose for
%               2H). This matches the amplitude-domain simulation (Fig 7), so the two
%               noise axes are directly comparable.
%   'rmsfid'  - the RMS of the FID at the max-signal timepoint (legacy behaviour).
addParameter(p, 'noise_reference', 'peakamp', @(x) ischar(x) || isstring(x));
addParameter(p, 'magnetic_field', 3.0, @(x) x > 0);
addParameter(p, 'n_timepoints', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'n_points', 2048, @(x) isnumeric(x) && x > 0);
addParameter(p, 'ppm_offset', 0, @isnumeric);
addParameter(p, 'default_tr', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'k_true', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'r1_pyr', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'r1_bic', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'pyr_0', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'FlipAngle', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'ConcStandardAmp', [],@(x) isempty(x) || (isnumeric(x) && x > 0) ); 
addParameter(p, 'fixed_linewidth',1e-1,@(x) isnumeric(x));
addParameter(p, 'flag_adjust_linewidths',true,@(x) islogical(x));


% DMI parameters
addParameter(p, 'a1', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'kp', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'v', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'kg', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'Vmax', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'kl', [], @(x) isempty(x) || (isnumeric(x) && x > 0));


parse(p, nucleus, varargin{:});
params = p.Results;

%% Set nucleus-specific defaults
if strcmp(params.nucleus, '13C')
    if isempty(params.n_timepoints)
        params.n_timepoints = 120;
    end
    if isempty(params.default_tr)
        default_TR = 1.0; % seconds
    else
        default_TR = params.default_tr;
    end
    spectral_width_ppm = 100;
    gamma = 10.7084e6; % Hz/T for 13C
    if isempty(params.ConcStandardAmp)
        ConcStandardAmp = 1e-2; 
    else
        ConcStandardAmp = params.ConcStandardAmp; 
    end 
else % 2H
    if isempty(params.n_timepoints)
        params.n_timepoints = 40;
    end
    if isempty(params.default_tr)
        default_TR = 2.0 * 60; % 2 minutes in seconds
    else
        default_TR = params.default_tr;
    end
    if isempty(params.ConcStandardAmp)
        ConcStandardAmp=1; 
    else
        ConcStandardAmp = params.ConcStandardAmp; 
    end
    spectral_width_ppm = 40;
    gamma = 6.536e6; % Hz/T for 2H
end

% Use existing generate_metabolite_signals function for HP-13C
k_true     = params.k_true;     if isempty(k_true),     k_true = 0.05;     end
r1_pyr     = params.r1_pyr;     if isempty(r1_pyr),     r1_pyr = 1/30;     end
r1_bic     = params.r1_bic;     if isempty(r1_bic),     r1_bic = 1/20;     end
pyr_0      = params.pyr_0;      if isempty(pyr_0),      pyr_0 = 1.0;       end
FlipAngle  = params.FlipAngle;  if isempty(FlipAngle),  FlipAngle = 15;    end

%DMI parameters
dmimodel_a1   = params.a1;   if isempty(dmimodel_a1),   dmimodel_a1 = 17.0;   end
dmimodel_kp   = params.kp;   if isempty(dmimodel_kp),   dmimodel_kp = 0.0416; end
dmimodel_v    = params.v;    if isempty(dmimodel_v),    dmimodel_v  = 0.72;   end
dmimodel_kg   = params.kg;   if isempty(dmimodel_kg),   dmimodel_kg = 0.141;  end
dmimodel_Vmax = params.Vmax; if isempty(dmimodel_Vmax), dmimodel_Vmax = 0.985; end
dmimodel_kl   = params.kl;   if isempty(dmimodel_kl),   dmimodel_kl = 0.124;  end


%% Calculate spectroscopic parameters
base_freq = gamma * params.magnetic_field; % Hz
spectral_width_hz = spectral_width_ppm * base_freq / 1e6;
dwell_time = 1 / spectral_width_hz; % seconds
time_axis = (0:params.n_points-1) * dwell_time;

%% Generate kinetic timecourses
time_points = (0:params.n_timepoints-1) * default_TR;

if strcmp(params.nucleus, '13C')
    flip_angles_pyr = deg2rad(FlipAngle) * ones(1, params.n_timepoints);
    flip_angles_bic = deg2rad(FlipAngle) * ones(1, params.n_timepoints);
    
    % Generate clean kinetic timecourses (no noise yet)
    amplitudes_matrix = zeros(6, params.n_timepoints); % 6 metabolites (added UreaRef)


    [pyr_signal, bic_signal, ~, ~] = generate_metabolite_signals(...
        k_true, r1_pyr, r1_bic, pyr_0, default_TR, params.n_timepoints, ...
        flip_angles_pyr, flip_angles_bic, 0, 'gaussian', 0);

    amplitudes_matrix(1, :) = pyr_signal; % Pyruvate
    amplitudes_matrix(2, :) = bic_signal; % Bicarbonate
    amplitudes_matrix(3, :) = bic_signal; % "Lactate"
    amplitudes_matrix(4, :) = 0.05 * pyr_signal; % fake PyrHydrate
    amplitudes_matrix(5, :) = 0.02 * bic_signal; %fake alanine
    amplitudes_matrix(6, :) = ConcStandardAmp; %5e-2; % UreaRef - CONSTANT concentration reference


    metabolite_names = {'1-13C Pyruvate', '1-13C Bicarbonate', '1-13C Lactate', ...
        '1-13C PyrHydrate', '1-13C Alanine', '1-13C UreaRef'};
    dominant_amp = max(amplitudes_matrix(1, :));  % pyruvate max amplitude (dominant peak; matches Fig 7)
else % 2H
    % Use existing generate_dmi_signals approach for DMI
    %dmi_params_true = struct('a1', 17.0, 'kp', 0.0416, 'v', 0.72, ...
    %                        'kg', 0.141, 'Vmax', 0.985, 'kl', 0.124);
    dmi_params_true = struct( ...
        'a1', dmimodel_a1, ...
        'kp', dmimodel_kp, ...
        'v',  dmimodel_v, ...
        'kg', dmimodel_kg, ...
        'Vmax', dmimodel_Vmax, ...
        'kl', dmimodel_kl ...
        );
    dmi_params_fixed = struct('f', 0.6, 'Km', 2.8);
    
    params_vec = [dmi_params_true.a1, dmi_params_true.kp, dmi_params_true.v, ...
        dmi_params_true.kg, dmi_params_true.Vmax, dmi_params_true.kl];
    
    [glucose_conc, lactate_conc] = solve_dmi_model(params_vec, time_points/60, dmi_params_fixed);
    
    % Normalize and create amplitude matrix
    max_signal = max([glucose_conc; lactate_conc]);
    amplitudes_matrix = zeros(6, params.n_timepoints);
    amplitudes_matrix(1, :) = 10; % HDO (assumed constant, large. I knoPK_3T_2H_6Peakspk_structw it's not really constant but it's not in Felix's model and I don't want to model it here!)
    amplitudes_matrix(2, :) = glucose_conc' / max_signal * 2; % Glucose (2 MR-visible 2H, matched to lactate so signal is proportional to concentration; was *3, which distorted glucose:lactate and biased Vmax)
    amplitudes_matrix(3, :) = 0.5; % Glx (constant, small)
    amplitudes_matrix(4, :) = lactate_conc' / max_signal * 2; % Lactate
    amplitudes_matrix(5, :) = 0.2; % Lipid (constant, small)
    amplitudes_matrix(6, :) = ConcStandardAmp; %1.0; % TSPPhantom - CONSTANT concentration reference at 0 ppm

    metabolite_names = {'D2O', 'Glucose', 'Glx', 'Lactate', 'Lipid', 'CDCl3'};
    dominant_amp = max(amplitudes_matrix(2, :));  % glucose max amplitude (dominant substrate peak; matches Fig 7)
end

%% Get prior knowledge for chemical shifts and linewidths
if strcmp(params.nucleus, '13C')
    pk_struct = PK_3T_1C_13C_6Peaks(params.ppm_offset, 1);
else
    pk_struct = PK_3T_2H_6Peaks(params.ppm_offset);
end

%% Generate clean FIDs for each timepoint
clean_fids = zeros(params.n_points, params.n_timepoints);

fprintf('Generating %d timepoints of %s spectra...\n', params.n_timepoints, params.nucleus);

% Add the fixed extra linewidth ONCE, before the timecourse loop. (Previously this
% lived inside the per-timepoint loop and mutated the shared pk_struct every
% iteration, so fixed_linewidth ACCUMULATED -- linewidths drifted wider across the
% acquisition, e.g. pyruvate 50->62 Hz over 120 timepoints. Linewidth is a fixed
% T2* property and must be constant in time; the drift biased both integration
% (windows sized to the final, widest snapshot) and AMARES (bounds exceeded late).)
if params.flag_adjust_linewidths
    for jx = 1:size(pk_struct.values.IVCellArray,1)
        pk_struct.values.IVCellArray{jx,3} = pk_struct.values.IVCellArray{jx,3} + params.fixed_linewidth;
    end
end

for t = 1:params.n_timepoints
    if mod(t, 10) == 0
        fprintf('  Timepoint %d/%d\n', t, params.n_timepoints);
    end

    % Generate clean FID
    current_amplitudes = amplitudes_matrix(:, t);
    clean_fid = make_model_fid(current_amplitudes, metabolite_names, pk_struct, ...
        time_axis, base_freq);
    clean_fids(:, t) = clean_fid;
end

%% Add noise AFTER generating all clean FIDs
% The reference against which noise_level is scaled is constant across the timecourse
% (noise power is constant, not time-varying).
switch lower(char(params.noise_reference))
    case 'peakamp'
        % Referenced to the maximum substrate peak amplitude, as in the amplitude-domain
        % simulation (Fig 7). The added complex FID noise then has RMS magnitude
        % = noise_level * dominant peak amplitude, so the two noise axes are comparable.
        max_signal_power = dominant_amp^2;
        ref_desc = 'max peak amplitude';
    case 'rmsfid'
        % Legacy: referenced to the RMS of the FID at the max-signal timepoint.
        max_signal_power = max(mean(abs(clean_fids).^2, 1));
        ref_desc = 'RMS FID';
    otherwise
        error('generate_fake_fids:noise_reference', ...
            'unknown noise_reference ''%s'' (use ''peakamp'' or ''rmsfid'')', char(params.noise_reference));
end

fprintf('Adding noise (level=%.4f, referenced to %s)...\n', params.noise_level, ref_desc);
noisy_fids = add_white_noise_timecourse(clean_fids, params.noise_level, max_signal_power);

%% Package results
fake_spectra = struct();
fake_spectra.nucleus = params.nucleus;
fake_spectra.parameters = params;
fake_spectra.metabolite_names = metabolite_names;
fake_spectra.true_amplitudes = amplitudes_matrix;
fake_spectra.time_points = time_points;
fake_spectra.time_axis = time_axis;
fake_spectra.dwell_time = dwell_time;
fake_spectra.base_freq = base_freq;
fake_spectra.clean_fids = clean_fids;
fake_spectra.noisy_fids = noisy_fids;
fake_spectra.prior_knowledge = pk_struct;
fake_spectra.n_samples = size(fake_spectra.time_axis,2); 

fake_spectra.freqAxis= ((-fake_spectra.n_samples/2):(fake_spectra.n_samples/2-1)).'/(fake_spectra.dwell_time*fake_spectra.n_samples);
fake_spectra.ppmAxis=fake_spectra.freqAxis /(fake_spectra.base_freq*1e-6); 

fprintf('Generated %d %s spectra successfully!\n', params.n_timepoints, params.nucleus);

end

function clean_fid = make_model_fid(amplitudes, metabolite_names, pk_struct, time_axis, base_freq)
%MAKE_MODEL_FID Generate a single clean model FID from metabolite parameters
%
% INPUTS:
%   amplitudes       - Vector of peak amplitudes
%   metabolite_names - Cell array of metabolite names
%   pk_struct        - Prior knowledge structure
%   time_axis        - Time points for FID (seconds)
%   base_freq        - Base frequency (Hz)

n_peaks = length(metabolite_names);
n_points = length(time_axis);
clean_fid = zeros(n_points, 1);

for i = 1:n_peaks
    if amplitudes(i) <= 0
        continue; % Skip zero amplitude peaks
    end
    
    % Find metabolite in prior knowledge
    pk_idx = find(strcmp(pk_struct.values.IVCellArray(:,1), metabolite_names{i}), 1);
    
    if isempty(pk_idx)
        warning('Metabolite %s not found in prior knowledge', metabolite_names{i});
        continue;
    end
    
    % Get parameters from prior knowledge
    chem_shift_ppm = pk_struct.values.IVCellArray{pk_idx, 2}; % ppm
    linewidth_hz = pk_struct.values.IVCellArray{pk_idx, 3}; % Hz
    phase_deg = pk_struct.values.IVCellArray{pk_idx, 5}; % degrees
    
    % Convert to appropriate units
    freq_hz = chem_shift_ppm * base_freq / 1e6; % Hz
    phase_rad = deg2rad(phase_deg);
    decay_rate = pi * linewidth_hz; % rad/s
    
    % Generate complex exponential signal
    signal = amplitudes(i) * exp(1i * phase_rad) * ...
        exp(1i * 2 * pi * freq_hz * time_axis(:)) .* ...
        exp(-decay_rate * time_axis(:));
    
    clean_fid = clean_fid + signal;
end

end

function noisy_fids = add_white_noise_timecourse(clean_fids, noise_level, reference_signal_power)
%ADD_WHITE_NOISE_TIMECOURSE Add constant white Gaussian noise to entire timecourse
%
% INPUTS:
%   clean_fids              - Clean complex FIDs [n_points x n_timepoints]
%   noise_level             - Noise level as fraction of reference signal
%   reference_signal_power  - Reference signal power (typically max across timecourse)
%
% This function adds noise with CONSTANT power across all timepoints,
% which is physically realistic. The noise power is calculated relative
% to a reference signal (typically the maximum signal).

if noise_level <= 0
    noisy_fids = clean_fids;
    return;
end

[n_points, n_timepoints] = size(clean_fids);

% Calculate noise power based on reference signal
noise_power = noise_level^2 * reference_signal_power;

% Standard deviation for real and imaginary parts
% (divide by 2 because total power splits equally between Re and Im)
noise_std = sqrt(noise_power / 2);

% Generate independent Gaussian noise for all timepoints at once (vectorized)
noise_real = noise_std * randn(n_points, n_timepoints);
noise_imag = noise_std * randn(n_points, n_timepoints);

% Add complex noise
noisy_fids = clean_fids + noise_real + 1i * noise_imag;

end


function [glucose_conc, lactate_conc] = solve_dmi_model(params, time_points, params_fixed)
%SOLVE_DMI_MODEL Solve DMI differential equations

a1 = params(1); kp = params(2); v = params(3);
kg = params(4); Vmax = params(5); kl = params(6);
f = params_fixed.f; Km = params_fixed.Km;

y0 = [0; 0];
ode_fun = @(t, y) dmi_odes(t, y, a1, kp, v, kg, Vmax, kl, f, Km);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, y] = ode45(ode_fun, time_points, y0, options);

glucose_conc = y(:, 1);
lactate_conc = y(:, 2);
end

function dydt = dmi_odes(t, y, a1, kp, v, kg, Vmax, kl, f, Km)
%DMI_ODES DMI differential equations
% NOTE ON PARAMETERISATION. This is Kreis et al. (Radiology 2020, Eqs S5-S6),
%   dG*/dt = kg*(a1*exp(-kp*t) - G*/v) - f*Vmax*G*/(f*v*(Km + G*/(f*v)))
%   dL*/dt = f*Vmax*G*/(f*v*(Km + G*/(f*v))) - kl*L*
% written in terms of the tissue-averaged pool: here a1, kg and Vmax stand for
% v*a1, kg/v and Vmax/v of Kreis et al. (Vmax per unit extracellular volume, so the
% tissue-averaged consumption is f*v*Vglc). The two forms give identical curves under
% that mapping. Note G_star_e below is the TOTAL extracellular glucose [Ge] = G*/(f*v)
% of Eq S3, not the labelled [G*e] = G*/v of Eq S2. The default parameter values are
% the Table E2 estimates, used directly in this parameterisation.

G_star = y(1); L_star = y(2);
G_star_p = a1 * exp(-kp * t);
G_star_e = G_star / (f * v);
Vglc = Vmax * G_star_e / (Km + G_star_e);

dG_dt = kg * (G_star_p - G_star) - f * v * Vglc;
dL_dt = f * v * Vglc - kl * L_star;

dydt = [dG_dt; dL_dt];
end

function outStruct = PK_3T_1C_13C_6Peaks(ppm_off, sgn)
%PK_3T_1C_13C_5PEAKS Prior knowledge for 13C metabolites (now 6 peaks with UreaRef)

fields.Bounds = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
values.boundsCellArray = {
    '1-13C Lactate',    [11.1,13.1],         [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C PyrHydrate', [7.2,9.2],           [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C Alanine',    [4.6,6.6],           [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C Pyruvate',   [-1,1],            [0,120],    [0,inf],     [0,360], [],               [];
    '1-13C Bicarbonate',[-11.1,-9.1],        [0,20],     [0,inf],     [0,360], [],               [];
    '1-13C UreaRef',    [4.0,5.1],           [5,15],     [0,inf],     [0,360], [],               [];
    };

fields.IV = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
values.IVCellArray = {
    '1-13C Lactate',      12.1,      10,          0.1,           0;
    '1-13C PyrHydrate',   8.2,       10,          0.1,           0;
    '1-13C Alanine',      5.6,       10,          0.1,           0;
    '1-13C Pyruvate',     0,         50,          10,            0;
    '1-13C Bicarbonate',  -10.1,     10,          0.1,           0;
    '1-13C UreaRef',      4.56,      10,          1.0,           0;
    };

fields.PK = {'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak'};
values.PKCellArray = {
    '1-13C Lactate',    [],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C PyrHydrate', [],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C Alanine',    [],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C Pyruvate',   [],          [],               [],               [],            [],            1,         [],         [],                 1;
    '1-13C Bicarbonate',[],          [],               [],               [],            [],            1,         [],         [],                 0;
    '1-13C UreaRef',    [],          [],               [],               [],            [],            1,         [],         [],                 0;
    };

if ((abs(ppm_off)>1d-10) || (sgn<0))
    for l=1:size(values.boundsCellArray,1)
        if sgn<0
            values.boundsCellArray{l,2} = sgn*(values.boundsCellArray{l,2}([2,1]) - ppm_off);
        else
            values.boundsCellArray{l,2} = values.boundsCellArray{l,2} - ppm_off;
        end
        values.IVCellArray{l,2} = sgn*(values.IVCellArray{l,2} - ppm_off);
    end
end

outStruct.fields = fields;
outStruct.values = values;
end

function outStruct = PK_3T_2H_6Peaks(ppm_offset)
%PK_3T_2H_5PEAKS Prior knowledge for 2H metabolites (now 6 peaks with TSPPhantom)

fields.Bounds = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
values.boundsCellArray = {
    'D2O',            [4.5,4.9],   [0,inf],     [0,inf],     [0,360], [],               [];
    'Glucose',        [3.6,4.0],   [9,11],      [0,inf],     [0,360], [],               [];
    'Glx',            [2.1,2.5],   [0,100],     [0,inf],     [0,360], [],               [];
    'Lactate',        [0.7,1.35],  [0,100],     [0,inf],     [0,360], [],               [];
    'Lipid',          [0,0.5],     [0,100],     [0,inf],     [0,360], [],               [];
    'CDCl3',          [7.1,7.3],   [5,15],      [0,inf],     [0,360], [],               [];
    };

fields.IV = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
values.IVCellArray = {
    'D2O',            4.7,         8,          0.5,          0;
    'Glucose',        3.8,         10,         0.5,           0;
    'Glx',            2.3,         8,          0.5,           0;
    'Lactate',        1.0,         8,          0.5,         0;
    'Lipid',          0,           8,          0.5,         0;
    'CDCl3',          7.2,           10,       0.5,         0;

%     'D2O',            4.7,         8,          10,          0;
%     'Glucose',        3.8,         10,         3,           0;
%     'Glx',            2.3,         8,          1,           0;
%     'Lactate',        1.0,         8,          0.5,         0;
%     'Lipid',          0,           8,          0.5,         0;
%     'CDCl3',          7.2,           10,         1.0,         0;
    };

fields.PK = {'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak'};
values.PKCellArray = {
    'D2O',             [],          [],               [],               [],           [],            1,         [],         [],                 1;
    'Glucose',         [],          [],               [],               [],           [],            1,         [],         [],                 0;
    'Glx',             [],          [],               [],               [],           [],            1,         [],         [],                 0;
    'Lactate',         [],          [],               [],               [],           [],            1,         [],         [],                 0;
    'Lipid',           [],          [],               [],               [],           [],            [],        [],         [],                 0;
    'CDCl3',           [],          [],               [],               [],           [],            1,         [],         [],                 0;
    };

for l=1:size(values.boundsCellArray,1)
    values.boundsCellArray{l,2} = values.boundsCellArray{l,2} + ppm_offset;
    values.IVCellArray{l,2} = values.IVCellArray{l,2} + ppm_offset;
end

outStruct.fields = fields;
outStruct.values = values;
end


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

normalise_sd = 1;

if strcmp(noise_type, 'gaussian')
    % Gaussian noise of mean zero and standard deviation (noise_level * max_pyr)
    pyr_noise = noise_level * max_pyr * randn(1, n_timepoints);
    bic_noise = noise_level * max_bic * randn(1, n_timepoints);
    pyr_signal = pyr_signal + pyr_noise;
    bic_signal = bic_signal + bic_noise;
    
elseif strcmp(noise_type, 'chi-squared')
    % Non-central chi-squared noise (just to see FTW / FTLolz what happens
    % if you put this in to the FID...)
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

function [glucose_noisy, lactate_noisy] = generate_dmi_signals(...
    glucose_true, lactate_true, noise_level, noise_type, censoring_threshold)
% Generate noisy DMI signals using the same noise distributions as HP-13C

global rate_V
global scale_delta
global lognormal_sigma
global lognormal_meanlog
global gamma_shape
global gamma_rate

n_timepoints = length(glucose_true);
max_glucose = max(glucose_true);
max_lactate = max(lactate_true);

normalise_sd = 1;

if strcmp(noise_type, 'gaussian')
    % Gaussian noise
    glucose_noise = noise_level * max_glucose * randn(1, n_timepoints);
    lactate_noise = noise_level * max_lactate * randn(1, n_timepoints);
    glucose_noisy = glucose_true + glucose_noise';
    lactate_noisy = lactate_true + lactate_noise';
    
elseif strcmp(noise_type, 'chi-squared')
    % Non-central chi-squared noise
    gn = ncx2rnd(rate_V, scale_delta, 1, n_timepoints);
    ln = ncx2rnd(rate_V, scale_delta, 1, n_timepoints);
    
    if normalise_sd
        std_gn = sqrt(2*(rate_V + 2*scale_delta));
        std_ln = sqrt(2*(rate_V + 2*scale_delta));
    else
        std_gn = 1;
        std_ln = 1;
    end
    
    glucose_noise = noise_level * max_glucose * (gn / std_gn);
    lactate_noise = noise_level * max_lactate * (ln / std_ln);
    
    glucose_noisy = glucose_true + glucose_noise';
    lactate_noisy = lactate_true + lactate_noise';
    
elseif strcmp(noise_type, 'lognormal')
    % Lognormal noise
    glucose_noise_raw = lognrnd(lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);
    lactate_noise_raw = lognrnd(lognormal_meanlog, lognormal_sigma, [1 n_timepoints]);
    
    if normalise_sd
        std_glucose_noise = std(glucose_noise_raw);
        std_lactate_noise = std(lactate_noise_raw);
    else
        std_glucose_noise = 1;
        std_lactate_noise = 1;
    end
    
    glucose_noise = noise_level * max_glucose * (glucose_noise_raw / std_glucose_noise);
    lactate_noise = noise_level * max_lactate * (lactate_noise_raw / std_lactate_noise);
    
    glucose_noisy = glucose_true + glucose_noise';
    lactate_noisy = lactate_true + lactate_noise';
    
elseif strcmp(noise_type, 'gamma') || strcmp(noise_type, 'censored')
    % Gamma noise
    glucose_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);
    lactate_noise_raw = gamrnd(gamma_shape, 1/gamma_rate, 1, n_timepoints);
    
    if normalise_sd
        std_glucose_noise = std(glucose_noise_raw);
        std_lactate_noise = std(lactate_noise_raw);
    else
        std_glucose_noise = 1;
        std_lactate_noise = 1;
    end
    
    glucose_noise = noise_level * max_glucose * (glucose_noise_raw / std_glucose_noise);
    lactate_noise = noise_level * max_lactate * (lactate_noise_raw / std_lactate_noise);
    
    glucose_noisy = glucose_true + glucose_noise';
    lactate_noisy = lactate_true + lactate_noise';
end

% Ensure non-negative signals
glucose_noisy = max(glucose_noisy, 0);
lactate_noisy = max(lactate_noisy, 0);

% Apply censoring if requested
if strcmp(noise_type, 'censored')
    glucose_threshold = censoring_threshold * max_glucose;
    lactate_threshold = censoring_threshold * max_lactate;
    
    glucose_noisy(glucose_noisy < glucose_threshold) = 0;
    lactate_noisy(lactate_noisy < lactate_threshold) = 0;
end

end
