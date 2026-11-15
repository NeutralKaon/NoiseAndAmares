function [glucose_conc, lactate_conc] = solve_dmi_model(params, time_points, params_fixed)
%SOLVE_DMI_MODEL Solve the DMI differential equation model
%
% Implements equations S5 and S6 from the Radiology supplementary material
% i.e. the SI for doi: 10.1148/radiol.2019191242

    % Extract parameters
    a1 = params(1);
    kp = params(2);
    v = params(3);
    kg = params(4);
    Vmax = params(5);
    kl = params(6);
    
    f = params_fixed.f;
    Km = params_fixed.Km;
    
    % Initial conditions: [G*, L*]
    y0 = [0; 0];  % Start with zero tissue concentrations
    
    % Define the ODE system
    ode_fun = @(t, y) dmi_odes(t, y, a1, kp, v, kg, Vmax, kl, f, Km);
    
    % Solve ODE system
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [~, y] = ode45(ode_fun, time_points, y0, options);
    
    % Extract concentrations
    glucose_conc = y(:, 1);
    lactate_conc = y(:, 2);
end

function dydt = dmi_odes(t, y, a1, kp, v, kg, Vmax, kl, f, Km)
%DMI_ODES System of ODEs for DMI kinetic model
%
% NOTE ON PARAMETERISATION. This is Kreis et al. (Radiology 2020, Eqs S5-S6),
%   dG*/dt = kg*(a1*exp(-kp*t) - G*/v) - f*Vmax*G*/(f*v*(Km + G*/(f*v)))
%   dL*/dt = f*Vmax*G*/(f*v*(Km + G*/(f*v))) - kl*L*
% written in terms of the tissue-averaged pool: here a1, kg and Vmax stand for
% v*a1, kg/v and Vmax/v of Kreis et al. (Vmax per unit extracellular volume, so the
% tissue-averaged consumption is f*v*Vglc). The two forms give identical curves under
% that mapping. Note G_star_e below is the TOTAL extracellular glucose [Ge] = G*/(f*v)
% of Eq S3, not the labelled [G*e] = G*/v of Eq S2. The default parameter values are
% the Table E2 estimates, used directly in this parameterisation.

    % Current state
    G_star = y(1);  % Tissue glucose concentration
    L_star = y(2);  % Tissue lactate concentration
    
    % Plasma glucose concentration (exponential decay)
    G_star_p = a1 * exp(-kp * t);
    
    % Total extracellular glucose, [Ge] = G*/(f v)  (Eq S3)
    G_star_e = G_star / (f * v);  % Adjusted by fractional enrichment and volume
    
    % Glucose consumption rate (Michaelis-Menten kinetics)
    Vglc = Vmax * G_star_e / (Km + G_star_e);
    
    % Differential equations
    dG_dt = kg * (G_star_p - G_star) - f * v * Vglc;
    dL_dt = f * v * Vglc - kl * L_star;
    
    dydt = [dG_dt; dL_dt];
end


