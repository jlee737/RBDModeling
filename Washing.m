function [Bound_NBs, tot_bound, Bound_molecules] = Washing(Init_NB, Tot_NB, Kd, density_func, wash_vol, num_wash)
    
    %% Function Information
    % Function models the process of washing both specific and non-specific binders off beads

    % Assumptions:
        % Specific and non-specific binding are independent
        % Non-specific binders exist as a proportion of specific binders in each Kd "bin"
        % Washing can be described by a first order exponential decay process
    
    % Inputs:
        % Init_NB - Initial number of NBs put into RBD
        % Tot_NB -  Number of NB-Ribosome complexes inputted into wash
        % Kd - Values of Kd for all nanobodies
        % density_func - Distribution describing Kd values for input NBs
        % wash_vol - Volume of wash buffer used per wash (L)
        % num_wash - Number of washes

    % Outputs:
        % Bound_NBs - Array containing prevalence of NBs and associated Kd
        % tot_bound - Total number of NBs recovered

    %% Code:
    
    % Define constants
    vol = (225 + 6.25) * 10^-6; % L, Volume of IVTT and beads
    wash_time = 30; % S
    k_on = 10^6; % 1/S Could potentially be measured via BLI?

    % Calculate wash function parameters
    dilution = vol / (vol + wash_vol); % Dilution factor during wash
    
    binding_frac = Tot_NB / Init_NB;
    spec_binding = 1 - 0.9 / (Kd(end)) .* Kd; % Specific binding fraction as function of Kd is unknown...
    non_spec_binding = 0.1; % Assumption

    % Calculate new distribution after washing
    Bound_NBs = density_func .* (non_spec_binding * dilution ^ num_wash + spec_binding .* exp(-k_on .*Kd * wash_time));
    Bound_molecules = Bound_NBs .* 6.022 * 10^23 * (vol + wash_vol);

    tot_bound = sum(Bound_molecules);

% end