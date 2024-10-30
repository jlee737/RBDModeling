<<<<<<< HEAD
function [Bound_NBs, tot_bound, Kd] = PosSelection(Targets, unbound_NB, pos_avg_Kd, pos_std_Kd)
    
    %% Function Information
    % Function models the equilibrium biopanning process in positive selection of NBs

    % Assumptions:
        % Enough time is taken for equilibrium to be reached between NBs and Targets
        % The number of targets is non-limiting (Remains constant)
        % The Kd distribution for NBs follows a normal distribution and can be modeled
        % Targets conjugated to bead are the equivalent of free-floating targets in solution

    % Inputs:
        % Targets - Number of Targets bound to all beads
        % Tot_NB - Number of NB-Ribosome complexes in total
        % avg_Kd - Average dissociation constant for NBs
        % std_Kd - Standard Deviation for NBs 

    % Outputs:
        % Bound_NBs - Array containing prevalence of NBs and associated Kd
        % tot_bound - Total number of NBs recovered
        % Kd - Values of Kd for all nanobodies

    %% Code
    
    % Define constants
    vol = (225 + 6.25) * 10^-6; % L, Volume of IVTT and beads
    Kd_bins = 10^7; % Number of possible Kd values
    
    % Convert to molarity
    c.Targets = Targets / (6.022 * 10^23) / vol; % Molar
    c.Tot_NB = unbound_NB / (6.022 * 10^23) / vol; % Molar
    
    % Define range for Kd values
    Kd = linspace(1e-20, 100 * pos_avg_Kd, Kd_bins);

    % Calculate number of nanobodies per bin
    NB_per_Kd = c.Tot_NB * (normpdf(Kd, pos_avg_Kd, pos_std_Kd) / sum(normpdf(Kd, pos_avg_Kd, pos_std_Kd))); % Molar
    
    % Calculate equilibrium to determine bound molecules
    Bound_NBs = NB_per_Kd .* c.Targets ./ (Kd + c.Targets); % Molar
    Bound_molecules = Bound_NBs .* 6.022 * 10^23 * vol;

    tot_bound = sum(Bound_molecules);
    
    % Plot Kd distribution
    figure;
    
    plot(Kd, NB_per_Kd * 6.022 * 10^23 * vol);
    hold on
    plot(Kd, Bound_molecules);
    hold off

    legend("All NB", 'Bound NB')

    title('Kd Distribution');
    xlabel('Kd (M)');
    ylabel('Number of Binders');

end
=======
function [Bound_NBs, tot_bound, Kd] = PosSelection(Targets, unbound_NB, pos_avg_Kd, pos_std_Kd)
    
    %% Function Information
    % Function models the equilibrium biopanning process in positive selection of NBs

    % Assumptions:
        % Enough time is taken for equilibrium to be reached between NBs and Targets
        % The number of targets is non-limiting (Remains constant)
        % The Kd distribution for NBs follows a normal distribution and can be modeled
        % Targets conjugated to bead are the equivalent of free-floating targets in solution

    % Inputs:
        % Targets - Number of Targets bound to all beads
        % Tot_NB - Number of NB-Ribosome complexes in total
        % avg_Kd - Average dissociation constant for NBs
        % std_Kd - Standard Deviation for NBs 

    % Outputs:
        % Bound_NBs - Array containing prevalence of NBs and associated Kd
        % tot_bound - Total number of NBs recovered
        % Kd - Values of Kd for all nanobodies

    %% Code
    
    % Define constants
    vol = (225 + 6.25) * 10^-6; % L, Volume of IVTT and beads
    Kd_bins = 10^8; % Number of possible Kd values
    
    % Convert to molarity
    c.Targets = Targets / (6.022 * 10^23) / vol; % Molar
    c.Tot_NB = unbound_NB / (6.022 * 10^23) / vol; % Molar
    
    % Define range for Kd values
    Kd = linspace(1e-20, 100 * pos_avg_Kd, Kd_bins);

    % Calculate number of nanobodies per bin
    NB_per_Kd = c.Tot_NB * (normpdf(Kd, pos_avg_Kd, pos_std_Kd) / sum(normpdf(Kd, pos_avg_Kd, pos_std_Kd))); % Molar
    
    % Calculate equilibrium to determine bound molecules
    Bound_NBs = NB_per_Kd .* c.Targets ./ (Kd + c.Targets); % Molar
    Bound_molecules = Bound_NBs .* 6.022 * 10^23 * vol;

    tot_bound = sum(Bound_molecules);
    
    % Plot Kd distribution
    figure;
    
    plot(Kd, NB_per_Kd * 6.022 * 10^23 * vol);
    hold on
    plot(Kd, Bound_molecules);
    hold off

    legend("All NB", 'Bound NB')

    title('Kd Distribution');
    xlabel('Kd (M)');
    ylabel('Number of Binders');

end
>>>>>>> 752a4959d60e9087716e6650951d4fc80afb62c3
