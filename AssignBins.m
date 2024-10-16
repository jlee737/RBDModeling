function [Kd, NB_per_bin, Tot_unique_NBs, Tot_NBs, new_dist] = AssignBins(Kd, Nb)
    
    %% Function Information
    % Function assigns different NBs with same affinities to different bins

    % Inputs:
        % Kd - Array of possible Kd values in model
        % Nb - Array of discrete nanobodies per Kd 
   
    % Outputs:
        % Kd - Array of possible Kd values with discrete bins per NB
        % NB_per_bin - Number of NBs per new bin
        % Tot_unique_NBs - Number of bins (unique binding NBs)
        % Tot_NBs - Total number of NBs counting replicates
        % new_dist - List containing Kd for each NB to be used to generate
        % new disib
    %% Code
    new_Kd = [];
    
    for i = 1:length(Kd)
        for j = 1:Nb(i) 
            new_Kd(end + 1) = Kd(i); 
        end
    end

    Kd = new_Kd;
    NB_per_bin = ones(1,length(Kd));
    Tot_unique_NBs = length(Kd);
    Tot_NBs = sum(NB_per_bin);

    new_dist = [];

    for i = 1:length(Kd)
        for j = 1:Nb(i)
        new_dist(end+1) = Kd(i);
        end
    end

end
