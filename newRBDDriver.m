close all;

%% Sensitivity Analysis:

% Metrics:
%   1. How many unique NBs
Unique = [];
%   2. Total NBs pulled down
Total = [];
%   2. Fraction enriched per bin
EnrichedFrac = [];
%   3. Average kd from pulldown
AvgKd = [];

% Parameters to perturb:
%   1. Kd_bins
%   2. Initial distribution
%   3. Number of targets
%       a. Positive selection
%       b. Negative selection
%   4. Number of initial nbs
%   5. Detection threshold
%   6. Number of rounds

vol = (225 + 6.25) * 10^-6; % L
avo = 6.022 * 10^23;

% Perturb:

Kd_bins = 10^5;
init_dist.a = 12;
init_dist.b = 1.2;

neg.a = 12;
neg.b = 1.2;

neg.targets = 10^13 / avo / vol;
pos.targets = 10^13/ avo / vol;

init_nb = 10^11;

threshold = 1000000;

rounds = 5;

% Kd_bins = [10^4, 10^5, 10^6, 10^7];
% 
% init_dist.a = [5, 10, 15, 20];
% init_dist.b = [1.1, 1.2, 1.3, 1.4];
% 
% neg.Targets = [10^11, 10^12, 10^13, 10^14] / avo / vol;
% targets = [10^11, 10^12, 10^13, 10^14] / avo / vol;
% 
% init_nb = [10^11, 10^12, 10^13, 10^14];
% 
% threshold = [1000, 10000, 100000];
% 
% rounds = [1, 2, 3, 4, 5];

%% Round 1 selection:

%% Initialize Parameters...
wash_vol = (225+6.25) * 10^-6;
wash_time = 30; 
k_on = 10^6;
num_wash = 4;

dilution = vol / (vol + wash_vol);

init_kds = linspace(1e-20, 0.001, Kd_bins);

%% Generate initial distribution pre-R1
r1dist = init_nb .* betapdf(init_kds, init_dist.a, init_dist.b) / sum(betapdf(init_kds, init_dist.a, init_dist.b)); % Number of nanobodies per affinity

Mr1dist = r1dist / avo / vol;

%% Round 1

% Negative Selection:

negNBpBin = init_nb .* (betapdf(init_kds, neg.a, neg.b) / sum(betapdf(init_kds, neg.a, neg.b)));

negBound = negNBpBin .* neg.targets ./ (init_kds + neg.targets);

neg_tot_bound = sum(negBound);
negEnrich = neg_tot_bound / init_nb; 

% Wash

spec_binding = 1 - 0.9 / (init_kds(end)) .* init_kds;
non_spec_binding = 0.1;

negBound = negBound.* (non_spec_binding * dilution ^ num_wash + spec_binding .* exp(-k_on .* init_kds * wash_time));
mNegBound = negBound;

neg_tot_bound = sum(mNegBound);

% Delete NBs lost to negative selection
bound = round(neg_tot_bound);

lostNBindex = randi([1, length(r1dist)], 1, bound);

lostNB = zeros(1,length(r1dist));
for i = 1:length(lostNBindex)
    lostNB(lostNBindex(i)) = lostNB(lostNBindex(i)) + 1;
end

r1dist = r1dist - lostNB;

% Positive Selection
Mr1bound = Mr1dist .* pos.targets ./ (init_kds + pos.targets); 

init_avgKd = sum(Mr1bound .* init_kds) / sum(Mr1bound)

prewash = sum(Mr1bound .* vol * avo);

% Washing
Bound = Mr1bound.* (non_spec_binding * dilution ^ num_wash + spec_binding .* exp(-k_on .* init_kds * wash_time));
mBound = Bound .* avo * (vol + wash_vol);

totmBound = sum(mBound);

[Kd, NB_per_bin, Tot_unique_NBs, Tot_NBs, new_dist] = AssignBins(init_kds, mBound);

% Post R1 Amplification
rand_amp = randi([10^3, 10^7], 1, Tot_unique_NBs);

r1amp = NB_per_bin .* rand_amp;
csum = sum(r1amp);
adj = init_nb / csum;

r1amp = r1amp * adj;

% Threshold for RT/PCR
overmin = r1amp >= threshold;

length(r1amp);
rAmp = r1amp(overmin);

length(rAmp);
Kd = Kd(overmin);

% if length(rAmp) == 0
%     disp("No binders")
%     return;
% end

Unique(end+1) = length(rAmp);
Total(end+1) = sum(rAmp);
AvgKd(end+1) = sum(rAmp .* Kd) / sum(rAmp);



%% Subsequent Rounds

for i = 1:(rounds - 1)

    [outputbins, outputNbpBin, uniqueNBs, totBound] = SelectionRound(Kd, rAmp, neg, init_nb, pos, threshold);
    
    % To next round
    Kd = outputbins;
    rAmp = outputNbpBin;
    
    % Metrics
    Unique(end+1) = uniqueNBs;
    Total(end+1) = totBound;
    AvgKd(end+1) = sum(rAmp .* Kd) / sum(rAmp);


end


avg_sens_kd(end+1) = AvgKd(end);
unique_sens_nb(end+1) = Unique(end);

