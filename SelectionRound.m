% Input: 
% 1. every single kd (repeated if multiple nbs per kd)
% 2. nbs per kd

% Output:
% 1. every single kd
% 2. nbs per kd

function [outputbins, outputNbpBin, uniqueNBs, totBound] = SelectionRound(inputbins, inputnbpbin, neg, init_nb, pos, threshold)

%% Negative Selection

% Define Parameters

avo = 6.022 * 10^23;
vol = (225 + 6.25) * 10^-6; % L

neg.a = 15;
neg.b = 1.5;

% Selection

negNBpBin = init_nb * (betapdf(inputbins, neg.a, neg.b) / sum(betapdf(inputbins, neg.a, neg.b)));

negBound = negNBpBin .* neg.targets ./ (inputbins + neg.targets);

neg_tot_bound = sum(negBound) .* avo * vol;
negEnrich = neg_tot_bound / init_nb; 


%% Wash

% Define Parameters

wash_vol = (225+6.25) * 10^-6;
wash_time = 30; 
k_on = 10^6;
num_wash = 3;

dilution = vol / (vol + wash_vol);

% Binding fractions
spec_binding = 1 - 0.9 / (inputbins(end)) .* inputbins;
non_spec_binding = 0.1;

% Washing
negBound = negBound.* (non_spec_binding * dilution ^ num_wash + spec_binding .* exp(-k_on .* inputbins * wash_time));
mNegBound = negBound .* avo * (vol + wash_vol);

tot_bound = sum(mNegBound);

%% Remove negative selected nb
unboundFrac = (init_nb * avo * vol - tot_bound) / (init_nb * avo * vol);

inputnbpbin = inputnbpbin .* unboundFrac / avo / vol;

%% Positive Selection

% Define Parameters

vol = (225 + 6.25) * 10^-6; % L
avo = 6.022 * 10^23;

mBound = inputnbpbin .* pos.targets ./ (inputbins + pos.targets);
Bound = mBound .* avo * vol;

sum(Bound);

% Washing
Bound = Bound.* (non_spec_binding * dilution ^ num_wash + spec_binding .* exp(-k_on .* inputbins * wash_time));
mBound = Bound .* avo * (vol + wash_vol);

sum(Bound);

% Amplification

totBound = sum(Bound);

rand_amp = randi([10^3, 10^7], 1, length(inputbins));
rAmp = Bound .* rand_amp;
csum = sum(rAmp);
adj = init_nb / csum;

rAmp = rAmp * adj;

% Delete nanobodies below detection threshold
overmin = rAmp >= threshold;

rAmp = rAmp(overmin);
inputbins = inputbins(overmin);


outputbins = inputbins;
outputNbpBin = rAmp;
uniqueNBs = length(inputbins);




end
