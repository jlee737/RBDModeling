% Assumes that the number of targets is not limiting

P = 10^11; % Molecules

c.P = P / (6.022 * 10^23) / ((225 + 6.25) * 10^-6); % M

Tot_NB = 10^11; % Molecules

c.NB = Tot_NB / (6.022 * 10^23) / ((225 + 6.25) * 10^-6); % M

avg = 10000 * 10^-9; % M
std = avg * 8;

Kd = linspace(0, 100 * avg, 10^8);

NB_per_Kd = c.NB * (normpdf(Kd, avg, std)) / sum(normpdf(Kd, avg, std)); % Molarity

Bound = NB_per_Kd .* (c.P ./ (Kd + c.P)); % M

Bound_molecules = Bound .* 6.022 * 10^23 * ((225 + 6.25) * 10^-6);

sum(Bound_molecules)

plot(Kd, normpdf(Kd, avg, std))


