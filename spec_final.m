clear; close all;

% this code puts down the parameters of the AR(2) process 
% we use to simulate the small sample bias and standard errors
% we select it in spec_findAR2coef.m

load('CleanData/US_data.mat');

phi1v = 1.42082559207046;
phi2v = -0.420883510141803;

phi1r = 0.289111306354948;
phi2r = -0.281919479948611;

% find residuals
demean = @(x) (x - mean(x));
horizon = 1;

Y_Beta = v_tilde(2 + 1:end, 2);
X2_Beta = [v_tilde(2:end - 1, 2) v_tilde(1:end - 1 - 1, 2)];
E_v = demean(Y_Beta) - demean(X2_Beta) * [phi1v, phi2v]';

Y_Beta = return_debt_tilde(2 + 1:end, 2);
X2_Beta = [v_tilde(2:end - 1, 2) v_tilde(1:end - 1 - 1, 2)];
E_r = demean(Y_Beta) - demean(X2_Beta) * [phi1r, phi2r]';

E_sim = [E_v E_r];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len_long = 1e5;

rng(1, 'twister');
E_rnd = datasample(E_sim, len_long, 'Replace', true);
E_v = E_rnd(:, 1);
E_r = E_rnd(:, 2);

v_long = zeros(len_long, 1);
r_long = zeros(len_long, 1);
s_long = zeros(len_long, 1);

for t = 1:len_long
    v_long(t + 2) = [v_long(t + 2 - 1) v_long(t + 2 - 2)] * [phi1v, phi2v]' + E_v(t);
    r_long(t + 2) = [v_long(t + 2 - 1) v_long(t + 2 - 2)] * [phi1r, phi2r]' + E_r(t);
    s_long(t + 2) = v_long(t + 1) - v_long(t + 2) + r_long(t + 2);
end

save('CleanData/simul.mat', 'phi1r', 'phi2r', 'phi1v', 'phi2v', ...
    'E_sim', 'v_long', 'r_long', 's_long', 'len_long');
