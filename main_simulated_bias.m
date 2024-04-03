clear; close all;
load('CleanData/US_data.mat');

% result from data
[coeff_RS, coeff_R, coeff_S, coeff_V, mdl_sest, mdl_rest, mdl_vest] = ...
    estimate_vsr(v_tilde(:, 2), return_debt_tilde(:, 2), surplus_tilde(:, 2), 10);

coeff_R_uni = coeff_R;
coeff_S_uni = coeff_S;
coeff_V_uni = coeff_V;
coeff_RS_uni = coeff_RS;
results_mat = [coeff_R(:, 1:10); coeff_S(:, 1:10); coeff_V(:, 1:10)];

% construct cumulative series
horizon_total = 10;
r_cum = zeros(T, horizon_total);
s_cum = zeros(T, horizon_total);

r_cum(1:end, 1) = return_debt_tilde(:, 2);
s_cum(1:end, 1) = surplus_tilde(:, 2);

for horizon = 2:horizon_total
    r_cum(1 + horizon:end, horizon) = r_cum(1 + horizon - 1:end - 1, horizon - 1) + return_debt_tilde(1 + horizon:end, 2);
    s_cum(1 + horizon:end, horizon) = s_cum(1 + horizon - 1:end - 1, horizon - 1) + surplus_tilde(1 + horizon:end, 2);
end

load('CleanData/simul.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% long sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[coeff_RS_sim, coeff_R_sim, coeff_S_sim, coeff_V_sim, mdl_sest_sim, mdl_rest_sim, mdl_vest_sim] = ...
    estimate_vsr(v_long, r_long, s_long, 10);

coeff_V_sim_long = coeff_V_sim(1, :);
coeff_S_sim_long = coeff_S_sim(1, :);
coeff_R_sim_long = coeff_R_sim(1, :);
coeff_RS_sim_long = coeff_RS_sim(1, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nboot_total = 1e4;
len_data = T;

coeff_V_sim_total = [];
coeff_S_sim_total = [];
coeff_R_sim_total = [];
coeff_RS_sim_total = [];

parfor nboot = 1:nboot_total
    rng(nboot);
    E_rnd = datasample(E_sim, len_data, 'Replace', true);
    E_v = E_rnd(:, 1);
    E_r = E_rnd(:, 2);

    v_sim = zeros(len_data, 1);
    r_sim = zeros(len_data, 1);
    s_sim = zeros(len_data, 1);

    idx = randi([1, len_long - 1], 1, 1);
    v_sim(1) = v_long(idx);
    r_sim(1) = r_long(idx);
    s_sim(1) = s_long(idx);

    v_sim(2) = v_long(idx + 1);
    r_sim(2) = r_long(idx + 1);
    s_sim(2) = s_long(idx + 1);

    for t = 1:(len_data - 2)
        v_sim(t + 2) = [v_sim(t + 2 - 1) v_sim(t + 2 - 2)] * [phi1v, phi2v]' + E_v(t);
        r_sim(t + 2) = [v_sim(t + 2 - 1) v_sim(t + 2 - 2)] * [phi1r, phi2r]' + E_r(t);
        s_sim(t + 2) = v_sim(t + 1) - v_sim(t + 2) + r_sim(t + 2);
    end

    [coeff_RS_sim, coeff_R_sim, coeff_S_sim, coeff_V_sim, mdl_sest_sim, mdl_rest_sim, mdl_vest_sim] = ...
        estimate_vsr(v_sim, r_sim, s_sim, 10);

    coeff_V_sim_total(nboot, :) = coeff_V_sim(1, :);
    coeff_S_sim_total(nboot, :) = coeff_S_sim(1, :);
    coeff_R_sim_total(nboot, :) = coeff_R_sim(1, :);
    coeff_RS_sim_total(nboot, :) = coeff_RS_sim(1, :);
end

for k = 1:10
    mean_coeff_V(k) = nanmean(coeff_V_sim_total(:, k));
    std_coeff_V(k) = nanstd(coeff_V_sim_total(:, k));

    mean_coeff_S(k) = nanmean(coeff_S_sim_total(:, k));
    std_coeff_S(k) = nanstd(coeff_S_sim_total(:, k));

    mean_coeff_R(k) = nanmean(coeff_R_sim_total(:, k));
    std_coeff_R(k) = nanstd(coeff_R_sim_total(:, k));

    mean_coeff_RS(k) = nanmean(coeff_RS_sim_total(:, k));
    std_coeff_RS(k) = nanstd(coeff_RS_sim_total(:, k));
end

coeff_R_uni(1, :) = mean_coeff_R;
coeff_V_uni(1, :) = mean_coeff_V;
coeff_S_uni(1, :) = mean_coeff_S;
coeff_RS_uni(1, :) = mean_coeff_RS;

coeff_R_uni(2, :) = std_coeff_R;
coeff_V_uni(2, :) = std_coeff_V;
coeff_S_uni(2, :) = std_coeff_S;
coeff_RS_uni(2, :) = std_coeff_RS;

Table_tex = round([coeff_R_uni(:, 1:10); coeff_S_uni(:, 1:10); coeff_V_uni(:, 1:10)], 2, 'significant')
rowLabels = [];
columnLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
newpath_tables = strcat(pwd, '\tables\');
matrix2latex(Table_tex, strcat(newpath_tables, 'tables_data_nobreak_UR', '.txt'));

horizon = 10;

coeff_S_total = coeff_S_uni(1, :);
coeff_V_total = coeff_V_uni(1, :);
coeff_R_total = coeff_R_uni(1, :);
coeff_RS_total = coeff_RS_uni(1, :);

std_coeff_S = coeff_S_uni(2, :);
std_coeff_V = coeff_V_uni(2, :);
std_coeff_R = coeff_R_uni(2, :);
std_coeff_RS = coeff_RS_uni(2, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

horizon = 10;

color = 'b';

f = figure;
subplot(2, 2, 1);
plot_fun(coeff_S(1, :) + coeff_S_sim_long - coeff_S_total, std_coeff_S, horizon, color);
plot(linspace(1, 10, 10), coeff_S_total, 'LineWidth', 3, 'Color', 'r');
title('Future $S$', 'interpreter', 'latex');

subplot(2, 2, 2);
plot_fun(coeff_R(1, :) + coeff_R_sim_long - coeff_R_total, std_coeff_R, horizon, color);
plot(linspace(1, 10, 10), coeff_R_total, 'LineWidth', 3, 'Color', 'r');
title('Future $-R$', 'interpreter', 'latex');

subplot(2, 2, 3);
plot_fun(coeff_V(1, :) + coeff_V_sim_long - coeff_V_total, std_coeff_V, horizon, color);
plot(linspace(1, 10, 10), coeff_V_total, 'LineWidth', 3, 'Color', 'r');
title('Future $V$', 'interpreter', 'latex');

subplot(2, 2, 4);
plot_fun(coeff_RS(1, :) + coeff_RS_sim_long - coeff_RS_total, std_coeff_RS, horizon, color);
plot(linspace(1, 10, 10), coeff_RS_total, 'LineWidth', 3, 'Color', 'r');
title('Future $S-R$', 'interpreter', 'latex');

filename = strcat('figs/boot_CSdecomp_AR2_direct_nobreak_simulation', '.pdf')
f.PaperSize = [6 6];
print(filename, '-dpdf', '-fillpage');
