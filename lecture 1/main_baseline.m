nboot_total = 10000;

if (bias_adjust == 1)
    label = 'nobreak_biasadj';
    color = 'b';
else
    label = 'nobreak';
    color = 'r';
end

if (break_debt == 1)
    label = label(3:end);
    load('CleanData/US_data_break.mat');
else
    load('CleanData/US_data.mat');
end

horizon_total = 10;

% result from data
[coeff_RS, coeff_R, coeff_S, coeff_V, mdl_sest, mdl_rest, mdl_vest] = ...
    estimate_vsr(v_tilde(:, 2), return_debt_tilde(:, 2), surplus_tilde(:, 2), horizon_total);

coeff_R_uni = coeff_R;
coeff_S_uni = coeff_S;
coeff_V_uni = coeff_V;
coeff_RS_uni = coeff_RS;
results_mat = [coeff_R(:, 1:horizon_total); coeff_S(:, 1:horizon_total); coeff_V(:, 1:horizon_total)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bias adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (bias_adjust == 1)
    T0 = size(v_tilde(:, 2), 1) - 1;
    phi_1 = coeff_V(1, 1);
    rho(1) = phi_1;
    E_v = mdl_vest{1}.Residuals.Raw;
    E_r = mdl_rest{1}.Residuals.Raw;
    E_s = mdl_sest{1}.Residuals.Raw;
    Res = [E_v E_r E_s];
    Cov = cov(Res);

    for j = 1:50
        bias(1, j) =- (1 / T0) * (j * (1 + rho(1)) + 2 * rho(1) * ((1 - rho(1) ^ j) / (1 - rho(1)))) * (Cov(1, 2) / Cov(1, 1));
        bias(2, j) =- (1 / T0) * (j * (1 + rho(1)) + 2 * rho(1) * ((1 - rho(1) ^ j) / (1 - rho(1)))) * (Cov(1, 3) / Cov(1, 1));
    end

    coeff_R_uni(1, :) = coeff_R(1, :) - bias(1, 1:horizon_total);
    coeff_S_uni(1, :) = coeff_S(1, :) - bias(2, 1:horizon_total);
    coeff_V_uni(1, :) = coeff_V(1, :) + bias(1, 1:horizon_total) + bias(2, 1:horizon_total);
    coeff_RS_uni(1, :) = coeff_R_uni(1, :) + coeff_S_uni(1, :);

    % report corr
    sum(abs(E_v + E_r + E_s))
    corr([E_v -E_r E_s])
    Mdl = arima(1, 0, 0);
    EstMdl = estimate(Mdl, v_tilde(:, 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len_data = T;

load('CleanData/simul.mat');

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
        estimate_vsr(v_sim, r_sim, s_sim, horizon_total);

    coeff_V_sim_total(nboot, :) = coeff_V_sim(1, :);
    coeff_S_sim_total(nboot, :) = coeff_S_sim(1, :);
    coeff_R_sim_total(nboot, :) = coeff_R_sim(1, :);
    coeff_RS_sim_total(nboot, :) = coeff_RS_sim(1, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% report results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:horizon_total
    mean_coeff_V(k) = nanmean(coeff_V_sim_total(:, k));
    std_coeff_V(k) = nanstd(coeff_V_sim_total(:, k));

    mean_coeff_S(k) = nanmean(coeff_S_sim_total(:, k));
    std_coeff_S(k) = nanstd(coeff_S_sim_total(:, k));

    mean_coeff_R(k) = nanmean(coeff_R_sim_total(:, k));
    std_coeff_R(k) = nanstd(coeff_R_sim_total(:, k));

    mean_coeff_RS(k) = nanmean(coeff_RS_sim_total(:, k));
    std_coeff_RS(k) = nanstd(coeff_RS_sim_total(:, k));
end

coeff_R_uni(2, :) = std_coeff_R;
coeff_V_uni(2, :) = std_coeff_V;
coeff_S_uni(2, :) = std_coeff_S;
coeff_RS_uni(2, :) = std_coeff_RS;

Table_tex = round([results_mat(1, 1:horizon_total); coeff_R_uni(2:3, 1:horizon_total); coeff_R_uni(1, 1:horizon_total); ...
                       results_mat(4, 1:horizon_total); coeff_S_uni(2:3, 1:horizon_total); coeff_S_uni(1, 1:horizon_total); ...
                       results_mat(7, 1:horizon_total); coeff_V_uni(2:3, 1:horizon_total); coeff_V_uni(1, 1:horizon_total)], 2);
columnLabels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
rowLabels = {'$-b_T^r$', '$s.e.$', '$R^2$', '$unbiased$', ...
                 '$b_T^s$', '$s.e.$', '$R^2$', '$unbiased$', ...
                 '$\phi$', '$s.e.$', '$R^2$', '$unbiased$'};
matrix2latex(Table_tex, strcat('tables/tables_data_', label, '.txt'), 'rowLabels', rowLabels);

coeff_RS_total = coeff_RS_uni(1, :);
coeff_S_total = coeff_S_uni(1, :);
coeff_V_total = coeff_V_uni(1, :);
coeff_R_total = coeff_R_uni(1, :);

std_coeff_S = coeff_S_uni(2, :);
std_coeff_V = coeff_V_uni(2, :);
std_coeff_R = coeff_R_uni(2, :);
std_coeff_RS = coeff_RS_uni(2, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
subplot(2, 2, 1);
plot_fun(coeff_S_total, std_coeff_S, horizon_total, color);
title('Future $S$', 'interpreter', 'latex');

subplot(2, 2, 2);
plot_fun(coeff_R_total, std_coeff_R, horizon_total, color);
title('Future $-R$', 'interpreter', 'latex');

subplot(2, 2, 3);
plot_fun(coeff_V_total, std_coeff_V, horizon_total, color);
title('Future $V$', 'interpreter', 'latex');

subplot(2, 2, 4);
plot_fun(coeff_RS_total, std_coeff_RS, horizon_total, color);
title('Future $S-R$', 'interpreter', 'latex');

newpath_figs = strcat(pwd, '\figs\');

filename = strcat('figs/boot_CSdecomp_AR2_direct_', label, '.pdf')
f.PaperSize = [6 6];
print(filename, '-dpdf', '-fillpage');

sprintf("%.2f ", [coeff_S_total(5), coeff_R_total(5), coeff_V_total(5)] * 100)
sprintf("%.2f ", [coeff_S_total(end), coeff_R_total(end), coeff_V_total(end)] * 100)
