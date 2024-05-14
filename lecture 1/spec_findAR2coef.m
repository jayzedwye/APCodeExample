clear; close all;
load('CleanData/US_data.mat');
horizon_total = 10;

% result from data
[coeff_RS, coeff_R, coeff_S, coeff_V, mdl_sest, mdl_rest, mdl_vest] = ...
    estimate_vsr(v_tilde(:, 2), return_debt_tilde(:, 2), surplus_tilde(:, 2), horizon_total);

coeff_R_uni = coeff_R;
coeff_S_uni = coeff_S;
coeff_V_uni = coeff_V;
coeff_RS_uni = coeff_RS;
results_mat = [coeff_R(:, 1:horizon_total); coeff_S(:, 1:horizon_total); coeff_V(:, 1:horizon_total)];

% construct cumulative series
r_cum = zeros(T, horizon_total);
s_cum = zeros(T, horizon_total);

r_cum(1:end, 1) = return_debt_tilde(:, 2);
s_cum(1:end, 1) = surplus_tilde(:, 2);
for horizon = 2:horizon_total
    r_cum(1 + horizon:end, horizon) = r_cum(1 + horizon - 1:end - 1, horizon - 1) + return_debt_tilde(1 + horizon:end, 2);
    s_cum(1 + horizon:end, horizon) = s_cum(1 + horizon - 1:end - 1, horizon - 1) + surplus_tilde(1 + horizon:end, 2);
end

% % guess the AR2 coefficients
% phi1s = .6:.05:1.4;
% phi2s = -.5:.05:1;
% 
% for i = 1:length(phi1s)
%     for j = 1:length(phi2s)
%         phi1 = phi1s(i);
%         phi2 = phi2s(j);
% 
%         if (max(abs(eig([phi1, phi2; 1, 0]))) >= 1)
%             res(i,j) = NaN;
%             continue;
%         end
% 
%         rho1 = phi1/(1-phi2);
%         if (rho1 <= .85)
%             res(i,j) = NaN;
%             continue;
%         end
% 
%         res(i,j) = mmtAR2vsData([phi1, phi2], T, v_tilde, coeff_V, horizon_total);
%     end
% 
%     save resZJ.mat res;
% end

phi1v = 1.42082559207046;
phi2v = -0.420883510141803;

% phi1v = 1.4542;
% phi2v = -0.48434;

% mmtAR2vsData([phi1v,phi2v], T, v_tilde, coeff_V)

param_init = [phi1v, phi2v];
options = optimset('DiffMinChange', 1e-2, ...
    'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1e2, 'MaxFunEval', 1e8, 'Display', 'iter');
[param, fval, exitflag, output] = fminsearch('mmtAR2vsData', param_init, options,...
    T, v_tilde, coeff_V, horizon_total);


demean = @(x) (x - mean(x));

Y_Beta = v_tilde(2 + 1:end, 2);
X2_Beta = [v_tilde(2:end - 1, 2) v_tilde(1:end - 1 - 1, 2)];
E_v = demean(Y_Beta) - demean(X2_Beta) * [phi1v, phi2v]';
autocorr(E_v)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nboot_total = 1e4;
len_data = T;

E_sim = [E_v];

coeff_V_sim_total = [];

parfor nboot = 1:nboot_total
    rng(nboot);
    coeff_V_sim = zeros(1, horizon_total);
    E_rnd = datasample(E_sim, len_data, 'Replace', true);
    E_v = E_rnd(:, 1);
    %     E_r = E_rnd(:, 2);

    v_sim = zeros(len_data, 1);
    %     r_sim = zeros(len_data, 1);
    %     s_sim = zeros(len_data, 1);

    for t = 1:len_data
        v_sim(t + 2) = [v_sim(t + 2 - 1) v_sim(t + 2 - 2)] * [phi1v, phi2v]' + E_v(t);
        %         r_sim(t + 2) = E_r(t);
        %         s_sim(t + 2) = v_sim(t + 1) - v_sim(t + 2) + r_sim(t + 2);
    end

    for horizon=1:horizon_total
        Y_Beta                 = v_sim(1+horizon:end);
        X_Beta                 = v_sim(1:end-horizon);
        mdl_vest               = LinearModel.fit([ X_Beta(:,:) ],Y_Beta);
        coeff_V_sim(1,horizon) = mdl_vest.Coefficients.Estimate(2);
    end

    coeff_V_sim_total(nboot, :) = coeff_V_sim(1, :);
end

mean_coeff_V = mean(coeff_V_sim_total);

err = sum((coeff_V(1, :) - mean_coeff_V).^2);

f= figure;
plot(1:horizon_total, coeff_V(1, :), 1:horizon_total, mean_coeff_V)
legend("data","best AR(2)")

f.PaperSize = [6 6];
print('figs/simul-bestAR2','-dpdf','-fillpage');
