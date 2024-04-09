clear all; close all; clc
loaddata;
% number of params
nparam = 20;

load('results.mat')

idx = 1;
cy = zeros(length(date), 1);
param = minparamval(:, idx);

[N, T, Psi, Sig, I_pi, I_gdp, I_y1, I_yspr, I_cy, inflpos, gdppos, y1pos, ...
     ysprpos, pi0, x0, ynom1q0, ...
     yspr0, cy0, X2, yielddata, yieldmaturity, eps2] = ...
    setup_model(yielddata, yieldmaturity, cy, ...
    [], [], infl, x, ...
    [], []);

mmt(param, N, T, Psi, Sig, I_pi, I_gdp, I_y1, I_yspr, I_cy, ...
    inflpos, gdppos, y1pos, ysprpos, pi0, x0, ynom1q0, yspr0, cy0, X2, yielddata, yieldmaturity, eps2)

%% No arbitrage restrictions on bond yields
x = param;
striphorizon = 1e2;

L0 = zeros(N, 1);
L1 = zeros(N, N);

L0(1:4) = x(1:4)';
tmp = zeros(4:4);
tmp(:) = x((4 + 1):(4 + 4 ^ 2));

L1(1:4, 1:4) = tmp ./ std(X2(:, 1:4));

Api = zeros(striphorizon, 1);
Bpi = zeros(N, striphorizon);

Api(1) = -ynom1q0 + cy0;
Bpi(:, 1) = -I_y1' + I_cy';

for j = 1:striphorizon
    Api(j + 1) =- ynom1q0 + Api(j) + .5 * Bpi(:, j)' * (Sig * Sig') * Bpi(:, j) - Bpi(:, j)' * Sig * L0;
    Bpi(:, j + 1) = (Bpi(:, j)' * Psi - I_y1' - Bpi(:, j)' * Sig * L1)';
end

predicted_yield = kron(ones(T, 1), -Api(yieldmaturity)' ./ yieldmaturity) - ((Bpi(:, yieldmaturity)' ./ kron(yieldmaturity', ones(1, N))) * X2')';

% 1Q
subplot(2, 4, 1)
plot(date, yielddata(:, 1), date, predicted_yield(:, 1), 'LineWidth', 2)
% 1Y
subplot(2, 4, 2)
plot(date, yielddata(:, 2), date, predicted_yield(:, 2), 'LineWidth', 2)
% 2Y
subplot(2, 4, 3)
plot(date, yielddata(:, 3), date, predicted_yield(:, 3), 'LineWidth', 2)
% 3Y
subplot(2, 4, 4)
plot(date, yielddata(:, 4), date, predicted_yield(:, 4), 'LineWidth', 2)
% 4Y
subplot(2, 4, 5)
plot(date, yielddata(:, 5), date, predicted_yield(:, 5), 'LineWidth', 2)
% 5Y
subplot(2, 4, 6)
plot(date, yielddata(:, 6), date, predicted_yield(:, 6), 'LineWidth', 2)
% 7Y
subplot(2, 4, 7)
plot(date, yielddata(:, 7), date, predicted_yield(:, 7), 'LineWidth', 2)
% 10Y
subplot(2, 4, 8)
plot(date, yielddata(:, 8), date, predicted_yield(:, 8), 'LineWidth', 2)

% optimizer
tic

optionssimplex = optimset('TolX', 1e-16, 'TolFun', 1e-16, 'MaxIter', 1e5, ...
    'MaxFunEval', 1e8, 'display', 'iter');

for ii = 1:10
    [param, fval] = ...
        fminsearch('mmt', param, optionssimplex, N, T, Psi, Sig, I_pi, ...
        I_gdp, I_y1, I_yspr, I_cy, inflpos, gdppos, y1pos, ...
        ysprpos, pi0, x0, ynom1q0, ...
        yspr0, cy0, X2, yielddata, yieldmaturity, eps2);
end

toc

mmt(param, N, T, Psi, Sig, I_pi, I_gdp, I_y1, I_yspr, I_cy, ...
    inflpos, gdppos, y1pos, ysprpos, pi0, x0, ynom1q0, yspr0, cy0, X2, yielddata, yieldmaturity, eps2)
