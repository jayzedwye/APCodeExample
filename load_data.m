clear
close all

newpath_tables = strcat(pwd, '\tables\');
newpath_figs = strcat(pwd, '\figs\');
[data_CRSP, txt_CRSP] = xlsread('CleanData/MasterData_JF.xlsx', 'CRSP');
[data_VAR, txt_VAR] = xlsread('CleanData/MasterData_JF.xlsx', 'Sheet1');


data_VAR_sample = data_VAR(2:end, :); % 1930 is 2nd row
data_t0 = data_VAR(1, :); % 1929
T = size(data_VAR_sample, 1);
dates = data_VAR_sample(:, 1);

x_tilde = zeros(T, 2);
x_tilde(:, 1) = dates;
x_tilde(:, 2) = data_VAR_sample(:, 5);

pi_tilde = zeros(T, 2);
pi_tilde(:, 1) = dates;
pi_tilde(:, 2) = data_VAR_sample(:, 6);

return_debt = zeros(T, 2);
return_debt(:, 1) = dates;
return_debt(:, 2) = log(1 + data_CRSP(2:end, 2)); % CRSP data
% return_debt0(:, 2) = log((data_CRSP(2:end, 2) + data_CRSP(2:end, 3) - data_CRSP(2:end, 5) + data_CRSP(2:end, 6)) ./ data_CRSP(1:end - 1, 6));
% corr(return_debt0(:, 2), return_debt(:, 2))

return_debt_tilde = zeros(T, 2);
return_debt_tilde(:, 1) = dates;
return_debt_tilde(:, 2) = return_debt(:, 2) - x_tilde(:, 2) - pi_tilde(:, 2);

v_tilde = zeros(T, 2);
v_tilde(:, 1) = dates;
v_tilde(:, 2) = log(data_VAR_sample(:, 7));
v_tilde0 = log(data_t0(7));

surplus_tilde = zeros(T, 2);
surplus_tilde(:, 1) = dates;
surplus_tilde(1, 2) = v_tilde0 - v_tilde(1, 2) + return_debt_tilde(1, 2);
surplus_tilde(2:end, 2) = v_tilde(1:end - 1, 2) - v_tilde(2:end, 2) + return_debt_tilde(2:end, 2);

% check corr with surplus/gdp
surplus_gdp = zeros(T, 2);
surplus_gdp(:, 1) = dates;
surplus_gdp(:, 2) = data_VAR_sample(:, 4);

% corr(surplus_tilde(:, 2), surplus_gdp(:, 2))
% corr(surplus_tilde(18:end, 2), surplus_gdp(18:end, 2))

% subsampling
t1947 = find(v_tilde(:, 1) == 1947);
surplus_tilde = surplus_tilde(t1947:end,:);
v_tilde = v_tilde(t1947:end,:);
return_debt_tilde = return_debt_tilde(t1947:end,:);
x_tilde = x_tilde(t1947:end,:);
pi_tilde = pi_tilde(t1947:end,:);
surplus_gdp = surplus_gdp(t1947:end,:);
dates = dates(t1947:end,:);
T = size(dates, 1);

% end result: 1947-2022
save("CleanData/US_data.mat");

