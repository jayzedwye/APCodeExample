function err = mmtAR2vsData(x, T, v_tilde, coeff_V, horizon_total)

    phi1 = x(1);
    phi2 = x(2);

    if (max(abs(eig([phi1, phi2; 1, 0]))) >= 0.9999)
        err = 1e10;
        return;
    end

    % find residuals
    demean = @(x) (x - mean(x));

    Y_Beta = v_tilde(2 + 1:end, 2);
    X2_Beta = [v_tilde(2:end - 1, 2) v_tilde(1:end - 1 - 1, 2)];
    E_v = demean(Y_Beta) - demean(X2_Beta) * [phi1, phi2]';

    % Y_Beta = r_cum(2 + 1:end, 1);
    % E_r = demean(Y_Beta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bootstrap
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
            v_sim(t + 2) = [v_sim(t + 2 - 1) v_sim(t + 2 - 2)] * [phi1, phi2]' + E_v(t);
            %         r_sim(t + 2) = E_r(t);
            %         s_sim(t + 2) = v_sim(t + 1) - v_sim(t + 2) + r_sim(t + 2);
        end

        for horizon = 1:horizon_total
            Y_Beta = v_sim(1 + horizon:end);
            X_Beta = v_sim(1:end - horizon);
            mdl_vest = LinearModel.fit([X_Beta(:, :)], Y_Beta);
            coeff_V_sim(1, horizon) = mdl_vest.Coefficients.Estimate(2);
        end

        coeff_V_sim_total(nboot, :) = coeff_V_sim(1, :);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % report results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mean_coeff_V = mean(coeff_V_sim_total);

    err = sum((coeff_V(1, :) - mean_coeff_V) .^ 2);
