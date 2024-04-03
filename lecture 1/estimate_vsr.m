function [coeff_RS, coeff_R, coeff_S, coeff_V, mdl_sest, mdl_rest, mdl_vest] = ...
        estimate_vsr(v, r_tilde, s, horizon_total)


    T = size(v, 1);
    r_cum = zeros(T, horizon_total);
    s_cum = zeros(T, horizon_total);

    r_cum(1:end, 1) = r_tilde;
    s_cum(1:end, 1) = s;

    for horizon = 2:horizon_total
        r_cum(1 + horizon:end, horizon) = r_cum(1 + horizon - 1:end - 1, horizon - 1) + r_tilde(1 + horizon:end);
        s_cum(1 + horizon:end, horizon) = s_cum(1 + horizon - 1:end - 1, horizon - 1) + s(1 + horizon:end);
    end

    for horizon = 1:horizon_total
        Y_Beta = s_cum(1 + horizon:end, horizon) - r_cum(1 + horizon:end, horizon);
        X_Beta = v(1:end - horizon);
        mdl{horizon} = LinearModel.fit([X_Beta(:, :)], Y_Beta);
        tmp_se = mdl{horizon}.Coefficients.SE(2);

        k = horizon;
        coeff_RS(1, k) = mdl{horizon}.Coefficients.Estimate(2);
        coeff_RS(2, k) = tmp_se;
        coeff_RS(3, k) = mdl{horizon}.Rsquared.Ordinary;
    end

    for horizon = 1:horizon_total
        Y_Beta = s_cum(1 + horizon:end, horizon);
        X_Beta = v(1:end - horizon);
        mdl_sest{horizon} = LinearModel.fit([X_Beta(:, :)], Y_Beta);
        tmp_se = mdl_sest{horizon}.Coefficients.SE(2);

        k = horizon;
        coeff_S(1, k) = mdl_sest{horizon}.Coefficients.Estimate(2);
        coeff_S(2, k) = tmp_se;
        coeff_S(3, k) = mdl_sest{horizon}.Rsquared.Ordinary;
    end

    for horizon = 1:horizon_total
        Y_Beta = -r_cum(1 + horizon:end, horizon);
        X_Beta = v(1:end - horizon);
        mdl_rest{horizon} = LinearModel.fit([X_Beta(:, :)], Y_Beta);
        tmp_se = mdl_rest{horizon}.Coefficients.SE(2);

        k = horizon;
        coeff_R(1, k) = mdl_rest{horizon}.Coefficients.Estimate(2);
        coeff_R(2, k) = tmp_se;
        coeff_R(3, k) = mdl_rest{horizon}.Rsquared.Ordinary;
    end

    for horizon = 1:horizon_total
        Y_Beta = v(1 + horizon:end);
        X_Beta = v(1:end - horizon);
        mdl_vest{horizon} = LinearModel.fit([X_Beta(:, :)], Y_Beta);
        tmp_se = mdl_vest{horizon}.Coefficients.SE(2);

        k = horizon;
        coeff_V(1, k) = mdl_vest{horizon}.Coefficients.Estimate(2);
        coeff_V(2, k) = tmp_se;
        coeff_V(3, k) = mdl_vest{horizon}.Rsquared.Ordinary;
    end

end
