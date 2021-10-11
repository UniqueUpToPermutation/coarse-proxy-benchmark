function [oracle, sample_space] = GuassianRTEHarderNonSym()
    
    % Random parameters are just center offsets
    np2 = 20;    
    pgrid = (0:np2) / np2;
    [xs,ys] = ndgrid(pgrid, pgrid);
    sample_space_temp = [xs(:)';    ys(:)'];

    % mu_t(x) and mu_s(x), transmission and scattering
    %sigs = [0.2, 0.4, 0.6];
    sigs = [0.2, 0.4, 0.6];
    %rhos = [10, 5, 3.33];
    rhos = [10, 7.5, 5];
    mu_s_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    mu_t_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    
    % Force function
    force_func = @(x1, x2) exp(-((x1-0.5).^2+(x2-0.5).^2)*256);
    
    % Create the oracle
    oracle = OracleRTENonSym(force_func, mu_t_func, mu_s_func);
    
    % Create the sample space
    num_sigs = length(sigs);
    num_temp_sample_space = size(sample_space_temp, 2);
    sample_space = zeros(4, num_temp_sample_space * num_sigs);
    for i = 1:num_sigs
        start_col = num_temp_sample_space * (i - 1) + 1;
        end_col = num_temp_sample_space * i;
        sample_space(1:2, start_col:end_col) = sample_space_temp;
        sample_space(3, start_col:end_col) = sigs(i);
        sample_space(4, start_col:end_col) = rhos(i);
    end
end