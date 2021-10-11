function [oracle, sample_space] = GaussianRTE()
    
    % Random parameters are just center offsets
    np2 = 20;    
    pgrid = (0:np2) / np2;
    [xs,ys] = ndgrid(pgrid, pgrid);
    sample_space = [xs(:)';    ys(:)'];

    % mu_t(x) and mu_s(x), transmission and scattering
    sig = 0.2;
    mu_s_func = @(x1, x2, param) 1 + 1.0 * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / sig^2);
    mu_t_func = @(x1, x2, param) 1 + 1.0 * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / sig^2);

    % Force function
    force_func = @(x1, x2) exp(-((x1-0.5).^2+(x2-0.5).^2)*256);
    
    % Create the oracle
    oracle = OracleRTE(force_func, mu_t_func, mu_s_func);
end

