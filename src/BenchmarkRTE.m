% Generates an oracle and a sample space for our radiative transport equation
%
% Arguments: 
%   benchmark (string) - Determines which benchmark to create, options are:
%       'test' : Spawns a simple benchmark for testing if the software works.
%       'final' : Spawns the benchmark used in the paper.
function [oracle, sample_space] = BenchmarkRTE(benchmark)
    if nargin == 0
        benchmark = 'test';
    end
    
    if strcmp(benchmark, 'test')
        [oracle, sample_space] = RTE_gaussian_simple();
    elseif strcmp(benchmark, 'final')
        [oracle, sample_space] = RTE_gaussian_final();
    else 
        error('benchmark is not recognized');
    end
end

function [oracle, sample_space] = RTE_gaussian_simple()
    
    % Random parameters are just center offsets
    np2 = 20;    
    pgrid = (0:np2) / np2;
    [xs,ys] = ndgrid(pgrid, pgrid);
    sample_space_temp = [xs(:)';    ys(:)'];

    % mu_t(x) and mu_s(x), transmission and scattering
    sigs = [0.2, 0.4, 0.6];
    rhos = [10, 7.5, 5];
    mu_s_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    mu_t_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    
    % Force function
    force_func = @(x1, x2) exp(-((x1-0.5).^2+(x2-0.5).^2)*256);
    
    % Create the oracle
    oracle = OracleRTE(force_func, mu_t_func, mu_s_func);
    
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

function [oracle, sample_space] = RTE_gaussian_final()
    
    % Random parameters are just center offsets
    np2 = 20;    
    pgrid = (0:np2) / np2;
    [xs,ys] = ndgrid(pgrid, pgrid);
    sample_space_temp = [xs(:)';    ys(:)'];

    % mu_t(x) and mu_s(x), transmission and scattering
    sigs = 0.2:0.1:0.6;
    rhos = 2:2:10;
    mu_s_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    mu_t_func = @(x1, x2, param) 1 + param(4) * ...
        exp(-((x1 - param(1)).^2 + (x2 - param(2)).^2) / param(3)^2);
    
    % Force function
    force_func = @(x1, x2) exp(-((x1-0.5).^2+(x2-0.5).^2)*256);
    
    % Create the oracle
    oracle = OracleRTE(force_func, mu_t_func, mu_s_func);
    
    % Create the sample space
    num_sigs = length(sigs);
    num_rhos = length(rhos);
    num_temp_sample_space = size(sample_space_temp, 2);
    sample_space = zeros(4, num_temp_sample_space * num_sigs * num_rhos);
    for i = 1:num_sigs
        for j = 1:num_rhos
            start_col = num_temp_sample_space * (num_rhos * (i - 1) + (j - 1)) + 1;
            end_col = num_temp_sample_space * (num_rhos * (i - 1) + j);
            sample_space(1:2, start_col:end_col) = sample_space_temp;
            sample_space(3, start_col:end_col) = sigs(i);
            sample_space(4, start_col:end_col) = rhos(j);
        end
    end
end