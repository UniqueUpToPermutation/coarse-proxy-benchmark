% Generates an oracle and a sample space for our boundary integral equation
% (i.e. Fredholm equation). 
function [oracle, sample_space] = MoreSamplesSingularBIE()
    rng(2142);
    sample_space_size = 4 * 1024 * 8;
    nK = 8;
    dev = 0.4;
    sample_space = 1 + dev*(2*rand(nK, sample_space_size)-1);
    n_coarse = 32;
    n_fine = 1024;
    oracle = OracleBIE(n_coarse, n_fine, @gaussianF);
end

function [result] = gaussianF(x, y)
    % x0 = 1;
    % y0 = 1;
    % sigma = 0.5;
    amp = 1;
    
    % result = amp * exp(-((x - x0)/sigma) .^ 2 - ((y - y0)/sigma) .^ 2);
    % result = amp * sin(50 * 2 * pi * x);
    
    x0 = 0.6;
    y0 = 0;
    result = amp ./ sqrt((x - x0).^2 + (y - y0).^2);
end