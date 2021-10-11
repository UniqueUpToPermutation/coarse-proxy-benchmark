% Generates an oracle and a sample space for our boundary integral equation
% (i.e. Fredholm equation). 
function [oracle, sample_space] = SingularBIE()
    rng(2142);
    sample_space_size = 1024 * 16;
    nK = 8;
    dev = 0.4;
    sample_space = 1 + dev*(2*rand(nK, sample_space_size)-1);
    n_coarse = 32;
    n_fine = 512;
    oracle = OracleBIE(n_coarse, n_fine, @singularF);
end

function [result] = singularF(x, y)
    amp = 1;
    x0 = 0.6;
    y0 = 0;
    result = amp ./ sqrt((x - x0).^2 + (y - y0).^2);
end