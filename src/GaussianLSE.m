% A simple hello world for LSE 
function [oracle, sample_space] = GaussianLSE()
    K_coarse = 4;
    K_fine = 2 * K_coarse;
    npw = 6;
    
    oracle = OracleLSE(K_coarse, K_fine, npw, @potential_gen, @force_gen);
    
    % Trivial sample space, for now
    sample_space = 1;
end

function [V] = potential_gen(N, ~)
    h = 1/(N + 1);   
    [x1, x2, x3] = ndgrid(h:h:1-h);
    sig = 0.2;

    V = exp(-((x1-0.5).^2 -(x2-0.5).^2 - (x3-0.5).^2) / sig^2);
end

function [f] = force_gen(N, ~)
    f = zeros(N, N, N);
    h = 1/(N + 1);
    pos = int32(N / 2);
    f(pos, pos, pos) = 1 / h^3;
end