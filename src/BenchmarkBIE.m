% Generates an oracle and a sample space for our boundary integral equation
% (i.e. Fredholm equation). 
% Arguments: 
%   benchmark (string) - Determines which benchmark to create, options are:
%       'test' : Spawns a simple benchmark for testing if the software works.
%       'final' : Spawns the benchmark used in the paper.
function [oracle, sample_space] = BenchmarkBIE(benchmark)
    if nargin == 0
        benchmark = 'test';
    end
    
    if strcmp(benchmark, 'test')
        [oracle, sample_space] = SingularBIE_simple();
    elseif strcmp(benchmark, 'final')
        [oracle, sample_space] = SingularBIE_final();
    else
        error('benchmark is not recognized');
    end
end

% Uses
function [oracle, sample_space] = SingularBIE_simple()
    sample_space_size = 1024 * 16;
    nK = 8;
    dev = 0.4;
    sample_space = 1 + dev*(2*rand(nK, sample_space_size)-1);
    n_coarse = 32;
    n_fine = 512;
    oracle = OracleBIE(n_coarse, n_fine, @singularF);
end 

% Uses significantly more samples than SingularBIE_simple
function [oracle, sample_space] = SingularBIE_final()
    sample_space_size = 4 * 1024 * 8;
    nK = 8;
    dev = 0.4;
    sample_space = 1 + dev*(2*rand(nK, sample_space_size)-1);
    n_coarse = 32;
    n_fine = 1024;
    oracle = OracleBIE(n_coarse, n_fine, @singularF);
end

function [result] = singularF(x, y)
    amp = 1;
    x0 = 0.6;
    y0 = 0;
    result = amp ./ sqrt((x - x0).^2 + (y - y0).^2);
end