function A = Bfunnonsym(i,j, param, b_add_identity)
    xs = param.xs;
    diagval = param.diagval;
    N = param.N;
    mus = param.mus;
    
    % Generates the operator K
    A = Kfun(xs(:,i),xs(:,j),param)  / N;
    
    % Generates A from K? 

    assert(size(mus, 2) == 1);
    A = mus(i) .* A; % Multiply rows of A by mu_s to form mu_s K
    
    [I,J] = ndgrid(i,j);
    A(I==J) = diagval(I(I==J));
    A(I==J) = A(I==J) + b_add_identity; % Add identity to form I + mu_s K
end
