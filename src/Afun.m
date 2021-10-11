function A = Afun(i,j,param)
    xs = param.xs;
    diagval = param.diagval;
    N = param.N;
    
    % Generates the operator K
    A = Kfun(xs(:,i),xs(:,j),param)  / N;
    
    % Generates A from K? 
    [I,J] = ndgrid(i,j);
    A(I==J) = diagval(I(I==J));
end
