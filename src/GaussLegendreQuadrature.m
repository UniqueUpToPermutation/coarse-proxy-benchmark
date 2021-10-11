function fval = GaussLegendreQuadrature(fun, Q,lvar,weight)
%Q points Gauss-Legendre quadrature: x = 0..1
    fval = 0;
    for i = 1:Q
        fval = fval + weight(i) .* fun(lvar(i));
    end
end

