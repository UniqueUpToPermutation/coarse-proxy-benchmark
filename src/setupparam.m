function [N,xs,param,mus] = setupparam(Q,nC,proxy,musfun,mutfun)
    
    [lvar, weight] = GaussLegendre_2(Q);%compute Gaussian quadrature weights
    lvar = (1+lvar)/2;
    weight = weight/2;
    
    [x1,x2] = ndgrid( (1/2:nC)/nC );        xs = [x1(:) x2(:)]';
    N = size(xs,2); %n^2
    
    mut = mutfun(xs(1,:), xs(2,:))';
    h = 1/nC;
    intgrl = zeros(N,1);
    [lvar2, weight2] = GaussLegendre_2(15);
    lvar2 = (1+lvar2) * h/4; weight2 = weight2 * h/4; %[-1,1] -> [0, h/2]
    [lvar2x, lvar2y] = ndgrid(lvar2);
    lvar2x = lvar2x(:); lvar2y = lvar2y(:);
    weight2 = kron(weight2, weight2');
    weight2 = weight2(:);
    length_d = sqrt(lvar2x.^2+lvar2y.^2);
    for k = 1:N
        intgrl(k) = sum(weight2 .* (exp(-length_d.*mut(k))-1)./length_d);
    end
    intgrl = intgrl + h*log(1+sqrt(2));
    intgrl = -2/pi * intgrl;
    
    mus = musfun(xs(1,:), xs(2,:))';
    
    diagval = 1./mus + intgrl(:);
    
    param.Q = Q;
    param.n = nC;
    param.xs = xs;
    param.N = N;
    param.proxy = proxy;
    param.dim = 2;
    param.diagval = diagval;
    param.lvar = lvar;
    param.weight = weight;
    param.musfun = musfun;
    param.mutfun = mutfun;
end

