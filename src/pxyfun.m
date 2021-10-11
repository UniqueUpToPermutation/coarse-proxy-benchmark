function [Kpxy,nbr] = pxyfun(x,slf,nbr,l,ctr, param)
% proxy function
    N = param.N;
    proxy = param.proxy;
    
    pxy = bsxfun(@plus,proxy*l,ctr');
    Kpxy = Kfun(pxy,x(:,slf), param)/N;
    dx = x(1,nbr) - ctr(1);
    dy = x(2,nbr) - ctr(2);
    dist = sqrt(dx.^2 + dy.^2);
    nbr = nbr(dist/l < 1.5);
end