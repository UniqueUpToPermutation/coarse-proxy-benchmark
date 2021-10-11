function K = Kfun(x,y, param)
    Q = param.Q;
    lvar = param.lvar;
    weight = param.weight;
    mutfun = param.mutfun;
    
    dx = bsxfun(@minus,x(1,:)',y(1,:));
    dy = bsxfun(@minus,x(2,:)',y(2,:));
    len1 = size(x,2);
    len2 = size(y,2);
    xinitial = repmat(x(1,:)',1,len2);
    yinitial = repmat(x(2,:)',1,len2);
    
    K = GaussLegendreQuadrature(@(l) mutfun(xinitial-l.*dx, yinitial-l.*dy), Q,lvar,weight);
    K = sqrt(dx.^2+dy.^2).*K;
    %LY:K = -1/(2*pi)*rdivide(exp(-K),sqrt( dx.^2+dy.^2 )) .* musfun(y(1,:), y(2,:));
    K = -1/(2*pi)*rdivide(exp(-K),sqrt( dx.^2+dy.^2 ));
end
