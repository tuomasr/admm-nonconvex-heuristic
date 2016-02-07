% See Takapoui et al. (2015)
% ADMM
function [x_best, y_best] = admm_sparse(x, y, m, n, Q, c, alpha, A, E, F, b, no_cont, lb, ub, rho, e_tol)
    iter = 1;
    max_iter = 300;
    
    x_best = x;
    y_best = y;
    u = zeros(n+m,1);
    Qspa = sparse(Q);
    F2 = sparse(F^2);
    temp = sparse([Q+rho*F2 (A')*E;E*A -(1/rho)*eye(m)]);
    cache1 = sparse([eye(n) zeros(n,m)]/temp);
    cache2 = sparse((A')*(E^2)*b);
    cache3 = sparse([(A')*E F]);
    cache4 = sparse([zeros(n,m) inv(F)]);
    cache5 = sparse([E*A;F]);
    cache6 = sparse([zeros(m,n);F]);
    cache7 = sparse([E*b;zeros(n,1)]);
    cspa = sparse(c);
    
    while iter < max_iter
        xhalf = cache1*[-c+rho*(F2*x+cache2-cache3*u);zeros(m,1)];
        xnew = proj(xhalf+cache4*u,no_cont,lb,ub); 
        u = u + cache5*xhalf - cache6*xnew - cache7;
        x = xnew;
        y = x'*Qspa*x+cspa'*x+alpha;
        
        if all(((A*x-b).^2)<e_tol) && y < y_best
            x_best = x;
            y_best = y;
        end
        
        iter = iter + 1;
    end
end