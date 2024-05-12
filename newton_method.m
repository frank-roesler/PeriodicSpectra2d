function [z_c,ctr] = newton_method(z_start, stepsize, maxiter, tol, theta_x, theta_y, k, z0, N, A)
    D = inf;
    discretization = 1e-12;
    z_c = z_start;
    ctr=1;
    while D>tol && ctr<maxiter
        [detK, grad] = numeric_gradient(z_c, discretization, theta_x, theta_y, z0, k, N, A);
        dDdz = 0.5*(grad(1)-1i*grad(2));
        z_c = z_c - detK/dDdz;
        D = abs(detK);
        ctr = ctr+1;
    end
end
    