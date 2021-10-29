function [z_c,ctr] = GD(z_start, stepsize, maxiter, tol, theta_x, theta_y, k, z0, N, A)
%   Performs gradient descent minimization  of |det(I-K)| from starting
%   point z_start until |det(I-K)|<tol.
    D = 1e+10;
    discretization = 1e-6;
    z_c = z_start;
    ctr=0;
    while D>tol && ctr<maxiter
        [detK, grad] = numeric_gradient(z_c, discretization, theta_x, theta_y, z0, k, N, A);
        
        z_c = z_c - grad(1)*stepsize - 1i*grad(2)*stepsize;
        if detK>=D
            stepsize = 0.5*stepsize;
            discretization = 0.5*discretization;
        end
        ctr = ctr+1;
        D = detK;
     end
end