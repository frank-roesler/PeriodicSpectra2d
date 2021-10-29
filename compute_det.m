function Determinant = compute_det(window, k, theta_x, theta_y, z0, Id, sqrt_H0_inv, V, N)
%     Computes |det(I-K(z))| for all z in the window.
    Determinant = zeros(size(window));
%     Constant term |θ|^2-1:
    qm = (theta_x^2 + theta_y^2 - z0) * Id;
%     First part of potential term:  H_0^(-1/2)(|θ|^2−1+V):
    A = sqrt_H0_inv * (qm+V);
    parfor n=1:length(window(:))
        z = window(n);
        diag_grad = 2*(theta_x*k(1,:) + theta_y*k(2,:))./(z - z0 - vecnorm(k).^2);%Gradient term-2iθ∇:
        grad = spdiags(diag_grad.', 0, N^2, N^2);
        diag_resolvent = sqrt(z0 + vecnorm(k).^2)./(z - z0 - vecnorm(k).^2);
        sqrt_H0_resolvent = spdiags(diag_resolvent.', 0, N^2, N^2); % Resolvent:
        K = grad + A * sqrt_H0_resolvent; % Compact operator K(z,θ):
        Determinant(n) = det(Id - K); %*exp(trace(K));
    end
end