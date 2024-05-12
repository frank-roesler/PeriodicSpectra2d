function [detK, grad] = numeric_gradient(z_c, discretization, theta_x, theta_y, z0, k, N, A)
% Computes the finite difference gradient of |det(I-K)| at the point z_c.
    Id = speye(N^2);
    stencil = [z_c, z_c+discretization, z_c+1i*discretization];
    dets = zeros(1,length(stencil));

    for j=1:length(stencil)
        z = stencil(j);
        diag_grad = 2*(theta_x*k(1,:) + theta_y*k(2,:))./(z - z0 - vecnorm(k).^2); %Gradient term-2iθ∇:
        grad = spdiags(diag_grad.', 0, N^2, N^2);
        diag_resolvent = sqrt(z0 + vecnorm(k).^2)./(z - z0 - vecnorm(k).^2);
        sqrt_H0_resolvent = spdiags(diag_resolvent.', 0, N^2, N^2); % Resolvent:
        K = grad + A * sqrt_H0_resolvent; % Compact operator K(z,θ):
        dets(j) = det(Id - K);
    end
    grad = [(dets(2)-dets(1)), (dets(3)-dets(1))];
    grad = grad./discretization;
    detK = dets(1);
end