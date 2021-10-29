 function a = fourier_coefficient(k, n)
%     Computes the Fourier coefficients of the potential function for a
%     vector k of frequencies and integral quadrature size n.
    a = zeros(1,size(k,2));
    x = linspace(0,1,n);
    y = linspace(0,1,n);
    midpoints_x = zeros(1,n-1);
    midpoints_y = zeros(1,n-1);
    for j=1:n-1
        midpoints_x(j) = (x(j)+x(j+1))/2;
        midpoints_y(j) = (y(j)+y(j+1))/2;
    end
    [xx,yy] = meshgrid(midpoints_x, midpoints_y);
    vals = potential(xx,yy);
    for j=1:size(k,2)
        exponential = exp(2*pi*1i*(k(1,j).*xx + k(2,j).*yy));
        a(j) = sum(vals.*exponential, 'all')/n^2;
    end
end