% MAIN SCRIPT: Computes spectrum of periodic Schrödinger operator with
% non-real potential in 2d.

%% Setup: 
clear;
tic

%  Define potential matrix:
N   = 10;  

disp('Building potential matrix...')
n_plot = 1000; 
k = build_k(2*N);
a = fourier_coefficient(k, n_plot);
V = compute_potential_matrix(k, a, N);
V = round(V);
disp('Done!')
pause(0.1)

%  Define the vector of k's containing integer multiples of 2π:
k = 2*pi*k;

resolution_z = 4.5;             % resolution in complex plane
resolution = 20;                % theta resolution
                                % Size of matrix K
n_t = round(2*pi*resolution);   % number of theta points
z0  = -4;                       % Spectral shift

r = 0.28;

t1 = linspace(0,2*pi,n_t); 
t2 = linspace(0,2*pi,n_t);
[tx,ty]  = meshgrid(t1, t2); % theta array
tx = tx.';
ty = ty.';

Id = speye(N^2);

% Define lattice in C:
r1 = -2;
r2 = 35;
i1 = -3;
i2 = 3;
M1 = round((r2-r1)*resolution_z);
M2 = round((i2-i1)*resolution_z); 

X = linspace(r1, r2, M1);
Y = linspace(i1, i2, M2);
[XX,YY] = meshgrid(X,Y);
L = XX+1i*YY;
h_L = 1/resolution_z;

% Inverse square root of H0:
sqrt_H0_inv = spdiags(1./sqrt(z0 + vecnorm(k).^2).', 0, N^2,N^2);

% Compute initial set of zeros:
disp('Computing starting points...')
Spectrum = [];
ctr=0;

while isempty(Spectrum)
    ctr = ctr+1;
    theta_x=tx(ctr);
    theta_y=ty(ctr);
    Determinant = compute_det(L, k, theta_x, theta_y, z0, Id, sqrt_H0_inv, V, N);
    mins_x = islocalmin(abs(Determinant),1, 'FlatSelection','all');
    mins_y = islocalmin(abs(Determinant),2, 'FlatSelection','all');
    mins_xy = mins_x & mins_y ; 
    Spectrum = unique([Spectrum, L(mins_xy).']);
end

Spectrum_temp = Spectrum;
Spectrum_temp_old = Spectrum_temp;
window = L(box(L, Spectrum_temp, r, 0));
disp('Done!')
% 

%% Main loop:

figure('Position',[100 600 1500 200])
tic

disp('Computing...')
save_ctr = 0;
for m=ctr:length(tx(:))
    theta_x = tx(m);
    theta_y = ty(m);

%     Constant term |θ|^2-1:
    qm = (theta_x^2 + theta_y^2 - z0) * Id;
%     First part of potential term:  H_0^(-1/2)(|θ|^2−1+V):
    A = sqrt_H0_inv * (qm+V);
    
    starting_points = window;
    
%     % Plot Results:
    plot_results(L, Spectrum, starting_points, window, [])
    
    endpoints = zeros(1,length(starting_points(:)));
    Determinant = zeros(length(starting_points(:)));
%     Start gradient descent:
    iters = length(endpoints);
    maxiter = 1000;
    stepsize = 0.04;
    tol = 1e-4;
    while max(abs(endpoints))==0
        parfor n=1:iters
            z_start = starting_points(n);
            [z_c,ctr] = GD(z_start, stepsize, maxiter, tol, theta_x, theta_y, k, z0, N, A);
            if ctr<maxiter
                endpoints(n) = z_c;
            end
        end
        maxiter = maxiter*5;
    end
    endpoints = endpoints(endpoints~=0);
    
    Spectrum_temp = endpoints(real(endpoints)<r2+1);
    
    if isempty(Spectrum_temp)
        window = L(box(L, bdry_pts(Spectrum_temp_old, h_L), r,0));
    else
        window = L(box(L, Spectrum_temp, r,0));
    end
    Spectrum_temp_old = Spectrum_temp;
    Spectrum = unique([Spectrum,Spectrum_temp]);
    
    if theta_x==0 && m>2
        plot(Spectrum+0.000001i,'.','MarkerEdgeColor',[0.3 0.3 1], 'MarkerSize',10);
        xlim([min(real(L(:))) max(real(L(:)))]);
        ylim([min(imag(L(:))) max(imag(L(:)))]);
        drawnow;
        
        save_ctr = save_ctr+1;
        save(['Spectrum_',num2str(save_ctr)], 'Spectrum','L','m','n_t','resolution','N')
        Spectrum_temp = [Spectrum_temp,Spectrum(1:10)];
        Spectrum = [];
        disp('-----------------------------------------')
        disp([num2str(round(theta_y/2/pi*100,1)), '% completed'])
        disp([num2str(toc/60),' ',' minutes'])
        tic
        disp('-----------------------------------------')
    end
end
disp('Done!')

