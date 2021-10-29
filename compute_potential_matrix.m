function V = compute_potential_matrix(k,a,N)
%     Assembles the NxN-matrix V_n^appr made up of the Fourier coefficients
%     given in the list a.
    W = zeros(N^2);
    for j=1:size(W,1)
        for m=1:size(W,2)
            W(j,m) = a(all(k==k(:,m)-k(:,j)));
        end
    end
    W(abs(W)<1e-6) = 0; % Band limit to make V sparse. Remove this for higher accuracy.
    W = sparse(W);
    V = W;
end
    